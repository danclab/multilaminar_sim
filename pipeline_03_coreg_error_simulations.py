import glob
import shutil
import sys
import time

import numpy as np
import nibabel as nib
import pickle
from scipy.spatial.transform import Rotation as R
from lameg.invert import coregister, invert_ebb
from lameg.laminar import model_comparison
from lameg.simulate import run_current_density_simulation
from lameg.util import spm_context
import json
import os
import os.path as op

from utilities.utils import gaussian, get_fiducial_coords


def run(subject_id, session_id, sim_vertex, err_level, win_size, json_file):

    with open(json_file) as pipeline_file:
        parameters = json.load(pipeline_file)

    out_path = os.path.join(parameters['output_path'], 'coreg_err_simulations')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    output_file = os.path.join(
        out_path,
        f"vx_{sim_vertex}_err_{err_level}_winsize_{win_size}.pickle"
    )

    if not os.path.exists(output_file):
        path = parameters["dataset_path"]
        der_path = op.join(path, "derivatives")
        proc_path = op.join(der_path, "processed")

        print("ID:", subject_id)

        sub_path = op.join(proc_path, subject_id)

        ses_path = op.join(sub_path, session_id)

        lpa, nas, rpa = get_fiducial_coords(subject_id, json_file)

        orig_fid = np.array([np.array(nas), np.array(lpa), np.array(rpa)])
        mean_fid = np.mean(orig_fid, axis=0)
        zero_mean_fid = np.hstack([(orig_fid - mean_fid), np.ones((3, 1))])

        # Patch size to use for inversion
        patch_size = 5
        # Number of temporal modes
        n_temp_modes = 4
        # Window of interest
        woi = [-int(.5*win_size), int(.5*win_size)]

        tmp_dir = op.join(out_path, f'sim_coregerr_{sim_vertex}_{err_level}_{win_size}')
        if not op.exists(tmp_dir):
            os.mkdir(tmp_dir)

        # Native space MRI to use for coregistration
        shutil.copy(os.path.join(sub_path, 't1w.nii'),
                    os.path.join(tmp_dir, 't1w.nii'))
        surf_files = glob.glob(os.path.join(sub_path, 't1w*.gii'))
        for surf_file in surf_files:
            shutil.copy(surf_file, tmp_dir)
        mri_fname = os.path.join(tmp_dir, 't1w.nii')

        surf_dir = os.path.join(sub_path, 'surf')

        n_layers = 11
        SNR = 0

        # Get name of each mesh that makes up the layers of the multilayer mesh - these will be used for the source
        # reconstruction
        orientation_method = 'link_vector.fixed'

        shutil.copy(os.path.join(surf_dir,'multilayer.11.ds.link_vector.fixed.gii'),
                    os.path.join(tmp_dir,'multilayer.11.ds.link_vector.fixed.gii'))
        multilayer_mesh_fname = os.path.join(tmp_dir, 'multilayer.11.ds.link_vector.fixed.gii')
        shutil.copy(os.path.join(surf_dir, 'FWHM5.00_multilayer.11.ds.link_vector.fixed.mat'),
                    os.path.join(tmp_dir, 'FWHM5.00_multilayer.11.ds.link_vector.fixed.mat'))

        layers = np.linspace(1, 0, n_layers)
        layer_fnames = []
        for layer in layers:

            if layer == 1:
                name = f'pial.ds.{orientation_method}'
            elif layer == 0:
                name = f'white.ds.{orientation_method}'
            else:
                name = f'{layer:.3f}.ds.{orientation_method}'

            shutil.copy(os.path.join(surf_dir, f'{name}.gii'),
                        os.path.join(tmp_dir, f'{name}.gii'))
            layer_fnames.append(os.path.join(tmp_dir, f'{name}.gii'))
            shutil.copy(os.path.join(surf_dir, f'FWHM5.00_{name}.mat'),
                        os.path.join(tmp_dir, f'FWHM5.00_{name}.mat'))

        data_file = os.path.join(ses_path,
            f'spm/pcspm_converted_autoreject-{subject_id}-{session_id}-motor-epo.mat'
        )
        data_path, data_file_name = os.path.split(data_file)
        data_base = os.path.splitext(data_file_name)[0]

        # Copy data files to tmp directory
        shutil.copy(
            os.path.join(data_path, f'{data_base}.mat'),
            os.path.join(tmp_dir, f'{data_base}.mat')
        )
        shutil.copy(
            os.path.join(data_path, f'{data_base}.dat'),
            os.path.join(tmp_dir, f'{data_base}.dat')
        )

        # Construct base file name for simulations
        base_fname = os.path.join(tmp_dir, f'{data_base}.mat')


        dipole_moment = 10

        # simulation signal
        time = np.linspace(-1, 1, num=1201)
        sim_signal = gaussian(time, 0.0, 0.4, dipole_moment)

        mesh = nib.load(multilayer_mesh_fname)
        verts_per_surf = int(mesh.darrays[0].data.shape[0] / n_layers)

        sim_vx_res = {}
        sim_vx_res["woi"] = woi
        sim_vx_res["sim_vertex"] = sim_vertex
        sim_vx_res["time"] = time
        sim_vx_res["err_level"] = err_level
        sim_vx_res["hann"] = True
        sim_vx_res["signal"] = sim_signal
        sim_vx_res['layerF'] = np.zeros((n_layers, n_layers))

        with spm_context() as spm:
            coregister(
                nas, lpa, rpa, mri_fname,
                multilayer_mesh_fname,
                base_fname, spm_instance=spm,
                viz=False
            )

            [_, _] = invert_ebb(
                multilayer_mesh_fname, base_fname,
                n_layers, patch_size=patch_size,
                n_temp_modes=n_temp_modes,
                spm_instance=spm,
                viz=False
            )

            while True:
                # Translation vector
                translation = err_level
                shift_vec = np.random.randn(3)
                shift_vec = shift_vec / np.linalg.norm(shift_vec) * np.random.randn() * translation

                # Rotation vector
                rotation_rad = err_level * np.pi / 180.0
                rot_vec = np.random.randn(3)
                rot_vec = rot_vec / np.linalg.norm(rot_vec) * np.random.randn() * rotation_rad

                # Apply transformation to fiducial locations
                P = np.concatenate((shift_vec, rot_vec))
                rotation = R.from_rotvec(P[3:])
                rotation_matrix = rotation.as_matrix()
                A = np.eye(4)
                A[:3, :3] = rotation_matrix
                A[:3, 3] = P[:3]

                # Transform zero_mean_fid
                new_fid_homogeneous = (A @ zero_mean_fid.T).T
                new_fid = new_fid_homogeneous[:, :3] + mean_fid
                max_dist = np.max(np.sqrt(np.sum((new_fid - orig_fid) ** 2, axis=-1)))
                if np.abs(err_level - max_dist) < .05:
                    break

            sim_vx_res['new_nas'] = new_fid[0, :]
            sim_vx_res['new_lpa'] = new_fid[1, :]
            sim_vx_res['new_rpa'] = new_fid[2, :]

            for l in range(n_layers):
                prefix = f'sim_{sim_vertex}_layer{str(l).zfill(2)}_'
                l_vertex = l * verts_per_surf + sim_vertex
                sim_fname = run_current_density_simulation(
                    base_fname,
                    prefix,
                    l_vertex,
                    sim_signal,
                    dipole_moment,
                    patch_size,
                    SNR,
                    average_trials=True,
                    spm_instance=spm
                )

                [layerF, _] = model_comparison(
                    new_fid[0, :],
                    new_fid[1, :],
                    new_fid[2, :],
                    mri_fname,
                    layer_fnames,
                    sim_fname,
                    method='MSP',
                    viz=False,
                    spm_instance=spm,
                    invert_kwargs={
                        'priors': [sim_vertex],
                        'woi': woi,
                        'hann_windowing': True,
                        'patch_size': patch_size,
                        'n_temp_modes': n_temp_modes
                    }
                )
                sim_vx_res['layerF'][:,l] = layerF

                extensions = ['mat', 'dat']
                for ext in extensions:
                    os.remove(
                        os.path.join(
                            tmp_dir, f'{prefix}pcspm_converted_autoreject-{subject_id}-{session_id}-motor-epo.{ext}'
                        )
                    )
                    os.remove(
                        os.path.join(
                            tmp_dir, f'm{prefix}pcspm_converted_autoreject-{subject_id}-{session_id}-motor-epo.{ext}'
                        )
                    )


        with open(output_file, "wb") as fp:
            pickle.dump(sim_vx_res, fp)

        shutil.rmtree(tmp_dir)

if __name__=='__main__':
    # parsing command line arguments
    try:
        sim_idx = int(sys.argv[1])
    except:
        print("incorrect simulation index")
        sys.exit()

    try:
        json_file = sys.argv[2]
        print("USING:", json_file)
    except:
        json_file = "settings.json"
        print("USING:", json_file)

    with open(json_file) as pipeline_file:
        parameters = json.load(pipeline_file)

    subject_id = 'sub-001'
    session_id = 'ses-01'
    path = parameters["dataset_path"]
    der_path = op.join(path, "derivatives")
    proc_path = op.join(der_path, "processed")
    sub_path = op.join(proc_path, subject_id)
    surf_dir = op.join(sub_path,'surf')
    mesh_fname = op.join(surf_dir, 'pial.ds.link_vector.fixed.gii')
    mesh = nib.load(mesh_fname)

    n_vertices = mesh.darrays[0].data.shape[0]
    np.random.seed(42)
    vertices = np.random.randint(0, n_vertices, 100)
    err_levels = [0, 0.5, 1, 2, 3, 4, 5]
    win_sizes = [50]

    all_verts = []
    all_err_levels = []
    all_win_sizes = []
    for vert in vertices:
        for err in err_levels:
            for win_size in win_sizes:
                all_verts.append(vert)
                all_err_levels.append(err)
                all_win_sizes.append(win_size)
    np.random.seed(int(time.time()))

    vertex_idx = all_verts[sim_idx]
    err_level = all_err_levels[sim_idx]
    win_size = all_win_sizes[sim_idx]

    run(subject_id, session_id, vertex_idx, err_level, win_size, json_file)
