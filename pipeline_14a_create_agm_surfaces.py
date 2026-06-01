import nump as np
from lameg.surf import postprocess_freesurfer_surfaces
import spm_standalone

postprocess_freesurfer_surfaces(
    'sub-001',
    '/home/common/bonaiuto/cued_action_meg/derivatives/processed/sub-001/agm_surf/',
    'multilayer.11.ds.link_vector.fixed.gii',
    n_surfaces=11,
    ds_factor=0.11,
    orientation='link_vector',
    fix_orientation=True,
    remove_deep=True,
    n_jobs=-1
)

spm = spm_standalone.initialize()

_ = spm.spm_eeg_smoothmesh_multilayer_mm(
    '/home/common/bonaiuto/cued_action_meg/derivatives/processed/sub-001/agm_surf/multilayer.11.ds.link_vector.fixed.gii',
    float(5.0),
    float(11),
    nargout=1
)

layers = np.linspace(1, 0, 11)
layer_fnames=[]
for layer in layers:
    if layer == 1:
        layer_fnames.append(
            '/home/common/bonaiuto/cued_action_meg/derivatives/processed/sub-001/agm_surf/pial.ds.link_vector.fixed.gii'
        )
    elif 0 < layer < 1:
        layer_fnames.append(
            f'/home/common/bonaiuto/cued_action_meg/derivatives/processed/sub-001/agm_surf/{layer:.3f}.ds.link_vector.fixed.gii'
        )
    elif layer == 0:
        layer_fnames.append(
            '/home/common/bonaiuto/cued_action_meg/derivatives/processed/sub-001/agm_surf/white.ds.link_vector.fixed.gii'
        )
print(layer_fnames)

for layer_fname in layer_fnames:
     _ = spm.spm_eeg_smoothmesh_multilayer_mm(
        layer_fname,
        float(5.0),
        float(1),
        nargout=1
    )

spm.terminate()