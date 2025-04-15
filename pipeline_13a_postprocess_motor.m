json_file = "settings.json";
parameters = read_json_file(json_file);

base_dir=fullfile(parameters.dataset_path, 'derivatives/processed');

spm('defaults','eeg');

subj_dirs=dir(base_dir);
for s=3:length(subj_dirs)
    subject=subj_dirs(s).name;
    subject_dir=fullfile(base_dir, subject);
    ses_dirs=dir(subject_dir);
    for se=3:length(ses_dirs)
        session=ses_dirs(se).name;
        if length(strfind(session,'ses'))
            session_dir=fullfile(subject_dir, session);
            disp(session_dir);
            
            spm_jobman('initcfg');

            clear jobs
            matlabbatch={};
            
            spm_dir=fullfile(session_dir, 'spm');
            run_files={};
            spm_files=dir(spm_dir);
            for f=3:length(spm_files)
                if strcmp(spm_files(f).name(1:3),'spm') && length(strfind(spm_files(f).name,'-motor-epo.mat'))
                    run_files{end+1}=fullfile(spm_dir, spm_files(f).name);
                    disp(spm_files(f).name);
                end
            end
            if length(run_files)>1
                matlabbatch{1}.spm.meeg.preproc.merge.D = run_files';
                matlabbatch{1}.spm.meeg.preproc.merge.recode.file = '.*';
                matlabbatch{1}.spm.meeg.preproc.merge.recode.labelorg = '.*';
                matlabbatch{1}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
                matlabbatch{1}.spm.meeg.preproc.merge.prefix = 'c';
                
                merged_file=fullfile(spm_dir, sprintf('cspm_converted_autoreject-%s-%s-motor-epo.mat',subject,session));
                matlabbatch{2}.spm.meeg.other.copy.D(1) = cfg_dep('Merging: Merged Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
                matlabbatch{2}.spm.meeg.other.copy.outfile = merged_file;

                matlabbatch{3}.spm.meeg.preproc.crop.D(1) = cfg_dep('Copy: Copied M/EEG datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
                matlabbatch{3}.spm.meeg.preproc.crop.timewin = [-1000 1000];
                matlabbatch{3}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
                matlabbatch{3}.spm.meeg.preproc.crop.channels{1}.all = 'all';
                matlabbatch{3}.spm.meeg.preproc.crop.prefix = 'p';
                
                matlabbatch{4}.spm.meeg.averaging.average.D(1) = cfg_dep('Copy: Copied M/EEG datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
                matlabbatch{4}.spm.meeg.averaging.average.userobust.standard = false;
                matlabbatch{4}.spm.meeg.averaging.average.plv = false;
                matlabbatch{4}.spm.meeg.averaging.average.prefix = 'm';
                
                matlabbatch{5}.spm.meeg.preproc.crop.D(1) = cfg_dep('Averaging: Averaged Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
                matlabbatch{5}.spm.meeg.preproc.crop.timewin = [-1000 1000];
                matlabbatch{5}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
                matlabbatch{5}.spm.meeg.preproc.crop.channels{1}.all = 'all';
                matlabbatch{5}.spm.meeg.preproc.crop.prefix = 'p';

                spm_jobman('run', matlabbatch);
            end
        end
    end
end