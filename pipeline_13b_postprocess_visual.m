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
        if length(strfind(session,'ses-'))
            session_dir=fullfile(subject_dir, session);
            disp(session_dir);
            
            ses_coherence={};
            ses_congruence={};
            beh_files=dir(fullfile(session_dir,'autoreject*-visual-beh.csv'));
            for b=1:length(beh_files)
                beh_file=beh_files(b).name;
                beh_df=readtable(fullfile(session_dir,beh_file));
                trial_coherence=beh_df.trial_coherence;
                trial_congruence=beh_df.trial_congruence;
                
                ses_coherence = [ses_coherence; trial_coherence];
                ses_congruence = [ses_congruence; trial_congruence];
            end
            
            spm_dir=fullfile(session_dir, 'spm');
            
            spm_file=fullfile(spm_dir, sprintf('cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            D=spm_eeg_load(spm_file);
            
            % Congruence
            cong_trials=find(strcmp(ses_congruence,'congruent'));
            incong_trials=find(strcmp(ses_congruence,'incongruent'));
            
            cong_file=fullfile(spm_dir, sprintf('congruent_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dcong = subset_trials(D, cong_file, cong_trials);
            
            incong_file=fullfile(spm_dir, sprintf('incongruent_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dincong = subset_trials(D, incong_file, incong_trials);
            
            % Coherence
            low_trials=find(strcmp(ses_coherence,'low'));
            med_trials=find(strcmp(ses_coherence,'med'));
            high_trials=find(strcmp(ses_coherence,'high'));
            
            low_file=fullfile(spm_dir, sprintf('coherence-low_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dlow = subset_trials(D, low_file, low_trials);
            
            med_file=fullfile(spm_dir, sprintf('coherence-med_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dmed = subset_trials(D, med_file, med_trials);
            
            high_file=fullfile(spm_dir, sprintf('coherence-high_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dhigh = subset_trials(D, high_file, high_trials);
            
            % Congruent coherence
            low_trials=find(strcmp(ses_congruence,'congruent') & strcmp(ses_coherence,'low'));
            med_trials=find(strcmp(ses_congruence,'congruent') & strcmp(ses_coherence,'med'));
            high_trials=find(strcmp(ses_congruence,'congruent') & strcmp(ses_coherence,'high'));
            
            cong_low_file=fullfile(spm_dir, sprintf('congruent_coherence-low_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dconglow = subset_trials(D, cong_low_file, low_trials);
            
            cong_med_file=fullfile(spm_dir, sprintf('congruent_coherence-med_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dcongmed = subset_trials(D, cong_med_file, med_trials);
            
            cong_high_file=fullfile(spm_dir, sprintf('congruent_coherence-high_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dconghigh = subset_trials(D, cong_high_file, high_trials);
            
            % Incongruent coherence
            low_trials=find(strcmp(ses_congruence,'incongruent') & strcmp(ses_coherence,'low'));
            med_trials=find(strcmp(ses_congruence,'incongruent') & strcmp(ses_coherence,'med'));
            high_trials=find(strcmp(ses_congruence,'incongruent') & strcmp(ses_coherence,'high'));
            
            incong_low_file=fullfile(spm_dir, sprintf('incongruent_coherence-low_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dinconglow = subset_trials(D, incong_low_file, low_trials);
            
            incong_med_file=fullfile(spm_dir, sprintf('incongruent_coherence-med_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dincongmed = subset_trials(D, incong_med_file, med_trials);
            
            incong_high_file=fullfile(spm_dir, sprintf('incongruent_coherence-high_cspm_converted_autoreject-%s-%s-visual-epo.mat', subject, session));
            Dinconghigh = subset_trials(D, incong_high_file, high_trials);
            
            
            condition_files={cong_file, incong_file, low_file, med_file, high_file, cong_low_file, cong_med_file, cong_high_file, incong_low_file, incong_med_file, incong_high_file};
            
            for f=1:length(condition_files)            
                spm_jobman('initcfg');

                clear jobs
                matlabbatch={};
                        
                matlabbatch{1}.spm.meeg.averaging.average.D(1) = {condition_files{f}};
                matlabbatch{1}.spm.meeg.averaging.average.userobust.standard = false;
                matlabbatch{1}.spm.meeg.averaging.average.plv = false;
                matlabbatch{1}.spm.meeg.averaging.average.prefix = 'm';

                matlabbatch{2}.spm.meeg.preproc.crop.D(1) = cfg_dep('Averaging: Averaged Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
                matlabbatch{2}.spm.meeg.preproc.crop.timewin = [-500 1000];
                matlabbatch{2}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
                matlabbatch{2}.spm.meeg.preproc.crop.channels{1}.all = 'all';
                matlabbatch{2}.spm.meeg.preproc.crop.prefix = 'rdk_p';

                matlabbatch{3}.spm.meeg.preproc.crop.D(1) = cfg_dep('Averaging: Averaged Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
                %matlabbatch{3}.spm.meeg.preproc.crop.timewin = [2250 2750];
                matlabbatch{3}.spm.meeg.preproc.crop.timewin = [1500 3000];
                matlabbatch{3}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
                matlabbatch{3}.spm.meeg.preproc.crop.channels{1}.all = 'all';
                matlabbatch{3}.spm.meeg.preproc.crop.prefix = 'instr_p';

                spm_jobman('run', matlabbatch);
            end
        end
    end
end
