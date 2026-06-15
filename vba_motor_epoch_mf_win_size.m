function vba_motor_epoch_mf_win_size()

win_sizes=[10,25,50];

for w_idx=1:length(win_sizes)
    win_size=win_sizes(w_idx);
    in_fname=sprintf('mf_aligned_all_layer_F_win_size_%d.mat', win_size);
    load(fullfile('./output', in_fname));

    F = all_layer_F


    % ========================
    % PARAMETERS
    % ========================
    dt = (0.5 - (-0.5)) / (601 - 1);  
    time = linspace(-0.5, 0.5, 601);

    options = struct();
    options.niter = 100
    options.DisplayWin = 0; 
    options.families = {[1, 2, 3], [4], [5, 6]} ;
    nSubjects   = size(F,1);   % 8 in ERF
    nModels     = size(F,2);   % 6 layers
    nTimepoints = size(F,3);   % times

    EP  = zeros(6, nTimepoints);
    famEP = zeros(3, nTimepoints);
    p_H0 = zeros(1, nTimepoints);           % BOR in time
    PEP  = zeros(3, nTimepoints);
    a_model  = zeros(6, nTimepoints); %


    for t = 1:nTimepoints

        % models × subjects
        F_t = squeeze(F(:,:,t)).';

        [posterior, out] = VBA_groupBMC(F_t, options);


        % Bayesian Omnibus
        p_H0(t) = out.bor; %null: p(H0|y) 

        % Protected Exceedance Probabilities
        K = length(out.families.ep); %3 family 
        EP(:,t) = out.ep; %ep for 6 layers 
        famEP(:,t) = out.families.ep; %ep for 3 family
        PEP(:,t) = (1 - out.bor) .* out.families.ep + out.bor / K; %pep for family bs our out.bor is against 3 family

    end


    bestModelOverTime = zeros(1, nTimepoints);
    for t = 1:nTimepoints

        if p_H0(t) < 0.5 % if it is more chance to have one dominant family 


            [~, bestModel] = max(PEP(:,t)); %  it can be er also here 
            bestModelOverTime(t) = bestModel;

        end

    end

    out_fname=sprintf('mf_aligned_overtime_family_win_size_%d.mat', win_size);
    save(fullfile('./output', out_fname), 'bestModelOverTime', 'time', 'p_H0', 'EP', 'famEP');


    %quick check 

    figure;
    plot(time, p_H0, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('p(H0|y)');
    title('Bayesian Omnibus Risk over Time');
    grid on;

    figure;
    imagesc(time, 1:3, PEP);
    set(gca,'YDir','normal');
    colorbar;
    xlabel('Time (s)');
    ylabel('Model');
    title('Protected Exceedance Probabilities (PEP)');

end