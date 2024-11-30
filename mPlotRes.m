%% mPlotRes: loads and plots the MMSlab2 results
% loads in saved data from MMSlab2 and shows combined plots of measured
% data. The script contains the basic structure to load an plot data, but
% requires additions and alterations to answer the questions to PA1 of the
% Human Controller correctly.
%
%
% intended for:      PA1 (MMSlab) of The Human Controller
% created by:        R.J. Kuiper
% last modified on:  29-4-2015
%-------------------------------------------------------------------------

clear all

%% define data files (edit appropriately)

% defining the relative path to the data files

fpath.sub{1} = './data_ours/data/maxdam/';      
fpath.sub{2} = './data_ours/data/lefteris_mourelatos/';  
fpath.sub{3} = './data_ours/data/tobiasknell/';  


% defining the file names for each task and repetition

fnames.sub(1).task(:,:) = {'Est_dataA1','Est_dataA2';...
                           'Est_dataB1','Est_dataB2';...
                           'Est_dataC1','Est_dataC2'};
fnames.sub(2).task(:,:) = {'Est_dataA1','Est_dataA2';...
                           'Est_dataB1','Est_dataB2';...
                           'Est_dataC1','Est_dataC2'};
fnames.sub(3).task(:,:) = {'Est_dataA1','Est_dataA2';...
                           'Est_dataB1','Est_dataB2';...
                           'Est_dataC1','Est_dataC2'};


% add here additional tasks and subjects to the name struct 'fnames'. You
% can add tasks by enlarging the written example with rows. Addind subjects
% is best done by copying the struct and repeat it for each subject and
% changing the index number 'sub(n)'
                       
%% read in data from files (no need to edit)

% define total number of subjects, tasks and repetitions per task
[nsub, dummy] = size(fnames.sub(:));                                    % number of subjects
[ntsk, nrep]  = size(fnames.sub(1).task);                               % number of tasks and repetitions

% read in all data files and extract variables into a single data struct
clear data
for s = 1:nsub
    for t = 1:ntsk
        for r = 1:nrep
            run([fpath.sub{s}, fnames.sub(s).task{t,r}])              % run selected data file in slected path

            data.sub(s).task(t).W(:,r)   = Wdata1_X0(:,1);              % frequency vector
            data.sub(s).task(t).Mp(:,r)  = Wdata1_X0(:,2);              % magnitude pilot vector
            data.sub(s).task(t).Ap(:,r)  = Wdata1_X0(:,3);              % phase angle pilot vector
            data.sub(s).task(t).Mol(:,r) = Wdata1_X0(:,4);              % magnitude open-loop system vector
            data.sub(s).task(t).Aol(:,r) = Wdata1_X0(:,5);              % phase angle open-loop system vector
            
            data.sub(s).task(t).Snum(1:length(num_Xs0),r) = num_Xs0;    % numerator system transfer function
            data.sub(s).task(t).Sden(1:length(den_Xs0),r) = den_Xs0;    % denumerator system transfer function
            
            data.sub(s).task(t).rms_err(r) = rms_errX0;                 % RMS of system error vector
            data.sub(s).task(t).std_inp(r) = std_inpX0;                 % STD of pilot input vector
            
            clear Wdata1_X0 num_Xs0 den_Xs0 rms_errX0 std_inpX0         % clear variables
        end
    end
end

%% calculate mean and standard deviation for bode plot data
% (edit when necessary for inter/between subject analysis)

% calculate mean and standard deviation of all repetitions
all_stats = [];

for s = 1:nsub
    for t = 1:ntsk
        Mp  = data.sub(s).task(t).Mp';                      % magnitude of pilot vector
        Mol = data.sub(s).task(t).Mol';                     % magnitude of system open-loop vector
        Ap  = data.sub(s).task(t).Ap';                      % phase angle of pilot vector
        Aol = data.sub(s).task(t).Aol';                     % phase angle of system open-loop vector

        data.sub(s).task(t).Mpmn   = mean( Mp(:,:) );       % mean of pilot magnitudes for all repetitions
        data.sub(s).task(t).Mpstd  = std(  Mp(:,:) );       % standard deviation of pilot magn. for all reps.
        data.sub(s).task(t).Molmn  = mean( Mol(:,:) );      % mean of open-loop magnitudes for all repetitions
        data.sub(s).task(t).Molstd = std(  Mol(:,:) );      % standard deviation of open-loop magn. for all reps.
        data.sub(s).task(t).Apmn   = mean( Ap(:,:) );       % mean of pilot phase angles for all repetitions
        data.sub(s).task(t).Apstd  = std(  Ap(:,:) );       % standard deviation of pilot phase angles for all reps.
        data.sub(s).task(t).Aolmn  = mean( Aol(:,:) );      % mean of open-loop phase angles for all repetitions
        data.sub(s).task(t).Aolstd = std(  Aol(:,:) );      % standard deviation of open-loop phase angles for all reps.

    end
end

task_means = struct();
for t=1:ntsk
    all_Mp = [];
    all_Mol = [];
    all_Ap = [];
    all_Aol = [];
    for s=1:nsub
        all_Mp = [all_Mp; data.sub(s).task(t).Mp'];
        all_Mol = [all_Mol; data.sub(s).task(t).Mol'];
        all_Ap = [all_Ap; data.sub(s).task(t).Ap'];
        all_Aol = [all_Aol; data.sub(s).task(t).Aol'];
    end
    task_means.task(t).all_Mp = all_Mp;
    task_means.task(t).all_Mol = all_Mol;
    task_means.task(t).all_Ap = all_Ap;
    task_means.task(t).all_Aol = all_Aol;
end



%% Select data for showing in bode plots
% alter this code when required, e.g. for between subject variability. But
% use this code as example of extracting data from the struct and plotting
% the data in bode plots

s = 1;                                        % select subject for bode plots
t = 1;                                        % select task for bode plots

% extract data for bode plot from data struct for the selected task and subject
W = data.sub(s).task(t).W(:,1)/(2*pi);   % frequency vector in Hz
Snum   = data.sub(s).task(t).Snum(:,1)';      % system transfer function numerator
Sden   = data.sub(s).task(t).Sden(:,1)';      % system transfer function denumerator

Mpmn = mean(task_means.task(t).all_Mp(:,:));
Mpstd = std(task_means.task(t).all_Mp(:,:));
Molmn = mean(task_means.task(t).all_Mol(:,:));
Molstd = std(task_means.task(t).all_Mol(:,:));
Apmn = mean(task_means.task(t).all_Ap(:,:));
Apstd = std(task_means.task(t).all_Ap(:,:));
Aolmn = mean(task_means.task(t).all_Aol(:,:));
Aolstd = std(task_means.task(t).all_Aol(:,:));

% Mpmn   = data.sub(s).task(t).Mpmn;            % pilot mean magniture
% Mpstd  = data.sub(s).task(t).Mpstd;           % pilot st. dev. magniture
% Molmn  = data.sub(s).task(t).Molmn;           % open-loop mean magniture
% Molstd = data.sub(s).task(t).Molstd;          % open-loop st, dev. magniture
% Apmn   = data.sub(s).task(t).Apmn;            % pilot mean phase angle
% Apstd  = data.sub(s).task(t).Apstd;           % pilot st. dev. phase angle
% Aolmn  = data.sub(s).task(t).Aolmn;           % open-loop mean phase angle
% Aolstd = data.sub(s).task(t).Aolstd;          % open-loop st. dev. phase angle

            
Hsys = tf(Snum,Sden);                         % transfer function of the measured system
% figure(1); clf; step(Hsys,10)               % uncomment this to show the step response of the selected system
[mag,phase] = bode(Hsys,W*2*pi);              % bodeplot of the system transfer function
Msys = squeeze(mag(1,1,:));                   % magnitude of the system bodeplot
Asys = squeeze(phase(1,1,:));                 % phase angle of the system bodeplot

% HINT: to determine the crossover frequency you can use the function
%   'find' to look for a specific point in a vector. Using the function
%   'interp' will allow to you to interpolate between the limited datya point
%   in the vector.
% Here is an example (so not giving you the answer...!):
%   W2   = interp(W,100);                     % interpolates the vector with 100 points between each sample
%   iw2  = find(W2 > 1.0, 1, 'first');        % finds the first index of the vector > 1.0
%   fprintf('%5.3f [Hz] \n',W2(iw2))          % prints the found vector value

% example of McRuer estimation of pilot model, change accordingly for the selected type of task
Hmcr = tf([4 0],[1 0])                        % McRuer estimate of pilot; zeroth order system with a gain of 4
% figure(2); clf; step(Hhum,10)               % uncomment this to show the step response of the selected system
[mag,phase] = bode(Hmcr,W);                   % bodeplot of the pilot transfer function
Mmcr = squeeze(mag(1,1,:));                   % magnitude of the pilot bodeplot
Amcr = squeeze(phase(1,1,:));                 % phase angle of the pilot bodeplot

s_fd = tf('s')                                % matlab transfer function s as building block for transfer functions

% Pilot models per task, as provided in McRuers table II

if t == 1
    Kp = 4;          % Pilot gain
    Te = 0.5;          % Time delay (tau_e)    

    Hmcr = Kp * exp(-Te * s_fd)

elseif t == 2
    Kp = 4;          % Pilot gain
    Te = 0.5;          % Time delay (tau_e)
    TI = 1;          % Time constant

    Hmcr = Kp * exp(-Te * s_fd) / (TI*s_fd + 1);  % Transfer function

    [mag,phase] = bode(Hmcr,W);                   % bodeplot of the pilot transfer function
    Mmcr = squeeze(mag(1,1,:));                   % magnitude of the pilot bodeplot
    Amcr = squeeze(phase(1,1,:));                 % phase angle of the pilot bodeplot

elseif t == 3
    Kp = 4;          % Pilot gain
    Te = 0.5;          % Time delay (tau_e)
    TL = 1;          % lag constant

    Hmcr = Kp * (TL*s_fd + 1) * exp(-Te * s_fd);  % Transfer function    

end

% Showing bode plot of selected subject and trial
figure(3); clf; 
ax31 = subplot(211); hold on
    pl31(1) = loglog(W, W*0+10^0,':k');             % straight line at 10^0 magnitude
    pl31(2) = errorbar(W,Mpmn,Mpstd,'Color','b');   % mean pilot magnitude with errorbars of the st. dev. of the magn.
    pl31(3) = errorbar(W,Molmn,Molstd,'Color','r'); % mean open-loop magnitude with errorbars of the st. dev. of the magn.
    pl31(4) = loglog(W,Msys,'.-g');                 % system magnitude for each frequency point
    pl31(5) = loglog(W,Mmcr,'--c');                 % uncomment for McRuer magnitude estimate
    % loglog(W2, W2*0+1,'.-m')                      % displays the interpolated vector when placed in the correct figure
    % loglog(W2(iw2), W2(iw2)*0+1,'ok')             % displays the found vector value with a circle when placed in the correct figure

    set(gca,'XScale','log','YScale','log')          % all shown in double logarithmic axis 
    legend(pl31(2:5),'Pilot','Open Loop','System', 'pilot model')  % legend three selected plots, extend this when uncommenting the McRuer plot
    ylabel('Magnitude')
    title(['Bode Plot of subject: ',num2str(s),' , task: ',num2str(t)])
    
ax32 = subplot(212); hold on
    errorbar(W,Apmn,Apstd,'Color','b')              % mean pilot phase angle with errorbars of the st. dev. of the angle
    errorbar(W,Aolmn,Aolstd,'Color','r')            % mean open-loop phase angle with errorbars of the st. dev. of the angle
    plot(W, Asys,'.-g')                             % system phase angle
    % plot(W, Amcr,'.--c')                          % uncomment for McRuer phase angle estimate
    
    set(gca,'XScale','log','YScale','linear')       % all shown in single logarithmic axis
    ylabel('Phase Angle [deg]')
    xlabel('frequency (Hz)')
    
linkaxes([ax31,ax32],'x')                           % linking x-axis of both subplots


%% calculate and show data of metric plots of RMS error and STD input
% show measured metrics for each subject and subtask with mean and standard
% deviation as errorbars for all repetitions

% calculate and plot RMS error and STD input for each subject and task
figure(4);clf
for s=1:nsub
    for t=1:ntsk
        
        % calculate mean and standard deviation of RMS error and STD input
        rms_err = data.sub(s).task(t).rms_err(:);
        std_inp = data.sub(s).task(t).std_inp(:);
        m_rms_err = mean( rms_err );     % mean RMS error of all repetitions
        s_rms_err = std(  rms_err );     % st. dev. RMS error of all repetitions
        m_std_inp = mean( std_inp );     % mean STD input of all repetitions
        s_std_inp = std(  std_inp );     % st. dev of STD inputs of all repetitions
        
        % plot mean and st.dev oi both RMS error and STD input with for
        % each subject a horizontal shift(s/10) over the x-axis per task
        subplot(211); hold on
            plot((t+s/10).*ones(length(rms_err(:)),1), rms_err(:),'.')     % individual data points
            plot(t+s/10,m_rms_err,'o')                                     % mean RMS error
            errorbar(t+s/10,m_rms_err,s_rms_err,'k')                       % errorbar of STD input
            title('RMS Error')                              
            set(gca,'XTick',[1 2 3]+2.5/10,'XTickLabel',{'A','B','C'})
        subplot(212); hold on
            plot((t+s/10).*ones(length(std_inp(:)),1), std_inp(:),'.')     % individual data points
            errorbar(t+s/10,m_std_inp,s_std_inp,'k')                       % errorbar of STD input
            plot(t+s/10,m_std_inp,'ob')                                    % mean STD input  
            title('STD Input')
            set(gca,'XTick',[1 2 3]+2.5/10,'XTickLabel',{'A','B','C'})
            
    end
end
