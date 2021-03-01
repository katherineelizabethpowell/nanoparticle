clear;
clc;

%Folder with the images.
folder='200nmbeads';

%Must name the folder variable after the folder that contains the oib files
files = dir([folder, '/*.oib']);

%Finds all filenames in the folder
filenames = {files.name};
Lfile = length(filenames);
newFolder = [folder, ' plots ', num2str(1), 'x', num2str(1)];
% mkdir(newFolder);
maxDye = ones(1,length(filenames));

%Loops through each oib file in the folder
for x=1:length(filenames)
    F = fullfile(folder, filenames(x));
    img = bfopen(F{1});

    
    %To just make plots
    %x = number of files being run
    %: = self defined
    Output(:,x) = greenBeadPlotter(img, folder, newFolder, filenames);
    
end
 
function Output = greenBeadPlotter(data, folder, newFolder, filenames)
[b,a] = butter(1,20/50);

%Change the number based on the step size for the images
step = 2

%extracts the cell array that contains all the zstacks
images = length(data{1,1}); 

%assigns the number of zstacks to rows
[rows, cols] = size(data{1,1}{1,1}); 

rBeadProf = []; %initializing 200nm Red Bead profiles
gBeadProf = []; %initializing 20nm or 500nm Green Bead profiles
DyeProf = [];
Dyenorm = [];

%Checks the number of channels and runs with 2 channels if no Dextran. Or
%runs with 3 channels if Dextran is present. 
num_channels = data{1,4}.getChannelCount(0);


    for slice = 1:(images/num_channels)

    %Channel 1 = {slice, 1}
    %Channel 2 = 
    
if num_channels == 1 
    Dye = mat2cell(data{1,1}{slice*num_channels,1}, rows, cols);
    
elseif num_channels == 2
    %Dye is in Channel 1? slice*num_channels-1
    %Dye is in Channel 2? slice*num_channels
    rBead = mat2cell(data{1,1}{slice*num_channels,1}, rows, cols);
    Dye = mat2cell(data{1,1}{slice*num_channels-1}, rows, cols);
    
    
elseif num_channels == 3
    %Dye is in Channel 1? slice*num_channels-2
    %Dye is in Channel 2? slice*num_channels-1
    %Dye is in Channel 3? slice*num_channels
    Dye = mat2cell(data{1,1}{slice*num_channels-1,1}, rows, cols);
end    

        for y=1:1
            %replaces them with the average of the intensity profile
            rBead{y} = mean(mean(rBead{y}));
            Dye{y}= mean(mean(Dye{y}));
        end   
    rBead = cell2mat(rBead);
    %assigns chunks(x) to out(x)
    out3 = rBead(y);
    rBeadProf = [out3,rBeadProf]; 
    
    Dye = cell2mat(Dye);
    out2 = Dye(y);
    DyeProf = [out2, DyeProf];
    end
    
    %Fix profiles so surface is to the left.
    DyeProf = fliplr(DyeProf);
    rBeadProf = fliplr(rBeadProf);

    Xaxis = [];
    for i = 1:(images/num_channels)
        Xaxis = [i Xaxis];
    end
    Xaxis = Xaxis/step; %Because there is 0.1um per slice
   % Xaxis = Xaxis*step; %For step as fraction of a um
%% FIGURE 1
    %To plot RAW Data:
    figure(1);
    hold on;
    plot(Xaxis,rBeadProf);
    %plot(Xaxis, DyeProf);
    hold off;
    
    xlabel('Zposition in stack [um]'); 
    ylabel('Raw Intensity');
    title('Raw Intensity Profiles')
    xlim([0 (images/num_channels)/step])
    %xlim([0 (images/num_channels)/2])
    %ylim([0 max(rBeadProf(:))])
    
 %   axis[0 images 0 max(rBeadProf(:))]; % [min_x max_x min_y max_y]
    %legend(); %General legend which says data1, data2, etc.
  %  legend(filenames());
    legend(filenames(:)); 
    box on
    set(gca,'LineWidth',2);
    
    Output = [Xaxis,999999999,rBeadProf];
    
    %% FIGURE 2
    %To plot NORMALIZED Data:
    rBeadnorm = (rBeadProf-rBeadProf(end))/(max(rBeadProf(:))-rBeadProf(end));
    
    figure(2);
    hold on;
    plot(Xaxis,rBeadnorm);
    hold off;

    xlabel('Zposition in stack [um]'); 
    ylabel('Normalized Intensity');
    %legend();
    legend(filenames(:));
    ylim([0 1.1])
    xlim([0 (images/num_channels)/step])
    %xlim([0 (images/num_channels)/2])
    title('Normalized Intensity Profiles')
    box on
    set(gca,'LineWidth',2);
    
    Output = [Xaxis,999999999,rBeadnorm];

%% FIGURE 3

 
    dyemax = max(DyeProf);
    [dmp,dpeak] = findpeaks(DyeProf, 'MinPeakHeight', dyemax-1)
    f = fit(Xaxis(dpeak-15:dpeak+15)',DyeProf(dpeak-15:dpeak+15)','gauss2')
     coefficients = coeffvalues(f);    
   
     num_fitpts = 1000; %Number of data points to make gaussian line out of
     sz = linspace(1,(images/num_channels)/step,num_fitpts); 
     f_line = zeros(1,num_fitpts);
     f_xaxis = [1:images/num_channels];  
         for m = 1:length(sz)
             f_line(m) = coefficients(1)*exp(-((sz(m)-coefficients(2))/coefficients(3))^2) + coefficients(4)*exp(-((sz(m)-coefficients(5))/coefficients(6))^2);
         end

    figure(3);
    hold on;
    plot (sz,f_line);
    plot (Xaxis, DyeProf);
    hold off;
    xlabel('Z slice'); 
    ylabel('Intensity');
    legend();
%    ylim([0 1.1])
%    xlim([0 max(reg_x)/2])
    title('Raw Data with Gaussian Fit')
    box on
     set(gca,'LineWidth',2);

     f_line = fliplr(f_line);
     DyeProf = fliplr(DyeProf);
     rBeadProf = fliplr(rBeadProf);
     Xaxis = fliplr(Xaxis);     
  
    %Find the peak in the fitted gaussian.
    %Returns the widths of the peaks as the vector w and 
    %the prominences of the peaks as the vector p.
     
    [pks,locs,w,p] = findpeaks(DyeProf);
    [Maxp,I_p] = max(p); %the peak with the max prominence  
    %Xaxis position of peak in the fitted curve
    DyeProfxlocus = locs(I_p);
    
    %Find the location of the peak in the DyeProfile curve and prepare to
    %plot only values on x-axis greater than this value for both curves.
    %sz(locs) is asking for the x-position of the peak like A(2).
    xthreshold = sz(DyeProfxlocus);
    moreThanThreshold = sz > xthreshold;
    f_line_truncated = f_line(locs:end);
    %x_axis_threshold goes from x@peak to end of zstack profile
    x_axis_threshold = sz(moreThanThreshold);
  
    
    %Rounding to nearest integer
    Threshold_inDyeProf = round(xthreshold);
    Offset = Threshold_inDyeProf - xthreshold;

    %Shift the FIT to align peak of gaussian at x = 0.
    x_axis_threshold_FITshift = xthreshold - xthreshold(1);
    
    %Shift the Dye and fit to realign around x = 0.
    x_axis_threshold_DYEshift = x_axis_threshold - Threshold_inDyeProf;
    
    DyeProf_truncated = DyeProf(DyeProfxlocus:end);
    rBeadProf_truncated = rBeadProf(DyeProfxlocus:end);
    Dye_x_axis = [0:length(DyeProf_truncated)-1];
    Dye_x_axis = Dye_x_axis/step
    
    figure(4);
    hold on;
    %plot(Dye_x_axis, DyeProf_truncated);
    plot (Dye_x_axis,rBeadProf_truncated);
    %plot(x_axis_threshold_FITshift,f_line_truncated);
    %legend('DyeProf','Fit');
    hold off;
    %plot(Dye_x_axis,rBeadProf_truncated);
    
    xlabel('Zposition above glass [um]'); 
    ylabel('Raw Intensity');
    legend(filenames(:));
    %legend();
%     ylim([0 1.1])
%     xlim([0 images/2])
    title('Truncated at Glass Profile')
    box on
    set(gca,'LineWidth',2);
    
% %% FIGURE 4
% 
%     %Normalize data from Figure 3.
%     %DyeProf_truncatednorm = (DyeProf_truncated-DyeProf_truncated(end))/(max(DyeProf_truncated(:))-DyeProf_truncated(end));
%     rBeadProf_truncatednorm = (rBeadProf_truncated-rBeadProf_truncated(end))/(max(rBeadProf_truncated(:))-rBeadProf_truncated(end));
%     
%     figure(5);
%     hold on;
%     plot(Dye_x_axis, rBeadProf_truncatednorm);
%     %plot(x_axis_threshold_FITshift,f_line_truncated);
%     %legend('DyeProf','Fit');
%     hold off;
% 
%     xlabel('Zposition above glass [um]'); 
%     ylabel('Normalized Intensity');
%     %legend();
%     legend(filenames(:));
%     ylim([0 1.1])
% %     xlim([0 images/2])
%     title('Normalized and Truncated at Glass Profile')
%     box on
%     set(gca,'LineWidth',2);

    
end