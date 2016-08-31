function Plot_Phoneme_ERPs(Phn, is_HG, appelation, yrange, varargin)
%% Plots depending on the number of inputs
ecog = Phn.ecog;
time_axis = Phn.time_axis;
fs_ecog = 250;

if isempty(yrange)
    scale_y = false;
else
    scale_y = true;
end
x_plot_range = [-100,250];
clip_data = true; %command to clip data to include only a specific time window
if clip_data
    twin(1) = find(time_axis <= x_plot_range(1),1,'last');
    twin(2) = find(time_axis >= x_plot_range(2),1);
    time_axis = time_axis(twin(1):twin(2));
    ecog = ecog(:,twin(1):twin(2),:);
end
%time_window = 1:(2*fs_ecog);

%time_axis = 1000*((1:375)/fs-0.2);
%time_axis = 1000*((1:size(ecog,2))/fs_ecog-(pre_event_duration));
num_chans = size(ecog,1);
num_tags = length(varargin);

%% Establish a nxinputs matrix of data_tags
if num_tags == 0
    Data_tag = true(length(Phn.phn_names),1);
else
    Data_tag = false(length(Phn.phn_names),num_tags);
    for i = 1:num_tags
        Data_tag(:,i) = varargin{i};
    end
end



%% Plot
%colorlist = {'k', 'r','b','c','y'};
 colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple

figure;
for i = 1:num_chans
    subplot(1,num_chans,i);
%     p = plotGridPosition(i, 8, 4); 
%     subplot('position',p);
    hold on;
    for j = 1:size(Data_tag,2)
        %shadedErrorBar(time_axis, mean(ecog(i,:,Data_tag(:,j)),3),nansem(ecog(i,:,Data_tag(:,j)),3),colorlist{j},0.8);
        shadedErrorBar(time_axis, mean(ecog(i,:,Data_tag(:,j)),3),nansem(ecog(i,:,Data_tag(:,j)),3),colorlist(j,:),0.8);
    end    
    xlim(x_plot_range);
    %ylim(get(gca,'ylim'))
    if scale_y
        ylim(yrange)
    end
    plot([0 0],get(gca,'ylim'),'k','LineWidth',1.5)
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Time (s)')
    %ylabel('Awesome Brain Stuff')
    if is_HG
        ylabel('High Gamma Intensity')
    else
        ylabel('ECoG Intensity')
    end
    %title(['Awesome Brain Stuff in Ch ',num2str(i),' For Phoneme ERPs'])
    %title(['High Gamma in Ch ',num2str(i),' For Phoneme ERPs']) %Fricative (blue) and Plosive (red) ERPs'])
    if is_HG
        title([appelation ' High Gamma at Phoneme Onset in Ch ',num2str(i)]) %Fricative (blue) and Plosive (red) ERPs'])
    else
        title([appelation ' Broadband at Phoneme Onset in Ch ',num2str(i)]) %Fricative (blue) and Plosive (red) ERPs'])
    end
        %axis('tight')
end
a=1;