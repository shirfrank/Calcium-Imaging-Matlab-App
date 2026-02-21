function [n_spike,f_spike,d_spike,p_spike,max_sig,p_spike_norm,amp_spike,mean_amp_spike,areaAmp_base,areaAmp_prom,onset_time,start_idx,offset_time,end_idx,rise_time30,decay_time30,timeLat,areaLat_30Prom,areaLat_0Base,duration,area_prom,area_base,inter_peak_intervals] = analys_oscillations_140124(sig,time_vec,glut_inj_idx,sanity_check,n_cells)

    %% Drift correction - optional
    %sig_corr = detrend(sig);
    %diff_drift=sig(1)-sig_corr(1);
    %sig=sig_corr+diff_drift;

    %% findpeaks()parameters - amplitude,index,half prominance width,prominance
    %Oligo peaks parameters
    max_sig=max(sig); %max signal amplitude
    normsig = (sig-min(sig))/(max(sig)-min(sig));
    O_duration_min = 21; %sec
    O_duration_max = 180; %sec
    O_prominence_min = 0.3*max_sig; 
    O_dis_min = 60; %sec
    O_amp_min = max(max_sig)*0.5;
    %,'WidthReference','halfheight'
    %find peaks parameters (unnormalized signal)
    [amp_spike,loc_spike,~,p_spike] = findpeaks(sig,time_vec,'MinPeakWidth',O_duration_min,'MaxPeakWidth',O_duration_max,'MinPeakProminence',O_prominence_min,'MinPeakHeight',O_amp_min,'MinPeakDistance',O_dis_min);
    %find peaks normalized prominence
    [~,~,~,p_spike_norm] = findpeaks(normsig,time_vec,'MinPeakWidth',O_duration_min,'MaxPeakWidth',O_duration_max,'MinPeakProminence',0.3,'MinPeakHeight',0.5,'MinPeakDistance',O_dis_min);
    n_spike=length(amp_spike);
    time_dur=time_vec(end)-time_vec(1);
    f_spike=(n_spike/time_dur)*60; %frequency - spike/min
    timeLat=loc_spike-(time_vec(glut_inj_idx)); %latency of spike relative to glutamate administration
    lat_idx=[];
    mean_amp_spike=[];
    for m=1:length(loc_spike)
        lat_idx(end+1)=find(time_vec==loc_spike(m));
        mean_amp_spike(end+1)=mean(sig(lat_idx(m)-1:lat_idx(m)+1));
    end       

    %% find duration btween peaks
    if n_spike ==1
       d_spike=NaN;
    else
       d_spike=diff(timeLat);
    end
    
    %% features
    %initialize features 
    rise_time30=[]; %rising time (30% to 90% of prominence)
    decay_time30=[]; %decay time (90% to 30% of prominence)
    area_prom=[]; %area under the curve (with prominence baseline)
    area_base=[]; %area under the curve (with baseline 0)
    areaLat_30Prom=[]; %time in which the area reaches 50% of total area (with prominence baseline)
    areaLat_0Base=[];
    areaAmp_base=[]; %area/area base (with baseline 0)
    areaAmp_prom=[]; %area/area base (with prominence baseline)
    duration=[]; %width of area base (bounded by 25% prominence threshold)
    onset_time=[]; %time to 30% of prominence before max amp
    start_idx=[];
    offset_time=[]; %time to 30% of prominence after max amp
    end_idx=[];
    inter_peak_intervals=[];
    baseline_prom=[];
    th30_val=[];
    th90_val=[];
    intx={};
    inty={};
    int90={};

    %% peaks intervals parameters
    for i=1:n_spike
        baseline_p=amp_spike(i)-p_spike(i); %set new baseline (prominence baseline)
        %th50=amp_spike(i)-p_spike(i)*(1-0.5); %50% prominence threshold
        %th50_win=find_thWin(sig,th50,time_vec,lat_idx(i)); %index of threshold crossing
        %th10=amp_spike(i)-p_spike(i)*(1-0.1); %10% prominence threshold
        %th10_win=find_thWin(sig,th10,time_vec,lat_idx(i)); %index of threshold crossing
        th30=amp_spike(i)-p_spike(i)*(1-0.30); %30% prominence threshold
        th30_win=find_thWin(sig,th30,time_vec,lat_idx(i)); %index of threshold crossing
        th90=amp_spike(i)-p_spike(i)*(1-0.9);%90% prominence threshold
        th90_win=find_thWin(sig,th90,time_vec,lat_idx(i)); %index of threshold crossing
        
        %time bound from 90% th
        bound=120;
        if time_vec(th90_win(1))-time_vec(th30_win(1))>bound
           th_start=th90_win(1)-bound;
        else
           th_start=th30_win(1);
        end

        if time_vec(th30_win(2))-time_vec(th90_win(2))>bound
           th_end=th90_win(2)+bound;
        else
            th_end=th30_win(2);
        end

        end_idx(end+1)=th_end;
        start_idx(end+1)=th_start;
        baseline_prom(end+1)=baseline_p;
        th30_val(end+1)=th30;
        th90_val(end+1)=th90;
        int90{end+1}=th90_win;

        %onset-offset intersection
        if i>1
            if th_start-end_idx(end-1)<0
               end_idx(end-1)=th_start;
            end
        end
    end

%% based on intervals parameters
    for i=1:n_spike
        th_end=end_idx(i);
        th_start=start_idx(i);
        baseline_p=baseline_prom(i);
        th90_win=int90{i};

        w_base=time_vec(th_end)-time_vec(th_start); %width of peak base
        rise_time=time_vec(th90_win(1))-time_vec(th_start); 
        decay_time=time_vec(th_end)-time_vec(th90_win(2));
        x_int=time_vec(th_start:th_end);
        y_int=sig(th_start:th_end);
        curr_area_base=trapz(x_int,y_int);
        if baseline_p>0
            curr_area_prom=trapz(x_int,y_int-baseline_p);   
        else
            curr_area_prom=trapz(x_int,y_int+baseline_p);
        end

        %areaLat_prom=find_areaLat(y_int,curr_area_prom,x_int,baseline_p);
        %areaLat_base=find_areaLat(y_int,curr_area_base,x_int,0);

        duration(end+1)=w_base;
        rise_time30(end+1)=rise_time;
        decay_time30(end+1)=decay_time;
        area_prom(end+1)=curr_area_prom;
        area_base(end+1)=curr_area_base;
        %areaLat_30Prom(end+1)=areaLat_prom-(time_vec(glut_inj_idx));
        %areaLat_0Base(end+1)=areaLat_base-(time_vec(glut_inj_idx));
        onset_time(end+1)=time_vec(th_start);
        offset_time(end+1)=time_vec(th_end);
        areaAmp_base(end+1)=curr_area_base/w_base;
        areaAmp_prom(end+1)=curr_area_prom/w_base;
        intx{end+1}=x_int;
        inty{end+1}=y_int;
    end

    if n_spike==1
        inter_peak_intervals=[];
    else
        for j=1:n_spike-1
            inter_peak_intervals(end+1)=onset_time(j+1)-offset_time(j);
        end
    end
    disp('Run completed successfully')

    %% sanity check
    if sanity_check==1
        for k=1:n_spike

                subplot(n_spike,1,k,'Color','w')
                hold on
                plot(time_vec,sig,'Color',[0,0,0],'LineWidth',2);
                hold on
                int_x=intx{k};
                int_y=inty{k};
                baseline_p=baseline_prom(k);

                if baseline_p>0
                   prominence_int=zeros(length(int_x),1)+baseline_p;
                    area([int_x;int_x],[prominence_int;int_y],0,'FaceColor',[0,0,0],'FaceAlpha',0.1,'EdgeColor',[1,0,0],'EdgeAlpha',0.1);
                    %area(int_x,int_y,0,'FaceColor',[0,0,0],'FaceAlpha',0.1,'EdgeColor',[1,0,0],'EdgeAlpha',0.1);
                else
                    prominence_int=zeros(length(int_x),1);
                    area([int_x;int_x],[prominence_int;int_y],baseline_p,'FaceColor',[0,0,0],'FaceAlpha',0.1,'EdgeColor',[1,0,0],'EdgeAlpha',0.1);
                    %area(int_x,int_y,baseline_prom(k),'FaceColor',[0,0,0],'FaceAlpha',0.1,'EdgeColor',[1,0,0],'EdgeAlpha',0.1);
                end
        
                %yl1=[-1 4];
                %ylim(yl1)
                scalebar('Border','LL','XUnit','sec','YUnit','\DeltaF/F_{0}');
                hold on
                yline(baseline_prom(k),'--','Prominence Baseline','LabelVerticalAlignment','bottom');
                yline(0,'-','Baseline','LabelVerticalAlignment','bottom');
                %yline(th50,'--','50% Prominence Threshold','LabelVerticalAlignment','bottom');
                yline(th30_val(k),'--','30% Prominence Threshold','LabelVerticalAlignment','bottom');
                %yline(th10,'--','10% Prominence Threshold','LabelVerticalAlignment','bottom');
                yline(th90_val(k),'--','90% Prominence Threshold','LabelVerticalAlignment','bottom');
                xline(onset_time(k),'--','Onset');
                xline(offset_time(k),'--','Offset');
                %xline(areaLat_0Base(k));
                %str = {'Peak Parameters:',['Area 25%Prom: ',num2str(curr_area_prom)],['Area 0Baseline: ',num2str(curr_area_base)],['Rise Time: ',num2str(rise_time),' [sec]'],['Decay time: ',num2str(decay_time),' [sec]']};
                %annotation('textbox','String',str,'BackgroundColor','white')\

        end
    end        

%% _____________ local functions _____________ %%
function[th_win] = find_thWin(sig,percentage_th10,time_vec,loc_peak)
% find intersection points of the signal with prominence percentage
% threshold.
    peakWidth = 0; %width of meaning window
    %std=0.05*percentage_th10;
    %first search for early crossing of baseline by going backwards from peak    
    start = NaN;
    for sample=loc_peak:-1:1+peakWidth
        currentAmp = mean(sig(sample-peakWidth:sample+peakWidth));
        if currentAmp <= percentage_th10
            start = sample;
            break;
        end
    end

    if sig(loc_peak)-sig(start)>120
        for sample=loc_peak:-1:1+peakWidth
            currentAmp = mean(sig(sample-peakWidth:sample+peakWidth));
            if currentAmp <= percentage_th10
                start = sample;
                break;
            end
        end
    end


    %now search for the later crossing of the threshold
    stop = NaN;
    for sample=loc_peak:length(time_vec)-peakWidth
        currentAmp = mean(sig(sample-peakWidth:sample+peakWidth)); 
        if currentAmp <= percentage_th10
            stop = sample;
            break;
        end
    end

    if sig(stop)-sig(loc_peak)>120
        for sample=loc_peak:-1:1+peakWidth
            currentAmp = mean(sig(sample-peakWidth:sample+peakWidth));
            if currentAmp <= percentage_th10
                stop = sample;
                break;
            end
        end
    end
    
    th_win=[start,stop];
end

%find percentage area latency
%function [areaLat]=find_areaLat(int_y,area,int_x,baseline)
%    
%    int_y=int_y-baseline;
%    precArea = area*0.5;
%        
%        sumArea = 0;
%        for j = 2:length(int_y)
%            currArea = trapz(int_x(j-1:j),int_y(j-1:j)); %area of 1 block
%            if sumArea < precArea % as long as we dornt reached the quantile area
%                sumArea = sumArea + currArea; %add the area of the current sample
%            else %if we reached the quantile area save the index
%                areaLat = int_x(j-1);
%                break
%            end
%        end 
% end
end

