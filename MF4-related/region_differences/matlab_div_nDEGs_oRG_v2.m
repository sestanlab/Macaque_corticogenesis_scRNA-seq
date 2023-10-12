
!Plot the trend on 3D surface plot
addpath('/Users/shaojie/Research/matlab/custom_scripts')
cd /Users/shaojie/Research/matlab/Region_div/
for jj = [ "oRG"]
    nrows = 4; ntimes = 100;
    file_name = strcat("RGC_diffs_nDEGs_", jj, "_log_v2")
    win_avg_data = readtable(file_name,'ReadRowNames',true,'HeaderLines',1)
    win_avg_data = table2array(win_avg_data)
    
    % switch the order of regions (rows): tc/oc
    new_order = [1:2, 4, 3];
    win_avg_data = win_avg_data(new_order, :)
    win_avg_data = win_avg_data(fliplr(1:nrows),:)
    win_avg_data = standardizeMissing(win_avg_data,-1)



    % Build the x-y axis data
    % xx = linspace(0,1,ntimes)
    xx = linspace(9.51, 18.65, 100)
    a = repmat(xx,nrows,1)

    
    % Do the ribbon plots
    s = ribboncoloredZ(a',win_avg_data');
    
    % Set the x y ticks & labels
    % yticks([0:0.2:1]);
    yticks([9:2:19]);
    xticks([1:nrows]);
    xticklabels(fliplr({"FC","MSC","OcC","TC"}));
    xlim([0, (nrows + 1)]); 
    ylim([9.5, 18.7])
    zlabel("Region Differences"); ylabel("Progression score"); xlabel("Region");

    % Set a red line at the each  end of the ribbon
    %for i = 1:nrows
    %    line([i,i],[0,0],[0,win_avg_data(i,1)],'Color','red','LineStyle','--')
    %end
    %for i = 1:nrows
    %    line([i,i],[1,1],[0,win_avg_data(i,length(xx))],'Color','red','LineStyle','--')
    %end
    set(gcf, 'Renderer', 'painters')
    clim([2.647, 8.558])
    colorbar
    
    view(131.4175,20.6);
end




