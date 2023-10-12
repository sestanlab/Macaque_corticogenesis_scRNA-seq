
!Plot the trend on 3D surface plot
addpath('/Users/shaojie/Research/matlab/custom_scripts')
cd /Users/shaojie/Research/matlab/IPC-ExN-3D_UMAP/
for jj = ["nonIT"]
    nrows = 4; ntimes = 100;
    file_name = strcat("ExN_sliding-window_", jj, "_nDEGs_log.txt")
    win_avg_data = readtable(file_name,'ReadRowNames',true,'HeaderLines',1)
    win_avg_data = table2array(win_avg_data)
    
    % switch the order of regions (rows): tc/oc
    new_order = [1:2, 4, 3];
    win_avg_data = win_avg_data(new_order, :)
    win_avg_data = win_avg_data(fliplr(1:nrows),:)
    win_avg_data = standardizeMissing(win_avg_data,-1)



    % Build the x-y axis data
    % xx = linspace(0,1,ntimes)
    xx = linspace(0.17, 25.9, 100)
    a = repmat(xx,nrows,1)

    
    % Do the ribbon plots
    s = ribboncoloredZ(a',win_avg_data');
    
    % Set the x y ticks & labels
    % yticks([0:0.2:1]);
    yticks([0:5:27]);
    xticks([1:nrows]);
    xticklabels(fliplr({"FC","MSC","OcC","TC"}));
    xlim([0, (nrows + 1)]); 
    ylim([0, 26])
    zlabel("Region Differences"); ylabel("Progression score"); xlabel("Region");

    % Set a red line at the each  end of the ribbon
    %for i = 1:nrows
    %    line([i,i],[0,0],[0,win_avg_data(i,1)],'Color','red','LineStyle','--')
    %end
    %for i = 1:nrows
    %    line([i,i],[1,1],[0,win_avg_data(i,length(xx))],'Color','red','LineStyle','--')
    %end
    set(gcf, 'Renderer', 'painters')
    colorbar
    clim([4 8.5])
    
    view(45.6574,28.3860);
end




