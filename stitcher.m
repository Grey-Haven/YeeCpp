clear;

Nxs = [256];

for g = Nxs
    path = pwd + "\results\" + num2str(g) + "x" + num2str(g);
    ExFiles = dir(path + "\Ex_*");
    EyFiles = dir(path + "\Ey_*");
    HzFiles = dir(path + "\Hz_*");
    dx = 50/(g-1);
    dy = 50/(g-1);
    x = -25:dx:25;
    y = -25:dy:25;
    for i = 1:length(HzFiles)
        file = HzFiles(i);
        fileName = strcat(file.folder, "\", file.name);
        disp(fileName);
        table = readmatrix(fileName);
        surf(x,y,table);
        xlim([-25,25]);
        ylim([-25,25]);
        zlim([-.1,.1])
        xlabel("x");
        ylabel("y");
        drawnow;
    end
end