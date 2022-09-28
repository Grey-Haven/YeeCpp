clear;

Nxs = [8,16,32,64,128.256];

vidPath = pwd + "/";

for g = Nxs
    
    gxg = num2str(g) + "x" + num2str(g);
    
    vidObj = VideoWriter(vidPath + gxg + "_vid", 'MPEG-4');
    open(vidObj);

    path = pwd + "/results/" + num2str(g) + "x" + num2str(g);
    ExFiles = dir(path + "/Ex_*");
    EyFiles = dir(path + "/Ey_*");
    HzFiles = dir(path + "/Hz_*");
    dx = 50/(g-1);
    dy = 50/(g-1);
    x = -25:dx:25;
    y = -25:dy:25;
    for i = 1:length(HzFiles)

        subplot(2,2,1);
        HzFile = HzFiles(i);
        HzFileName = strcat(HzFile.folder, "/", HzFile.name);
        disp(HzFileName);
        table = readmatrix(HzFileName);
        subplot(2,2,1);
        surf(x,y,table);
        xlim([-25,25]);
        ylim([-25,25]);
        title("Hz");
        xlabel("x");
        ylabel("y");

        subplot(2,2,2);
        ExFile = ExFiles(i);
        ExFileName = strcat(ExFile.folder, "/", ExFile.name);
        table = readmatrix(ExFileName);
        surf(x,y,table);
        xlim([-25,25]);
        ylim([-25,25]);
        title("Ex");
        xlabel("x");
        ylabel("y");
        
        subplot(2,2,3);
        EyFile = EyFiles(i);
        EyFileName = strcat(EyFile.folder, "/", EyFile.name);
        table = readmatrix(EyFileName);
        surf(x,y,table);
        xlim([-25,25]);
        ylim([-25,25]);
        title("Ey");
        xlabel("x");
        ylabel("y");

        drawnow;
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end
    close(vidObj);
end