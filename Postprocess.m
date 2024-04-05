
%ICL, meshing limit line and curvature interference limit line
%ICL of second envelope
figure;
fid=fopen('resuICLOSE.txt');
colo=rand(1,3);
numbTemp=999;
dataXY=zeros(numbTemp,2);
eid=0;
lid=0;
while ~feof(fid)
    stri=fgets(fid);
    if stri(1) == '#' || stri(1) == '@'
        lid=lid+1;
        plot(dataXY(1:eid,1),dataXY(1:eid,2),'-','Color',colo);hold on;
        eid=0;
        if stri(1) == '@'
            colo=rand(1,3);
        end
    else
        eid=eid+1;
        coor=sscanf(stri,'%e %e %e');
        dataXY(eid,1)=-coor(3);
        dataXY(eid,2)=sqrt(coor(1)^2+coor(2)^2);
    end
end
fclose(fid);
%meshing limit line of second envelope
data=load('resuMLLOSE.txt');
numbTemp=size(data,1);
if numbTemp > 0
    dataXY=zeros(numbTemp,2);
    dataXY(:,1)=-data(1:numbTemp,3);
    dataXY(:,2)=sqrt(data(1:numbTemp,1).^2+data(1:numbTemp,2).^2);
    plot(dataXY(:,1),dataXY(:,2),'b-','LineWidth',1.5);hold on;
end
%curvature interference limit line of second envelope
data=load('resuCILLOSE.txt');
numbTemp=size(data,1);
if numbTemp > 0
    dataXY=zeros(numbTemp,2);
    dataXY(:,1)=-data(1:numbTemp,3);
    dataXY(:,2)=sqrt(data(1:numbTemp,1).^2+data(1:numbTemp,2).^2);
    plot(dataXY(1:numbTemp,1),dataXY(1:numbTemp,2),'ko','MarkerSize',1.5);hold on;
end
xlim([-0.5,0.5]);
ylim([0.0,0.6]);

%Tooth surface grid discretization
a=0.25;%working center distance
%worm
figure;
data = load('resuWOTSGRID.txt');
for ti=1:1:size(data,1)
    data(ti,1:3)=([1.0,0.0,0.0;0.0,0.0,1.0;0.0,-1.0,0.0] * (data(ti,1:3)).' +[-a;0.0;0.0]).';
end
scatter3(data(:,1),data(:,2),data(:,3),5,[0.1,0.1,0.5],'filled');
hold on;
%worm wheel
data = load('resuWHTSGRID.txt');
data_f = load('resuWHTSPHAS.txt');
data_p=zeros(size(data,1),3);
for tj=1:1:5
    num_tj=0;
    for ti=1:1:size(data,1)
        if data_f(ti)==tj
            num_tj = num_tj+1;
            data_p(num_tj,1:3) = data(ti,1:3);
        end
    end
    if tj == 1
        scatter3(data_p(1:num_tj,1),data_p(1:num_tj,2),data_p(1:num_tj,3), 15,'c','filled');
        hold on;
    elseif tj == 2
        scatter3(data_p(1:num_tj,1),data_p(1:num_tj,2),data_p(1:num_tj,3), 15,'r','filled');
        hold on;
    elseif tj == 3
        scatter3(data_p(1:num_tj,1),data_p(1:num_tj,2),data_p(1:num_tj,3), 15,'b','filled');
        hold on;
    elseif tj == 4
        scatter3(data_p(1:num_tj,1),data_p(1:num_tj,2),data_p(1:num_tj,3), 15,'g','filled');
        hold on;
    elseif tj == 5
        scatter3(data_p(1:num_tj,1),data_p(1:num_tj,2),data_p(1:num_tj,3), 15,'m','filled');
        hold on;
    end
end
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

%Displacement
node = load('resuNode_0.txt');
elem = load('resuElem_0.txt')+1;
disp = load('resuDisp_0.txt');
figure;
patch('Vertices', node(:,:), 'Faces', elem(:,1:3), 'FaceVertexCData', disp(:,3), 'FaceColor', 'interp');
hold on;
axis equal;
colorbar;
colormap(jet);
view(30,40);
xlabel('x');
ylabel('y');
zlabel('z');

%Contact pressure
nodeSlav = load('resuNode_0.txt');
contForc=load('resuCont.txt');
contForc=contForc(contForc(:,2)>0.0,:);
figure;
scatter3(nodeSlav(1+contForc(:,1),1),nodeSlav(1+contForc(:,1),2),nodeSlav(1+contForc(:,1),3),15,contForc(:,2),'filled');
colorbar;
colormap(jet);
xlabel('x');
ylabel('y');
zlabel('z');

%Initial contact gap (very time consuming)
contElem=load('resuContElem.txt');
figure;
scatter3(contElem(:,1),contElem(:,2),contElem(:,3),15,contElem(:,7),'filled');hold on;
colorbar;
