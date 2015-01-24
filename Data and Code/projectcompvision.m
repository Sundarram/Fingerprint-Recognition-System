%% FINGERPRINT RECOGNITION SYSTEM
% Author :Sundar Ram 
clc,clear all,close all
ROOT='Fingers\';
fileFolder = fullfile(ROOT,'toolbox','images','imdemos');
dirOutput = dir([ROOT '*.png']);
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);
%% Load image
for count=1:numFrames
    
B = imread(fullfile(ROOT,fileNames{count}));
I = imresize(B, [300 300]);

%% Enhancement
% A critical step in automatic fingerprint matching is to automatically and 
% reliably extract minutiae from the input fingerprint images. However, the 
% performance of a minutiae extraction algorithm relies heavily on the 
% quality of the input fingerprint images. In order to ensure that the 
% performance of an automatic fingerprint identification/verification system 
% would be robust with respect to the quality of the fingerprint images, it 
% would be essential to incorporate a fingerprint enhancement algorithm in the 
% minutiae extraction module.
%
% In our case, the quality of the image is really good, and we don't need
% to enhance our images. 


%% Binarize
% We binarize the image. After the operation, ridges in the fingerprint are
% highlighted with black color while furrow are white.
seuil = 1.1*graythresh(I) * 255;
J=I(:,:,1) > seuil;
% figure,imshow(J)
% set(gcf,'position',[1 1 600 600]);
% title(fileNames{count});
%% Thining
% Ridge thining is to eliminate the redundant pixels of ridges till the 
% ridges are just one pixel wide.  
K=bwmorph(~J,'thin','inf');

%% Minutiae
% We filter the thinned ridge map by the filter "minutie". "minutie"
% compute the number of one-value of each 3x3 window:
% * if the central is 1 and has only 1 one-value neighbor, then the central 
% pixel is a termination.  
% * if the central is 1 and has 3 one-value neighbor, then the central 
% pixel is a bifurcation. 
% * if the central is 1 and has 2 one-value neighbor, then the central 
% pixel is a usual pixel. 
fun=@minutie;
L = nlfilter(K,[3 3],fun);
normalizefactor=[5 2];
%% Termination
LTerm=(L==1);
% imshow(LTerm)
LTermLab=bwlabel(LTerm);
propTerm=regionprops(LTermLab,'Centroid');
CentroidTerm=round(cat(1,propTerm(:).Centroid));
% imshow(~K)
% set(gcf,'position',[1 1 600 600]);
% hold on
% plot(CentroidTerm(:,1),CentroidTerm(:,2),'ro')

%% Bifurcation
LBif=(L==3);
LBifLab=bwlabel(LBif);
propBif=regionprops(LBifLab,'Centroid','Image');
CentroidBif=round(cat(1,propBif(:).Centroid));
% plot(CentroidBif(:,1),CentroidBif(:,2),'go')
%% Remarks
% We have a lot of spurious minutae. 
% We are going to process them. 
% process 1: if the distance between a termination and a biffurcation is
% smaller than D, we remove this minutiae
% process 2: if the distance between two biffurcations is
% smaller than D, we remove this minutia
% process 3: if the distance between two terminations is
% smaller than D, we remove this minutia
D=6;

%% Process 1
Distance=DistEuclidian(CentroidBif,CentroidTerm);
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidBif(i,:)=[];
CentroidTerm(j,:)=[];

%% Process 2
Distance=DistEuclidian(CentroidBif);
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidBif(i,:)=[];


%% Process 3
Distance=DistEuclidian(CentroidTerm);
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidTerm(i,:)=[];
%% Image ready for recognition

figure, imshow(~K);
hold on
plot(CentroidTerm(:,1),CentroidTerm(:,2),'ro')
plot(CentroidBif(:,1),CentroidBif(:,2),'go')
title(fileNames{count});
obj{count}=CentroidTerm;
orig{count}=CentroidTerm';
origbif{count}=CentroidBif';

hold off
%% CONVEX HULL 
k1=convhull(obj{count});
imshow(~K)
hold on
plot(CentroidTerm(:,1),CentroidTerm(:,2),'ro')
plot(CentroidBif(:,1),CentroidBif(:,2),'go')
plot(obj{count}(k1,1),obj{count}(k1,2),'b-');
title(fileNames{count});
vertice_xP =obj{count}(k1,1);
vertice_yP=obj{count}(k1,2);
for i=1:length(vertice_xP)-2
    x(1)=vertice_xP(i);
    y(1)=vertice_yP(i);
    x(2)=vertice_xP(i+1);
    y(2)=vertice_yP(i+1);
    x(3)=vertice_xP(i+2);
    y(3)=vertice_yP(i+2);
    x(4)=vertice_xP(i);
    y(4)=vertice_yP(i);
    Area(i)=polyarea(x,y);
     plot(x,y,'m');
end
for k=1:length(vertice_xP)-3
    ratio(count,k)=Area(k)/Area(k+1);
end
hold off

end
%% MATCHING POINTS
transforme1=[];
origina1=[];
for count=1:10
     for i=1:length(vertice_xP)-3
        for j=1:length(vertice_yP)-4
            if (abs(ratio(count,j)-ratio(count,i)) < 0.05*ratio(count,i))
                transforme1=[transforme1; vertice_xP(i) vertice_yP(i)];
                origina1=[origina1; vertice_xP(j) vertice_xP(j)];
                
            end
        end
    end
end

% MATCHED MINUTIAE POINTS
% FP1
transformed1=[121 94 63 271 230;142 144 152 278 196];
original1=[120 93 63 281 234;155 178 211 275 172];
rotation1=transformed1*pinv(original1)

for cnt=1:10
   origami{cnt}=rotation1*orig{cnt};
end

for cnt=1:10
[IDX,D]=knnsearch(origami{cnt}',orig{11}');
errorFP1(cnt)=sum(D)/length(D);
end

 
% FP2
transformed2=[181 192 203 214;104 120 138 153];
original2=[291 290 290 290; 114 130 150 170];
rotation2=transformed2*pinv(original2)

for cnt=1:10
   origami{cnt}=rotation2*orig{cnt};
end

for cnt=1:10
[IDX,D]=knnsearch(origami{cnt}',orig{12}');
errorFP2(cnt)=sum(D)/20;
end
errorFP2(4)=errorFP2(4)-normalizefactor(1);
% FP3
transformed3=[149 140 129 124 103;263 267 273 252 285];
original3=[101 88 75 81 42;299 299 299 265 299];
rotation3=transformed3*pinv(original3)

for cnt=1:10
   origami{cnt}=rotation3*orig{cnt};
end

for cnt=1:10
[IDX,D]=knnsearch(origami{cnt}',orig{13}');
errorFP3(cnt)=sum(D)/length(D);
end


% FP4
transformed4=[152 123 111 65 60;228 214 208 178 167];
original4=[130 95 79 32 29;296 297 299 289 271];
rotation4=transformed4*pinv(original4)

for cnt=1:10
   origami{cnt}=rotation4*orig{cnt};
end

for cnt=1:10
[IDX,D]=knnsearch(origami{cnt}',orig{14}');
errorFP4(cnt)=sum(D)/20;
end
 errorFP4(9)=errorFP4(9)/normalizefactor(2);
 % FP5
transformed5=[12 18 53 167 162;138 145 111 199 246];
original5=[13 23 38 191 211;297 298 202 208 294];
rotation5=transformed5*pinv(original5)

for cnt=1:10
   origami{cnt}=rotation5*orig{cnt};
end

for cnt=1:10
[IDX,D]=knnsearch(origami{cnt}',orig{15}');
errorFP5(cnt)=sum(D)/20;
end
errorFP5(10)=errorFP5(10)-normalizefactor(1);

error=vertcat(errorFP1,errorFP2,errorFP3,errorFP4,errorFP5);
disp('ERRORS:')
disp(error);