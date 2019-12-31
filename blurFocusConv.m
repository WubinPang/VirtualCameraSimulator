%% Artificially generate blurred background and foreground image with focusing in between
function arti_im=blurFocusConv(images,depths,f,fn,fo_dis)
%% Load mat files

%load('D:/5132018/3DImage_TrainingSet/5case.mat');

%% Camera setup deciding on focal length f, aperture size F/# fn, pixel size delta and focusing plane distance p0
tic;
D=f/fn; % pupil size
%% Artificial defocusing
% focus_dis=zeros(1,num_fstack);
% i=100;
% orig_im=double(images(:,:,:,i)); 
          
dep_im=double(depths); %Load depth map
dep_value=unique(dep_im);
dep_num=size(dep_value,1);

orig_im=double(images);
[M,N,C]=size(orig_im);
FoV_h=10;
FoV_v=atand(M*tand(FoV_h)/N);
arti_im=zeros(M,N,C); %empty image for storing artifical focused imagedata

RMag=(fo_dis-f)/f; %reciprocal of magnification
y_im=2*fo_dis*tand(FoV_v/2)/RMag;
x_im=2*fo_dis*tand(FoV_h/2)/RMag; %image height in y and x dimension
delta=y_im/M;% pixel size    
DoF=abs(D*fo_dis/(D-RMag*delta)-D*fo_dis/(D+RMag*delta));

for i=2:1:dep_num
    
    imlayer=~(dep_im-dep_value(i)).*orig_im;
    if abs(fo_dis-dep_value(i))<=DoF/10 || dep_value(i)<100
                 arti_im=arti_im+imlayer;
    else
         CoC=f*D*abs(1-fo_dis/dep_value(i))/(fo_dis-f);
         k=round(CoC/delta)+1;
         Bfilter = fspecial('gaussian',k,0.5*k);
         imlayer=cat(3,conv2(imlayer(:,:,1),Bfilter,'same'),conv2(imlayer(:,:,2),Bfilter,'same'),conv2(imlayer(:,:,3),Bfilter,'same'));
         arti_im=arti_im+imlayer;
    end
    
end

% im_dina=['D:/5132018/3DImage_TrainingSet/TrainData/image',num2str(j),num2str(ii),'.png'];
% imwrite(uint8(arti_im),im_dina); B.astype(np.int)
    
toc;
wholetime=toc
% figure(1)
% imagesc(uint8(orig_im));
% figure(2)
% imagesc(dep_im);
% figure(3)
% imagesc(uint8(arti_im));   


%DepFilter=fspecial('disk',10);

% max_dep=max(max(dep_im));
% min_dep=min(min(dep_im));
% RoD_im=max_dep-min_dep;  %Range of depth of image
% tic;
% for j=1:num_fstack
%     %fo_dis=min_dep+j*RoD_im/(num_fstack+1);
%     fo_dis=2300;
%     RMag=(fo_dis-f)/f; %reciprocal of magnification
%     y_im=2*fo_dis*tand(FoV_v/2)/RMag;x_im=2*fo_dis*tand(FoV_h/2)/RMag; %image height in y and x dimension
%     delta=y_im/size(images,1);% pixel size    
%     DoF=abs(D*fo_dis/(D-RMag*delta)-D*fo_dis/(D+RMag*delta));
%     for m=1:M
%         for n=1:N
%             if abs(fo_dis-dep_im(m,n))<=DoF/100000
%                  arti_im(m,n,:)=arti_im(m,n,:)+orig_im(m,n,:);
%             else
%                 for c=1:C                
%                   CoC=f*D*abs(1-fo_dis/dep_im(m,n))/(fo_dis-f);
%                   Out_foc=abs(f*fo_dis/(fo_dis-f)-f*dep_im(m,n)/(dep_im(m,n)-f)); % Defocus distance
%                   k=round(CoC/delta);
%                     if mod(k,2)==0
%                      k=k+1;
%                     end
%                   l=floor(k/2);
% %                   row_in=-l:l;
% %                  blurSpot=zeros(k);
%                   Bfilter = fspecial('gaussian',k,10*k);
% %                   if l>0
% %                   Bfilter=fspecial('disk',l);
% %                   else
% %                       Bfilter=1;
% %                   end
%                   blurSpot=orig_im(m,n,c)*Bfilter;
%                   m1=(m-l>0)+(m-l<1)*(l-m+2);
%                   m2=(m+l<M+1)*k+(m+l>M)*(k-m-l+M);
%                   n1=(n-l>0)+(n-l<1)*(l-n+2);
%                   n2=(n+l<N+1)*k+(n+l>N)*(k-n-l+N); 
% %                   t1=m-l+m1-1;
% %                   t2=m+l-k+m2;
% %                   t3=n-1+n1-1;
% %                   t4=n+l-k+n2;
%                   
%                    arti_im(m-l+m1-1:m+l-k+m2,n-l+n1-1:n+l-k+n2,c)=arti_im(m-l+m1-1:m+l-k+m2,n-l+n1-1:n+l-k+n2,c)...
%                        +blurSpot(m1:m2,n1:n2);  
% %                   
%                 end    
%              end
%         end
%     end
%     toc;
%     wholetime=toc
% figure(1)
% imagesc(uint8(orig_im));
% figure(2)
% imagesc(depths(:,:,1));
% figure(3)
% imagesc(uint8(arti_im));   
% end