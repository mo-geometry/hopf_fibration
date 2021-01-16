addpath('src') 
nFrames = 80;
nCircles = 20;
cayBasis = 'left';%[ 'left', 'right' ]
videoString = 'noVideo'; % ['writeVideo','noVideo'];
% The cayley matrices
cay = cayley(cayBasis);
%% Sequence #0
hopf0 = sequence0(2*nCircles);
ctr(1)=plot2sphere1a(hopf0,'frames0',videoString); % ctr = 520;
%% Sequence #1
hopf1 = sequence1(2*nCircles,nFrames);
ctr(2)=plot2sphere1a(hopf1,'frames1',videoString); % ctr = 320;
%% Sequence #2
hopf2 = sequence2(2*nCircles,nFrames);
ctr(3)=plot2sphere1a(hopf2,'frames2',videoString); % ctr = 80;
%% Sequence #3
hopf3 = sequence3(2*nCircles,nFrames);
ctr(4)=plot2sphere1a(hopf3,'frames3',videoString); % ctr = 400;
%% Sequence #4
[hopf4a,hopf4b,hopf4c] = sequence4(nCircles,2*nFrames);
ctr(5)=plot2sphere3a(hopf4a,hopf4b,hopf4c,'frames4',videoString); % ctr = 160;
%% Sequence #5
[hopf5a,hopf5b,hopf5c] = sequence5(nCircles,nFrames);
ctr(6)=plot2sphere3a(hopf5a,hopf5b,hopf5c,'frames5',videoString); % ctr = 80;
%% Sequence #6
[hopf6a,hopf6b,hopf6c] = sequence6(nCircles,4*nFrames);
ctr(7)=plot2sphere3a(hopf6a,hopf6b,hopf6c,'frames6',videoString);
%% Sequence #7
[hopf7a,hopf7b,hopf7c] = sequence7(nCircles,2*nFrames);
ctr(8)=plot2sphere3a(hopf7a,hopf7b,hopf7c,'frames7',videoString);
% %% 
% for jj=1:8     
%      images = cell(ctr(jj),1);
%      for ii=1:ctr(jj)
%         images{ii} = imread(['frames',num2str(jj-1),'/img',num2str(ii,'%05.f'),'.png']);
%      end
%      % create the video writer with 1 fps
%      writerObj = VideoWriter(['frames',num2str(jj-1),'/hopf',num2str(jj-1),'.avi']);
%      % open the video writer
%      open(writerObj);
%      % write the frames to the video
%      for u=1:length(images)
%             frame = im2frame(images{u});
%             writeVideo(writerObj, frame);
%      end
%      % close the writer object
%      close(writerObj);
% end
