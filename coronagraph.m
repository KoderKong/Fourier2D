clear
close all
[im1,im0] = makeimage([100 300],16,uint8(100));
width = [112 400]; % Annulus circles

rng default
k = [0.5 1 2 4 8 16]; % Row/col step
[MN,secs,bytes] = complexity(im1,width,k);
addfigures('secsvspels',MN,secs)
addfigures('bytesvspels',MN,bytes)
clear global

rng default
[Wb,im2] = occulter(im0,width);
[Xa,Xt] = opticalsystem(im1,Wb);
makefigure(im2,Xt.*Xt)
Ph = zeros(size(Xa));
Ph = optimize(Ph,Xa,Wb,100);

function [im1,im0] = makeimage(exopos,exodia,bright)
im0 = imread('Airy_disk_D65.png');
pos = [flip(exopos) exodia/2; exopos exodia/2]; % Twins
col = repmat(bright,1,3);
im1 = insertShape(im0,'FilledCircle',pos,'Color',col);
im1 = rgb2gray(im1); % Ignore color
im0 = rgb2gray(im0); % Source image
end

function [Wb,im] = occulter(im,width)
dims = size(im);
centre = round(flip(dims)/2);
radius = round(sort(width)/2);
mask = cell(1,2);
im = jetcolor(im); % False color
for k = 1:2
    mask{k} = zeros(dims,'uint8');
    pos = [centre radius(k)];
    mask{k} = insertShape(mask{k},'FilledCircle',pos);
    mask{k} = logical(rgb2gray(mask{k}));
    im = insertShape(im,'Circle',pos,'Color','cyan');
end
Wb = mask{1} | ~mask{2};
end

function im = jetcolor(im)
dims = size(im);
cmap = jet(256);
im = cmap(double(im)+1,:);
im = uint8(im*255);
im = reshape(im,[dims 3]);
end

function [Xa,Xt,Ph] = opticalsystem(im,Wb,Ph)
if nargin < 3
    Xr = rand(size(im));
    Ph = angle(fft2(Xr));
    assert(isvalid(Ph))
end
Xt = sqrt(double(im)/255);
Xt(Wb) = 0; % Blackout
Yt = fft2(Xt);
Ya = Yt.*complex(cos(-Ph),sin(-Ph));
Xa = real(ifft2(Ya));
end

function is = isvalid(Ph)
[M,N] = size(Ph);
[~,~,Ru] = dftmat(M);
[~,~,Rv] = dftmat(N);
Ph = mod(Ph+Ru*Ph*Rv,2*pi);
is = all(Ph == 0,'all');
end

function [W,K,R] = dftmat(N)
k = 0:N-1;
K = mod(k'*k,N);
Ph = K*(-360/N);
W = complex(cosd(Ph),sind(Ph));
if nargout > 2
    i = 1:N;
    j = [1 N:-1:2];
    R = sparse(i,j,1,N,N,N);
end
end

function makefigure(varargin)
frame = cell(1,nargin);
for k = 1:nargin
    switch ndims(varargin{k})
        case 3
            image(varargin{k}) % False (jet) color
        case 2
            image(varargin{k},'CDataMapping','scaled')
            set(gca,'CLim',[0 1])
    end
    colormap(jet(256))
    colorbar
    axis image off
    frame{k} = getframe(gcf);
    frame{k} = border(frame{k}.cdata);
    close
end
try
    imwrite([frame{:}],'groundtruth.jpg')
catch
    for k = 1:numel(frame)
        file = 'groundtruth%d.jpg';
        imwrite(frame{k},sprintf(file,k));
    end
end
end

function im = border(im)
white = intmax(class(im));
im([1 end],:,:) = white;
im(:,[1 end],:) = white;
end

function Ph = optimize(Ph,Xa,Wb,maxiter)
outfun = @(Ph,optval,state) oneiter(Ph,Xa,optval,state,maxiter);
options = optimoptions('fminunc',...
    'Display','iter','MaxIterations',maxiter-1,'OutputFcn',outfun,...
    'Algorithm','trust-region','SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(Xt,DPh) hessmfun(Xt,DPh,Wb));
Ph = fminunc(@(Ph) ssefun(Ph,Xa,Wb),Ph,options);
end

function [SSE,GSSE,Xt] = ssefun(Ph,Xa,Wb)
global BYTES %#ok<GVMIS>
Ya = fft2(Xa);
Yt = Ya.*complex(cos(Ph),sin(Ph));
Xt = real(ifft2(Yt));
Wf = ~Wb & Xt < 0;
Xe = (Wb | Wf).*Xt;
if nargout > 1
    MN = numel(Xt);
    Ye = fft2(Xe);
    GSSE = 2/MN*imag(conj(Yt).*Ye);
end
Xe = Xe(:);
SSE = Xe'*Xe;
if isscalar(BYTES)
    S = whos;
    BYTES = sum([S.bytes]);
end
end

function HMF = hessmfun(Xt,DPh,Wb)
global BYTES %#ok<GVMIS>
Yt = fft2(Xt);
W = Wb | (~Wb & Xt < 0);
Ye = fft2(W.*Xt);
[M,N] = size(Xt);
DPh = reshape(DPh,M,N,[]);
MN = numel(Xt);
DYt = 1j*(Yt.*DPh);
DXt = real(ifft2(DYt));
DYe = fft2(W.*DXt);
DYY = imag(conj(DYt).*Ye);
YDY = imag(conj(Yt).*DYe);
HMF = 2/MN*(DYY+YDY);
HMF = reshape(HMF,MN,[]);
if isscalar(BYTES)
    S = whos;
    BYTES = sum([S.bytes]);
end
end

function stop = oneiter(Ph,Xa,optval,state,maxiter)
persistent frame SSE
stop = false;
switch state
    case 'init'
        frame = struct('cdata',[],'colormap',[]);
        frame = repmat(frame,maxiter+1,1);
        SSE = NaN(maxiter+1,1);
    case 'iter'
        k = optval.iteration+1;
        SSE(k) = optval.fval;
        Yt = fft2(Xa).*complex(cos(Ph),sin(Ph));
        Xt = real(ifft2(Yt));
        frame(k) = oneframe(Xt.*Xt,SSE);
    case 'done'
        stop = true;
        k = optval.iteration+1;
        makevideo('coronagraph',frame(1:k))
        close
end
end

function frame = oneframe(Xtt,SSE)
[M,N] = size(Xtt);
maxiter = numel(SSE)-1;
ratio = [M/5 N/maxiter 1];
xlims = [0 maxiter];
ylims = [0 -5]+ceil(max(log10(SSE)));
image(xlims,ylims,Xtt,'CDataMapping','scaled')
set(gca,'DataAspectRatio',ratio)
set(gca,'YDir','normal')
set(gca,'CLim',[0 1])
colormap(jet(256))
colorbar
hold on
plot(0:maxiter,log10(SSE),'cyan-')
hold off
xlabel('Iteration')
ylabel('log_{10} SSE')
frame = getframe(gcf);
end

function makevideo(file,frame)
numFrame = numel(frame);
audio = load('gong');
audioLen = length(audio.y);
Fs = numFrame*audio.Fs/audioLen;
audioFrameLen = floor(audioLen/numFrame);
try
    videoFWriter = vision.VideoFileWriter(strcat(file,'.avi'),...
        'FrameRate',Fs,'VideoCompressor','DV Video Encoder',...
        'AudioInputPort',true); % Ignore audio compressor
catch
    videoFWriter = vision.VideoFileWriter(strcat(file,'.avi'),...
        'FrameRate',Fs,'AudioInputPort',true);
end
for k = 1:numFrame
    a = 1+(k-1)*audioFrameLen; % Start of frame
    b = a+audioFrameLen-1; % End of frame
    videoFWriter(frame(k).cdata,audio.y(a:b))
end
release(videoFWriter) % Close file
im1 = border(frame(1).cdata);
im2 = border(frame(numFrame).cdata);
imwrite([im1 im2],strcat(file,'.jpg'))
end

function [MN,secs,bytes] = complexity(im,width,step)
dotdot('analyze complexity',10)
num = numel(step);
MN = zeros(num,1);
secs = zeros(num,3);
bytes = zeros(num,3);
for i = 1:num
    dotdot(true)
    k = step(i);
    imk = imresize(im,1/k);
    Wb = occulter(imk,width/k);
    Xa = opticalsystem(imk,Wb);
    Ph = zeros(size(Xa));
    MN(i) = numel(imk);
    secs(i,:) = timeit(Ph,Xa,Wb,k);
    bytes(i,:) = spaceit(Ph,Xa,Wb);
end
dotdot(false)
end

function dotdot(varargin)
persistent count
persistent width
switch nargin
    case 1
        if varargin{1}
            count = count+1;
            fprintf('.')
            if count >= width
                count = 0;
                fprintf('\n')
            end
        else
            if count > 0
                count = 0;
                fprintf('\n')
            end
        end
    case 2
        fprintf('%s\n',varargin{1})
        count = 0;
        width = varargin{2};
end
end

function secs = timeit(Ph,Xa,Wb,step)
global BYTES %#ok<GVMIS>
BYTES = [];
secs = zeros(1,3);
reps = 20*step*step;
tstart = tic;
for j = 1:reps
    [~] = ssefun(Ph,Xa,Wb);
end
secs(1) = toc(tstart)/reps;
tstart = tic;
for j = 1:reps
    [~,~] = ssefun(Ph,Xa,Wb);
end
secs(2) = toc(tstart)/reps;
[~,~,Xt] = ssefun(Ph,Xa,Wb);
DPh = cat(3,Ph,Ph); % P == 2
tstart = tic;
for j = 1:reps
    [~] = hessmfun(Xt,DPh,Wb);
end
secs(3) = toc(tstart)/reps;
end

function bytes = spaceit(Ph,Xa,Wb)
global BYTES %#ok<GVMIS>
BYTES = 8; % Scalar 'double'
bytes = zeros(1,3);
[~] = ssefun(Ph,Xa,Wb);
bytes(1) = BYTES-8;
[~,~] = ssefun(Ph,Xa,Wb);
bytes(2) = BYTES-8;
[~,~,Xt] = ssefun(Ph,Xa,Wb);
DPh = cat(3,Ph,Ph); % P == 2
[~] = hessmfun(Xt,DPh,Wb);
bytes(3) = BYTES-8;
end

function addfigures(file,N,value)
loglog(N,value,'.-')
axis tight
grid on
xlabel('Number of pixels')
switch file
    case 'secsvspels'
        ylabel('Compute time (s)')
    case 'bytesvspels'
        ylabel('Memory space (B)')
    otherwise
        error('Invalid filename!')
end
legend('SSE only','SSE & gradient',...
    'HMF (2 pages)','Location','NW')
set(gcf,'PaperSize',[6 4])
set(gcf,'PaperPosition',[0 0 6 4])
print('-dpdf',file)
close
end
