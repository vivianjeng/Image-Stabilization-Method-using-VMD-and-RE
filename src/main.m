clear all;clc;
octaves = 4;
intervals = 2;
cache = 1;
file_name = 'frames'; % dir of output frame
v = VideoReader('data/IMG_3714_resize.mp4'); % input video
v_out = VideoWriter('output/3714_out1.avi'); % output video
adjust_direction = 'vertical';

idx = 1;

videos = {};
videos_original = {};

%% SIFT
while hasFrame(v)
    disp(idx);
    video = readFrame(v);
    videos_original{idx} = video;
    s1 = '/';
    s2 = [ int2str(idx) ];
    s3 = '.pgm';
    pgms = [file_name s1 s2 s3];
    video = imresize(video,0.5);
    imwrite(video, pgms, 'pgm');
    videos{idx} = pgmread( [pgms] )./255;
    [ pos{idx}, scale{idx}, orient{idx}, desc{idx} ] = SIFT( videos{idx}, octaves, intervals, ones(size(videos{idx})),  0.0001, 10.0, 2 ,idx);
    idx = idx+1;
end

picNum =length(pos);

fpos = {};

for i=1:picNum-1
    fpos{i} = [];
end

for i=1:(picNum-1)
    [fpos{i}] = matching(fpos{i},pos{i},desc{i},pos{i+1},desc{i+1});
end


%% Calculating GMV
vertical_shift = GMV(picNum, fpos,adjust_direction);

%%  VMD

f2 = vertical_shift;

% VMD time domain
T = picNum-1;
fs = 1/T;
t = (1:T)/T;
freqs = 2*pi*(t-0.5-1/T)/(fs);

% composite signal, including noise
f_hat = fftshift((fft(f2)));

% some sample parameters for VMD
alpha = 100;        % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)
K = 5;              % 5 modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;

%--------------- Run actual VMD code

[u, u_hat, omega] = VMD(f2, alpha, tau, K, DC, init, tol);


% For convenience here: Order omegas increasingly and reindex u/u_hat
[~, sortIndex] = sort(omega(end,:));
omega = omega(:,sortIndex);
u_hat = u_hat(:,sortIndex);
u = u(sortIndex,:);
linestyles = {'b', 'g', 'm', 'c', 'c', 'r', 'k'};

for k = 1:K
    figure('Name', ['Reconstructed mode ' num2str(K)]);
    plot(t,u(k,:), linestyles{k});   hold on;
    %if ~isempty(fsub)
    %    plot(t, fsub{min(k,length(fsub))}, 'k:');
    %end
    set(gca, 'XLim', [0 1]);
end

plot(u(1,:)+u(2,:)+u(3,:), 'b'); hold on;
plot(vertical_shift,'r');

%%  Relative Entropy

% plot the histogram of each mode and normalize
subhist = {};
hist_values = {};
for i=1:K
    subhist{i} = histogram((u(i,:)-min(u(i,:)))/(max(u(i,:))-min(u(i,:))),20);
    hist_values{i} = subhist{i}.Values/picNum;
end

% calculating relative entropy
a = zeros(5,5);
for i=1:5
    for j=1:5
        a(i,j) = 0;
        for k=1:length(hist_values{i})
            if hist_values{j}(k) ~= 0 && hist_values{i}(k) ~= 0
                a(i,j) = a(i,j) - hist_values{j}(k)*(log(hist_values{i}(k))-log(hist_values{j}(k)));
            end
        end
    end
end

fprintf('relative entropy:\n\n');
disp(a);

%% Adjust original images &  Write Video

output_videos = {};

open(v_out);
for k=2:picNum-1
    init = videos{k}*255; % input video
    [R, C] = size(init); 
    output_videos{k-1} = zeros(R, C); 
    delX = -(u(4,k)+u(5,k)); 
    delY = 0; 
    tras = [1 0 delX; 0 1 delY; 0 0 1]; % translate transform 

    for i = 1 : R
        for j = 1 : C
            temp = [i; j; 1];
            temp = tras * temp; 
            x = floor(temp(1, 1));
            y = floor(temp(2, 1));
            if (x <= R) & (y <= C) & (x >= 1) & (y >= 1)
                output_videos{k-1}(x, y) = init(i, j);
            end
        end
    end

    writeVideo(v_out,output_videos{k-1});
end
close(v_out);
