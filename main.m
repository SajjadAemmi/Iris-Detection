clc
close all
clear all
tic

img = imread('dataset/(1).bmp');

if ndims(img) == 3
    img = rgb2gray(img);
end
[W, L] = size(img);

%% Iris inner contour - Pupil

H = imhist(img);
%plot(H);
[~,indexes] = findpeaks(H);
omega = indexes(1);

s0 = L - 0.25 * W;

xs = 0.25 * W : 30 : s0;
ys = [0.25 * W, 0.5 * W, 0.75 * W];

k = 1 : 0.5 * omega;

index = 1;

for i2 = ys
    for i1 = xs
        
        SW = img(round(i2 - 0.25 * W + 1 : i2 + 0.25 * W)  ,round(i1 - 0.25 * W + 1 : i1 + 0.25 * W));
        
        for row = 1 : size(SW, 1)
            for col = 1 : size(SW, 2)
                
                Landa = omega + k;
                if SW(row, col) < Landa
                    BW(row, col) = 1;
                else
                    BW(row, col) = 0;
                end
            end
        end
        
        BW = logical(BW);
        BW = imfill(BW, 'holes');
        se = strel('disk',3);
        BW = imopen(BW, se);
        
        %         figure;
        %         imshow(a);
        %         rectangle('Position', [i1 - 0.25 * W,i2 - 0.25 * W , size(SW, 1), size(SW, 2)],'EdgeColor','red','LineWidth',2);
        %
%         figure;
%         imshow(BW);
        
        % get regions
        rp = regionprops(BW, 'BoundingBox', 'Area' , 'PixelList');
        
        if size(rp,1) > 0
            for j = 1 : size([rp.Area],2)
                
                x_min = min(rp(j).PixelList(:,1));
                y_min = min(rp(j).PixelList(:,2));
                x_max = max(rp(j).PixelList(:,1));
                y_max = max(rp(j).PixelList(:,2));
                
                b1 = max(x_max - x_min, y_max - y_min);
                b2 = min(x_max - x_min, y_max - y_min);
                
                % roundness test
                if (0.6 * b1 <= b2 && b2 <= b1)
                    
                    %rectangle('Position', rp(j).BoundingBox,'EdgeColor','red','LineWidth',2);
                    
                    % Radius detected object
                    round_object(index).r0 = abs(round(0.5 * (0.5 * (x_max - x_min) + 0.5 * (y_max - y_min))));
                    
                    % Number of pixel in detected object
                    K0 = rp(j).Area;
                    
                    x_offset = i1 - 0.25 * W;
                    y_offset = i2 - 0.25 * W;
                    
                    % center of detected object
                    round_object(index).x0 = round((sum(rp(j).PixelList(:,1) + x_offset)) / K0);
                    round_object(index).y0 = round((sum(rp(j).PixelList(:,2) + y_offset)) / K0);
                    index = index + 1;
                end
            end
        end
        close all
    end
end

[~, i] = max([round_object.r0]);
rp = round_object(i).r0;
xp = round_object(i).x0;
yp = round_object(i).y0;

% figure;
% imshow(img);
% viscircles([xp, yp], rp, 'Color','w');

%% Iris outer contour Detection

% Extract a sub image
xp2 = xp;
yp2 = yp;
rp2 = rp;

if(xp > rp * 4)
    x_start = xp - rp * 4;
    xp2 = rp * 4;
else
    x_start = 1;
end

if(L - xp >= rp * 4)
    x_end = xp + rp * 4;
else
    x_end = L;
end

if(yp > rp * 4)
    y_start = yp - rp * 4;
    yp2 = rp * 4;
else
    y_start = 1;
end

if(W - yp >= rp * 4)
    y_end = yp + rp * 4;
else
    y_end = W;
end

s = img(y_start : y_end, x_start : x_end);
% imshow(s);

S_o = imcomplement(s);
S_o = imfill(S_o, 'holes');
S_o = imcomplement(S_o);
% imshow(S_o);

s_0_hat = medfilt2(S_o,[15,15]);
% imshow(s_0_hat);

d1 = double(min(min(s_0_hat)));
d2 = double(max(max(s_0_hat)));

s_m = round(255 * ((double(s_0_hat) - d1) / (d2 - d1)));
s_m = uint8(s_m);
% imshow(s_m);

% figure;
% imshow(s_m);
% viscircles([xp2, yp2], rp2, 'Color','w');

%%
Gray_level_Intensity = smooth(double(s_m(xp2, yp2 : -1 : 1)));

figure;
subplot(2,2,1);
imshow(s_m);

subplot(2,2,2);
plot(1:size(Gray_level_Intensity,1), Gray_level_Intensity);
xlabel('Pixels along the radial segment');
ylabel('Gray Level Intensity');
title('Gray Level Transition along the Pupil to Irris Region')

Gradient = abs(smooth(gradient(double(Gray_level_Intensity))));
subplot(2,2,3);
plot(1:size(Gradient,1),Gradient);
xlabel('Pixels taken along the radial segment');
ylabel('Gradient');
title('Peak Reprensenting a Point Located at the Irris Inner Contour');

s_sld = imresize(s_m, 0.5);

xp_sld = xp2 / 2;
yp_sld = yp2 / 2;
rp_sld = rp2 / 2;

% Search for the outer (limbus) boundary of the iris, using the result for the inner (pupil) boundary

[W, L] = size(s_sld);
s_sld = im2double(s_sld);
min_r = round(rp_sld * 1.5);
max_r = round(rp_sld * 3.5);
h = [];
angels = 0 : 0.05 : 2 * pi;
R = min_r : max_r;
for r = R
    
    Sum = 0;
    for ang = angels
        y = round(yp_sld - cos(ang) * r);
        x = round(xp_sld + sin(ang) * r);
        if y < 1
            y = 1;
        elseif y > W
            y = W;
        end
        if x < 1
            x = 1;
        elseif x > L
            x = L;
        end
        Sum = Sum + s_sld(y,x);
    end
    h(end + 1) = Sum;
end

h = diff(h);

% BLURRING
f = fspecial('gaussian', 3, 2);
IDO = abs(convn(h, f, 'same'));
[~ , i] = max(IDO(:));
r_out = R(i);

r_out = 2 * r_out;

% plot
figure
imshow(s);
title('Iris Detection');
viscircles([xp2, yp2], rp2,'EdgeColor','w');
viscircles([xp2, yp2], r_out,'EdgeColor','w');
toc