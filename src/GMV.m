function motion_vector = GMV(picNum, fpos, mode)
% Global Motion Vector
% Input Parameters:
% picNum  -number of video frames
% fpos    -SIFT matching points
%         -fpos{i}(j,k)- the ith frame matches the (i+1)th frame
%                        the jth points
%                        k=1 x-axis of ith frame, k=2 y-axis of ith frame
%                        k=3 x-axis of (i+1)th frame, k=4 y-axis of (i+1)th frame
% mode    - 'vertical' or 'horizontal'
%
% Output:
% motion_vector - average pixel of movement


motion_vector = [];
accumulate = 0;
if mode == 'vertical'
    for i=1:picNum-1
        count = 0;
        shift_sum = 0;
        for j=1:length(fpos{i})
            if fpos{i}(j,4) ~= fpos{i}(j,2)
                shift_sum = shift_sum + (fpos{i}(j,4)-fpos{i}(j,2));
                count = count+1;
            end
        end
        if count ~= 0
            accumulate = accumulate + shift_sum/count;
        end
            motion_vector = [motion_vector accumulate];
    end
elseif mode == 'horizontal'
    for i=1:picNum-1
        count = 0;
        shift_sum = 0;
        for j=1:length(fpos{i})
            if fpos{i}(j,3) ~= fpos{i}(j,1)
                shift_sum = shift_sum + (fpos{i}(j,3)-fpos{i}(j,1));
                count = count+1;
            end
        end
        if count ~= 0
            accumulate = accumulate + shift_sum/count;
        end
            motion_vector = [motion_vector accumulate];
    end
end

% plot the motion_vector
figure('Name', "Global Motion Vector");
plot(motion_vector);