function [fpos1] = matching(fpos1,pos1,desc1,pos2,desc2)

size(desc1(:,1),1);

distRatio = 0.6;

desc2t = desc2';
for i = 1 : size(desc1,1)
   dotprods = desc1(i,:) * desc2t;        % Computes vector of dot products
   [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

   % Check if nearest neighbor has angle less than distRatio times 2nd.
   if (vals(1) < distRatio * vals(2))
      match(i) = indx(1);
   else
      match(i) = 0;
   end
end

%s1='matchpoint';
%s2 = [ int2str(number) ];
%s3 = '.txt';
%s = [s1 s2 s3];
%TempName = 'temp3.txt' ;
%fid = fopen(TempName, 'w');

for i = 1: size(desc1,1)
  if (match(i) > 0)
         fpos1 = [ fpos1 ; pos1(i,1) , pos1(i,2) , pos2(match(i),1), pos2(match(i),2)];
         %fpos2 = [fpos2; pos2(match(i),1), pos2(match(i),2)];
         %fprintf( fid , '%d %d %d %d \n', pos1(i,1) , pos1(i,2) , pos2(match(i),1), pos2(match(i),2)); 
  end
%fclose(fid); 
end