% script to output the smooth W for GFP case
% fid = fopen('W_GFP_smooth','w');
fid = fopen('W_rsaL_peak_BF_smooth','w');

fprintf(fid,'%12.10e ',ttmp);

%WS = W(:,1:5:1000);

% for i = 1:
%  
%    tmp = WS(i,:);  
%    fprintf(fid, '%12.10e ', tmp);
%    fprintf(fid, '\n');
%    
% end

fclose(fid);