function [neighb] = neighMGen( params, max_dist)
%  NEIGHMGEN generate adjacency matrix automatically from standard electrode cap location files
%  supported by eeglab.
%  Parameter  Value:
%  - max_dist is the radius of neighbourhood
%   Example:  
%       params{1} = '.\NMP_BrainVision_32.loc';
%       params{2} = 'filetype';
%       params{3} = 'autodetect';
%       neighbM = neighMGen(params, 0.6);       

[chans] = readlocs(params{:});
%[tmp tmp2 chans]  = eeg_checkchanlocs(chans);
chaninfo = [];
chaninfo.filename = params{1};
% backup file content etc...
tmptext         = loadtxt( chaninfo.filename, 'delim', [], 'verbose', 'off', 'convert', 'off');
chaninfo.filecontent = strvcat(tmptext{:});

% set urchan structure
urchans = chans;
for index = 1:length(chans)
    chans(index).urchan = index;
end;                        

if ~isfield(chans, 'datachan')
    chans(1).datachan = [];
end;
for index = 1:length(chans)
    if isempty(chans(index).datachan)
        chans(index).datachan = 1;
    end;
end;                       
          
chans2 = convertlocs(chans, 'topo2all');
p_num = length(chans2);
xyz = zeros( p_num, 3 );
xyz(:,1) = [chans2(:).X]';  xyz(:,2) = [chans2(:).Y]';  xyz(:,3) = [chans2(:).Z]';

dist_mat = zeros(p_num, p_num);
for ind = 1: p_num 
    for ind2 = ind: p_num 
        vec = xyz(ind,:)-xyz(ind2,:);
        dist = sum(vec.^2);
        dist_mat(ind,ind2) = dist;
        dist_mat(ind2,ind) = dist;
    end
end
[sort_dist, sort_ind] = sort(dist_mat,  2, 'ascend');

neighb = zeros(p_num, p_num);
for ind  = 1:p_num 
    for ind2  = 1:p_num
    if   dist_mat(ind,ind2) <max_dist && ind~=ind2
        neighb(ind,ind2) = 1;   
    end
    end
end

%figure; imshow(conn_mat);