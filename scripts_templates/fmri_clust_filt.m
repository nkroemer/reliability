function [Y,XYZ,C] = fmri_clust_filt(Y,k)
% identification of voxel cluster size k in Y

 XYZ = {};
 
 % Find all non-NaN or non-zero voxels in the volume
 if any(isnan(Y))
   good_vox_idxs = find(~isnan(Y));
   multiplier = NaN;
 else
  good_vox_idxs = find(Y);
   multiplier = 0;
 end

fprintf('Finding clusters of size >=%d in %d good voxels\n', k, length(good_vox_idxs));

 % Create an XYZ list
 [XYZ{1:3}] = ind2sub(size(Y),good_vox_idxs);
 XYZ = cat(2,XYZ{:})';
 
 % Get a list of cluster memberships
 A = spm_clusters(XYZ);
 
 % Create a list of indices that exceed cluster threshold
 Q = [];
 C = [];
 nc = 0;
 for i = 1:max(A)
   j = find(A == i);
   if length(j) >= k
     Q = [Q j]; 
    nc = nc+1;
     C(nc).size = length(j);
     C(nc).XYZ = XYZ(:,j);
     curr_good_idxs = good_vox_idxs(j);
     [C(nc).maxvoxval maxidxs] = max(Y(curr_good_idxs));
     C(nc).maxvoxidx_in_vol = curr_good_idxs(maxidxs);
     C(nc).maxvoxidx_in_clust = find(Y(curr_good_idxs)==C(nc).maxvoxval);
     C(nc).voxidx_in_vol = good_vox_idxs(j);
 end
 end
 
% Weed out voxels that didn't pass the threshold
fprintf('Found %d voxels belonging to %d clusters\n', length(Q), nc);
XYZ = XYZ(:,Q);
 good_vox_idxs = good_vox_idxs(Q);
 
% Filter the data volume
Y1 = zeros(size(Y))*multiplier;
Y1(good_vox_idxs) = 1;
 
Y = Y.*Y1;

return