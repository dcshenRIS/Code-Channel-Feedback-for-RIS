function [on_gird_value,on_gird_index] = On_gird2Off_gird(off_gird_value_vector,on_gird_sequence)
%% descriptions of this function
% This function is to trans. the off-gird angles to the on-gird angles with
% unavoidable errors of course.
% ---------------- input descriptions -------------------------------------
%   "off_gird_value_vector" is the real angles off-grids.
%   "on_gird_sequence" is the grid sequence vector.
% 
% ---------------- output descriptions ------------------------------------
%   "on_gird_value" is is the  angles on-grids.
%   "on_gird_index" is the angle grid indexes.
%% Note that there are some spelling errors (i.e., 'grid' and 'gird')
length=max(size(off_gird_value_vector,2),size(off_gird_value_vector,1));
for i=1:length
[~,on_gird_index(i)]=min(abs(on_gird_sequence-off_gird_value_vector(i)));
on_gird_value(i)=on_gird_sequence(on_gird_index(i));
end
end

