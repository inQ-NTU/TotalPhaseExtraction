function outmat = remove_singleton(inmat,singleton_pos)
% function only removing the singleton-dimensions specified in singleton_pos
%
% TS 2016-06

size_in = size(inmat);

%check if the specified dimensions are really singleton dimensions
size_to_cut = size_in(singleton_pos);
if sum(size_to_cut) ~= length(singleton_pos)
    error('Wrong specification of singleton-pos! Not all specified dimensions are singleton-dimensions!')
end

no_dim = length(size(inmat));

size_ind_out = ones(no_dim,1);
size_ind_out(singleton_pos) = 0; 
size_ind_out = logical(size_ind_out);
size_out = size_in(size_ind_out);

outmat = reshape(inmat,size_out);

end