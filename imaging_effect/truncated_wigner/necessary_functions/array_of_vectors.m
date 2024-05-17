function vec_cloned = array_of_vectors(vec,array_size,vec_dim)

reshape_mat = ones(1,length(array_size));
reshape_mat(vec_dim)=length(vec);
vec_at_correct_dim = reshape(vec,reshape_mat);

clone_mat = array_size;
clone_mat(vec_dim) = 1;
vec_cloned = repmat(vec_at_correct_dim,clone_mat);

end