function [B, src_indices, dest_indices] = noncircshift(A, offsets)
%Like circshift, but shifts are not circulant. Missing data are filled with
%zeros.
%
%  [B,src_indices,dest_indices]=noncirchift(A,offsets)
%
%B is the resulting array and the other outputs are such that
%
%  B(dest_indices{:})=A(src_indices{:})
%
    siz=size(A);
    N=length(siz);
    if length(offsets)<N
       offsets(N)=0; 
    end
    B=zeros(siz);
    indices=cell(3,N);
    for ii=1:N
          for ss=[1,3]
           idx=(1:siz(ii))+(ss-2)*offsets(ii);
            idx(idx<1)=[];
            idx(idx>siz(ii))=[];
           indices{ss,ii}=idx;
          end
      end
    src_indices=indices(1,:);
    dest_indices=indices(3,:);
    B(dest_indices{:})=A(src_indices{:});
end