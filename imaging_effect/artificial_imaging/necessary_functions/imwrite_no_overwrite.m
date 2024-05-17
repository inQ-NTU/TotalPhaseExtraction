function imwrite_no_overwrite(data,pic_name)

if ~exist(pic_name,'file')
    imwrite(data,pic_name,'tif','Compression','deflate');
else
    error(['Picture ',pic_name,' already exists!']);
end

end