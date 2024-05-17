function pool_size = initialize_parpool(parallel)

% check if a pool object already exists
pool_obj = gcp('nocreate');

% initialize parpool
if parallel == 1
    % create a pool object if it does not exist
    if isempty(pool_obj)
        pool_obj = parpool;
    end
    % get pool size
    pool_size = pool_obj.NumWorkers;
    
    fprintf('parallel mode (size of worker pool: %d)\n', pool_size);
else
%     % delete pool object if it exists
%     if ~isempty(pool_obj)
%         delete(pool_obj);
%     end
    pool_size = 0;
    
    fprintf('serial mode (size of worker pool: %d)\n', pool_size);
end

end