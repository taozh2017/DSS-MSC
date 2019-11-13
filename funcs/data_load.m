

function [X,gt] = data_load(data_id)

% ------------------------------- 
data_path = 'data/';

load([data_path,data_id, '.mat']);


