current_path = pwd;
addpath(strcat(current_path, '/core'), ...
        strcat(current_path, '/scripts'));


welcome_message = fileread('welcome_message.txt');
disp(welcome_message)