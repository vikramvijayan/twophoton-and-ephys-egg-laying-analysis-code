list = dir('*/*.tif');

for i = 1:1:length(list)
    i
    name = list(i).name;
    path = list(i).folder;
    new_registration_script_nomem(name,path);
     % new_registration_script_nomem2(name,path);

end


