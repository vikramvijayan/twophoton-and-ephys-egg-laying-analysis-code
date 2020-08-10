%script to process new DLC and conver PWM lasr to power

for i = 1:1:length(files)
    load(files(i).name);
    [recording] = convert_PWMlaser_to_power(recording);
    [recording] = processes_more_DLC_variables(recording);
    save(files(i).name,'recording');
end
