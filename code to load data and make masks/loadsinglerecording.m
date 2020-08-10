function y = loadsinglerecording(filename)

        
    
foo = load(filename);
whichVariables = fieldnames(foo);

y = foo.recording;

end