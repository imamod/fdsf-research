function [data] = readFile(filename)
    f = fopen(filename,'r');
    data = fscanf(f, '%f');
    data = data';
    fclose(f); 
end