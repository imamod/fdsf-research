function writeFile(filename, data, format)
    f = fopen(filename , 'w');
    fprintf(f, format, data);
    fclose(f);
end