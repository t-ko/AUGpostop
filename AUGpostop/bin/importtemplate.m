function filesplit = importtemplate(filename, delimiter)
    fid = fopen([ filename ], 'r');
    tline = fgets(fid);
    tlines = '';
    while ischar(tline)
        tlines = [tlines tline];
        tline = fgets(fid);
    end
    fclose(fid);
    filesplit = strsplit(tlines,delimiter);
end