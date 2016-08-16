function MPToVTK(NURBS, d, GNum, ParaPts, filename, fieldname)
% MPToVTK(NURBS, d, GNum, ParaPts, filename, fieldname)
str1 = cat (2,'<?xml version="1.0"?> \n', ...
    '<VTKFile type="Collection" version="0.1"> \n', ...
    '<Collection> \n');

str2 = cat (2, '<DataSet part="%d" file="%s.vts"/> \n');

str3 = cat (2, ...
    '</Collection>\n', ...
    '</VTKFile> \n');

if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
    pvd_filename = cat (2, filename, '.pvd');
else
    pvd_filename = filename;
    filename = filename (1:end-4);
end

fid = fopen (pvd_filename, 'w');
if (fid < 0)
    error ('MPToVTK: could not open file %s', pvd_filename);
end

fprintf (fid, str1);
for iptc = 1 : numel(NURBS)
    filename_patch = cat (2, filename, '_', num2str (iptc));
    fprintf (fid, str2, iptc, filename_patch);
    SPToVTK(NURBS{iptc}, d(GNum{iptc}), ParaPts, filename_patch, fieldname);
end
fprintf (fid, str3);

fclose (fid);
end