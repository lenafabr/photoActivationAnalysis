function boundregion = polygongeo(pt,meshsize,outfile)
% give a polygon defined by a 2xn list of points
% create a .geo file to be read in by gmsh for generation of a mesh
% returns the number for the boundary region

npt = size(pt,2);

% ---- write out a geometry file for gmesh --------
% heading info for the geo file
fid = fopen(outfile,'w');
fprintf(fid,'// Geometry file generated with matlab\n\n');
fprintf(fid,'lc = %f;\n\n',meshsize);

% put points in the geo file
for pc = 1:npt
    fprintf(fid,'Point(%d) = {%f, %f, %f, lc};\n',pc,pt(1,pc),pt(2,pc),0)
end
    
% put lines in the geo file
ptwrap = [1:npt,1];

for pc = 1:npt
    fprintf(fid,'Line(%d) = {%d, %d};\n',pc,ptwrap(pc),ptwrap(pc+1))
end

% set up polygon surface in geo file
str = sprintf('\nLine Loop(%d) = {',npt+1);
str = [str sprintf('%d,',1:npt)]
str = [str(1:end-1) '};\n']
fprintf(fid,str);

fprintf(fid,'\nPlane Surface(%d) = {%d} ;\n',npt+2,npt+1);

% Physical boundary
boundregion = npt+100;
fprintf(fid,'\nBoundNum = %d;\n',boundregion);
str = 'Physical Line(BoundNum) = {';
str = [str sprintf('%d,',1:npt)];
str = [str(1:end-1) '};\n'];
fprintf(fid,str);

fprintf(fid,'Physical Surface(%d) = {%d};\n',npt+101,npt+2);


% Recombine into quadrangles
fprintf(fid,'\nRecombine Surface{%d};\n', npt+2);
fprintf(fid,'Mesh.Algorithm = 8;\n');

fclose(fid);