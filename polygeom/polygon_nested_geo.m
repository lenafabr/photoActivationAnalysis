function [boundregion,boundregionB] = polygon_nested_geo(pt,ptB, indle,indbody,meshsize,outfile)
% give a polygon defined by a 2xn list of points (ptA)
% and another nested polygon given by ptB
% create a .geo file to be read in by gmsh for generation of a mesh
% returns the number for the boundary region
% indle is a list of indices into pt defining the leading edge
% indbody is a list of indices into pt defining the cell body boundary

npt = size(pt,2);
nptB= size(ptB,2);

% ---- write out a geometry file for gmesh --------
% heading info for the geo file
fid = fopen(outfile,'w');
fprintf(fid,'// Geometry file generated with matlab\n\n');
fprintf(fid,'lc = %f;\n\n',meshsize*2);
%fprintf(fid,'lcb = %f;\n\n',meshsize);

% put points in the geo file
for pc = 1:npt
    fprintf(fid,'Point(%d) = {%f, %f, %f, lc};\n',pc,pt(1,pc),pt(2,pc),0);
end
for pc = 1:nptB
    fprintf(fid,'Point(%d) = {%f, %f, %f, lc};\n',pc+npt,ptB(1,pc),ptB(2,pc),0);
end    

% put lines in the geo file
ptwrap = [1:npt,1];
ptBwrap = [npt+1:npt+nptB,npt+1];

for pc = 1:npt
    fprintf(fid,'Line(%d) = {%d, %d};\n',pc,ptwrap(pc),ptwrap(pc+1));
end
for pc = 1:nptB
    fprintf(fid,'Line(%d) = {%d, %d};\n',pc+npt,ptBwrap(pc),ptBwrap(pc+1));
end

ind = npt+nptB+1;

% set up polygon surface in geo file
str = sprintf('\nLine Loop(%d) = {',ind);
str = [str sprintf('%d,',[1:npt -(npt+nptB):-(npt+1)])];
str = [str(1:end-1) '};\n'];
fprintf(fid,str);

ind = ind +1;
fprintf(fid,'\nPlane Surface(%d) = {%d} ;\n',ind,ind-1);

% Physical boundary
boundregion = ind+100;
fprintf(fid,'\nBoundNum = %d;\n',boundregion);
str = 'Physical Line(BoundNum) = {';
str = [str sprintf('%d,',indle)];
str = [str(1:end-1) '};\n'];
fprintf(fid,str);

fprintf(fid,'\nBoundNum2 = %d;\n',boundregion+1);
str = 'Physical Line(BoundNum2) = {';
str = [str sprintf('%d,',indbody)];
str = [str(1:end-1) '};\n'];
fprintf(fid,str);


fprintf(fid,'Physical Surface(%d) = {%d};\n',boundregion+1,ind-2);

boundregionB = ind+200;
fprintf(fid,'\nBoundNumB = %d;\n',boundregionB);
str = 'Physical Line(BoundNumB) = {';
str = [str sprintf('%d,',npt+1:npt+nptB)];
str = [str(1:end-1) '};\n'];
fprintf(fid,str);

fprintf(fid,'Physical Surface(%d) = {%d};\n',boundregionB+1,ind);

% Recombine into quadrangles
fprintf(fid,'\nRecombine Surface{%d};\n', ind);
fprintf(fid,'Mesh.Algorithm = 8;\n');

fclose(fid);

boundregions = [boundregion, boundregion+1, boundregionB];