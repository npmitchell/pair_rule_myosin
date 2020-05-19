function [faces, faceIDgrid] = defineFacesRectilinearGrid(uv, nU, nV) 
%DEFINEFACESRECTILINEARGRID(sp, nU, nV) Given rectilinear grid of points
%uv that is nU x nV, define a regular triangulation of faces. 
%
% Define faces with normals facing in (u x v) direction 
% Note that this is out of page when drawn as follows:
%
% V ^ 
%   |
%   |
%   |
%   o---------------> U
%
% but in imagesc() or imshow() this will be INTO the page:
%  
%   o---------------> U
%   |
%   |
%   |
% V v
%
% Parameters
% ----------
% uv : nU*nV x 2 float array
%   The positions of the mesh vertices
% nU : int
%   Number of vertices along the first dimension
% nV : int
%   Number of vertices along the second dimension
% 
% Returns
% -------
% faces : 2 * (nU-1) * (nV-1) x 3 int array
%   A simple triangulation of the rectilinear mesh (connectivityList)
%
% NPMitchell 2020

% check that uv has increasing u, then increasing v
assert(uv(1, 1) ~= uv(2, 1))
assert(uv(1, 2) == uv(2, 2))

% Define faces with normals facing in (u x v) direction 
% Note that this is out of page when drawn as follows:
%
% V ^ 
%   |
%   |
%   |
%   o---------------> U
%
% but in imagesc() or imshow() this will be INTO the page:
%  
%   o---------------> U
%   |
%   |
%   |
% V v
%
faces = zeros(2 * (nU-1) * (nV-1), 3) ;
kk = 1 ;
for i = 1:(nU-1)
    for j = 1:(nV-1)
        % first triangle goes from bottom row (right hugging bottom)
        faces(kk, :) = [(j-1)*nU + i, (j-1)*nU + i+1, j*nU + i+1] ;
        kk = kk + 1 ;
        % next to left of same row (left hugging top)
        faces(kk, :) = [(j-1)*nU + i, j*nU + i+1, j*nU + i] ;
        kk = kk + 1 ;
    end
end

% Check it
% triplot(tri, uv(:, 1), uv(:, 2))
% hold on;
% plot(uv(:, 1), uv(:, 2), 'o')

% Create faceIDgrid, indexing faces as a rectilinear grid
if nargout > 1
    % num faces in each row, col is nfU, nfV
    % faces are arranged as (nU-1)*(nV-1) * 2
    nfU = (nU - 1) ;
    nfV = (nV - 1) * 2;
    [xx, yy] = meshgrid(1:nfV, 1:nfU) ;
    faceIDgrid = yy + (xx-1) * nfU ;
end
