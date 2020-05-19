function [divs, rots, harms, glueMesh] = ...
    helmHodgeDECRectGridPullback(cutM, facevf, varargin)
% helmHodgeDECRectGridPullback(cutM, facevf, varargin)
%   Perform Hodge decomposition using discrete exterior calculus on
%   vector field defined on faces. Clip, denoise and/or smooth the results 
%   based on varargin.
%
% Parameters
% ----------
% cutM : struct with fields v,f
% facevf : #faces x 3 float array
%   vector field defined on faces
% varargin : keyword arguments (optional)
%   niterSmoothing : int or two ints 
%       how many smoothing steps to perform. 
%       If two values given, applies to div, rot separately
%       If method is 'both', applies to 
%       div-denoise, div-smooth, rot-denoise, rot-smooth
%   filterMethod : 'smooth' 'denoise' 'both'
%       What smoothing method to apply
%           denoise : average local extrema with neighboring face values
%           smooth : average all values with neighboring face values
%           both : first denoise, then smooth
%   Epsilon : float
%       small value setting threshold for local extremum for method == denoise
%   clipDiv : list of two floats
%       values to clip the rotation field
%   clipRot : list of two floats
%       values to clip the curl field
%
% Returns
% -------
% divs : struct
% rots : struct 
% harms : struct
% glueMesh : struct
%   glued rectilinear mesh
%
% NPMitchell 2020


% Method options
niter_div = 0 ;
niter_rot = 0 ;
niterU2d_div = 0 ;
niterU2d_rot = 0 ;
niterU2d_harm = 0 ;
clipDiv = [-Inf, Inf] ;  
clipRot = [-Inf, Inf] ;  
method = 'smooth' ;     % options: smooth, denoise
eps = 1e-16 ;
preview = false ;
for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double')
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
        
    if ~isempty(regexp(varargin{i}, '^[Nn]iter[Ss]moothing', 'match'))
        niter = varargin{i+1} ;
        % Allow for different niters for divergence and rotation
        if length(niter) > 1
            niter_div = niter(1) ;
            niter_rot = niter(2) ;
        else
            niter_div = niter ;
            niter_rot = niter ;
        end
    end    
    
    if ~isempty(regexp(varargin{i}, '^[Nn]iter[Ss]moothing[Uu]2[Dd]', 'match'))
        niter = varargin{i+1} ;
        % Allow for different niters for divergence and rotation
        if length(niter) > 1
            niterU2d_div = niter(1) ;
            niterU2d_rot = niter(2) ;
            niterU2d_harm = niter(3) ;
        else
            niterU2d_div = niter ;
            niterU2d_rot = niter ;
            niterU2d_harm = niter ;
        end
    end    
    
    if ~isempty(regexp(varargin{i}, '^[Cc]lip[Dd]iv', 'match'))
        clipDiv = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Cc]lip[Rr]ot', 'match'))
        clipRot = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Mm]ethod', 'match'))
        method = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Ee]psilon', 'match'))
        eps = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Pp]review', 'match'))
        preview = varargin{i+1} ;
    end    
end

% Unpack the cutMesh
FF = cutM.f ;
V2D = cutM.u ;
v3drs = cutM.v ;
nU = cutM.nU ;
nV = cutM.nV ;

% % compute COM for each triangle in 2d and 3d --> note that we do this
% % before gluing so that the 2d case is not messed up, and the 3d case is
% % the same either way
% bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
% bc = mean( bc, 3 ) ;
% bc2d = cat( 3, V2D(FF(:,1), :), V2D(FF(:,2), :), V2D(FF(:,3), :) );
% bc2d = mean( bc2d, 3 ) ;

% num faces in each row, col is nfU, nfV
% faces are arranged as (nU-1)*(nV-1) * 2
% [~, faceIDgrid] = defineFacesRectilinearGrid(V2D, nU, nV) ;

% % Check that this makes sense
% trisurf(FF, v3drs(:, 1), v3drs(:, 2), v3drs(:, 3), divv, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none')

% hold on;
% inds = 1:21:length(bc) ;
% quiver3(bc(inds, 1), bc(inds, 2), bc(inds, 3), ...
%     vfsm(inds, 1), vfsm(inds, 2), vfsm(inds, 3), 1) 
% % quiver3(bc(inds, 1), bc(inds, 2), bc(inds, 3), ...
% %     v0t(inds, 1), v0t(inds, 2), v0t(inds, 3), 1) 
% axis equal

% Take divergence and curl
        
% TODO: glue the mesh back together, FF will change
% cutMC is cutM that is Closed at the seam
[glueMesh, glue2cut] = glueCylinderCutMeshSeam(cutM) ;  
DEC = DiscreteExteriorCalculus( glueMesh.f, glueMesh.v ) ;

% Now resolve the vector field for decomposition
[v0n, v0t, v0t2d, jac3d_to_2d, ~, ~, dilation] = ...
    resolveTangentNormalVelocities(FF, v3drs, facevf, V2D, 1:length(FF)) ;

divv = DEC.divergence(v0t) ;
rotv = DEC.curl(v0t) ;

% Clip the fields using supplied bounds for valid values
divv(divv < clipDiv(1)) = clipDiv(1) ;
divv(divv > clipDiv(2)) = clipDiv(2) ;
rotv(rotv < clipRot(1)) = clipRot(1) ;
rotv(rotv > clipRot(2)) = clipRot(2) ;

% Perform Helmholtz-Hodge decomposition
[divU, rotU, harmU, scalarP, vectorP] = ...
    DEC.helmholtzHodgeDecomposition(facevf) ;

% Preview divergence
if preview 
    tmp = reshape(divv, [nU, (nV-1)]) ;
    imagesc(tmp)
    colormap bwr 
    colorbar
    caxis([-0.4, 0.4])
    title('Divergence, before denoising/smoothing')
    pause(1)
end

% LAPLACIAN SMOOTHING on vertices (divergence field)
fixed_verts = [] ;  % note: could use boundaries here, seems unnecessary
% Options are: laplacian_smooth(V,F,L_method,b,lambda,method,S,max_iter)
%
% todo: denoise on vertices here
%
divvsm = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
    0.05, 'explicit', divv, niter_div) ;

% View results on divergence
if preview 
    tmp = reshape(divvsm, [nU, (nV-1)]) ;
    imagesc(tmp)
    colormap bwr 
    colorbar
    caxis([-0.4, 0.4])
    title('Divergence, after denoising/smoothing')
    pause(1)
end

% LAPLACIAN SMOOTHING on faces (rotational field)
% weight curl by face area
% Note: hodge dual (star) is applied at end of Curl(), so no area exists in
% the curl/rotational field, so just average with any weighting we like (we
% think).

% First inspect the field
if preview 
    tmp = reshape(rotv, [nU-1, (nV-1)*2]) ;
    imagesc(tmp)
    colormap bwr 
    colorbar
    caxis([-0.4, 0.4])
    title('Curl, before denoising/smoothing')
    pause(1)
end

rotvsm = rotv ;
% SMOOTHING -- DENOISE, SMOOTH, or ITERATE BOTH
if strcmp(method, 'denoise')
    % apply harmonic smoothing several times to extrema (DENOISE)
    for qq = 1:niter_rot
        rotvsm = laplacian_smooth_faces(rotvsm, glueMesh, ...
            'Weight', 'Area', 'Method', 'denoise', 'Epsilon', eps) ;

        if preview 
            tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
            imagesc(tmp)
            colormap bwr 
            colorbar
            caxis([-0.4, 0.4])
            title(['Curl, denoising #' num2str(qq)])
            pause(3)
        end
    end
elseif strcmp(method, 'smooth')
    % apply harmonic smoothing several times (SMOOTHING)
    for qq = 1:niter_rot
        rotvsm = laplacian_smooth_faces(rotvsm, glueMesh,...
            'Weight', 'Area', 'Method', 'smooth', 'Epsilon', eps) ;

        % Check it
        if preview
            tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
            imagesc(tmp)
            colormap bwr 
            colorbar
            caxis([-0.4, 0.4])
            title(['Curl, smoothing #' num2str(qq)])
            pause(1)
        end
    end
elseif strcmp(method, 'both')
    % apply denoising and harmonic smoothing several times (BOTH)
    for qq = 1:niter_rot
        rotvsm = laplacian_smooth_faces(rotvsm, glueMesh,...
            'Weight', 'Area', 'Method', 'denoise', 'Epsilon', eps) ;

        % Check it
        if preview
            tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
            imagesc(tmp)
            colormap bwr 
            colorbar
            caxis([-0.4, 0.4])
            title(['Curl, denoising #' num2str(qq)])
            pause(1)
        end
        
        rotvsm = laplacian_smooth_faces(rotvsm, glueMesh,...
            'Weight', 'Area', 'Method', 'smooth', 'Epsilon', eps) ;
        
        % Check it
        if preview
            tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
            imagesc(tmp)
            colormap bwr 
            colorbar
            caxis([-0.4, 0.4])
            title(['Curl, smoothing #' num2str(qq)])
            pause(1)
        end
        
    end
else
    error('SmoothingMethod not recognized')
end

% tmp = reshape(rotv, [nU-1, (nV-1)*2])

% add nU points back to divv from phi=0 by duplication --> convert back to
% cut mesh indices
divvCut = divvsm(glue2cut) ;
% note: rot is on faces, and faces are preserved
rotvCut = rotvsm ; 

% Note: rotU, divU, and harmU are on faces, and faces are preserved

% Coarse grain
% bcx = bc(:, 1) ;
% bcy = bc(:, 2) ;
% bcz = bc(:, 3) ;
% bcgridx = imresize(bcx(faceIDgrid), 1/qsub) ;
% bcgridy = imresize(bcy(faceIDgrid), 1/qsub) ;
% bcgridz = imresize(bcz(faceIDgrid), 1/qsub) ;
% divUgrid = divU(faceIDgrid) ;
% rotUgrid = rotU(faceIDgrid) ;
% harmUgrid = harmU(faceIDgrid) ;

% Pullback Vector Fields to Domain of Parameterization ====================
divU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    divU2d(f,:) = jac3d_to_2d{FF(f)} * divU(f,:)' / dilation(f) ;
end
rotU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    rotU2d(f,:) = jac3d_to_2d{FF(f)} * rotU(f,:)'  / dilation(f) ;
end
harmU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    harmU2d(f,:) = jac3d_to_2d{FF(f)} * harmU(f,:)' / dilation(f) ;
end

% Smooth output 2d fields
for dim=1:2
    divU2d(:, dim) = laplacian_smooth_faces(divU2d(:, dim), glueMesh, ...
        'Weight', 'Area', 'Method', method, 'Epsilon', eps, ...
        'niter', niterU2d_div) ;
    rotU2d(:, dim) = laplacian_smooth_faces(rotU2d(:, dim), glueMesh, ...
        'Weight', 'Area', 'Method', method, 'Epsilon', eps, ...
        'niter', niterU2d_rot) ;
    harmU2d(:, dim) = laplacian_smooth_faces(rotU2d(:, dim), glueMesh, ...
        'Weight', 'Area', 'Method', method, 'Epsilon', eps, ...
        'niter', niterU2d_harm) ;
end
    
% Store fields in structs
divs.divv = divvCut ;
divs.divU = divU ;
divs.divU2d = divU2d ;
rots.rotv = rotvCut ;
rots.rotU = rotU ;
rots.rotU2d = rotU2d ;
harms.harmU = harmU ;
harms.harmU2d = harmU2d ;

return

