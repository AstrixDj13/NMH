function [ h, hu, hv ] = apply_boundary_conditions(constants, grid, bconds, h, hu, hv)
%APPLY_BOUNDARY_CONDITIONS

%% ---- Set boundary conditions -------------------------------------------
for k = 1:length(bconds.bwest)
    if strcmp( bconds.bwest{k}, 'PER' )
        % ---- Periodic boundary condition on west boundary ---------------
        h(1,:)  = h(end-1,:);
        hu(1,:) = hu(end-1,:);
        hv(1,:) = hv(end-1,:);
        
    elseif strcmp( bconds.bwest{k}, 'HFIX' )
        % ---- Fixed flow depth on west boundary --------------------------
        h(1,:) = bconds.hwest; %2 * bconds.hwest - h(2,:);
        
    elseif strcmp( bconds.bwest{k}, 'HUFIX' )
        % ---- Fixed discharge on west boundary ---------------------------
        hu(1,:) = bconds.huwest; %2 * bconds.huwest - hu(2,:);
        hv(1,:) = 0;
        
    elseif strcmp( bconds.bwest{k}, 'WALL' )
        % ---- Wall boundary condition ------------------------------------
        h(1,:) = h(2,:);
        hu(1,:) = -hu(2,:);
        hv(1,:) = hv(2,:);
    else
        error('No known boundary condition')
    end
end

for k = 1:length(bconds.beast)
    if strcmp( bconds.beast, 'PER' )
        % ---- Periodic boundary condition on east boundary ---------------
        h(end,:) = h(2,:);
        hu(end,:) = hu(2,:);
        hv(end,:) = hv(2,:);
        
    elseif strcmp( bconds.beast{k}, 'HFIX' )
        % ---- Fixed flow depth on east boundary --------------------------
        h(end,:) = 2 * bconds.heast - h(end-1,:);
        
    elseif strcmp( bconds.beast{k}, 'HUFIX' )
        % ---- Fixed discharge on east boundary ---------------------------
        hu(end,:) = 2 * bconds.hueast - hu(end-1,:);
        hv(end,:) = 0;
        
    elseif strcmp( bconds.beast{k}, 'HUEXT' )
        % ---- Extrapolate value to east boundary -------------------------
        % Compute Froude number of the cell to the west
        Fr = hu(end-1,:) ./ sqrt( constants.g * h(end-1,:).^3 );
        
        for j = 2:length(grid.y)-1
            if Fr(j) > -1
                % linear extrapolation
                hu(end,j) = 2 * hu(end-1,j) - hu(end-2,j);
                hv(end,j) = 2 * hv(end-1,j) - hv(end-2,j);
            else
                % if supercritical flow from east to west, then we cannot
                % extrapolate
                warning('Cannot extrapolate')
            end
        end
        
    elseif strcmp( bconds.beast{k}, 'WALL' )
        % ---- Wall boundary condition ------------------------------------
        h(end,:) = h(end-1,:);
        hu(end,:) = -hu(end-1,:);
        hv(end,:) = hv(end-1,:);
        
    else
        error('No known boundary condition')
    end
end


for k = 1:length(bconds.bsouth)
    % ---- Periodic boundary condition on south boundary ------------------
    if strcmp( bconds.bsouth{k}, 'PER' )
        h(:,1) = h(:,end-1);
        hu(:,1) = hu(:,end-1);
        hv(:,1) = hv(:,end-1);
        
    elseif strcmp( bconds.bsouth{k}, 'HFIX' )
        % ---- Fixed flow depth on south boundary -------------------------
        h(:,1) = 2 * bconds.hsouth - h(:,2);
        
    elseif strcmp( bconds.bsouth{k}, 'HVFIX' )
        % ---- Fixed discharge on south boundary --------------------------
        hv(:,1) = 2 * bconds.hvsouth - hv(:,2);
        
    elseif strcmp( bconds.bsouth{k}, 'WALL' )
        % ---- Wall boundary condition ------------------------------------
        h(:,1) = h(:,2);
        hu(:,1) = hu(:,2);
        hv(:,1) = -hv(:,2);
        
    else
        error('No known boundary condition')
    end
end

for k = 1:length(bconds.bnorth)
    % ---- Periodic boundary condition on north boundary ------------------
    if strcmp( bconds.bnorth{k}, 'PER' )
        h(:,end)   = h(:,2);
        hu(:,end)   = hu(:,2);
        hv(:,end)   = hv(:,2);
        
    elseif strcmp( bconds.bnorth{k}, 'HFIX' )
        % ---- Fixed flow depth on north boundary -------------------------
        h(:,end) = 2 * bconds.hnorth - h(:,end-1);
        
    elseif strcmp( bconds.bnorth{k}, 'HVFIX' )
        % ---- Fixed discharge on north boundary --------------------------
        hv(:,end) = 2 * bconds.hvnorth - hv(:,end-1);
        
    elseif strcmp( bconds.bnorth{k}, 'HVEXT' )
        % ---- Extrapolate value to north boundary ------------------------
        % Compute Froude number of the cell to the south
        Fr = hv(:,end-1) ./ sqrt( constants.g * h(:,end-1).^3 );
        
        for i = 2:length(grid.x)-1
            if Fr(i) > -1
                % linear extrapolation
                hv(i,end) = 2 * hv(i,end-1) - hv(i,end-2);
            else
                % if supercritical flow from north to sourth, then we cannot
                % extrapolate
                warning('Cannot extrapolate')
            end
        end
        
    elseif strcmp( bconds.bnorth{k}, 'WALL' )
        % ---- Wall boundary condition ------------------------------------
        h(:,end) = h(:,end-1);
        hu(:,end) = hu(:,end-1);
        hv(:,end) = -hv(:,end-1);
        
    else
        error('No known boundary condition')
    end
end

% ---- Set corner values to NaN -------------------------------------------
h(1,1) = NaN;
h(1,end) = NaN;
h(end,1) = NaN;
h(end,end) = NaN;

hu(1,1) = NaN;
hu(1,end) = NaN;
hu(end,1) = NaN;
hu(end,end) = NaN;

hv(1,1) = NaN;
hv(1,end) = NaN;
hv(end,1) = NaN;
hv(end,end) = NaN;
end


