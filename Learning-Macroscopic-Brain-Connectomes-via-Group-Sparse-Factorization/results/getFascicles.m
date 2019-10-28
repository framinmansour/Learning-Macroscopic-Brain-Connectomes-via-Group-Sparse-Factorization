function fascicles = getFascicles(Phi, coords, start, size)
    % Returns a sptensor of only the fascicles that pass through the cube
    % between start and start + [size, size, size]
    % Arguments:
    % Phi - the full tensor Phi containing all fascicles
    % coords - fe.roi.coords, map of indicies in Phi to actual coordinates
    % start - a 3 dimensional vector containing the corner of the cube
    %         which is minimal in the x, y, and z coordinates.
    % size - the side length of the cube being checked
    finish = start + [size, size, size];
    fascicleIndices = [];
    % Loop over all non zeroes in Phi
    for i = Phi.subs().'
        % Get XYZ coordinates of voxel of current non zero
        coord = coords(i(2),:);
        if coord(1) >= start(1) && coord(1) <= finish(1) && ...
            coord(2) >= start(2) && coord(2) <= finish(2) && ...
            coord(3) >= start(3) && coord(3) <= finish(3)
            % If coordinates between start and end points, add fascicle
            % index of current non zero to list
            fascicleIndices = [fascicleIndices, i(3)];
        end
    end
    % Remove all duplicates
    fascicleIndices = unique(fascicleIndices);
    % Return copy of Phi with only selected fascicles
    fascicles = Phi(:,:,fascicleIndices);
end