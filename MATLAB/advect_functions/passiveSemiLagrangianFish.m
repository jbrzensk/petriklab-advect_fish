function [Flux] = passiveSemiLagrangianFish(conc_matrix, Flux, idx, dir, current, dt, dx_m, dy_m, neighbors, grid_mask, area)
%% semiLagrangianFish.m
% -------------------------------------------------------------------------
% Simulates fish advection in a 2D field using a 
% semi-Lagrangian approach.
%
% This function calculates the change in predator population (Flux). 
% It uses semi-Lagrangian advection to trace fish backward along the 
% flow field and compute their new location.
%
% INPUTS:
%   conc_matrix       : [n x m] array of fish biomass  [g/m^3]
%   Flux              : [n x m] array of flux of biomass due to movement
%   idx               : index of the grid cell [m_i, n_i]
%   dir               : [4 x 1] array of direction of more food [logical]
%   current           : [4 x 1] apparent current (swim speed taken
%                       into account ) at idx [m/s]
%   dt                : timestep [s]
%   dx_m, dy_m        : grid spacing in x and y directions [m]
%
% OUTPUTS:
%   Flux              : updated flux of fish array for grid cell idx.
%
% AUTHOR: [JARED BRZENSKI
% DATE  : 30-06-2025
% -------------------------------------------------------------------------
    % Directions matrix
    % directions = [
    %             -1, 0;  % Up
    %              1, 0;  % Down
    %              0, 1;  % Right
    %              0, -1];% Left  
    
    % Only go through loop if there is a current
    if all( current == 0 )
        return
    end

    tolerance = zeros(4,1);
    tolerance( current ~= 0 ) = 1e-3;

    [ n, m] = size( conc_matrix );

    % Current cell coordinates
    i = idx(1);
    j = idx(2);
    if ( size(dx_m,1) ~= 1 )
        distance = [
                dy_m(i,j);
                dy_m(i,j);
                dx_m(i,j);
                dx_m(i,j);];
    else
       distance = [
                dy_m;
                dy_m;
                dx_m;
                dx_m;];
    end

    valid_neighbor = check_valid_neighbor(neighbors, grid_mask);

    % Current cell concentration ( saves lookups later )
    cell_concentration = conc_matrix(i,j) .* area(i,j) ;

    proportion_out = valid_neighbor .* abs(current) .* dt .* flipud(distance) ./ area(i,j);
    %
    conc_moving_out = proportion_out .* cell_concentration;
    %
    total_leaving = sum(conc_moving_out);
    %
    % If too much wants to leave, cap at cell capacity - tolerance
    if total_leaving > cell_concentration
        % Scale amounts moving out to available concentration
        conc_moving_out = conc_moving_out .* (cell_concentration / total_leaving ) - tolerance;
        total_leaving = sum(conc_moving_out);
    end

    % Remove movers from the current cell
    Flux(i,j) = Flux(i,j) - total_leaving; % update for fish leaving this timestep
    %
    % Update fish positions around based on flux
    for k = 1:4
        ni = neighbors( k, 1);
        nj = neighbors( k, 2);
        Flux(ni, nj) = Flux(ni, nj) + conc_moving_out(k);
    end
end