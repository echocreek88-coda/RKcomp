function [state, mass] = read3bodyCSV(filename)
% READ3BODYCSV Reads a 3-body setup CSV and outputs state & mass vectors.
%
% Input:
%   filename - path to CSV file formatted like:
%       BODY, Body1, Body2, Body3
%       MASS, m1, m2, m3
%       X, x1, x2, x3
%       Y, y1, y2, y3
%       Z, z1, z2, z3
%       Vx, vx1, vx2, vx3
%       Vy, vy1, vy2, vy3
%       Vz, vz1, vz2, vz3
%
% Output:
%   state - 1x6 cell vector: {[r1], [v1], [r2], [v2], [r3], [v3]}
%   mass  - 1x3 numeric vector [m1 m2 m3]

    % Read table, auto-detect numeric content
    T = readtable(filename, 'ReadVariableNames', true);

    % Extract numeric values ignoring first column (labels)
    data = T{:, 2:end};
    rows = string(T{:, 1});

    % Row lookup
    mIdx  = find(rows == "MASS");
    xIdx  = find(rows == "X");
    yIdx  = find(rows == "Y");
    zIdx  = find(rows == "Z");
    vxIdx = find(rows == "Vx");
    vyIdx = find(rows == "Vy");
    vzIdx = find(rows == "Vz");

    % Masses
    mass = data(mIdx, :);

    % Body 1
    r1 = [data(xIdx,1), data(yIdx,1), data(zIdx,1)];
    v1 = [data(vxIdx,1), data(vyIdx,1), data(vzIdx,1)];

    % Body 2
    r2 = [data(xIdx,2), data(yIdx,2), data(zIdx,2)];
    v2 = [data(vxIdx,2), data(vyIdx,2), data(vzIdx,2)];

    % Body 3
    r3 = [data(xIdx,3), data(yIdx,3), data(zIdx,3)];
    v3 = [data(vxIdx,3), data(vyIdx,3), data(vzIdx,3)];

    % Output nested row vector
    state = {r1, v1, r2, v2, r3, v3};
end
