function RadGy = RadiusGyration(current_xyz,current_mol)

     centerOfMass = sum(current_xyz .* current_mol(:, 2)) / sum(current_mol(:, 2));
     relativePositions = current_xyz(:, 1:3) - centerOfMass;
%% 2. 计算每个原子相对于质心的质量乘积
    massProducts = current_mol(:, 2) .* sum(relativePositions.^2, 2);

% 3. 计算关于旋转轴的惯性张量
    inertiaTensor = zeros(3, 3);
    for i = 1:size(current_mol, 1)
        inertiaTensor = inertiaTensor + current_mol(i, 2) * ...
            [relativePositions(i, 2)^2 + relativePositions(i, 3)^2, -relativePositions(i, 1) * relativePositions(i, 2), -relativePositions(i, 1) * relativePositions(i, 3);
            -relativePositions(i, 1) * relativePositions(i, 2), relativePositions(i, 1)^2 + relativePositions(i, 3)^2, -relativePositions(i, 2) * relativePositions(i, 3);
            -relativePositions(i, 1) * relativePositions(i, 3), -relativePositions(i, 2) * relativePositions(i, 3), relativePositions(i, 1)^2 + relativePositions(i, 2)^2];
    end

% 4. 对惯性张量进行对角化
    [eigenVectors, eigenValues] = eig(inertiaTensor);

% 5. 获取主轴和对应的主惯性矩
    principalAxes = diag(eigenValues);
    [~, idx] = sort(principalAxes, 'descend');
    principalAxes = principalAxes(idx);
    eigenVectors = eigenVectors(:, idx);

% 6. 计算分子的旋转半径
    totalMass = sum(current_mol(:, 2));
    RadGy = sqrt((principalAxes(1) + principalAxes(2)) / totalMass);

end