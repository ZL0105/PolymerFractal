function fracdim = FractalDimension(current_xyz)

        box_x = max(current_xyz(:,1))-min(current_xyz(:,1));
        box_y = max(current_xyz(:,2))-min(current_xyz(:,2));
        box_z = max(current_xyz(:,3))-min(current_xyz(:,3));
        ini_cube_size = max([box_x box_y box_z]);
        ini_center = mean(current_xyz);

        clear box_sizes box_counts

        for i6 = 1:1:10
            num_divisions = 2^(i6-1); % Break every dimension into 2 parts
            subcube_size = ini_cube_size/num_divisions;
            if subcube_size<1
                break;
            end
            cube_flags = zeros(num_divisions,num_divisions,num_divisions);
            for ix = 1:num_divisions                
                for iy = 1:num_divisions                    
                    for iz = 1:num_divisions
                        subcube_centerX = ini_center(1) - 0.5*ini_cube_size + (ix-0.5)*subcube_size;
                        subcube_centerY = ini_center(2) - 0.5*ini_cube_size + (iy-0.5)*subcube_size;
                        subcube_centerZ = ini_center(3) - 0.5*ini_cube_size + (iz-0.5)*subcube_size;
                        subcube_center = [subcube_centerX subcube_centerY subcube_centerZ];

                        subcube_min = subcube_center - 0.5*subcube_size;
                        subcube_max = subcube_center + 0.5*subcube_size;

                        i_count = 0;
                        for itemp = 1:1:length(current_xyz)
                            if current_xyz(itemp,1)>=subcube_min(1) && current_xyz(itemp,1)<=subcube_max(1)
                                if current_xyz(itemp,2)>=subcube_min(2) && current_xyz(itemp,2)<=subcube_max(2)
                                    if current_xyz(itemp,3)>=subcube_min(3) && current_xyz(itemp,3)<=subcube_max(3)
                                        i_count = i_count + 1;
                                    end
                                end
                            end
                        end
                        if i_count > 0 
                            cube_flags(ix,iy,iz) = 1;
                        end
                    end
                end
            end 
            box_sizes(i6) = subcube_size;
            box_counts(i6) = sum(cube_flags(:));
        end


        coefficients = polyfit(log(box_sizes),log(box_counts),1);
        fracdim = -coefficients(1);
end