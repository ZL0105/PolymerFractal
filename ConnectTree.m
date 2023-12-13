function coordinate = ConnectTree(data)
global border_length

imod = 1;
visited = zeros(size(data, 1), 1);
visited(1) = data(1,1);

for i = 1:1:length(visited)
    
    ref_atom = visited(i);
    iref = find(data(:,1)==ref_atom);

    ref_x = data(iref,2);
    ref_y = data(iref,3);
    ref_z = data(iref,4);

    for icn = 5:1:8
        if data(iref,icn) ~= 0
           if ~any(visited == data(iref,icn))

                icrt = find(data(:,1)==data(iref,icn));
                dist = abs(data(iref,2:4)-data(icrt,2:4));

                if dist(1) > 8 && data(iref,2) > data(icrt,2)
                    data(icrt,2) = data(icrt,2) + border_length;
                elseif dist(1) > 8 && data(iref,2) < data(icrt,2)
                    data(icrt,2) = data(icrt,2) - border_length;
                end
                if dist(2) > 8 && data(iref,3) > data(icrt,3)
                    data(icrt,3) = data(icrt,3) + border_length;
                elseif dist(2) > 8 && data(iref,3) < data(icrt,3)
                    data(icrt,3) = data(icrt,3) - border_length;
                end
                if dist(3) > 8 && data(iref,4) > data(icrt,4)
                    data(icrt,4) = data(icrt,4) + border_length;
                elseif dist(3) > 8 && data(iref,4) < data(icrt,4)
                    data(icrt,4) = data(icrt,4) - border_length;
                end
                
                imod = imod + 1;
                visited(imod) = data(iref,icn);

            end
        end
    end
end

coordinate = data(:,2:4);
