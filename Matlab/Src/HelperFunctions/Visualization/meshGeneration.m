function meshGeneration(gridObject, displacement, application)
    [rows, cols] = size(displacement.x);
    
    x = gridObject.grid.x;
    y = gridObject.grid.y;
    
    Ux = displacement.x;
    Uy = displacement.y;
    
    xMesh = zeros(rows, cols);
    yMesh = zeros(rows, cols);
    for i = 1: rows
        for j = 1: cols
            if(application == "Source")
                 xpos = ceil(x(i, j) - ceil(Ux(i,j)));
                 ypos = ceil(y(i,j) - ceil(Uy(i,j)));
            else
                 xpos = ceil(x(i, j) + ceil(Ux(i,j)));
                 ypos = ceil(y(i,j) + ceil(Uy(i,j)));
            end
           
            if(xpos > 0 && xpos <= rows && ypos >0 && ypos <= cols)
                xMesh(i, j) = xpos;
                yMesh(i, j) = ypos;
            end
        end
    end
    figure;
    mesh(xMesh, yMesh, xMesh * 0,'FaceLighting','gouraud','LineWidth',0.8);
    grid on; colormap gray;
    xlabel("x"); ylabel("y"); zlabel("z");
    title("Deformation Fields Applied To " + application + " Coordinate Plane");
end