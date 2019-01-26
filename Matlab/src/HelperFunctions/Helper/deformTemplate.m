function Template = deformTemplate(Template, U, gridObject)

x = gridObject.x;
y = gridObject.y;

for i = 2: length(U.x)
    for j = 2:length(U.y)
        xpos = ceil(x(i) - U.x(i,j));
        ypos = ceil(y(j) - U.y(i,j));
        if(xpos > 0 && ypos > 0 && xpos < x(end)+1 && ypos < y(end) + 1)
           Template(x(i), y(j)) = Template(xpos, ypos);
        end
    end
end

end

