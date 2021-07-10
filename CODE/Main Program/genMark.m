function [ mark ] = genMark( iline, icolor, ishape )

line = {'-' '-.' '--'};
color = {'r' 'b' 'g' 'c' 'm' 'y' 'k'}; 
shape = {'s' '*' 'o'};

if iline > length(line)
   iline = 1;
end
if icolor > length(color)
   icolor = 1;
end
if ishape > length(shape)
   ishape = 1;
end

mark = char(strcat(line(iline), color(icolor), shape(ishape)));

end

