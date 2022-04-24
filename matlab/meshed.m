% nump = 100;
% t = 0:pi/nump:2*pi;
% matrixMeshed = sin(t') *cos(t);
% surface(matrixMeshed)

t = 0:pi/50:10*pi;
plot3(sin(t),cos(t),t);

for i = 1:length(t)
    for j=1:length(t) 
        m2(i, j) = sin(t(i)) * cos(t(j));
    end
end
