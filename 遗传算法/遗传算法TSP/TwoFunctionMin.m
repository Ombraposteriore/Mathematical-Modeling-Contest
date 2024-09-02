%采用GA工具，求解2元函数的最大值
%hzj
%2023-08-11


xi = linspace(-6,2,300);
yi = linspace(-4,4,300);
[X,Y] = meshgrid(xi,yi);
Z = ps_example([X(:),Y(:)]);
Z = reshape(Z,size(X));
surf(X,Y,Z,'MeshStyle','none');hold on
colormap 'jet'
view(-26,43)
xlabel('x(1)')
ylabel('x(2)')
title('ps\_example(x)')

x = ga(@ps_example,2);
z = ps_example(x);
plot3(x(1),x(2),z,'r*');
hold off