clc;
clear;

%constantes
k = 8.987*10^9;

%número de cargas
qn = input("¿Cuantas cargas habran? ");
while qn <= 0
    qn = input("Inserte un valor válido ");
end

%listas en 0 cargas
Q = zeros(1, qn);
%se llenan de nan para permitir posiciones en x = 0 y y = 0
qx = nan(1, qn);
qy = nan(1, qn);

%función preguntar posición
function [x] = ask_pos(L, i)
    prompt = ['Posición en ', L, ' de la carga ', num2str(i), ': '];
    x = str2double(input(prompt, 's'));
    while isnan(x)
        x = str2double(input("Introduce un valor válido: ", 's'));
    end
end

%preguntar carga y posición
for i = 1:qn
    prompt = ['Valor de la carga ', num2str(i), ': '];
    q = str2double(input(prompt, 's'));
    while isnan(q) || q == 0
        q = str2double(input("Introduce un valor válido: ", 's'));
    end
    x = ask_pos('x', i);
    y = ask_pos('y', i);

    %cargar datos a lista
    Q(i) = q;

    %verificar disponibilidad de posición
    x_qx = find(qx == x);
    y_qy = find(qy == y);
    while ~isempty(intersect(x_qx, y_qy))
        fprintf("\nLa posición ya esta ocupada por otra carga\n")
        x = ask_pos('x', i);
        y = ask_pos('y', i);
        x_qx = find(qx == x);
        y_qy = find(qy == y);
    end
    qx(i) = x;
    qy(i) = y;
end

%limites del gráfico

% evitar 0 como mínimo
function [n] = min0(x)
    if min(x) == 0
        n = -0.85;
    else
        n = min(x);
    end
end

% evitar 0 como máximo
function [n] = max0(x)
    if max(x) == 0
        n = 0.85;
    else
        n = max(x);
    end
end

min_y = min0(qy) - abs(min0(qy))*0.25;
min_x = min0(qx) - abs(min0(qx))*0.25;
max_y = max0(qy) + abs(max0(qy))*0.25;
max_x = max0(qx) + abs(max0(qx))*0.25;

%plano xy
n = 50;
x = linspace(min_x, max_x, n);
y = linspace(min_y, max_y, n);
[X, Y] = meshgrid(x, y);

%zeros campos y potencial
Ex = zeros(length(x), length(y));
Ey = zeros(length(x), length(y));
V = zeros(length(x), length(y));
grad_Vx = zeros(length(x), length(y));
grad_Vy = zeros(length(x), length(y));

%funciones

%función distancia entre dos puntos
function [d] = r(p, q)
    d = sqrt(dif(p(1),q(1))^2 + dif(p(2),q(2))^2);
    if d == 0 
       d = 10^-6;
    end
end

% función gradiente Potencial eléctrico total en x
function [Vxf] = grad_vx (p, q, k, n, x, y)
        Vxf_lista = zeros(1, n);
        for i = 1:n
            a = -k*q(i);
            b = ((dif(p(1),x(i))^2 + dif(p(2),y(i))^2)^(3/2));
            Vxf_lista(i) = a*(dif(p(1),x(i))/b);
        end
        Vxf = sum(Vxf_lista);
end

% función gradiente Potencial eléctrico total en y
function [Vyf] = grad_vy (p, q, k, n, x, y)
        Vyf_lista = zeros(1, n);
        for i = 1:n
            a = -k*q(i);
            b = ((dif(p(1),x(i))^2 + dif(p(2),y(i))^2)^(3/2));
            Vyf_lista(i) = a*(dif(p(2),y(i))/b);
        end
        Vyf = sum(Vyf_lista);
end

% función delta x
function [c] = dif(xf, xi)
    c = xf - xi;
end

% función potencial eléctrico total
function [C] = v (p,q,k,n,x,y)
    C_lista = zeros(1, n);
    for i = 1:n
        pos_q = [x(i), y(i)];
        C_lista(i) = k*q(i)/r(p,pos_q);
    end
    C = sum(C_lista);
end

% función campo eléctrico total en x
function [c] = ex (p,q,k,n,x,y)
    c = -1*grad_vx (p, q, k, n, x, y);
end

% función campo eléctrico total en y
function [c] = ey (p,q,k,n,x,y)
    c = -1*grad_vy (p, q, k, n, x, y);
end

%ciclo
for i = 1:length(x)
    for j = 1:length(y)
        p = [X(i, j), Y(i, j)];
        
        grad_Vx(i, j) = grad_vx (p, Q, k, qn, qx, qy);
        grad_Vy(i, j) = grad_vy (p, Q, k, qn, qx, qy);
        V(i, j) = v(p,Q,k,qn,qx,qy);
        Ex(i, j) = ex (p, Q, k, qn, qx, qy);
        Ey(i, j) = ey (p, Q, k, qn, qx, qy);
    end
end

% lineas de campo no conectadas
quiver(X,Y,Ex,Ey)
% mapeo de color por potencial de voltaje
hold on
pcolor(X, Y, V);
shading interp;
colormap('jet');
% encuadre
axis([min_x, max_x, min_y, max_y])
% lineas de campo conectadas
streamslice(X, Y, Ex, Ey);
hold off
