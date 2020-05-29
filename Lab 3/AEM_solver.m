function [T,Y] = AEM_solver(f, t0, tN, y0, h)
    tol = 1e-8;
    T = [t0];
    Y = [y0];
    j = 2;
    while (T(end)<tN)
        while 1
            y = f(T(j-1), Y(j-1))*h + Y(j-1);
            Z_MID = f(T(j-1), Y(j-1))*(h/2) + Y(j-1);
            Z = f(T(j-1)+(h/2), Z_MID)*(h/2) + Z_MID;
            if (abs(Z-y)<tol)
                Y = [Y Z+(Z-y)];
                T = [T T(j-1)+h];
                j = j + 1;
                break
            else
                h = 0.9*h*min(max(tol/abs(Z-y),0.3),2);
            end
        end
    end
