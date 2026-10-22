function T_m = rfTransEPG(alpha,phi)
        % rfTransEPG returns the tranistion matrix describing an RF pulse of flip 
        % angle alpha (degrees) with a phase of phi (radians); phi is relative to
        % the x-axis; T_m is 3x3 matrix
        %
        % % Example Inversion Pulse
        % alpha = 180
        % phi = 0
        % T_m = [
        %         Fn      Fn*     Zn
        %         0       1       0     Fn
        %         1       0       0     Fn*
        %         0       0       -1    Zn
        % ]
        %
        % alpha = 180
        % phi = 90
        % T_m = [
        %         Fn      Fn*     Zn
        %         0       -1      0     Fn
        %         -1      0       0     Fn*
        %         0       0       -1    Zn
        % ]
        %

        T_m = [ cosd(alpha/2)^2                 exp(2i*phi)*sind(alpha/2)^2     -1i*exp(1i*phi)*sind(alpha); ...
                exp(-2i*phi)*sind(alpha/2)^2    cosd(alpha/2)^2                 1i*exp(-1i*phi)*sind(alpha); ...
                -1i/2*exp(-1i*phi)*sind(alpha)  1i/2*exp(1i*phi)*sind(alpha)    cosd(alpha)];

end
