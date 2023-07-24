%% Nylander's Parametric Rose

parametric_rose(false)

%%
%Animate rose blooming
opening = linspace(2*pi,pi,50);
density = linspace(5,12,50);

s = parametric_rose(true);

for i = 1:length(scale)

    % Calculate surface
    [x,y,z] = rose_vals(1,density(i));

    % Assign values
    s.XData = x;
    s.YData = y;
    s.ZData = z;

    pause(0.05)
end

%% Nylander's Parametric Rose function

function [rose_surface] = parametric_rose(custom_colormap)
    
    % Parameters
    n = 800;
    
    A = 1.995653;
    B = 1.27689;
    C = 8;
    petalNum=3.6;
   
    % Grid
    r_ = linspace(0,1,n);
    theta_ = linspace(-2,20*pi,n);
    
    [r,theta] = ndgrid(r_,theta_);
    
    % Helper variables
    phi = (pi/2)*exp(-theta/(C*pi));
    corte_petalo = 1 - (1/2)*((5/4)*(1 - mod(petalNum*theta, 2*pi)/pi).^2 - 1/4).^2;
    curva_petalo = A*(r.^2).*(B*r - 1).^2.*sin(phi);
    petalo = corte_petalo.*(r.*sin(phi) + curva_petalo.*cos(phi));
    
    % Assign 3d coordinates
    X = petalo.*sin(theta);
    Y = petalo.*cos(theta);
    Z = corte_petalo.*(r.*cos(phi)-curva_petalo.*sin(phi));

    %Create custom colormap
    % Define the number of color shades
    num_shades = 128;
    % Create a custom colormap with shades of yellow for the rose flower effect
    rose_colormap = [linspace(1, 1, num_shades); linspace(0.6, 0, num_shades); linspace(0.4, 0, num_shades)]';

    rose_surface = surf(X,Y,Z,'LineStyle','none');
    if custom_colormap
        colormap(flipud(rose_colormap))
    end

    % Configure plot
    xlim([-1 1]) 
    ylim([-1 1])
    zlim([-0.4 1])
    % Set equal aspect ratios for the axes
    daspect([1 1 1])
    hold on
end

function [X,Y,Z] = rose_vals(opening, density)
    n = 800;
    A = 1.995653;
    B = 1.27689;
    C = 8;
    r=linspace(0,1,n);
    theta=linspace(-2,20*pi,n);
    [R,THETA]=ndgrid(r,theta);
    % define the number of petals we want per cycle. Roses have 3 and a bit.
    petalNum=3.6;
    x = 1 - (1/2)*((5/4)*(1 - mod(petalNum*THETA, 2*pi)/pi).^2 - 1/4).^2;
    phi = 1/opening*(pi/2)*exp(-THETA/(density*pi));
    y = A*(R.^2).*(B*R - 1).^2.*sin(phi);
    R2 = x.*(R.*sin(phi) + y.*cos(phi));
    
    X=R2.*sin(THETA);
    Y=R2.*cos(THETA);
    Z=x.*(R.*cos(phi)-y.*sin(phi));

end