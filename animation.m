% Plot an algorithmic Dahlia flower
%
% The math for this was found as a rose in Ned's community blog:
%   https://blogs.mathworks.com/community/2021/02/12/happy-valentines-day/
% but Ned found it at Paul Nylander's math art site:
%   http://www.bugman123.com/Math/index.html
% Updated for labeled parameters, rounded petals, new colors, and simplified some expressions. 
% newbee attempt to do animation

ppr=12.6; % # petals per 1 revolution
nr = 30; % radius resolution
pr = 10; % petal resolution
np = 140; % total number of petals
pf = -1.2; % How much the ends of the petals tilt up or down


%% Animation
%Opening motion
for k = 0:0.1:1.1
    ol = [ 0.11 0.11+k ]; % How open is it? [ inner outer ]
    M = motion(ppr,nr,pr,np,pf,ol);
end
%rest before closing motion 
for freezeFrames = 0:6
    M = motion(ppr,nr,pr,np,pf,ol);
end
%Closing motion 
for k = 0:0.1:1.1
    ol = [ 0.11 1.1-k ]; % How open is it? [ inner outer ]
    M = motion(ppr,nr,pr,np,pf,ol);
end

%play the frames M as a movie 
movie(M,1,70)

%% Functions 
%draw the flower 
function [X,Y,Z,C]=drawFlower(ppr,nr,pr,np,pf,ol)
    pt = (1/ppr) * pi * 2;
    theta=linspace(0, np*pt,np*pr+1);
    [R,THETA]=ndgrid(linspace(0,1,nr),theta);
    x = 1-(abs(1-mod(ppr*THETA, 2*pi)/pi).^2)*.7;
    phi = (pi/2)*linspace(ol(1),ol(2),np*pr+1).^2;
    y = pf*(R.^2).*(1.27689*R-1).^2.*sin(phi);
    R2 = x.*(R.*sin(phi) + y.*cos(phi));
    
    X=R2.*sin(THETA);
    Y=R2.*cos(THETA);
    Z=x.*(R.*cos(phi)-y.*sin(phi));
    C=hypot(hypot(X,Y), Z);
end
%motion function
function [M] = motion(ppr,nr,pr,np,pf,ol)
        [X,Y,Z,C] = drawFlower(ppr,nr,pr,np,pf,ol);
        surf(X,Y,Z,C,FaceColor='interp', EdgeColor='none');
        colormap(gca,[ linspace(.6,1,256); linspace(.1,.8,256); linspace(.7,1,256); ]');
    
        %resolve axis
        axis equal      % make figure look square 
        %keep axis the same for each frame
        xlim([-1 1])    
        ylim([-1 1])
        zlim([-0.8 1])

        %save figure as a movie frame 
        M=getframe();
end 