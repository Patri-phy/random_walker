%2D-SAW run until the distance N is reached
clear; clc;
%1.***data settings***
copies=80; %numbers of repetitions of the trail
maxdis=20; %maximum among distance values to consider
% maxdis-1 is the number of distance values to consider
av=zeros(1,maxdis-1); %averages
sd=zeros(1,maxdis-1); %standard deviation (eval as sqrt of variance)
%2.***execution and data analysis***
[av,sd]=firstpassage(maxdis,copies)
%3***plotting***
vect=2:1:maxdis;
%---first plot---
figure(1); plot(vect,av);
% title(['Monte Carlo simulation of ', num2str(copies),' repetitions of 2D-SAW stopped when the distance is reached'])
xlabel("Euclidean distance"); ylabel("Average amount of time elapsed"); grid on; hold on;
%---second plot---
figure(2); loglog(vect,av);
% title(['Monte Carlo simulation of ', num2str(copies),' repetitions of 2D-SAW stopped when the distance is reached'])
xlabel("Euclidean distance"); ylabel("Average amount of time elapsed"); grid on; hold on;
%4.***functions definitions***
function [av,sd]=firstpassage(maxdis,copies)
    M=zeros(maxdis-1,copies); %matrix for averages
    av=zeros(1,maxdis-1); %averages
    delta2=zeros(1,maxdis-1); %quadratic deviations from average
    sd=zeros(1,maxdis-1); %standard deviation (eval as sqrt of variance)
    for i=1:maxdis-1
        alls=fmax(i);
        for j=1:copies
            [a,b,flag1,t]=dist2DSAW(i,alls)
            while(flag1==0) %enter this cycle if the preceding SAW trail is aborted
                [a,b,flag1,t]=dist2DSAW(i,alls);
            end
            M(i,j)=t
        end
    end
    %average
    for i=1:maxdis-1
        av(i)=mean(M(i,:))
    end
    %variance
    for i=1:maxdis-1
        del=0;
        for j=1:copies
            delta2(j)=(av(i)-M(i,j))^2;
            del=del+delta2(j)
        end
        sd(i)=sqrt(del/copies)
    end
end
% The number of allowed sites (alls) for the SAWr the with distance from origin less than N
function alls=fmax(N)
    alls=0; o=N+1
    for i=1:2*N+2
        for j=1:2*N+2
            %circle of radius N centered in (N+1,N+1)
            if sqrt((i-o)^2+(j-o)^2) < N
                alls=alls+1
            end 
        end
    end
end
function [a,b,flag1,t]=dist2DSAW(N,alls) %N must be greater than 1
    a=zeros(1,alls); b=zeros(1,alls); %spacial coordinates time evolution
    x=0; y=0; t=0; mem=zeros(2*N+2,2*N+2); o=N+1; mem(o,o)=1;
    while(sqrt(x^2+y^2)<N)
        flag0=1; flag1=1;
        while flag0==1
            r=mem(o+x+1,o+y); l=mem(o+x-1,o+y); u=mem(o+x,o+y+1); d=mem(o+x,o+y-1);
            if (r==1 & l==1 & u==1 & d==1) %verifies if the SAWr cannot move anymore
                flag1=0; break; % break implies that the statements to execute
                % are those that follow the end of the while/for cycle (in this case the while)
            end
            ind=randi(4); x1=x; y1=y;
            switch ind
                case 1; x1=x1+1; case 2; x1=x1-1; case 3; y1=y1+1; case 4; y1=y1-1;
            end
            flag0=mem(o+x1,o+y1); %now establish if in this while the flag==1 is fulfilled 
        end
        if flag1 == 0
            break %exit from while cycle
        end
        x=x1; y=y1; t=t+1; a(t)=x; b(t)=y;
        mem(o+x,o+y)=1
    end
    %plot(a,b)
    % figure;
    % xlim([-10,10]); ylim([-10,10]);
    % hold on;
    % pause
    % for i=1:t
    %     line(a(i:i+1),b(i:i+1))
    %     pause
    % end
end

