%2D-RW run until the distance N is reached
clear; clc;
%1.***data settings***
copies=80; %numbers of repetitions of the trail
maxdis=20; %maximum among distance values to consider
% maxdis-1 is the number of distance values to consider
av=zeros(1,maxdis-1); %averages
sd=zeros(1,maxdis-1); %standard deviation (eval as sqrt of variance)
%2.***execution and data analysis***
[av,sd]=firstpassage(maxdis,copies);
%3 ***plotting***
vect=2:1:maxdis;
%---first plot---
figure(1); plot(vect,av);
% title(['Monte Carlo simulation of ', num2str(copies),' repetitions of 2D-RW stopped when a distance is reached'])
xlabel("Euclidean distance"); ylabel("Average amount of time elapsed"); grid on; hold on;
%---second plot---
figure(2); loglog(vect,av);
% title(['Monte Carlo simulation of ', num2str(copies),' repetitions of 2D-RW stopped when a distance is reached'])
xlabel("Euclidean distance"); ylabel("Average amount of time elapsed"); grid on; hold on;
%4.***functions definitions***
function [av,sd]=firstpassage(maxdis,copies)
    M=zeros(maxdis-1,copies); %matrix for averages
    av=zeros(1,maxdis-1); %averages
    delta2=zeros(1,maxdis-1); %quadratic deviations from average
    sd=zeros(1,maxdis-1); %standard deviation (eval as sqrt of variance)
    for i=1:maxdis-1
        for j=1:copies
            t=dist2DRW(i)
            M(i,j)=t
        end
    end
    %average
    for i=1:maxdis-1
        av(i)=mean(M(i,:))
    end
    %standard deviation (eval as sqrt of variance)
    for i=1:maxdis-1
        del=0;
        for j=1:copies
            delta2(j)=(av(i)-M(i,j))^2;
            del=del+delta2(j);
        end
        sd(i)=sqrt(del/copies) %sqrt of variance
    end
end
function t=dist2DRW(N) %N must be greater than 1
    x=0; y=0; t=0;
    while(sqrt(x^2+y^2)<N)
        ind=randi(4);
        switch ind
            case 1; x=x+1; case 2; x=x-1; case 3; y=y+1; case 4; y=y-1;
        end
        t=t+1;
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

