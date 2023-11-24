%Aim: analysis of the statistic of the records of a 1-D random walk race

%terms explanations
%rw(r) means random walk(er)
%trial refers to the execution of a single rw
%rw time: time to run a single movement in a rw

%normalizations
%rw time unit and space unit are both normalized to 1; the origin of position is 0.
%trial time unit is normalized to 1;

%variable meaning
% d     distance that, once reached, determines the end of the rw trial
% p     position of the rwr
% t     amount of time elapsed in a trial
% ib    index of beated record (not to update)
% inbt  index of nonbeated trial (updated for every new record)
% ammt  ammount of time between two subsequent records
% ttv   trial time (between two subsequent records) vector
% trv   rw time (between two subsequent records) vector
% l     number of records

%algorithm of [ttv,trv,l]=record_rw_1D(d)
%1. a randow walk is run until the rwr reaches the distance of d (input)
%2. the amount of time elapsed in the trial is calculated through the function
%   t=trw_1D(d); if it is different from d then the record can be beated; otherwise it cannot and
%   the program is ended (the subsequent while cycle cannot be entered)
%3. a sequence of random walks is run until the randow walk with the minimum time value has achieved;
%   at the end of each random walk one of the two following conditions takes place: a or b.
%       a. The current random walk beats the current record;
%       b. The current randow walk does not beat the current record.

clc; clear;
d=15; copies=60;
lim=1000; %upper limit of the value t-d
sd=num2str(d); sc=num2str(copies);
ttv=zeros(1,d); trv=zeros(1,d); l=0;
avt=zeros(1,lim); avr=zeros(1,lim);
% ---single process simulation ---
[ttv,trv,l]=record_rw_1D(d)
% adaptation for the plots; l is the number of records
vect=1:l
ttv0=zeros(1,l); trv0=zeros(1,l);
ttv0(1:l)=ttv(1:l); trv0(1:l)=trv(1:l);
% plotting
figure(1);
% plot(vect,ttv0);
semilogy(vect,ttv0);
title(['1.Updates history of time records (1-D rwrs trials within the distance ', sd,' )'])
xlabel("Index of record"); ylabel("Number of trials between this record and the preceding"); hold on;
figure(2);
% plot(vect,trv0);
semilogy(vect,trv0)
title(['2.Updates history of time records (1-D rwrs trials within the distance ', sd,' )'])
xlabel("Index of record"); ylabel("Time between this record and the preceding"); hold on;
% --- Monte Carlo simulation ---
[avt,avr,lm]=record_rw_1D_mc(d,copies,lim)
lmvect=1:lm;
avt=avt(lmvect); avr=avr(lmvect);
% zeros(1,lm); avr0=zeros(1,lm);
% avt0(1:lm)=avt(1:lm); avr0(1:lm)=avr(1:lm);
figure(3)
semilogy(lmvect,avt);
title(['3.MC simulation of ', sc,' repetitions of a 1-D rwrs race within the distance ', sd])
xlabel("Index of record");
ylabel("Average over the n. of trials betw. this record and the preceding"); hold on;
figure(4)
semilogy(lmvect,avr);
title(['4.MC simulation of ', sc,' repetitions of a 1-D rwrs race within the distance ', sd])
xlabel("Index of record")
ylabel("Average over the t. coor.tes between this record and the preceding")
function [avt,avr,lm]=record_rw_1D_mc(d,copies,lim)
    Mt=zeros(lim, copies); Mr=zeros(lim, copies); %matrices for ttv and trv copies storage
    avt=zeros(1,lim); avr=zeros(1,lim); %vectors for averages ttv and trv storage
    ll=zeros(copies); %vector for the length of ttv and trv copies
    % 1.data storage in Mt and Mr
    for i=1:copies
        [ttv,trv,l]=record_rw_1D(d)
        while(l==[1,1] | lim < l)
            [ttv,trv,l]=record_rw_1D(d)
        end
        %this race is nontrivial and the data can be stored
        for j=1:l
            Mt(j,i)=ttv(j); Mr(j,i)=trv(j);
        end
        ll(i)=l
    end
    lm=max(ll);
    % 2.averages calculations
    for i=1:lm
        c1=0; c2=0; %numbers of data for ttv and trv
        for j=1:copies
            if(Mt(i,j)~=0)
                c1=c1+1
            end
            if(Mr(i,j)~=0)
                c2=c2+1
            end
        end
        % averages of trial time values
        if (c1~=0)
            avt(i)=sum(Mt(i,:))/c1;
        else
            avt(i)=0;
        end
        % averages of random walk time values
        if (c2~=0)
            avr(i)=sum(Mr(i,:))/c2;
        else
            avr(i)=0;
        end
    end
end
%single trial generation; evaluation of the amount of time for the trial
function t=trw_1D(d)
    p=0; t=0;
    while(p~=d & p~=-d)
    %while(-d<p&p<d)
        f=rand;
        if (f<0.5)
            p=p+1;
        else
            p=p-1;
        end
        t=t+1;
    end
end
%function for records storage (ttv and trv) and numbers of records evalutation (l)
%t is the amount of time for the first trail
%ttv and trv are created with d-t column because the maximum number of records that the
%rwr can beat after a trial of t unit of time is less than or equal to d-t
function [ttv,trv,l]=record_rw_1D(d)
    ib=1; inbt=0; ammt=0; l=0;
    t=trw_1D(d); ttv=zeros(1,d-t); trv=zeros(1,d-t);
    cur_r=t; 
    % i=0
    while(cur_r~=d)
        % i=i+1;
        t=trw_1D(d);
        if (t<cur_r)
            cur_r=t;
            ttv(1,ib)=inbt; trv(1,ib)=ammt;
            ib=ib+1; inbt=0; ammt=0;
        else
            inbt=inbt+1;
            ammt=ammt+t;
        end
    end
    l=ib-1;
end