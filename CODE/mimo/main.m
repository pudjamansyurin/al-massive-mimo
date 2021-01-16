clc
clear all
close all

for x=[30 -30]
  load(sprintf("%d/xt.mat",x));
  fc = size(xt,3);
  data = ifft(xt, fc, 3)
  range = 2*pi;
  point = range*(180/pi);
  t = linspace(0, range, point);

  if x == 30
    col = '-r*'
  else
    col = '-b*'
  endif
  
  for f=1
    AF = af(data(:,1,f), range, point);
    
    figure(1);
    plot(AF, col);
    hold on;
    
    figure(2);
    polar(t, AF, col);
    hold on;
  endfor  
endfor

for x=[3030]
  load(sprintf("%d/xt.mat",x));
  fc = size(xt,3);
  data = ifft(xt, fc, 3)
  range = 2*pi;
  point = range*(180/pi);
  t = linspace(0, range, point);

  for f=1:fc
    AF = af(data(:,1,f), range, point);
    
    figure(1);
    plot(AF);
    hold on;
    
    figure(2);
    polar(t, AF);
    hold on;
  endfor  
endfor