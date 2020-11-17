clear all; clc;
%==============================Modulasi 4 QAM==============================
hasil=[];
simbol=[-1 1; 1 1; 1 -1 ;-1 -1]/sqrt(2); 
 
%Konversi biner integer ke desimal integer
jBit = 4*10;
dataSource = randint(1,jBit);
biner=dataSource;
b=[];
n=2; %Jumlah bit biner
binerString=int2str(biner);
for i=1:3:length(binerString)
    a=binerString(1,i);
    b=[b,a];
end
b; %Biner dalam format string
desimal=[];
e=[];
for t=1:n:length(b),
    d=[];
    for k=1:n;
        c=b(1,t+(k-1));
        d=[d c];
    end
    e=[e; d];
end
e; % 2 Biner matriks kolom
for no=1:size(e,1),
    des=bin2dec(e(no,:));
    desimal=[desimal des];
end
desimal; % Hasil konversi ke desimal
dataDesimal=desimal;
kompleks=simbol(dataDesimal+1,:); %Bentuk simbol dari data desimal
aData=[];
for n = 1:size(kompleks,1),
    hasil=complex(kompleks(n,1),kompleks(n,2));
    aData=[aData hasil]; % Reshape matriks hasil ke kolom
end
%Hasil modulasi 4 QAM
aData

