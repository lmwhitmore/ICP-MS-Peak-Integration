%% Specify parameters


path='c:\Scripts\Rob12Apr2016\';
%define output file name
output_filename=[path 'Integrals.dat'];   

time_limit=27;
time_limit_end=70;

%% Do calculations

listdir=dir([path,'*.TXT']);
fid2=fopen(output_filename,'w');

for ifile=1:length(listdir)
filename=[path, char(listdir(ifile).name)]
fid=fopen(filename,'r');
%read header lines and first record
headerline=fgetl(fid);
for i=2:7
 cline=fgetl(fid);
end;
space_position=find(isspace(headerline)==1);

nst=0;
cdata=[];


while isstr(cline)
 cdata=[cdata; str2num(cline)];
 cline=fgetl(fid);
end;

time_ind=find(cdata(:,1)>=time_limit,1);
time_ind_end=find(cdata(:,1)>=time_limit_end,1);
disp(' ');
disp(['Data are from: ' filename]);

if (~isempty(time_ind)) && (time_ind>1) && (~isempty(time_ind_end)) && (time_ind_end>time_ind)
 
 disp(['Time interval:    ' num2str(cdata(1,1)) '-' num2str(time_limit) '          '...
    num2str(time_limit) '-' num2str(time_limit_end)]);   
disp('================================================================')

%Write isotopes headers
if ifile==1
fprintf(fid2,'%s','Filename ');    
for num_param=2:length(cdata(1,:))   
    if num_param+2<=length(space_position)
     var_header_end= space_position(num_param+2)-1;
    else
     var_header_end=length(headerline);   
    end; 
    fprintf(fid2,'%s',['INT(' headerline(space_position(num_param+1)+1:var_header_end) ')  ']);
end;
fprintf(fid2,'%s\n','');
end;

fprintf(fid2,'%s',[char(listdir(ifile).name) ' ']);
for num_param=2:length(cdata(1,:))
    param_value=    interp1q(cdata(:,1),cdata(:,num_param),time_limit); 
    param_value_end=interp1q(cdata(:,1),cdata(:,num_param),time_limit_end);
    int_t0_tlim  =trapz([cdata(1:time_ind-1,1); time_limit],[cdata(1:time_ind-1,num_param); param_value]);
    int_tlim_tend=trapz([time_limit; cdata(time_ind:time_ind_end-1,1); time_limit_end],...
                        [param_value; cdata(time_ind:time_ind_end-1,num_param); param_value_end]);
    if num_param+2<=length(space_position)
     var_header_end= space_position(num_param+2)-1;
    else
     var_header_end=length(headerline);   
    end; 
    disp(['INT(' headerline(space_position(num_param+1)+1:var_header_end) ')   ' num2str(int_t0_tlim) '       ' num2str(int_tlim_tend)]);
    fprintf(fid2,'%s',[num2str(int_tlim_tend) '     ']);
end;  
fprintf(fid2,'%s\n','');

fclose(fid);
else
  disp('Check the limits of integration')
end;
end;
fclose('all');

