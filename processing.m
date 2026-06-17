%一发一收距离像显示程序
%可显示对消前后整体距离像及指定周期距离像
%可设置显示区域及线缆补偿
%李虎泉
% 设置默认字符编码为 UTF-8
feature('DefaultCharacterSet', 'UTF-8');
clear
clc
close all
%% =====================信号参数设置====================%
vc=3e8;                   %光速
startf=1.6e9;             %起始频率
stopf=2.2e9;              %截止频率
delta_f=2e6;              %频率步进
T0=100e-6;                %脉宽
NFFT=4096;                %FFT点数
Rmax=vc/(delta_f)/2;      %最大测量距离 75米
K=(stopf-startf)/delta_f; %频点数   
%% =====================显示参数设置====================%
%电缆补偿长度（电磁波传播实际距离）单位米
wire_length_1 = 0;  %一通道补偿值     
wire_length_2 = 0;  %二通道补偿值
%处理周期设置
delete_Period=10;        %删掉的周期
deal_period=300;         %实际处理的周期数
%检测区域 单位米
detectRmin = 0;         %起始距离
detectRmax = 5;         %结束距离
%显示设置
% OriginalProfile=figure;%原始距离像
% OriginalImage=figure;%原始单周期距离像
% AfterCutProfile=figure;%对消后距离像
% AfterCutImage=figure;%对消后单周期距离像
OriginalProfile=figure('Visible','off');%原始距离像
OriginalImage=figure('Visible','off');%原始单周期距离像
AfterCutProfile=figure('Visible','off');%对消后距离像
AfterCutImage=figure('Visible','off');%对消后单周期距离像
SN=80; %观察周期序号
%数据文件名
% filename='F:\行为识别数据\raw_radar_data\cc\box\'; 
filename='E:\zjx\jmh_pick_Radar1\jmh_pick_Radar1\0\'; 
% open_filename='C:\Users\Administrator\Desktop\1.8\no4\cpy\box\0\data_20190108_105415.data'; 
%% =====================数据读取与预处理====================%
filePattern=fullfile(filename,'*.data');
dataFiles=dir(filePattern);
  
for k=1:length(dataFiles)
    lastname=dataFiles(k).name;
    open_filename=fullfile(filename,dataFiles(k).name);
    fprintf('正在处理文件：%s\n',open_filename);
    fid = fopen(open_filename, 'r');
    if fid==-1
        error('无法打开文件：%s',open_filename);
    end
% 数据文件预处理 去除前后不完整周期
oridata = fread(fid,'int8');
index = find(oridata == 64);
if (index(1) == 2 && oridata(index(1)-1) == 3) 
    fid = fopen(open_filename, 'r');
    rdata = fread(fid,'int16');
else
    for i = 1:length(index)
        if(oridata(index(i)-1) == 3)
            oridata(1:index(i)-2) = [];
            break
        end
    end
    fid = fopen(open_filename, 'w');
    fwrite(fid,oridata);
    fid = fopen(open_filename, 'r');
    rdata = fread(fid,'int16');
end
index = find(rdata == 16387);
rdata(index(length(index)):end) = [];%再次找到特定值的索引，移除这些索引之后的所有数据，去除文件末尾的不完整周期
% 按通道存储
p=length(rdata);
if mod(p, 4) ~= 0
    rdata(p-p/4*4+1:end) = []; % 去除不能被4整除的部分
    p = length(rdata); % 更新rdata的长度
end
m=zeros(p/4,4);
for i=1:p/4
    m(i,1)=rdata((i-1)*4+1);          %  通道一 I
    m(i,2)=rdata((i-1)*4+2);          %  通道一 Q 
    m(i,3)=rdata((i-1)*4+3);          %  通道二 I 
    m(i,4)=rdata((i-1)*4+4);          %  通道二 Q 
end
I1=m(:,1);

Q1=m(:,2);
I2=m(:,3);
Q2=m(:,4);
I1=I1';
Q1=Q1';
I2=I2';
Q2=Q2';
%把数据按通道存储，并转换为复数形式，以便进行后续处理

fclose(fid);
clear oridata rdata p m fid i index
%% =====================数据完整性检测=====================%
index1=seekflag(I1,16387);   
NumFreqTest=index1(2:end)-index1(1:end-1);%NumFreqTest每个元素为对应周期的成功采样的频点数目
badperiod=find(NumFreqTest~=K+1);         %数据有缺失或过多的周期 (K+1:频点数加标志位)
Nbadperiod=length(badperiod);             %统计不正常周期数目,输出到命令窗口
disp(['不正常周期数目：',num2str(Nbadperiod)]);
if Nbadperiod/length(index1)>0.1          %如果不正常周期超过10/100，输出警告窗口
    msgbox('数据不正常周期超过10%！','警告','warn');
end
for k=Nbadperiod:-1:1                   %将周期采样频点数不正常的去掉
    nbadperiod=badperiod(k);
    I1(index1(nbadperiod):index1(nbadperiod)+NumFreqTest(nbadperiod)-1)=[];
    Q1(index1(nbadperiod):index1(nbadperiod)+NumFreqTest(nbadperiod)-1)=[];
    I2(index1(nbadperiod):index1(nbadperiod)+NumFreqTest(nbadperiod)-1)=[];
    Q2(index1(nbadperiod):index1(nbadperiod)+NumFreqTest(nbadperiod)-1)=[];
end
index1=seekflag(I1,16387);    %重新寻找周期索引
%实际处理的周期数,去掉最后一个周期，确保对消后还有deal_period个周期
NS=min(deal_period+1,length(index1)-delete_Period-1); 
index1=index1(delete_Period+1:delete_Period+NS);%周期索引
%需要处理的周期
Idata1=I1(index1(1):index1(end)+K);
Qdata1=Q1(index1(1):index1(end)+K);
Idata2=I2(index1(1):index1(end)+K);
Qdata2=Q2(index1(1):index1(end)+K);
%转矩阵
Idata1=reshape(Idata1,[K+1,NS]);
Qdata1=reshape(Qdata1,[K+1,NS]);
Idata2=reshape(Idata2,[K+1,NS]);
Qdata2=reshape(Qdata2,[K+1,NS]);
%======================去标志位=====================%
Idata1=Idata1(2:K+1,:);
Qdata1=Qdata1(2:K+1,:);
Idata2=Idata2(2:K+1,:);
Qdata2=Qdata2(2:K+1,:);
%======================去直流=====================%
Idata1=Idata1-repmat(mean(Idata1),K,1);
Qdata1=Qdata1-repmat(mean(Qdata1),K,1);
Idata2=Idata2-repmat(mean(Idata2),K,1);
Qdata2=Qdata2-repmat(mean(Qdata2),K,1);
%======================IQ路合成=====================%
IQresample1=(Idata1+1j*Qdata1).';
IQresample2=(Idata2+1j*Qdata2).';

clear Idata1 Qdata1 Idata2 Qdata2 index1 NumFreqTest badperiod Nbadperiod

%% =====================回波数据预处理====================%
% 1，2为加hamming窗后数据 11，22为不加窗数据
%======================加窗=====================%
IQdatawin1=IQresample1.*repmat(hamming(K).',NS,1);
IQdatawin2=IQresample2.*repmat(hamming(K).',NS,1);
IQdatawin11=IQresample1;
IQdatawin22=IQresample2;
%======================补零=====================%
IQdatawin1=[zeros(NS,startf/delta_f),IQdatawin1];
IQdatawin2=[zeros(NS,startf/delta_f),IQdatawin2];
IQdatawin11=[zeros(NS,startf/delta_f),IQdatawin11];
IQdatawin22=[zeros(NS,startf/delta_f),IQdatawin22];
%======================IFFT=====================%
[~,Nn]=size(IQdatawin1);
IQdata_FFT1=[IQdatawin1,zeros(NS,NFFT-Nn)];
IQdata_FFT2=[IQdatawin2,zeros(NS,NFFT-Nn)];

[~,Nn]=size(IQdatawin11);
IQdata_FFT11=[IQdatawin11,zeros(NS,NFFT-Nn)];
IQdata_FFT22=[IQdatawin22,zeros(NS,NFFT-Nn)];
%为IFFT变换准备信号，如果信号长度不足NFFT，会补零
clear IQdatawin1 IQdatawin2 IQdatawin11 IQdatawin22 IQresample1 IQresample2

%对每个周期的信号进行快速傅里叶变换（FFT），转换到频域
for k=1:NS
    IQdata_FFT1(k,:)=fft(IQdata_FFT1(k,:),NFFT);
    IQdata_FFT2(k,:)=fft(IQdata_FFT2(k,:),NFFT); 
    IQdata_FFT11(k,:)=fft(IQdata_FFT11(k,:),NFFT);
    IQdata_FFT22(k,:)=fft(IQdata_FFT22(k,:),NFFT);
end
rr=linspace(0,Rmax,NFFT);%显示距离


clear k
%===============电缆补偿，检测区域确定=============%
[select_data1,rlabel_small]=extractDetectZone(IQdata_FFT1,rr,wire_length_1,detectRmin,detectRmax,Rmax,NFFT); 
[select_data2,rlabel_small]=extractDetectZone(IQdata_FFT2,rr,wire_length_2,detectRmin,detectRmax,Rmax,NFFT);

[select_data11,rlabel_small]=extractDetectZone(IQdata_FFT11,rr,wire_length_1,detectRmin,detectRmax,Rmax,NFFT);
[select_data22,rlabel_small]=extractDetectZone(IQdata_FFT22,rr,wire_length_2,detectRmin,detectRmax,Rmax,NFFT);

rlabel_small_length=length(rlabel_small);


clear IQdata_FFT1 IQdata_FFT2 IQdata_FFT11 IQdata_FFT22


%% 对消前距离像
% if exist('OriginalProfile','var')
%     figure(OriginalProfile)
%     imagesc(rlabel_small,1:NS,20*log10(abs(select_data1)));
%     title('原始距离像');xlabel('目标距离(m)');ylabel('周期序号');
%     colormap('jet')
% end   
% %% 单周期对消前距离像
% if exist('OriginalImage','var')
%     figure(OriginalImage)
%     plot(rlabel_small,20*log10(abs(select_data1(SN,:))./max(abs(select_data1(SN,:)))),'r');
%     hold on
%     plot(rlabel_small,20*log10(abs(select_data11(SN,:))./max(abs(select_data11(SN,:)))),'k');
%     hold off
%     grid
%     title(['第',num2str(SN),'个周期原始距离像']);xlabel('Range(m)');ylabel('dB');
% end



%====================脉冲对消=====================%
MTI_data1=select_data1(2:NS,:)-select_data1(1:NS-1,:);
MTI_data2=select_data2(2:NS,:)-select_data2(1:NS-1,:);
MTI_data11=select_data11(2:NS,:)-select_data11(1:NS-1,:);
MTI_data22=select_data22(2:NS,:)-select_data22(1:NS-1,:);

NS=NS-1;
clear select_data1 select_data2 select_data11 select_data22

%% 对消后距离像
if exist('AfterCutProfile','var')
    % figure('Visible','off')
    % figure(AfterCutProfile);
    savemat =20*log10(abs(MTI_data1));
    max_r_db=max(savemat(:));
    min_r_db=min(savemat(:));
    disp(max_r_db);
    disp(min_r_db);
    r_dB_normalized = normalizeData(MTI_data1, max_r_db, min_r_db);
    r_dB_real=real(r_dB_normalized);
    r_dB_real(r_dB_real < 0) = 0;
    disp(r_dB_real)
    imagesc(rlabel_small,1:NS,20*log10(abs(MTI_data1)));
    xlabel('Range(m)');
    ylabel('Slow Time');
    colormap('jet')
    %%axis off;
    clim = get(gca,'CLim');
    set(gca, 'CLim', clim(2)+[-20,0]);%%添加这步后卡了阈值，但不明白什么意思
    %%savefig('E:\桌面\mat\first.fig')
    %%set(gca,'xtick',[],'ytick',[]);%去除坐标轴
    set(gca,'LooseInset', get(gca,'TightInset'))    %去除白边
    % title(lastname)
    for w=1:length(dataFiles)
        filename1 = sprintf('%03d', w); % 使用sprintf生成文件名
        savepath=['E:\zjx\jmh_pick_Radar1\jmh_pick_Radar1\0\After\',lastname(1:end-4),'jpg'];%保存地址根据文件名变化
        % savepath1=['E:\桌面\11Motion\TR\MAT\run\',lastname(1:end-4),'mat'];
        saveas(gca,savepath);%图窗保存
        E=0.7;
        V=(E-min_r_db)/(max_r_db-min_r_db);
        r_dB_real(r_dB_real < V) = V;
        % save(savepath1,'r_dB_real');
    end
    % savepath=['F:\行为识别数据\multi-view-data\radar1\fwq\1\box\135\PNG\',lastname(1:end-4),'.jpg'];
    % saveas(gca,savepath,'jpg');
end

% % 单周期对消后距离像
% if exist('AfterCutImage','var')
%     disp(size(MTI_data1));
%     disp(MTI_data1(1:5)); % 显示前5个值作为示例
%     disp(size(rlabel_small));
%     figure(AfterCutImage)
%     plot(rlabel_small,20*log10(abs(MTI_data1(SN,:))./max(abs(MTI_data1(SN,:)))),'r');
%     plot(rlabel_small,20*log10(abs(MTI_data11(SN,:))./max(abs(MTI_data11(SN,:)))),'k','LineWidth',2);
%     grid
%     title(['第',num2str(SN),'个周期对消后距离像']);xlabel('Range(m)');ylabel('dB');
% end
end

function normalizedData = normalizeData(data, maxVal, minVal)  
    normalizedData = (data - minVal) / (maxVal - minVal);  
end
