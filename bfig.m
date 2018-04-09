% Generate a figure on the Desktop
%
function bfig(ext,dirname,fname)

if nargin<1, ext = 'png'; end

set(gcf,'PaperPositionMode','Auto')

c = computer;
ismac = strncmp(c,'MAC',3);

if ~exist('dirname','var'), dirname = 'figs'; end

if strcmp(c,'PCWIN')
    dirname = [getenv('USERPROFILE'),'/Desktop/',dirname,'/'];
elseif (ismac) || strcmp(c,'GLNX86') || strcmp(c,'GLNA64')
    dirname = [getenv('HOME'),'/Desktop/',dirname,'/'];
else
    fprintf('Platform not recognized.\n'); return;    
end
if ~exist(dirname,'dir')
    mkdir(dirname);
end

if ~exist('fname','var')
	tmp = dir([dirname,'fig*.',ext]);
	tmp = {tmp.name}';

	if isempty(tmp)
		fname = 'fig_1';
	else
		num = regexp(tmp,['(?<=fig_)\d+(?=(.',ext,'))'],'match');
		num = str2double([num{:}]);
		num = max(num)+1;
		fname = ['fig_',num2str(num)];
	end
end

% fprintf('dirname=%s .. fname=%s\n',dirname,fname);
switch ext
    case 'eps'
        print('-depsc','-noui',[dirname,fname])
    case 'png'
        print('-dpng','-r0','-noui',[dirname,fname])
    case 'pdf'
        if (ismac)
            print('-depsc','-noui',[dirname,fname])
            eval(['!pstopdf ',dirname,fname,'.eps'])
        else
            fprintf('I only know how to make pdf in Mac OS\n')
            return;
        end
    otherwise
        fprintf('Not understood!\n')
end
