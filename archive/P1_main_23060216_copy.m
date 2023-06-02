clear all; close all;
dates.start= datestr(now,'mmdd,HHMMSS')

addpath(genpath('../src')); addpath(genpath('../data'));

[~,mlx_name,fnex] = fileparts(matlab.desktop.editor.getActiveFilename);
if ~isfolder('../archive')
    mkdir('../archive');
end
if strcmp(fnex,'.mlx')
    matlab.internal.liveeditor.openAndConvert([mlx_name '.mlx'],['../archive/',mlx_name,'_', datestr(now,'yymmddHH') '_copy.m'])
end

prob.randseed = rng(1,'twister');

%! choose algorithms to activate from
%! algs = {'sdpD','greedyD','gradgreedyD','crslcrD','sdpD_lmi'}
algs = {'greedyD','gradgreedyD','sdpD','crslcrD'}%

prob.algs = algs;
algs{2,1} = [];
obj     = struct(algs{:}); ctime   = struct(algs{:});
err     = struct(algs{:}); sens    = struct(algs{:});
rlx     = struct();
itrnum  = [];
algs = algs(1,:);

F_obj = @(Wo,maxrank) F_calc_det(Wo,maxrank);
% F_obj = @(Wo,maxrank) F_calc_traceinv(Wo,maxrank);
%% sensor selection
% Problem set

dtype = 'random'; %! no need to download public dataset
% dtype = 'PIV'; %! dtype = {'PIV','cylinder','SST'}; for real-world example with open data

caption = [
    'clone from github repo'
    ];
disp(caption);

pathwork = ['../work/' dtype '/']
if ~isfolder(pathwork)
    mkdir(pathwork)
end

filen = [pathwork,'diary_tmp.txt'];
diary(filen) %! empty output if .mlx
%% 
% Set matrices

h = waitbar(0,'Sensor selection...','Name',sprintf('%s.m',mlx_name));

prob = F_problem_spec(prob,dtype)
for idata = 1:prob.itr_data
    [u, loc] = F_data_load(prob,idata); 
    for j = 1:size(prob.n,2)
        waitbar((j+(idata-1)*size(prob.n,2))/size(prob.n,2)/prob.itr_data,h)
        ptmp    = prob.p(j);
        rtmp    = prob.r(j);
        ntmp    = prob.n(j);
        fprintf('Current problem: p=%04i, r=%04i, n=%06i\n',ptmp,rtmp,ntmp)
        [A,C]   = F_set_matrices(prob,rtmp,ntmp,u);
        
            [evecs,eigsA] = eig(A); eigsA = diag(eigsA);
            a = 0.99; prob.ampeig = a;
            eigsA = (abs(eigsA) > 1) .* eigsA .* 0.99 ./ abs(eigsA+(abs(eigsA)==0)) + (abs(eigsA) <= 1) .* eigsA;
            A = (a/max(abs(eigsA)))*A;

        for itime = 1:prob.itr_time
%% 
% random

            alg = 'randomD';
                tic
                for itr = 1:10
                    [sensors]=F_sensor_random( C, ptmp, 1);
                    [Wo] = F_calc_Gram_dlyap(A,C,sensors,[]);
                    objtmp(itr) = F_obj(Wo,[]);
                end
                ctime.(alg)( (idata-1)*prob.itr_time+itime, j ) = toc;
                obj.(alg)( idata, j, 1 ) = max(objtmp);
                fprintf('rank of %s ... %1.3d\n',alg,rank(Wo))
%% 
% normal greedy    

            alg = 'greedyD';
            if ismember(alg,algs)
                disp(alg)
                tic
                [sensors,pfull] = F_sensor_dg_gram(A,C,[],ptmp,F_obj);
                ctime.(alg)( (idata-1)*prob.itr_time+itime, j ) = toc;

                sens(j).(alg)(:,idata) = sensors;
                [Wo] = F_calc_Gram_dlyap(A,C,sensors,[]);
                obj.(alg)( idata, j, 1 ) = F_obj(Wo,[]);
                fprintf('rank of %s ... %1.3d\n',alg,rank(Wo))
            end
%% 
% gradient greedy

            alg = 'gradgreedyD';
            if ismember(alg,algs)
                disp(alg)
                for k = 1:size(prob.amp_reg,2)
                    tic
                    [sensors] = F_sensor_dg_gramgrad(A,C,[], ...
                        ptmp,prob.amp_reg(k));
                    ctime.(alg)( (idata-1)*prob.itr_time+itime, j, k) = toc;

                    sens(j).(alg)(:,idata,k) = sensors;
                    [Wo] =  F_calc_Gram_dlyap(A,C,sensors,[]);
                    obj.(alg)( idata, j, k ) = F_obj(Wo,[]);
                end
                fprintf('rank of %s ... %1.3d\n',alg,rank(Wo))
            end
%% 
% customized randomized subspace convex relaxation

            alg = 'crslcrD';
            if ismember(alg,algs)
                disp(alg)
                for k = 1:size(prob.dnm,2)
                    tic
                    [zhat, ~, z, ~,rslt,iter] = F_sensor_gramianCRSNC(C,...
                        ptmp, prob.maxiteration, A, prob.dnm(k));
                    ctime.(alg) ( (idata-1)*prob.itr_time+itime, j, k) = toc;
                    itrnum.(alg)( (idata-1)*prob.itr_time+itime, j, k) = iter;

                    [~,sensors] = maxk(z,ptmp);
                    sens(j).(alg)(:,idata,k) = sensors;
                    [Wo] = F_calc_Gram_dlyap(A,C,sensors,[]);
                    obj.(alg)( idata, j, k ) = F_obj(Wo,[]);
                    [Wo] = F_calc_Gram_vi_rlx(A,C,z);
                    obj.([alg '_rlx'])( idata, j, k ) = F_obj(Wo,[]);
                    %     rlx(ntmp).(alg)(:,idata) = z;
                end
                fprintf('rank of %s ... %1.3d\n',alg,rank(Wo))
            end
%% 
% sdp

            alg = 'sdpD';
            if ismember(alg,algs)
                try                    
                    cvx_where;
                catch ME
                    txterror = [mfilename ': Cannot use CVX; ',...
                        'Install from <a href="http://cvxr.com/cvx/">CVX RESEARCH</a> and get a license'];
                    error(txterror)
                end

                disp(alg)
                tic
                [z, cvx_cputime, cvx_optbnd, cvx_optval, cvx_slvitr,...
                    cvx_slvtol, cvx_status] = F_sensor_cvx_gram(A,C,ptmp);
                ctime.(alg)((idata-1)*prob.itr_time+itime, j, 1) = toc;
                echo off

                [~,sensors] = maxk(z,ptmp);
                sens(j).(alg)(:,idata) = sensors;
                [Wo] = F_calc_Gram_dlyap(A,C,sensors,[]);
                obj.(string(alg))( idata, j, 1 ) = F_obj(Wo,[]);
                [Wo] = F_calc_Gram_vi_rlx(A,C,z);
                obj.(string([alg '_rlx']))( idata, j, 1 ) = F_obj(Wo,[]);
                rlx(ntmp).(alg)(:,idata) = z;
                fprintf('rank of %s ... %1.3d\n',alg,rank(Wo))
            end

            alg = 'sdpD_lmi'; %! heavy but better approximation
            if ismember(alg,algs)
                disp(alg)
                tic
                [z, cvx_cputime, cvx_optbnd, cvx_optval, cvx_slvitr,...
                    cvx_slvtol, cvx_status] = F_sensor_cvx_gram_sdp_lmi(A,C,ptmp);
                ctime.(alg)((idata-1)*prob.itr_time+itime, j, 1) = toc;
                echo off

                [~,sensors] = maxk(z,ptmp);
                sens(j).(alg)(:,idata) = sensors;
                [Wo] = F_calc_Gram_dlyap(A,C,sensors,[]);
                obj.(string(alg))( idata, j, 1 ) = F_obj(Wo,[]);
                [Wo] = F_calc_Gram_vi_rlx(A,C,z);
                obj.(string([alg '_rlx']))( idata, j, 1 ) = F_obj(Wo,[]);
                rlx(ntmp).(alg)(:,idata) = z;
                fprintf('rank of %s ... %1.3d\n',alg,rank(Wo))
            end
        end
    end
    waitbar(idata/prob.itr_data,h)
end
close(h)
%%
dates.end = datestr(now,'mmdd,HHMMSS')
filen = [pathwork 'gramian_sensors_' prob.data_type '_' dates.end];
cname = getenv('COMPUTERNAME');
% mfn = mfilename('fullpath');
save(filen, 'algs','sens', 'obj', 'ctime','filen','err','itrnum',...'rlx',...'obj_map',...
    'prob','dates','caption','cname','mlx_name','pathwork'...
    );
who(matfile(filen))
%%
diary off
if strcmp(fnex,'.m')
    movefile([pathwork,'diary_tmp.txt'], [pathwork,'diary_',dates.start,'.txt']);
end
%%
function prob = F_problem_spec(prob,dtype)
    prob.p = [1:9,10:10:100];%[2,4,8,12,16];
        prob.pmax = max(prob.p); prob.pmin = min(prob.p);
    prob.r = 10:10:60;%2.^[4:8];%10:10:60;%[2,4,8,12,16];
        prob.rmax = max(prob.r); prob.rmin =  min(prob.r);
    prob.maxiteration = 5000;
    prob.amp_reg = 10.^[-10];%[-20 -10];
    prob.data_type = dtype;
    if strcmp(dtype,'PIV')
        prob.itr_data = 5;  prob.itr_time = 1;
        prob.cs = {'18deg'}; prob.case_num = char(prob.cs);
        prob.num_data = 1*10^4;
        prob.num_cut = prob.itr_data;
        prob.num_int = fix(prob.num_data/prob.num_cut);
        for icut = 1:prob.num_cut
            prob.usedind(icut) = mat2cell(colon(prob.num_int*(icut-1)+1,...
                prob.num_int*icut),1);
        end
        prob.pathdata = ['../data/PIV/work/' prob.case_num '/PIV_flowfield/ufield_allt.mat'];
        if ~isfile(prob.pathdata)
            txterror = [mfilename ': Cannot locate data <' dtype,'>; ',...
                'download from <a href="http://www.aero.mech.tohoku.ac.jp/rom/Nonomura2021-Airfoil-PIV-data-for-linear-ROM-18.zip">Our web page</a>, ',...
                'also see <a href="https://github.com/Aerodynamics-Lab/Airfoil-PIV-data-for-linear-ROM">Our github repo</a> for visualization'];            
            error(txterror)
        end
        [utest, ~] = F_PIV_svd(1, prob.pathdata, prob.case_num);
        prob.n = size(utest.Uorg,1);
        prob.nmax = size(utest.Uorg,1); prob.nmin = size(utest.Uorg,1);
        prob.m = size(cell2mat(prob.usedind(1)),2);
    
    elseif strcmp(dtype,'cylinder')
        prob.itr_data = 1;  prob.itr_time = 1;
        
        prob.pathdata = '../data/cylinder/DATA/FLUIDS/CYLINDER_ALL.mat';
        if ~isfile(prob.pathdata)
            txterror = [mfilename ': Cannot locate data <' dtype,'>; ',...
                'download from <a href="http://dmdbook.com/DATA.zip">dmdbook.com</a>'];
            error(txterror)
        end
        utest = matfile(prob.pathdata);
        prob.n = size(utest,'VORTALL',1);
        prob.nmax = max(prob.n); prob.nmin = min(prob.n);
    
        prob.num_data = size(utest,'VORTALL',2);
        prob.num_cut = prob.itr_data;
        prob.num_int = fix(prob.num_data/prob.num_cut);
        for icut = 1:prob.num_cut
            prob.usedind(icut) = mat2cell(colon(prob.num_int*(icut-1)+1,...
                prob.num_int*icut),1);
        end
        prob.m = size(cell2mat(prob.usedind(1)),2);
    elseif strcmp(dtype,'SST')
        prob.itr_data = 1;  prob.itr_time = 1;
        prob.pathdata = '../data/SST/';
        if ~isfile([prob.pathdata,'sst.wkmean.1990-present.nc'])
            txterror = [mfilename ': Cannot locate data <' dtype,'>; ',...
                'download [''sst.wkmean.1990-present.nc''] and [''lsmask.nc''] under data/SST/',...
                ' from <a href="https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html/">noaa.oisst.v2</a>, ',...
                ' or its high-reso version <a href="https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html/">noaa.oisst.v2 high</a>'];
            error(txterror)
        end
        [~, ~, time, mask, sst]...
            ...= F_pre_read_NOAA_SST( prob.pathdata, [prob.pathdata,'/../lsmask.nc'] );
            = F_pre_read_NOAA_SST( [prob.pathdata,'sst.wkmean.1990-present.nc'], [prob.pathdata,'lsmask.nc'] );
        save([prob.pathdata,'/sst.wkmean.1990-present.mat'],'sst','time','mask','-v7.3')

        prob.n = size(find(mask),1);
        prob.nmax = max(prob.n); prob.nmin = min(prob.n);

        prob.m = size(time,1);
    
        prob.num_data = prob.m;
        prob.num_cut = prob.itr_data;
        prob.num_int = fix(prob.num_data/prob.num_cut);
        for icut = 1:prob.num_cut
            prob.usedind(icut) = mat2cell(colon(prob.num_int*(icut-1)+1,...
                prob.num_int*icut),1);
        end
    elseif strcmp(dtype,'random')
        prob.itr_data = 1;  prob.itr_time = 1;
%         prob.n = 2.^(15:16);
        prob.n = 2.^(10:16);
        prob.nmax = max(prob.n); prob.nmin =  min(prob.n);
        prob.f = [ rand(1,prob.rmax).*10 ];
        prob.g = [-rand(1,prob.rmax).*1e-2     ];
        prob.num_data = 1*10^1;
        prob.num_cut = 1;
        prob.num_int = prob.num_data/prob.num_cut;
        for icut = 1:prob.num_cut
            prob.usedind(icut) = mat2cell(colon(prob.num_int*(icut-1)+1,...
                prob.num_int*icut),1);
        end
        prob.m = size(cell2mat(prob.usedind(1)),2);
    else
        error([mfilename ': Cannot identify ' dtype])
    end
    prob.dnm =  2^4; %! randomized subspace size for crslcr
    %_____________________________
%     prob.n = prob.n                ; prob.p = repmat(prob.pmin,1,size(prob.n,2)); prob.r = repmat(prob.rmin,1,size(prob.n,2));
%     prob.r = prob.r                ; prob.p = repmat(prob.pmin,1,size(prob.r,2)); prob.n = repmat(prob.nmin,1,size(prob.r,2));
    prob.p = prob.p                ; prob.r = repmat(prob.rmin,1,size(prob.p,2)); prob.n = repmat(prob.nmin,1,size(prob.p,2));
    %_____________________________
end
function [u, loc] = F_data_load(prob,idata)
    dtype = prob.data_type;
    if strcmp(dtype,'PIV')
        [u, loc] = F_PIV_svd(cell2mat(prob.usedind(idata)), prob.pathdata, prob.case_num);
    elseif strcmp(dtype,'cylinder')
        utest       = load(prob.pathdata);
        h           = utest.m; 
        w           = utest.n;
        [ix_tmp,iy_tmp] = meshgrid(1:w,1:h);
        Rtmp        = 25;
        cyl_mask    = ~((ix_tmp-49).^2 + (iy_tmp-99).^2 <= Rtmp^2);
        cyl_mask_vec    = reshape(cyl_mask,[],1);
        ind_nonmask = find(cyl_mask);
        N           = size(find(cyl_mask),1);
        loc.normal  = cyl_mask_vec;
        loc.mask    = cyl_mask;        
        
        u.vector    = [utest.VORTALL, utest.VORTEXTRA(:,end)];
        u.mean      = mean(u.vector,2,'omitnan');
        u.fluc      = u.vector - u.mean;
        [u.Uorg,u.Sorg,u.Vorg] = svd(u.fluc(loc.normal,:),'econ');
    elseif strcmp(dtype,'SST')
        load([prob.pathdata,'/sst.wkmean.1990-present.mat'])
        [u.Uorg, u.Sorg, u.Vorg, u.fluc, u.mean, ~] = F_pre_SVD_NOAA_SST(prob.m, time, mask, sst);
        loc.mask    = mask;
        loc.normal  = mask(:);
    elseif strcmp(dtype,'random')
        u.Sorg  = [];
        u.Uorg  = [];
        u.Vorg  = [];
        loc     = [];
    else
        error([mfilename ': Cannot identify ' dtype])
    end
end
function [A,C] = F_set_matrices(prob,rtmp,ntmp,u)
    dtype = prob.data_type;
    if strcmp(dtype,'PIV')||strcmp(dtype,'cylinder')||strcmp(dtype,'SST')
        A = (u.Sorg(1:rtmp,1:rtmp)*u.Vorg(2:prob.m,1:rtmp)')...
            *pinv(u.Sorg(1:rtmp,1:rtmp)*u.Vorg(1:prob.m-1,1:rtmp)');
        C = u.Uorg(:,1:rtmp);
    
    elseif strcmp(dtype,'random')
        assert(length(prob.f(1:rtmp/2))==length(prob.g(1:rtmp/2)))
    
        %! construct low-rank continuous-time operator (rank=k)
        k = 2*length(prob.f);  %! to perform DMD/TDMD with correct rank, set r=k
        A1 = [];
        for ii = 1:length(prob.f(1:rtmp/2))
            A2 = [[prob.g(ii) 2*pi*prob.f(ii); -2*pi*prob.f(ii) prob.g(ii)]];
            A1 = [A1 A2];
        end
        Alowrank = [];
        for ii = 1:length(prob.f(1:rtmp/2))
            Alowrank = blkdiag(Alowrank,A1(:,(ii-1)*2+1:2*ii));
        end
        dt = 0.01;
        A = expm(dt*Alowrank);
        [C,~,~] = svd(randn(ntmp,rtmp),'econ'); prob.Ctype = 'orthonormal';
    else
        error([mfilename ': Cannot identify ' dtype])
    end
end