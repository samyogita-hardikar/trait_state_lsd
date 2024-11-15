
%% loading data, surface, parcel indeces

lh_pial = SurfStatReadSurf1('schaefer_atlas/parcellations/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/surf/lh.pial');
rh_pial = SurfStatReadSurf1('schaefer_atlas/parcellations/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/surf/rh.pial');

SC = combine_surfaces(lh_pial, rh_pial);

parcel_index=importdata('schaefer_atlas/fsaverage5_400parcels.txt');

%multivariate gradient matrix

n_grads={1,2,3,4,5};


for n = 1: length(n_grads)
    grad = importdata(sprintf('/all_lsd/gradients_400_parcels/surfstat/groupmat_parcels_grad%d.txt', n));
    grad_mult(:,:,n)= grad;
end
  


    
var_mat= readtable('/all_lsd/gradients_400_parcels/LSD_beh_vars.csv');
var_mat= var_mat(:,[1 12:36]);


trait1=term(var_mat.FAC1_2);
trait2=term(var_mat.FAC2_2);
trait3=term(var_mat.FAC3_2);
trait4=term(var_mat.FAC4_2);
trait5=term(var_mat.FAC5_2);

thought1=term(var_mat.Thought_C1);
thought2=term(var_mat.Thought_C2);
thought3=term(var_mat.Thought_C3);
thought4=term(var_mat.Thought_C4);
thought5=term(var_mat.Thought_C5);

gender=term(var_mat.gender);
age=term(var_mat.age_years);
motion=term(var_mat.motion);


savepath= '/data/p_02260/all_lsd/results/'


%%



M=1+ trait1+ trait2+ trait3+ trait4+ trait5+ thought1+ thought2+ thought3+ thought4+ thought5+ age+ motion+ gender;
%M0_trait1=1+  trait2+ trait3+ trait4+trait5+ thought1+ thought2+ thought3+ thought4+ thought5+ age+ motion+ gender;

var_list= {1 trait1 trait2 trait3 trait4 trait5 thought1 thought2 thought3 thought4 thought5 age motion gender};
var_names= {'intercept' 'trait1' 'trait2' 'trait3' 'trait4' 'trait5' 'thought1' 'thought2' 'thought3' 'thought4' 'thought5' 'age' 'motion' 'gender'};

% M=1+ thought1+ thought2+ thought3+ thought4+ thought5+ age+ motion+ gender;
% var_list= {1 thought1 thought2 thought3 thought4 thought5 age motion gender};
% var_names= {'intercept' 'thought1' 'thought2' 'thought3' 'thought4' 'thought5' 'age' 'motion' 'gender'};

%slm1=SurfStatLinMod(grad_mult(:,:,1:3),M);
%slm0 = SurfStatLinMod(grad_mult(:,:,1:3), M0_trait4);
%slm = SurfStatF( slm1, slm0 ); 

%SurfStatView( slm.t', SC, 'F statistic, 3,143 df' );

slm=SurfStatLinMod(grad_mult(:,:,1:3),M);


for v = 7: 9
    
    slm =SurfStatT(slm, var_list{v});
    
    FDR_mask=zeros(1, 400);
    for k = 1:length(SurfStatQ(slm).Q)
        if SurfStatQ(slm).Q(:,k) < 1
            FDR_mask(:,k)=1;
        end
    end 
    bonf = sum(FDR_mask==1);
    if any(FDR_mask(:)==1)
        for g = 1:3
            uni_slm=SurfStatLinMod(grad_mult(:,:,g),M);
            uni_slm =SurfStatT(uni_slm, var_list{v});
            uniBonf_mask=zeros(1, 400);
            pp = 1 - tcdf(uni_slm.t,uni_slm.df);
            pn = 1 - tcdf(-uni_slm.t,uni_slm.df);
            
            p= zeros(size(pp));
            p = pp<pn;
            p_all= zeros(size(p));
            p_all(p==1) = pp(p==1);
            p_all(p==0) = pn(p==0);

            for j = 1:length(pp)
                if pp(:,j) < 1
                    uniBonf_mask(:,j)=1;
                elseif pn(:,j) < 1
                    uniBonf_mask(:,j)=1;
                end
                                    
                if pp(:,j) < 1
                    uniUncorr_mask(:,j)=1;
                    %p(:,j)=pp(:,j)
                elseif pn(:,j) < 1
                    uniUncorr_mask(:,j)=1;
                    %p(:,j)=pn(:,j)
                end
            end
            p=[pp ; pn];
             dlmwrite([savepath, 'p_Uncorr_univar_' var_names{v} '_grad' num2str(g) '.txt'], p_all(:));
%             dlmwrite([savepath, 'p_Bonf_univar_' var_names{v} '_grad' num2str(g) '.txt'], p_all(:).*bonf);
             dlmwrite([savepath, 't_univar_' var_names{v} '_grad' num2str(g) '.txt'], uni_slm.t');
%             dlmwrite([savepath, 'coef_univar_' var_names{v} '_grad' num2str(g) '.txt'], uni_slm.coef');
%             dlmwrite([savepath, 'ef_univar_' var_names{v} '_grad' num2str(g) '.txt'], uni_slm.ef');
%             dlmwrite([savepath, 'sd_univar_' var_names{v} '_grad' num2str(g) '.txt'], uni_slm.sd');
%  


                    
            var2plot= uni_slm.coef(v,:)'.*uniBonf_mask(:);
            var2plot=var2plot.*FDR_mask(:);
            COnSurf = parcel2full(var2plot, parcel_index);
            f = figure;
            Rs = COnSurf';
            %SurfStatView(Rs(10243:20484),rh_pial, ['sept_Ef ' var_names{v} ' G' num2str(g)  ' RH (p<' num2str(0.025/bonf) ')']);
            SurfStatViewData(Rs,SC, ['sept_coef ' var_names{v} ' G' num2str(g)  ' p<' num2str(0.025/bonf) ]);
            colormap(redblue);
            SurfStatColLim([-0.01 0.01]);
 %           exportfigbo(f,[savepath, 'sept_Reg_Coef_Bonf_RH' var_names{v} '_grad' num2str(g) '.png'],'png', 9);
%             close(f)
            
        end
        
        %if v>1
%         parcels2plot= find(FDR_mask==1);
%         for p=1:length(parcels2plot)
%             for g= 1:3
%                 f=figure;
%                 SurfStatPlot( var_list{v}, grad_mult(:,parcels2plot(p),g), M);
%                 title(var_names{v});
%                 %scatter(var_list{v}, grad_mult(:,parcels2plot(p),g));
%                 ylabel(['G' num2str(g) '-' num2str(parcels2plot(p)) '- corrected']);
%                 saveas(f, [savepath, 'scatter_plots/T-stat_' var_names{v} '_grad' num2str(g) '-' num2str(parcels2plot(p)) '.png'],'png');
%                 %close(f)
%             end
%         end
        
    
    %end
        
    end
%     dlmwrite([savepath, 'pFDR_' var_names{v} '.txt'], SurfStatQ(slm).Q');
end
    
