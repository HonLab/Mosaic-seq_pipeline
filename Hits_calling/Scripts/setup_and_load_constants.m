mkdir('pvals');
mkdir('plots');
mkdir('cpm');

regions = {'CD52.SE3.HS2',...
             'CD52.SE3.HS4',...
             'CD52.SE3.HS5',...
             'CD52.SE3.HS6',...
             'Control',...
             'PIM1.SE2.HS2',...
             'PIM1.promoter',...
             'PIM1.SE1.HS1',...
             'PIM1.SE1.HS10',...
             'PIM1.SE1.HS2',...
             'PIM1.SE1.HS3',...
             'PIM1.SE1.HS6',...
             'PIM1.SE1.HS7',...
             'PIM1.SE1.HS8',...
             'PIM1.SE1.HS9',...
             'PIM1.SE2.HS1',...
             'PIM1.SE2.HS10',...
             'PIM1.SE2.HS11',...
             'PIM1.SE2.HS12',...
             'PIM1.SE2.HS3',...
             'PIM1.SE2.HS4',...
             'PIM1.SE2.HS5',...
             'PIM1.SE2.HS6',...
             'PIM1.SE2.HS9',...
             'beta_globin.HBE1_positive_control',...
             'beta_globin.HBE1.promoter',...
             'beta_globin.HS2_positive_control',...
             'beta_globin.SE1.HS1',...
             'beta_globin.SE1.HS2',...
             'beta_globin.SE1.HS3',...
             'beta_globin.SE1.HS4',...
             'CD52.promoter',...
             'CD52.SE1.HS1',...
             'CD52.SE1.HS2',...
             'CD52.SE1.HS3',...
             'CD52.SE2.HS2',...
             'CD52.SE2.HS3',...
             'CD52.SE2.HS4',...
             'CD52.SE2.HS5',...
             'CD52.SE2.HS6',...
             'CD52.SE3.HS1',...
             'GFP',...
             'FADS1.promoter',...
             'FADS1.SE1.HS1',...
             'FADS1.SE2.HS1',...
             'FADS1.SE2.HS2',...
             'FADS1.SE3.HS3',...
             'FADS1.SE3.HS4',...
             'FADS1.SE4.HS1',...
             'FADS1.SE4.HS2',...
             'FADS1.SE4.HS3',...
             'PRKAR2B.promoter',...
             'PRKAR2B.SE1.HS1',...
             'PRKAR2B.SE1.HS2',...
             'PRKAR2B.SE1.HS3',...
             'PRKAR2B.SE2.HS1',...
             'PRKAR2B.SE2.HS2',...
             'PRKAR2B.SE3.HS1',...
             'PRKAR2B.SE3.HS2',...
             'PRKAR2B.SE3.HS4',...
             'PRKAR2B.SE3.HS5',...
             'PRKAR2B.SE3.HS8',...
             'SMYD3.promoter',...
             'SMYD3.SE2.HS1',...
             'SMYD3.SE2.HS10',...
             'SMYD3.SE2.HS11',...
             'SMYD3.SE2.HS12',...
             'SMYD3.SE2.HS2',...
             'SMYD3.SE2.HS3',...
             'SMYD3.SE2.HS4',...
             'SMYD3.SE2.HS6',...
             'SMYD3.SE2.HS7',...
             'SMYD3.SE2.HS8',...
             'SMYD3.SE3.HS1',...
             'SMYD3.SE3.HS3',...
             'SMYD3.SE3.HS4',...
             'SMYD3.SE3.HS5',...
             'SMYD3.SE3.HS6',...
             'SMYD3.SE3.HS7',...
             'UFD1L.promoter',...
             'UFD1L.SE1.HS1',...
             'UFD1L.SE1.HS4'};
