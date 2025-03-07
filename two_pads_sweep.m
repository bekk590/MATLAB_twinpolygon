function [Freqs, Q ,m_eff, S_F, eta, rl2_match, Q_match] = ...
         two_pads_sweep(stress, h_mbr, l0, w0, N,...
                           N_sol,rl1, rl2, rw1, rw2, ...
                           l_trans, l_pad, w_pad, values,...
                           plot_flag, plot_op_flag, pad_trigger)
    % Description
    % COMSOL MATLAB Livelink function for simulating binary-tree mechanical
    % resonators. This function was originally made for simulating the
    % engineered binary tree beams for the integrated samples project
    % Inputs:
    % stress: Depostion stress
    % h_mbr: Thickness of the membrane
    % l0: Half length of the central branch
    % w0: Width of the central branch
    % N: Number of branching generations
    % theta: Branching anlge (radians)
    % rl: Branches length contraction ratio
    % rw: Branches width contraction ratio
    % l_pad: Length of the interaction pad in the center
    % w_pad: Width of the interaction pad in the center
    % N_sol: Number of desired eigenmodes
    % plot_flag: Boolean - 1 for plotting the geometry and the mode shapes
    % Outputs:
    % Freqs: Frequency  of the eigenmodes (in units of Hz)
    % Q: Quality factor of the eigenmodes (calculated based on the
    %    dissipation-dilution formula)
    % m_eff: Effective mass of the eigenmode
    % S_F: Thermal force power spectral density at room temperature =
    %      S_F = 4 * K_B * T* m_eff * 2*pi*f_m / Q
    % eta: Out-of-plane/in-plane indicator - 1 for OP, 0 for IP
    %%                   
    import com.comsol.model.*;
    import com.comsol.model.util.*;

    % Loading the previously created model
    % model = mphload('SiN_PhC_cavity.mph');

    model = ModelUtil.create('Model');
    ModelUtil.showProgress(1);
    %% defining the global parameters:
    Qint = 12500 * (h_mbr/100e-9);
    f0 = 100e3;
    %l_trans = 5e-6;
    %N_seg = 2^(N+1) - 1 + 2^N; % number of the segments 
    rotangle = (2*pi/N);
    %rotvalue = (2*pi/N);
    Diameter = sqrt((l0*rl1)^2+(l0*rl2)^2);
    rad = Diameter/2;
    theta = (pi-rotangle)/2;
    Radx = rad*sin(atan(max(rl1/rl2,rl2/rl1))-pi/4);
    Rady = rad*cos(atan(max(rl1/rl2,rl2/rl1))-pi/4);
    segments = cell(3,1);
    segmentsString = cell(3,1);
    x_seg = cell(3,1);
    y_seg = cell(3,1);
    theta_seg = cell(3,1);
    l_seg = cell(3,1);
    w_seg = cell(3,1);
    theta_seg{1} = 'pi/2';
    disp(theta_seg{1});

   
    model.param.set('stress', sprintf('%d[Pa]',stress), 'Initial stress');
    model.param.set('h_mbr', sprintf('%d[m]', h_mbr), 'Thickness of the resonator');
    model.param.set('l0', sprintf('%d[m]', l0), 'Length of the fundamental segment');
    model.param.set('w0', sprintf('%d[m]', w0), 'Width of the fundamental segment');
    model.param.set('N', sprintf('%i', N), 'Number of the tops');
    model.param.set('rl1', sprintf('%d',rl1), 'Length contraction ratio left');
    model.param.set('rw1', sprintf('%d',rw1), 'Width expansion ratio left');
    model.param.set('rl2', sprintf('%d',rl2), 'Length contraction ratio right');
    model.param.set('rw2', sprintf('%d',rw2), 'Width expansion ratio right');
    model.param.set('l_pad', sprintf('%d[m]', l_pad), 'Length of the pad');
    model.param.set('w_pad', sprintf('%d[m]', w_pad), 'Width of the pad');
    model.param.set('l_trans', sprintf('%d[m]', l_trans), 'Length of the spline for the pad');
    model.param.set('Qint', '12500 * (h_mbr/100e-9[m])');
    model.param.descr('Qint', 'Intrinsic Q');    
    model.param.set('rotangle', '(2*pi/N)');
    model.param.descr('rotangle', 'Angle of rotation');
    model.param.set('theta', '((pi-rotangle)/2)');
    model.param.descr('theta', 'Branching angle (radians)');
    model.param.set('Diameter', 'sqrt((l0*rl1)^2+(l0*rl2)^2)');
    model.param.descr('Diameter', 'Diameter of circumscribed circle');
    model.param.set('rad', 'Diameter/2');
    model.param.descr('rad', 'Radius of circumscribed circle');
    model.param.set('Radx', 'rad*sin(atan(max(rl1/rl2,rl2/rl1))-pi/4)');
    model.param.descr('Radx', 'x-coordinate of polygon vertex');
    model.param.set('Rady', 'rad*cos(atan(max(rl1/rl2,rl2/rl1))-pi/4)');
    model.param.descr('Rady', 'y-coordinate of polygon vertex');

    
    
    %disp(rotvale);
    %rw = 1/rW1;

    l_seg{1} = 'l0';
    w_seg{1} = 'w0';

    if rl1>=rl2
        x_seg{1} = 'Radx';
    else
        x_seg{1} = '-Radx';
    end
    y_seg{1} = sprintf('Rady + (%s/2)',l_seg{1});
    
    l_clamp = '50 * h_mbr';
    n_clamp_x = 10;
    n_clamp_y = 30;
    mesh_max = 1e-6;
    mesh_min = mesh_max/2;
    
    segmentsString{1} = 'r1';
    
    clamp_lines = cell(2^N,1);
    clamp_linesString = cell(2^N,1);

    for i = 2:3
        ind = i;
        ind_p = 1;
        segmentsString{ind} = sprintf('r%i',ind);
        disp(theta_seg{ind_p});
        if i == 2
            theta_seg{ind} = sprintf('(%s - theta)', theta_seg{ind_p});
            l_seg{ind} = 'l0*rl1/2';
            w_seg{ind} = 'w0*rw1';
            dl_rw = sprintf('0.5 * %s * (cos(%s) - %s) / sin(%s)', w_seg{ind_p}, 'theta', 'rw1', 'theta');
            if dl_rw<0
                dl_rw = 0;
            end
        disp(theta_seg{ind});
        end
        if i == 3
            theta_seg{ind} = sprintf('(%s + theta)', theta_seg{ind_p});
            l_seg{ind} = 'l0*rl2/2';
            w_seg{ind} = 'w0*rw2';
            dl_rw = sprintf('0.5 * %s * (cos(%s) - %s) / sin(%s)', w_seg{ind_p},'theta','rw2','theta');
            if dl_rw<0
                dl_rw = 0;
            end
        end

%        if i == N
%            dl_seg = 0.5*w_seg(ind) * min(abs(cot(theta_seg(ind))),...
%                                          abs(tan(theta_seg(ind))));
%            l_seg(ind) = l_seg(ind) + dl_seg;
%        end
%        dl_rw = 0.5 * w_seg(ind_p) * (cos(theta1) - rw) / sin(theta1);
%        if dl_rw<0
%            dl_rw = 0;
%        end
%        x_seg(ind) = x_seg(ind_p) - ...
%                     (l_seg(ind_p)/2-dl_rw) * cos(theta_seg(ind_p)) - ...
%                     l_seg(ind)/2 * cos(theta_seg(ind));
%        y_seg(ind) = y_seg(ind_p) - ...
%                     (l_seg(ind_p)/2-dl_rw) * sin(theta_seg(ind_p)) - ...
%                     l_seg(ind)/2 * sin(theta_seg(ind));


        x_seg{ind} = sprintf('%s - (%s/2-%s) * cos(%s) - %s/2 * cos(%s)', x_seg{ind_p}, l_seg{ind_p}, dl_rw, theta_seg{ind_p}, l_seg{ind},  theta_seg{ind} );
        y_seg{ind} = sprintf('%s - (%s/2-%s) * sin(%s) - %s/2 * sin(%s)', y_seg{ind_p}, l_seg{ind_p}, dl_rw, theta_seg{ind_p}, l_seg{ind},  theta_seg{ind} );
        disp(x_seg{ind});
    end

%    model.param.set('Qint', sprintf('%d',Qint), 'Intrinsic Q');
%    model.param.set('stress', sprintf('%d[Pa]',stress), 'Initial stress');
%    model.param.set('h_mbr', sprintf('%d[m]', h_mbr), 'Thickness of the resonator');
%    model.param.set('l0', sprintf('%d[m]', l0), 'Length of the fundamental segment');
%    model.param.set('w0', sprintf('%d[m]', w0), 'Width of the fundamental segment');
%    model.param.set('N', sprintf('%i', N), 'Number of the tops');
%    model.param.set('theta', sprintf('%d', theta), 'Branching angle (radians)');
%    model.param.set('rl1', sprintf('%d',rl1), 'Length contraction ratio left');
%    model.param.set('rw1', sprintf('%d',rw1), 'Width expansion ratio left');
%    model.param.set('rl2', sprintf('%d',rl2), 'Length contraction ratio right');
%    model.param.set('rw2', sprintf('%d',rw2), 'Width expansion ratio right');
%    model.param.set('l_pad', sprintf('%d[m]', l_pad), 'Length of the pad');
%    model.param.set('w_pad', sprintf('%d[m]', w_pad), 'Width of the pad');
%    model.param.set('l_trans', sprintf('%d[m]', l_trans), 'Length of the spline for the pad');



    %% Geometry
    geom1 = model.geom.create('geom1', 3);
    wp1 = geom1.feature.create('wp1', 'WorkPlane');
    wp1.set('quickplane', 'xy');

    disp(x_seg{1});
    disp(l_seg{1});
    
    for i = 1:3
        segments{i} = wp1.geom.feature.create(segmentsString{i}, 'Rectangle');
        segments{i}.set('size', {l_seg{i} w_seg{i}});
        segments{i}.set('base', 'center');
        segments{i}.set('rot', sprintf('%s*180/pi', theta_seg{i}));
        segments{i}.set('pos', {x_seg{i} y_seg{i}});


        
        

        if i == 1
            l_cut =  sprintf('%s / max(abs(cos(%s)), abs(sin(%s)))',w_seg{i}, theta_seg{i},theta_seg{i});
            theta_mod = sprintf('mod(%s, 2*pi)', theta_seg{i});
            dl_seg = sprintf('0.5*%s * min(abs(cot(%s)), abs(tan(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
            dx_cut = sprintf('%s * fix(sqrt(2)*cos(%s))/2', l_cut, theta_mod);
            dy_cut = sprintf('%s * fix(sqrt(2)*sin(%s))/2', l_cut, theta_mod);

            if mod(eval(theta_mod),pi/4) == 0 && mod(eval(theta_mod), pi/2) ~= 0
                if cos(theta_mod) < 1e-10
                    dx_cut = '0';
                else
                    dx_cut = sprintf('sign(cos(%s)) * %s/2', theta_mod, l_cut);
                end
                dy_cut = '0';
            end
            x_cut = sprintf('%s + (%s - 2 * %s) * cos(%s) / 2 + %s', ...
                 x_seg{i}, l_seg{i}, dl_seg, theta_seg{i}, dx_cut);
            y_cut = sprintf('%s + (%s - 2 * %s) * sin(%s) / 2 + %s', ...
                 y_seg{i}, l_seg{i}, dl_seg, theta_seg{i}, dy_cut);
            
            i_cut = i+3;
            segmentsString{i_cut} = sprintf('r%i',i_cut);
            segments{i_cut} = wp1.geom.feature.create(segmentsString{i_cut}, 'Rectangle');
            segments{i_cut}.set('size', {l_cut l_cut});
            segments{i_cut}.set('base', 'center');
            segments{i_cut}.set('pos', {x_cut y_cut});
            
            i_cl = 1;
            if i_cl>0
                x_cl = sprintf('%s + (%s - 2 * %s - %s) * cos(%s) / 2', ...
                    x_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i});
                y_cl = sprintf('%s + (%s - 2 * %s - %s) * sin(%s) / 2', ...
                    y_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i});
                l_cl = l_cut;
                theta_cl = sprintf('fix(sqrt(2)*cos(%s)) * pi / 2', theta_mod);
                if mod(eval(theta_mod),pi/4) < 1e-10
                    theta_cl = 'pi/2';
                end
                
                x1_cl = sprintf('%s - %s * sin(%s) / 2', x_cl, l_cl, theta_cl);
                x2_cl = sprintf('%s + %s * sin(%s) / 2', x_cl, l_cl, theta_cl);
                y1_cl = sprintf('%s - %s * cos(%s) / 2', y_cl, l_cl, theta_cl);
                y2_cl = sprintf('%s + %s * cos(%s) / 2', y_cl, l_cl, theta_cl);
                clamp_linesString{i_cl} = sprintf('ls%i',i_cl);
                clamp_lines{i_cl} = wp1.geom.create(clamp_linesString{i_cl}, 'LineSegment');
                clamp_lines{i_cl}.set('specify1', 'coord');
                clamp_lines{i_cl}.set('coord1', {x1_cl y1_cl});
                clamp_lines{i_cl}.set('specify2', 'coord');
                clamp_lines{i_cl}.set('coord2', {x2_cl y2_cl});

                x_bl = sprintf('%s', x_seg{i});
                y_bl = sprintf('%s-%s/2', y_seg{i}, l_seg{i});
                x1_bl = sprintf('%s - w0/2', x_bl);
                y1_bl = sprintf('%s', y_bl);
                x2_bl = sprintf('%s + w0/2', x_bl);
                y2_bl = sprintf('%s', y_bl);
                clamp_linesString{2} = sprintf('ls%i',2);
                clamp_lines{2} = wp1.geom.create(clamp_linesString{2}, 'LineSegment');
                clamp_lines{2}.set('specify1', 'coord');
                clamp_lines{2}.set('coord1', {x1_bl y1_bl});
                clamp_lines{2}.set('specify2', 'coord');
                clamp_lines{2}.set('coord2', {x2_bl y2_bl});
            end
        end

    end
    uni1 = wp1.geom.create('uni1', 'Union');
    uni1.selection('input').set(segmentsString(1:3));
    uni1.set('intbnd', false);

    uni5 = wp1.geom.create('uni5', 'Union');
    uni5.selection('input').set({'ls1', 'ls2'});
    uni5.set('intbnd', false);
%     
    dif1 = wp1.geom.create('dif1', 'Difference');
    dif1.selection('input').set('uni1');
    dif1.selection('input2').set(segmentsString{4});

    rotates = cell(2*N,1);
    rotatesString = cell(2*N+2,1);
%    rotatesString = cell(fix(N),1);


    mir1 = wp1.geom.create('mir1', 'Mirror');
    mir1.selection('input').set('dif1');
    mir1.set('pos', {'0' '0'});
    mir1.set('axis', {'1' '0'});
    mir1.set('keep', true);
    
    for r = 2:N
        rotatesString{r} = sprintf('rot%i', r);
        rotates{r} = wp1.geom.create(rotatesString{r}, 'Rotate');
        if mod(r,2) == 1
            rotates{r}.selection('input').set('dif1');
        else
            rotates{r}.selection('input').set('mir1');
        end
        
        rotates{r}.set('rot', sprintf('rotangle*180/pi * (%d-1)', r));
        rotates{r}.set('keep', true)
        rotates{r}.set('pos', {'0' '0'})

    end

    rotatesString{1} = 'rot1';
    rotates{1} = wp1.geom.create(rotatesString{1}, 'Rotate');
    rotates{1}.selection('input').set('dif1');
    rotates{1}.set('rot', '0');
    rotates{1}.set('keep', false)
    rotates{1}.set('pos', {'0' '0'})

    del1 = wp1.geom.create('del1', 'Delete');
    del1.selection('input').init;
    del1.selection('input').set('mir1');

    
    mir2 = wp1.geom.create('mir2', 'Mirror');
    mir2.selection('input').set('uni5');
    mir2.set('pos', {'0' '0'});
    mir2.set('axis', {'1' '0'});
    mir2.set('keep', true);

    for r = N+2:2*N
        rotatesString{r} = sprintf('rot%i', r);
        rotates{r} = wp1.geom.create(rotatesString{r}, 'Rotate');
        if mod(r,2) == 1
            rotates{r}.selection('input').set('uni5');
        else
            rotates{r}.selection('input').set('mir2');
        end
        
        rotates{r}.set('rot', sprintf('rotangle * 180/pi * (%d - N - 1)', r));
        rotates{r}.set('keep', true)
        rotates{r}.set('pos', {'0' '0'})
    end

    rotatesString{N+1} = sprintf('rot%i', N+1);
    rotates{N+1} = wp1.geom.create(rotatesString{N+1}, 'Rotate');
    rotates{N+1}.selection('input').set('uni5');
    rotates{N+1}.set('rot', 0);
    rotates{N+1}.set('keep', false)
    rotates{N+1}.set('pos', {'0' '0'})

    del2 = wp1.geom.create('del2', 'Delete');
    del2.selection('input').init;
    del2.selection('input').set('mir2');

    uni2 = wp1.geom.create('uni2', 'Union');
    uni2.selection('input').set(rotatesString(1:N));
    uni2.set('intbnd', false);

    %{
    mir1 = wp1.geom.create('mir1', 'Mirror');
    mir1.selection('input').set('dif1');
    mir1.set('pos', [0 0]);
    mir1.set('axis', [1 0]);
    mir1.set('keep', true);
    uni2 = wp1.geom.create('uni2', 'Union');
    uni2.selection('input').set({'mir1','dif1'});
    uni2.set('intbnd', false);
    mir2 = wp1.geom.create('mir2', 'Mirror');
    mir2.selection('input').set(clamp_linesString);
    mir2.set('pos', [0 0]);
    mir2.set('axis', [1 0]);
    mir2.set('keep', true);
    %}

    if (w_pad>w0) && (l_pad>0)
        
        paramcurve1 = wp1.geom.create('pc1', 'ParametricCurve');
        paramcurve1.set('parmin', 0);
        paramcurve1.set('parmax', 1);
        paramcurve1.set('pos', {'-l_trans-l_pad/2' 'w0*rw1/2'});
        paramcurve1.set('coord', {'s*l_trans' '0.5* (w_pad-w0*rw1)*(-2*s^3+3*s^2)'});
        ls_pad = wp1.geom.create(sprintf('ls%i',3), 'LineSegment');
        ls_pad.set('specify1', 'coord');
        ls_pad.set('coord1', [-l_pad/2 w_pad/2]);
        ls_pad.set('specify2', 'coord');
        ls_pad.set('coord2', [l_pad/2 w_pad/2]);
        paramcurve2 = wp1.geom.create('pc2', 'ParametricCurve');
        paramcurve2.set('parmin', 0);
        paramcurve2.set('parmax', 1);
        paramcurve2.set('pos', [l_trans+l_pad/2 w0*rw1/2]);
        paramcurve2.set('coord', {'-s*l_trans' '0.5* (w_pad-w0*rw1)*(-2*s^3+3*s^2)'});
        mir3 = wp1.geom.create('mir3', 'Mirror');
        mir3.selection('input').set({'pc1', 'pc2', sprintf('ls%i',2)});
        mir3.set('pos', [0 0]);
        mir3.set('axis', [0 1]);
        mir3.set('keep', true);
        
        rec_pad = wp1.geom.create(sprintf('r%i',5), 'Rectangle');
        rec_pad.set('size', [l_pad+2*l_trans w0*rw1]);
        rec_pad.set('base', 'center');
        rec_pad.set('pos', [0 0]);       
        csol1 = wp1.geom.create('csol1', 'ConvertToSolid');
        csol1.selection('input').set({'pc1', 'pc2', sprintf('ls%i',3), 'mir3',sprintf('r%i',5)});

        rotatesString{2*N+1} = sprintf('rot%i', 2*N+1);
        rotates{2*N+1} = wp1.geom.create(rotatesString{2*N+1}, 'Rotate');
        rotates{2*N+1}.selection('input').set('csol1');
        rotates{2*N+1}.set('rot', sprintf('%s*180/pi',theta_seg{2}));
        rotates{2*N+1}.set('keep', false)
        rotates{2*N+1}.set('pos', [0 0])
    
    
        mov1 = wp1.geom.create('mov1', 'Move');
        mov1.selection('input').set(rotatesString{2*N+1});
        
        mov1.set('displ', {sprintf('%s - (l0 * rl1 * cos(%s) / 2)',x_seg{1}, theta_seg{2}) sprintf('Rady-(l0*rl1*sin(%s)/2)', theta_seg{2})});


        %mov1.set('displ', sprintf('{%s - (l0 * rl1 * cos(%s) / 2) Rady - (l0 * rl1 * sin(%s) / 2)}', ...
        %    x_seg{1}, theta_seg{2},  theta_seg{2}));


        mov1.set('keep', false);
    
        mir4 = wp1.geom.create('mir4', 'Mirror');
        mir4.selection('input').set('mov1');
        mir4.set('pos', [0 0]);
        mir4.set('axis', [1 -1]);
        mir4.set('keep', true);
    
    
        uni3 = wp1.geom.create('uni3', 'Union');
        uni3.selection('input').set({'uni2','mov1','mir4'});
        uni3.set('intbnd', false);


        if pad_trigger

            paramcurve3 = wp1.geom.create('pc3', 'ParametricCurve');
            paramcurve3.set('parmin', 0);
            paramcurve3.set('parmax', 1);
            paramcurve3.set('pos', [-l_trans-l_pad/2 w0*rw2/2]);
            paramcurve3.set('coord', {'s*l_trans' '0.5* (w_pad-w0*rw2)*(-2*s^3+3*s^2)'});
            ls_pad = wp1.geom.create(sprintf('ls%i',4), 'LineSegment');
            ls_pad.set('specify1', 'coord');
            ls_pad.set('coord1', [-l_pad/2 w_pad/2]);
            ls_pad.set('specify2', 'coord');
            ls_pad.set('coord2', [l_pad/2 w_pad/2]);
            paramcurve4 = wp1.geom.create('pc4', 'ParametricCurve');
            paramcurve4.set('parmin', 0);
            paramcurve4.set('parmax', 1);
            paramcurve4.set('pos', [l_trans+l_pad/2 w0*rw2/2]);
            paramcurve4.set('coord', {'-s*l_trans' '0.5* (w_pad-w0*rw2)*(-2*s^3+3*s^2)'});
            mir5 = wp1.geom.create('mir5', 'Mirror');
            mir5.selection('input').set({'pc3', 'pc4', sprintf('ls%i',3)});
            mir5.set('pos', [0 0]);
            mir5.set('axis', [0 1]);
            mir5.set('keep', true);

            rec_pad = wp1.geom.create(sprintf('r%i',6), 'Rectangle');
            rec_pad.set('size', [l_pad+2*l_trans w0*rw2]);
            rec_pad.set('base', 'center');
            rec_pad.set('pos', [0 0]);       
            csol1 = wp1.geom.create('csol2', 'ConvertToSolid');
            csol1.selection('input').set({'pc3', 'pc4', sprintf('ls%i',4), 'mir5',sprintf('r%i',6)});


            rotatesString{2*N+2} = sprintf('rot%i', 2*N+2);
            rotates{2*N+2} = wp1.geom.create(rotatesString{2*N+2}, 'Rotate');
            rotates{2*N+2}.selection('input').set('csol2');
            rotates{2*N+2}.set('rot', sprintf('%s*180/pi', theta_seg{3}));
            rotates{2*N+2}.set('keep', false)
            rotates{2*N+2}.set('pos', [0 0])
        
        
            mov2 = wp1.geom.create('mov2', 'Move');
            mov2.selection('input').set(rotatesString{2*N+2});
            mov2.set('displ', {sprintf('%s-(l0*rl2*cos(%s)/2)',x_seg{1},theta_seg{3}) sprintf('Rady-(l0*rl2*sin(%s)/2)',theta_seg{3})});
            mov2.set('keep', false);
        
            mir6 = wp1.geom.create('mir6', 'Mirror');
            mir6.selection('input').set('mov2');
            mir6.set('pos', [0 0]);
            mir6.set('axis', [1 1]);
            mir6.set('keep', true);
        
        
            uni4 = wp1.geom.create('uni4', 'Union');
            uni4.selection('input').set({'uni3','mov2','mir6'});
            uni4.set('intbnd', false);
        end



    end
    
    



    geom1.run;
    
    if plot_flag
        mphgeom(model,'geom1', 'vertexlabels', 'off','facemode','off')
        zlim([-10*h_mbr,10*h_mbr])
        view(0,90)
        xlim([-2*l0,2*l0])
        ylim([-2*l0,2*l0])
    end


       %% Selection
    %     hold on
        idx_fixed = [];
        
        fprintf('%s-1\n','rotangle');
        
        disp(rotangle);
        
        for i = 1:1
            if i > 0
                l_cut =  sprintf('%s / max(abs(cos(%s)), abs(sin(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
                theta_mod = sprintf('mod(%s, 2*pi)', theta_seg{i});
                dl_seg = sprintf('0.5 * %s * min(abs(cot(%s)), abs(tan(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
                dx_cut = sprintf('%s * fix(sqrt(2) * cos(%s)) / 2', l_cut, theta_mod);
                dy_cut = sprintf('%s * fix(sqrt(2) * sin(%s)) / 2', l_cut, theta_mod);
                if mod(eval(theta_mod),pi/4) == 0 && mod(eval(theta_mod), pi/2) ~= 0
                    if cos(theta_mod) < 1e-10
                        dx_cut = '0';
                    else
                        dx_cut = sprintf('sign(cos(%s)) * %s / 2', theta_mod, l_cut);
                    end
                dy_cut = 0;
                end
                x_cut = sprintf('%s + (%s - 2 * %s) * cos(%s) / 2 + %s', ...
                    x_seg{i}, l_seg{i}, dl_seg, theta_seg{i}, dx_cut);
                y_cut = sprintf('%s + (%s - 2 * %s) * sin(%s) / 2 + %s', ...
                    y_seg{i}, l_seg{i}, dl_seg, theta_seg{i}, dy_cut);
                l_select = sprintf('%s + 1e-9', l_cut);
                x_select_1 = sprintf('%s - %s / 2', x_cut, l_select);
                y_select_1 = sprintf('%s - %s / 2', y_cut, l_select);
                x_select_2 = sprintf('%s + %s / 2', x_cut, l_select);
                y_select_2 = sprintf('%s + %s / 2', y_cut, l_select);
                
                
                for n = 1:N
                    R= ...
                    [...
                        cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                        sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                        0, 0, 1 ...
                    ];
                    if mod(n,2)==1
                        box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
                    else
                        box_coordinates = [(-1)*eval(x_select_1), eval(y_select_1), 0; (-1)*eval(x_select_2), eval(y_select_2), 0]';
                    end
                    rotatedCoords = (R*box_coordinates);
    
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');
    
                    idx_fixed = [idx_fixed, idx];
                end
                
                %{
                for n = 1:N
                    
                    angle_str = sprintf('(%d-1) * %s', n, rotangle);
                    disp(angle_str);
                    R = {...
                        sprintf('cos(%s)', angle_str), sprintf('-sin(%s)', angle_str), '0';...
                        sprintf('sin(%s)', angle_str), sprintf('cos(%s)', angle_str), '0';...
                        '0', '0', '1'...
                    };
                    
                    if mod(n, 2) == 1
                        box_coordinates = {x_select_1, y_select_1, '0'; x_select_2, y_select_2, '0'}';
                    else
                        box_coordinates = {-x_select_1, y_select_1, '0'; -x_select_2, y_select_2, '0'}';
                    end
                    
                    rotatedCoords_str = {
                        sprintf('%s * %s + %s * %s + %s * %s', R{1, 1}, box_coordinates{1, 1}, R{1, 2}, box_coordinates{2, 1}, R{1, 3}, box_coordinates{3, 1}),...
                        sprintf('%s * %s + %s * %s + %s * %s', R{1, 1}, box_coordinates{1, 2}, R{1, 2}, box_coordinates{2, 2}, R{1, 3}, box_coordinates{3, 2});...
    
                        sprintf('%s * %s + %s * %s + %s * %s', R{2, 1}, box_coordinates{1, 1}, R{2, 2}, box_coordinates{2, 1}, R{2, 3}, box_coordinates{3, 1}),...
                        sprintf('%s * %s + %s * %s + %s * %s', R{2, 1}, box_coordinates{1, 2}, R{2, 2}, box_coordinates{2, 2}, R{2, 3}, box_coordinates{3, 2});...
    
                        sprintf('%s * %s + %s * %s + %s * %s', R{3, 1}, box_coordinates{1, 1}, R{3, 2}, box_coordinates{2, 1}, R{3, 3}, box_coordinates{3, 1}),...
                        sprintf('%s * %s + %s * %s + %s * %s', R{3, 1}, box_coordinates{1, 2}, R{3, 2}, box_coordinates{2, 2}, R{3, 3}, box_coordinates{3, 2})...
                    };
                    
                    rotatedCoords = cellfun(@(c) eval(c), rotatedCoords_str);
                    idx = mphselectbox(model, 'geom1', rotatedCoords,'edge');
                    
                    idx_fixed = {idx_fixed, idx};
                end
                %}
    
    
    
    
    
    
    
    
            end
        end
        sel1 = model.selection.create('sel1').geom(1);
        sel1.set(idx_fixed);
        sel1.label('clamping edges');
        idx_clamps = [];
        idx_clamps_alledges = [];
    
    
        for i = 1:1
            if i > 0
                l_cut =  sprintf('%s / max(abs(cos(%s)), abs(sin(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
                theta_mod = sprintf('mod(%s, 2*pi)', theta_seg{i});
                dl_seg = sprintf('0.5 * %s * min(abs(cot(%s)), abs(tan(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
                dx_cut = sprintf('%s * fix(sqrt(2) * cos(%s)) / 2', l_cut, theta_mod);
                dy_cut = sprintf('%s * fix(sqrt(2) * sin(%s)) / 2', l_cut, theta_mod);
                if mod(eval(theta_mod),pi/4) == 0 && mod(eval(theta_mod), pi/2) ~= 0
                    if cos(theta_mod) < 1e-10
                        dx_cut = '0';
                    else
                        dx_cut = sprintf('sign(cos(%s)) * %s / 2', theta_mod, l_cut);
                    end
                dy_cut = 0;
                end
                x_cut = sprintf('%s + (%s - 2 * %s - %s / 2) * cos(%s) / 2 ', ...
                    x_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i});
                y_cut = sprintf('%s + (%s - 2 * %s - %s / 2) * sin(%s) / 2 ', ...
                    y_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i});
                l_select = sprintf('%s + %s + 1e-9', l_clamp, l_cut);
                x_select_1 = sprintf('%s - %s / 2', x_cut, l_select);
                y_select_1 = sprintf('%s - %s / 2', y_cut, l_select);
                x_select_2 = sprintf('%s + %s / 2', x_cut, l_select);
                y_select_2 = sprintf('%s + %s / 2', y_cut, l_select);
    
    
                for n = 1:N
                    R= ...
                    [...
                        cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                        sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                        0, 0, 1 ...
                    ];    
                    if mod(n,2)==1
                        box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
                    else
                        box_coordinates = [(-1)*eval(x_select_1), eval(y_select_1), 0; (-1)*eval(x_select_2), eval(y_select_2), 0]';
                    end    
                    rotatedCoords = (R*box_coordinates);    
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'boundary');    
                    idx_clamps = [idx_clamps, idx];
                end
    
                for n = 1:N
                    R= ...
                    [...
                        cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                        sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                        0, 0, 1 ...
                    ];
                    if mod(n,2)==1
                        box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
                    else
                        box_coordinates = [(-1)*eval(x_select_1), eval(y_select_1), 0; (-1)*eval(x_select_2), eval(y_select_2), 0]';
                    end
                    rotatedCoords = (R*box_coordinates);    
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');    
                    idx_clamps_alledges = [idx_clamps_alledges, idx];
                end
    
            end
        end
        sel2 = model.selection.create('sel2').geom(2);
        sel2.set(idx_clamps);
        sel2.label('clamping pads');
        
        idx_clamps_x_inner = [];
        for i = 1:1
            if i > 0
                l_cut =  sprintf('%s / max(abs(cos(%s)), abs(sin(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
                theta_mod = sprintf('mod(%s, 2*pi)', theta_seg{i});
                dl_seg = sprintf('0.5 * %s * min(abs(cot(%s)), abs(tan(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
                dx_cut = sprintf('%s * fix(sqrt(2) * cos(%s)) / 2', l_cut, theta_mod);
                dy_cut = sprintf('%s * fix(sqrt(2) * sin(%s)) / 2', l_cut, theta_mod);
                if mod(eval(theta_mod),pi/4) == 0 && mod(eval(theta_mod), pi/2) ~= 0
                    if cos(theta_mod) < 1e-10
                        dx_cut = '0';
                    else
                        dx_cut = sprintf('sign(cos(%s)) * %s / 2', theta_mod, l_cut);
                    end
                dy_cut = 0;
                end
                x_cut = sprintf('%s + (%s - 2 * %s - %s) * cos(%s) / 2 - %s', ...
                    x_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i}, dx_cut);
                y_cut = sprintf('%s + (%s - 2 * %s - %s) * sin(%s) / 2 - %s', ...
                    y_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i}, dy_cut);
                l_select = sprintf('%s + 1e-9', l_cut);
                x_select_1 = sprintf('%s - %s / 2', x_cut, l_select);
                y_select_1 = sprintf('%s - %s / 2', y_cut, l_select);
                x_select_2 = sprintf('%s + %s / 2', x_cut, l_select);
                y_select_2 = sprintf('%s + %s / 2', y_cut, l_select);
    
                for n = 1:N
                    R= ...
                    [...
                        cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                        sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                        0, 0, 1 ...
                    ];
                    if mod(n,2)==1
                        box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
                    else
                        box_coordinates = [(-1)*eval(x_select_1), eval(y_select_1), 0; (-1)*eval(x_select_2), eval(y_select_2), 0]';
                    end
                    rotatedCoords = (R*box_coordinates);
    
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');
    
                    idx_clamps_x_inner = [idx_clamps_x_inner, idx];
                end
                
            end
        end
        idx_clamps_y = setdiff(idx_clamps_alledges,[idx_fixed, idx_clamps_x_inner]);
        sel3 = model.selection.create('sel3').geom(1);
        sel3.set(idx_clamps_y);
        sel3.label('clamping edges y');
        idx_all = mphselectbox(model, 'geom1',...
                    [-((2+rl1+rl2)*l0),-(2*l0+(rl1+rl2)*l0),0;...
                     (2*l0+(rl1+rl2)*l0),(2*l0+(rl1+rl2)*l0),0]', 'edge');
        idx_segments = setdiff(idx_all, idx_clamps_alledges);
        
        sel4 = model.selection.create('sel4').geom(1);
        sel4.set(idx_segments);
        sel4.label('segments edges');
        
        idx_all_boundaries = mphselectbox(model, 'geom1',...
                    [-((2+rl1+rl2)*l0),-(2*l0+(rl1+rl2)*l0),0;...
                     (2*l0+(rl1+rl2)*l0),(2*l0+(rl1+rl2)*l0),0]', 'boundary');
        idx_segments_area = setdiff(idx_all_boundaries, idx_clamps);
        
        sel5 = model.selection.create('sel5').geom(2);
        sel5.set(idx_segments_area);
        sel5.label('beam area');

        idx_polygon_boundary = mphselectbox(model, 'geom1',...
                    [-((rw1+rw2)*w0+Rady),-((rw1+rw2)*w0+Rady),0;...
                     ((rw1+rw2)*w0+Rady),((rw1+rw2)*w0+Rady),0]', 'boundary');
        idx_polygon_area = setdiff(idx_polygon_boundary, idx_clamps);
        
        sel6 = model.selection.create('sel6').geom(2);
        sel6.set(idx_polygon_area);
        sel6.label('polygon area'); 

        idx_support_area = setdiff(idx_segments_area, idx_polygon_area);
        sel7 = model.selection.create('sel7').geom(2);
        sel7.set(idx_support_area);
        sel7.label('support area');



%% Material
    mat1 = model.material.create('mat1');
    mat1.propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
    mat1.label('Si3N4 - Silicon nitride');
    mat1.propertyGroup('def').set('density', '3100[kg/m^3]');
    mat1.propertyGroup('Enu').set('E', '250e9[Pa]');
    mat1.propertyGroup('Enu').set('nu', '0.23');
    mat1.set('family', 'plastic');
    mat1.selection.all;

%% Probes
    bnd1 = model.component('mod1').probe.create('bnd1', 'Boundary');
    bnd1.set('intsurface', true);
    bnd1.label('OOP bending energy');
    bnd1.set('probename', 'bend_energy_op');
    bnd1.set('type', 'integral');
    bnd1.set('expr', 'shell.E*h_mbr^3*((dtang(dtang(w,x),x)+dtang(dtang(w,y),y))^2+2*(1-shell.nu)*(dtang(dtang(w,x),y)^2-dtang(dtang(w,x),x)*dtang(dtang(w,y),y)))/(12*(1-shell.nu^2))');

    bnd2 = model.component('mod1').probe.create('bnd2', 'Boundary');
    bnd2.set('intsurface', true);
    bnd2.label('IP bending energy');
    bnd2.set('probename', 'bend_energy_ip');
    bnd2.set('type', 'integral');
    bnd2.set('expr', 'shell.E*h_mbr*((1-shell.nu)*((dtang(u,x)^2)+(dtang(v,y)^2))+2*shell.nu*dtang(u,x)*dtang(v,y)+0.5*(1-2*shell.nu)*((dtang(u,y)+dtang(v,x))^2))/((1+shell.nu)*(1-2*shell.nu))');

    bnd3 = model.component('mod1').probe.create('bnd3', 'Boundary');
    bnd3.set('intsurface', true);
    bnd3.label('Kinetic energy');
    bnd3.set('probename', 'kin_energy');
    bnd3.set('type', 'integral');
    bnd3.set('expr', 'h_mbr*shell.rho*shell.omega^2*shell.disp^2');

    Maximum = model.component('mod1').cpl.create('maxop1', 'Maximum');
    Maximum.selection.geom('geom1', 2);
    Maximum.selection.named('sel5');

    Average = model.component('mod1').cpl.create('aveop1', 'Average');
    Average.selection.geom('geom1', 2);
    Average.selection.named('sel5');

    Maxpolygon = model.component('mod1').cpl.create('maxop2', 'Maximum');
    Maxpolygon.selection.geom('geom1', 2);
    Maxpolygon.selection.named('sel6');

    Maxsupport = model.component('mod1').cpl.create('maxop3', 'Maximum');
    Maxsupport.selection.geom('geom1', 2);
    Maxsupport.selection.named('sel7');
    
    var1 = model.component('mod1').variable.create('var1');
    var1.set('Qm', 'Qint * kin_energy / max(bend_energy_op, bend_energy_ip)');
    var1.set('DQ', 'kin_energy/max(bend_energy_op, bend_energy_ip)')
    var1.set('m_eff', 'kin_energy / ((2*pi*shell.freq)^2 * maxop1(shell.disp)^2)');
    var1.set('eta_op', 'sqrt(maxop1(w^2)/maxop1(u^2+v^2))>1');
    var1.set('amp_peri', 'sqrt(maxop2(w^2)/maxop3(w^2))>10');
%% Physics
    shell = model.physics.create('shell', 'Shell', 'geom1');
    shell.feature('to1').set('d', 'h_mbr');
    fix1 = shell.create('fix1', 'Fixed', 1);
    fix1.selection.named('sel1');
    iss1 = shell.feature('emm1').create('iss1', 'InitialStressandStrain', 2);
    iss1.set('Ni', {'stress*h_mbr' '0' '0' 'stress*h_mbr'});



    %% Mesh
    mesh = model.mesh.create('mesh', 'geom1');
%         mesh.autoMeshSize(4);
    mesh_clamp = mesh.create('map1', 'Map');
    mesh_clamp.selection.set(idx_clamps);
    dis1 = mesh_clamp.create('dis1', 'Distribution');
    dis1.selection.set(idx_fixed);
    dis1.set('numelem', n_clamp_x); 
    dis2 = mesh_clamp.create('dis2', 'Distribution');
    dis2.selection.set(idx_clamps_y);
    dis2.set('numelem', n_clamp_y);
    
    ftri1 = mesh.feature.create('ftri1', 'FreeTri');
    size1 = ftri1.create('size', 'Size');
    size1.set('hmax', mesh_max);
    size1.set('hmin', mesh_min);
    ftri1.selection.set(idx_segments_area);
   
    
    mesh.run;

%% Study
    
    std1 = model.study.create('std1');    % study  node

    %{
    param = std1.create('param', "Parametric");
    param.set('sweeptype', 'sparse');       
    param.setIndex('pname', 'rl2', 0);
    param.setIndex('plistarr', '{range(2.45,0.1,2.55)}', 0);
    %}

    stat = std1.create('stat', 'Stationary');
    stat.set('geometricNonlinearity', true);
    
    eig = std1.feature.create('eig', 'Eigenfrequency');  % Mode analysis solver
    eig.set('geometricNonlinearity', true);
    eig.set('neigsactive', true);   % Enable entering number of solutions
    eig.set('neigs', N_sol);  % number of solutions
    eig.set('eigunit', 'kHz'); % solutions unit
    eig.set('shift', num2str(f0/1e3)); % Solve for frequency around 194THz

%{
    model.sol.create('sol3');
    model.sol('sol3').study('std1');
    model.sol('sol3').label('Parametric Solutions 1');

    model.result.dataset.create('dset4', 'Solution');
    model.result.dataset('dset4').set('solution', 'sol3');
    model.result.numerical.create('gev1', 'EvalGlobal');
    model.result.numerical('gev1').set('data', 'dset4');

    model.result.table.create('tbl2', 'Table');
    model.result.table('tbl2').comments('Global Evaluation 1');

    model.result.numerical('gev1').set('table', 'tbl2');
    model.result.numerical('gev1').set('expr', {'DQ' 'eta_op' 'Qm' 'm_eff'});
    model.result.numerical('gev1').set('unit', {'1' '' '1' 'kg'});
    model.result.numerical('gev1').set('descr', {'' '' '' ''});
    model.result.numerical('gev1').set('const', {'shell.refpntx' '0' 'Reference point for moment computation, x-coordinate';  ...
    'shell.refpnty' '0' 'Reference point for moment computation, y-coordinate';  ...
    'shell.refpntz' '0' 'Reference point for moment computation, z-coordinate';  ...
    'shell.z' '1' 'Local z-coordinate [-1,1] for thickness-dependent results'});
%}
    steps = 0;
    pgsegments = cell(length(values), 1);

    index_match = [];
    freqs_match = [];
    eta_match = [];
    amp_peri_match = [];
    Q_match = [];
    m_eff_match = [];
    rl2_match = [];
    S_F_match = [];
    D_Q_match = [];

    for setparamator = values
        model.param.set('rl2', setparamator);
        model.study('std1').run;

        %% Global evaluation
    %{
        eta = mphglobal(model, 'sqrt(maxop1(w^2)/maxop1(u^2+v^2))')>1;
        % ind_ip = find(eta<1);
        ind_ip = eta < 1;
        ind_op = find(eta == 1);
        Freqs = mphglobal(model, 'shell.freq');
        disp(Freqs);
        W_bend = mphglobal(model, 'bend_energy_op');
        W_bend_ip = mphglobal(model, 'bend_energy_ip');
        W_bend(ind_ip) = W_bend_ip(ind_ip);
        W_kin = mphglobal(model, 'kin_energy');
        Q = Qint * W_kin ./ W_bend;
        m_eff = W_kin ./ ((2*pi*Freqs).^2 .* mphglobal(model,'maxop1(shell.disp)').^2);
        S_F = 8*pi * 1.38e-23 * 300 * m_eff .* Freqs ./ Q;
    %}
        eta = mphglobal(model, 'sqrt(maxop1(w^2)/maxop1(u^2+v^2))') > 1;
        amp_peri = mphglobal(model, 'sqrt(maxop2(w^2)/maxop3(w^2))') > 10;
        ind_ip = eta < 1;
        %ind_op = find(eta == 1);
        ind_op = find(eta == 1 & amp_peri == 1);
        Freqs = mphglobal(model, 'shell.freq');
        W_bend = mphglobal(model, 'bend_energy_op');
        W_bend_ip = mphglobal(model, 'bend_energy_ip');
        W_bend(ind_ip) = W_bend_ip(ind_ip);
        W_kin = mphglobal(model, 'kin_energy');
        Q = Qint * W_kin ./ W_bend;
        m_eff = W_kin ./ ((2*pi*Freqs).^2 .* mphglobal(model, 'maxop1(shell.disp)').^2);
        S_F = 8 * pi * 1.38e-23 * 300 * m_eff .* Freqs ./ Q;
        DQ = W_kin ./ W_bend;



        for j = 1:length(ind_op)
            i = ind_op(j);
            %index_match = [inex_match; i];
            freqs_match = [freqs_match; Freqs(i)];
            eta_match = [eta_match; eta(i)];
            amp_peri_match = [amp_peri_match; amp_peri(i)];
            Q_match = [Q_match; Q(i)];
            m_eff_match = [m_eff_match; m_eff(i)];
            rl2_match = [rl2_match; setparamator];
            S_F_match = [S_F_match; S_F(i)];
            D_Q_match = [D_Q_match; DQ(i)];
        end
    



        
    %% creating a 3D plot group

        steps = steps + 1;
        if (plot_flag && plot_op_flag) & steps ==1 
            pgsegments{steps} = model.result.create(sprintf('pg%d', steps), 'PlotGroup3D');
            pgsegments{steps}.create('surf1', 'Surface');
            pgsegments{steps}.feature('surf1').set('wireframe', true);
            pgsegments{steps}.feature('surf1').create('def1', 'Deform');
            pgsegments{steps}.run;
            fig1 = figure;
            fig1.Position = [0,0,500,500];
            for j = 1:length(ind_op)
                i = ind_op(j);
                subplot(2,ceil(length(ind_op)/2),j)
                pgsegments{steps}.setIndex('looplevel', i, 0);
                mphplot(model, sprintf('pg%d', steps));
                title({['f = ', num2str(Freqs(i)/1e6), ' MHz, ', '\eta = ', num2str(eta(i))], ...
                       ['PerimeterMode = ', num2str(amp_peri(i)), ', Q = ', num2str(Q(i)/1e6,3), ' x 10^6'], ...
                       ['m_{eff} = ', num2str(m_eff(i)*1e15,3), ' pg, rl2 = ', num2str(rl2)], ...
                       ['S_F^{th} = ', num2str(S_F(i)*1e36,3), ' aN^2/Hz, DQ = ', num2str(DQ(i))]})
            end
        end
    end


        %{
        
    if plot_flag && plot_op_flag
        pgsegments{steps} = model.result.create(sprintf('pg%d', steps), 'PlotGroup3D');
        pgsegments{steps}.create('surf1', 'Surface');
        pgsegments{steps}.feature('surf1').set('wireframe', true);
        pgsegments{steps}.feature('surf1').create('def1', 'Deform');
        pgsegments{steps}.run;
        fig1 = figure;
        fig1.Position = [0,0,500,500];        
        for i = 1:length(freqs_match) % 
            %i = ind_op(j);
            subplot(2,ceil(length(freqs_match)/2),i)
            %pgsegments{steps}.setIndex('looplevel', i, 0);
            mphplot(model, sprintf('pg%d', steps));
            %title({['f = ',num2str(freqs_match(i)/1e6),' MHz'],...
            %       ['\eta = ', num2str(eta_match(i))],...
            %       ['PerimeterMode = ', num2str(amp_peri_match(i))],...
            %       ['Q = ',num2str(Q_match(i)/1e6,3), ' x 10^6'],...
            %       ['m_{eff} = ',num2str(m_eff_match(i)*1e15,3), ' pg'],...
            %       ['rl2 = ', num2str(rl2_match(i))],...
            %       ['S_F^{th} = ', num2str(S_F_match(i)*1e36,3), ' aN^2/Hz']})
            title({['f = ', num2str(freqs_match(i)/1e6), ' MHz, ', '\eta = ', num2str(eta_match(i))], ...
                   ['PerimeterMode = ', num2str(amp_peri_match(i)), ', Q = ', num2str(Q_match(i)/1e6,3), ' x 10^6'], ...
                   ['m_{eff} = ', num2str(m_eff_match(i)*1e15,3), ' pg, rl2 = ', num2str(rl2_match(i))], ...
                   ['S_F^{th} = ', num2str(S_F_match(i)*1e36,3), ' aN^2/Hz']})
        end
    end
        %}
        


    %end
%}

Periresults = table (freqs_match, eta_match, amp_peri_match, Q_match, D_Q_match, m_eff_match, rl2_match, S_F_match);
figTable = figure;
figTable.Name = sprintf('Periresults Table PadTrigger = %d', pad_trigger); 
uitable('Data',Periresults{:,:},'ColumnName',Periresults.Properties.VariableNames,...
    'RowName',Periresults.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

values_str = sprintf('%.2f_', values);
values_str = values_str(1:end-1);
filename = sprintf('periresult_%s.jpg', values_str);
frame = getframe(figTable);
imwrite(frame.cdata, filename);
mphsave(model, 'Practice_Resonator_four_pads_sweep.mph')
end
