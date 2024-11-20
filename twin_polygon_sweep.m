function [Freqs, Q ,m_eff, S_F, eta, rl2_match, Q_match] = ...
         twin_polygon_sweep(stress, h_mbr, l0, w0, N,...
                           N_sol,rl1, rl2, rw1, rw2, lc, wc,...
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
    rotangle = (pi);
    %rotvalue = (2*pi/N);
    Diameter = sqrt((l0*rl1)^2+(l0*rl2)^2);
    rad = Diameter/2;
    theta = (pi-rotangle)/2;
    Radx = rad*sin(atan(max(rl1/rl2,rl2/rl1))-pi/4);
    Rady = rad*cos(atan(max(rl1/rl2,rl2/rl1))-pi/4);
    segments = cell(7,1);
    segmentsString = cell(7,1);
    x_seg = cell(5,1);
    y_seg = cell(5,1);
    theta_seg = cell(5,1);
    l_seg = cell(5,1);
    w_seg = cell(5,1);
    %disp(theta_seg{1});

   
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
    model.param.set('wc', sprintf('%d[m]', wc), 'Width of the coupler segment');
    model.param.set('lc', sprintf('%d[m]', lc), 'Length of the coupler segment');
    model.param.set('Qint', '12500 * (h_mbr/100e-9[m])');
    model.param.descr('Qint', 'Intrinsic Q');    
    model.param.set('rotangle', '(pi)');
    model.param.descr('rotangle', 'Angle of rotation');
    model.param.set('theta', '((pi-rotangle)/2)');
    model.param.descr('theta', 'Branching angle (radians)');
    model.param.set('Diameter', 'sqrt((l0*rl1)^2+(l0*rl2)^2)');
    model.param.descr('Diameter', 'Diameter of circumscribed circle');
    model.param.set('rad', '(Diameter/2)');
    model.param.descr('rad', 'Radius of circumscribed circle');
    model.param.set('Radx', '(rad*sin(atan(max(rl1/rl2,rl2/rl1))-pi/4))');
    model.param.descr('Radx', 'x-coordinate of polygon vertex');
    model.param.set('Rady', '(rad*cos(atan(max(rl1/rl2,rl2/rl1))-pi/4))');
    model.param.descr('Rady', 'y-coordinate of polygon vertex');

    
    
    %disp(rotvale);
    %rw = 1/rW1;
    
    l_seg{1} = '(lc/2)';
    w_seg{1} = 'wc';
    theta_seg{1} = '0';
    segmentsString{1} = 'r1';

    x_seg{1} = '(lc/4)';
    y_seg{1} = '0';
    
    l_clamp = '(50 * h_mbr)';
    n_clamp_x = 10;
    n_clamp_y = 30;
    mesh_max = 1e-6;
    mesh_min = mesh_max/2;
    
    clamp_lines = cell(5,1);
    clamp_linesString = cell(5,1);

    for i = 2:5
        %ind = i;
        %ind_p = 1;
        segmentsString{i} = sprintf('r%i',i);
%        j = round(i/2);
%        if mod(j, 2) == 1 & j ~= 1
%            j = j + 1;
%        end
        j = i - 1;
        %dl_rw = '0';
        dl_rl = '0';
        fprintf('j = %d\n', j);
        %disp(theta_seg{ind_p});
        
        if i == 2            
            l_seg{i} = '(l0*rl1)';
            w_seg{i} = '(w0*rw1)';            
            theta_seg{i} = sprintf('(%s + (pi/4))', theta_seg{j});
            disp(eval(theta_seg{i}));
        end
        if i == 3 || i == 5   
            l_seg{i} = 'l0';
            w_seg{i} = 'w0';            
        
            theta_seg{i} = sprintf('(%s + (pi/4))', theta_seg{j});
            disp(eval(theta_seg{i}));
                        
            dl_seg = sprintf('(0.5* %s * min(abs(cot(%s)), abs(tan(%s))))', w_seg{i}, theta_seg{i}, theta_seg{i});
            l_seg{i} = sprintf('(%s - %s)', l_seg{i}, dl_seg);
            
        end
        if i == 4          
            l_seg{i} = '(l0*rl2)';
            w_seg{i} = '(w0*rw2)';            
        
            theta_seg{i} = sprintf('(%s - pi/4 * 3)', theta_seg{j});
            disp(eval(theta_seg{i}));
        end

        dl_rw = sprintf('(0.5 * %s * (cos(pi/4) - (%s/%s)) / sin(pi/4))', w_seg{j}, w_seg{i}, w_seg{j});
        %dl_rw = '0';   

        if eval(dl_rw)<0
            dl_rw = '0';
            fprintf('i = %d\n', i);
            dl_rl = sprintf('(0.5 * %s * (cos(pi/4) - (%s/%s)) / sin(pi/4))', w_seg{i}, w_seg{j}, w_seg{i});
        end
       
        x_seg{i} = sprintf('%s + (%s/2-%s) * cos(%s) + (%s/2-%s) * cos(%s)', x_seg{j}, l_seg{j}, dl_rw, theta_seg{j}, l_seg{i}, dl_rl, theta_seg{i} );
        
        
        if i == 4
            y_seg{i} = sprintf('%s - (%s/2-%s) * sin(%s) + (%s/2-%s) * sin(%s)', y_seg{j}, l_seg{j}, dl_rw, theta_seg{j}, l_seg{i}, dl_rl,  theta_seg{i} );
        else
            y_seg{i} = sprintf('%s + (%s/2-%s) * sin(%s) + (%s/2-%s) * sin(%s)', y_seg{j}, l_seg{j}, dl_rw, theta_seg{j}, l_seg{i}, dl_rl,  theta_seg{i} );
        end
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
    disp(theta_seg{2});
    disp(theta_seg{3});
    disp(theta_seg{4});
    disp(theta_seg{5});
    
    for i = 1:5
        segments{i} = wp1.geom.feature.create(segmentsString{i}, 'Rectangle');
        segments{i}.set('size', {l_seg{i} w_seg{i}});
        segments{i}.set('base', 'center');
        segments{i}.set('rot', sprintf('%s*180/pi', theta_seg{i}));
        segments{i}.set('pos', {x_seg{i} y_seg{i}});


        
        

        if i == 3 || i == 5
            l_cut =  sprintf('%s / max(abs(cos(%s)), abs(sin(%s)))',w_seg{i}, theta_seg{i},theta_seg{i});
            theta_mod = sprintf('mod(%s, 2*pi)', theta_seg{i});
            dl_seg = sprintf('0.5*%s * min(abs(cot(%s)), abs(tan(%s)))', w_seg{i}, theta_seg{i}, theta_seg{i});
            dx_cut = sprintf('%s * fix(sqrt(2)*cos(%s))/2', l_cut, theta_mod);
            dy_cut = sprintf('%s * fix(sqrt(2)*sin(%s))/2', l_cut, theta_mod);
            fprintf('theta_mod = %d\n', eval(theta_mod));
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
            
            if i == 3
                i_cut = 6;
                i_cl = 1;
            else
                i_cut = 7;
                i_cl = 3;
            end
            %i_cut = i+3;
            segmentsString{i_cut} = sprintf('r%i',i_cut);
            segments{i_cut} = wp1.geom.feature.create(segmentsString{i_cut}, 'Rectangle');
            segments{i_cut}.set('size', {l_cut l_cut});
            segments{i_cut}.set('base', 'center');
            segments{i_cut}.set('pos', {x_cut y_cut});
            
            %i_cl = 1;
            %if i_cl>0
            x_cl = sprintf('%s + (%s - 2 * %s - %s) * cos(%s) / 2', ...
                x_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i});
            y_cl = sprintf('%s + (%s - 2 * %s - %s) * sin(%s) / 2', ...
                y_seg{i}, l_seg{i}, dl_seg, l_clamp, theta_seg{i});
            l_cl = l_cut;

            theta_cl = sprintf('fix(sqrt(2)*cos(%s)) * pi / 2', theta_mod);
            if mod(eval(theta_mod),pi/4) < 1e-10 && mod(eval(theta_mod),pi/2) > 1e-5
                theta_cl = 'pi/2';
            end
            fprintf('theta_cl = %d\n', eval(theta_cl));
            x1_cl = sprintf('%s - %s * cos(%s) / 2', x_cl, l_cl, theta_cl);
            x2_cl = sprintf('%s + %s * cos(%s) / 2', x_cl, l_cl, theta_cl);
            y1_cl = sprintf('%s - %s * sin(%s) / 2', y_cl, l_cl, theta_cl);
            y2_cl = sprintf('%s + %s * sin(%s) / 2', y_cl, l_cl, theta_cl);
            clamp_linesString{i_cl} = sprintf('ls%i',i_cl);
            clamp_lines{i_cl} = wp1.geom.create(clamp_linesString{i_cl}, 'LineSegment');
            clamp_lines{i_cl}.set('specify1', 'coord');
            clamp_lines{i_cl}.set('coord1', {x1_cl y1_cl});
            clamp_lines{i_cl}.set('specify2', 'coord');
            clamp_lines{i_cl}.set('coord2', {x2_cl y2_cl});

            x_bl = sprintf('%s-%s/2 * cos(%s)', x_seg{i}, l_seg{i}, theta_seg{i});
            y_bl = sprintf('%s-%s/2 * sin(%s)', y_seg{i}, l_seg{i}, theta_seg{i});
            x1_bl = sprintf('%s - w0/2 * sin(%s)', x_bl, theta_seg{i});
            y1_bl = sprintf('%s - w0/2 * cos(%s)', y_bl, theta_seg{i});
            x2_bl = sprintf('%s + w0/2 * sin(%s)', x_bl, theta_seg{i});
            y2_bl = sprintf('%s + w0/2 * cos(%s)', y_bl, theta_seg{i});
            clamp_linesString{i_cl+1} = sprintf('ls%i',i_cl+1);
            clamp_lines{i_cl+1} = wp1.geom.create(clamp_linesString{i_cl+1}, 'LineSegment');
            clamp_lines{i_cl+1}.set('specify1', 'coord');
            clamp_lines{i_cl+1}.set('coord1', {x1_bl y1_bl});
            clamp_lines{i_cl+1}.set('specify2', 'coord');
            clamp_lines{i_cl+1}.set('coord2', {x2_bl y2_bl});
        end

    end

    x1_bl = 'lc/2';
    y1_bl = 'wc/2';
    x2_bl = 'lc/2';
    y2_bl = '-wc/2';
    clamp_linesString{5} = sprintf('ls%i',5);
    clamp_lines{5} = wp1.geom.create(clamp_linesString{5}, 'LineSegment');
    clamp_lines{5}.set('specify1', 'coord');
    clamp_lines{5}.set('coord1', {x1_bl y1_bl});
    clamp_lines{5}.set('specify2', 'coord');
    clamp_lines{5}.set('coord2', {x2_bl y2_bl});


    dif1 = wp1.geom.create('dif1', 'Difference');
    dif1.selection('input').set(segmentsString(3));
    dif1.selection('input2').set(segmentsString(6));

    dif2 = wp1.geom.create('dif2', 'Difference');
    dif2.selection('input').set(segmentsString(5));
    dif2.selection('input2').set(segmentsString(7));

    uni1 = wp1.geom.create('uni1', 'Union');
    uni1.selection('input').set({segmentsString{2}, segmentsString{4}, 'dif1'});
    uni1.set('intbnd', false);



    if (w_pad>w0) && (l_pad>0)
        
        paramcurve1 = wp1.geom.create('pc1', 'ParametricCurve');
        paramcurve1.set('parmin', 0);
        paramcurve1.set('parmax', 1);
        paramcurve1.set('pos', {'-l_trans-l_pad/2' 'w0*rw1/2'});
        paramcurve1.set('coord', {'s*l_trans' '0.5* (w_pad-w0*rw1)*(-2*s^3+3*s^2)'});
        ls_pad = wp1.geom.create(sprintf('ls%i',6), 'LineSegment');
        ls_pad.set('specify1', 'coord');
        ls_pad.set('coord1', {'-l_pad/2' 'w_pad/2'});
        ls_pad.set('specify2', 'coord');
        ls_pad.set('coord2', {'l_pad/2' 'w_pad/2'});
        paramcurve2 = wp1.geom.create('pc2', 'ParametricCurve');
        paramcurve2.set('parmin', 0);
        paramcurve2.set('parmax', 1);
        paramcurve2.set('pos', [l_trans+l_pad/2 w0*rw1/2]);
        paramcurve2.set('coord', {'-s*l_trans' '0.5* (w_pad-w0*rw1)*(-2*s^3+3*s^2)'});
        mir1 = wp1.geom.create('mir1', 'Mirror');
        mir1.selection('input').set({'pc1', 'pc2', sprintf('ls%i',6)});
        mir1.set('pos', [0 0]);
        mir1.set('axis', [0 1]);
        mir1.set('keep', true);
        
        rec_pad = wp1.geom.create(sprintf('r%i',8), 'Rectangle');
        rec_pad.set('size', [l_pad+2*l_trans w0*rw1]);
        rec_pad.set('base', 'center');
        rec_pad.set('pos', [0 0]);       
        csol1 = wp1.geom.create('csol1', 'ConvertToSolid');
        csol1.selection('input').set({'pc1', 'pc2', sprintf('ls%i',6), 'mir1',sprintf('r%i',8)});

        rotatesString{1} = 'rot1';
        rotates{1} = wp1.geom.create(rotatesString{1}, 'Rotate');
        rotates{1}.selection('input').set('csol1');
        rotates{1}.set('rot', sprintf('%s*180/pi',theta_seg{2}));
        rotates{1}.set('keep', false)
        rotates{1}.set('pos', [0 0])
    
    
        mov1 = wp1.geom.create('mov1', 'Move');
        mov1.selection('input').set(rotatesString{1});
        mov1.set('displ', {x_seg{2} y_seg{2}});
        mov1.set('keep', false);


        if pad_trigger

            paramcurve3 = wp1.geom.create('pc3', 'ParametricCurve');
            paramcurve3.set('parmin', 0);
            paramcurve3.set('parmax', 1);
            paramcurve3.set('pos', {'-l_trans-l_pad/2' 'w0*rw2/2'});
            paramcurve3.set('coord', {'s*l_trans' '0.5* (w_pad-w0*rw2)*(-2*s^3+3*s^2)'});
            ls_pad = wp1.geom.create(sprintf('ls%i',7), 'LineSegment');
            ls_pad.set('specify1', 'coord');
            ls_pad.set('coord1', {'-l_pad/2' 'w_pad/2'});
            ls_pad.set('specify2', 'coord');
            ls_pad.set('coord2', {'l_pad/2' 'w_pad/2'});
            paramcurve4 = wp1.geom.create('pc4', 'ParametricCurve');
            paramcurve4.set('parmin', 0);
            paramcurve4.set('parmax', 1);
            paramcurve4.set('pos', {'l_trans+l_pad/2' 'w0*rw2/2'});
            paramcurve4.set('coord', {'-s*l_trans' '0.5* (w_pad-w0*rw2)*(-2*s^3+3*s^2)'});
            mir2 = wp1.geom.create('mir2', 'Mirror');
            mir2.selection('input').set({'pc3', 'pc4', sprintf('ls%i',7)});
            mir2.set('pos', [0 0]);
            mir2.set('axis', [0 1]);
            mir2.set('keep', true);

            rec_pad = wp1.geom.create(sprintf('r%i',9), 'Rectangle');
            rec_pad.set('size', [l_pad+2*l_trans w0*rw2]);
            rec_pad.set('base', 'center');
            rec_pad.set('pos', [0 0]);       
            csol2 = wp1.geom.create('csol2', 'ConvertToSolid');
            csol2.selection('input').set({'pc3', 'pc4', sprintf('ls%i',7), 'mir2',sprintf('r%i',9)});


            rotatesString{2} = sprintf('rot%i', 2);
            rotates{2} = wp1.geom.create(rotatesString{2}, 'Rotate');
            rotates{2}.selection('input').set('csol2');
            rotates{2}.set('rot', sprintf('%s*180/pi', theta_seg{4}));
            rotates{2}.set('keep', false)
            rotates{2}.set('pos', [0 0])
        
        
            mov2 = wp1.geom.create('mov2', 'Move');
            mov2.selection('input').set(rotatesString{2});
            mov2.set('displ', {x_seg{4} y_seg{4}});
            mov2.set('keep', false);

            uni2 = wp1.geom.create('uni2', 'Union');
            uni2.selection('input').set({'uni1', 'mov1', 'mov2'});
            uni2.set('intbnd', false);
        else
            uni2 = wp1.geom.create('uni2', 'Union');
            uni2.selection('input').set({'uni1', 'mov1'});
            uni2.set('intbnd', false);
        end
    end

    posx = sprintf('(abs((%s + %s * cos(%s) / 2) + (%s - %s * cos(%s) / 2)) / 2)', ...
        x_seg{4}, l_seg{4}, theta_seg{4}, x_seg{2}, l_seg{2}, theta_seg{2});
    
    posy = sprintf('(abs((%s + %s * sin(%s) / 2) + (%s - %s * sin(%s) / 2)) / 2)', ...
        y_seg{4}, l_seg{4}, theta_seg{4}, y_seg{2}, l_seg{2}, theta_seg{2});

    rotatesString{3} = sprintf('rot%i', 3);
    rotates{3} = wp1.geom.create(rotatesString{3}, 'Rotate');
    rotates{3}.selection('input').set('uni2');
    rotates{3}.set('rot', '180');
    rotates{3}.set('keep', true);
    rotates{3}.set('pos', {posx posy});

    
    uni3 = wp1.geom.create('uni3', 'Union');
    uni3.selection('input').set({rotatesString{3}, 'uni2', segmentsString{1}, 'dif2'})
    uni3.set('intbnd', false);

    rotatesString{4} = sprintf('rot%i', 4);
    rotates{4} = wp1.geom.create(rotatesString{4}, 'Rotate');
    rotates{4}.selection('input').set('uni3');
    rotates{4}.set('rot', '180');
    rotates{4}.set('keep', true);
    rotates{4}.set('pos', {'0' '0'});

    uni4 = wp1.geom.create('uni4', 'Union');
    uni4.selection('input').set({rotatesString{4}, 'uni3'})
    uni4.set('intbnd', false);

    uni5 = wp1.geom.create('uni5', 'Union');
    uni5.selection('input').set({clamp_linesString{1:2}})
    uni5.set('intbnd', false);

    rotatesString{5} = sprintf('rot%i', 5);
    rotates{5} = wp1.geom.create(rotatesString{5}, 'Rotate');
    rotates{5}.selection('input').set({'uni5'});
    rotates{5}.set('rot', '180');
    rotates{5}.set('keep', true);
    rotates{5}.set('pos', {posx posy});

    uni6 = wp1.geom.create('uni6', 'Union');
    uni6.selection('input').set({clamp_linesString{3:5}, 'uni5', rotatesString{5}})
    uni6.set('intbnd', false);

    rotatesString{6} = sprintf('rot%i', 6);
    rotates{6} = wp1.geom.create(rotatesString{6}, 'Rotate');
    rotates{6}.selection('input').set({'uni6'});
    rotates{6}.set('rot', '180');
    rotates{6}.set('keep', true);
    rotates{6}.set('pos', {'0' '0'});
    
    geom1.run;
    
    if plot_flag
        mphgeom(model,'geom1', 'vertexlabels', 'off','facemode','off')
        zlim([-10*h_mbr,10*h_mbr])
        view(0,90)
        xlim([-6*l0,6*l0])
        ylim([-4*l0,4*l0])
    end


       %% Selection
    %     hold on
        idx_fixed = [];
        disp('test');
        fprintf('%s-1\n','rotangle');
        
        disp(rotangle);
        
        for i = [3,5]
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
            
            for n = 1:2
                R= ...
                [...
                    cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                    sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                    0, 0, 1 ...
                ];
                box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
                rotatedCoords = (R*box_coordinates);
                idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');
                idx_fixed = [idx_fixed, idx];

                if i == 3
                    box_coordinates = [eval(sprintf('2*(%s) - (%s)', posx, x_select_1)), eval(sprintf('2*(%s) - (%s)', posy, y_select_1)), 0;...
                        eval(sprintf('2*(%s) - (%s)', posx, x_select_2)), eval(sprintf('2*(%s) - (%s)', posy, y_select_2)), 0]';
                    rotatedCoords = (R*box_coordinates);
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');
                    idx_fixed = [idx_fixed, idx];
                end
            end
        end
        sel1 = model.selection.create('sel1').geom(1);
        sel1.set(idx_fixed);
        sel1.label('clamping edges');
        idx_clamps = [];
        idx_clamps_alledges = [];
    
    
        for i = [3,5]
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


            for n = 1:2
                R= ...
                [...
                    cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                    sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                    0, 0, 1 ...
                ];    

                box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
  
                rotatedCoords = (R*box_coordinates);    
                idx = mphselectbox(model, 'geom1',rotatedCoords, 'boundary');    
                idx_clamps = [idx_clamps, idx];
                idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');    
                idx_clamps_alledges = [idx_clamps_alledges, idx];

                if i == 3
                    box_coordinates = [eval(sprintf('2*(%s) - (%s)', posx, x_select_1)), eval(sprintf('2*(%s) - (%s)', posy, y_select_1)), 0;...
                        eval(sprintf('2*(%s) - (%s)', posx, x_select_2)), eval(sprintf('2*(%s) - (%s)', posy, y_select_2)), 0]';
                    rotatedCoords = (R*box_coordinates);
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'boundary');
                    idx_clamps = [idx_clamps, idx];
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');    
                    idx_clamps_alledges = [idx_clamps_alledges, idx];
                end
             end
        end
        sel2 = model.selection.create('sel2').geom(2);
        sel2.set(idx_clamps);
        sel2.label('clamping pads');
        
        idx_clamps_x_inner = [];

        for i = [3,5]

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

            for n = 1:2
                R= ...
                [...
                    cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                    sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                    0, 0, 1 ...
                ];
                box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
                rotatedCoords = (R*box_coordinates);
                idx = mphselectbox(model, 'geom1',rotatedCoords, 'edge');
                idx_clamps_x_inner = [idx_clamps_x_inner, idx];            

                if i == 3
                    box_coordinates = [eval(sprintf('2*(%s) - (%s)', posx, x_select_1)), eval(sprintf('2*(%s) - (%s)', posy, y_select_1)), 0;...
                        eval(sprintf('2*(%s) - (%s)', posx, x_select_2)), eval(sprintf('2*(%s) - (%s)', posy, y_select_2)), 0]';
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

        idx_clamps_support = [];
        for i = [3,5]
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
            x_cut = sprintf('%s', ...
                x_seg{i});
            y_cut = sprintf('%s', ...
                y_seg{i});
            l_select = sprintf('%s + %s + 1e-9', l_seg{i}, l_cut);
            x_select_1 = sprintf('%s - %s / 2', x_cut, l_select);
            y_select_1 = sprintf('%s - %s / 2', y_cut, l_select);
            x_select_2 = sprintf('%s + %s / 2', x_cut, l_select);
            y_select_2 = sprintf('%s + %s / 2', y_cut, l_select);


            for n = 1:2
                R= ...
                [...
                    cos((n-1)*rotangle), -sin((n-1)*rotangle), 0;...
                    sin((n-1)*rotangle), cos((n-1)*rotangle),0;...
                    0, 0, 1 ...
                ];    

                box_coordinates = [eval(x_select_1), eval(y_select_1), 0; eval(x_select_2), eval(y_select_2), 0]';
  
                rotatedCoords = (R*box_coordinates);    
                idx = mphselectbox(model, 'geom1',rotatedCoords, 'boundary');    
                idx_clamps_support = [idx_clamps_support, idx];

                if i == 3
                    box_coordinates = [eval(sprintf('2*(%s) - (%s)', posx, x_select_1)), eval(sprintf('2*(%s) - (%s)', posy, y_select_1)), 0;...
                        eval(sprintf('2*(%s) - (%s)', posx, x_select_2)), eval(sprintf('2*(%s) - (%s)', posy, y_select_2)), 0]';
                    rotatedCoords = (R*box_coordinates);
                    idx = mphselectbox(model, 'geom1',rotatedCoords, 'boundary');
                    idx_clamps_support = [idx_clamps_support, idx];

                end
             end
        end

        idx_all = mphselectbox(model, 'geom1',...
                    [-((2+rl1+rl2)*l0)-lc,-(2*l0+(rl1+rl2)*l0),0;...
                     (2*l0+(rl1+rl2)*l0)+lc,(2*l0+(rl1+rl2)*l0),0]', 'edge');
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

        idx_support = setdiff(idx_clamps_support,idx_clamps);
        sel6 = model.selection.create('sel6').geom(2);
        sel6.set(idx_support);
        sel6.label('support area');

        idx_center_area = mphselectbox(model, 'geom1',...
                    [-((rw1+rw2)*w0+Rady),-((rw1+rw2)*w0+Rady),0;...
                     ((rw1+rw2)*w0+Rady),((rw1+rw2)*w0+Rady),0]', 'boundary');
        sel7 = model.selection.create('sel7').geom(2);
        sel7.set(idx_center_area);
        sel7.label('center area'); 

        idx_polygon_area = setdiff(idx_all_boundaries, [idx_clamps_support, idx_center_area]);
        
        sel8 = model.selection.create('sel8').geom(2);
        sel8.set(idx_polygon_area);
        sel8.label('polygon area'); 



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
    Maxpolygon.selection.named('sel8');

    Maxsupport = model.component('mod1').cpl.create('maxop3', 'Maximum');
    Maxsupport.selection.geom('geom1', 2);
    Maxsupport.selection.named('sel6');
    
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


mphsave(model, 'Practice_Resonator_twin_practice.mph')

%[Freqs, Q ,m_eff, S_F, eta, rl2_match, Q_match] = [1,1,1,1,1,1,1];
%Freqs = 1;
%Q = 1;
%m_eff =1;
%S_F = 1;
%eta = 1;
%rl2_match = 1;
%Q_match = 1;

end
