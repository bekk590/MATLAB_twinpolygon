function out = model
%
% variables.m
%
% Model exported on Oct 29 2024, 16:07 by COMSOL 6.2.0.415.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath(['C:\home\bessho\' native2unicode(hex2dec({'75' '59'}), 'unicode')  native2unicode(hex2dec({'5b' '66'}), 'unicode')  native2unicode(hex2dec({'95' 'a2'}), 'unicode')  native2unicode(hex2dec({'90' '23'}), 'unicode') '\Research\20241029_codeproducefromCOMSOL']);

model.param.set('Qint', '2500');
model.param.descr('Qint', 'Intrinsic Q');
model.param.set('stress', '1100000000[Pa]');
model.param.descr('stress', 'Initial stress');
model.param.set('h_mbr', '2.000000e-08[m]');
model.param.descr('h_mbr', 'Thickness of the resonator');
model.param.set('l0', '2.800000e-04[m]');
model.param.descr('l0', 'Length of the fundamental segment');
model.param.set('w0', '2.000000e-07[m]');
model.param.descr('w0', 'Width of the fundamental segment');
model.param.set('N', '4');
model.param.descr('N', 'Number of the tops');
model.param.set('theta', '7.853982e-01');
model.param.descr('theta', 'Branching angle (radians)');
model.param.set('rl1', '2.500000e+00');
model.param.descr('rl1', 'Length contraction ratio left');
model.param.set('rw1', '7.071068e-01');
model.param.descr('rw1', 'Width expansion ratio left');
model.param.set('rl2', '2.658750e+00');
model.param.descr('rl2', 'Length contraction ratio right');
model.param.set('rw2', '7.071068e-01');
model.param.descr('rw2', 'Width expansion ratio right');
model.param.set('l_pad', '2.000000e-05[m]');
model.param.descr('l_pad', 'Length of the pad');
model.param.set('w_pad', '4.300000e-07[m]');
model.param.descr('w_pad', 'Width of the pad');
model.param.set('l_trans', '5.000000e-06[m]');
model.param.descr('l_trans', 'Length of the spline for the pad');
model.param.set('Qint', '12500 * (h_mbr/100e-9[m])');
model.param.set('angle', '2*pi/N');
model.param.set('theta', '(pi-angle)/2');
model.param.set('Diameter', 'sqrt((l0*rl1)^2+(l0*rl2)^2)');
model.param.set('rad', 'Diameter/2');
model.param.set('Radx', 'rad*sin(atan(max(rl1/rl2,rl2/rl1))-pi/4');
model.param.set('Rady', 'rad*cos(atan(max(rl1/rl2,rl2/rl1))-pi/4)');
model.param.set('Radx', 'rad*sin(atan(max(rl1/rl2,rl2/rl1))-pi/4)');
model.param.descr('Diameter', 'Diameter of circumscribed circle');
model.param.descr('rad', 'Radius of circumscribed circle');
model.param.descr('Radx', 'x-coordinate of polygon vertex');
model.param.descr('Rady', 'x-coordinate of polygon vertex');
model.param.descr('Rady', 'y-coordinate of polygon vertex');
model.param.descr('angle', 'Angle of rotation');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('shell', 'Shell', 'geom1');

model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'w0' 'l0'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', {'Radx' 'Rady'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('size', {'l0*rl2/2' 'w0*rw2'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('pos', {'Radx' '0'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', {'Radx' 'Rady + (l0/2)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('pos', {'Radx+ cos(theta)' '0'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r2');

out = model;
