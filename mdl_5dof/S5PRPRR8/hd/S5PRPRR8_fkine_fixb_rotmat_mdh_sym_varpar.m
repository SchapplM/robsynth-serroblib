% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRPRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:16
% EndTime: 2019-12-05 16:02:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (137->54), mult. (294->61), div. (0->0), fcn. (405->10), ass. (0->38)
t26 = sin(pkin(9));
t27 = sin(pkin(5));
t20 = t26 * t27;
t32 = sin(qJ(2));
t21 = t27 * t32;
t35 = cos(qJ(2));
t50 = t27 * t35;
t28 = cos(pkin(9));
t49 = t28 * t27;
t29 = cos(pkin(5));
t48 = t29 * t32;
t47 = t29 * t35;
t46 = pkin(6) * t49;
t45 = t26 * pkin(1) + 0;
t44 = qJ(1) + 0;
t43 = t28 * pkin(1) + pkin(6) * t20 + 0;
t42 = t29 * pkin(6) + t44;
t8 = t26 * t32 - t28 * t47;
t9 = t26 * t35 + t28 * t48;
t41 = t9 * pkin(2) + t8 * qJ(3) + t45;
t10 = t26 * t47 + t28 * t32;
t11 = -t26 * t48 + t28 * t35;
t40 = t11 * pkin(2) + t10 * qJ(3) + t43;
t39 = pkin(2) * t21 - qJ(3) * t50 + t42;
t38 = t29 * pkin(3) + pkin(7) * t21 + t39;
t37 = pkin(3) * t20 + t11 * pkin(7) + t40;
t36 = t9 * pkin(7) + (-pkin(3) - pkin(6)) * t49 + t41;
t34 = cos(qJ(4));
t33 = cos(qJ(5));
t31 = sin(qJ(4));
t30 = sin(qJ(5));
t13 = t29 * t34 - t31 * t50;
t12 = t29 * t31 + t34 * t50;
t4 = t8 * t31 - t34 * t49;
t3 = t31 * t49 + t8 * t34;
t2 = t10 * t31 + t34 * t20;
t1 = -t10 * t34 + t31 * t20;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t26, 0, 0; t26, t28, 0, 0; 0, 0, 1, t44; 0, 0, 0, 1; t11, -t10, t20, t43; t9, -t8, -t49, t45 - t46; t21, t50, t29, t42; 0, 0, 0, 1; t20, -t11, t10, t40; -t49, -t9, t8, t41 - t46; t29, -t21, -t50, t39; 0, 0, 0, 1; t2, -t1, t11, t37; t4, t3, t9, t36; t13, -t12, t21, t38; 0, 0, 0, 1; t11 * t30 + t2 * t33, t11 * t33 - t2 * t30, t1, t2 * pkin(4) + t1 * pkin(8) + t37; t9 * t30 + t4 * t33, -t4 * t30 + t9 * t33, -t3, t4 * pkin(4) - t3 * pkin(8) + t36; t13 * t33 + t30 * t21, -t13 * t30 + t33 * t21, t12, t13 * pkin(4) + t12 * pkin(8) + t38; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
