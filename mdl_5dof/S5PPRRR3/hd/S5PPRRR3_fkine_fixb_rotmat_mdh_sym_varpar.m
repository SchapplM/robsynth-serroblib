% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PPRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:14
% EndTime: 2019-12-05 15:16:14
% DurationCPUTime: 0.16s
% Computational Cost: add. (114->54), mult. (191->67), div. (0->0), fcn. (266->10), ass. (0->37)
t20 = sin(pkin(9));
t27 = cos(qJ(3));
t45 = t20 * t27;
t21 = sin(pkin(8));
t44 = t21 * t20;
t22 = cos(pkin(9));
t43 = t21 * t22;
t25 = sin(qJ(3));
t42 = t21 * t25;
t41 = t21 * t27;
t23 = cos(pkin(8));
t40 = t23 * t20;
t39 = t23 * t22;
t38 = t23 * t25;
t37 = t23 * t27;
t24 = sin(qJ(4));
t36 = t24 * t44;
t35 = t24 * t40;
t18 = qJ(1) + 0;
t34 = t23 * pkin(1) + t21 * qJ(2) + 0;
t33 = t20 * pkin(2) + t18;
t32 = pkin(2) * t39 + pkin(5) * t40 + t34;
t31 = t21 * pkin(1) - t23 * qJ(2) + 0;
t30 = -t22 * pkin(5) + t33;
t29 = pkin(2) * t43 + pkin(5) * t44 + t31;
t28 = -pkin(7) - pkin(6);
t26 = cos(qJ(4));
t19 = qJ(4) + qJ(5);
t14 = cos(t19);
t13 = sin(t19);
t11 = t26 * pkin(4) + pkin(3);
t9 = t20 * t25;
t4 = t22 * t37 + t42;
t3 = t22 * t38 - t41;
t2 = t22 * t41 - t38;
t1 = t22 * t42 + t37;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t21, 0, 0; t21, t23, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t39, -t40, t21, t34; t43, -t44, -t23, t31; t20, t22, 0, t18; 0, 0, 0, 1; t4, -t3, t40, t32; t2, -t1, t44, t29; t45, -t9, -t22, t30; 0, 0, 0, 1; t4 * t26 + t35, -t4 * t24 + t26 * t40, t3, t4 * pkin(3) + t3 * pkin(6) + t32; t2 * t26 + t36, -t2 * t24 + t26 * t44, t1, t2 * pkin(3) + t1 * pkin(6) + t29; -t22 * t24 + t26 * t45, -t22 * t26 - t24 * t45, t9, (pkin(3) * t27 + pkin(6) * t25) * t20 + t30; 0, 0, 0, 1; t13 * t40 + t4 * t14, -t4 * t13 + t14 * t40, t3, pkin(4) * t35 + t4 * t11 - t3 * t28 + t32; t13 * t44 + t2 * t14, -t2 * t13 + t14 * t44, t1, pkin(4) * t36 - t1 * t28 + t2 * t11 + t29; -t22 * t13 + t14 * t45, -t13 * t45 - t22 * t14, t9, (-pkin(4) * t24 - pkin(5)) * t22 + (t11 * t27 - t25 * t28) * t20 + t33; 0, 0, 0, 1;];
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
