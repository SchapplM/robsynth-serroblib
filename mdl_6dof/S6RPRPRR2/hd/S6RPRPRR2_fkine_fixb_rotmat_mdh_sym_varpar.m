% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:17
% EndTime: 2018-11-23 16:03:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->53), mult. (102->56), div. (0->0), fcn. (152->12), ass. (0->40)
t20 = qJ(5) + qJ(6);
t13 = sin(t20);
t19 = qJ(1) + pkin(10);
t9 = sin(t19);
t45 = t9 * t13;
t14 = cos(t20);
t44 = t9 * t14;
t22 = sin(qJ(5));
t43 = t9 * t22;
t25 = cos(qJ(5));
t42 = t9 * t25;
t11 = cos(t19);
t41 = t11 * t13;
t40 = t11 * t14;
t39 = t11 * t22;
t38 = t11 * t25;
t37 = pkin(6) + 0;
t24 = sin(qJ(1));
t36 = pkin(1) * t24 + 0;
t27 = cos(qJ(1));
t35 = pkin(1) * t27 + 0;
t26 = cos(qJ(3));
t7 = pkin(3) * t26 + pkin(2);
t34 = t11 * t7 + t35;
t12 = qJ(2) + t37;
t21 = -qJ(4) - pkin(7);
t33 = t11 * t21 + t7 * t9 + t36;
t23 = sin(qJ(3));
t32 = pkin(3) * t23 + t12;
t18 = qJ(3) + pkin(11);
t10 = cos(t18);
t8 = sin(t18);
t31 = pkin(4) * t10 + pkin(8) * t8;
t28 = -pkin(9) - pkin(8);
t6 = pkin(5) * t25 + pkin(4);
t30 = t10 * t6 - t28 * t8;
t29 = -t9 * t21 + t34;
t4 = t11 * t8;
t3 = t9 * t8;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t24, 0, 0; t24, t27, 0, 0; 0, 0, 1, t37; 0, 0, 0, 1; t11, -t9, 0, t35; t9, t11, 0, t36; 0, 0, 1, t12; 0, 0, 0, 1; t11 * t26, -t11 * t23, t9, pkin(2) * t11 + pkin(7) * t9 + t35; t9 * t26, -t9 * t23, -t11, pkin(2) * t9 - pkin(7) * t11 + t36; t23, t26, 0, t12; 0, 0, 0, 1; t11 * t10, -t4, t9, t29; t9 * t10, -t3, -t11, t33; t8, t10, 0, t32; 0, 0, 0, 1; t10 * t38 + t43, -t10 * t39 + t42, t4, t11 * t31 + t29; t10 * t42 - t39, -t10 * t43 - t38, t3, t31 * t9 + t33; t8 * t25, -t8 * t22, -t10, pkin(4) * t8 - pkin(8) * t10 + t32; 0, 0, 0, 1; t10 * t40 + t45, -t10 * t41 + t44, t4 (pkin(5) * t22 - t21) * t9 + t30 * t11 + t34; t10 * t44 - t41, -t10 * t45 - t40, t3, -pkin(5) * t39 + t30 * t9 + t33; t8 * t14, -t8 * t13, -t10, t10 * t28 + t6 * t8 + t32; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
