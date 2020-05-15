% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2018-11-23 15:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:52:05
% EndTime: 2018-11-23 15:52:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (200->53), mult. (102->56), div. (0->0), fcn. (152->12), ass. (0->38)
t20 = qJ(1) + pkin(9);
t10 = sin(t20);
t19 = qJ(3) + pkin(10);
t12 = cos(t19);
t43 = t10 * t12;
t21 = sin(pkin(11));
t42 = t10 * t21;
t22 = cos(pkin(11));
t41 = t10 * t22;
t13 = cos(t20);
t40 = t13 * t12;
t39 = t13 * t21;
t38 = t13 * t22;
t37 = pkin(6) + 0;
t26 = sin(qJ(1));
t36 = t26 * pkin(1) + 0;
t28 = cos(qJ(1));
t35 = t28 * pkin(1) + 0;
t27 = cos(qJ(3));
t7 = t27 * pkin(3) + pkin(2);
t34 = t13 * t7 + t35;
t14 = qJ(2) + t37;
t23 = -qJ(4) - pkin(7);
t33 = t10 * t7 + t13 * t23 + t36;
t25 = sin(qJ(3));
t32 = t25 * pkin(3) + t14;
t24 = -pkin(8) - qJ(5);
t6 = t22 * pkin(5) + pkin(4);
t9 = sin(t19);
t31 = t12 * t6 - t24 * t9;
t30 = pkin(4) * t12 + qJ(5) * t9;
t29 = -t10 * t23 + t34;
t18 = pkin(11) + qJ(6);
t11 = cos(t18);
t8 = sin(t18);
t4 = t13 * t9;
t3 = t10 * t9;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t26, 0, 0; t26, t28, 0, 0; 0, 0, 1, t37; 0, 0, 0, 1; t13, -t10, 0, t35; t10, t13, 0, t36; 0, 0, 1, t14; 0, 0, 0, 1; t13 * t27, -t13 * t25, t10, t13 * pkin(2) + t10 * pkin(7) + t35; t10 * t27, -t10 * t25, -t13, t10 * pkin(2) - t13 * pkin(7) + t36; t25, t27, 0, t14; 0, 0, 0, 1; t40, -t4, t10, t29; t43, -t3, -t13, t33; t9, t12, 0, t32; 0, 0, 0, 1; t12 * t38 + t42, -t12 * t39 + t41, t4, t30 * t13 + t29; t12 * t41 - t39, -t12 * t42 - t38, t3, t30 * t10 + t33; t9 * t22, -t9 * t21, -t12, t9 * pkin(4) - t12 * qJ(5) + t32; 0, 0, 0, 1; t10 * t8 + t11 * t40, t10 * t11 - t8 * t40, t4, t31 * t13 + (pkin(5) * t21 - t23) * t10 + t34; t11 * t43 - t13 * t8, -t13 * t11 - t8 * t43, t3, -pkin(5) * t39 + t31 * t10 + t33; t9 * t11, -t9 * t8, -t12, t12 * t24 + t9 * t6 + t32; 0, 0, 0, 1;];
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
