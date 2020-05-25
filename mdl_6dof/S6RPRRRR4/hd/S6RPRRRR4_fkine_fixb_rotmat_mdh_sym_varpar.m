% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:34:09
% EndTime: 2018-11-23 16:34:09
% DurationCPUTime: 0.14s
% Computational Cost: add. (198->53), mult. (89->54), div. (0->0), fcn. (138->12), ass. (0->35)
t28 = sin(qJ(1));
t22 = pkin(11) + qJ(3);
t17 = qJ(4) + t22;
t14 = qJ(5) + t17;
t7 = sin(t14);
t42 = t28 * t7;
t30 = cos(qJ(1));
t41 = t30 * t7;
t25 = cos(pkin(11));
t13 = t25 * pkin(2) + pkin(1);
t27 = sin(qJ(6));
t40 = t28 * t27;
t29 = cos(qJ(6));
t39 = t28 * t29;
t38 = t30 * t27;
t37 = t30 * t29;
t26 = -pkin(7) - qJ(2);
t23 = pkin(6) + 0;
t21 = -pkin(8) + t26;
t16 = cos(t22);
t4 = pkin(3) * t16 + t13;
t18 = -pkin(9) + t21;
t12 = cos(t17);
t3 = pkin(4) * t12 + t4;
t36 = t30 * t18 + t28 * t3 + 0;
t24 = sin(pkin(11));
t35 = t24 * pkin(2) + t23;
t15 = sin(t22);
t34 = pkin(3) * t15 + t35;
t8 = cos(t14);
t33 = pkin(5) * t8 + pkin(10) * t7;
t11 = sin(t17);
t32 = pkin(4) * t11 + t34;
t31 = -t28 * t18 + t30 * t3 + 0;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t28, 0, 0; t28, t30, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t30 * t25, -t30 * t24, t28, t30 * pkin(1) + t28 * qJ(2) + 0; t28 * t25, -t28 * t24, -t30, t28 * pkin(1) - t30 * qJ(2) + 0; t24, t25, 0, t23; 0, 0, 0, 1; t30 * t16, -t30 * t15, t28, t30 * t13 - t28 * t26 + 0; t28 * t16, -t28 * t15, -t30, t28 * t13 + t30 * t26 + 0; t15, t16, 0, t35; 0, 0, 0, 1; t30 * t12, -t30 * t11, t28, -t28 * t21 + t30 * t4 + 0; t28 * t12, -t28 * t11, -t30, t30 * t21 + t28 * t4 + 0; t11, t12, 0, t34; 0, 0, 0, 1; t30 * t8, -t41, t28, t31; t28 * t8, -t42, -t30, t36; t7, t8, 0, t32; 0, 0, 0, 1; t8 * t37 + t40, -t8 * t38 + t39, t41, t33 * t30 + t31; t8 * t39 - t38, -t8 * t40 - t37, t42, t33 * t28 + t36; t7 * t29, -t7 * t27, -t8, t7 * pkin(5) - t8 * pkin(10) + t32; 0, 0, 0, 1;];
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
