% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 15:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:40:13
% EndTime: 2018-11-23 15:40:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (160->43), mult. (76->38), div. (0->0), fcn. (117->10), ass. (0->34)
t14 = qJ(1) + pkin(9);
t6 = sin(t14);
t13 = qJ(4) + pkin(10);
t7 = cos(t13);
t40 = t6 * t7;
t8 = cos(t14);
t39 = t8 * t7;
t16 = sin(qJ(6));
t38 = t6 * t16;
t17 = sin(qJ(4));
t37 = t6 * t17;
t19 = cos(qJ(6));
t36 = t6 * t19;
t35 = t8 * t16;
t34 = t8 * t19;
t33 = pkin(6) + 0;
t18 = sin(qJ(1));
t32 = t18 * pkin(1) + 0;
t21 = cos(qJ(1));
t31 = t21 * pkin(1) + 0;
t30 = t6 * pkin(2) + t32;
t29 = -pkin(4) * t17 - qJ(3);
t9 = qJ(2) + t33;
t28 = t8 * pkin(2) + t6 * qJ(3) + t31;
t27 = pkin(3) + t9;
t5 = sin(t13);
t26 = pkin(5) * t5 - pkin(8) * t7;
t20 = cos(qJ(4));
t25 = t20 * pkin(4) + t27;
t15 = -qJ(5) - pkin(7);
t24 = -t6 * t15 + t30;
t23 = -t8 * qJ(3) + t30;
t22 = pkin(4) * t37 - t8 * t15 + t28;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t18, 0, 0; t18, t21, 0, 0; 0, 0, 1, t33; 0, 0, 0, 1; t8, -t6, 0, t31; t6, t8, 0, t32; 0, 0, 1, t9; 0, 0, 0, 1; 0, -t8, t6, t28; 0, -t6, -t8, t23; 1, 0, 0, t9; 0, 0, 0, 1; t37, t6 * t20, t8, t8 * pkin(7) + t28; -t8 * t17, -t8 * t20, t6, t6 * pkin(7) + t23; t20, -t17, 0, t27; 0, 0, 0, 1; t6 * t5, t40, t8, t22; -t8 * t5, -t39, t6, t29 * t8 + t24; t7, -t5, 0, t25; 0, 0, 0, 1; t5 * t36 + t35, -t5 * t38 + t34, -t40, t26 * t6 + t22; -t5 * t34 + t38, t5 * t35 + t36, t39 (-t26 + t29) * t8 + t24; t7 * t19, -t7 * t16, t5, t7 * pkin(5) + t5 * pkin(8) + t25; 0, 0, 0, 1;];
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
