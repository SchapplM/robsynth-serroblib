% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 15:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:48:35
% EndTime: 2018-11-23 15:48:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (162->50), mult. (95->51), div. (0->0), fcn. (141->10), ass. (0->36)
t17 = sin(qJ(5));
t15 = qJ(1) + pkin(10);
t8 = sin(t15);
t41 = t8 * t17;
t21 = cos(qJ(4));
t40 = t8 * t21;
t9 = cos(t15);
t39 = t9 * t17;
t16 = qJ(5) + qJ(6);
t11 = sin(t16);
t18 = sin(qJ(4));
t38 = t11 * t18;
t12 = cos(t16);
t37 = t12 * t18;
t36 = t17 * t18;
t20 = cos(qJ(5));
t35 = t18 * t20;
t34 = pkin(6) + 0;
t19 = sin(qJ(1));
t33 = t19 * pkin(1) + 0;
t22 = cos(qJ(1));
t32 = t22 * pkin(1) + 0;
t31 = t8 * pkin(2) + t33;
t10 = qJ(2) + t34;
t3 = t8 * pkin(7);
t30 = t3 + t31;
t29 = t9 * pkin(2) + t8 * qJ(3) + t32;
t28 = pkin(3) + t10;
t27 = pkin(4) * t18 - pkin(8) * t21;
t23 = -pkin(9) - pkin(8);
t7 = t20 * pkin(5) + pkin(4);
t26 = t18 * t7 + t21 * t23;
t25 = t9 * pkin(7) + t29;
t24 = -t9 * qJ(3) + t31;
t1 = t9 * t21;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t34; 0, 0, 0, 1; t9, -t8, 0, t32; t8, t9, 0, t33; 0, 0, 1, t10; 0, 0, 0, 1; 0, -t9, t8, t29; 0, -t8, -t9, t24; 1, 0, 0, t10; 0, 0, 0, 1; t8 * t18, t40, t9, t25; -t9 * t18, -t1, t8, t24 + t3; t21, -t18, 0, t28; 0, 0, 0, 1; t8 * t35 + t39, t9 * t20 - t8 * t36, -t40, t27 * t8 + t25; -t9 * t35 + t41, t8 * t20 + t9 * t36, t1 (-qJ(3) - t27) * t9 + t30; t21 * t20, -t21 * t17, t18, t21 * pkin(4) + t18 * pkin(8) + t28; 0, 0, 0, 1; t9 * t11 + t8 * t37, t9 * t12 - t8 * t38, -t40, pkin(5) * t39 + t26 * t8 + t25; t8 * t11 - t9 * t37, t8 * t12 + t9 * t38, t1, pkin(5) * t41 + (-qJ(3) - t26) * t9 + t30; t21 * t12, -t21 * t11, t18, -t18 * t23 + t21 * t7 + t28; 0, 0, 0, 1;];
T_ges = t2;
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
