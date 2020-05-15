% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 16:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:38:00
% EndTime: 2018-11-23 16:38:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (140->63), mult. (129->66), div. (0->0), fcn. (184->10), ass. (0->37)
t25 = -pkin(9) - pkin(8);
t22 = cos(qJ(4));
t6 = t22 * pkin(4) + pkin(3);
t19 = sin(qJ(4));
t21 = sin(qJ(1));
t41 = t21 * t19;
t20 = sin(qJ(3));
t40 = t21 * t20;
t39 = t21 * t22;
t23 = cos(qJ(3));
t38 = t21 * t23;
t24 = cos(qJ(1));
t37 = t24 * t19;
t36 = t24 * t20;
t35 = t24 * t22;
t18 = qJ(4) + qJ(5);
t16 = pkin(6) + 0;
t34 = t21 * pkin(1) + 0;
t33 = pkin(2) + t16;
t32 = t24 * pkin(1) + t21 * qJ(2) + 0;
t11 = t21 * pkin(7);
t31 = t11 + t34;
t30 = t24 * pkin(7) + t32;
t29 = pkin(3) * t20 - pkin(8) * t23;
t8 = cos(t18);
t1 = pkin(5) * t8 + t6;
t17 = -pkin(10) + t25;
t28 = t1 * t20 + t17 * t23;
t27 = t20 * t6 + t23 * t25;
t26 = -t24 * qJ(2) + t34;
t10 = qJ(6) + t18;
t7 = sin(t18);
t5 = cos(t10);
t4 = sin(t10);
t3 = t24 * t23;
t2 = t19 * pkin(4) + pkin(5) * t7;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t21, 0, 0; t21, t24, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; 0, -t24, t21, t32; 0, -t21, -t24, t26; 1, 0, 0, t16; 0, 0, 0, 1; t40, t38, t24, t30; -t36, -t3, t21, t11 + t26; t23, -t20, 0, t33; 0, 0, 0, 1; t20 * t39 + t37, -t19 * t40 + t35, -t38, t29 * t21 + t30; -t20 * t35 + t41, t19 * t36 + t39, t3 (-qJ(2) - t29) * t24 + t31; t23 * t22, -t23 * t19, t20, t23 * pkin(3) + t20 * pkin(8) + t33; 0, 0, 0, 1; t24 * t7 + t8 * t40, t24 * t8 - t7 * t40, -t38, pkin(4) * t37 + t27 * t21 + t30; t21 * t7 - t8 * t36, t21 * t8 + t7 * t36, t3, pkin(4) * t41 + (-qJ(2) - t27) * t24 + t31; t23 * t8, -t23 * t7, t20, -t20 * t25 + t23 * t6 + t33; 0, 0, 0, 1; t24 * t4 + t5 * t40, t24 * t5 - t4 * t40, -t38, t24 * t2 + t28 * t21 + t30; t21 * t4 - t5 * t36, t21 * t5 + t4 * t36, t3, t21 * t2 + (-qJ(2) - t28) * t24 + t31; t23 * t5, -t23 * t4, t20, t23 * t1 - t20 * t17 + t33; 0, 0, 0, 1;];
T_ges = t9;
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
