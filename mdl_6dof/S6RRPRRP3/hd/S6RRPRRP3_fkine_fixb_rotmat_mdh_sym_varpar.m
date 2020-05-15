% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:11:36
% EndTime: 2018-11-23 17:11:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (189->61), mult. (137->64), div. (0->0), fcn. (196->10), ass. (0->45)
t32 = -pkin(9) - pkin(8);
t26 = sin(qJ(4));
t49 = t26 * pkin(4);
t29 = cos(qJ(4));
t13 = t29 * pkin(4) + pkin(3);
t22 = qJ(2) + pkin(10);
t15 = sin(t22);
t24 = qJ(4) + qJ(5);
t17 = sin(t24);
t48 = t15 * t17;
t28 = sin(qJ(1));
t47 = t28 * t17;
t18 = cos(t24);
t46 = t28 * t18;
t45 = t28 * t26;
t44 = t28 * t29;
t31 = cos(qJ(1));
t43 = t31 * t17;
t42 = t31 * t18;
t41 = t31 * t26;
t40 = t31 * t29;
t23 = pkin(6) + 0;
t30 = cos(qJ(2));
t14 = t30 * pkin(2) + pkin(1);
t39 = t31 * t14 + 0;
t27 = sin(qJ(2));
t38 = t27 * pkin(2) + t23;
t25 = -qJ(3) - pkin(7);
t37 = t28 * t14 + t31 * t25 + 0;
t16 = cos(t22);
t36 = pkin(3) * t16 + pkin(8) * t15;
t21 = -qJ(6) + t32;
t5 = pkin(5) * t18 + t13;
t35 = -t15 * t21 + t16 * t5;
t34 = t13 * t16 - t15 * t32;
t33 = -t28 * t25 + t39;
t11 = t31 * t15;
t10 = t28 * t15;
t7 = t15 * t18;
t6 = pkin(5) * t17 + t49;
t4 = t16 * t42 + t47;
t3 = -t16 * t43 + t46;
t2 = t16 * t46 - t43;
t1 = -t16 * t47 - t42;
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t28, 0, 0; t28, t31, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t31 * t30, -t31 * t27, t28, t31 * pkin(1) + t28 * pkin(7) + 0; t28 * t30, -t28 * t27, -t31, t28 * pkin(1) - t31 * pkin(7) + 0; t27, t30, 0, t23; 0, 0, 0, 1; t31 * t16, -t11, t28, t33; t28 * t16, -t10, -t31, t37; t15, t16, 0, t38; 0, 0, 0, 1; t16 * t40 + t45, -t16 * t41 + t44, t11, t36 * t31 + t33; t16 * t44 - t41, -t16 * t45 - t40, t10, t36 * t28 + t37; t15 * t29, -t15 * t26, -t16, t15 * pkin(3) - t16 * pkin(8) + t38; 0, 0, 0, 1; t4, t3, t11, t34 * t31 + (-t25 + t49) * t28 + t39; t2, t1, t10, -pkin(4) * t41 + t34 * t28 + t37; t7, -t48, -t16, t15 * t13 + t16 * t32 + t38; 0, 0, 0, 1; t4, t3, t11, t35 * t31 + (-t25 + t6) * t28 + t39; t2, t1, t10, t35 * t28 - t31 * t6 + t37; t7, -t48, -t16, t15 * t5 + t16 * t21 + t38; 0, 0, 0, 1;];
T_ges = t8;
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
