% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:23:53
% EndTime: 2018-11-23 16:23:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (203->46), mult. (115->45), div. (0->0), fcn. (169->10), ass. (0->37)
t26 = qJ(1) + pkin(10);
t18 = sin(t26);
t27 = qJ(3) + qJ(4);
t21 = sin(t27);
t11 = t18 * t21;
t22 = cos(t27);
t48 = t18 * t22;
t19 = cos(t26);
t12 = t19 * t21;
t47 = t19 * t22;
t28 = sin(qJ(5));
t46 = t21 * t28;
t45 = t22 * t28;
t31 = cos(qJ(5));
t44 = t22 * t31;
t43 = pkin(6) + 0;
t30 = sin(qJ(1));
t42 = t30 * pkin(1) + 0;
t33 = cos(qJ(1));
t41 = t33 * pkin(1) + 0;
t20 = qJ(2) + t43;
t32 = cos(qJ(3));
t17 = t32 * pkin(3) + pkin(2);
t34 = -pkin(8) - pkin(7);
t40 = t18 * t17 + t19 * t34 + t42;
t29 = sin(qJ(3));
t39 = t29 * pkin(3) + t20;
t38 = pkin(4) * t48 + pkin(9) * t11 + t40;
t37 = t19 * t17 - t18 * t34 + t41;
t36 = pkin(4) * t47 + pkin(9) * t12 + t37;
t35 = t21 * pkin(4) - t22 * pkin(9) + t39;
t15 = t21 * t31;
t4 = t18 * t28 + t19 * t44;
t3 = -t18 * t31 + t19 * t45;
t2 = t18 * t44 - t19 * t28;
t1 = t18 * t45 + t19 * t31;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t30, 0, 0; t30, t33, 0, 0; 0, 0, 1, t43; 0, 0, 0, 1; t19, -t18, 0, t41; t18, t19, 0, t42; 0, 0, 1, t20; 0, 0, 0, 1; t19 * t32, -t19 * t29, t18, t19 * pkin(2) + t18 * pkin(7) + t41; t18 * t32, -t18 * t29, -t19, t18 * pkin(2) - t19 * pkin(7) + t42; t29, t32, 0, t20; 0, 0, 0, 1; t47, -t12, t18, t37; t48, -t11, -t19, t40; t21, t22, 0, t39; 0, 0, 0, 1; t4, -t3, t12, t36; t2, -t1, t11, t38; t15, -t46, -t22, t35; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(5) + t3 * qJ(6) + t36; t2, t11, t1, t2 * pkin(5) + t1 * qJ(6) + t38; t15, -t22, t46 (pkin(5) * t31 + qJ(6) * t28) * t21 + t35; 0, 0, 0, 1;];
T_ges = t5;
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
