% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2018-11-23 15:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:56:16
% EndTime: 2018-11-23 15:56:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (190->48), mult. (102->46), div. (0->0), fcn. (152->10), ass. (0->39)
t21 = qJ(3) + pkin(10);
t13 = sin(t21);
t25 = sin(qJ(5));
t44 = t13 * t25;
t22 = qJ(1) + pkin(9);
t14 = sin(t22);
t43 = t14 * t25;
t28 = cos(qJ(5));
t42 = t14 * t28;
t16 = cos(t22);
t41 = t16 * t25;
t40 = t16 * t28;
t39 = pkin(6) + 0;
t27 = sin(qJ(1));
t38 = t27 * pkin(1) + 0;
t30 = cos(qJ(1));
t37 = t30 * pkin(1) + 0;
t29 = cos(qJ(3));
t12 = t29 * pkin(3) + pkin(2);
t36 = t16 * t12 + t37;
t17 = qJ(2) + t39;
t24 = -qJ(4) - pkin(7);
t35 = t14 * t12 + t16 * t24 + t38;
t26 = sin(qJ(3));
t34 = t26 * pkin(3) + t17;
t15 = cos(t21);
t33 = pkin(4) * t15 + pkin(8) * t13;
t11 = t28 * pkin(5) + pkin(4);
t23 = -qJ(6) - pkin(8);
t32 = t11 * t15 - t13 * t23;
t31 = -t14 * t24 + t36;
t10 = t13 * t28;
t8 = t16 * t13;
t7 = t14 * t13;
t4 = t15 * t40 + t43;
t3 = -t15 * t41 + t42;
t2 = t15 * t42 - t41;
t1 = -t15 * t43 - t40;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t39; 0, 0, 0, 1; t16, -t14, 0, t37; t14, t16, 0, t38; 0, 0, 1, t17; 0, 0, 0, 1; t16 * t29, -t16 * t26, t14, t16 * pkin(2) + t14 * pkin(7) + t37; t14 * t29, -t14 * t26, -t16, t14 * pkin(2) - t16 * pkin(7) + t38; t26, t29, 0, t17; 0, 0, 0, 1; t16 * t15, -t8, t14, t31; t14 * t15, -t7, -t16, t35; t13, t15, 0, t34; 0, 0, 0, 1; t4, t3, t8, t16 * t33 + t31; t2, t1, t7, t14 * t33 + t35; t10, -t44, -t15, t13 * pkin(4) - t15 * pkin(8) + t34; 0, 0, 0, 1; t4, t3, t8, t32 * t16 + (pkin(5) * t25 - t24) * t14 + t36; t2, t1, t7, -pkin(5) * t41 + t14 * t32 + t35; t10, -t44, -t15, t13 * t11 + t15 * t23 + t34; 0, 0, 0, 1;];
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
