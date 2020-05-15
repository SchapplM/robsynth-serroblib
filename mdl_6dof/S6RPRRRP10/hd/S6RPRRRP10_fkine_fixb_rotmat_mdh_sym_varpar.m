% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2018-11-23 16:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRP10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:29:36
% EndTime: 2018-11-23 16:29:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (136->55), mult. (144->53), div. (0->0), fcn. (203->8), ass. (0->39)
t23 = sin(qJ(4));
t25 = sin(qJ(1));
t48 = t25 * t23;
t24 = sin(qJ(3));
t47 = t25 * t24;
t26 = cos(qJ(4));
t46 = t25 * t26;
t27 = cos(qJ(3));
t45 = t25 * t27;
t22 = qJ(4) + qJ(5);
t14 = sin(t22);
t44 = t27 * t14;
t29 = -pkin(9) - pkin(8);
t43 = t27 * t29;
t28 = cos(qJ(1));
t42 = t28 * t23;
t41 = t28 * t24;
t40 = t28 * t26;
t21 = pkin(6) + 0;
t39 = t25 * pkin(1) + 0;
t38 = pkin(2) + t21;
t17 = t25 * pkin(7);
t37 = t17 + t39;
t36 = t28 * pkin(1) + t25 * qJ(2) + 0;
t35 = t28 * pkin(7) + t36;
t34 = pkin(3) * t24 - pkin(8) * t27;
t33 = -t28 * qJ(2) + t39;
t13 = t26 * pkin(4) + pkin(3);
t32 = t27 * t13 - t24 * t29 + t38;
t31 = pkin(4) * t42 + t13 * t47 + t25 * t43 + t35;
t30 = pkin(4) * t48 + (-t13 * t24 - qJ(2) - t43) * t28 + t37;
t15 = cos(t22);
t12 = t28 * t27;
t8 = t27 * t15;
t4 = t25 * t14 - t15 * t41;
t3 = t14 * t41 + t25 * t15;
t2 = t28 * t14 + t15 * t47;
t1 = t14 * t47 - t28 * t15;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; 0, -t28, t25, t36; 0, -t25, -t28, t33; 1, 0, 0, t21; 0, 0, 0, 1; t47, t45, t28, t35; -t41, -t12, t25, t17 + t33; t27, -t24, 0, t38; 0, 0, 0, 1; t24 * t46 + t42, -t23 * t47 + t40, -t45, t34 * t25 + t35; -t24 * t40 + t48, t23 * t41 + t46, t12 (-qJ(2) - t34) * t28 + t37; t27 * t26, -t27 * t23, t24, t27 * pkin(3) + t24 * pkin(8) + t38; 0, 0, 0, 1; t2, -t1, -t45, t31; t4, t3, t12, t30; t8, -t44, t24, t32; 0, 0, 0, 1; t2, -t45, t1, t2 * pkin(5) + t1 * qJ(6) + t31; t4, t12, -t3, t4 * pkin(5) - t3 * qJ(6) + t30; t8, t24, t44 (pkin(5) * t15 + qJ(6) * t14) * t27 + t32; 0, 0, 0, 1;];
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
