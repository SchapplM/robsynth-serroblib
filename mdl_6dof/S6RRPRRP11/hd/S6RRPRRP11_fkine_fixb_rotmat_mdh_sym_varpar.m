% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:18:21
% EndTime: 2018-11-23 17:18:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (141->65), mult. (162->60), div. (0->0), fcn. (221->8), ass. (0->40)
t26 = sin(qJ(1));
t28 = cos(qJ(2));
t11 = t26 * t28;
t25 = sin(qJ(2));
t40 = qJ(3) * t25;
t48 = pkin(2) * t11 + t26 * t40;
t30 = -pkin(9) - pkin(8);
t24 = sin(qJ(4));
t47 = t24 * pkin(4);
t27 = cos(qJ(4));
t13 = pkin(4) * t27 + pkin(3);
t46 = t26 * t25;
t45 = t26 * t27;
t23 = qJ(4) + qJ(5);
t14 = sin(t23);
t44 = t28 * t14;
t15 = cos(t23);
t43 = t28 * t15;
t29 = cos(qJ(1));
t42 = t29 * t25;
t41 = t29 * t27;
t12 = t29 * t28;
t22 = pkin(6) + 0;
t39 = pkin(1) * t26 + 0;
t38 = pkin(2) * t25 + t22;
t37 = pkin(1) * t29 + pkin(7) * t26 + 0;
t36 = t39 + t48;
t21 = -qJ(6) + t30;
t6 = pkin(5) * t14 + t47;
t35 = -t21 * t28 + t25 * t6;
t34 = -t29 * pkin(7) + t39;
t33 = pkin(2) * t12 + t29 * t40 + t37;
t32 = t25 * t47 - t28 * t30;
t31 = -qJ(3) * t28 + t38;
t5 = pkin(5) * t15 + t13;
t4 = t14 * t46 - t15 * t29;
t3 = t14 * t29 + t15 * t46;
t2 = t14 * t42 + t15 * t26;
t1 = -t14 * t26 + t15 * t42;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t12, -t42, t26, t37; t11, -t46, -t29, t34; t25, t28, 0, t22; 0, 0, 0, 1; t26, -t12, t42, t33; -t29, -t11, t46, t34 + t48; 0, -t25, -t28, t31; 0, 0, 0, 1; t24 * t42 + t45, -t24 * t26 + t25 * t41, t12, pkin(3) * t26 + pkin(8) * t12 + t33; t24 * t46 - t41, t24 * t29 + t25 * t45, t11, pkin(8) * t11 + (-pkin(3) - pkin(7)) * t29 + t36; -t28 * t24, -t28 * t27, t25, pkin(8) * t25 + t31; 0, 0, 0, 1; t2, t1, t12, t26 * t13 + t29 * t32 + t33; t4, t3, t11 (-pkin(7) - t13) * t29 + t32 * t26 + t36; -t44, -t43, t25, -t25 * t30 + (-qJ(3) - t47) * t28 + t38; 0, 0, 0, 1; t2, t1, t12, t26 * t5 + t29 * t35 + t33; t4, t3, t11 (-pkin(7) - t5) * t29 + t35 * t26 + t36; -t44, -t43, t25, -t25 * t21 + (-qJ(3) - t6) * t28 + t38; 0, 0, 0, 1;];
T_ges = t7;
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
