% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:53:16
% EndTime: 2018-11-23 17:53:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (173->62), mult. (130->60), div. (0->0), fcn. (184->10), ass. (0->40)
t26 = cos(qJ(1));
t20 = qJ(2) + qJ(3);
t14 = sin(t20);
t38 = qJ(4) * t14;
t16 = cos(t20);
t8 = t26 * t16;
t48 = pkin(3) * t8 + t26 * t38;
t21 = sin(qJ(5));
t47 = pkin(5) * t21;
t23 = sin(qJ(1));
t46 = t23 * t14;
t19 = qJ(5) + qJ(6);
t15 = cos(t19);
t45 = t23 * t15;
t7 = t23 * t16;
t44 = t23 * t21;
t24 = cos(qJ(5));
t43 = t23 * t24;
t42 = t26 * t14;
t41 = t26 * t15;
t40 = t26 * t21;
t39 = t26 * t24;
t18 = pkin(6) + 0;
t25 = cos(qJ(2));
t11 = t25 * pkin(2) + pkin(1);
t37 = t26 * t11 + 0;
t22 = sin(qJ(2));
t36 = t22 * pkin(2) + t18;
t28 = -pkin(8) - pkin(7);
t35 = t23 * t11 + t26 * t28 + 0;
t34 = t37 + t48;
t33 = t14 * pkin(3) + t36;
t32 = pkin(3) * t7 + t23 * t38 + t35;
t31 = -t23 * t28 + t37;
t27 = -pkin(10) - pkin(9);
t30 = t14 * t47 - t16 * t27;
t29 = -t16 * qJ(4) + t33;
t13 = sin(t19);
t10 = t24 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t23, 0, 0; t23, t26, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t26 * t25, -t26 * t22, t23, t26 * pkin(1) + t23 * pkin(7) + 0; t23 * t25, -t23 * t22, -t26, t23 * pkin(1) - t26 * pkin(7) + 0; t22, t25, 0, t18; 0, 0, 0, 1; t8, -t42, t23, t31; t7, -t46, -t26, t35; t14, t16, 0, t36; 0, 0, 0, 1; t23, -t8, t42, t31 + t48; -t26, -t7, t46, t32; 0, -t14, -t16, t29; 0, 0, 0, 1; t14 * t40 + t43, t14 * t39 - t44, t8, pkin(9) * t8 + (pkin(4) - t28) * t23 + t34; t14 * t44 - t39, t14 * t43 + t40, t7, -t26 * pkin(4) + pkin(9) * t7 + t32; -t16 * t21, -t16 * t24, t14, t14 * pkin(9) + t29; 0, 0, 0, 1; t13 * t42 + t45, -t23 * t13 + t14 * t41, t8, t30 * t26 + (t10 - t28) * t23 + t34; t13 * t46 - t41, t26 * t13 + t14 * t45, t7, -t26 * t10 + t30 * t23 + t32; -t16 * t13, -t16 * t15, t14, -t14 * t27 + (-qJ(4) - t47) * t16 + t33; 0, 0, 0, 1;];
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
