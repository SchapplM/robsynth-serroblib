% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2018-11-23 17:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:37:30
% EndTime: 2018-11-23 17:37:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (166->68), mult. (286->77), div. (0->0), fcn. (388->10), ass. (0->40)
t30 = sin(qJ(2));
t32 = cos(qJ(3));
t14 = t30 * t32;
t29 = sin(qJ(3));
t50 = t30 * t29;
t52 = pkin(3) * t14 + qJ(4) * t50;
t26 = sin(pkin(10));
t51 = t26 * t29;
t31 = sin(qJ(1));
t15 = t31 * t30;
t33 = cos(qJ(2));
t49 = t31 * t33;
t34 = cos(qJ(1));
t17 = t34 * t30;
t48 = t34 * t33;
t47 = qJ(5) * t30;
t25 = pkin(6) + 0;
t46 = t30 * pkin(2) + t25;
t45 = pkin(5) * t26 + qJ(4);
t44 = t34 * pkin(1) + t31 * pkin(7) + 0;
t43 = t46 + t52;
t42 = t31 * pkin(1) - t34 * pkin(7) + 0;
t41 = pkin(2) * t48 + pkin(8) * t17 + t44;
t40 = -t33 * pkin(8) + t46;
t6 = t31 * t29 + t32 * t48;
t39 = t6 * pkin(3) + t41;
t38 = pkin(2) * t49 + pkin(8) * t15 + t42;
t4 = -t34 * t29 + t32 * t49;
t37 = t4 * pkin(3) + t38;
t5 = t29 * t48 - t31 * t32;
t36 = t5 * qJ(4) + t39;
t3 = t29 * t49 + t34 * t32;
t35 = t3 * qJ(4) + t37;
t28 = -pkin(9) - qJ(5);
t27 = cos(pkin(10));
t24 = pkin(10) + qJ(6);
t19 = cos(t24);
t18 = sin(t24);
t13 = t27 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t48, -t17, t31, t44; t49, -t15, -t34, t42; t30, t33, 0, t25; 0, 0, 0, 1; t6, -t5, t17, t41; t4, -t3, t15, t38; t14, -t50, -t33, t40; 0, 0, 0, 1; t6, t17, t5, t36; t4, t15, t3, t35; t14, -t33, t50, t40 + t52; 0, 0, 0, 1; t5 * t26 + t6 * t27, -t6 * t26 + t5 * t27, -t17, t6 * pkin(4) - t34 * t47 + t36; t3 * t26 + t4 * t27, -t4 * t26 + t3 * t27, -t15, t4 * pkin(4) - t31 * t47 + t35; (t27 * t32 + t51) * t30 (-t26 * t32 + t27 * t29) * t30, t33, pkin(4) * t14 + (-pkin(8) + qJ(5)) * t33 + t43; 0, 0, 0, 1; t5 * t18 + t6 * t19, -t6 * t18 + t5 * t19, -t17, t6 * t13 + t28 * t17 + t45 * t5 + t39; t3 * t18 + t4 * t19, -t4 * t18 + t3 * t19, -t15, t4 * t13 + t28 * t15 + t45 * t3 + t37; (t18 * t29 + t19 * t32) * t30 (-t18 * t32 + t19 * t29) * t30, t33 (-pkin(8) - t28) * t33 + (pkin(5) * t51 + t13 * t32) * t30 + t43; 0, 0, 0, 1;];
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
