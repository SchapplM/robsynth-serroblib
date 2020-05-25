% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:42:26
% EndTime: 2018-11-23 18:42:26
% DurationCPUTime: 0.21s
% Computational Cost: add. (609->86), mult. (635->104), div. (0->0), fcn. (717->18), ass. (0->56)
t38 = sin(qJ(4));
t61 = pkin(4) * t38 + pkin(9);
t45 = -pkin(11) - pkin(10);
t36 = qJ(4) + qJ(5);
t28 = sin(t36);
t69 = pkin(5) * t28 + t61;
t42 = cos(qJ(4));
t27 = pkin(4) * t42 + pkin(3);
t68 = cos(qJ(3));
t37 = sin(pkin(6));
t41 = sin(qJ(1));
t67 = t41 * t37;
t44 = cos(qJ(1));
t66 = t44 * t37;
t65 = cos(pkin(6));
t64 = pkin(6) - qJ(2);
t63 = pkin(6) + qJ(2);
t62 = pkin(7) + 0;
t60 = t37 * t68;
t59 = pkin(8) * t65 + t62;
t58 = pkin(1) * t44 + pkin(8) * t67 + 0;
t57 = cos(t63);
t56 = sin(t64);
t53 = cos(t64) / 0.2e1;
t18 = t53 - t57 / 0.2e1;
t55 = pkin(2) * t18 + t59;
t52 = sin(t63) / 0.2e1;
t17 = t52 - t56 / 0.2e1;
t43 = cos(qJ(2));
t12 = -t17 * t41 + t43 * t44;
t54 = pkin(2) * t12 + t58;
t51 = pkin(1) * t41 - pkin(8) * t66 + 0;
t10 = t17 * t44 + t41 * t43;
t50 = pkin(2) * t10 + t51;
t16 = t52 + t56 / 0.2e1;
t49 = -pkin(9) * t16 + t55;
t40 = sin(qJ(2));
t46 = t53 + t57 / 0.2e1;
t11 = t40 * t44 + t41 * t46;
t48 = pkin(9) * t11 + t54;
t9 = t40 * t41 - t44 * t46;
t47 = pkin(9) * t9 + t50;
t39 = sin(qJ(3));
t35 = -pkin(12) + t45;
t31 = qJ(6) + t36;
t29 = cos(t36);
t26 = cos(t31);
t25 = sin(t31);
t15 = pkin(5) * t29 + t27;
t8 = t18 * t68 + t39 * t65;
t7 = t18 * t39 - t65 * t68;
t4 = t12 * t68 + t39 * t67;
t3 = t12 * t39 - t41 * t60;
t2 = t10 * t68 - t39 * t66;
t1 = t10 * t39 + t44 * t60;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t44, -t41, 0, 0; t41, t44, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t12, -t11, t67, t58; t10, -t9, -t66, t51; t18, t16, t65, t59; 0, 0, 0, 1; t4, -t3, t11, t48; t2, -t1, t9, t47; t8, -t7, -t16, t49; 0, 0, 0, 1; t11 * t38 + t4 * t42, t11 * t42 - t38 * t4, t3, pkin(3) * t4 + pkin(10) * t3 + t48; t2 * t42 + t38 * t9, -t2 * t38 + t42 * t9, t1, pkin(3) * t2 + pkin(10) * t1 + t47; -t16 * t38 + t42 * t8, -t16 * t42 - t38 * t8, t7, pkin(3) * t8 + pkin(10) * t7 + t49; 0, 0, 0, 1; t11 * t28 + t29 * t4, t11 * t29 - t28 * t4, t3, t11 * t61 + t4 * t27 - t3 * t45 + t54; t2 * t29 + t28 * t9, -t2 * t28 + t29 * t9, t1, -t1 * t45 + t2 * t27 + t61 * t9 + t50; -t16 * t28 + t29 * t8, -t16 * t29 - t28 * t8, t7, -t16 * t61 + t8 * t27 - t7 * t45 + t55; 0, 0, 0, 1; t11 * t25 + t26 * t4, t11 * t26 - t25 * t4, t3, t11 * t69 + t4 * t15 - t3 * t35 + t54; t2 * t26 + t25 * t9, -t2 * t25 + t26 * t9, t1, -t1 * t35 + t2 * t15 + t69 * t9 + t50; -t16 * t25 + t26 * t8, -t16 * t26 - t25 * t8, t7, t8 * t15 - t16 * t69 - t7 * t35 + t55; 0, 0, 0, 1;];
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
