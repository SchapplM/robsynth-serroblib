% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRRPR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:20:19
% EndTime: 2018-11-23 18:20:19
% DurationCPUTime: 0.23s
% Computational Cost: add. (609->86), mult. (635->104), div. (0->0), fcn. (717->18), ass. (0->56)
t39 = sin(qJ(4));
t61 = t39 * pkin(4) + pkin(9);
t36 = qJ(4) + pkin(12);
t28 = sin(t36);
t69 = pkin(5) * t28 + t61;
t43 = cos(qJ(4));
t27 = t43 * pkin(4) + pkin(3);
t68 = cos(qJ(3));
t37 = sin(pkin(6));
t42 = sin(qJ(1));
t67 = t42 * t37;
t45 = cos(qJ(1));
t66 = t45 * t37;
t38 = -qJ(5) - pkin(10);
t65 = cos(pkin(6));
t64 = pkin(6) - qJ(2);
t63 = pkin(6) + qJ(2);
t62 = pkin(7) + 0;
t60 = t37 * t68;
t59 = t65 * pkin(8) + t62;
t58 = t45 * pkin(1) + pkin(8) * t67 + 0;
t57 = cos(t63);
t56 = sin(t64);
t53 = cos(t64) / 0.2e1;
t18 = t53 - t57 / 0.2e1;
t55 = t18 * pkin(2) + t59;
t52 = sin(t63) / 0.2e1;
t17 = t52 - t56 / 0.2e1;
t44 = cos(qJ(2));
t12 = -t42 * t17 + t45 * t44;
t54 = t12 * pkin(2) + t58;
t51 = t42 * pkin(1) - pkin(8) * t66 + 0;
t10 = t45 * t17 + t42 * t44;
t50 = t10 * pkin(2) + t51;
t16 = t52 + t56 / 0.2e1;
t49 = -t16 * pkin(9) + t55;
t41 = sin(qJ(2));
t46 = t53 + t57 / 0.2e1;
t11 = t45 * t41 + t42 * t46;
t48 = t11 * pkin(9) + t54;
t9 = t42 * t41 - t45 * t46;
t47 = t9 * pkin(9) + t50;
t40 = sin(qJ(3));
t35 = -pkin(11) + t38;
t30 = qJ(6) + t36;
t29 = cos(t36);
t26 = cos(t30);
t25 = sin(t30);
t15 = pkin(5) * t29 + t27;
t8 = t18 * t68 + t65 * t40;
t7 = t18 * t40 - t65 * t68;
t4 = t12 * t68 + t40 * t67;
t3 = t12 * t40 - t42 * t60;
t2 = t10 * t68 - t40 * t66;
t1 = t10 * t40 + t45 * t60;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t45, -t42, 0, 0; t42, t45, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t12, -t11, t67, t58; t10, -t9, -t66, t51; t18, t16, t65, t59; 0, 0, 0, 1; t4, -t3, t11, t48; t2, -t1, t9, t47; t8, -t7, -t16, t49; 0, 0, 0, 1; t11 * t39 + t4 * t43, t11 * t43 - t4 * t39, t3, t4 * pkin(3) + t3 * pkin(10) + t48; t2 * t43 + t9 * t39, -t2 * t39 + t9 * t43, t1, t2 * pkin(3) + t1 * pkin(10) + t47; -t16 * t39 + t8 * t43, -t16 * t43 - t8 * t39, t7, t8 * pkin(3) + t7 * pkin(10) + t49; 0, 0, 0, 1; t11 * t28 + t4 * t29, t11 * t29 - t4 * t28, t3, t61 * t11 + t4 * t27 - t3 * t38 + t54; t2 * t29 + t9 * t28, -t2 * t28 + t9 * t29, t1, -t1 * t38 + t2 * t27 + t61 * t9 + t50; -t16 * t28 + t8 * t29, -t16 * t29 - t8 * t28, t7, -t61 * t16 + t8 * t27 - t7 * t38 + t55; 0, 0, 0, 1; t11 * t25 + t4 * t26, t11 * t26 - t4 * t25, t3, t69 * t11 + t4 * t15 - t3 * t35 + t54; t2 * t26 + t9 * t25, -t2 * t25 + t9 * t26, t1, -t1 * t35 + t2 * t15 + t69 * t9 + t50; -t16 * t25 + t8 * t26, -t16 * t26 - t8 * t25, t7, t8 * t15 - t69 * t16 - t7 * t35 + t55; 0, 0, 0, 1;];
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
