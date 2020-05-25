% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:55:35
% EndTime: 2018-11-23 14:55:36
% DurationCPUTime: 0.25s
% Computational Cost: add. (717->88), mult. (534->103), div. (0->0), fcn. (579->22), ass. (0->64)
t47 = pkin(6) - qJ(2);
t32 = cos(t47) / 0.2e1;
t46 = pkin(6) + qJ(2);
t42 = cos(t46);
t81 = t32 - t42 / 0.2e1;
t31 = sin(t46) / 0.2e1;
t39 = sin(t47);
t20 = t31 - t39 / 0.2e1;
t49 = sin(pkin(10));
t50 = sin(pkin(6));
t29 = t49 * t50;
t52 = cos(pkin(10));
t78 = t52 * t50;
t45 = qJ(2) + pkin(11);
t77 = qJ(1) + 0;
t48 = sin(pkin(12));
t76 = pkin(5) * t48 + pkin(8);
t55 = pkin(7) + qJ(3);
t14 = pkin(2) * t20 - t50 * t55;
t59 = cos(qJ(2));
t35 = pkin(2) * t59 + pkin(1);
t75 = t14 * t52 + t35 * t49 + 0;
t74 = pkin(6) - t45;
t73 = pkin(6) + t45;
t65 = sin(t73) / 0.2e1;
t69 = sin(t74);
t17 = t65 - t69 / 0.2e1;
t41 = cos(t45);
t8 = t17 * t52 + t41 * t49;
t72 = pkin(3) * t8 + t75;
t71 = -t49 * t14 + t35 * t52 + 0;
t70 = cos(t73);
t53 = cos(pkin(6));
t68 = pkin(2) * t81 + t53 * t55 + t77;
t10 = -t17 * t49 + t41 * t52;
t67 = pkin(3) * t10 + t71;
t66 = cos(t74) / 0.2e1;
t19 = t66 - t70 / 0.2e1;
t64 = pkin(3) * t19 + t68;
t37 = sin(t45);
t60 = t70 / 0.2e1 + t66;
t7 = t37 * t49 - t52 * t60;
t63 = pkin(8) * t7 + t72;
t9 = t37 * t52 + t49 * t60;
t62 = pkin(8) * t9 + t67;
t18 = t69 / 0.2e1 + t65;
t61 = -pkin(8) * t18 + t64;
t58 = cos(qJ(4));
t57 = sin(qJ(2));
t56 = sin(qJ(4));
t54 = -pkin(9) - qJ(5);
t51 = cos(pkin(12));
t44 = pkin(12) + qJ(6);
t40 = cos(t44);
t36 = sin(t44);
t34 = pkin(5) * t51 + pkin(4);
t21 = t32 + t42 / 0.2e1;
t12 = t19 * t58 + t53 * t56;
t11 = t19 * t56 - t53 * t58;
t4 = t10 * t58 + t29 * t56;
t3 = t10 * t56 - t29 * t58;
t2 = -t56 * t78 + t58 * t8;
t1 = t56 * t8 + t58 * t78;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t52, -t49, 0, 0; t49, t52, 0, 0; 0, 0, 1, t77; 0, 0, 0, 1; -t20 * t49 + t52 * t59, -t21 * t49 - t52 * t57, t29, pkin(1) * t52 + pkin(7) * t29 + 0; t20 * t52 + t49 * t59, t21 * t52 - t49 * t57, -t78, pkin(1) * t49 - pkin(7) * t78 + 0; t81, t31 + t39 / 0.2e1, t53, pkin(7) * t53 + t77; 0, 0, 0, 1; t10, -t9, t29, t71; t8, -t7, -t78, t75; t19, t18, t53, t68; 0, 0, 0, 1; t4, -t3, t9, t62; t2, -t1, t7, t63; t12, -t11, -t18, t61; 0, 0, 0, 1; t4 * t51 + t48 * t9, -t4 * t48 + t51 * t9, t3, pkin(4) * t4 + qJ(5) * t3 + t62; t2 * t51 + t48 * t7, -t2 * t48 + t51 * t7, t1, pkin(4) * t2 + qJ(5) * t1 + t63; t12 * t51 - t18 * t48, -t12 * t48 - t18 * t51, t11, pkin(4) * t12 + qJ(5) * t11 + t61; 0, 0, 0, 1; t36 * t9 + t4 * t40, -t36 * t4 + t40 * t9, t3, -t3 * t54 + t4 * t34 + t76 * t9 + t67; t2 * t40 + t36 * t7, -t2 * t36 + t40 * t7, t1, -t1 * t54 + t2 * t34 + t7 * t76 + t72; t12 * t40 - t18 * t36, -t12 * t36 - t18 * t40, t11, -t11 * t54 + t12 * t34 - t18 * t76 + t64; 0, 0, 0, 1;];
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
