% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:02:35
% EndTime: 2018-11-23 17:02:35
% DurationCPUTime: 0.23s
% Computational Cost: add. (717->88), mult. (534->103), div. (0->0), fcn. (579->22), ass. (0->64)
t47 = pkin(6) - qJ(2);
t31 = cos(t47) / 0.2e1;
t46 = pkin(6) + qJ(2);
t42 = cos(t46);
t81 = t31 - t42 / 0.2e1;
t30 = sin(t46) / 0.2e1;
t39 = sin(t47);
t20 = t30 - t39 / 0.2e1;
t49 = sin(pkin(6));
t56 = sin(qJ(1));
t32 = t56 * t49;
t59 = cos(qJ(1));
t78 = t59 * t49;
t77 = pkin(7) + 0;
t45 = qJ(2) + pkin(11);
t48 = sin(pkin(12));
t76 = pkin(5) * t48 + pkin(9);
t53 = pkin(8) + qJ(3);
t14 = t20 * pkin(2) - t49 * t53;
t58 = cos(qJ(2));
t35 = t58 * pkin(2) + pkin(1);
t75 = t59 * t14 + t56 * t35 + 0;
t74 = pkin(6) - t45;
t73 = pkin(6) + t45;
t64 = sin(t73) / 0.2e1;
t69 = sin(t74);
t17 = t64 - t69 / 0.2e1;
t41 = cos(t45);
t8 = t59 * t17 + t56 * t41;
t72 = t8 * pkin(3) + t75;
t71 = -t56 * t14 + t59 * t35 + 0;
t70 = cos(t73);
t51 = cos(pkin(6));
t68 = t81 * pkin(2) + t51 * t53 + t77;
t10 = -t56 * t17 + t59 * t41;
t67 = t10 * pkin(3) + t71;
t65 = cos(t74) / 0.2e1;
t19 = t65 - t70 / 0.2e1;
t66 = t19 * pkin(3) + t68;
t37 = sin(t45);
t60 = t70 / 0.2e1 + t65;
t7 = t56 * t37 - t59 * t60;
t63 = t7 * pkin(9) + t72;
t9 = t59 * t37 + t56 * t60;
t62 = t9 * pkin(9) + t67;
t18 = t69 / 0.2e1 + t64;
t61 = -t18 * pkin(9) + t66;
t57 = cos(qJ(4));
t55 = sin(qJ(2));
t54 = sin(qJ(4));
t52 = -pkin(10) - qJ(5);
t50 = cos(pkin(12));
t44 = pkin(12) + qJ(6);
t40 = cos(t44);
t36 = sin(t44);
t34 = t50 * pkin(5) + pkin(4);
t21 = t31 + t42 / 0.2e1;
t12 = t19 * t57 + t51 * t54;
t11 = t19 * t54 - t51 * t57;
t4 = t10 * t57 + t54 * t32;
t3 = t10 * t54 - t57 * t32;
t2 = -t54 * t78 + t8 * t57;
t1 = t8 * t54 + t57 * t78;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t59, -t56, 0, 0; t56, t59, 0, 0; 0, 0, 1, t77; 0, 0, 0, 1; -t56 * t20 + t59 * t58, -t56 * t21 - t59 * t55, t32, t59 * pkin(1) + pkin(8) * t32 + 0; t59 * t20 + t56 * t58, t59 * t21 - t56 * t55, -t78, t56 * pkin(1) - pkin(8) * t78 + 0; t81, t30 + t39 / 0.2e1, t51, t51 * pkin(8) + t77; 0, 0, 0, 1; t10, -t9, t32, t71; t8, -t7, -t78, t75; t19, t18, t51, t68; 0, 0, 0, 1; t4, -t3, t9, t62; t2, -t1, t7, t63; t12, -t11, -t18, t61; 0, 0, 0, 1; t4 * t50 + t9 * t48, -t4 * t48 + t9 * t50, t3, t4 * pkin(4) + t3 * qJ(5) + t62; t2 * t50 + t7 * t48, -t2 * t48 + t7 * t50, t1, t2 * pkin(4) + t1 * qJ(5) + t63; t12 * t50 - t18 * t48, -t12 * t48 - t18 * t50, t11, t12 * pkin(4) + t11 * qJ(5) + t61; 0, 0, 0, 1; t9 * t36 + t4 * t40, -t4 * t36 + t9 * t40, t3, -t3 * t52 + t4 * t34 + t76 * t9 + t67; t2 * t40 + t7 * t36, -t2 * t36 + t7 * t40, t1, -t1 * t52 + t2 * t34 + t76 * t7 + t72; t12 * t40 - t18 * t36, -t12 * t36 - t18 * t40, t11, -t11 * t52 + t12 * t34 - t76 * t18 + t66; 0, 0, 0, 1;];
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
