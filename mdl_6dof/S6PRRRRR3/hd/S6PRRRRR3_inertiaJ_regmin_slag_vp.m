% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t62 = sin(qJ(5));
t63 = sin(qJ(4));
t67 = cos(qJ(5));
t68 = cos(qJ(4));
t40 = t62 * t63 - t67 * t68;
t51 = -t68 * pkin(4) - pkin(3);
t31 = t40 * pkin(5) + t51;
t94 = 0.2e1 * t31;
t93 = 0.2e1 * t51;
t64 = sin(qJ(3));
t92 = -0.2e1 * t64;
t69 = cos(qJ(3));
t91 = -0.2e1 * t69;
t90 = 0.2e1 * t69;
t89 = pkin(9) + pkin(10);
t88 = pkin(3) * t68;
t87 = pkin(8) * t63;
t61 = sin(qJ(6));
t86 = t61 * pkin(5);
t85 = t62 * pkin(4);
t66 = cos(qJ(6));
t44 = -t69 * pkin(3) - t64 * pkin(9) - pkin(2);
t39 = t68 * t44;
t75 = t68 * t64;
t21 = -pkin(10) * t75 + t39 + (-pkin(4) - t87) * t69;
t74 = t68 * t69;
t71 = pkin(8) * t74;
t25 = t71 + (-pkin(10) * t64 + t44) * t63;
t76 = t67 * t25;
t11 = t62 * t21 + t76;
t41 = t62 * t68 + t67 * t63;
t32 = t41 * t64;
t9 = -t32 * pkin(11) + t11;
t84 = t66 * t9;
t83 = t69 * pkin(4);
t82 = t69 * pkin(5);
t59 = sin(pkin(6));
t81 = t59 * sin(qJ(2));
t80 = t59 * cos(qJ(2));
t79 = t63 * t64;
t78 = t63 * t68;
t77 = t63 * t69;
t52 = t64 * pkin(8);
t43 = pkin(4) * t79 + t52;
t73 = t64 * t90;
t72 = t66 * t85;
t10 = t67 * t21 - t62 * t25;
t33 = t40 * t64;
t6 = t33 * pkin(11) + t10 - t82;
t1 = t66 * t6 - t61 * t9;
t45 = t89 * t63;
t46 = t89 * t68;
t26 = -t67 * t45 - t62 * t46;
t54 = t67 * pkin(4);
t50 = t54 + pkin(5);
t34 = t66 * t50 - t61 * t85;
t2 = t61 * t6 + t84;
t27 = -t62 * t45 + t67 * t46;
t60 = cos(pkin(6));
t58 = t69 ^ 2;
t57 = t68 ^ 2;
t56 = t64 ^ 2;
t55 = t63 ^ 2;
t53 = t66 * pkin(5);
t37 = t60 * t64 + t69 * t81;
t36 = -t60 * t69 + t64 * t81;
t35 = t61 * t50 + t72;
t30 = t63 * t44 + t71;
t29 = -pkin(8) * t77 + t39;
t24 = t37 * t68 - t63 * t80;
t23 = -t37 * t63 - t68 * t80;
t22 = t32 * pkin(5) + t43;
t20 = -t61 * t40 + t66 * t41;
t19 = t66 * t40 + t61 * t41;
t17 = -t40 * pkin(11) + t27;
t16 = -t41 * pkin(11) + t26;
t15 = -t61 * t32 - t66 * t33;
t14 = t66 * t32 - t61 * t33;
t13 = t62 * t23 + t67 * t24;
t12 = t67 * t23 - t62 * t24;
t8 = t61 * t16 + t66 * t17;
t7 = t66 * t16 - t61 * t17;
t4 = t61 * t12 + t66 * t13;
t3 = t66 * t12 - t61 * t13;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t80, -t81, 0, 0, 0, 0, 0, t69 * t80, -t64 * t80, 0, 0, 0, 0, 0, -t23 * t69 + t36 * t79, t24 * t69 + t36 * t75, 0, 0, 0, 0, 0, -t12 * t69 + t36 * t32, t13 * t69 - t36 * t33, 0, 0, 0, 0, 0, t36 * t14 - t3 * t69, t36 * t15 + t4 * t69; 0, 1, 0, 0, t56, t73, 0, 0, 0, pkin(2) * t90, pkin(2) * t92, t57 * t56, -0.2e1 * t56 * t78, t74 * t92, t63 * t73, t58, -0.2e1 * t29 * t69 + 0.2e1 * t56 * t87, 0.2e1 * t56 * pkin(8) * t68 + 0.2e1 * t30 * t69, t33 ^ 2, 0.2e1 * t33 * t32, -t33 * t91, t32 * t90, t58, -0.2e1 * t10 * t69 + 0.2e1 * t43 * t32, 0.2e1 * t11 * t69 - 0.2e1 * t43 * t33, t15 ^ 2, -0.2e1 * t15 * t14, t15 * t91, t14 * t90, t58, -0.2e1 * t1 * t69 + 0.2e1 * t22 * t14, 0.2e1 * t22 * t15 + 0.2e1 * t2 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, 0, 0, 0, 0, -t36 * t68, t36 * t63, 0, 0, 0, 0, 0, t36 * t40, t36 * t41, 0, 0, 0, 0, 0, t36 * t19, t36 * t20; 0, 0, 0, 0, 0, 0, t64, t69, 0, -t52, -t69 * pkin(8), t63 * t75 (-t55 + t57) * t64, -t77, -t74, 0, -pkin(8) * t75 + (-pkin(3) * t64 + pkin(9) * t69) * t63, pkin(9) * t74 + (t87 - t88) * t64, -t33 * t41, -t41 * t32 + t33 * t40, -t41 * t69, t40 * t69, 0, -t26 * t69 + t51 * t32 + t43 * t40, t27 * t69 - t51 * t33 + t43 * t41, t15 * t20, -t20 * t14 - t15 * t19, -t20 * t69, t19 * t69, 0, t31 * t14 + t22 * t19 - t7 * t69, t31 * t15 + t22 * t20 + t8 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t55, 0.2e1 * t78, 0, 0, 0, 0.2e1 * t88, -0.2e1 * pkin(3) * t63, t41 ^ 2, -0.2e1 * t41 * t40, 0, 0, 0, t40 * t93, t41 * t93, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t94, t20 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, 0, 0, 0, 0, t12, -t13, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t79, -t69, t29, -t30, 0, 0, -t33, -t32, -t69, -t67 * t83 + t10, -t76 + (-t21 + t83) * t62, 0, 0, t15, -t14, -t69, -t34 * t69 + t1, t35 * t69 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t68, 0, -t63 * pkin(9), -t68 * pkin(9), 0, 0, t41, -t40, 0, t26, -t27, 0, 0, t20, -t19, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t85, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32, -t69, t10, -t11, 0, 0, t15, -t14, -t69, -t66 * t82 + t1, -t84 + (-t6 + t82) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t26, -t27, 0, 0, t20, -t19, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t54, -t85, 0, 0, 0, 0, 1, t34 + t53, -t72 + (-pkin(5) - t50) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, -t69, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
