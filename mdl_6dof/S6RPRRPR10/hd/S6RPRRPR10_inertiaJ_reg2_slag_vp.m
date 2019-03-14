% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t69 = sin(qJ(4));
t63 = t69 ^ 2;
t72 = cos(qJ(4));
t65 = t72 ^ 2;
t88 = t63 + t65;
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t24 = t68 * t69 + t71 * t72;
t73 = cos(qJ(3));
t19 = t73 * t24;
t110 = -0.2e1 * t19;
t104 = pkin(4) + pkin(5);
t87 = t69 * qJ(5);
t21 = t104 * t72 + pkin(3) + t87;
t109 = 0.2e1 * t21;
t108 = -0.2e1 * t69;
t107 = 0.2e1 * t72;
t106 = 0.2e1 * t73;
t105 = 2 * qJ(2);
t70 = sin(qJ(3));
t103 = pkin(8) * t70;
t102 = t69 * pkin(8);
t101 = t70 * pkin(4);
t100 = t73 * pkin(3);
t16 = t24 * t70;
t66 = t73 ^ 2;
t74 = -pkin(1) - pkin(7);
t99 = t66 * t74;
t98 = t69 * t72;
t52 = t69 * t73;
t97 = t70 * t74;
t96 = t71 * t69;
t95 = t71 * t70;
t54 = t72 * t70;
t55 = t72 * t73;
t94 = t72 * t74;
t79 = -t72 * pkin(4) - t87;
t38 = -pkin(3) + t79;
t93 = t73 * t38;
t92 = t73 * t70;
t58 = t73 * t74;
t37 = t70 * pkin(3) - t73 * pkin(8) + qJ(2);
t91 = -t72 * t37 + t69 * t97;
t12 = t69 * t37 + t70 * t94;
t90 = t88 * t103;
t89 = t88 * pkin(8) ^ 2;
t64 = t70 ^ 2;
t46 = t64 + t66;
t86 = t72 * qJ(5);
t85 = t69 * t92;
t84 = t66 * t98;
t59 = t70 * qJ(5);
t6 = t59 + t12;
t83 = (pkin(8) - pkin(9)) * t69;
t82 = -t100 - t103;
t3 = -pkin(9) * t55 - t104 * t70 + t91;
t4 = pkin(9) * t52 + t6;
t1 = t71 * t3 - t68 * t4;
t2 = t68 * t3 + t71 * t4;
t7 = t91 - t101;
t81 = t6 * t72 + t7 * t69;
t80 = -t93 + t103;
t78 = pkin(4) * t69 - t86;
t77 = t12 * t72 + t69 * t91;
t75 = qJ(2) ^ 2;
t67 = t74 ^ 2;
t62 = t72 * pkin(8);
t57 = t66 * t67;
t53 = t65 * t66;
t51 = t69 * t70;
t50 = t63 * t66;
t42 = t69 * t55;
t40 = -t72 * pkin(9) + t62;
t39 = t92 * t107;
t36 = t71 * qJ(5) - t68 * t104;
t34 = t68 * qJ(5) + t71 * t104;
t33 = 0.2e1 * t88 * pkin(8);
t32 = t46 * t74;
t30 = t46 * t72;
t29 = t88 * t70;
t28 = -t68 * t72 + t96;
t27 = t46 * t69;
t26 = (-t63 + t65) * t73;
t20 = t88 * t64 + t66;
t17 = t68 * t55 - t73 * t96;
t14 = t68 * t54 - t69 * t95;
t13 = t73 * t78 - t58;
t10 = t71 * t40 + t68 * t83;
t8 = t68 * t40 - t71 * t83;
t5 = t58 + (-t104 * t69 + t86) * t73;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t105, pkin(1) ^ 2 + t75, t66, -0.2e1 * t92, 0, t64, 0, 0, t70 * t105, t73 * t105, -0.2e1 * t32, t64 * t67 + t57 + t75, t53, -0.2e1 * t84, t39, t50, -0.2e1 * t85, t64, -0.2e1 * t69 * t99 - 0.2e1 * t70 * t91, -0.2e1 * t12 * t70 - 0.2e1 * t66 * t94 (-t12 * t69 + t72 * t91) * t106, t12 ^ 2 + t91 ^ 2 + t57, t53, t39, 0.2e1 * t84, t64, 0.2e1 * t85, t50, 0.2e1 * t13 * t52 - 0.2e1 * t7 * t70 (-t6 * t69 + t7 * t72) * t106, -0.2e1 * t13 * t55 + 0.2e1 * t6 * t70, t13 ^ 2 + t6 ^ 2 + t7 ^ 2, t19 ^ 2, t17 * t110, t70 * t110, t17 ^ 2, 0.2e1 * t17 * t70, t64, -0.2e1 * t1 * t70 + 0.2e1 * t5 * t17, 0.2e1 * t5 * t19 + 0.2e1 * t2 * t70, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t46, t32, 0, 0, 0, 0, 0, 0, -t27, -t30, 0, t70 * t77 + t99, 0, 0, 0, 0, 0, 0, -t27, 0, t30, -t13 * t73 + t70 * t81, 0, 0, 0, 0, 0, 0, t14 * t70 + t73 * t17, t16 * t70 + t73 * t19, t14 * t19 - t16 * t17, -t1 * t14 + t2 * t16 + t5 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t16 ^ 2 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t70, 0, t58, -t97, 0, 0, t42, t26, t51, -t42, t54, 0, t72 * t58 + t82 * t69, -t69 * t58 + t82 * t72, t77, pkin(3) * t58 + pkin(8) * t77, t42, t51, -t26, 0, -t54, -t42, -t13 * t72 - t69 * t80, t81, -t13 * t69 + t72 * t80, pkin(8) * t81 + t13 * t38, t19 * t28, -t28 * t17 - t19 * t24, -t28 * t70, t17 * t24, t16, 0, t21 * t17 + t5 * t24 + t8 * t70, t10 * t70 + t21 * t19 + t5 * t28, -t1 * t28 - t10 * t17 + t8 * t19 - t2 * t24, -t1 * t8 + t2 * t10 + t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t70, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t52, t29, t90 + t100, 0, 0, 0, 0, 0, 0, t55, t29, t52, t90 - t93, 0, 0, 0, 0, 0, 0, t19, t73 * t28, t14 * t28 - t16 * t24, t16 * t10 + t14 * t8 + t73 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t63, 0.2e1 * t98, 0, t65, 0, 0, pkin(3) * t107, pkin(3) * t108, t33, pkin(3) ^ 2 + t89, t63, 0, -0.2e1 * t98, 0, 0, t65, -0.2e1 * t38 * t72, t33, t38 * t108, t38 ^ 2 + t89, t28 ^ 2, -0.2e1 * t28 * t24, 0, t24 ^ 2, 0, 0, t24 * t109, t28 * t109, -0.2e1 * t10 * t24 + 0.2e1 * t8 * t28, t10 ^ 2 + t21 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t52, t70, -t91, -t12, 0, 0, 0, t55, 0, t70, t52, 0, -t91 + 0.2e1 * t101, t79 * t73, 0.2e1 * t59 + t12, -t7 * pkin(4) + t6 * qJ(5), 0, 0, -t19, 0, t17, t70, t34 * t70 - t1, t36 * t70 + t2, -t36 * t17 + t34 * t19, -t1 * t34 + t2 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t54, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, t54, -t78 * t70, 0, 0, 0, 0, 0, 0, t14, t16, 0, t14 * t34 + t16 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t72, 0, -t102, -t62, 0, 0, 0, t69, 0, 0, -t72, 0, -t102, -t78, t62, -t78 * pkin(8), 0, 0, -t28, 0, t24, 0, t8, t10, -t36 * t24 + t34 * t28, t10 * t36 + t8 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, 0.2e1 * t36, 0, t34 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t55, 0, t7, 0, 0, 0, 0, 0, 0, -t95, t68 * t70, -t68 * t17 - t71 * t19, t1 * t71 + t2 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t71 + t16 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, -t68 * t24 - t71 * t28, t10 * t68 - t8 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t71, t68, 0, -t34 * t71 + t36 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 ^ 2 + t71 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, -t70, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t24, 0, -t8, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t34, -t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t9;