% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t71 = cos(qJ(4));
t57 = t71 * pkin(3);
t53 = t57 + pkin(4);
t66 = sin(qJ(5));
t67 = sin(qJ(4));
t100 = t67 * pkin(3);
t70 = cos(qJ(5));
t86 = t70 * t100;
t35 = t66 * t53 + t86;
t33 = pkin(10) + t35;
t65 = sin(qJ(6));
t59 = t65 ^ 2;
t69 = cos(qJ(6));
t61 = t69 ^ 2;
t89 = t59 + t61;
t115 = t89 * t33;
t101 = t66 * pkin(4);
t51 = pkin(10) + t101;
t114 = t89 * t51;
t72 = cos(qJ(3));
t63 = sin(pkin(11));
t103 = t63 * pkin(1);
t49 = pkin(7) + t103;
t99 = pkin(8) + t49;
t84 = t99 * t72;
t68 = sin(qJ(3));
t85 = t99 * t68;
t15 = -t67 * t85 + t71 * t84;
t39 = -t67 * t68 + t71 * t72;
t12 = t39 * pkin(9) + t15;
t14 = -t67 * t84 - t71 * t85;
t40 = t67 * t72 + t71 * t68;
t77 = -t40 * pkin(9) + t14;
t5 = t66 * t12 - t70 * t77;
t113 = t5 ^ 2;
t24 = -t70 * t39 + t66 * t40;
t112 = t24 ^ 2;
t64 = cos(pkin(11));
t102 = t64 * pkin(1);
t50 = -pkin(2) - t102;
t41 = -t72 * pkin(3) + t50;
t27 = -t39 * pkin(4) + t41;
t111 = 0.2e1 * t27;
t110 = 0.2e1 * t40;
t109 = 0.2e1 * t68;
t108 = pkin(5) * t65;
t107 = t24 * pkin(5);
t26 = t66 * t39 + t70 * t40;
t106 = t26 * pkin(10);
t105 = t5 * t24;
t104 = t5 * t69;
t81 = t66 * t100 - t70 * t53;
t32 = -pkin(5) + t81;
t98 = t32 * t69;
t56 = t70 * pkin(4);
t52 = -t56 - pkin(5);
t97 = t52 * t69;
t96 = t59 * t26;
t95 = t65 * t26;
t94 = t65 * t69;
t93 = t69 * t26;
t90 = pkin(10) * t89;
t60 = t68 ^ 2;
t62 = t72 ^ 2;
t88 = t60 + t62;
t87 = -0.2e1 * t26 * t24;
t80 = -pkin(5) * t26 - pkin(10) * t24;
t7 = t70 * t12 + t66 * t77;
t8 = -t106 + t27 + t107;
t2 = -t65 * t7 + t69 * t8;
t3 = t65 * t8 + t69 * t7;
t1 = -t2 * t65 + t3 * t69;
t79 = -t24 * t33 + t26 * t32;
t78 = -t24 * t51 + t26 * t52;
t58 = pkin(5) * t69;
t47 = 0.2e1 * t94;
t45 = t52 * t65;
t38 = t40 ^ 2;
t37 = t39 ^ 2;
t30 = t32 * t65;
t23 = t26 ^ 2;
t20 = t69 * t24;
t19 = t61 * t26;
t18 = t61 * t23;
t17 = t65 * t24;
t16 = t59 * t23;
t13 = t65 * t93;
t10 = t19 + t96;
t9 = t19 - t96;
t4 = t5 * t65;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t102, -0.2e1 * t103, 0 (t63 ^ 2 + t64 ^ 2) * pkin(1) ^ 2, t60, t72 * t109, 0, t62, 0, 0, -0.2e1 * t50 * t72, t50 * t109, 0.2e1 * t88 * t49, t88 * t49 ^ 2 + t50 ^ 2, t38, t39 * t110, 0, t37, 0, 0, -0.2e1 * t41 * t39, t41 * t110, -0.2e1 * t14 * t40 + 0.2e1 * t15 * t39, t14 ^ 2 + t15 ^ 2 + t41 ^ 2, t23, t87, 0, t112, 0, 0, t24 * t111, t26 * t111, -0.2e1 * t7 * t24 + 0.2e1 * t5 * t26, t27 ^ 2 + t7 ^ 2 + t113, t18, -0.2e1 * t23 * t94, 0.2e1 * t24 * t93, t16, t65 * t87, t112, 0.2e1 * t2 * t24 + 0.2e1 * t5 * t95, -0.2e1 * t3 * t24 + 0.2e1 * t5 * t93, 0.2e1 * (-t2 * t69 - t3 * t65) * t26, t2 ^ 2 + t3 ^ 2 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t39 + t15 * t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t26 + t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t26 + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t112, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 + t16 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t72, 0, -t68 * t49, -t72 * t49, 0, 0, 0, 0, t40, 0, t39, 0, t14, -t15 (t39 * t67 - t40 * t71) * pkin(3) (t14 * t71 + t15 * t67) * pkin(3), 0, 0, t26, 0, -t24, 0, -t5, -t7, -t24 * t35 + t26 * t81, t35 * t7 + t5 * t81, t13, t9, t17, -t13, t20, 0, t79 * t65 - t104, t79 * t69 + t4, t1, t1 * t33 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t68, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, 0 (t39 * t71 + t40 * t67) * pkin(3), 0, 0, 0, 0, 0, 0, -t24, -t26, 0, t24 * t81 + t26 * t35, 0, 0, 0, 0, 0, 0, -t20, t17, t10, t115 * t26 + t24 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t100, 0 (t67 ^ 2 + t71 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t81, -0.2e1 * t35, 0, t35 ^ 2 + t81 ^ 2, t59, t47, 0, t61, 0, 0, -0.2e1 * t98, 0.2e1 * t30, 0.2e1 * t115, t89 * t33 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t39, 0, t14, -t15, 0, 0, 0, 0, t26, 0, -t24, 0, -t5, -t7 (-t24 * t66 - t26 * t70) * pkin(4) (-t5 * t70 + t66 * t7) * pkin(4), t13, t9, t17, -t13, t20, 0, t78 * t65 - t104, t78 * t69 + t4, t1, t1 * t51 + t5 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0 (-t24 * t70 + t26 * t66) * pkin(4), 0, 0, 0, 0, 0, 0, -t20, t17, t10, t114 * t26 + t24 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t100, 0, 0, 0, 0, 0, 0, 0, 1, t56 - t81, -t86 + (-pkin(4) - t53) * t66, 0 (t35 * t66 - t70 * t81) * pkin(4), t59, t47, 0, t61, 0, 0 (-t32 - t52) * t69, t45 + t30, t114 + t115, t114 * t33 + t32 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t101, 0 (t66 ^ 2 + t70 ^ 2) * pkin(4) ^ 2, t59, t47, 0, t61, 0, 0, -0.2e1 * t97, 0.2e1 * t45, 0.2e1 * t114, t89 * t51 ^ 2 + t52 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, -t5, -t7, 0, 0, t13, t9, t17, -t13, t20, 0, t80 * t65 - t104, t80 * t69 + t4, t1, -t5 * pkin(5) + t1 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t17, t10, t89 * t106 - t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t81, -t35, 0, 0, t59, t47, 0, t61, 0, 0, t58 - t98, t30 - t108, t90 + t115, -t32 * pkin(5) + pkin(10) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t101, 0, 0, t59, t47, 0, t61, 0, 0, t58 - t97, t45 - t108, t90 + t114, -t52 * pkin(5) + pkin(10) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t59, t47, 0, t61, 0, 0, 0.2e1 * t58, -0.2e1 * t108, 0.2e1 * t90, t89 * pkin(10) ^ 2 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, -t95, t24, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t93, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t69, 0, -t65 * t33, -t69 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t69, 0, -t65 * t51, -t69 * t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t69, 0, -t65 * pkin(10), -t69 * pkin(10), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t6;