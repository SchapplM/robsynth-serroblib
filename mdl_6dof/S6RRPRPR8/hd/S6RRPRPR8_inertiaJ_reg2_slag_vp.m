% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:58:59
% EndTime: 2019-05-06 14:59:06
% DurationCPUTime: 1.98s
% Computational Cost: add. (1687->176), mult. (3393->313), div. (0->0), fcn. (3753->8), ass. (0->93)
t115 = cos(qJ(4));
t72 = sin(pkin(10));
t73 = cos(pkin(10));
t75 = sin(qJ(4));
t55 = t115 * t72 + t75 * t73;
t76 = sin(qJ(2));
t45 = t55 * t76;
t53 = -t115 * t73 + t75 * t72;
t47 = t53 * t76;
t56 = (pkin(3) * t72 + pkin(7)) * t76;
t14 = t45 * pkin(4) + t47 * qJ(5) + t56;
t124 = t45 ^ 2;
t123 = t53 ^ 2;
t118 = -pkin(4) - pkin(5);
t64 = -t73 * pkin(3) - pkin(2);
t82 = t55 * qJ(5) - t64;
t15 = t118 * t53 + t82;
t122 = 0.2e1 * t15;
t121 = 0.2e1 * t64;
t78 = cos(qJ(2));
t120 = -0.2e1 * t78;
t119 = 0.2e1 * t78;
t117 = pkin(7) * t72;
t116 = t76 * pkin(7);
t99 = pkin(8) + qJ(3);
t88 = t99 * t73;
t89 = t99 * t72;
t35 = t115 * t89 + t75 * t88;
t114 = t35 * t78;
t37 = t115 * t88 - t75 * t89;
t113 = t37 * t78;
t112 = t45 * t53;
t111 = t47 * t45;
t110 = t55 * t53;
t109 = t55 * t78;
t108 = t72 * t73;
t107 = t72 * t76;
t106 = t72 * t78;
t105 = t73 * t76;
t104 = t73 * t78;
t102 = t76 * t78;
t101 = t78 * t45;
t100 = t78 * t53;
t60 = -t78 * pkin(2) - t76 * qJ(3) - pkin(1);
t50 = t73 * t60;
t29 = -pkin(8) * t105 + t50 + (-pkin(3) - t117) * t78;
t40 = pkin(7) * t104 + t72 * t60;
t33 = -pkin(8) * t107 + t40;
t98 = -t115 * t29 + t75 * t33;
t13 = t115 * t33 + t75 * t29;
t68 = t72 ^ 2;
t69 = t73 ^ 2;
t96 = t68 + t69;
t95 = t78 * qJ(5);
t94 = 0.2e1 * t102;
t93 = t35 ^ 2 + t37 ^ 2;
t92 = t72 * t105;
t67 = t78 * pkin(4);
t9 = t67 + t98;
t87 = -t35 * t47 - t37 * t45;
t8 = -t95 + t13;
t3 = t78 * pkin(5) + t47 * pkin(9) + t9;
t4 = t45 * pkin(9) + t8;
t74 = sin(qJ(6));
t77 = cos(qJ(6));
t1 = t77 * t3 - t74 * t4;
t2 = t74 * t3 + t77 * t4;
t86 = -pkin(2) * t76 + qJ(3) * t78;
t39 = -pkin(7) * t106 + t50;
t85 = -t39 * t72 + t40 * t73;
t84 = -t55 * t45 + t47 * t53;
t83 = 0.2e1 * t35 * t55 - 0.2e1 * t37 * t53;
t81 = -t55 * pkin(9) + t35;
t80 = pkin(7) ^ 2;
t71 = t78 ^ 2;
t70 = t76 ^ 2;
t66 = t70 * t80;
t59 = t77 * qJ(5) + t74 * t118;
t57 = t74 * qJ(5) - t77 * t118;
t51 = t55 ^ 2;
t44 = t47 ^ 2;
t38 = t47 * t120;
t30 = t47 * t55;
t28 = t74 * t53 + t77 * t55;
t26 = -t77 * t53 + t74 * t55;
t25 = t53 * pkin(4) - t82;
t21 = t53 * pkin(9) + t37;
t20 = t74 * t45 - t77 * t47;
t18 = -t77 * t45 - t74 * t47;
t10 = t45 * pkin(5) + t14;
t7 = t77 * t21 + t74 * t81;
t5 = t74 * t21 - t77 * t81;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t70, t94, 0, t71, 0, 0, pkin(1) * t119, -0.2e1 * pkin(1) * t76, 0.2e1 * (t70 + t71) * pkin(7), pkin(1) ^ 2 + t71 * t80 + t66, t69 * t70, -0.2e1 * t70 * t108, -0.2e1 * t73 * t102, t68 * t70, t72 * t94, t71, 0.2e1 * t70 * t117 - 0.2e1 * t39 * t78, 0.2e1 * t70 * pkin(7) * t73 + 0.2e1 * t40 * t78, 0.2e1 * (-t39 * t73 - t40 * t72) * t76, t39 ^ 2 + t40 ^ 2 + t66, t44, 0.2e1 * t111, -t38, t124, 0.2e1 * t101, t71, 0.2e1 * t56 * t45 + 0.2e1 * t78 * t98, 0.2e1 * t13 * t78 - 0.2e1 * t56 * t47, -0.2e1 * t13 * t45 - 0.2e1 * t47 * t98, t13 ^ 2 + t56 ^ 2 + t98 ^ 2, t44, -t38, -0.2e1 * t111, t71, -0.2e1 * t101, t124, 0.2e1 * t14 * t45 + 0.2e1 * t9 * t78, -0.2e1 * t8 * t45 - 0.2e1 * t9 * t47, 0.2e1 * t14 * t47 - 0.2e1 * t8 * t78, t14 ^ 2 + t8 ^ 2 + t9 ^ 2, t20 ^ 2, -0.2e1 * t20 * t18, t20 * t119, t18 ^ 2, t18 * t120, t71, 0.2e1 * t1 * t78 - 0.2e1 * t10 * t18, -0.2e1 * t10 * t20 - 0.2e1 * t2 * t78, -0.2e1 * t1 * t20 - 0.2e1 * t2 * t18, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, t78, 0, -t116, -t78 * pkin(7), 0, 0, t92 (-t68 + t69) * t76, -t106, -t92, -t104, 0, -pkin(7) * t105 + t86 * t72, pkin(7) * t107 + t86 * t73, t85, -pkin(2) * t116 + t85 * qJ(3), -t30, t84, -t109, t112, t100, 0, t64 * t45 + t56 * t53 + t114, -t64 * t47 + t56 * t55 + t113, -t13 * t53 + t55 * t98 + t87, t13 * t37 + t35 * t98 + t56 * t64, -t30, -t109, -t84, 0, -t100, t112, t14 * t53 + t25 * t45 + t114, -t8 * t53 + t9 * t55 + t87, -t14 * t55 + t25 * t47 - t113, t14 * t25 + t9 * t35 + t8 * t37, t20 * t28, -t28 * t18 - t20 * t26, t28 * t78, t18 * t26, -t26 * t78, 0, -t10 * t26 + t15 * t18 - t5 * t78, -t10 * t28 + t15 * t20 - t7 * t78, -t1 * t28 - t7 * t18 - t2 * t26 + t5 * t20, -t1 * t5 - t10 * t15 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t68, 0.2e1 * t108, 0, t69, 0, 0, 0.2e1 * pkin(2) * t73, -0.2e1 * pkin(2) * t72, 0.2e1 * t96 * qJ(3), t96 * qJ(3) ^ 2 + pkin(2) ^ 2, t51, -0.2e1 * t110, 0, t123, 0, 0, t53 * t121, t55 * t121, t83, t64 ^ 2 + t93, t51, 0, 0.2e1 * t110, 0, 0, t123, 0.2e1 * t25 * t53, t83, -0.2e1 * t25 * t55, t25 ^ 2 + t93, t28 ^ 2, -0.2e1 * t28 * t26, 0, t26 ^ 2, 0, 0, t26 * t122, t28 * t122, -0.2e1 * t7 * t26 + 0.2e1 * t5 * t28, t15 ^ 2 + t5 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t105, 0, t116, 0, 0, 0, 0, 0, 0, t45, -t47, 0, t56, 0, 0, 0, 0, 0, 0, t45, 0, t47, t14, 0, 0, 0, 0, 0, 0, -t18, -t20, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t72, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t53, t55, 0, t64, 0, 0, 0, 0, 0, 0, t53, 0, -t55, t25, 0, 0, 0, 0, 0, 0, -t26, -t28, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, -t45, -t78, -t98, -t13, 0, 0, 0, -t47, 0, -t78, t45, 0, -0.2e1 * t67 - t98, pkin(4) * t47 - t45 * qJ(5), -0.2e1 * t95 + t13, -t9 * pkin(4) + t8 * qJ(5), 0, 0, -t20, 0, t18, -t78, -t57 * t78 - t1, -t59 * t78 + t2, -t59 * t18 + t57 * t20, -t1 * t57 + t2 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t53, 0, -t35, -t37, 0, 0, 0, t55, 0, 0, t53, 0, -t35, -pkin(4) * t55 - t53 * qJ(5), t37, -t35 * pkin(4) + t37 * qJ(5), 0, 0, -t28, 0, t26, 0, t5, t7, -t59 * t26 + t57 * t28, t5 * t57 + t7 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, 0.2e1 * t59, 0, t57 ^ 2 + t59 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t47, 0, t9, 0, 0, 0, 0, 0, 0, t77 * t78, -t74 * t78, -t74 * t18 - t77 * t20, t1 * t77 + t2 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t26 - t77 * t28, -t5 * t77 + t7 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t77, t74, 0, -t57 * t77 + t59 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 ^ 2 + t77 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, t78, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, 0, -t5, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t57, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t6;
