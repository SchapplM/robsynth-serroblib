% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:47:43
% EndTime: 2019-05-06 20:47:47
% DurationCPUTime: 1.11s
% Computational Cost: add. (1795->142), mult. (4451->274), div. (0->0), fcn. (5302->12), ass. (0->112)
t113 = cos(qJ(5));
t74 = sin(pkin(12));
t75 = sin(pkin(6));
t76 = cos(pkin(12));
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t50 = (t74 * t84 + t76 * t81) * t75;
t77 = cos(pkin(6));
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t40 = t50 * t83 + t77 * t80;
t79 = sin(qJ(5));
t88 = -t50 * t80 + t77 * t83;
t24 = -t113 * t88 + t79 * t40;
t123 = -0.2e1 * t24;
t107 = t75 * t84;
t108 = t75 * t81;
t49 = -t76 * t107 + t74 * t108;
t122 = 0.2e1 * t49;
t67 = -t76 * pkin(2) - pkin(3);
t62 = -t83 * pkin(4) + t67;
t121 = 0.2e1 * t62;
t120 = 0.2e1 * t80;
t119 = pkin(1) * t81;
t117 = t49 * pkin(4);
t64 = t77 * t84 * pkin(1);
t97 = pkin(8) + qJ(3);
t41 = t77 * pkin(2) - t97 * t108 + t64;
t94 = t77 * t119;
t44 = t97 * t107 + t94;
t29 = t74 * t41 + t76 * t44;
t27 = t77 * pkin(9) + t29;
t61 = (-pkin(2) * t84 - pkin(1)) * t75;
t31 = t49 * pkin(3) - t50 * pkin(9) + t61;
t13 = -t80 * t27 + t83 * t31;
t11 = -t40 * pkin(10) + t117 + t13;
t14 = t83 * t27 + t80 * t31;
t12 = t88 * pkin(10) + t14;
t6 = t113 * t11 - t79 * t12;
t4 = -t49 * pkin(5) - t6;
t82 = cos(qJ(6));
t118 = t4 * t82;
t116 = t79 * pkin(4);
t91 = t113 * pkin(4);
t70 = -t91 - pkin(5);
t115 = pkin(5) - t70;
t66 = t74 * pkin(2) + pkin(9);
t114 = pkin(10) + t66;
t25 = t113 * t40 + t79 * t88;
t78 = sin(qJ(6));
t19 = t82 * t25 + t49 * t78;
t16 = t19 * t78;
t57 = t114 * t83;
t89 = t113 * t80;
t34 = t114 * t89 + t79 * t57;
t112 = t34 * t82;
t18 = t78 * t25 - t49 * t82;
t103 = t79 * t80;
t59 = -t113 * t83 + t103;
t111 = t59 * t18;
t60 = t79 * t83 + t89;
t110 = t60 * t49;
t71 = t75 ^ 2;
t109 = t71 * t84;
t22 = t78 * t24;
t106 = t78 * t60;
t69 = pkin(11) + t116;
t105 = t78 * t69;
t104 = t78 * t82;
t102 = t80 * t49;
t101 = t80 * t66;
t23 = t82 * t24;
t100 = t82 * t60;
t99 = t82 * t69;
t98 = t83 * t66;
t96 = -0.2e1 * t60 * t59;
t95 = 0.2e1 * t75 * t77;
t93 = t24 * t106;
t92 = t24 * t100;
t90 = t113 * t12;
t28 = t76 * t41 - t74 * t44;
t87 = -pkin(5) * t60 - pkin(11) * t59;
t86 = -t59 * t69 + t60 * t70;
t26 = -t77 * pkin(3) - t28;
t7 = t79 * t11 + t90;
t17 = -t88 * pkin(4) + t26;
t73 = t82 ^ 2;
t72 = t78 ^ 2;
t65 = 0.2e1 * t104;
t58 = t60 ^ 2;
t56 = pkin(8) * t107 + t94;
t55 = -pkin(8) * t108 + t64;
t54 = t82 * t59;
t53 = t78 * t59;
t51 = t78 * t100;
t48 = t49 ^ 2;
t45 = t83 * t49;
t38 = t59 * t49;
t37 = (-t72 + t73) * t60;
t35 = -t114 * t103 + t113 * t57;
t33 = t59 * pkin(5) - t60 * pkin(11) + t62;
t32 = t34 * t78;
t21 = t78 * t33 + t82 * t35;
t20 = t82 * t33 - t78 * t35;
t15 = t19 * t59;
t9 = -t78 * t18 + t19 * t82;
t8 = t24 * pkin(5) - t25 * pkin(11) + t17;
t5 = t49 * pkin(11) + t7;
t3 = t4 * t78;
t2 = t82 * t5 + t78 * t8;
t1 = -t78 * t5 + t82 * t8;
t10 = [1, 0, 0, t71 * t81 ^ 2, 0.2e1 * t81 * t109, t81 * t95, t84 * t95, t77 ^ 2, 0.2e1 * pkin(1) * t109 + 0.2e1 * t55 * t77, -0.2e1 * t71 * t119 - 0.2e1 * t56 * t77, -0.2e1 * t28 * t50 - 0.2e1 * t29 * t49, t28 ^ 2 + t29 ^ 2 + t61 ^ 2, t40 ^ 2, 0.2e1 * t40 * t88, t40 * t122, t88 * t122, t48, 0.2e1 * t13 * t49 - 0.2e1 * t26 * t88, -0.2e1 * t14 * t49 + 0.2e1 * t26 * t40, t25 ^ 2, t25 * t123, t25 * t122, t49 * t123, t48, 0.2e1 * t17 * t24 + 0.2e1 * t6 * t49, 0.2e1 * t17 * t25 - 0.2e1 * t7 * t49, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t24, t18 * t123, t24 ^ 2, 0.2e1 * t1 * t24 + 0.2e1 * t4 * t18, 0.2e1 * t4 * t19 - 0.2e1 * t2 * t24; 0, 0, 0, 0, 0, t108, t107, t77, t55, -t56 (-t49 * t74 - t50 * t76) * pkin(2) (t28 * t76 + t29 * t74) * pkin(2), t40 * t80, t40 * t83 + t80 * t88, t102, t45, 0, -t49 * t101 - t26 * t83 - t67 * t88, t26 * t80 + t67 * t40 - t49 * t98, t25 * t60, -t60 * t24 - t25 * t59, t110, -t38, 0, t17 * t59 + t62 * t24 - t34 * t49, t17 * t60 + t62 * t25 - t35 * t49, t19 * t100 (-t18 * t82 - t16) * t60, t15 + t92, -t93 - t111, t24 * t59, t1 * t59 + t4 * t106 + t34 * t18 + t20 * t24, t4 * t100 + t34 * t19 - t2 * t59 - t21 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t74 ^ 2 + t76 ^ 2) * pkin(2) ^ 2, t80 ^ 2, t83 * t120, 0, 0, 0, -0.2e1 * t67 * t83, t67 * t120, t58, t96, 0, 0, 0, t59 * t121, t60 * t121, t73 * t58, -0.2e1 * t58 * t104, 0.2e1 * t59 * t100, t78 * t96, t59 ^ 2, 0.2e1 * t34 * t106 + 0.2e1 * t20 * t59, 0.2e1 * t34 * t100 - 0.2e1 * t21 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, t45, -t102, 0, 0, 0, 0, 0, -t38, -t110, 0, 0, 0, 0, 0, -t93 + t111, t15 - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t88, t49, t13, -t14, 0, 0, t25, -t24, t49, t49 * t91 + t6, -t90 + (-t11 - t117) * t79, t16, t9, t22, t23, 0, -t24 * t105 + t70 * t18 - t118, t70 * t19 - t24 * t99 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t83, 0, -t101, -t98, 0, 0, t60, -t59, 0, -t34, -t35, t51, t37, t53, t54, 0, t86 * t78 - t112, t82 * t86 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t80, 0, 0, 0, 0, 0, -t59, -t60, 0, 0, 0, 0, 0, -t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t91, -0.2e1 * t116, t72, t65, 0, 0, 0, -0.2e1 * t70 * t82, 0.2e1 * t70 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, t49, t6, -t7, t16, t9, t22, t23, 0, -pkin(5) * t18 - pkin(11) * t22 - t118, -pkin(5) * t19 - pkin(11) * t23 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t59, 0, -t34, -t35, t51, t37, t53, t54, 0, t87 * t78 - t112, t82 * t87 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t60, 0, 0, 0, 0, 0, -t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t91, -t116, t72, t65, 0, 0, 0, t115 * t82, -t115 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t72, t65, 0, 0, 0, 0.2e1 * pkin(5) * t82, -0.2e1 * pkin(5) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t24, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t106, t59, t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82, 0, -t105, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82, 0, -t78 * pkin(11), -t82 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
