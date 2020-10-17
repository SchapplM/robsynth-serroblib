% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:49:25
% EndTime: 2019-05-06 04:49:30
% DurationCPUTime: 1.50s
% Computational Cost: add. (2993->164), mult. (8126->324), div. (0->0), fcn. (9777->14), ass. (0->119)
t86 = sin(pkin(6));
t106 = cos(pkin(13));
t107 = cos(pkin(7));
t96 = t107 * t106;
t108 = cos(pkin(6));
t85 = sin(pkin(7));
t99 = t108 * t85;
t135 = t86 * t96 + t99;
t100 = t86 * t106;
t84 = sin(pkin(13));
t127 = pkin(1) * t84;
t60 = qJ(2) * t100 + t108 * t127;
t43 = t135 * pkin(9) + t60;
t90 = sin(qJ(3));
t93 = cos(qJ(3));
t119 = t84 * t86;
t101 = pkin(1) * t106;
t74 = t108 * t101;
t49 = t108 * pkin(2) + t74 + (-t107 * pkin(9) - qJ(2)) * t119;
t54 = (-pkin(9) * t84 * t85 - t106 * pkin(2) - pkin(1)) * t86;
t94 = t107 * t49 + t54 * t85;
t27 = -t90 * t43 + t94 * t93;
t122 = cos(qJ(5));
t48 = t90 * t99 + (t84 * t93 + t90 * t96) * t86;
t57 = t85 * t100 - t108 * t107;
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t33 = t48 * t89 + t57 * t92;
t34 = t48 * t92 - t57 * t89;
t88 = sin(qJ(5));
t25 = t122 * t33 + t88 * t34;
t134 = -0.2e1 * t25;
t47 = t90 * t119 - t135 * t93;
t133 = -0.2e1 * t47;
t132 = 0.2e1 * t47;
t131 = -0.2e1 * t48;
t80 = -t92 * pkin(4) - pkin(3);
t130 = 0.2e1 * t80;
t129 = 0.2e1 * t92;
t128 = pkin(10) + pkin(11);
t125 = t47 * pkin(4);
t31 = t107 * t54 - t85 * t49;
t20 = t47 * pkin(3) - t48 * pkin(10) + t31;
t28 = t93 * t43 + t94 * t90;
t22 = -t57 * pkin(10) + t28;
t13 = t92 * t20 - t89 * t22;
t10 = -t34 * pkin(11) + t125 + t13;
t14 = t89 * t20 + t92 * t22;
t12 = -t33 * pkin(11) + t14;
t6 = t122 * t10 - t88 * t12;
t4 = -t47 * pkin(5) - t6;
t91 = cos(qJ(6));
t126 = t4 * t91;
t124 = t88 * pkin(4);
t104 = t122 * pkin(4);
t79 = -t104 - pkin(5);
t123 = pkin(5) - t79;
t26 = t122 * t34 - t88 * t33;
t87 = sin(qJ(6));
t18 = t91 * t26 + t47 * t87;
t16 = t18 * t87;
t118 = t85 * t90;
t61 = t92 * t107 - t89 * t118;
t62 = t89 * t107 + t92 * t118;
t38 = -t122 * t61 + t88 * t62;
t121 = t38 * t91;
t102 = t122 * t89;
t71 = t128 * t92;
t52 = t128 * t102 + t88 * t71;
t120 = t52 * t91;
t117 = t85 * t93;
t23 = t87 * t25;
t67 = t88 * t92 + t102;
t116 = t87 * t67;
t78 = pkin(12) + t124;
t115 = t87 * t78;
t114 = t87 * t91;
t113 = t88 * t89;
t112 = t89 * t47;
t24 = t91 * t25;
t111 = t91 * t67;
t110 = t91 * t78;
t109 = t92 * t47;
t66 = -t122 * t92 + t113;
t105 = -0.2e1 * t67 * t66;
t103 = t122 * t12;
t98 = -pkin(5) * t67 - pkin(12) * t66;
t97 = -t66 * t78 + t67 * t79;
t7 = t88 * t10 + t103;
t21 = t57 * pkin(3) - t27;
t15 = t33 * pkin(4) + t21;
t83 = t91 ^ 2;
t82 = t87 ^ 2;
t81 = t86 ^ 2;
t75 = 0.2e1 * t114;
t65 = t67 ^ 2;
t64 = t91 * t66;
t63 = t87 * t66;
t59 = -qJ(2) * t119 + t74;
t56 = t87 * t111;
t53 = -t128 * t113 + t122 * t71;
t50 = t52 * t87;
t46 = t47 ^ 2;
t45 = (-t82 + t83) * t67;
t44 = t66 * pkin(5) - t67 * pkin(12) + t80;
t39 = t122 * t62 + t88 * t61;
t37 = t38 * t87;
t36 = -t87 * t117 + t91 * t39;
t35 = -t91 * t117 - t87 * t39;
t30 = t87 * t44 + t91 * t53;
t29 = t91 * t44 - t87 * t53;
t17 = t87 * t26 - t47 * t91;
t11 = -t87 * t17 + t18 * t91;
t8 = t25 * pkin(5) - t26 * pkin(12) + t15;
t5 = t47 * pkin(12) + t7;
t3 = t4 * t87;
t2 = t91 * t5 + t87 * t8;
t1 = -t87 * t5 + t91 * t8;
t9 = [1, 0, 0, 0.2e1 * t81 * t101 + 0.2e1 * t59 * t108, -0.2e1 * t60 * t108 - 0.2e1 * t81 * t127, 0.2e1 * (t106 * t60 - t59 * t84) * t86, t81 * pkin(1) ^ 2 + t59 ^ 2 + t60 ^ 2, t48 ^ 2, t47 * t131, t57 * t131, t57 * t132, t57 ^ 2, -0.2e1 * t27 * t57 + 0.2e1 * t31 * t47, 0.2e1 * t28 * t57 + 0.2e1 * t31 * t48, t34 ^ 2, -0.2e1 * t34 * t33, t34 * t132, t33 * t133, t46, 0.2e1 * t13 * t47 + 0.2e1 * t21 * t33, -0.2e1 * t14 * t47 + 0.2e1 * t21 * t34, t26 ^ 2, t26 * t134, t26 * t132, t25 * t133, t46, 0.2e1 * t15 * t25 + 0.2e1 * t6 * t47, 0.2e1 * t15 * t26 - 0.2e1 * t7 * t47, t18 ^ 2, -0.2e1 * t18 * t17, 0.2e1 * t18 * t25, t17 * t134, t25 ^ 2, 0.2e1 * t1 * t25 + 0.2e1 * t4 * t17, 0.2e1 * t4 * t18 - 0.2e1 * t2 * t25; 0, 0, 0, -t100, t119, 0, -t86 * pkin(1), 0, 0, 0, 0, 0, t107 * t47 - t57 * t117, t107 * t48 + t57 * t118, 0, 0, 0, 0, 0, -t33 * t117 + t61 * t47, -t34 * t117 - t62 * t47, 0, 0, 0, 0, 0, -t25 * t117 - t38 * t47, -t26 * t117 - t39 * t47, 0, 0, 0, 0, 0, t38 * t17 + t35 * t25, t38 * t18 - t36 * t25; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, -t57, t27, -t28, t34 * t89, -t89 * t33 + t34 * t92, t112, t109, 0, -pkin(3) * t33 - pkin(10) * t112 - t21 * t92, -pkin(3) * t34 - pkin(10) * t109 + t21 * t89, t26 * t67, -t67 * t25 - t26 * t66, t67 * t47, -t66 * t47, 0, t15 * t66 + t80 * t25 - t52 * t47, t15 * t67 + t80 * t26 - t53 * t47, t18 * t111 (-t17 * t91 - t16) * t67, t111 * t25 + t18 * t66, -t116 * t25 - t17 * t66, t25 * t66, t1 * t66 + t116 * t4 + t52 * t17 + t29 * t25, t111 * t4 + t52 * t18 - t2 * t66 - t30 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, -t118, 0, 0, 0, 0, 0, t92 * t117, -t89 * t117, 0, 0, 0, 0, 0, -t66 * t117, -t67 * t117, 0, 0, 0, 0, 0, t116 * t38 + t35 * t66, t111 * t38 - t36 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t89 ^ 2, t89 * t129, 0, 0, 0, pkin(3) * t129, -0.2e1 * pkin(3) * t89, t65, t105, 0, 0, 0, t66 * t130, t67 * t130, t83 * t65, -0.2e1 * t65 * t114, 0.2e1 * t66 * t111, t87 * t105, t66 ^ 2, 0.2e1 * t116 * t52 + 0.2e1 * t29 * t66, 0.2e1 * t111 * t52 - 0.2e1 * t30 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, t47, t13, -t14, 0, 0, t26, -t25, t47, t47 * t104 + t6, -t103 + (-t10 - t125) * t88, t16, t11, t23, t24, 0, -t115 * t25 + t79 * t17 - t126, -t110 * t25 + t79 * t18 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62, 0, 0, 0, 0, 0, -t38, -t39, 0, 0, 0, 0, 0, -t121, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t92, 0, -t89 * pkin(10), -t92 * pkin(10), 0, 0, t67, -t66, 0, -t52, -t53, t56, t45, t63, t64, 0, t87 * t97 - t120, t91 * t97 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t104, -0.2e1 * t124, t82, t75, 0, 0, 0, -0.2e1 * t79 * t91, 0.2e1 * t79 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t47, t6, -t7, t16, t11, t23, t24, 0, -pkin(5) * t17 - pkin(12) * t23 - t126, -pkin(5) * t18 - pkin(12) * t24 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t39, 0, 0, 0, 0, 0, -t121, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t66, 0, -t52, -t53, t56, t45, t63, t64, 0, t87 * t98 - t120, t91 * t98 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t104, -t124, t82, t75, 0, 0, 0, t123 * t91, -t123 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t82, t75, 0, 0, 0, 0.2e1 * pkin(5) * t91, -0.2e1 * pkin(5) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t116, t66, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t91, 0, -t115, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t91, 0, -t87 * pkin(12), -t91 * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
