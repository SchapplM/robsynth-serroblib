% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR11
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 05:39:01
% EndTime: 2019-05-06 05:39:06
% DurationCPUTime: 1.66s
% Computational Cost: add. (3038->190), mult. (8310->371), div. (0->0), fcn. (9941->14), ass. (0->120)
t85 = cos(pkin(6));
t126 = pkin(1) * t85;
t80 = sin(pkin(13));
t83 = cos(pkin(13));
t82 = sin(pkin(6));
t98 = qJ(2) * t82;
t56 = t80 * t126 + t83 * t98;
t81 = sin(pkin(7));
t84 = cos(pkin(7));
t110 = t83 * t84;
t95 = t82 * t110;
t37 = (t81 * t85 + t95) * pkin(9) + t56;
t89 = sin(qJ(3));
t93 = cos(qJ(3));
t114 = t80 * t82;
t71 = t83 * t126;
t41 = t85 * pkin(2) + t71 + (-pkin(9) * t84 - qJ(2)) * t114;
t48 = (-pkin(9) * t80 * t81 - pkin(2) * t83 - pkin(1)) * t82;
t94 = t41 * t84 + t48 * t81;
t21 = -t89 * t37 + t94 * t93;
t113 = t81 * t89;
t40 = t85 * t113 + (t89 * t110 + t80 * t93) * t82;
t111 = t82 * t83;
t54 = -t81 * t111 + t85 * t84;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t32 = t40 * t88 - t54 * t92;
t135 = -0.2e1 * t32;
t134 = 0.2e1 * t32;
t112 = t81 * t93;
t39 = -t85 * t112 + t89 * t114 - t93 * t95;
t133 = -0.2e1 * t39;
t86 = sin(qJ(6));
t87 = sin(qJ(5));
t90 = cos(qJ(6));
t91 = cos(qJ(5));
t60 = t86 * t87 - t90 * t91;
t53 = t60 * t88;
t132 = 0.2e1 * t53;
t74 = -t91 * pkin(5) - pkin(4);
t131 = 0.2e1 * t74;
t130 = -0.2e1 * t88;
t129 = 0.2e1 * t92;
t128 = pkin(11) + pkin(12);
t75 = t82 ^ 2;
t127 = pkin(1) * t75;
t125 = pkin(4) * t91;
t124 = pkin(10) * t87;
t123 = t32 * pkin(5);
t122 = t86 * pkin(5);
t30 = -t81 * t41 + t84 * t48;
t17 = t39 * pkin(3) - t40 * pkin(10) + t30;
t22 = t93 * t37 + t94 * t89;
t20 = t54 * pkin(10) + t22;
t11 = t92 * t17 - t88 * t20;
t9 = -t39 * pkin(4) - t11;
t121 = t9 * t87;
t120 = t9 * t91;
t119 = t90 * pkin(5);
t33 = t40 * t92 + t54 * t88;
t24 = t33 * t87 - t39 * t91;
t12 = t88 * t17 + t92 * t20;
t10 = t39 * pkin(11) + t12;
t19 = -t54 * pkin(3) - t21;
t14 = t32 * pkin(4) - t33 * pkin(11) + t19;
t7 = t91 * t10 + t87 * t14;
t5 = -t24 * pkin(12) + t7;
t118 = t90 * t5;
t117 = t92 * pkin(5);
t25 = t33 * t91 + t39 * t87;
t116 = t25 * t87;
t115 = t32 * t92;
t109 = t87 * t32;
t108 = t87 * t88;
t107 = t87 * t91;
t106 = t87 * t92;
t105 = t88 * t32;
t104 = t88 * t39;
t64 = -t92 * pkin(4) - t88 * pkin(11) - pkin(3);
t100 = t91 * t92;
t96 = pkin(10) * t100;
t44 = t96 + (-pkin(12) * t88 + t64) * t87;
t103 = t90 * t44;
t102 = t91 * t32;
t101 = t91 * t88;
t99 = t92 * t39;
t97 = t88 * t129;
t6 = -t87 * t10 + t91 * t14;
t4 = -t25 * pkin(12) + t123 + t6;
t1 = t90 * t4 - t86 * t5;
t59 = t91 * t64;
t38 = -pkin(12) * t101 + t59 + (-pkin(5) - t124) * t92;
t26 = t90 * t38 - t86 * t44;
t61 = t86 * t91 + t90 * t87;
t79 = t92 ^ 2;
t78 = t91 ^ 2;
t77 = t88 ^ 2;
t76 = t87 ^ 2;
t68 = t128 * t91;
t67 = t128 * t87;
t63 = (pkin(5) * t87 + pkin(10)) * t88;
t58 = t92 * t113 + t88 * t84;
t57 = t88 * t113 - t92 * t84;
t55 = -t80 * t98 + t71;
t52 = t61 * t88;
t50 = t87 * t64 + t96;
t49 = -pkin(10) * t106 + t59;
t47 = -t86 * t67 + t90 * t68;
t46 = -t90 * t67 - t86 * t68;
t43 = -t87 * t112 + t91 * t58;
t42 = -t91 * t112 - t87 * t58;
t31 = t32 ^ 2;
t29 = t86 * t42 + t90 * t43;
t28 = t90 * t42 - t86 * t43;
t27 = t86 * t38 + t103;
t16 = -t86 * t24 + t90 * t25;
t15 = t90 * t24 + t86 * t25;
t8 = t24 * pkin(5) + t9;
t2 = t86 * t4 + t118;
t3 = [1, 0, 0, 0.2e1 * t83 * t127 + 0.2e1 * t55 * t85, -0.2e1 * t80 * t127 - 0.2e1 * t56 * t85, 0.2e1 * (-t55 * t80 + t56 * t83) * t82, t75 * pkin(1) ^ 2 + t55 ^ 2 + t56 ^ 2, t40 ^ 2, t40 * t133, 0.2e1 * t40 * t54, t54 * t133, t54 ^ 2, 0.2e1 * t21 * t54 + 0.2e1 * t30 * t39, -0.2e1 * t22 * t54 + 0.2e1 * t30 * t40, t33 ^ 2, t33 * t135, 0.2e1 * t33 * t39, t32 * t133, t39 ^ 2, 0.2e1 * t11 * t39 + 0.2e1 * t19 * t32, -0.2e1 * t12 * t39 + 0.2e1 * t19 * t33, t25 ^ 2, -0.2e1 * t25 * t24, t25 * t134, t24 * t135, t31, 0.2e1 * t9 * t24 + 0.2e1 * t6 * t32, 0.2e1 * t9 * t25 - 0.2e1 * t7 * t32, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t134, t15 * t135, t31, 0.2e1 * t1 * t32 + 0.2e1 * t8 * t15, 0.2e1 * t8 * t16 - 0.2e1 * t2 * t32; 0, 0, 0, -t111, t114, 0, -t82 * pkin(1), 0, 0, 0, 0, 0, t54 * t112 + t84 * t39, -t54 * t113 + t84 * t40, 0, 0, 0, 0, 0, -t32 * t112 - t57 * t39, -t33 * t112 - t58 * t39, 0, 0, 0, 0, 0, t57 * t24 + t42 * t32, t57 * t25 - t43 * t32, 0, 0, 0, 0, 0, t57 * t15 + t28 * t32, t57 * t16 - t29 * t32; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t39, t54, t21, -t22, t33 * t88, t33 * t92 - t105, t104, t99, 0, -pkin(3) * t32 - pkin(10) * t104 - t19 * t92, -pkin(3) * t33 - pkin(10) * t99 + t19 * t88, t25 * t101 (-t24 * t91 - t116) * t88, t32 * t101 - t25 * t92, -t87 * t105 + t24 * t92, -t115, t49 * t32 - t6 * t92 + (pkin(10) * t24 + t121) * t88, -t50 * t32 + t7 * t92 + (pkin(10) * t25 + t120) * t88, -t16 * t53, t53 * t15 - t16 * t52, -t16 * t92 - t53 * t32, t15 * t92 - t52 * t32, -t115, -t1 * t92 + t63 * t15 + t26 * t32 + t8 * t52, t63 * t16 + t2 * t92 - t27 * t32 - t8 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t113, 0, 0, 0, 0, 0, t92 * t112, -t88 * t112, 0, 0, 0, 0, 0, t57 * t108 - t42 * t92, t57 * t101 + t43 * t92, 0, 0, 0, 0, 0, -t28 * t92 + t57 * t52, t29 * t92 - t57 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t77, t97, 0, 0, 0, pkin(3) * t129, pkin(3) * t130, t78 * t77, -0.2e1 * t77 * t107, t100 * t130, t87 * t97, t79, 0.2e1 * t77 * t124 - 0.2e1 * t49 * t92, 0.2e1 * t77 * pkin(10) * t91 + 0.2e1 * t50 * t92, t53 ^ 2, t52 * t132, t92 * t132, t52 * t129, t79, -0.2e1 * t26 * t92 + 0.2e1 * t63 * t52, 0.2e1 * t27 * t92 - 0.2e1 * t63 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t39, t11, -t12, t116, -t87 * t24 + t25 * t91, t109, t102, 0, -pkin(4) * t24 - pkin(11) * t109 - t120, -pkin(4) * t25 - pkin(11) * t102 + t121, t16 * t61, -t61 * t15 - t16 * t60, t61 * t32, -t60 * t32, 0, t74 * t15 + t46 * t32 + t8 * t60, t74 * t16 - t47 * t32 + t8 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t58, 0, 0, 0, 0, 0, -t57 * t91, t57 * t87, 0, 0, 0, 0, 0, t57 * t60, t57 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t92, 0, -t88 * pkin(10), -t92 * pkin(10), t87 * t101 (-t76 + t78) * t88, -t106, -t100, 0, -pkin(10) * t101 + (-pkin(4) * t88 + pkin(11) * t92) * t87, pkin(11) * t100 + (t124 - t125) * t88, -t53 * t61, -t61 * t52 + t53 * t60, -t61 * t92, t60 * t92, 0, -t46 * t92 + t74 * t52 + t63 * t60, t47 * t92 - t74 * t53 + t63 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t76, 0.2e1 * t107, 0, 0, 0, 0.2e1 * t125, -0.2e1 * pkin(4) * t87, t61 ^ 2, -0.2e1 * t61 * t60, 0, 0, 0, t60 * t131, t61 * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, t32, t6, -t7, 0, 0, t16, -t15, t32, t32 * t119 + t1, -t118 + (-t4 - t123) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t43, 0, 0, 0, 0, 0, t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t108, -t92, t49, -t50, 0, 0, -t53, -t52, -t92, -t117 * t90 + t26, -t103 + (-t38 + t117) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t91, 0, -t87 * pkin(11), -t91 * pkin(11), 0, 0, t61, -t60, 0, t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t119, -0.2e1 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t32, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, -t92, t26, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t60, 0, t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t119, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
