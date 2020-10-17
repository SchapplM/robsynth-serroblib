% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:32:27
% EndTime: 2019-05-08 06:32:35
% DurationCPUTime: 1.62s
% Computational Cost: add. (2712->228), mult. (6148->413), div. (0->0), fcn. (6989->10), ass. (0->118)
t82 = sin(qJ(4));
t83 = sin(qJ(3));
t106 = t82 * t83;
t81 = sin(qJ(5));
t121 = cos(qJ(5));
t85 = cos(qJ(4));
t90 = t121 * t85;
t46 = -t106 * t81 + t83 * t90;
t135 = -0.2e1 * t46;
t79 = sin(pkin(6));
t84 = sin(qJ(2));
t111 = t79 * t84;
t80 = cos(pkin(6));
t86 = cos(qJ(3));
t50 = t111 * t83 - t80 * t86;
t134 = -0.2e1 * t50;
t133 = 0.2e1 * t50;
t51 = t111 * t86 + t80 * t83;
t132 = -0.2e1 * t51;
t91 = t121 * t82;
t56 = t81 * t85 + t91;
t131 = -0.2e1 * t56;
t71 = -t85 * pkin(4) - pkin(3);
t130 = 0.2e1 * t71;
t129 = -0.2e1 * t83;
t128 = 0.2e1 * t86;
t127 = -pkin(11) - pkin(10);
t87 = cos(qJ(2));
t125 = pkin(1) * t87;
t63 = pkin(8) * t111;
t41 = t63 + (-pkin(2) - t125) * t80;
t21 = t50 * pkin(3) - t51 * pkin(10) + t41;
t110 = t79 * t87;
t126 = pkin(1) * t84;
t95 = pkin(8) * t110;
t42 = t95 + (pkin(9) + t126) * t80;
t43 = (-pkin(2) * t87 - pkin(9) * t84 - pkin(1)) * t79;
t26 = t86 * t42 + t83 * t43;
t23 = -pkin(10) * t110 + t26;
t12 = t82 * t21 + t85 * t23;
t33 = t110 * t85 + t51 * t82;
t10 = -t33 * pkin(11) + t12;
t11 = t85 * t21 - t82 * t23;
t34 = -t110 * t82 + t51 * t85;
t8 = t50 * pkin(4) - t34 * pkin(11) + t11;
t4 = t121 * t10 + t81 * t8;
t124 = pkin(3) * t85;
t123 = pkin(9) * t82;
t48 = t50 * pkin(5);
t72 = t81 * pkin(4);
t122 = t86 * pkin(5);
t25 = -t83 * t42 + t86 * t43;
t22 = pkin(3) * t110 - t25;
t120 = t22 * t82;
t119 = t22 * t85;
t118 = t34 * t82;
t60 = t127 * t85;
t36 = -t127 * t91 - t81 * t60;
t117 = t36 * t50;
t116 = t36 * t86;
t108 = t81 * t82;
t37 = t108 * t127 - t121 * t60;
t115 = t37 * t50;
t114 = t37 * t86;
t113 = t50 * t86;
t74 = t79 ^ 2;
t112 = t74 * t87;
t109 = t80 * t84;
t107 = t82 * t50;
t105 = t82 * t85;
t104 = t82 * t86;
t103 = t83 * t50;
t102 = t85 * t50;
t101 = t85 * t83;
t100 = t85 * t86;
t59 = -t86 * pkin(3) - t83 * pkin(10) - pkin(2);
t54 = t85 * t59;
t31 = -pkin(11) * t101 + t54 + (-pkin(4) - t123) * t86;
t96 = pkin(9) * t100;
t35 = t96 + (-pkin(11) * t83 + t59) * t82;
t17 = t121 * t35 + t81 * t31;
t73 = t83 * pkin(9);
t58 = pkin(4) * t106 + t73;
t99 = t86 * qJ(6);
t98 = 0.2e1 * t110;
t97 = t83 * t128;
t47 = t50 * qJ(6);
t1 = t47 + t4;
t94 = t83 * t110;
t93 = t86 * t110;
t92 = t121 * pkin(4);
t3 = -t81 * t10 + t121 * t8;
t16 = t121 * t31 - t81 * t35;
t2 = -t3 - t48;
t13 = t33 * pkin(4) + t22;
t89 = 0.2e1 * pkin(5);
t88 = 0.2e1 * qJ(6);
t78 = t86 ^ 2;
t77 = t85 ^ 2;
t76 = t83 ^ 2;
t75 = t82 ^ 2;
t69 = t92 + pkin(5);
t67 = t72 + qJ(6);
t55 = -t90 + t108;
t53 = pkin(1) * t109 + t95;
t52 = t125 * t80 - t63;
t49 = t50 ^ 2;
t45 = t56 * t83;
t40 = t82 * t59 + t96;
t39 = -pkin(9) * t104 + t54;
t29 = t55 * pkin(5) - t56 * qJ(6) + t71;
t24 = t45 * pkin(5) - t46 * qJ(6) + t58;
t19 = t121 * t34 - t81 * t33;
t18 = t121 * t33 + t81 * t34;
t15 = -t16 + t122;
t14 = -t99 + t17;
t5 = t18 * pkin(5) - t19 * qJ(6) + t13;
t6 = [1, 0, 0, t74 * t84 ^ 2, 0.2e1 * t84 * t112, 0.2e1 * t79 * t109, t80 * t98, t80 ^ 2, 0.2e1 * pkin(1) * t112 + 0.2e1 * t52 * t80, -0.2e1 * t126 * t74 - 0.2e1 * t53 * t80, t51 ^ 2, t50 * t132, t110 * t132, t50 * t98, t74 * t87 ^ 2, -0.2e1 * t110 * t25 + 0.2e1 * t41 * t50, 0.2e1 * t110 * t26 + 0.2e1 * t41 * t51, t34 ^ 2, -0.2e1 * t34 * t33, t34 * t133, t33 * t134, t49, 0.2e1 * t11 * t50 + 0.2e1 * t22 * t33, -0.2e1 * t12 * t50 + 0.2e1 * t22 * t34, t19 ^ 2, -0.2e1 * t19 * t18, t19 * t133, t18 * t134, t49, 0.2e1 * t13 * t18 + 0.2e1 * t3 * t50, 0.2e1 * t13 * t19 - 0.2e1 * t4 * t50, 0.2e1 * t5 * t18 - 0.2e1 * t2 * t50, -0.2e1 * t1 * t18 + 0.2e1 * t2 * t19, 0.2e1 * t1 * t50 - 0.2e1 * t5 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t111, t110, t80, t52, -t53, t51 * t83, t51 * t86 - t103, -t94, -t93, 0, -pkin(2) * t50 + pkin(9) * t94 - t41 * t86, -pkin(2) * t51 + pkin(9) * t93 + t41 * t83, t34 * t101 (-t33 * t85 - t118) * t83, t101 * t50 - t34 * t86, -t103 * t82 + t33 * t86, -t113, -t11 * t86 + t39 * t50 + (pkin(9) * t33 + t120) * t83, t12 * t86 - t40 * t50 + (pkin(9) * t34 + t119) * t83, t19 * t46, -t46 * t18 - t19 * t45, -t19 * t86 + t46 * t50, t18 * t86 - t45 * t50, -t113, t13 * t45 + t16 * t50 + t58 * t18 - t3 * t86, t13 * t46 - t17 * t50 + t58 * t19 + t4 * t86, -t15 * t50 + t24 * t18 + t2 * t86 + t5 * t45, -t1 * t45 - t14 * t18 + t15 * t19 + t2 * t46, -t1 * t86 + t14 * t50 - t24 * t19 - t5 * t46, t1 * t14 + t2 * t15 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t76, t97, 0, 0, 0, pkin(2) * t128, pkin(2) * t129, t77 * t76, -0.2e1 * t76 * t105, t100 * t129, t82 * t97, t78, 0.2e1 * t123 * t76 - 0.2e1 * t39 * t86, 0.2e1 * t76 * pkin(9) * t85 + 0.2e1 * t40 * t86, t46 ^ 2, t45 * t135, t86 * t135, t45 * t128, t78, -0.2e1 * t16 * t86 + 0.2e1 * t58 * t45, 0.2e1 * t17 * t86 + 0.2e1 * t58 * t46, 0.2e1 * t15 * t86 + 0.2e1 * t24 * t45, -0.2e1 * t14 * t45 + 0.2e1 * t15 * t46, -0.2e1 * t14 * t86 - 0.2e1 * t24 * t46, t14 ^ 2 + t15 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, -t110, t25, -t26, t118, -t82 * t33 + t34 * t85, t107, t102, 0, -pkin(3) * t33 - pkin(10) * t107 - t119, -pkin(3) * t34 - pkin(10) * t102 + t120, t19 * t56, -t56 * t18 - t19 * t55, t56 * t50, -t55 * t50, 0, t13 * t55 + t71 * t18 - t117, t13 * t56 + t71 * t19 - t115, t29 * t18 + t5 * t55 - t117, -t1 * t55 - t37 * t18 + t36 * t19 + t2 * t56, -t29 * t19 - t5 * t56 + t115, t1 * t37 + t2 * t36 + t5 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t86, 0, -t73, -t86 * pkin(9), t82 * t101 (-t75 + t77) * t83, -t104, -t100, 0, -pkin(9) * t101 + (-pkin(3) * t83 + pkin(10) * t86) * t82, pkin(10) * t100 + (t123 - t124) * t83, t46 * t56, -t56 * t45 - t46 * t55, -t56 * t86, t55 * t86, 0, t71 * t45 + t58 * t55 + t116, t71 * t46 + t58 * t56 + t114, t24 * t55 + t29 * t45 + t116, -t14 * t55 + t15 * t56 + t36 * t46 - t37 * t45, -t24 * t56 - t29 * t46 - t114, t14 * t37 + t15 * t36 + t24 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t75, 0.2e1 * t105, 0, 0, 0, 0.2e1 * t124, -0.2e1 * pkin(3) * t82, t56 ^ 2, t55 * t131, 0, 0, 0, t55 * t130, t56 * t130, 0.2e1 * t29 * t55, 0.2e1 * t36 * t56 - 0.2e1 * t37 * t55, t29 * t131, t29 ^ 2 + t36 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, t50, t11, -t12, 0, 0, t19, -t18, t50, t50 * t92 + t3, -t50 * t72 - t4, t69 * t50 - t2, -t67 * t18 - t69 * t19, t67 * t50 + t1, t1 * t67 - t2 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t106, -t86, t39, -t40, 0, 0, t46, -t45, -t86, -t86 * t92 + t16, t72 * t86 - t17 (-pkin(5) - t69) * t86 + t16, -t67 * t45 - t69 * t46 (-qJ(6) - t67) * t86 + t17, t14 * t67 - t15 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t85, 0, -t82 * pkin(10), -t85 * pkin(10), 0, 0, t56, -t55, 0, -t36, -t37, -t36, -t67 * t55 - t69 * t56, t37, -t36 * t69 + t37 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t92, -0.2e1 * t72, 0.2e1 * t69, 0, 0.2e1 * t67, t67 ^ 2 + t69 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t50, t3, -t4, -t2 + t48, -pkin(5) * t19 - t18 * qJ(6), 0.2e1 * t47 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, -t86, t16, -t17, t16 - 0.2e1 * t122, -pkin(5) * t46 - t45 * qJ(6), -0.2e1 * t99 + t17, -t15 * pkin(5) + t14 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t55, 0, -t36, -t37, -t36, -pkin(5) * t56 - t55 * qJ(6), t37, -t36 * pkin(5) + t37 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t92, -t72, t89 + t92, 0, t88 + t72, t69 * pkin(5) + t67 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t89, 0, t88, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t46, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
