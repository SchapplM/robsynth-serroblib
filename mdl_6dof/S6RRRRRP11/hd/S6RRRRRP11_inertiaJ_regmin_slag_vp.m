% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:07:50
% EndTime: 2019-05-08 07:07:59
% DurationCPUTime: 2.08s
% Computational Cost: add. (3588->235), mult. (9523->480), div. (0->0), fcn. (10953->12), ass. (0->127)
t81 = sin(pkin(6));
t91 = cos(qJ(2));
t117 = t81 * t91;
t87 = sin(qJ(2));
t134 = pkin(1) * t87;
t83 = cos(pkin(6));
t60 = pkin(9) * t117 + t83 * t134;
t80 = sin(pkin(7));
t82 = cos(pkin(7));
t115 = t82 * t91;
t95 = t81 * t115;
t36 = (t80 * t83 + t95) * pkin(10) + t60;
t86 = sin(qJ(3));
t90 = cos(qJ(3));
t118 = t81 * t87;
t71 = t83 * t91 * pkin(1);
t40 = t83 * pkin(2) + t71 + (-pkin(10) * t82 - pkin(9)) * t118;
t46 = (-pkin(10) * t80 * t87 - pkin(2) * t91 - pkin(1)) * t81;
t92 = t40 * t82 + t46 * t80;
t20 = -t86 * t36 + t92 * t90;
t120 = t80 * t86;
t39 = t83 * t120 + (t86 * t115 + t87 * t90) * t81;
t53 = t80 * t117 - t83 * t82;
t85 = sin(qJ(4));
t89 = cos(qJ(4));
t27 = t39 * t85 + t53 * t89;
t141 = -0.2e1 * t27;
t119 = t80 * t90;
t38 = t86 * t118 - t83 * t119 - t90 * t95;
t140 = 0.2e1 * t38;
t139 = -0.2e1 * t39;
t56 = t89 * t120 + t85 * t82;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t41 = t88 * t119 + t84 * t56;
t138 = -0.2e1 * t41;
t137 = -0.2e1 * t56;
t136 = -0.2e1 * t85;
t135 = 0.2e1 * t89;
t133 = pkin(2) * t86;
t132 = pkin(2) * t90;
t131 = pkin(4) * t88;
t130 = pkin(11) * t84;
t25 = -t80 * t40 + t82 * t46;
t14 = t38 * pkin(3) - t39 * pkin(11) + t25;
t21 = t90 * t36 + t92 * t86;
t17 = -t53 * pkin(11) + t21;
t8 = t89 * t14 - t85 * t17;
t6 = -t38 * pkin(4) - t8;
t129 = t6 * t84;
t128 = t6 * t88;
t127 = t84 * pkin(5);
t28 = t39 * t89 - t53 * t85;
t23 = t28 * t88 + t38 * t84;
t126 = t23 * t84;
t96 = pkin(10) * t119;
t51 = t96 + (pkin(11) + t133) * t82;
t52 = (-pkin(3) * t90 - pkin(11) * t86 - pkin(2)) * t80;
t32 = -t85 * t51 + t89 * t52;
t30 = pkin(4) * t119 - t32;
t125 = t30 * t84;
t124 = t30 * t88;
t42 = -t84 * t119 + t88 * t56;
t123 = t42 * t84;
t75 = t80 ^ 2;
t122 = t75 * t90;
t76 = t81 ^ 2;
t121 = t76 * t91;
t116 = t82 * t86;
t114 = t84 * t27;
t55 = t85 * t120 - t89 * t82;
t113 = t84 * t55;
t112 = t84 * t88;
t111 = t84 * t89;
t110 = t85 * t27;
t109 = t85 * t38;
t108 = t85 * t55;
t107 = t88 * t27;
t106 = t88 * t55;
t105 = t88 * t85;
t104 = t88 * t89;
t103 = t89 * t38;
t102 = -qJ(6) - pkin(12);
t101 = qJ(6) * t85;
t100 = 0.2e1 * t119;
t99 = 0.2e1 * t81 * t83;
t98 = t85 * t135;
t97 = pkin(11) * t104;
t94 = t85 * t119;
t93 = t89 * t119;
t16 = t53 * pkin(3) - t20;
t11 = t27 * pkin(4) - t28 * pkin(12) + t16;
t9 = t85 * t14 + t89 * t17;
t7 = t38 * pkin(12) + t9;
t3 = t88 * t11 - t84 * t7;
t69 = pkin(10) * t120;
t50 = t69 + (-pkin(3) - t132) * t82;
t29 = t55 * pkin(4) - t56 * pkin(12) + t50;
t33 = t89 * t51 + t85 * t52;
t31 = -pkin(12) * t119 + t33;
t18 = t88 * t29 - t84 * t31;
t4 = t84 * t11 + t88 * t7;
t19 = t84 * t29 + t88 * t31;
t79 = t88 ^ 2;
t78 = t85 ^ 2;
t77 = t84 ^ 2;
t74 = -t88 * pkin(5) - pkin(4);
t66 = t102 * t88;
t65 = t102 * t84;
t64 = -t89 * pkin(4) - t85 * pkin(12) - pkin(3);
t63 = (pkin(11) + t127) * t85;
t61 = t88 * t64;
t59 = pkin(2) * t116 + t96;
t58 = -pkin(9) * t118 + t71;
t57 = t82 * t132 - t69;
t48 = t84 * t64 + t97;
t47 = -pkin(11) * t111 + t61;
t43 = t97 + (t64 - t101) * t84;
t37 = -t88 * t101 + t61 + (-pkin(5) - t130) * t89;
t24 = t41 * pkin(5) + t30;
t22 = t28 * t84 - t38 * t88;
t13 = -t41 * qJ(6) + t19;
t12 = t55 * pkin(5) - t42 * qJ(6) + t18;
t5 = t22 * pkin(5) + t6;
t2 = -t22 * qJ(6) + t4;
t1 = t27 * pkin(5) - t23 * qJ(6) + t3;
t10 = [1, 0, 0, t76 * t87 ^ 2, 0.2e1 * t87 * t121, t87 * t99, t91 * t99, t83 ^ 2, 0.2e1 * pkin(1) * t121 + 0.2e1 * t58 * t83, -0.2e1 * t76 * t134 - 0.2e1 * t60 * t83, t39 ^ 2, t38 * t139, t53 * t139, t53 * t140, t53 ^ 2, -0.2e1 * t20 * t53 + 0.2e1 * t25 * t38, 0.2e1 * t21 * t53 + 0.2e1 * t25 * t39, t28 ^ 2, t28 * t141, t28 * t140, t38 * t141, t38 ^ 2, 0.2e1 * t16 * t27 + 0.2e1 * t8 * t38, 0.2e1 * t16 * t28 - 0.2e1 * t9 * t38, t23 ^ 2, -0.2e1 * t23 * t22, 0.2e1 * t23 * t27, t22 * t141, t27 ^ 2, 0.2e1 * t6 * t22 + 0.2e1 * t3 * t27, 0.2e1 * t6 * t23 - 0.2e1 * t4 * t27, -0.2e1 * t1 * t23 - 0.2e1 * t2 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t118, t117, t83, t58, -t60, t39 * t120 (-t38 * t86 + t39 * t90) * t80, -t53 * t120 + t39 * t82, -t53 * t119 - t38 * t82, -t53 * t82, t20 * t82 - t57 * t53 + (-pkin(2) * t38 - t25 * t90) * t80, -t21 * t82 + t59 * t53 + (-pkin(2) * t39 + t25 * t86) * t80, t28 * t56, -t56 * t27 - t28 * t55, -t28 * t119 + t56 * t38, t27 * t119 - t55 * t38, -t38 * t119, -t8 * t119 + t16 * t55 + t50 * t27 + t32 * t38, t9 * t119 + t16 * t56 + t50 * t28 - t33 * t38, t23 * t42, -t42 * t22 - t23 * t41, t23 * t55 + t42 * t27, -t22 * t55 - t41 * t27, t27 * t55, t18 * t27 + t30 * t22 + t3 * t55 + t6 * t41, -t19 * t27 + t30 * t23 - t4 * t55 + t6 * t42, -t1 * t42 - t12 * t23 - t13 * t22 - t2 * t41, t1 * t12 + t2 * t13 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t75 * t86 ^ 2, 0.2e1 * t86 * t122, 0.2e1 * t80 * t116, t82 * t100, t82 ^ 2, 0.2e1 * pkin(2) * t122 + 0.2e1 * t57 * t82, -0.2e1 * t75 * t133 - 0.2e1 * t59 * t82, t56 ^ 2, t55 * t137, t119 * t137, t55 * t100, t75 * t90 ^ 2, -0.2e1 * t32 * t119 + 0.2e1 * t50 * t55, 0.2e1 * t33 * t119 + 0.2e1 * t50 * t56, t42 ^ 2, t42 * t138, 0.2e1 * t42 * t55, t55 * t138, t55 ^ 2, 0.2e1 * t18 * t55 + 0.2e1 * t30 * t41, -0.2e1 * t19 * t55 + 0.2e1 * t30 * t42, -0.2e1 * t12 * t42 - 0.2e1 * t13 * t41, t12 ^ 2 + t13 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, -t53, t20, -t21, t28 * t85, t28 * t89 - t110, t109, t103, 0, -pkin(3) * t27 - pkin(11) * t109 - t16 * t89, -pkin(3) * t28 - pkin(11) * t103 + t16 * t85, t23 * t105 (-t22 * t88 - t126) * t85, t27 * t105 - t23 * t89, -t110 * t84 + t22 * t89, -t27 * t89, t47 * t27 - t3 * t89 + (pkin(11) * t22 + t129) * t85, -t48 * t27 + t4 * t89 + (pkin(11) * t23 + t128) * t85, -t43 * t22 - t37 * t23 + (-t1 * t88 - t2 * t84) * t85, t1 * t37 + t2 * t43 + t5 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, t82, t57, -t59, t56 * t85, t56 * t89 - t108, -t94, -t93, 0, -pkin(3) * t55 + pkin(11) * t94 - t50 * t89, -pkin(3) * t56 + pkin(11) * t93 + t50 * t85, t42 * t105 (-t41 * t88 - t123) * t85, t55 * t105 - t42 * t89, -t108 * t84 + t41 * t89, -t55 * t89, -t18 * t89 + t47 * t55 + (pkin(11) * t41 + t125) * t85, t19 * t89 - t48 * t55 + (pkin(11) * t42 + t124) * t85, -t37 * t42 - t43 * t41 + (-t12 * t88 - t13 * t84) * t85, t12 * t37 + t13 * t43 + t24 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t78, t98, 0, 0, 0, pkin(3) * t135, pkin(3) * t136, t79 * t78, -0.2e1 * t78 * t112, t104 * t136, t84 * t98, t89 ^ 2, 0.2e1 * t130 * t78 - 0.2e1 * t47 * t89, 0.2e1 * t78 * pkin(11) * t88 + 0.2e1 * t48 * t89, 0.2e1 * (-t37 * t88 - t43 * t84) * t85, t37 ^ 2 + t43 ^ 2 + t63 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t38, t8, -t9, t126, -t84 * t22 + t23 * t88, t114, t107, 0, -pkin(4) * t22 - pkin(12) * t114 - t128, -pkin(4) * t23 - pkin(12) * t107 + t129, -t1 * t84 + t2 * t88 + t66 * t22 - t65 * t23, t1 * t65 - t2 * t66 + t5 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t55, -t119, t32, -t33, t123, -t84 * t41 + t42 * t88, t113, t106, 0, -pkin(4) * t41 - pkin(12) * t113 - t124, -pkin(4) * t42 - pkin(12) * t106 + t125, -t12 * t84 + t13 * t88 + t66 * t41 - t65 * t42, t12 * t65 - t13 * t66 + t24 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t89, 0, -t85 * pkin(11), -t89 * pkin(11), t84 * t105 (-t77 + t79) * t85, -t111, -t104, 0, -pkin(11) * t105 + (-pkin(4) * t85 + pkin(12) * t89) * t84, pkin(12) * t104 + (t130 - t131) * t85 (-t65 * t85 + t43) * t88 + (t66 * t85 - t37) * t84, t37 * t65 - t43 * t66 + t63 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t77, 0.2e1 * t112, 0, 0, 0, 0.2e1 * t131, -0.2e1 * pkin(4) * t84, -0.2e1 * t65 * t84 - 0.2e1 * t66 * t88, t65 ^ 2 + t66 ^ 2 + t74 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t27, t3, -t4, -pkin(5) * t23, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, t55, t18, -t19, -pkin(5) * t42, t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t84 * t85, -t89, t47, -t48, -pkin(5) * t105, t37 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t88, 0, -t84 * pkin(12), -t88 * pkin(12), -t127, t65 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
