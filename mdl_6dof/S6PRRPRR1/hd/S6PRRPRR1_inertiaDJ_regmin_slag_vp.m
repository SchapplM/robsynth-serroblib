% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:07
% EndTime: 2019-03-08 21:54:10
% DurationCPUTime: 1.08s
% Computational Cost: add. (1921->164), mult. (4731->310), div. (0->0), fcn. (4927->12), ass. (0->117)
t85 = cos(pkin(12));
t114 = pkin(3) * t85 + pkin(4);
t140 = cos(qJ(5));
t83 = sin(pkin(12));
t142 = pkin(3) * t83;
t88 = sin(qJ(5));
t143 = -t114 * t140 + t142 * t88;
t91 = cos(qJ(6));
t82 = t91 ^ 2;
t87 = sin(qJ(6));
t129 = t87 ^ 2 - t82;
t109 = t129 * qJD(6);
t111 = qJD(5) * t140;
t126 = qJD(5) * t88;
t131 = -qJ(4) - pkin(8);
t89 = sin(qJ(3));
t74 = t131 * t89;
t92 = cos(qJ(3));
t75 = t131 * t92;
t47 = t74 * t85 + t75 * t83;
t68 = t83 * t92 + t85 * t89;
t38 = -pkin(9) * t68 + t47;
t48 = t74 * t83 - t75 * t85;
t67 = -t83 * t89 + t85 * t92;
t39 = pkin(9) * t67 + t48;
t110 = qJD(3) * t131;
t61 = t92 * qJD(4) + t110 * t89;
t62 = -t89 * qJD(4) + t110 * t92;
t37 = t61 * t85 + t62 * t83;
t98 = t68 * qJD(3);
t94 = -pkin(9) * t98 + t37;
t36 = -t61 * t83 + t62 * t85;
t123 = t92 * qJD(3);
t124 = t89 * qJD(3);
t64 = t123 * t85 - t124 * t83;
t97 = pkin(9) * t64 - t36;
t11 = t111 * t39 + t126 * t38 + t140 * t97 + t88 * t94;
t20 = -t140 * t38 + t39 * t88;
t80 = qJD(6) * t91;
t141 = t11 * t87 + t20 * t80;
t99 = t140 * t67 - t68 * t88;
t24 = qJD(5) * t99 + t140 * t64 - t88 * t98;
t44 = t140 * t68 + t67 * t88;
t139 = t44 * t24;
t138 = t44 * t87;
t137 = t44 * t91;
t84 = sin(pkin(6));
t90 = sin(qJ(2));
t136 = t84 * t90;
t93 = cos(qJ(2));
t135 = t84 * t93;
t25 = qJD(5) * t44 + t140 * t98 + t88 * t64;
t134 = t87 * t25;
t133 = t91 * t24;
t132 = t91 * t25;
t96 = t114 * t88 + t140 * t142;
t57 = t96 * qJD(5);
t59 = -pkin(5) + t143;
t130 = t57 * t87 + t59 * t80;
t128 = qJD(2) * t90;
t127 = qJD(2) * t93;
t125 = qJD(6) * t87;
t121 = -0.2e1 * pkin(2) * qJD(3);
t120 = pkin(5) * t125;
t119 = pkin(5) * t80;
t79 = pkin(3) * t124;
t118 = t93 * t124;
t117 = t84 * t128;
t116 = t84 * t127;
t115 = t87 * t80;
t78 = -pkin(3) * t92 - pkin(2);
t113 = -0.4e1 * t87 * t137;
t112 = t125 * t59 - t57 * t91;
t53 = -pkin(4) * t67 + t78;
t19 = -pkin(5) * t99 - pkin(10) * t44 + t53;
t21 = t140 * t39 + t38 * t88;
t108 = t19 * t91 - t21 * t87;
t107 = t19 * t87 + t21 * t91;
t60 = pkin(10) + t96;
t106 = -t44 * t59 - t60 * t99;
t86 = cos(pkin(6));
t65 = -t136 * t89 + t86 * t92;
t66 = t136 * t92 + t86 * t89;
t40 = t65 * t85 - t66 * t83;
t41 = t65 * t83 + t66 * t85;
t23 = t140 * t41 + t40 * t88;
t105 = t135 * t91 + t23 * t87;
t104 = t135 * t87 - t23 * t91;
t102 = t24 * t87 + t44 * t80;
t101 = t125 * t44 - t133;
t16 = -t80 * t99 + t134;
t100 = -t125 * t99 - t132;
t49 = pkin(4) * t98 + t79;
t56 = t143 * qJD(5);
t95 = t24 * t59 - t25 * t60 + t44 * t57 - t56 * t99;
t76 = 0.2e1 * t115;
t72 = -0.2e1 * t109;
t46 = -qJD(3) * t66 - t116 * t89;
t45 = qJD(3) * t65 + t116 * t92;
t42 = t44 ^ 2;
t30 = t45 * t85 + t46 * t83;
t29 = -t45 * t83 + t46 * t85;
t22 = -t140 * t40 + t41 * t88;
t17 = t20 * t125;
t14 = -t109 * t44 + t133 * t87;
t13 = t25 * pkin(5) - t24 * pkin(10) + t49;
t12 = qJD(6) * t113 - t129 * t24;
t10 = -t111 * t38 + t126 * t39 - t140 * t94 + t88 * t97;
t8 = qJD(5) * t23 - t140 * t29 + t88 * t30;
t7 = -t111 * t40 + t126 * t41 - t140 * t30 - t29 * t88;
t6 = t125 * t22 - t8 * t91;
t5 = t22 * t80 + t8 * t87;
t4 = qJD(6) * t104 + t117 * t91 + t87 * t7;
t3 = qJD(6) * t105 - t117 * t87 + t91 * t7;
t2 = -qJD(6) * t107 + t87 * t10 + t91 * t13;
t1 = -qJD(6) * t108 + t91 * t10 - t87 * t13;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t127 * t84 ^ 2 * t90 + 0.2e1 * t29 * t40 + 0.2e1 * t30 * t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t117, -t116, 0, 0, 0, 0, 0 (-t128 * t92 - t118) * t84 (-t123 * t93 + t128 * t89) * t84, -t29 * t68 + t30 * t67 - t40 * t64 - t41 * t98, t29 * t47 + t30 * t48 + t40 * t36 + t41 * t37 + (-pkin(3) * t118 + t128 * t78) * t84, 0, 0, 0, 0, 0 (-t128 * t99 - t25 * t93) * t84 (t128 * t44 - t24 * t93) * t84, 0, 0, 0, 0, 0, t102 * t22 - t105 * t25 + t138 * t8 - t4 * t99, -t101 * t22 + t104 * t25 + t137 * t8 - t3 * t99; 0, 0, 0, 0, 0.2e1 * t89 * t123, 0.2e1 * (-t89 ^ 2 + t92 ^ 2) * qJD(3), 0, 0, 0, t89 * t121, t92 * t121, -0.2e1 * t36 * t68 + 0.2e1 * t37 * t67 - 0.2e1 * t47 * t64 - 0.2e1 * t48 * t98, 0.2e1 * t36 * t47 + 0.2e1 * t37 * t48 + 0.2e1 * t78 * t79, 0.2e1 * t139, 0.2e1 * t24 * t99 - 0.2e1 * t25 * t44, 0, 0, 0, 0.2e1 * t25 * t53 - 0.2e1 * t49 * t99, 0.2e1 * t24 * t53 + 0.2e1 * t44 * t49, -0.2e1 * t115 * t42 + 0.2e1 * t139 * t82, 0.2e1 * t109 * t42 + t113 * t24, 0.2e1 * t101 * t99 + 0.2e1 * t132 * t44, 0.2e1 * t102 * t99 - 0.2e1 * t134 * t44, -0.2e1 * t99 * t25, 0.2e1 * t102 * t20 + 0.2e1 * t108 * t25 + 0.2e1 * t11 * t138 - 0.2e1 * t2 * t99, -0.2e1 * t1 * t99 - 0.2e1 * t101 * t20 - 0.2e1 * t107 * t25 + 0.2e1 * t11 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0 (t29 * t85 + t30 * t83) * pkin(3), 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, t123, -t124, 0, -pkin(8) * t123, pkin(8) * t124 (-t85 * t64 - t83 * t98) * pkin(3) (t36 * t85 + t37 * t83) * pkin(3), 0, 0, t24, -t25, 0, -t11, t10, t14, t12, t16, -t100, 0, t17 + (-qJD(6) * t106 - t11) * t91 + t95 * t87, t106 * t125 + t91 * t95 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t57, 0.2e1 * t56, t76, t72, 0, 0, 0, 0.2e1 * t112, 0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, t25, t24, 0, 0, 0, 0, 0, -t100, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25, 0, -t11, t10, t14, t12, t16, -t100, 0, t17 + (-pkin(5) * t24 - pkin(10) * t25) * t87 + (-t11 + (-pkin(5) * t44 + pkin(10) * t99) * qJD(6)) * t91, pkin(5) * t101 + pkin(10) * t100 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, t76, t72, 0, 0, 0, t112 - t120, -t119 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t72, 0, 0, 0, -0.2e1 * t120, -0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t102, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t125, 0, t56 * t87 - t60 * t80, t125 * t60 + t56 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t125, 0, -pkin(10) * t80, pkin(10) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
