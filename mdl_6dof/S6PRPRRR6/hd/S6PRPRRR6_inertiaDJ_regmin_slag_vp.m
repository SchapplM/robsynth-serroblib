% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:46
% EndTime: 2019-03-08 20:48:50
% DurationCPUTime: 1.30s
% Computational Cost: add. (893->173), mult. (2404->327), div. (0->0), fcn. (2195->10), ass. (0->128)
t137 = qJD(5) + qJD(6);
t67 = sin(qJ(6));
t68 = sin(qJ(5));
t71 = cos(qJ(6));
t72 = cos(qJ(5));
t42 = t67 * t68 - t71 * t72;
t19 = t137 * t42;
t69 = sin(qJ(4));
t140 = t19 * t69;
t43 = t67 * t72 + t71 * t68;
t73 = cos(qJ(4));
t139 = t43 * t73;
t86 = t69 * pkin(4) - t73 * pkin(9);
t46 = qJ(3) + t86;
t38 = t68 * t46;
t127 = t69 * t72;
t75 = -pkin(2) - pkin(8);
t55 = t75 * t127;
t138 = -t38 - t55;
t115 = qJD(5) * t73;
t105 = t68 * t115;
t112 = t69 * qJD(4);
t98 = t72 * t112;
t34 = -t98 - t105;
t63 = t72 ^ 2;
t123 = t68 ^ 2 - t63;
t91 = t123 * qJD(5);
t136 = 2 * qJD(3);
t135 = pkin(9) + pkin(10);
t114 = qJD(5) * t75;
t104 = t68 * t114;
t116 = qJD(5) * t72;
t87 = pkin(4) * t73 + pkin(9) * t69;
t40 = t87 * qJD(4) + qJD(3);
t60 = t73 * qJD(4);
t96 = t75 * t60;
t16 = t69 * t104 - t46 * t116 - t68 * t40 - t72 * t96;
t101 = t68 * t112;
t103 = t72 * t115;
t36 = t101 - t103;
t9 = t36 * pkin(10) - t16;
t134 = t71 * t9;
t20 = t137 * t43;
t133 = t20 * t69;
t65 = sin(pkin(6));
t70 = sin(qJ(2));
t132 = t65 * t70;
t74 = cos(qJ(2));
t131 = t65 * t74;
t129 = t68 * t73;
t21 = -pkin(10) * t129 - t138;
t130 = t67 * t21;
t128 = t68 * t75;
t126 = t71 * t21;
t125 = t72 * t73;
t124 = t73 * t75;
t62 = t69 ^ 2;
t64 = t73 ^ 2;
t122 = t62 - t64;
t121 = t62 + t64;
t120 = qJD(2) * t74;
t66 = cos(pkin(6));
t31 = t73 * t131 + t66 * t69;
t119 = qJD(4) * t31;
t118 = qJD(4) * t72;
t117 = qJD(5) * t68;
t113 = qJD(6) * t67;
t111 = qJ(3) * qJD(4);
t110 = -0.2e1 * pkin(4) * qJD(5);
t109 = pkin(5) * t117;
t108 = pkin(5) * t60;
t107 = pkin(5) * t113;
t106 = qJD(6) * t71 * pkin(5);
t57 = qJD(2) * t132;
t102 = t65 * t120;
t100 = t68 * t116;
t99 = t75 * t112;
t97 = t69 * t60;
t30 = t72 * t40;
t93 = pkin(5) - t128;
t8 = t30 + (-t55 + (pkin(10) * t73 - t46) * t68) * qJD(5) + (pkin(10) * t127 + t93 * t73) * qJD(4);
t95 = -t67 * t9 + t71 * t8;
t39 = t72 * t46;
t18 = -pkin(10) * t125 + t93 * t69 + t39;
t94 = -t69 * pkin(5) - t18;
t92 = qJD(5) * t135;
t90 = t122 * qJD(4);
t89 = t68 * t96;
t88 = t68 * t98;
t85 = t71 * t18 - t130;
t84 = t67 * t18 + t126;
t32 = -t69 * t131 + t66 * t73;
t24 = t72 * t132 - t32 * t68;
t25 = t68 * t132 + t32 * t72;
t83 = t71 * t24 - t67 * t25;
t82 = t67 * t24 + t71 * t25;
t53 = t135 * t68;
t54 = t135 * t72;
t81 = -t71 * t53 - t67 * t54;
t80 = -t67 * t53 + t71 * t54;
t23 = t32 * qJD(4) - t73 * t57;
t79 = t31 * t116 + t23 * t68;
t78 = t31 * t117 - t23 * t72;
t28 = t42 * t73;
t76 = qJD(4) * t42;
t59 = -t72 * pkin(5) - pkin(4);
t56 = 0.2e1 * t97;
t45 = t72 * t92;
t44 = t68 * t92;
t41 = (pkin(5) * t68 - t75) * t73;
t35 = t69 * t116 + t68 * t60;
t33 = t69 * t117 - t72 * t60;
t26 = -t36 * pkin(5) + t99;
t22 = -t69 * t57 + t119;
t17 = t138 * qJD(5) + t30 - t89;
t15 = -t80 * qJD(6) + t67 * t44 - t71 * t45;
t14 = -t81 * qJD(6) + t71 * t44 + t67 * t45;
t13 = -t113 * t129 + (t137 * t125 - t101) * t71 + t34 * t67;
t12 = -qJD(4) * t139 + t140;
t11 = -t137 * t139 + t69 * t76;
t10 = t73 * t76 + t133;
t7 = t24 * qJD(5) + t68 * t102 - t22 * t72;
t6 = -t25 * qJD(5) + t72 * t102 + t22 * t68;
t4 = -t84 * qJD(6) + t95;
t3 = -t85 * qJD(6) - t67 * t8 - t134;
t2 = -t82 * qJD(6) + t71 * t6 - t67 * t7;
t1 = -t83 * qJD(6) - t67 * t6 - t71 * t7;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t57, -t102, t57, t102 (qJD(3) * t70 + (-pkin(2) * t70 + qJ(3) * t74) * qJD(2)) * t65, 0, 0, 0, 0, 0 (t69 * t120 + t70 * t60) * t65 (-t70 * t112 + t73 * t120) * t65, 0, 0, 0, 0, 0 (-t68 * t119 + t6) * t69 + (qJD(4) * t24 + t79) * t73 (-t31 * t118 - t7) * t69 + (-qJD(4) * t25 - t78) * t73, 0, 0, 0, 0, 0, t31 * t13 + t139 * t23 + t2 * t69 + t83 * t60, t1 * t69 + t31 * t11 - t23 * t28 - t82 * t60; 0, 0, 0, 0, 0, t136, qJ(3) * t136, -0.2e1 * t97, 0.2e1 * t90, 0, 0, 0, 0.2e1 * qJD(3) * t69 + 0.2e1 * t73 * t111, 0.2e1 * qJD(3) * t73 - 0.2e1 * t69 * t111, -0.2e1 * t64 * t100 - 0.2e1 * t63 * t97, 0.2e1 * t64 * t91 + 0.4e1 * t73 * t88, -0.2e1 * t69 * t105 - 0.2e1 * t122 * t118, -0.2e1 * t69 * t103 + 0.2e1 * t68 * t90, t56, -0.2e1 * t64 * t72 * t114 + 0.2e1 * t39 * t60 + 0.2e1 * (t17 + t89) * t69, 0.2e1 * t64 * t104 + 0.2e1 * t16 * t69 + 0.2e1 * (-t38 + t55) * t60, -0.2e1 * t28 * t11, -0.2e1 * t11 * t139 + 0.2e1 * t28 * t13, 0.2e1 * t11 * t69 - 0.2e1 * t28 * t60, -0.2e1 * t13 * t69 - 0.2e1 * t139 * t60, t56, 0.2e1 * t41 * t13 + 0.2e1 * t139 * t26 + 0.2e1 * t4 * t69 + 0.2e1 * t85 * t60, 0.2e1 * t41 * t11 - 0.2e1 * t26 * t28 + 0.2e1 * t3 * t69 - 0.2e1 * t84 * t60; 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t116, t121 * t117, 0, 0, 0, 0, 0, t12 * t69 - t73 * t13, t10 * t69 - t73 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0, 0, 0, 0, 0, t78, t79, 0, 0, 0, 0, 0, t31 * t20 + t23 * t42, -t31 * t19 + t23 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t60, 0, -t99, -t96, -t73 * t91 - t88, -0.4e1 * t73 * t100 + t123 * t112, t35, -t33, 0 (-t68 * t124 - t87 * t72) * qJD(5) + (t86 * t68 - t55) * qJD(4) (-t72 * t124 + t87 * t68) * qJD(5) + (-pkin(9) * t125 + (pkin(4) * t72 + t128) * t69) * qJD(4), t11 * t43 + t28 * t19, -t11 * t42 - t43 * t13 + t139 * t19 + t28 * t20, t43 * t60 - t140, -t42 * t60 - t133, 0, t109 * t139 + t59 * t13 + t15 * t69 + t41 * t20 + t26 * t42 + t81 * t60, -t28 * t109 + t59 * t11 + t14 * t69 - t41 * t19 + t26 * t43 - t80 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t60, 0, 0, 0, 0, 0, t34, t36, 0, 0, 0, 0, 0, t42 * t112 - t73 * t20, t43 * t112 + t73 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100, -0.2e1 * t91, 0, 0, 0, t68 * t110, t72 * t110, -0.2e1 * t43 * t19, 0.2e1 * t19 * t42 - 0.2e1 * t43 * t20, 0, 0, 0, 0.2e1 * t42 * t109 + 0.2e1 * t59 * t20, 0.2e1 * t43 * t109 - 0.2e1 * t59 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, t60, t17, t16, 0, 0, t11, -t13, t60, t71 * t108 + (t67 * t94 - t126) * qJD(6) + t95, -t134 + (-t8 - t108) * t67 + (t71 * t94 + t130) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t33, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, -t117, 0, -pkin(9) * t116, pkin(9) * t117, 0, 0, -t19, -t20, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t107, -0.2e1 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t13, t60, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
