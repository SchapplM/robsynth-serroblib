% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:52
% EndTime: 2019-03-09 04:55:57
% DurationCPUTime: 1.61s
% Computational Cost: add. (1198->241), mult. (2508->375), div. (0->0), fcn. (1765->4), ass. (0->132)
t120 = qJ(5) * qJD(3);
t69 = cos(qJ(3));
t100 = t69 * t120;
t67 = sin(qJ(3));
t153 = t67 * qJD(5) + t100;
t65 = -pkin(4) - qJ(6);
t152 = t65 * t69;
t151 = 0.2e1 * t153;
t68 = cos(qJ(4));
t59 = qJD(5) * t68;
t66 = sin(qJ(4));
t131 = t66 * qJ(5);
t83 = -pkin(4) * t68 - t131;
t150 = qJD(4) * t83 + t59;
t149 = t65 * t68 - t131;
t61 = t66 ^ 2;
t63 = t68 ^ 2;
t135 = t61 - t63;
t95 = t135 * qJD(4);
t145 = pkin(8) * t67;
t35 = -pkin(3) + t83;
t124 = t67 * qJD(3);
t107 = t66 * t124;
t119 = qJ(5) * qJD(4);
t101 = t66 * t119;
t102 = t67 * t120;
t127 = qJD(4) * t69;
t109 = t68 * t127;
t70 = -pkin(1) - pkin(7);
t52 = t70 * t124;
t91 = pkin(4) * t109 + t101 * t69 + t102 * t68 + t52;
t8 = -pkin(4) * t107 - t59 * t69 + t91;
t148 = (-t35 * t69 + t145) * qJD(4) + t8;
t147 = 2 * qJD(2);
t146 = pkin(5) + pkin(8);
t144 = t69 * pkin(8);
t143 = t35 * t67;
t141 = t66 * t69;
t140 = t67 * t70;
t139 = t68 * t69;
t125 = t66 * qJD(5);
t128 = qJD(4) * t66;
t55 = pkin(4) * t128;
t17 = -t119 * t68 - t125 + t55;
t138 = t69 * t17;
t137 = t69 * t70;
t136 = t153 * t68;
t134 = t61 + t63;
t62 = t67 ^ 2;
t64 = t69 ^ 2;
t133 = t62 - t64;
t132 = qJ(5) * t68;
t130 = qJD(3) * t66;
t129 = qJD(3) * t68;
t60 = qJD(4) * t68;
t126 = qJD(4) * t70;
t122 = t67 * qJD(6);
t58 = t69 * qJD(3);
t121 = qJ(2) * qJD(3);
t118 = qJ(5) * qJD(5);
t117 = pkin(5) * t67 * t68;
t116 = -0.2e1 * pkin(3) * qJD(4);
t45 = t66 * t140;
t46 = t68 * t140;
t103 = t70 * t58;
t88 = pkin(3) * t69 + t145;
t31 = qJD(3) * t88 + qJD(2);
t87 = pkin(3) * t67 - t144;
t34 = qJ(2) + t87;
t115 = t103 * t68 + t31 * t66 + t34 * t60;
t114 = pkin(4) * t58;
t113 = pkin(8) * t128;
t112 = t66 * t127;
t111 = t66 * t126;
t110 = t67 * t60;
t108 = t68 * t126;
t106 = t66 * t60;
t105 = t68 * t124;
t104 = t67 * t58;
t99 = t65 * t67;
t98 = t134 * t69;
t97 = -t68 * t34 + t45;
t96 = qJD(4) * (t62 + t64);
t94 = t133 * qJD(3);
t93 = 0.2e1 * t104;
t92 = t103 * t66 + t108 * t67 + t128 * t34 - t31 * t68;
t90 = t66 * t105;
t29 = -pkin(3) + t149;
t3 = (qJ(6) * qJD(4) - qJD(5)) * t139 + (qJD(3) * t99 + qJD(6) * t69) * t66 + t91;
t89 = t127 * t29 - t3;
t30 = t66 * t34;
t14 = -qJ(5) * t67 - t30 - t46;
t10 = -pkin(5) * t141 - t14;
t9 = t45 + (pkin(5) * t69 - t34) * t68 + t99;
t85 = t10 * t68 + t66 * t9;
t84 = t10 * t66 - t68 * t9;
t15 = -t67 * pkin(4) + t97;
t82 = t14 * t68 - t15 * t66;
t81 = t14 * t66 + t15 * t68;
t80 = qJ(6) * t66 - t132;
t79 = -t68 * qJD(6) - t125;
t25 = t58 * t66 + t110;
t6 = t111 * t67 - t115;
t77 = -pkin(5) * t112 + t92;
t12 = qJD(4) * t80 + t55 + t79;
t53 = pkin(4) * t141;
t13 = t53 + (-t70 + t80) * t69;
t76 = -qJD(4) * t13 - t69 * t12 + t124 * t29;
t4 = -t100 + (-qJD(5) + t111) * t67 - t115;
t5 = t92 - t114;
t75 = qJD(4) * t81 - t4 * t68 + t5 * t66;
t32 = t146 * t128;
t56 = pkin(8) * t60;
t33 = pkin(5) * t60 + t56;
t39 = t146 * t66;
t40 = t146 * t68;
t74 = -t32 * t68 + t33 * t66 + (t39 * t68 - t40 * t66) * qJD(4);
t16 = t53 + (-t70 - t132) * t69;
t73 = -qJD(4) * t16 - t138 + (t143 + t144) * qJD(3);
t72 = pkin(5) * t107 + (-pkin(5) * t139 - t45) * qJD(4) + t115;
t71 = 0.2e1 * qJD(5);
t44 = t66 * t102;
t26 = -t107 + t109;
t24 = t105 + t112;
t23 = t128 * t67 - t58 * t68;
t22 = t68 * t96;
t21 = t66 * t96;
t20 = qJD(3) * t98;
t11 = (-0.1e1 + t134) * t93;
t2 = t72 + t153;
t1 = -t122 + (-t117 + t152) * qJD(3) + t77;
t7 = [0, 0, 0, 0, t147, qJ(2) * t147, -0.2e1 * t104, 0.2e1 * t94, 0, 0, 0, 0.2e1 * qJD(2) * t67 + 0.2e1 * t121 * t69, 0.2e1 * qJD(2) * t69 - 0.2e1 * t121 * t67, -0.2e1 * t104 * t63 - 0.2e1 * t106 * t64, 0.2e1 * t64 * t95 + 0.4e1 * t69 * t90, -0.2e1 * t112 * t67 - 0.2e1 * t129 * t133, -0.2e1 * t109 * t67 + 0.2e1 * t66 * t94, t93, -0.2e1 * t64 * t108 - 0.2e1 * t92 * t67 + 0.2e1 * (-t97 + 0.2e1 * t45) * t58, 0.2e1 * t64 * t111 + 0.2e1 * t6 * t67 + 0.2e1 * (-t30 + t46) * t58, -0.2e1 * t81 * t124 + 0.2e1 * (qJD(4) * t82 + t4 * t66 + t5 * t68) * t69, 0.2e1 * (t130 * t16 + t5) * t67 + 0.2e1 * (qJD(3) * t15 - t16 * t60 - t8 * t66) * t69, 0.2e1 * (t129 * t16 - t4) * t67 + 0.2e1 * (-qJD(3) * t14 + t128 * t16 - t8 * t68) * t69, 0.2e1 * t14 * t4 + 0.2e1 * t15 * t5 + 0.2e1 * t16 * t8, 0.2e1 * t84 * t124 + 0.2e1 * (-qJD(4) * t85 + t1 * t68 - t2 * t66) * t69, 0.2e1 * (t129 * t13 + t2) * t67 + 0.2e1 * (qJD(3) * t10 + t128 * t13 - t3 * t68) * t69, 0.2e1 * (-t13 * t130 - t1) * t67 + 0.2e1 * (-qJD(3) * t9 + t13 * t60 + t3 * t66) * t69, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t13 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, t22, -t21 (-qJD(3) * t82 - t8) * t69 + (qJD(3) * t16 + t75) * t67, 0, -t21, -t22 (qJD(3) * t85 - t3) * t69 + (qJD(3) * t13 - qJD(4) * t84 + t1 * t66 + t2 * t68) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t58, 0, -t52, -t103, -t69 * t95 - t90, -0.4e1 * t106 * t69 + t124 * t135, t25, -t23, 0 (-t137 * t66 - t68 * t88) * qJD(4) + (t66 * t87 - t46) * qJD(3) (-t137 * t68 + t66 * t88) * qJD(4) + (-pkin(8) * t139 + (pkin(3) * t68 + t66 * t70) * t67) * qJD(3), t75, t148 * t68 + t66 * t73, -t148 * t66 + t68 * t73, pkin(8) * t75 + t16 * t17 + t8 * t35 (-t39 * t124 + t33 * t69 + t2 + (-t40 * t69 + t9) * qJD(4)) * t68 + (t40 * t124 + t32 * t69 + t1 + (-t39 * t69 - t10) * qJD(4)) * t66, -t32 * t67 + t40 * t58 + t66 * t89 + t68 * t76, -t33 * t67 - t39 * t58 - t66 * t76 + t68 * t89, t1 * t39 - t10 * t32 + t12 * t13 + t2 * t40 + t29 * t3 + t33 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t58, 0, 0, 0, 0, 0, -t24, -t26, t20, t24, t26, -t138 + (pkin(8) * t98 + t143) * qJD(3), t20, t26, -t24 (-t12 + (t39 * t66 + t40 * t68) * qJD(3)) * t69 + (qJD(3) * t29 + t74) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t106, -0.2e1 * t95, 0, 0, 0, t66 * t116, t68 * t116, 0, -0.2e1 * t128 * t35 + 0.2e1 * t17 * t68, -0.2e1 * t17 * t66 - 0.2e1 * t35 * t60, 0.2e1 * t35 * t17, 0.2e1 * t74, -0.2e1 * t12 * t66 - 0.2e1 * t29 * t60, -0.2e1 * t12 * t68 + 0.2e1 * t128 * t29, 0.2e1 * t12 * t29 - 0.2e1 * t32 * t40 + 0.2e1 * t33 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, t58, -t92, t6, pkin(4) * t105 + t44 + (-t125 + (pkin(4) * t66 - t132) * qJD(4)) * t69, t92 - 0.2e1 * t114, -t6 + t151, -pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t14, -t65 * t105 + t44 + ((-t65 * t66 - t132) * qJD(4) + t79) * t69, t72 + t151, 0.2e1 * t122 + (t117 - 0.2e1 * t152) * qJD(3) - t77, qJ(5) * t2 + qJD(5) * t10 - qJD(6) * t9 + t1 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t23, 0, t25, -t23, -pkin(4) * t25 - t101 * t67 + t136, 0, -t23, -t25, t65 * t110 + (t65 * t58 + (-qJD(6) - t119) * t67) * t66 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t128, 0, -t56, t113, t150, t56, -t113, t150 * pkin(8), qJD(4) * t149 - qJD(6) * t66 + t59, -t32, -t33, -qJ(5) * t32 + qJD(5) * t40 - qJD(6) * t39 + t33 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0.2e1 * t118, 0, t71, 0.2e1 * qJD(6), -0.2e1 * qJD(6) * t65 + 0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t58, 0, t5, -t24, 0, -t58, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, t56, t60, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t58, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, 0, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
