% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:21
% EndTime: 2019-12-31 22:12:28
% DurationCPUTime: 1.79s
% Computational Cost: add. (2085->248), mult. (5840->478), div. (0->0), fcn. (5308->8), ass. (0->128)
t88 = sin(qJ(3));
t160 = -0.4e1 * t88;
t85 = sin(pkin(5));
t89 = sin(qJ(2));
t153 = t85 * t89;
t86 = cos(pkin(5));
t92 = cos(qJ(2));
t50 = pkin(7) * t153 + (-pkin(1) * t92 - pkin(2)) * t86;
t91 = cos(qJ(3));
t58 = t88 * t153 - t86 * t91;
t59 = t91 * t153 + t86 * t88;
t26 = t58 * pkin(3) - t59 * pkin(9) + t50;
t152 = t85 * t92;
t129 = pkin(7) * t152;
t157 = pkin(1) * t89;
t51 = t129 + (pkin(8) + t157) * t86;
t52 = (-pkin(2) * t92 - pkin(8) * t89 - pkin(1)) * t85;
t148 = t91 * t51 + t88 * t52;
t28 = -pkin(9) * t152 + t148;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t9 = t87 * t26 + t90 * t28;
t107 = -pkin(3) * t91 - pkin(9) * t88;
t69 = -pkin(2) + t107;
t150 = t90 * t91;
t77 = pkin(8) * t150;
t48 = t87 * t69 + t77;
t140 = qJD(2) * t92;
t117 = t85 * t140;
t139 = qJD(3) * t58;
t33 = t91 * t117 - t139;
t159 = -qJD(4) * t152 + t33;
t83 = t90 ^ 2;
t145 = t87 ^ 2 - t83;
t112 = t145 * qJD(4);
t158 = 0.2e1 * t85;
t156 = pkin(8) * t85;
t155 = pkin(8) * t87;
t141 = qJD(2) * t89;
t118 = t85 * t141;
t134 = qJD(4) * t87;
t20 = t87 * t118 - t59 * t134 + t159 * t90;
t154 = t20 * t87;
t151 = t88 * t90;
t149 = -qJ(5) - pkin(9);
t133 = qJD(4) * t90;
t106 = pkin(3) * t88 - pkin(9) * t91;
t99 = t106 * qJD(3);
t147 = -t69 * t133 - t87 * t99;
t138 = qJD(3) * t88;
t123 = t87 * t138;
t146 = pkin(8) * t123 + t90 * t99;
t82 = t88 ^ 2;
t144 = -t91 ^ 2 + t82;
t143 = qJ(5) * t88;
t142 = qJ(5) * t90;
t137 = qJD(3) * t90;
t136 = qJD(3) * t91;
t135 = qJD(3) * t92;
t132 = qJD(4) * t91;
t131 = t90 * qJD(5);
t130 = t91 * t155;
t128 = -0.2e1 * pkin(2) * qJD(3);
t127 = -0.2e1 * pkin(3) * qJD(4);
t126 = pkin(8) * t136;
t125 = pkin(4) * t134;
t80 = t85 ^ 2;
t124 = t80 * t140;
t122 = t90 * t136;
t120 = t87 * t132;
t119 = t90 * t132;
t116 = t87 * t133;
t115 = t88 * t136;
t8 = t90 * t26 - t28 * t87;
t114 = -t88 * t51 + t52 * t91;
t113 = qJD(4) * t149;
t111 = t144 * qJD(3);
t110 = 0.2e1 * t115;
t109 = t89 * t124;
t108 = t87 * t122;
t27 = pkin(3) * t152 - t114;
t35 = -t87 * t152 + t59 * t90;
t6 = pkin(4) * t58 - qJ(5) * t35 + t8;
t34 = t90 * t152 + t59 * t87;
t7 = -qJ(5) * t34 + t9;
t105 = -t6 * t90 - t7 * t87;
t104 = -t34 * t90 - t35 * t87;
t53 = (pkin(2) * t89 - pkin(8) * t92) * t85 * qJD(2);
t54 = -t86 * pkin(1) * t140 + pkin(7) * t118;
t14 = -t51 * t136 - t52 * t138 + t53 * t91 + t88 * t54;
t12 = -pkin(3) * t118 - t14;
t103 = t12 * t87 + t27 * t133;
t102 = -t12 * t90 + t27 * t134;
t32 = t59 * qJD(3) + t88 * t117;
t101 = t58 * t133 + t32 * t87;
t100 = t58 * t134 - t32 * t90;
t13 = -t52 * t136 + t51 * t138 - t88 * t53 + t91 * t54;
t11 = pkin(9) * t118 - t13;
t55 = (t86 * t157 + t129) * qJD(2);
t93 = t32 * pkin(3) - t33 * pkin(9) + t55;
t3 = -t90 * t11 - t26 * t133 + t28 * t134 - t87 * t93;
t98 = t88 * t135 + t91 * t141;
t97 = -t91 * t135 + t88 * t141;
t96 = -t88 * t134 + t122;
t95 = t88 * t137 + t120;
t94 = t88 * t133 + t87 * t136;
t4 = -t9 * qJD(4) - t87 * t11 + t90 * t93;
t79 = -pkin(4) * t90 - pkin(3);
t71 = t149 * t90;
t70 = t149 * t87;
t65 = (pkin(4) * t87 + pkin(8)) * t88;
t64 = t90 * t69;
t57 = -t87 * qJD(5) + t90 * t113;
t56 = t87 * t113 + t131;
t47 = t64 - t130;
t43 = t94 * pkin(4) + t126;
t36 = -t87 * t143 + t48;
t31 = -t88 * t142 + t64 + (-pkin(4) - t155) * t91;
t30 = -t48 * qJD(4) + t146;
t29 = t95 * pkin(8) + t147;
t21 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t151 + (-qJD(5) * t88 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t91) * t87 - t147;
t19 = -t90 * t118 + t59 * t133 + t159 * t87;
t18 = -t88 * t131 + (pkin(4) * t88 - t91 * t142) * qJD(3) + (-t77 + (-t69 + t143) * t87) * qJD(4) + t146;
t15 = pkin(4) * t34 + t27;
t5 = pkin(4) * t19 + t12;
t2 = -qJ(5) * t19 - qJD(5) * t34 - t3;
t1 = t32 * pkin(4) - t20 * qJ(5) - t35 * qJD(5) + t4;
t10 = [0, 0, 0, 0.2e1 * t109, 0.2e1 * (-t89 ^ 2 + t92 ^ 2) * t80 * qJD(2), 0.2e1 * t86 * t117, -0.2e1 * t86 * t118, 0, -0.2e1 * pkin(1) * t80 * t141 - 0.2e1 * t55 * t86, -0.2e1 * pkin(1) * t124 + 0.2e1 * t54 * t86, 0.2e1 * t59 * t33, -0.2e1 * t32 * t59 - 0.2e1 * t33 * t58, (t59 * t141 - t33 * t92) * t158, (-t58 * t141 + t32 * t92) * t158, -0.2e1 * t109, 0.2e1 * t50 * t32 + 0.2e1 * t55 * t58 + 0.2e1 * (t114 * t141 - t14 * t92) * t85, 0.2e1 * t50 * t33 + 0.2e1 * t55 * t59 + 0.2e1 * (-t13 * t92 - t148 * t141) * t85, 0.2e1 * t35 * t20, -0.2e1 * t19 * t35 - 0.2e1 * t20 * t34, 0.2e1 * t20 * t58 + 0.2e1 * t32 * t35, -0.2e1 * t19 * t58 - 0.2e1 * t32 * t34, 0.2e1 * t58 * t32, 0.2e1 * t12 * t34 + 0.2e1 * t19 * t27 + 0.2e1 * t32 * t8 + 0.2e1 * t4 * t58, 0.2e1 * t12 * t35 + 0.2e1 * t20 * t27 + 0.2e1 * t3 * t58 - 0.2e1 * t32 * t9, -0.2e1 * t1 * t35 - 0.2e1 * t19 * t7 - 0.2e1 * t2 * t34 - 0.2e1 * t20 * t6, 0.2e1 * t1 * t6 + 0.2e1 * t15 * t5 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, t117, -t118, 0, -t55, t54, t59 * t136 + t33 * t88, -t88 * t32 + t33 * t91 + (-t58 * t91 - t59 * t88) * qJD(3), t97 * t85, t98 * t85, 0, -pkin(2) * t32 + t50 * t138 - t97 * t156 - t55 * t91, -pkin(2) * t33 + t50 * t136 - t98 * t156 + t55 * t88, t20 * t151 + t96 * t35, t104 * t136 + (-t19 * t90 - t154 + (t34 * t87 - t35 * t90) * qJD(4)) * t88, (t58 * t137 - t20) * t91 + (qJD(3) * t35 - t100) * t88, (-t87 * t139 + t19) * t91 + (-qJD(3) * t34 - t101) * t88, t58 * t138 - t32 * t91, t30 * t58 + t32 * t47 + (-t4 + (pkin(8) * t34 + t27 * t87) * qJD(3)) * t91 + (pkin(8) * t19 + qJD(3) * t8 + t103) * t88, t29 * t58 - t32 * t48 + (-t3 + (pkin(8) * t35 + t27 * t90) * qJD(3)) * t91 + (pkin(8) * t20 - qJD(3) * t9 - t102) * t88, -t18 * t35 - t19 * t36 - t20 * t31 - t21 * t34 + t105 * t136 + (-t1 * t90 - t2 * t87 + (t6 * t87 - t7 * t90) * qJD(4)) * t88, t1 * t31 + t15 * t43 + t18 * t6 + t2 * t36 + t21 * t7 + t5 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -0.2e1 * t111, 0, 0, 0, t88 * t128, t91 * t128, 0.2e1 * t83 * t115 - 0.2e1 * t82 * t116, t108 * t160 + 0.2e1 * t82 * t112, 0.2e1 * t88 * t120 + 0.2e1 * t144 * t137, -0.2e1 * t87 * t111 + 0.2e1 * t88 * t119, -0.2e1 * t115, 0.2e1 * t47 * t138 - 0.2e1 * t30 * t91 + 0.2e1 * (t110 * t87 + t133 * t82) * pkin(8), -0.2e1 * t48 * t138 - 0.2e1 * t29 * t91 + 0.2e1 * (t110 * t90 - t134 * t82) * pkin(8), 0.2e1 * (-t31 * t90 - t36 * t87) * t136 + 0.2e1 * (-t18 * t90 - t21 * t87 + (t31 * t87 - t36 * t90) * qJD(4)) * t88, 0.2e1 * t18 * t31 + 0.2e1 * t21 * t36 + 0.2e1 * t43 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t118, t14, t13, t35 * t133 + t154, t104 * qJD(4) - t87 * t19 + t20 * t90, t101, -t100, 0, -pkin(3) * t19 - pkin(9) * t101 + t102, -pkin(3) * t20 + pkin(9) * t100 + t103, qJD(4) * t105 - t1 * t87 + t19 * t71 + t2 * t90 - t20 * t70 - t34 * t56 - t35 * t57, t1 * t70 + t125 * t15 - t2 * t71 + t5 * t79 + t56 * t7 + t57 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t138, 0, -t126, pkin(8) * t138, -t88 * t112 + t108, t116 * t160 - t145 * t136, -t119 + t123, t95, 0, (pkin(9) * t150 + (-pkin(3) * t90 + t155) * t88) * qJD(4) + (t107 * t87 - t77) * qJD(3), (pkin(8) * t151 + t106 * t87) * qJD(4) + (t107 * t90 + t130) * qJD(3), (-t70 * t136 - t57 * t88 + t21 + (t71 * t88 - t31) * qJD(4)) * t90 + (t71 * t136 - t56 * t88 - t18 + (t70 * t88 - t36) * qJD(4)) * t87, t125 * t65 + t18 * t70 - t21 * t71 + t31 * t57 + t36 * t56 + t43 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t116, -0.2e1 * t112, 0, 0, 0, t87 * t127, t90 * t127, 0.2e1 * t56 * t90 - 0.2e1 * t57 * t87 + 0.2e1 * (-t70 * t90 + t71 * t87) * qJD(4), 0.2e1 * t125 * t79 - 0.2e1 * t56 * t71 + 0.2e1 * t57 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t32, t4, t3, -pkin(4) * t20, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t94, t138, t30, t29, -t96 * pkin(4), t18 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, -t134, 0, -pkin(9) * t133, pkin(9) * t134, -pkin(4) * t133, t57 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
