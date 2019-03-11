% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:36
% EndTime: 2019-03-09 12:48:41
% DurationCPUTime: 1.55s
% Computational Cost: add. (2162->204), mult. (4678->362), div. (0->0), fcn. (3937->6), ass. (0->127)
t160 = pkin(3) + pkin(7);
t131 = qJD(4) + qJD(5);
t92 = sin(qJ(5));
t93 = sin(qJ(4));
t95 = cos(qJ(5));
t96 = cos(qJ(4));
t58 = t92 * t96 + t93 * t95;
t32 = t131 * t58;
t137 = qJD(5) * t92;
t141 = qJD(4) * t93;
t151 = t95 * t96;
t33 = t131 * t151 - t93 * t137 - t92 * t141;
t59 = -t92 * t93 + t151;
t159 = pkin(4) * t95;
t80 = pkin(5) + t159;
t165 = pkin(4) * (t33 * t92 + (t58 * t95 - t59 * t92) * qJD(5)) - t80 * t32;
t94 = sin(qJ(2));
t143 = t94 * qJ(3);
t97 = cos(qJ(2));
t98 = -pkin(2) - pkin(8);
t106 = -t97 * t98 + t143;
t51 = -pkin(1) - t106;
t119 = pkin(9) * t97 - t51;
t157 = t94 * pkin(4);
t70 = t160 * t94;
t61 = t96 * t70;
t27 = t119 * t93 + t157 + t61;
t60 = t93 * t70;
t148 = t96 * t51 + t60;
t150 = t96 * t97;
t30 = -pkin(9) * t150 + t148;
t28 = t95 * t30;
t149 = t92 * t27 + t28;
t156 = pkin(9) - t98;
t64 = t156 * t93;
t65 = t156 * t96;
t147 = -t95 * t64 - t92 * t65;
t139 = qJD(4) * t97;
t123 = t96 * t139;
t134 = t94 * qJD(2);
t125 = t93 * t134;
t164 = -t123 + t125;
t91 = t97 ^ 2;
t114 = qJD(2) * (t94 ^ 2 - t91);
t88 = t93 ^ 2;
t146 = -t96 ^ 2 + t88;
t113 = t146 * qJD(4);
t133 = t97 * qJD(3);
t71 = t160 * t97;
t135 = t71 * qJD(4);
t163 = t106 * qJD(2) - t133 - t135;
t161 = 0.2e1 * qJD(3);
t158 = pkin(5) * t32;
t154 = t59 * t32;
t152 = t93 * t97;
t144 = qJ(3) * t97;
t78 = t93 * pkin(4) + qJ(3);
t142 = qJD(2) * t71;
t140 = qJD(4) * t96;
t138 = qJD(4) * t98;
t136 = qJD(5) * t95;
t84 = t97 * qJD(2);
t74 = pkin(4) * t140 + qJD(3);
t132 = qJ(3) * qJD(4);
t130 = -0.2e1 * pkin(1) * qJD(2);
t129 = t95 * t150;
t45 = pkin(4) * t150 + t71;
t128 = pkin(4) * t137;
t127 = pkin(4) * t136;
t126 = pkin(7) * t134;
t124 = t93 * t139;
t122 = t94 * t84;
t121 = t96 * t134;
t120 = t93 * t140;
t112 = pkin(2) * t134 - t94 * qJD(3);
t36 = (pkin(8) * t94 - t144) * qJD(2) + t112;
t81 = pkin(7) * t84;
t63 = pkin(3) * t84 + t81;
t116 = -t93 * t36 + t96 * t63;
t11 = (-pkin(9) * t93 * t94 + pkin(4) * t97) * qJD(2) + (t119 * t96 - t60) * qJD(4) + t116;
t103 = t121 + t124;
t15 = -t70 * t140 + t51 * t141 - t96 * t36 - t93 * t63;
t13 = t103 * pkin(9) - t15;
t118 = t95 * t11 - t92 * t13;
t117 = t95 * t27 - t30 * t92;
t115 = t64 * t92 - t95 * t65;
t111 = t93 * t121;
t110 = -pkin(2) * t97 - t143;
t19 = -qJD(5) * t129 + (t131 * t152 + t121) * t92 + t164 * t95;
t39 = t58 * t97;
t109 = t19 * t59 + t32 * t39;
t108 = -t33 * t58 + t154;
t104 = -t33 * t94 - t58 * t84;
t3 = -t92 * t11 - t95 * t13 - t27 * t136 + t30 * t137;
t56 = t156 * t141;
t57 = qJD(4) * t65;
t17 = t65 * t136 - t64 * t137 - t92 * t56 + t95 * t57;
t62 = t160 * t134;
t102 = -t62 + (-t94 * t98 - t144) * qJD(4);
t4 = -t149 * qJD(5) + t118;
t1 = pkin(5) * t84 - t19 * qJ(6) + t39 * qJD(6) + t4;
t20 = t95 * t121 - t92 * t125 + t32 * t97;
t38 = t92 * t152 - t129;
t2 = qJ(6) * t20 + qJD(6) * t38 - t3;
t7 = pkin(5) * t94 + qJ(6) * t39 + t117;
t8 = qJ(6) * t38 + t149;
t101 = t1 * t59 + t2 * t58 - t32 * t7 + t33 * t8;
t25 = -qJ(6) * t59 + t115;
t26 = -qJ(6) * t58 + t147;
t5 = -qJ(6) * t33 - qJD(6) * t58 - t17;
t18 = -t147 * qJD(5) + t95 * t56 + t92 * t57;
t6 = t32 * qJ(6) - t59 * qJD(6) + t18;
t100 = t25 * t32 - t26 * t33 - t5 * t58 - t59 * t6;
t99 = t110 * qJD(2) + t133;
t34 = -pkin(4) * t124 + (-pkin(4) * t96 - t160) * t134;
t73 = 0.2e1 * t122;
t66 = -pkin(1) + t110;
t47 = -t94 * t140 - t93 * t84;
t46 = -t94 * t141 + t96 * t84;
t41 = -qJ(3) * t84 + t112;
t37 = pkin(5) * t58 + t78;
t31 = -pkin(5) * t38 + t45;
t29 = pkin(5) * t33 + t74;
t21 = -t32 * t94 + t59 * t84;
t16 = -qJD(4) * t148 + t116;
t14 = -t20 * pkin(5) + t34;
t9 = [0, 0, 0, t73, -0.2e1 * t114, 0, 0, 0, t94 * t130, t97 * t130, 0, -0.2e1 * t66 * t134 + 0.2e1 * t41 * t97, -0.2e1 * t41 * t94 - 0.2e1 * t66 * t84, 0.2e1 * t66 * t41, 0.2e1 * t91 * t120 - 0.2e1 * t88 * t122, -0.4e1 * t97 * t111 - 0.2e1 * t91 * t113, 0.2e1 * t93 * t114 - 0.2e1 * t94 * t123, 0.2e1 * t96 * t114 + 0.2e1 * t94 * t124, t73, 0.2e1 * (-t96 * t142 + t16) * t94 + 0.2e1 * ((-t51 * t93 + t61) * qJD(2) - t62 * t96 - t93 * t135) * t97, 0.2e1 * (t93 * t142 + t15) * t94 + 0.2e1 * (-t148 * qJD(2) - t96 * t135 + t62 * t93) * t97, -0.2e1 * t39 * t19, 0.2e1 * t19 * t38 - 0.2e1 * t20 * t39, 0.2e1 * t19 * t94 - 0.2e1 * t39 * t84, 0.2e1 * t20 * t94 + 0.2e1 * t38 * t84, t73, 0.2e1 * t117 * t84 - 0.2e1 * t45 * t20 - 0.2e1 * t34 * t38 + 0.2e1 * t4 * t94, -0.2e1 * t149 * t84 + 0.2e1 * t45 * t19 + 0.2e1 * t3 * t94 - 0.2e1 * t34 * t39, 0.2e1 * t1 * t39 - 0.2e1 * t19 * t7 + 0.2e1 * t2 * t38 + 0.2e1 * t20 * t8, 0.2e1 * t1 * t7 + 0.2e1 * t14 * t31 + 0.2e1 * t2 * t8; 0, 0, 0, 0, 0, t84, -t134, 0, -t81, t126, t99, t81, -t126, t99 * pkin(7), t97 * t113 + t111, 0.4e1 * t97 * t120 - t146 * t134, t46, t47, 0, t102 * t93 - t163 * t96, t102 * t96 + t163 * t93, t109, -t19 * t58 + t20 * t59 - t32 * t38 + t33 * t39, t21, t104, 0, t115 * t84 + t18 * t94 - t78 * t20 + t45 * t33 + t34 * t58 - t74 * t38, -t147 * t84 + t17 * t94 + t78 * t19 - t45 * t32 + t34 * t59 - t74 * t39, -t19 * t25 + t20 * t26 + t38 * t5 + t39 * t6 - t101, t1 * t25 + t14 * t37 + t2 * t26 + t29 * t31 + t5 * t8 + t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, qJ(3) * t161, -0.2e1 * t120, 0.2e1 * t113, 0, 0, 0, 0.2e1 * qJD(3) * t93 + 0.2e1 * t96 * t132, 0.2e1 * qJD(3) * t96 - 0.2e1 * t93 * t132, -0.2e1 * t154, 0.2e1 * t32 * t58 - 0.2e1 * t33 * t59, 0, 0, 0, 0.2e1 * t33 * t78 + 0.2e1 * t58 * t74, -0.2e1 * t32 * t78 + 0.2e1 * t59 * t74, 0.2e1 * t100, 0.2e1 * t25 * t6 + 0.2e1 * t26 * t5 + 0.2e1 * t29 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, t81, 0, 0, 0, 0, 0, t46, t47, 0, 0, 0, 0, 0, t21, t104, t20 * t58 + t33 * t38 - t109, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t108, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t103, t84, t16, t15, 0, 0, t19, t20, t84, t84 * t159 + (-t28 + (-t27 - t157) * t92) * qJD(5) + t118 (-t94 * t136 - t92 * t84) * pkin(4) + t3, -t80 * t19 + (t20 * t92 + (t38 * t95 - t39 * t92) * qJD(5)) * pkin(4), t1 * t80 + (t2 * t92 + (-t7 * t92 + t8 * t95) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t140, 0, -t93 * t138, -t96 * t138, 0, 0, -t32, -t33, 0, t18, t17, -t165, t6 * t80 + (t5 * t92 + (-t25 * t92 + t26 * t95) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t140, 0, 0, 0, 0, 0, -t32, -t33, 0, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t128, -0.2e1 * t127, 0, 0.2e1 * (-t80 + t159) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, t84, t4, t3, -pkin(5) * t19, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t18, t17, t158, t6 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, -t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t127, 0, -pkin(5) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
