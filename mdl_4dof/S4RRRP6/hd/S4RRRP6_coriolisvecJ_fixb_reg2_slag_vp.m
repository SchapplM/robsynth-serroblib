% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:10
% EndTime: 2019-12-31 17:19:15
% DurationCPUTime: 1.36s
% Computational Cost: add. (1328->238), mult. (3537->344), div. (0->0), fcn. (2100->4), ass. (0->139)
t120 = qJD(1) * qJD(2);
t88 = sin(qJ(2));
t79 = t88 * t120;
t102 = pkin(5) * t79;
t89 = cos(qJ(3));
t125 = qJD(3) * t89;
t87 = sin(qJ(3));
t126 = qJD(3) * t87;
t90 = cos(qJ(2));
t66 = -t90 * pkin(2) - t88 * pkin(6) - pkin(1);
t51 = t66 * qJD(1);
t99 = pkin(2) * t88 - pkin(6) * t90;
t64 = t99 * qJD(2);
t52 = qJD(1) * t64;
t121 = t90 * qJD(1);
t83 = pkin(5) * t121;
t72 = qJD(2) * pkin(6) + t83;
t106 = t89 * t102 - t51 * t125 + t72 * t126 - t87 * t52;
t29 = t89 * t51 - t87 * t72;
t77 = -qJD(3) + t121;
t167 = t29 * t77 - t106;
t166 = -0.2e1 * t120;
t142 = t89 * t90;
t78 = pkin(5) * t142;
t40 = t87 * t66 + t78;
t30 = t87 * t51 + t89 * t72;
t14 = -qJD(3) * t30 + t87 * t102 + t89 * t52;
t165 = t30 * t77 - t14;
t128 = qJD(1) * t88;
t112 = t87 * t128;
t123 = t89 * qJD(2);
t59 = t112 - t123;
t152 = t59 * t77;
t110 = t90 * t120;
t119 = qJD(2) * qJD(3);
t33 = qJD(3) * t112 + (-t110 - t119) * t89;
t164 = -t33 + t152;
t111 = t89 * t128;
t124 = t87 * qJD(2);
t61 = t111 + t124;
t150 = t61 * t77;
t94 = t90 * t124 + t88 * t125;
t34 = t94 * qJD(1) + t87 * t119;
t163 = t34 - t150;
t162 = t61 ^ 2;
t161 = pkin(3) * t87;
t160 = t59 * pkin(3);
t21 = -t59 * qJ(4) + t30;
t159 = t21 * t77;
t26 = t34 * pkin(3) + pkin(5) * t110;
t158 = t26 * t87;
t157 = t26 * t89;
t154 = t33 * t87;
t153 = t34 * t89;
t151 = t61 * t59;
t149 = t61 * t87;
t129 = qJD(2) * pkin(2);
t82 = pkin(5) * t128;
t71 = t82 - t129;
t148 = t71 * t87;
t147 = t71 * t89;
t146 = t77 * t87;
t145 = t77 * t89;
t144 = t87 * t90;
t143 = t88 * t89;
t92 = qJD(1) ^ 2;
t141 = t90 * t92;
t91 = qJD(2) ^ 2;
t140 = t91 * t88;
t139 = t91 * t90;
t138 = -qJ(4) - pkin(6);
t20 = -t61 * qJ(4) + t29;
t15 = -t77 * pkin(3) + t20;
t137 = t15 - t20;
t107 = qJD(3) * t138;
t63 = t99 * qJD(1);
t35 = pkin(5) * t112 + t89 * t63;
t130 = qJ(4) * t90;
t96 = pkin(3) * t88 - t89 * t130;
t136 = t96 * qJD(1) + t87 * qJD(4) - t89 * t107 + t35;
t122 = t89 * qJD(4);
t47 = t87 * t63;
t135 = t47 + (-pkin(5) * t143 - t87 * t130) * qJD(1) - t87 * t107 - t122;
t134 = t66 * t125 + t87 * t64;
t133 = t88 * pkin(5) * t124 + t89 * t64;
t85 = t88 ^ 2;
t132 = -t90 ^ 2 + t85;
t131 = qJ(4) * t88;
t127 = qJD(2) * t90;
t118 = pkin(5) * t144;
t117 = t88 * t141;
t116 = pkin(5) * t127;
t115 = t88 * t126;
t114 = t77 * t125;
t113 = t77 * t128;
t108 = -qJD(4) - t160;
t105 = t59 + t123;
t104 = -t61 + t124;
t103 = pkin(1) * t166;
t101 = pkin(3) * t79;
t100 = t90 * t79;
t98 = -t29 * t89 - t30 * t87;
t97 = qJD(1) * t85 - t77 * t90;
t95 = t34 * qJ(4) + t106;
t93 = t33 * qJ(4) + t14;
t81 = -t89 * pkin(3) - pkin(2);
t70 = t138 * t89;
t69 = t138 * t87;
t65 = (pkin(5) + t161) * t88;
t58 = t89 * t66;
t56 = t59 ^ 2;
t54 = t121 * t161 + t83;
t39 = t58 - t118;
t38 = (-t77 - t121) * t88 * qJD(2);
t37 = t94 * pkin(3) + t116;
t36 = -pkin(5) * t111 + t47;
t32 = -t108 + t71;
t31 = -t87 * t131 + t40;
t28 = -t89 * t131 + t58 + (-pkin(5) * t87 - pkin(3)) * t90;
t24 = -t56 + t162;
t23 = -t150 - t34;
t22 = -t33 - t152;
t19 = -t40 * qJD(3) + t133;
t18 = (-t88 * t123 - t90 * t126) * pkin(5) + t134;
t17 = -t114 + (t104 * t88 + t77 * t142) * qJD(1);
t16 = t77 * t126 + (t105 * t88 - t77 * t144) * qJD(1);
t12 = -t59 * t146 - t153;
t11 = -t61 * t145 - t154;
t10 = t34 * t87 * t88 + t94 * t59;
t9 = -t33 * t143 + (t90 * t123 - t115) * t61;
t8 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t143 + (-qJD(4) * t88 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t90) * t87 + t134;
t7 = -t88 * t122 + t96 * qJD(2) + (-t78 + (-t66 + t131) * t87) * qJD(3) + t133;
t6 = t88 * t114 + t34 * t90 + (-t59 * t88 - t97 * t87) * qJD(2);
t5 = t77 * t115 + t33 * t90 + (t61 * t88 + t97 * t89) * qJD(2);
t4 = -t59 * qJD(4) - t95;
t3 = -t61 * qJD(4) + t101 + t93;
t2 = -t163 * t87 + t164 * t89;
t1 = (-t59 * t89 - t149) * t127 + (t154 - t153 + (t59 * t87 - t61 * t89) * qJD(3)) * t88;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100, t132 * t166, t139, -0.2e1 * t100, -t140, 0, -pkin(5) * t139 + t88 * t103, pkin(5) * t140 + t90 * t103, 0, 0, t9, t1, t5, t10, t6, t38, -t14 * t90 - t19 * t77 + (pkin(5) * t34 + t71 * t125) * t88 + ((pkin(5) * t59 + t148) * t90 + (t29 + (t39 + t118) * qJD(1)) * t88) * qJD(2), -t106 * t90 + t18 * t77 + (-pkin(5) * t33 - t71 * t126) * t88 + ((pkin(5) * t61 + t147) * t90 + (-t30 + (-t40 + t78) * qJD(1)) * t88) * qJD(2), -t18 * t59 - t19 * t61 + t39 * t33 - t40 * t34 + t98 * t127 + (t106 * t87 - t14 * t89 + (t29 * t87 - t30 * t89) * qJD(3)) * t88, -t106 * t40 + t14 * t39 + t30 * t18 + t29 * t19 + (t71 + t82) * t116, t9, t1, t5, t10, t6, t38, t65 * t34 + t37 * t59 - t7 * t77 + (t32 * t124 - t3) * t90 + (t32 * t125 + t158 + (qJD(1) * t28 + t15) * qJD(2)) * t88, -t65 * t33 + t37 * t61 + t8 * t77 + (t32 * t123 + t4) * t90 + (-t32 * t126 + t157 + (-qJD(1) * t31 - t21) * qJD(2)) * t88, t28 * t33 - t31 * t34 - t8 * t59 - t7 * t61 + (-t15 * t89 - t21 * t87) * t127 + (-t3 * t89 - t4 * t87 + (t15 * t87 - t21 * t89) * qJD(3)) * t88, t15 * t7 + t21 * t8 + t26 * t65 + t3 * t28 + t4 * t31 + t32 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t132 * t92, 0, t117, 0, 0, t92 * pkin(1) * t88, pkin(1) * t141, 0, 0, t11, t2, t17, t12, t16, t113, -pkin(2) * t34 + t35 * t77 + (pkin(6) * t145 + t148) * qJD(3) + ((-pkin(6) * t124 - t29) * t88 + (-t105 * pkin(5) - t148) * t90) * qJD(1), pkin(2) * t33 - t36 * t77 + (-pkin(6) * t146 + t147) * qJD(3) + ((-pkin(6) * t123 + t30) * t88 + (t104 * pkin(5) - t147) * t90) * qJD(1), t35 * t61 + t36 * t59 + ((qJD(3) * t61 - t34) * pkin(6) + t167) * t89 + ((qJD(3) * t59 - t33) * pkin(6) + t165) * t87, -t29 * t35 - t30 * t36 + (-t71 - t129) * t83 + (t98 * qJD(3) - t106 * t89 - t14 * t87) * pkin(6), t11, t2, t17, t12, t16, t113, -t157 + t81 * t34 - t54 * t59 + t136 * t77 + (t32 + t160) * t126 + (-t32 * t144 + (qJD(2) * t69 - t15) * t88) * qJD(1), t158 - t81 * t33 - t54 * t61 - t135 * t77 + (pkin(3) * t149 + t32 * t89) * qJD(3) + (-t32 * t142 + (qJD(2) * t70 + t21) * t88) * qJD(1), t69 * t33 + t70 * t34 + t136 * t61 + t135 * t59 + (t15 * t77 + t4) * t89 + (-t3 + t159) * t87, t26 * t81 + t3 * t69 - t4 * t70 + (pkin(3) * t126 - t54) * t32 - t135 * t21 - t136 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t24, t22, -t151, t23, t79, -t71 * t61 - t165, t71 * t59 - t167, 0, 0, t151, t24, t22, -t151, t23, t79, 0.2e1 * t101 - t159 + (t108 - t32) * t61 + t93, -t162 * pkin(3) - t20 * t77 + (qJD(4) + t32) * t59 + t95, t33 * pkin(3) - t137 * t59, t137 * t21 + (-t32 * t61 + t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t164, -t56 - t162, t15 * t61 + t21 * t59 + t26;];
tauc_reg = t13;
