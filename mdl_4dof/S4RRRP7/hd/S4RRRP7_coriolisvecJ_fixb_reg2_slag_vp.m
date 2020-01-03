% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:09
% DurationCPUTime: 1.59s
% Computational Cost: add. (1319->247), mult. (3477->341), div. (0->0), fcn. (2028->4), ass. (0->134)
t71 = sin(qJ(2));
t132 = qJD(1) * t71;
t72 = cos(qJ(3));
t109 = t72 * t132;
t70 = sin(qJ(3));
t126 = t70 * qJD(2);
t50 = t109 + t126;
t73 = cos(qJ(2));
t124 = t73 * qJD(1);
t60 = -qJD(3) + t124;
t145 = t50 * t60;
t110 = t70 * t132;
t125 = t72 * qJD(2);
t48 = t110 - t125;
t148 = t48 * t60;
t122 = qJD(1) * qJD(2);
t106 = t73 * t122;
t121 = qJD(2) * qJD(3);
t27 = qJD(3) * t110 + (-t106 - t121) * t72;
t128 = qJD(3) * t72;
t79 = t73 * t126 + t71 * t128;
t28 = t79 * qJD(1) + t70 * t121;
t172 = (t27 - t148) * t72 + (t28 - t145) * t70;
t63 = t71 * t122;
t100 = pkin(5) * t63;
t129 = qJD(3) * t70;
t55 = -t73 * pkin(2) - t71 * pkin(6) - pkin(1);
t43 = t55 * qJD(1);
t97 = pkin(2) * t71 - pkin(6) * t73;
t53 = t97 * qJD(2);
t44 = qJD(1) * t53;
t66 = pkin(5) * t124;
t58 = qJD(2) * pkin(6) + t66;
t104 = -t70 * t100 + t58 * t128 + t43 * t129 - t72 * t44;
t20 = t70 * t43 + t72 * t58;
t81 = -t20 * t60 - t104;
t130 = qJD(2) * t73;
t149 = t28 * t72;
t150 = t27 * t70;
t169 = ((t48 * t70 - t50 * t72) * qJD(3) - t149 + t150) * t71 - (t48 * t72 + t50 * t70) * t130;
t112 = t60 * t128;
t68 = t71 ^ 2;
t88 = qJD(1) * t68 - t60 * t73;
t168 = (t48 * t71 + t88 * t70) * qJD(2) - t71 * t112 - t73 * t28;
t103 = t48 + t125;
t142 = t70 * t73;
t167 = (t103 * t71 - t60 * t142) * qJD(1) + t60 * t129;
t166 = -0.2e1 * t122;
t17 = -t60 * qJ(4) + t20;
t99 = pkin(3) * t63;
t5 = -t99 + t104;
t165 = t17 * t60 + t5;
t163 = t28 + t145;
t127 = qJD(3) * t73;
t14 = (t71 * t126 - t72 * t127) * pkin(5) - t55 * t129 + t72 * t53;
t159 = t50 ^ 2;
t158 = pkin(5) * t72;
t157 = pkin(6) * t50;
t156 = pkin(6) * t60;
t3 = t28 * pkin(3) + pkin(5) * t106 + t27 * qJ(4) - t50 * qJD(4);
t155 = t3 * t70;
t154 = t3 * t72;
t133 = qJD(2) * pkin(2);
t65 = pkin(5) * t132;
t57 = t65 - t133;
t18 = t48 * pkin(3) - t50 * qJ(4) + t57;
t152 = t18 * t50;
t146 = t50 * t48;
t144 = t57 * t70;
t143 = t57 * t72;
t52 = t97 * qJD(1);
t141 = t72 * t52;
t140 = t72 * t55;
t139 = t72 * t73;
t75 = qJD(1) ^ 2;
t138 = t73 * t75;
t74 = qJD(2) ^ 2;
t137 = t74 * t71;
t136 = t74 * t73;
t95 = pkin(3) * t70 - qJ(4) * t72;
t135 = t70 * qJD(4) + t60 * t95 + t66;
t62 = pkin(5) * t139;
t34 = t70 * t55 + t62;
t134 = -t73 ^ 2 + t68;
t131 = qJD(2) * t71;
t19 = t72 * t43 - t70 * t58;
t123 = qJD(4) - t19;
t120 = t70 * t156;
t119 = t72 * t156;
t118 = pkin(5) * t142;
t116 = t71 * t138;
t115 = pkin(6) * t125;
t113 = t71 * t129;
t111 = t60 * t132;
t108 = t48 ^ 2 - t159;
t107 = pkin(5) * t70 + pkin(3);
t102 = -t50 + t126;
t101 = pkin(1) * t166;
t98 = t73 * t63;
t96 = t72 * pkin(3) + t70 * qJ(4);
t16 = t60 * pkin(3) + t123;
t94 = t16 * t72 - t17 * t70;
t93 = -t19 * t72 - t20 * t70;
t89 = (qJD(3) * t48 - t27) * pkin(6);
t85 = t43 * t128 - t58 * t129 + t70 * t44;
t84 = pkin(5) + t95;
t83 = (qJ(4) - t158) * t132;
t77 = -t148 * t70 - t149;
t13 = t70 * t53 + t55 * t128 + (-t71 * t125 - t70 * t127) * pkin(5);
t7 = -t72 * t100 + t85;
t76 = t28 * t70 * t71 + t79 * t48;
t54 = -pkin(2) - t96;
t41 = t70 * t52;
t37 = t84 * t71;
t33 = -t118 + t140;
t32 = (-t60 - t124) * t131;
t30 = -pkin(5) * t109 + t41;
t29 = pkin(5) * t110 + t141;
t26 = t107 * t73 - t140;
t25 = -t73 * qJ(4) + t34;
t24 = pkin(6) * t149;
t23 = t50 * pkin(3) + t48 * qJ(4);
t22 = -t107 * t132 - t141;
t21 = t41 + t83;
t15 = -t27 - t148;
t12 = (t96 * qJD(3) - qJD(4) * t72) * t71 + t84 * t130;
t11 = -t112 + (t102 * t71 + t60 * t139) * qJD(1);
t10 = -pkin(3) * t131 - t14;
t9 = qJ(4) * t131 - t73 * qJD(4) + t13;
t6 = -t145 * t72 - t150;
t4 = -t27 * t72 * t71 + (t73 * t125 - t113) * t50;
t2 = qJD(2) * t83 - t60 * qJD(4) + t85;
t1 = t60 * t113 + t27 * t73 + (t50 * t71 + t88 * t72) * qJD(2);
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t98, t134 * t166, t136, -0.2e1 * t98, -t137, 0, -pkin(5) * t136 + t71 * t101, pkin(5) * t137 + t73 * t101, 0, 0, t4, t169, t1, t76, -t168, t32, -t14 * t60 + t104 * t73 + (pkin(5) * t28 + t57 * t128) * t71 + ((pkin(5) * t48 + t144) * t73 + (t19 + (t33 + t118) * qJD(1)) * t71) * qJD(2), t13 * t60 + t7 * t73 + (-pkin(5) * t27 - t57 * t129) * t71 + ((pkin(5) * t50 + t143) * t73 + (-t20 + (-t34 + t62) * qJD(1)) * t71) * qJD(2), -t13 * t48 - t14 * t50 + t33 * t27 - t34 * t28 + t93 * t130 + (-t7 * t70 + t72 * t104 + (t19 * t70 - t20 * t72) * qJD(3)) * t71, t20 * t13 + t19 * t14 - t104 * t33 + t7 * t34 + (t57 + t65) * pkin(5) * t130, t4, t1, -t169, t32, t168, t76, t10 * t60 + t12 * t48 + t37 * t28 + (t18 * t126 + t5) * t73 + (t18 * t128 + t155 + (-qJD(1) * t26 - t16) * qJD(2)) * t71, t10 * t50 - t25 * t28 - t26 * t27 - t9 * t48 + t94 * t130 + (-t2 * t70 + t5 * t72 + (-t16 * t70 - t17 * t72) * qJD(3)) * t71, -t12 * t50 + t37 * t27 - t9 * t60 + (-t18 * t125 - t2) * t73 + (t18 * t129 - t154 + (qJD(1) * t25 + t17) * qJD(2)) * t71, t16 * t10 + t18 * t12 + t17 * t9 + t2 * t25 + t5 * t26 + t3 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t134 * t75, 0, t116, 0, 0, t75 * pkin(1) * t71, pkin(1) * t138, 0, 0, t6, -t172, t11, t77, t167, t111, -pkin(2) * t28 + t29 * t60 + (t119 + t144) * qJD(3) + ((-pkin(6) * t126 - t19) * t71 + (-t103 * pkin(5) - t144) * t73) * qJD(1), pkin(2) * t27 - t30 * t60 + (-t120 + t143) * qJD(3) + ((t20 - t115) * t71 + (t102 * pkin(5) - t143) * t73) * qJD(1), t29 * t50 + t30 * t48 - t24 + (t19 * t124 + t7 + (-t19 + t157) * qJD(3)) * t72 + (t89 - t81) * t70, -t19 * t29 - t20 * t30 + (-t57 - t133) * t66 + (t93 * qJD(3) + t104 * t70 + t7 * t72) * pkin(6), t6, t11, t172, t111, -t167, t77, -t22 * t60 + t54 * t28 - t154 - t135 * t48 + (t18 * t70 + t119) * qJD(3) + (t16 * t71 + (-pkin(6) * t131 - t18 * t73) * t70) * qJD(1), t21 * t48 - t22 * t50 - t24 + (-t16 * t124 + t2 + (t16 + t157) * qJD(3)) * t72 + (t89 + t165) * t70, t21 * t60 + t54 * t27 - t155 + t135 * t50 + (-t18 * t72 + t120) * qJD(3) + (t18 * t139 + (-t17 + t115) * t71) * qJD(1), -t16 * t22 - t17 * t21 + t3 * t54 - t135 * t18 + (qJD(3) * t94 + t2 * t72 + t5 * t70) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, -t108, t15, -t146, -t163, t63, -t57 * t50 + t81, -t19 * t60 + t57 * t48 - t7, 0, 0, t146, t15, t108, t63, t163, -t146, -t23 * t48 - t152 + t81 + 0.2e1 * t99, pkin(3) * t27 - t28 * qJ(4) + (t17 - t20) * t50 + (t16 - t123) * t48, -t18 * t48 + t23 * t50 + (-0.2e1 * qJD(4) + t19) * t60 + (0.2e1 * qJ(4) - t158) * t63 + t85, -t5 * pkin(3) + t2 * qJ(4) + t123 * t17 - t16 * t20 - t18 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63 + t146, t15, -t60 ^ 2 - t159, t152 + t165;];
tauc_reg = t8;
