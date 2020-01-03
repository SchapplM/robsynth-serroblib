% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:58
% EndTime: 2020-01-03 12:04:01
% DurationCPUTime: 0.95s
% Computational Cost: add. (1692->175), mult. (2975->239), div. (0->0), fcn. (2238->8), ass. (0->128)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t101 = cos(pkin(9));
t106 = cos(qJ(4));
t144 = t106 * t101;
t99 = qJD(1) + qJD(2);
t138 = t99 * t144;
t100 = sin(pkin(9));
t103 = sin(qJ(4));
t145 = t103 * t100;
t140 = t99 * t145;
t58 = -t138 + t140;
t50 = t105 * t58;
t79 = t106 * t100 + t103 * t101;
t60 = t79 * t99;
t29 = t102 * t60 + t50;
t98 = qJD(4) + qJD(5);
t155 = t29 * t98;
t143 = qJD(5) * t102;
t82 = qJD(4) * t138;
t54 = -qJD(4) * t140 + t82;
t71 = t79 * qJD(4);
t55 = t99 * t71;
t9 = -qJD(5) * t50 - t102 * t55 + t105 * t54 - t60 * t143;
t177 = t9 + t155;
t150 = t100 ^ 2 + t101 ^ 2;
t107 = cos(qJ(2));
t148 = pkin(1) * qJD(2);
t136 = qJD(1) * t148;
t77 = t99 * qJD(3) + t107 * t136;
t176 = t150 * t77;
t119 = t102 * t58 - t105 * t60;
t110 = t119 * qJD(5) - t102 * t54 - t105 * t55;
t156 = t119 * t98;
t175 = t110 - t156;
t149 = pkin(1) * qJD(1);
t137 = t107 * t149;
t121 = qJD(3) - t137;
t174 = t119 * t29;
t173 = t119 ^ 2 - t29 ^ 2;
t104 = sin(qJ(2));
t139 = t104 * t149;
t81 = t99 * qJ(3) + t139;
t134 = pkin(7) * t99 + t81;
t47 = t134 * t100;
t48 = t134 * t101;
t118 = t103 * t47 - t106 * t48;
t18 = -t58 * pkin(8) - t118;
t93 = -t101 * pkin(3) - pkin(2);
t57 = t93 * t99 + t121;
t33 = t58 * pkin(4) + t57;
t112 = t79 * t77;
t6 = -t54 * pkin(8) + t118 * qJD(4) - t112;
t172 = t33 * t29 + t18 * t143 + (-t18 * t98 - t6) * t102;
t169 = -t103 * t48 - t106 * t47;
t78 = -t144 + t145;
t170 = t78 * t77;
t5 = -t55 * pkin(8) + t169 * qJD(4) - t170;
t171 = -t102 * t5 + t105 * t6 + t33 * t119;
t84 = (-pkin(7) - qJ(3)) * t100;
t94 = t101 * pkin(7);
t85 = t101 * qJ(3) + t94;
t116 = -t103 * t84 - t106 * t85;
t168 = t116 * qJD(4) - t121 * t79;
t167 = t101 * t104;
t166 = qJD(5) - t98;
t142 = t106 * qJD(4);
t165 = (qJD(3) * t100 + qJD(4) * t85) * t103 - t78 * t137 - qJD(3) * t144 - t84 * t142;
t164 = pkin(4) * t60;
t70 = t78 * qJD(4);
t163 = t70 * pkin(8);
t162 = t71 * pkin(4);
t161 = t78 * pkin(4);
t160 = t79 * pkin(8);
t39 = -t102 * t78 + t105 * t79;
t20 = t39 * qJD(5) - t102 * t70 + t105 * t71;
t38 = t102 * t79 + t105 * t78;
t90 = t104 * t136;
t40 = t55 * pkin(4) + t90;
t159 = t33 * t20 + t40 * t38;
t19 = -t38 * qJD(5) - t102 * t71 - t105 * t70;
t158 = t33 * t19 + t40 * t39;
t157 = t107 * pkin(1);
t154 = -t57 * t70 + t79 * t90;
t153 = t57 * t71 + t78 * t90;
t146 = t105 * t18;
t141 = t104 * t148;
t17 = -t60 * pkin(8) + t169;
t14 = qJD(4) * pkin(4) + t17;
t135 = -pkin(4) * t98 - t14;
t88 = t107 * t148 + qJD(3);
t131 = t150 * t88;
t130 = t150 * t107;
t128 = t150 * qJD(3);
t127 = t99 * t141;
t126 = t99 * t139;
t67 = t71 * pkin(8);
t125 = -qJD(5) * (-t103 * t85 + t106 * t84 - t160) + t67 + t165;
t75 = t78 * pkin(8);
t124 = qJD(5) * (-t116 - t75) - t163 - t168;
t123 = (-qJD(2) + t99) * t149;
t122 = (-qJD(1) - t99) * t148;
t92 = t104 * pkin(1) + qJ(3);
t72 = (-pkin(7) - t92) * t100;
t73 = t101 * t92 + t94;
t117 = -t103 * t72 - t106 * t73;
t83 = t93 - t157;
t114 = -t139 + t162;
t111 = t72 * t142 + t88 * t144 + (-qJD(4) * t73 - t100 * t88) * t103;
t109 = t117 * qJD(4) - t79 * t88;
t86 = t100 * t90;
t80 = -t99 * pkin(2) + t121;
t66 = t71 * qJD(4);
t65 = t70 * qJD(4);
t56 = t93 + t161;
t49 = t141 + t162;
t45 = t83 + t161;
t27 = -t117 - t75;
t26 = -t103 * t73 + t106 * t72 - t160;
t21 = t54 * t79 - t60 * t70;
t16 = t20 * t98;
t15 = t19 * t98;
t12 = t109 + t163;
t11 = -t67 + t111;
t3 = -t54 * t78 - t79 * t55 + t70 * t58 - t60 * t71;
t2 = -t119 * t19 + t9 * t39;
t1 = t110 * t39 + t119 * t20 - t19 * t29 - t9 * t38;
t4 = [0, 0, 0, 0, -t90 - t127, t107 * t122, t122 * t167, t100 * t127 + t86, t99 * t131 + t176, t81 * t131 + t92 * t176 + (t80 + (-pkin(2) - t157) * qJD(1)) * t141, t21, t3, -t65, -t66, 0, t109 * qJD(4) + t58 * t141 + t83 * t55 + t153, -t111 * qJD(4) + t60 * t141 + t83 * t54 + t154, t2, t1, t15, -t16, 0, t49 * t29 - t45 * t110 + (-t102 * t11 + t105 * t12 + (-t102 * t26 - t105 * t27) * qJD(5)) * t98 + t159, -t49 * t119 + t45 * t9 - (t102 * t12 + t105 * t11 + (-t102 * t27 + t105 * t26) * qJD(5)) * t98 + t158; 0, 0, 0, 0, -t90 + t126, t107 * t123, t123 * t167, -t100 * t126 + t86, (-t130 * t149 + t128) * t99 + t176, t81 * t128 + qJ(3) * t176 + (-t81 * t130 + (-pkin(2) * qJD(2) - t80) * t104) * t149, t21, t3, -t65, -t66, 0, t168 * qJD(4) - t58 * t139 + t93 * t55 + t153, t165 * qJD(4) - t60 * t139 + t93 * t54 + t154, t2, t1, t15, -t16, 0, -t56 * t110 + (t125 * t102 - t124 * t105) * t98 + t114 * t29 + t159, t56 * t9 + (t124 * t102 + t125 * t105) * t98 - t114 * t119 + t158; 0, 0, 0, 0, 0, 0, 0, 0, -t150 * t99 ^ 2, -t150 * t99 * t81 + t90, 0, 0, 0, 0, 0, 0.2e1 * t60 * qJD(4), t82 + (-t58 - t140) * qJD(4), 0, 0, 0, 0, 0, -t110 - t156, t9 - t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, t82 + (t58 - t140) * qJD(4), 0, 0, -t57 * t60 - t112, t57 * t58 + t170, -t174, t173, t177, t175, 0, -t29 * t164 - (-t102 * t17 - t146) * t98 + (t135 * t102 - t146) * qJD(5) + t171, t119 * t164 + (t135 * qJD(5) + t17 * t98 - t5) * t105 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t173, t177, t175, 0, t166 * (-t102 * t14 - t146) + t171, (-t166 * t14 - t5) * t105 + t172;];
tauc_reg = t4;
