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
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:03:00
% DurationCPUTime: 1.04s
% Computational Cost: add. (1688->171), mult. (2959->236), div. (0->0), fcn. (2230->8), ass. (0->125)
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t100 = cos(pkin(9));
t105 = cos(qJ(4));
t142 = t105 * t100;
t98 = qJD(1) + qJD(2);
t136 = t98 * t142;
t102 = sin(qJ(4));
t99 = sin(pkin(9));
t144 = t102 * t99;
t139 = t98 * t144;
t58 = -t136 + t139;
t50 = t104 * t58;
t79 = t102 * t100 + t105 * t99;
t60 = t79 * t98;
t29 = t101 * t60 + t50;
t97 = qJD(4) + qJD(5);
t153 = t29 * t97;
t141 = qJD(5) * t101;
t82 = qJD(4) * t136;
t54 = -qJD(4) * t139 + t82;
t71 = t79 * qJD(4);
t55 = t98 * t71;
t9 = -qJD(5) * t50 - t101 * t55 + t104 * t54 - t60 * t141;
t175 = t9 + t153;
t148 = t100 ^ 2 + t99 ^ 2;
t106 = cos(qJ(2));
t146 = pkin(1) * qJD(2);
t134 = qJD(1) * t146;
t77 = t98 * qJD(3) + t106 * t134;
t174 = t148 * t77;
t119 = t101 * t58 - t104 * t60;
t109 = t119 * qJD(5) - t101 * t54 - t104 * t55;
t154 = t119 * t97;
t173 = t109 - t154;
t147 = pkin(1) * qJD(1);
t135 = t106 * t147;
t121 = qJD(3) - t135;
t172 = t119 * t29;
t171 = t119 ^ 2 - t29 ^ 2;
t103 = sin(qJ(2));
t137 = t103 * t147;
t81 = t98 * qJ(3) + t137;
t132 = pkin(7) * t98 + t81;
t47 = t132 * t99;
t48 = t132 * t100;
t118 = t102 * t47 - t105 * t48;
t18 = -t58 * pkin(8) - t118;
t92 = -t100 * pkin(3) - pkin(2);
t57 = t92 * t98 + t121;
t33 = t58 * pkin(4) + t57;
t111 = t79 * t77;
t6 = -t54 * pkin(8) + t118 * qJD(4) - t111;
t170 = t33 * t29 + t18 * t141 + (-t18 * t97 - t6) * t101;
t167 = -t102 * t48 - t105 * t47;
t115 = t142 - t144;
t168 = t115 * t77;
t5 = -t55 * pkin(8) + t167 * qJD(4) + t168;
t169 = -t101 * t5 + t104 * t6 + t33 * t119;
t84 = (-pkin(7) - qJ(3)) * t99;
t93 = t100 * pkin(7);
t85 = t100 * qJ(3) + t93;
t116 = -t102 * t84 - t105 * t85;
t166 = t116 * qJD(4) - t121 * t79;
t165 = t100 * t103;
t164 = qJD(5) - t97;
t140 = t105 * qJD(4);
t163 = (qJD(3) * t99 + qJD(4) * t85) * t102 + t115 * t135 - qJD(3) * t142 - t84 * t140;
t162 = pkin(4) * t60;
t70 = t115 * qJD(4);
t161 = t70 * pkin(8);
t160 = t71 * pkin(4);
t159 = t115 * pkin(4);
t158 = t79 * pkin(8);
t39 = t101 * t115 + t104 * t79;
t20 = t39 * qJD(5) + t101 * t70 + t104 * t71;
t38 = t101 * t79 - t104 * t115;
t89 = t103 * t134;
t40 = t55 * pkin(4) + t89;
t157 = t33 * t20 + t40 * t38;
t19 = -t38 * qJD(5) - t101 * t71 + t104 * t70;
t156 = t33 * t19 + t40 * t39;
t155 = t106 * pkin(1);
t152 = t57 * t70 + t79 * t89;
t151 = -t115 * t89 + t57 * t71;
t143 = t104 * t18;
t138 = t103 * t146;
t17 = -t60 * pkin(8) + t167;
t14 = qJD(4) * pkin(4) + t17;
t133 = -pkin(4) * t97 - t14;
t87 = t106 * t146 + qJD(3);
t129 = t148 * t87;
t128 = t148 * t106;
t126 = t148 * qJD(3);
t67 = t71 * pkin(8);
t125 = -qJD(5) * (-t102 * t85 + t105 * t84 - t158) + t67 + t163;
t75 = t115 * pkin(8);
t124 = qJD(5) * (-t116 + t75) + t161 - t166;
t123 = (-qJD(2) + t98) * t147;
t122 = (-qJD(1) - t98) * t146;
t91 = t103 * pkin(1) + qJ(3);
t72 = (-pkin(7) - t91) * t99;
t73 = t100 * t91 + t93;
t117 = -t102 * t72 - t105 * t73;
t83 = t92 - t155;
t113 = -t137 + t160;
t110 = t72 * t140 + t87 * t142 + (-qJD(4) * t73 - t87 * t99) * t102;
t108 = t117 * qJD(4) - t79 * t87;
t80 = -t98 * pkin(2) + t121;
t66 = t71 * qJD(4);
t65 = t70 * qJD(4);
t56 = t92 - t159;
t49 = t138 + t160;
t45 = t83 - t159;
t27 = -t117 + t75;
t26 = -t102 * t73 + t105 * t72 - t158;
t21 = t54 * t79 + t60 * t70;
t16 = t20 * t97;
t15 = t19 * t97;
t12 = t108 - t161;
t11 = -t67 + t110;
t3 = t115 * t54 - t79 * t55 - t70 * t58 - t60 * t71;
t2 = -t119 * t19 + t9 * t39;
t1 = t109 * t39 + t119 * t20 - t19 * t29 - t9 * t38;
t4 = [0, 0, 0, 0, -t98 * t138 - t89, t106 * t122, t122 * t165, t98 * t129 + t174, t81 * t129 + t91 * t174 + (t80 + (-pkin(2) - t155) * qJD(1)) * t138, t21, t3, t65, -t66, 0, t108 * qJD(4) + t58 * t138 + t83 * t55 + t151, -t110 * qJD(4) + t60 * t138 + t83 * t54 + t152, t2, t1, t15, -t16, 0, t49 * t29 - t45 * t109 + (-t101 * t11 + t104 * t12 + (-t101 * t26 - t104 * t27) * qJD(5)) * t97 + t157, -t49 * t119 + t45 * t9 - (t101 * t12 + t104 * t11 + (-t101 * t27 + t104 * t26) * qJD(5)) * t97 + t156; 0, 0, 0, 0, t98 * t137 - t89, t106 * t123, t123 * t165, (-t128 * t147 + t126) * t98 + t174, t81 * t126 + qJ(3) * t174 + (-t81 * t128 + (-pkin(2) * qJD(2) - t80) * t103) * t147, t21, t3, t65, -t66, 0, qJD(4) * t166 - t58 * t137 + t92 * t55 + t151, qJD(4) * t163 - t60 * t137 + t92 * t54 + t152, t2, t1, t15, -t16, 0, -t56 * t109 + (t125 * t101 - t124 * t104) * t97 + t113 * t29 + t157, t56 * t9 + (t124 * t101 + t125 * t104) * t97 - t113 * t119 + t156; 0, 0, 0, 0, 0, 0, 0, -t148 * t98 ^ 2, -t148 * t98 * t81 + t89, 0, 0, 0, 0, 0, 0.2e1 * t60 * qJD(4), t82 + (-t58 - t139) * qJD(4), 0, 0, 0, 0, 0, -t109 - t154, t9 - t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, t82 + (t58 - t139) * qJD(4), 0, 0, -t57 * t60 - t111, t57 * t58 - t168, -t172, t171, t175, t173, 0, -t29 * t162 - (-t101 * t17 - t143) * t97 + (t133 * t101 - t143) * qJD(5) + t169, t119 * t162 + (t133 * qJD(5) + t17 * t97 - t5) * t104 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t171, t175, t173, 0, t164 * (-t101 * t14 - t143) + t169, (-t14 * t164 - t5) * t104 + t170;];
tauc_reg = t4;
