% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:43
% EndTime: 2019-03-08 20:16:47
% DurationCPUTime: 1.25s
% Computational Cost: add. (1640->255), mult. (3696->371), div. (0->0), fcn. (2466->8), ass. (0->147)
t78 = sin(pkin(6));
t140 = qJD(1) * t78;
t80 = sin(qJ(5));
t82 = sin(qJ(2));
t156 = t80 * t82;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t101 = pkin(4) * t84 + pkin(9) * t81;
t55 = t101 * qJD(4) + qJD(3);
t83 = cos(qJ(5));
t85 = cos(qJ(2));
t180 = t83 * t55 - (-t81 * t156 + t83 * t85) * t140;
t86 = -pkin(2) - pkin(8);
t153 = t81 * t86;
t64 = pkin(4) * t81 - pkin(9) * t84 + qJ(3);
t145 = t83 * t153 + t80 * t64;
t129 = t83 * qJD(4);
t131 = qJD(5) * t83;
t152 = t82 * t83;
t179 = (t81 * t152 + t80 * t85) * t140 - t84 * t86 * t129 - t64 * t131 - t80 * t55;
t139 = qJD(1) * t81;
t113 = t85 * t140;
t97 = qJD(3) - t113;
t48 = t86 * qJD(2) + t97;
t79 = cos(pkin(6));
t31 = -t79 * t139 + t48 * t84;
t130 = qJD(5) * t84;
t116 = t80 * t130;
t117 = t81 * t129;
t92 = t116 + t117;
t35 = t92 * qJD(2) - qJD(5) * t129;
t137 = qJD(2) * t81;
t119 = t80 * t137;
t136 = qJD(2) * t84;
t118 = t83 * t136;
t135 = qJD(4) * t80;
t60 = t118 + t135;
t36 = -qJD(4) * t119 + t60 * qJD(5);
t178 = t60 ^ 2;
t158 = t79 * t84;
t72 = qJD(1) * t158;
t32 = t81 * t48 + t72;
t24 = qJD(4) * pkin(9) + t32;
t114 = t82 * t140;
t39 = t64 * qJD(2) + t114;
t12 = -t24 * t80 + t83 * t39;
t6 = -qJ(6) * t60 + t12;
t73 = qJD(5) + t137;
t3 = pkin(5) * t73 + t6;
t177 = t3 - t6;
t155 = t80 * t86;
t112 = pkin(5) - t155;
t128 = t83 * qJD(6);
t132 = qJD(5) * t80;
t176 = qJ(6) * t117 - t145 * qJD(5) + (qJ(6) * t132 + t112 * qJD(4) - t128) * t84 + t180;
t115 = t83 * t130;
t175 = -qJ(6) * t115 + (-qJD(6) * t84 + (qJ(6) * qJD(4) - qJD(5) * t86) * t81) * t80 - t179;
t138 = qJD(2) * t78;
t121 = t82 * t138;
t104 = t84 * t121;
t134 = qJD(4) * t81;
t148 = -qJD(4) * t72 - t48 * t134;
t18 = -qJD(1) * t104 - t148;
t174 = t18 * t80;
t173 = t18 * t83;
t23 = -qJD(4) * pkin(4) - t31;
t172 = t23 * t80;
t171 = t23 * t83;
t170 = t35 * t80;
t169 = t35 * t81;
t168 = t36 * t81;
t58 = t80 * t136 - t129;
t166 = t58 * t73;
t165 = t60 * t73;
t127 = qJD(2) * qJ(3);
t63 = t114 + t127;
t164 = t63 * t85;
t163 = t73 * t80;
t162 = t73 * t81;
t161 = t73 * t83;
t160 = t78 * t85;
t88 = qJD(2) ^ 2;
t159 = t78 * t88;
t157 = t80 * t39;
t154 = t81 * t83;
t151 = t84 * t35;
t150 = -qJ(6) - pkin(9);
t62 = t101 * qJD(2);
t149 = t83 * t31 + t80 * t62;
t106 = qJD(5) * t150;
t147 = -qJ(6) * t119 + t80 * t106 + t128 - t149;
t107 = -t80 * t31 + t83 * t62;
t146 = -t80 * qJD(6) + t83 * t106 - (pkin(5) * t84 + qJ(6) * t154) * qJD(2) - t107;
t77 = t84 ^ 2;
t144 = t81 ^ 2 - t77;
t87 = qJD(4) ^ 2;
t143 = -t87 - t88;
t142 = qJ(6) * t84;
t141 = qJD(2) * pkin(2);
t133 = qJD(4) * t84;
t126 = qJD(2) * qJD(4);
t125 = t82 * t159;
t124 = t85 * t159;
t17 = t48 * t133 + (-qJD(4) * t79 + t121) * t139;
t33 = (t55 + t113) * qJD(2);
t123 = -t39 * t131 - t83 * t17 - t80 * t33;
t120 = t85 * t138;
t111 = pkin(5) * t80 - t86;
t110 = t84 * t126;
t109 = -t80 * t17 + t83 * t33;
t108 = t73 * t86 + t24;
t105 = t73 + t137;
t103 = t81 * t121;
t102 = qJD(5) * t81 + qJD(2);
t100 = -t63 + t114;
t13 = t24 * t83 + t157;
t7 = -qJ(6) * t58 + t13;
t99 = t3 * t83 + t7 * t80;
t98 = t3 * t80 - t7 * t83;
t96 = qJD(2) * t77 - t162;
t47 = -t81 * t160 + t158;
t28 = t78 * t152 - t47 * t80;
t29 = t78 * t156 + t47 * t83;
t46 = t84 * t160 + t79 * t81;
t95 = -pkin(9) * t133 + t23 * t81;
t94 = t100 * qJD(2);
t93 = -t24 * t132 - t123;
t91 = t100 - t127;
t10 = pkin(5) * t36 + t18;
t90 = -t13 * qJD(5) + t109;
t56 = (qJD(3) + t113) * qJD(2);
t89 = t97 * qJD(2) - t86 * t87 + t56;
t67 = t150 * t83;
t66 = t150 * t80;
t57 = t97 - t141;
t54 = t58 ^ 2;
t53 = t83 * t64;
t27 = t47 * qJD(4) - t104;
t26 = -t46 * qJD(4) + t103;
t25 = -t80 * t142 + t145;
t19 = t112 * t81 - t83 * t142 + t53;
t15 = pkin(5) * t58 + qJD(6) + t23;
t9 = t28 * qJD(5) + t80 * t120 + t26 * t83;
t8 = -t29 * qJD(5) + t83 * t120 - t26 * t80;
t2 = -qJ(6) * t36 - qJD(6) * t58 + t93;
t1 = pkin(5) * t110 + t35 * qJ(6) - t60 * qJD(6) + t90;
t4 = [0, 0, -t125, -t124, t125, t124 (t56 * t82 + (t164 + (t57 - t113) * t82) * qJD(2)) * t78, 0, 0, 0, 0, 0, t81 * t124 + (-t27 + t104) * qJD(4), t84 * t124 + (-t26 - t103) * qJD(4), 0, 0, 0, 0, 0, t110 * t28 + t27 * t58 + t36 * t46 + t73 * t8, -t110 * t29 + t27 * t60 - t35 * t46 - t73 * t9, t28 * t35 - t29 * t36 - t58 * t9 - t60 * t8, t1 * t28 + t10 * t46 + t15 * t27 + t2 * t29 + t3 * t8 + t7 * t9; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t56 * qJ(3) + t63 * qJD(3) + (-t164 + (-t57 - t141) * t82) * t140, -0.2e1 * t81 * t110, 0.2e1 * t144 * t126, -t87 * t81, -t87 * t84, 0, -t91 * t133 + t89 * t81, t91 * t134 + t89 * t84, -t83 * t151 - t60 * t92 (t58 * t83 + t60 * t80) * t134 + (t170 - t36 * t83 + (t58 * t80 - t60 * t83) * qJD(5)) * t84, -t73 * t116 - t169 + (t60 * t84 + t83 * t96) * qJD(4), -t73 * t115 - t168 + (-t58 * t84 - t80 * t96) * qJD(4), t105 * t133 (-t132 * t64 + t180) * t73 + ((t58 * t86 - t172) * qJD(4) + (-t108 * t83 - t157) * qJD(5) + t109) * t81 + (t58 * t114 + t23 * t131 + t174 - t86 * t36 + (-t73 * t155 + (-t153 * t80 + t53) * qJD(2) + t12) * qJD(4)) * t84, t179 * t73 + (t108 * t132 + (t60 * t86 - t171) * qJD(4) + t123) * t81 + (t60 * t114 - t23 * t132 + t173 + t86 * t35 + (-qJD(2) * t145 - t13) * qJD(4)) * t84, t19 * t35 - t25 * t36 - t176 * t60 - t175 * t58 + t99 * t134 + (qJD(5) * t98 - t1 * t83 - t2 * t80) * t84, t1 * t19 + t2 * t25 + t175 * t7 + t176 * t3 - t15 * t111 * t134 + (t10 * t111 + (pkin(5) * t131 + t114) * t15) * t84; 0, 0, 0, 0, 0, -t88, t94, 0, 0, 0, 0, 0, t143 * t81, t143 * t84, 0, 0, 0, 0, 0, -t84 * t36 - t102 * t161 + (-t105 * t80 * t84 + t58 * t81) * qJD(4), t151 + t102 * t163 + (-t84 * t161 + (t60 - t118) * t81) * qJD(4) (t102 * t60 - t133 * t58 - t168) * t83 + (t102 * t58 + t133 * t60 - t169) * t80, -t99 * qJD(2) + (-qJD(4) * t98 - t10) * t84 + (qJD(4) * t15 - qJD(5) * t99 - t1 * t80 + t2 * t83) * t81; 0, 0, 0, 0, 0, 0, 0, t84 * t88 * t81, -t144 * t88, 0, 0, 0, t32 * qJD(4) + t84 * t94 + t148, -t100 * t137, t60 * t161 - t170 (-t35 - t166) * t83 + (-t36 - t165) * t80, t73 * t131 + (t73 * t154 + (-t60 + t135) * t84) * qJD(2), -t73 * t132 + (-t80 * t162 + (t58 + t129) * t84) * qJD(2), -t73 * t136, -pkin(4) * t36 - t173 - t107 * t73 - t32 * t58 + (-pkin(9) * t161 + t172) * qJD(5) + (-t12 * t84 + t80 * t95) * qJD(2), pkin(4) * t35 + t174 + t149 * t73 - t32 * t60 + (pkin(9) * t163 + t171) * qJD(5) + (t13 * t84 + t83 * t95) * qJD(2), t66 * t35 + t67 * t36 - t146 * t60 - t147 * t58 + (-t3 * t73 + t2) * t83 + (-t7 * t73 - t1) * t80, -t2 * t67 + t1 * t66 + t10 * (-pkin(5) * t83 - pkin(4)) + t147 * t7 + t146 * t3 + (pkin(5) * t163 - t32) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t54 + t178, t166 - t35, t165 - t36, t110, t13 * t73 - t23 * t60 + t90, t12 * t73 + t23 * t58 - t93, pkin(5) * t35 - t177 * t58, t177 * t7 + (-t15 * t60 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 - t178, t3 * t60 + t58 * t7 + t10;];
tauc_reg  = t4;
