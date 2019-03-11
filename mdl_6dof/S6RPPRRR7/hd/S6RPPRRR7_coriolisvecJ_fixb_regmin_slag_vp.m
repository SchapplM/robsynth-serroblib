% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:48
% EndTime: 2019-03-09 02:33:55
% DurationCPUTime: 2.34s
% Computational Cost: add. (3754->262), mult. (8693->354), div. (0->0), fcn. (6816->8), ass. (0->157)
t123 = cos(qJ(6));
t166 = qJD(6) * t123;
t121 = sin(qJ(5));
t124 = cos(qJ(5));
t117 = sin(pkin(10));
t118 = cos(pkin(10));
t125 = cos(qJ(4));
t170 = qJD(1) * t125;
t122 = sin(qJ(4));
t171 = qJD(1) * t122;
t87 = -t117 * t170 - t118 * t171;
t158 = t118 * t170;
t159 = t117 * t171;
t88 = t158 - t159;
t59 = t121 * t88 - t124 * t87;
t220 = t123 * t59;
t231 = t166 + t220;
t92 = -t122 * t117 + t125 * t118;
t114 = qJD(4) + qJD(5);
t91 = t125 * t117 + t122 * t118;
t139 = t121 * t92 + t124 * t91;
t89 = t91 * qJD(4);
t164 = t125 * qJD(4);
t165 = t122 * qJD(4);
t90 = -t117 * t165 + t118 * t164;
t41 = t139 * qJD(5) + t121 * t90 + t124 * t89;
t230 = t41 * t114;
t180 = t59 * t114;
t168 = qJD(5) * t124;
t169 = qJD(5) * t121;
t80 = qJD(1) * t89;
t97 = qJD(4) * t159;
t81 = qJD(4) * t158 - t97;
t33 = -t121 * t81 - t124 * t80 + t87 * t168 - t88 * t169;
t229 = t33 + t180;
t216 = -qJD(6) - t59;
t228 = qJD(6) + t216;
t140 = t121 * t87 + t124 * t88;
t120 = sin(qJ(6));
t167 = qJD(6) * t120;
t16 = t114 * t166 + t123 * t33 - t140 * t167;
t53 = t120 * t114 + t123 * t140;
t63 = t121 * t91 - t124 * t92;
t227 = -t63 * t16 - t41 * t53;
t17 = qJD(6) * t53 + t120 * t33;
t51 = -t123 * t114 + t120 * t140;
t226 = -t120 * t17 + t16 * t123 - t231 * t51;
t12 = t16 * t120;
t225 = t231 * t53 + t12;
t34 = t140 * qJD(5) - t121 * t80 + t124 * t81;
t29 = t120 * t34;
t56 = t216 * t166;
t190 = t29 - t56;
t195 = t53 * t140;
t224 = -t216 * t220 + t190 - t195;
t119 = -pkin(1) - qJ(3);
t205 = t119 * qJD(1);
t98 = qJD(2) + t205;
t155 = -pkin(7) * qJD(1) + t98;
t82 = t155 * t117;
t83 = t155 * t118;
t138 = -t122 * t83 - t125 * t82;
t48 = t87 * pkin(8) - t138;
t186 = t121 * t48;
t207 = -t122 * t82 + t125 * t83;
t47 = -t88 * pkin(8) + t207;
t46 = qJD(4) * pkin(4) + t47;
t20 = t124 * t46 - t186;
t18 = -t114 * pkin(5) - t20;
t223 = t18 * t59;
t116 = qJD(1) * qJ(2);
t110 = qJD(3) + t116;
t111 = t117 * pkin(3);
t95 = qJD(1) * t111 + t110;
t67 = -t87 * pkin(4) + t95;
t221 = t67 * t59;
t219 = t140 * t59;
t179 = qJD(1) * t91;
t218 = t87 - t179;
t217 = t120 * t216;
t215 = t92 * qJD(3);
t181 = t140 * t114;
t214 = -t34 + t181;
t212 = t140 ^ 2 - t59 ^ 2;
t38 = pkin(5) * t140 + t59 * pkin(9);
t194 = t140 * t51;
t209 = t216 * t140;
t31 = t123 * t34;
t208 = -t167 * t216 - t31;
t127 = -t63 * qJD(5) - t121 * t89 + t124 * t90;
t206 = t127 * t114;
t172 = t117 ^ 2 + t118 ^ 2;
t204 = t172 * qJD(3);
t182 = t124 * t48;
t21 = t121 * t46 + t182;
t19 = t114 * pkin(9) + t21;
t24 = t59 * pkin(5) - pkin(9) * t140 + t67;
t141 = t120 * t19 - t123 * t24;
t203 = t140 * t141 + t18 * t167;
t132 = t91 * qJD(3);
t35 = -t81 * pkin(8) - qJD(1) * t132 + t207 * qJD(4);
t135 = t215 * qJD(1);
t36 = t80 * pkin(8) + t138 * qJD(4) - t135;
t152 = t121 * t35 - t124 * t36;
t3 = t21 * qJD(5) + t152;
t5 = t120 * t24 + t123 * t19;
t202 = t3 * t120 + t5 * t140 + t18 * t166;
t201 = -t67 * t140 - t152;
t151 = t121 * t36 - t48 * t169;
t2 = (qJD(5) * t46 + t35) * t124 + t151;
t192 = -pkin(7) + t119;
t93 = t192 * t117;
t94 = t192 * t118;
t54 = -t92 * pkin(8) - t122 * t93 + t125 * t94;
t137 = -t122 * t94 - t125 * t93;
t55 = -t91 * pkin(8) - t137;
t27 = t121 * t54 + t124 * t55;
t106 = qJ(2) + t111;
t73 = t91 * pkin(4) + t106;
t28 = pkin(5) * t139 + pkin(9) * t63 + t73;
t26 = t121 * t55 - t124 * t54;
t129 = t94 * t164 - t93 * t165 - t132;
t44 = -t90 * pkin(8) + t129;
t128 = t137 * qJD(4) - t215;
t45 = t89 * pkin(8) + t128;
t6 = -t26 * qJD(5) + t121 * t45 + t124 * t44;
t200 = -t139 * (qJD(6) * t24 + t2) + (qJD(6) * t28 + t6) * t216 - t18 * t41 - t27 * t34 - t3 * t63;
t115 = qJD(1) * qJD(2);
t162 = 0.2e1 * t115;
t199 = t88 * pkin(4);
t197 = t18 * t63;
t196 = t28 * t34;
t193 = t63 * t34;
t188 = t120 * t53;
t175 = t90 * qJD(4);
t160 = t63 * t167;
t156 = -pkin(4) * t114 - t46;
t68 = t81 * pkin(4) + t115;
t75 = t90 * pkin(4) + qJD(2);
t108 = t121 * pkin(4) + pkin(9);
t146 = qJD(6) * t108 + t199 + t38;
t145 = qJD(1) * t172;
t144 = qJD(6) * t139 + qJD(1);
t22 = t121 * t47 + t182;
t143 = pkin(4) * t169 - t22;
t142 = t216 * t41 - t193;
t136 = t217 * t59 - t208;
t23 = t124 * t47 - t186;
t130 = -t108 * t34 + t223 - (-pkin(4) * t168 + t23) * t216;
t126 = qJD(1) ^ 2;
t109 = -t124 * pkin(4) - pkin(5);
t84 = t89 * qJD(4);
t10 = pkin(5) * t127 + pkin(9) * t41 + t75;
t9 = t34 * pkin(5) - t33 * pkin(9) + t68;
t8 = t123 * t9;
t7 = t27 * qJD(5) + t121 * t44 - t124 * t45;
t1 = [0, 0, 0, 0, t162, qJ(2) * t162, t117 * t162, t118 * t162, 0.2e1 * qJD(3) * t145 (t110 + t116) * qJD(2) + (-t98 - t205) * t204, -t80 * t92 - t88 * t89, t80 * t91 - t92 * t81 - t89 * t87 - t88 * t90, -t84, -t175, 0, -t218 * qJD(2) + t128 * qJD(4) + t106 * t81 + t95 * t90, -t106 * t80 - t95 * t89 - t129 * qJD(4) + (qJD(1) * t92 + t88) * qJD(2), -t140 * t41 - t33 * t63, -t127 * t140 - t139 * t33 + t41 * t59 + t193, -t230, -t206, 0, -t7 * t114 + t127 * t67 + t139 * t68 + t73 * t34 + t75 * t59, -t6 * t114 + t140 * t75 + t73 * t33 - t41 * t67 - t63 * t68, t123 * t227 + t53 * t160 -(-t123 * t51 - t188) * t41 - (-t12 - t123 * t17 + (t120 * t51 - t123 * t53) * qJD(6)) * t63, t123 * t142 + t127 * t53 + t139 * t16 - t160 * t216, -t120 * t142 - t127 * t51 - t139 * t17 - t56 * t63, -t127 * t216 + t139 * t34, t26 * t17 - t141 * t127 + t7 * t51 + t8 * t139 + (-t10 * t216 + t196 + (-t139 * t19 + t216 * t27 - t197) * qJD(6)) * t123 + t200 * t120, t26 * t16 - t5 * t127 + t7 * t53 + ((-qJD(6) * t27 + t10) * t216 - t196 - (-qJD(6) * t19 + t9) * t139 + qJD(6) * t197) * t120 + t200 * t123; 0, 0, 0, 0, -t126, -t126 * qJ(2), -t126 * t117, -t126 * t118, 0 (-t110 - t204) * qJD(1), 0, 0, 0, 0, 0, qJD(1) * t87 - t84, -qJD(1) * t88 - t175, 0, 0, 0, 0, 0, -qJD(1) * t59 - t230, -qJD(1) * t140 - t206, 0, 0, 0, 0, 0, -t139 * t29 + t63 * t17 + t41 * t51 - (-t120 * t127 - t123 * t144) * t216, -t139 * t31 - (t120 * t144 - t123 * t127) * t216 - t227; 0, 0, 0, 0, 0, 0, 0, 0, -t172 * t126, t98 * t145 + t115, 0, 0, 0, 0, 0, -t97 + (t88 + t158) * qJD(4), t218 * qJD(4), 0, 0, 0, 0, 0, t34 + t181, t33 - t180, 0, 0, 0, 0, 0, t136 - t194, -t123 * t216 ^ 2 - t195 - t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t87, -t87 ^ 2 + t88 ^ 2 (-t87 - t179) * qJD(4), t97 + (t88 - t158) * qJD(4), 0, -t95 * t88 - t135, qJD(3) * t179 - t95 * t87, t219, t212, t229, t214, 0, -t59 * t199 + t22 * t114 + (t156 * t121 - t182) * qJD(5) + t201, -t140 * t199 + t23 * t114 + t221 + (t156 * qJD(5) - t35) * t124 - t151, t225, t188 * t216 + t226, t224, t136 + t194, t209, t109 * t17 + t143 * t51 + (t146 * t216 - t3) * t123 + t130 * t120 + t203, t109 * t16 + t123 * t130 + t143 * t53 - t146 * t217 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, t212, t229, t214, 0 (-qJD(5) + t114) * t21 + t201, t20 * t114 - t2 + t221, t225, t217 * t53 + t226, t224, -t216 * t217 + t194 + t31, t209, -pkin(5) * t17 - t3 * t123 + (-t120 * t20 + t123 * t38) * t216 - t21 * t51 + t120 * t223 - t190 * pkin(9) + t203, -pkin(5) * t16 - (t120 * t38 + t123 * t20) * t216 - t21 * t53 + t18 * t220 + t208 * pkin(9) + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t51, -t51 ^ 2 + t53 ^ 2, -t216 * t51 + t16, -t216 * t53 - t17, t34, -t120 * t2 - t18 * t53 - t228 * t5 + t8, -t120 * t9 - t123 * t2 + t141 * t228 + t18 * t51;];
tauc_reg  = t1;
