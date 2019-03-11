% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:15
% EndTime: 2019-03-08 21:32:23
% DurationCPUTime: 2.45s
% Computational Cost: add. (4134->337), mult. (10590->467), div. (0->0), fcn. (8044->10), ass. (0->171)
t125 = sin(pkin(11));
t127 = cos(pkin(11));
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t212 = -qJ(4) - pkin(8);
t165 = qJD(3) * t212;
t95 = qJD(4) * t133 + t130 * t165;
t96 = -t130 * qJD(4) + t133 * t165;
t55 = t125 * t96 + t127 * t95;
t192 = t127 * t133;
t106 = t125 * t130 - t192;
t134 = cos(qJ(2));
t126 = sin(pkin(6));
t187 = qJD(1) * t126;
t169 = t134 * t187;
t79 = t106 * t169;
t207 = t55 + t79;
t185 = qJD(2) * t130;
t99 = qJD(2) * t192 - t125 * t185;
t94 = qJD(5) - t99;
t107 = t125 * t133 + t127 * t130;
t208 = -t107 * t169 + t125 * t95 - t127 * t96;
t132 = cos(qJ(5));
t181 = qJD(5) * t132;
t129 = sin(qJ(5));
t100 = t107 * qJD(3);
t92 = qJD(2) * t100;
t87 = t129 * t92;
t227 = -t181 * t94 - t87;
t131 = sin(qJ(2));
t170 = t131 * t187;
t182 = qJD(5) * t129;
t103 = t106 * qJD(3);
t180 = t130 * qJD(3);
t177 = pkin(3) * t180;
t57 = pkin(4) * t100 + pkin(9) * t103 + t177;
t173 = -pkin(3) * t133 - pkin(2);
t66 = pkin(4) * t106 - pkin(9) * t107 + t173;
t111 = t212 * t130;
t112 = t212 * t133;
t81 = t111 * t125 - t112 * t127;
t226 = -t66 * t181 + t81 * t182 - t207 * t132 + (t170 - t57) * t129;
t88 = t132 * t92;
t225 = -t182 * t94 + t88;
t109 = qJD(2) * pkin(8) + t170;
t152 = qJD(4) + t169;
t128 = cos(pkin(6));
t186 = qJD(1) * t128;
t168 = t130 * t186;
t224 = (-t109 * t133 - t168) * qJD(3) + (-qJD(3) * t133 * qJ(4) - t130 * t152) * qJD(2);
t161 = qJ(4) * qJD(2) + t109;
t77 = t133 * t161 + t168;
t67 = t125 * t77;
t116 = t133 * t186;
t76 = -t130 * t161 + t116;
t71 = qJD(3) * pkin(3) + t76;
t30 = t127 * t71 - t67;
t27 = -qJD(3) * pkin(4) - t30;
t101 = t107 * qJD(2);
t179 = t132 * qJD(3);
t82 = t101 * t129 - t179;
t84 = qJD(3) * t129 + t101 * t132;
t15 = pkin(5) * t82 - qJ(6) * t84 + t27;
t119 = pkin(3) * t125 + pkin(9);
t203 = t119 * t92;
t223 = t15 * t94 - t203;
t178 = qJD(2) * qJD(3);
t166 = t133 * t178;
t167 = t130 * t178;
t148 = -t125 * t167 + t127 * t166;
t50 = qJD(5) * t84 + t129 * t148;
t222 = t84 ^ 2;
t221 = t94 ^ 2;
t202 = t127 * t77;
t31 = t125 * t71 + t202;
t28 = qJD(3) * pkin(9) + t31;
t93 = qJD(2) * t173 + qJD(4) - t169;
t44 = -pkin(4) * t99 - pkin(9) * t101 + t93;
t13 = t129 * t44 + t132 * t28;
t7 = qJ(6) * t94 + t13;
t220 = t7 * t94;
t219 = t92 * pkin(5);
t218 = qJ(6) * t100 + qJD(6) * t106 - t226;
t206 = t129 * t66 + t132 * t81;
t58 = -t129 * t79 - t132 * t170;
t217 = -t100 * pkin(5) + qJD(5) * t206 + t129 * t55 - t132 * t57 - t58;
t153 = pkin(5) * t129 - qJ(6) * t132;
t154 = pkin(5) * t132 + qJ(6) * t129;
t216 = t153 * t103 - (qJD(5) * t154 - qJD(6) * t132) * t107 - t208;
t215 = t13 * t94;
t214 = t82 * t99;
t213 = t84 * t82;
t164 = t84 * t94;
t33 = t125 * t76 + t202;
t211 = t129 * qJD(6) - t153 * t94 + t33;
t35 = t127 * t76 - t67;
t56 = pkin(3) * t185 + pkin(4) * t101 - pkin(9) * t99;
t210 = t129 * t56 + t132 * t35;
t209 = -t129 * t50 - t181 * t82;
t205 = qJD(2) * pkin(2);
t204 = t101 * t82;
t201 = t129 * t94;
t49 = -qJD(5) * t179 + t101 * t182 - t132 * t148;
t200 = t49 * t129;
t199 = t84 * t101;
t198 = t92 * qJ(6);
t197 = t103 * t129;
t196 = t103 * t132;
t195 = t126 * t131;
t194 = t126 * t134;
t136 = qJD(2) ^ 2;
t193 = t126 * t136;
t135 = qJD(3) ^ 2;
t191 = t135 * t130;
t190 = t135 * t133;
t12 = -t129 * t28 + t132 * t44;
t189 = qJD(6) - t12;
t98 = pkin(3) * t167 + qJD(2) * t170;
t188 = t130 ^ 2 - t133 ^ 2;
t184 = qJD(2) * t131;
t183 = qJD(5) * t119;
t176 = t94 * t183;
t174 = t131 * t193;
t120 = -pkin(3) * t127 - pkin(4);
t172 = t126 * t184;
t171 = qJD(2) * t194;
t48 = (-t109 * t130 + t116) * qJD(3) + (-qJ(4) * t180 + t133 * t152) * qJD(2);
t18 = t125 * t48 - t127 * t224;
t19 = t125 * t224 + t127 * t48;
t41 = pkin(4) * t92 - pkin(9) * t148 + t98;
t163 = t129 * t19 - t132 * t41 + t181 * t28 + t182 * t44;
t80 = -t111 * t127 - t125 * t112;
t160 = t130 * t171;
t159 = t133 * t171;
t3 = pkin(5) * t50 + qJ(6) * t49 - qJD(6) * t84 + t18;
t158 = -t3 - t176;
t6 = -pkin(5) * t94 + t189;
t157 = -t129 * t7 + t132 * t6;
t155 = t18 * t107 - t81 * t92;
t151 = -t132 * t94 * t99 - t227;
t150 = t201 * t99 + t225;
t149 = t15 * t84 + t163;
t104 = t128 * t130 + t133 * t195;
t147 = t128 * t133 - t130 * t195;
t62 = t104 * t127 + t125 * t147;
t46 = t129 * t62 + t132 * t194;
t47 = -t129 * t194 + t132 * t62;
t146 = t129 * t41 + t132 * t19 + t181 * t44 - t182 * t28;
t74 = qJD(3) * t147 + t159;
t75 = -qJD(3) * t104 - t160;
t34 = t125 * t75 + t127 * t74;
t10 = qJD(5) * t47 + t129 * t34 - t132 * t172;
t32 = t125 * t74 - t127 * t75;
t61 = t104 * t125 - t127 * t147;
t144 = -t10 * t94 + t32 * t82 - t46 * t92 + t50 * t61;
t143 = t205 * qJD(2);
t142 = t27 * t94 - t203;
t139 = -0.2e1 * qJD(3) * t205;
t9 = -qJD(5) * t46 + t129 * t172 + t132 * t34;
t138 = t32 * t84 - t47 * t92 - t49 * t61 - t9 * t94;
t105 = t120 - t154;
t43 = pkin(5) * t84 + qJ(6) * t82;
t36 = t107 * t153 + t80;
t25 = -pkin(5) * t106 + t129 * t81 - t132 * t66;
t24 = qJ(6) * t106 + t206;
t21 = t82 * t94 - t49;
t14 = -pkin(5) * t101 + t129 * t35 - t132 * t56;
t11 = qJ(6) * t101 + t210;
t2 = t163 - t219;
t1 = qJD(6) * t94 + t146 + t198;
t4 = [0, 0, -t174, -t134 * t193, 0, 0, 0, 0, 0, -t133 * t174 + (t75 - t160) * qJD(3), t130 * t174 + (-t74 - t159) * qJD(3), t101 * t32 + t148 * t61 + t34 * t99 - t62 * t92, t18 * t61 + t19 * t62 - t30 * t32 + t31 * t34 + (-t134 * t98 + t184 * t93) * t126, 0, 0, 0, 0, 0, t144, t138, t144, t10 * t84 - t46 * t49 - t47 * t50 - t82 * t9, -t138, t1 * t47 + t10 * t6 + t15 * t32 + t2 * t46 + t3 * t61 + t7 * t9; 0, 0, 0, 0, 0.2e1 * t130 * t166, -0.2e1 * t188 * t178, t190, -t191, 0, -pkin(8) * t190 + t130 * t139, pkin(8) * t191 + t133 * t139, -t31 * t100 + t101 * t208 + t30 * t103 - t19 * t106 + t148 * t80 + t207 * t99 + t155, t19 * t81 + t18 * t80 + t98 * t173 + (-t170 + t177) * t93 + t207 * t31 - t208 * t30, -t84 * t196 + (-t49 * t132 - t182 * t84) * t107 -(-t129 * t84 - t132 * t82) * t103 + (t200 - t132 * t50 + (t129 * t82 - t132 * t84) * qJD(5)) * t107, t84 * t100 - t49 * t106 + t107 * t225 - t196 * t94, -t82 * t100 - t50 * t106 + t107 * t227 + t197 * t94, t100 * t94 + t106 * t92, -t163 * t106 + t12 * t100 + t80 * t50 + t58 * t94 + t208 * t82 + ((-qJD(5) * t81 + t57) * t94 + t66 * t92 + t27 * qJD(5) * t107) * t132 + ((-qJD(5) * t66 - t55) * t94 - t27 * t103 + t155) * t129, -t206 * t92 - t146 * t106 - t13 * t100 - t80 * t49 - t27 * t196 + t226 * t94 + t208 * t84 + (t18 * t132 - t182 * t27) * t107, -t15 * t197 - t6 * t100 - t2 * t106 - t25 * t92 + t36 * t50 - t217 * t94 - t216 * t82 + (t3 * t129 + t15 * t181) * t107, -t24 * t50 - t25 * t49 + t217 * t84 - t218 * t82 - t157 * t103 + (-t1 * t129 + t2 * t132 + (-t129 * t6 - t132 * t7) * qJD(5)) * t107, t15 * t196 + t1 * t106 + t7 * t100 + t24 * t92 + t36 * t49 + t218 * t94 + t216 * t84 + (-t3 * t132 + t15 * t182) * t107, t1 * t24 - t15 * t216 + t2 * t25 + t217 * t6 + t218 * t7 + t3 * t36; 0, 0, 0, 0, -t130 * t136 * t133, t188 * t136, 0, 0, 0, t130 * t143, t133 * t143 (-t35 + t30) * t99 + (t31 - t33) * t101 + (-t125 * t92 - t127 * t148) * pkin(3), t30 * t33 - t31 * t35 + (t125 * t19 - t127 * t18 - t185 * t93) * pkin(3), t132 * t164 - t200 (-t49 + t214) * t132 - t84 * t201 + t209, t151 - t199, t150 + t204, -t94 * t101, -t12 * t101 + t120 * t50 - t33 * t82 + (-t18 + (-t56 - t183) * t94) * t132 + (t35 * t94 + t142) * t129, -t120 * t49 + t210 * t94 + t13 * t101 - t33 * t84 + (t18 + t176) * t129 + t142 * t132, t6 * t101 + t105 * t50 + t129 * t223 + t158 * t132 + t14 * t94 - t211 * t82, t11 * t82 - t14 * t84 + (-t119 * t50 - t6 * t99 + t1 + (t119 * t84 + t6) * qJD(5)) * t132 + (-t119 * t49 + t7 * t99 + t2 + (t119 * t82 - t7) * qJD(5)) * t129, -t7 * t101 + t105 * t49 - t11 * t94 + t158 * t129 - t132 * t223 + t211 * t84, t3 * t105 - t7 * t11 - t6 * t14 - t211 * t15 + (qJD(5) * t157 + t1 * t132 + t2 * t129) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 ^ 2 - t99 ^ 2, t101 * t30 - t31 * t99 + t98, 0, 0, 0, 0, 0, t150 - t204, -t132 * t221 - t199 - t87, -t201 * t94 - t204 + t88 (t49 + t214) * t132 + t129 * t164 + t209, t151 + t199, -t15 * t101 + (-t2 + t220) * t132 + (t6 * t94 + t1) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t82 ^ 2 + t222, t21, t164 - t50, t92, -t27 * t84 - t163 + t215, t12 * t94 + t27 * t82 - t146, -t43 * t82 - t149 + t215 + 0.2e1 * t219, pkin(5) * t49 - t50 * qJ(6) + (-t13 + t7) * t84 + (t6 - t189) * t82, 0.2e1 * t198 - t15 * t82 + t43 * t84 + (0.2e1 * qJD(6) - t12) * t94 + t146, -t2 * pkin(5) + t1 * qJ(6) - t6 * t13 - t15 * t43 + t189 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t101 + t213, t21, -t221 - t222, t149 - t219 - t220;];
tauc_reg  = t4;
