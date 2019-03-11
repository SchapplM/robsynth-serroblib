% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:52
% EndTime: 2019-03-09 04:56:00
% DurationCPUTime: 2.52s
% Computational Cost: add. (2942->372), mult. (5964->484), div. (0->0), fcn. (3312->4), ass. (0->183)
t114 = sin(qJ(3));
t113 = sin(qJ(4));
t116 = cos(qJ(3));
t184 = qJD(4) * t116;
t169 = t113 * t184;
t115 = cos(qJ(4));
t182 = t115 * qJD(3);
t241 = t114 * t182 + t169;
t83 = pkin(3) * t114 - pkin(8) * t116 + qJ(2);
t61 = t83 * qJD(1);
t117 = -pkin(1) - pkin(7);
t98 = t117 * qJD(1) + qJD(2);
t82 = t114 * t98;
t63 = qJD(3) * pkin(8) + t82;
t26 = t113 * t63 - t115 * t61;
t190 = qJD(3) * t113;
t191 = qJD(1) * t116;
t75 = t115 * t191 + t190;
t140 = pkin(5) * t75 + t26;
t198 = qJD(5) + t140;
t193 = qJD(1) * t114;
t103 = qJD(4) + t193;
t153 = qJD(4) * t114 + qJD(1);
t171 = t116 * t182;
t126 = t153 * t113 - t171;
t189 = qJD(3) * t114;
t39 = t241 * qJD(1) - qJD(4) * t182;
t131 = -t116 * t39 - t75 * t189;
t181 = qJD(1) * qJD(3);
t106 = t116 * t181;
t151 = t114 * t106;
t122 = -t103 * t126 + t115 * t151 + t131;
t179 = 0.2e1 * qJD(1);
t27 = t113 * t61 + t115 * t63;
t20 = -qJ(5) * t103 - t27;
t152 = pkin(4) * t106;
t188 = qJD(3) * t116;
t173 = t113 * t188;
t185 = qJD(4) * t115;
t186 = qJD(4) * t113;
t149 = pkin(3) * t116 + pkin(8) * t114;
t72 = t149 * qJD(3) + qJD(2);
t53 = t72 * qJD(1);
t158 = -t115 * t53 + t98 * t173 + t63 * t185 + t61 * t186;
t6 = -t152 + t158;
t240 = t103 * t20 + t6;
t101 = qJ(5) * t106;
t90 = t103 * qJD(5);
t239 = t101 + t90;
t100 = t103 ^ 2;
t71 = t75 ^ 2;
t238 = -t71 - t100;
t208 = qJD(4) * t75;
t174 = t113 * t189;
t91 = qJD(1) * t174;
t40 = -t91 + t208;
t192 = qJD(1) * t115;
t210 = qJ(5) * t114;
t236 = -pkin(4) * t186 + qJD(5) * t113 + t192 * t210 + t82;
t73 = t113 * t191 - t182;
t235 = t73 ^ 2;
t234 = pkin(5) + pkin(8);
t233 = pkin(5) * t73;
t133 = qJ(5) * t39 - qJD(5) * t75 + t98 * t189;
t7 = pkin(4) * t40 + t133;
t232 = t113 * t7;
t231 = t115 * t7;
t213 = t116 * t98;
t64 = -qJD(3) * pkin(3) - t213;
t129 = -qJ(5) * t75 + t64;
t227 = pkin(4) + qJ(6);
t13 = t227 * t73 + t129;
t230 = t13 * t75;
t24 = pkin(4) * t73 + t129;
t229 = t24 * t75;
t228 = t75 * t73;
t204 = t114 * t115;
t130 = -pkin(5) * t204 - t227 * t116;
t205 = t113 * t116;
t78 = t149 * qJD(1);
t163 = t115 * t78 - t98 * t205;
t89 = t234 * t115;
t226 = t130 * qJD(1) - qJD(4) * t89 - t163;
t206 = t113 * t114;
t202 = t115 * t116;
t222 = t113 * t78 + t98 * t202;
t225 = (pkin(5) * t206 + qJ(5) * t116) * qJD(1) + t222 + t234 * t186;
t209 = qJ(5) * t115;
t143 = qJ(6) * t113 - t209;
t166 = t227 * t114;
t224 = -t113 * qJD(1) * t166 - t143 * qJD(4) + qJD(6) * t115 + t236;
t223 = -pkin(4) * t113 * t193 + qJ(5) * t185 + t236;
t203 = t114 * t117;
t220 = t113 * t83 + t115 * t203;
t219 = qJ(5) * t40;
t10 = -t227 * t103 + t198;
t218 = t10 * t103;
t216 = t103 * t73;
t215 = t103 * t75;
t212 = t39 * t113;
t211 = t73 * qJ(5);
t207 = t103 * t115;
t118 = qJD(3) ^ 2;
t201 = t118 * t114;
t200 = t118 * t116;
t119 = qJD(1) ^ 2;
t199 = t119 * qJ(2);
t197 = qJD(5) + t26;
t16 = t27 - t233;
t196 = -qJD(6) - t16;
t111 = t116 ^ 2;
t195 = t114 ^ 2 - t111;
t194 = -t118 - t119;
t187 = qJD(3) * t117;
t183 = qJD(4) * t117;
t170 = t116 * t187;
t180 = t113 * t72 + t115 * t170 + t83 * t185;
t178 = pkin(8) * t103 * t113;
t177 = pkin(8) * t207;
t176 = pkin(8) * t188;
t175 = qJD(2) * t179;
t168 = t114 * t183;
t167 = t115 * t184;
t165 = -t113 * qJ(5) - pkin(3);
t164 = -t106 + t228;
t95 = t113 * t203;
t162 = t115 * t83 - t95;
t161 = -t64 + t213;
t160 = qJD(3) * t227;
t159 = -t113 * t53 - t98 * t171 - t61 * t185 + t63 * t186;
t157 = pkin(4) * t167 + t241 * qJ(5) + t114 * t187;
t156 = -t75 + t190;
t155 = t73 + t182;
t150 = -t90 + t159;
t37 = -t210 - t220;
t12 = qJD(6) - t20 - t233;
t148 = t10 * t115 - t113 * t12;
t147 = t10 * t113 + t115 * t12;
t19 = -pkin(4) * t103 + t197;
t146 = t113 * t20 + t115 * t19;
t145 = t113 * t19 - t115 * t20;
t144 = -t113 * t170 - t83 * t186 + (-t168 + t72) * t115;
t142 = (t103 + t193) * t116;
t141 = qJD(1) * t111 - t103 * t114;
t139 = -pkin(5) * t39 + t158;
t138 = -pkin(5) * t40 - t159;
t137 = -t114 * t24 + t176;
t136 = t114 * t64 - t176;
t2 = qJD(6) * t73 + t227 * t40 + t133;
t135 = t113 * t2 + t13 * t185;
t134 = -t115 * t2 + t13 * t186;
t132 = t103 * t27 - t158;
t127 = -t13 * t73 + t138;
t4 = -t101 + t150;
t125 = t146 * qJD(4) + t113 * t6 - t115 * t4;
t124 = -t160 * t191 + t139;
t123 = t215 - t40;
t21 = t216 - t39;
t121 = t126 * t73 - t40 * t204 - t39 * t206 + (t114 * t185 + t173 + t192) * t75;
t120 = t73 * t189 + (-t40 - t91) * t116 + (-t153 * t115 - t173) * t103;
t105 = pkin(4) * t205;
t99 = 0.2e1 * t101;
t88 = t234 * t113;
t84 = -pkin(4) * t115 + t165;
t67 = -t227 * t115 + t165;
t43 = t105 + (-t117 - t209) * t116;
t38 = -pkin(4) * t114 - t162;
t36 = t105 + (-t117 + t143) * t116;
t33 = pkin(4) * t75 + t211;
t31 = -pkin(5) * t205 - t37;
t30 = -pkin(4) * t191 - t163;
t29 = -qJ(5) * t191 - t222;
t25 = t95 + (pkin(5) * t116 - t83) * t115 - t166;
t22 = t227 * t75 + t211;
t18 = -pkin(4) * t174 - qJD(5) * t202 + t157;
t14 = -pkin(4) * t188 - t144;
t11 = -qJ(5) * t188 + (t113 * t183 - qJD(5)) * t114 - t180;
t9 = (qJ(6) * qJD(4) - qJD(5)) * t202 + (qJD(6) * t116 - t114 * t160) * t113 + t157;
t8 = (-pkin(5) * t185 + qJ(5) * qJD(3)) * t116 + (qJD(5) + (pkin(5) * qJD(3) - t183) * t113) * t114 + t180;
t5 = -pkin(5) * t169 + t130 * qJD(3) - qJD(6) * t114 - t144;
t3 = t138 + t239;
t1 = -qJD(6) * t103 + t124;
t15 = [0, 0, 0, 0, t175, qJ(2) * t175, -0.2e1 * t151, 0.2e1 * t195 * t181, -t201, -t200, 0, -t117 * t201 + (qJ(2) * t188 + qJD(2) * t114) * t179, -t117 * t200 + (-qJ(2) * t189 + qJD(2) * t116) * t179, t131 * t115 - t75 * t169 (t113 * t75 + t115 * t73) * t189 + (t212 - t115 * t40 + (t113 * t73 - t115 * t75) * qJD(4)) * t116, -t103 * t169 - t114 * t39 + (t141 * t115 + t116 * t75) * qJD(3), -t103 * t167 - t114 * t40 + (-t141 * t113 - t116 * t73) * qJD(3), qJD(3) * t142, t144 * t103 - t158 * t114 + (-t117 * t40 + t185 * t64) * t116 + ((qJD(1) * t162 - t26) * t116 + (t113 * t161 + t117 * t73) * t114) * qJD(3) -(-t113 * t168 + t180) * t103 + t159 * t114 + (t117 * t39 - t186 * t64) * t116 + ((-qJD(1) * t220 - t27) * t116 + (t115 * t161 + t117 * t75) * t114) * qJD(3), t11 * t73 + t14 * t75 + t37 * t40 - t38 * t39 - t146 * t189 + (-qJD(4) * t145 + t113 * t4 + t115 * t6) * t116, t103 * t14 - t18 * t73 - t40 * t43 + (t24 * t190 + t6) * t114 + (-t24 * t185 - t232 + (qJD(1) * t38 + t19) * qJD(3)) * t116, -t103 * t11 - t18 * t75 + t39 * t43 + (t24 * t182 - t4) * t114 + (t24 * t186 - t231 + (-qJD(1) * t37 - t20) * qJD(3)) * t116, t11 * t20 + t14 * t19 + t18 * t24 + t37 * t4 + t38 * t6 + t43 * t7, -t25 * t39 - t31 * t40 + t5 * t75 - t73 * t8 - t148 * t189 + (-qJD(4) * t147 + t1 * t115 - t113 * t3) * t116, t103 * t8 + t36 * t39 - t75 * t9 + (t13 * t182 + t3) * t114 + ((qJD(1) * t31 + t12) * qJD(3) + t134) * t116, -t103 * t5 + t36 * t40 + t73 * t9 + (-t13 * t190 - t1) * t114 + ((-qJD(1) * t25 - t10) * qJD(3) + t135) * t116, t1 * t25 + t10 * t5 + t12 * t8 + t13 * t9 + t2 * t36 + t3 * t31; 0, 0, 0, 0, -t119, -t199, 0, 0, 0, 0, 0, t194 * t114, t194 * t116, 0, 0, 0, 0, 0, t120, -t122, t121, t116 * t40 + t153 * t207 + (t113 * t142 - t114 * t73) * qJD(3), t122, t146 * qJD(1) + (qJD(3) * t145 - t7) * t116 + (qJD(3) * t24 + t125) * t114, t121, t122, t120, t148 * qJD(1) + (qJD(3) * t147 - t2) * t116 + (qJD(3) * t13 + qJD(4) * t148 + t1 * t113 + t115 * t3) * t114; 0, 0, 0, 0, 0, 0, t116 * t119 * t114, -t195 * t119, 0, 0, 0, -t116 * t199, t114 * t199, t75 * t207 - t212 (-t39 - t216) * t115 + (-t40 - t215) * t113, t103 * t185 + (t103 * t204 + t156 * t116) * qJD(1), -t103 * t186 + (-t103 * t206 + t155 * t116) * qJD(1), -t103 * t191, -pkin(3) * t40 - t163 * t103 - t155 * t82 + (t113 * t64 - t177) * qJD(4) + (t113 * t136 + t26 * t116) * qJD(1), pkin(3) * t39 + t222 * t103 + t156 * t82 + (t115 * t64 + t178) * qJD(4) + (t115 * t136 + t27 * t116) * qJD(1), -t29 * t73 - t30 * t75 + (-t4 + t103 * t19 + (-t40 + t208) * pkin(8)) * t115 + ((qJD(4) * t73 - t39) * pkin(8) + t240) * t113, -t103 * t30 + t231 - t40 * t84 + t223 * t73 + (-t113 * t24 + t177) * qJD(4) + (t113 * t137 - t116 * t19) * qJD(1), t103 * t29 - t232 + t39 * t84 + t223 * t75 + (-t115 * t24 - t178) * qJD(4) + (t115 * t137 + t116 * t20) * qJD(1), pkin(8) * t125 - t19 * t30 - t20 * t29 - t223 * t24 + t7 * t84, -t39 * t88 - t40 * t89 - t226 * t75 + t225 * t73 + (t3 + t218) * t115 + (-t103 * t12 + t1) * t113, t39 * t67 + t224 * t75 - t225 * t103 + (-t13 * t204 + (qJD(3) * t89 - t12) * t116) * qJD(1) - t135, t40 * t67 - t224 * t73 + t226 * t103 + (t13 * t206 + (-qJD(3) * t88 + t10) * t116) * qJD(1) + t134, t1 * t88 - t10 * t226 - t12 * t225 - t13 * t224 + t2 * t67 + t3 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t71 - t235, t21, t123, t106, -t64 * t75 + t132, -t103 * t26 + t64 * t73 + t159, pkin(4) * t39 - t219 + (-t20 - t27) * t75 + (t19 - t197) * t73, t33 * t73 - t132 - 0.2e1 * t152 + t229, t103 * t197 - t24 * t73 + t33 * t75 - t150 + t99, -pkin(4) * t6 - qJ(5) * t4 - t19 * t27 - t197 * t20 - t24 * t33, -t219 + t227 * t39 + (t12 + t196) * t75 + (t10 - t198) * t73, t103 * t140 + t22 * t75 + t127 + 0.2e1 * t90 + t99, -t230 - t22 * t73 + (0.2e1 * qJD(6) + t16) * t103 + 0.2e1 * t227 * t106 - t139, qJ(5) * t3 - t1 * t227 + t10 * t196 + t12 * t198 - t13 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t164, t238, t229 + t240, t21, t238, t164, t230 + (-qJD(6) - t12) * t103 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t106 + t228, -t100 - t235, t127 + t218 + t239;];
tauc_reg  = t15;
