% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:28
% EndTime: 2019-03-09 04:48:36
% DurationCPUTime: 2.97s
% Computational Cost: add. (4364->353), mult. (9381->490), div. (0->0), fcn. (5837->6), ass. (0->175)
t148 = sin(pkin(9));
t149 = cos(pkin(9));
t150 = sin(qJ(4));
t152 = cos(qJ(4));
t248 = -t148 * t150 + t149 * t152;
t114 = t148 * t152 + t149 * t150;
t100 = t114 * qJD(1);
t151 = sin(qJ(3));
t99 = t114 * qJD(4);
t235 = t151 * t100 + t99;
t200 = t151 * qJD(1);
t139 = qJD(4) + t200;
t199 = t152 * qJD(3);
t153 = cos(qJ(3));
t208 = qJD(1) * t153;
t116 = t150 * t208 - t199;
t179 = t152 * t208;
t201 = t150 * qJD(3);
t118 = t179 + t201;
t69 = t149 * t116 + t148 * t118;
t247 = t139 * t69;
t204 = qJD(4) * t151;
t206 = qJD(3) * t153;
t246 = t114 * t204 - t206 * t248 + t100;
t203 = qJD(4) * t152;
t205 = qJD(4) * t150;
t101 = -t148 * t205 + t149 * t203;
t243 = t248 * qJD(1);
t231 = t243 * t151 + t101;
t166 = -t148 * t116 + t149 * t118;
t245 = t166 ^ 2;
t196 = 0.2e1 * qJD(1);
t241 = -qJ(5) - pkin(8);
t175 = qJD(4) * t241;
t159 = -t150 * qJD(5) + t152 * t175;
t168 = pkin(3) * t153 + pkin(8) * t151;
t120 = t168 * qJD(1);
t105 = t152 * t120;
t154 = -pkin(1) - pkin(7);
t134 = t154 * qJD(1) + qJD(2);
t218 = t153 * t134;
t221 = t151 * t152;
t53 = -t150 * t218 + t105 + (pkin(4) * t153 + qJ(5) * t221) * qJD(1);
t190 = t150 * t200;
t212 = t150 * t120 + t152 * t218;
t62 = qJ(5) * t190 + t212;
t198 = t152 * qJD(5);
t93 = t150 * t175 + t198;
t236 = (t159 - t53) * t149 + (t62 - t93) * t148;
t107 = -qJD(3) * pkin(3) - t218;
t74 = t116 * pkin(4) + qJD(5) + t107;
t23 = t69 * pkin(5) - qJ(6) * t166 + t74;
t244 = t23 * t166;
t122 = t151 * pkin(3) - t153 * pkin(8) + qJ(2);
t220 = t151 * t154;
t211 = t150 * t122 + t152 * t220;
t189 = t151 * t201;
t80 = -qJD(1) * t189 + qJD(4) * t118;
t103 = t122 * qJD(1);
t121 = t151 * t134;
t106 = qJD(3) * pkin(8) + t121;
t65 = t150 * t103 + t152 * t106;
t51 = -t116 * qJ(5) + t65;
t47 = t149 * t51;
t64 = t152 * t103 - t150 * t106;
t50 = -t118 * qJ(5) + t64;
t21 = t148 * t50 + t47;
t242 = t21 * t166;
t112 = t168 * qJD(3) + qJD(2);
t92 = t112 * qJD(1);
t82 = t152 * t92;
t158 = -qJD(4) * t65 + t82;
t180 = t151 * t199;
t202 = qJD(4) * t153;
t186 = t150 * t202;
t160 = -t180 - t186;
t79 = qJD(1) * t160 + qJD(4) * t199;
t15 = -t79 * qJ(5) - t118 * qJD(5) + (pkin(4) * qJD(1) - t134 * t150) * t206 + t158;
t187 = t153 * t199;
t161 = t103 * t203 - t106 * t205 + t134 * t187 + t150 * t92;
t20 = -t80 * qJ(5) - t116 * qJD(5) + t161;
t3 = -t148 * t20 + t149 * t15;
t4 = t148 * t15 + t149 * t20;
t27 = t148 * t53 + t149 * t62;
t24 = qJ(6) * t208 + t27;
t60 = t148 * t159 + t149 * t93;
t240 = -t24 + t60;
t239 = -pkin(5) * t208 + t236;
t222 = t150 * t154;
t176 = pkin(4) - t222;
t95 = t152 * t112;
t31 = qJ(5) * t180 + t95 - t211 * qJD(4) + (qJ(5) * t205 + t176 * qJD(3) - t198) * t153;
t185 = t152 * t202;
t193 = t150 * t112 + t122 * t203 + t154 * t187;
t34 = -qJ(5) * t185 + (-qJD(5) * t153 + (qJ(5) * qJD(3) - qJD(4) * t154) * t151) * t150 + t193;
t9 = t148 * t31 + t149 * t34;
t238 = t114 * qJD(6) + t121 + t231 * qJ(6) - t235 * pkin(5) + (-t190 - t205) * pkin(4);
t43 = t139 * pkin(4) + t50;
t19 = t148 * t43 + t47;
t188 = t153 * t201;
t237 = -t101 * t151 - t148 * t187 - t149 * t188 - t243;
t111 = t152 * t122;
t219 = t152 * t153;
t67 = -qJ(5) * t219 + t176 * t151 + t111;
t223 = t150 * t153;
t75 = -qJ(5) * t223 + t211;
t40 = t148 * t67 + t149 * t75;
t234 = t148 * t51;
t233 = t79 * t150;
t230 = t107 * t150;
t229 = t116 * t139;
t228 = t118 * t139;
t227 = t139 * t152;
t224 = t150 * t139;
t217 = t153 * t154;
t155 = qJD(3) ^ 2;
t216 = t155 * t151;
t215 = t155 * t153;
t156 = qJD(1) ^ 2;
t214 = t156 * qJ(2);
t22 = t149 * t50 - t234;
t213 = qJD(6) - t22;
t147 = t153 ^ 2;
t210 = t151 ^ 2 - t147;
t209 = -t155 - t156;
t207 = qJD(3) * t151;
t197 = qJD(1) * qJD(3);
t177 = t153 * t197;
t195 = qJ(6) * t177 + t4;
t194 = qJD(2) * t196;
t192 = t150 * t220;
t191 = -t152 * pkin(4) - pkin(3);
t61 = t80 * pkin(4) + t134 * t207;
t178 = t241 * t150;
t45 = t148 * t79 + t149 * t80;
t174 = pkin(4) * t223 - t217;
t173 = t139 + t200;
t172 = t116 + t199;
t171 = -t118 + t201;
t170 = qJD(1) + t204;
t46 = -t148 * t80 + t149 * t79;
t129 = t241 * t152;
t77 = -t148 * t129 - t149 * t178;
t78 = -t149 * t129 + t148 * t178;
t169 = -t78 * t45 + t77 * t46 - t60 * t69;
t167 = -t69 ^ 2 - t245;
t8 = -t148 * t34 + t149 * t31;
t18 = t149 * t43 - t234;
t39 = -t148 * t75 + t149 * t67;
t165 = qJD(1) * t147 - t139 * t151;
t164 = -pkin(8) * t206 + t107 * t151;
t2 = -pkin(5) * t177 - t3;
t162 = t154 * t207 + (t185 - t189) * pkin(4);
t87 = t114 * t151;
t89 = t248 * t151;
t157 = -t166 * t237 + t246 * t69 - t89 * t45 + t87 * t46;
t5 = t45 * pkin(5) - t46 * qJ(6) - qJD(6) * t166 + t61;
t143 = -t149 * pkin(4) - pkin(5);
t141 = t148 * pkin(4) + qJ(6);
t90 = -t148 * t223 + t149 * t219;
t88 = t114 * t153;
t66 = -pkin(5) * t248 - t114 * qJ(6) + t191;
t58 = -t148 * t189 + t149 * t180 + t153 * t99;
t56 = qJD(3) * t87 - t202 * t248;
t49 = t88 * pkin(5) - t90 * qJ(6) + t174;
t35 = -t151 * pkin(5) - t39;
t33 = t151 * qJ(6) + t40;
t28 = t118 * pkin(4) + pkin(5) * t166 + qJ(6) * t69;
t12 = t139 * qJ(6) + t19;
t11 = -t139 * pkin(5) + qJD(6) - t18;
t10 = -t56 * pkin(5) + t58 * qJ(6) - t90 * qJD(6) + t162;
t7 = -pkin(5) * t206 - t8;
t6 = qJ(6) * t206 + t151 * qJD(6) + t9;
t1 = t139 * qJD(6) + t195;
t13 = [0, 0, 0, 0, t194, qJ(2) * t194, -0.2e1 * t151 * t177, 0.2e1 * t210 * t197, -t216, -t215, 0, -t154 * t216 + (qJ(2) * t206 + qJD(2) * t151) * t196, -t154 * t215 + (-qJ(2) * t207 + qJD(2) * t153) * t196, t118 * t160 + t219 * t79 (t116 * t152 + t118 * t150) * t207 + (-t233 - t152 * t80 + (t116 * t150 - t118 * t152) * qJD(4)) * t153, -t139 * t186 + t79 * t151 + (t118 * t153 + t152 * t165) * qJD(3), -t139 * t185 - t80 * t151 + (-t116 * t153 - t150 * t165) * qJD(3), t173 * t206, -t80 * t217 + t95 * t139 + t82 * t151 + (t107 * t219 - t139 * t211 - t151 * t65) * qJD(4) + ((t154 * t116 - t230) * t151 + (-t139 * t222 + (t111 - t192) * qJD(1) + t64) * t153) * qJD(3) -(-qJD(4) * t192 + t193) * t139 - t161 * t151 + (-t107 * t205 - t154 * t79) * t153 + ((-qJD(1) * t211 - t65) * t153 + (t154 * t118 + (-t107 + t218) * t152) * t151) * qJD(3), -t166 * t8 + t18 * t58 + t19 * t56 - t3 * t90 - t39 * t46 - t4 * t88 - t40 * t45 - t9 * t69, t162 * t74 + t174 * t61 + t18 * t8 + t19 * t9 + t3 * t39 + t4 * t40, t10 * t69 - t7 * t139 - t2 * t151 - t23 * t56 + t49 * t45 + t5 * t88 + (-qJD(1) * t35 - t11) * t206, -t1 * t88 - t11 * t58 + t12 * t56 + t166 * t7 + t2 * t90 - t33 * t45 + t35 * t46 - t6 * t69, t1 * t151 - t10 * t166 + t6 * t139 + t23 * t58 - t49 * t46 - t5 * t90 + (qJD(1) * t33 + t12) * t206, t1 * t33 + t23 * t10 + t11 * t7 + t12 * t6 + t2 * t35 + t5 * t49; 0, 0, 0, 0, -t156, -t214, 0, 0, 0, 0, 0, t209 * t151, t209 * t153, 0, 0, 0, 0, 0, -t153 * t80 - t170 * t227 + (t151 * t116 - t173 * t223) * qJD(3), -t153 * t79 + t170 * t224 + (-t139 * t219 + (t118 - t179) * t151) * qJD(3), t157, -t61 * t153 + t237 * t18 - t19 * t246 + t74 * t207 - t3 * t87 + t4 * t89, -t153 * t45 + t237 * t139 + (t151 * t69 - t208 * t87) * qJD(3), t157, t153 * t46 - t246 * t139 + (-t151 * t166 + t208 * t89) * qJD(3), t1 * t89 - t237 * t11 - t12 * t246 - t5 * t153 + t2 * t87 + t23 * t207; 0, 0, 0, 0, 0, 0, t153 * t156 * t151, -t210 * t156, 0, 0, 0, -t153 * t214, t151 * t214, t118 * t227 + t233 (t79 - t229) * t152 + (-t80 - t228) * t150, t139 * t203 + (t139 * t221 + t153 * t171) * qJD(1), -t139 * t205 + (-t151 * t224 + t153 * t172) * qJD(1), -t139 * t208, -pkin(3) * t80 - t105 * t139 + (t139 * t223 - t151 * t172) * t134 + (-pkin(8) * t227 + t230) * qJD(4) + (t150 * t164 - t64 * t153) * qJD(1), -pkin(3) * t79 + t212 * t139 + t171 * t121 + (pkin(8) * t224 + t107 * t152) * qJD(4) + (t152 * t164 + t65 * t153) * qJD(1), -t3 * t114 - t166 * t236 - t231 * t18 - t235 * t19 + t248 * t4 + t27 * t69 + t169, t4 * t78 - t3 * t77 + t61 * t191 + (pkin(4) * t224 - t121) * t74 + (t60 - t27) * t19 + t236 * t18, -t5 * t248 + t66 * t45 - t238 * t69 + t235 * t23 + t239 * t139 + (-qJD(3) * t77 + t11) * t208, t1 * t248 + t231 * t11 + t2 * t114 - t235 * t12 - t166 * t239 + t24 * t69 + t169, -t5 * t114 - t66 * t46 + t238 * t166 - t231 * t23 + t240 * t139 + (qJD(3) * t78 - t12) * t208, t1 * t78 - t239 * t11 + t240 * t12 + t2 * t77 - t238 * t23 + t5 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t116, -t116 ^ 2 + t118 ^ 2, t79 + t229, t228 - t80, t177, -t107 * t118 - t134 * t188 + t65 * t139 + t158, t107 * t116 + t64 * t139 - t161, t19 * t166 - t242 + (-t148 * t45 - t149 * t46) * pkin(4) + (-t18 + t22) * t69, t18 * t21 - t19 * t22 + (-t118 * t74 + t148 * t4 + t149 * t3) * pkin(4), t21 * t139 - t244 - t28 * t69 + (pkin(5) - t143) * t177 + t3, t12 * t166 - t141 * t45 + t143 * t46 - t242 + (t11 - t213) * t69, t141 * t177 - t23 * t69 + t28 * t166 + (0.2e1 * qJD(6) - t22) * t139 + t195, t1 * t141 - t11 * t21 + t12 * t213 + t2 * t143 - t23 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t166 * t18 + t19 * t69 + t61, t139 * t166 + t45, t167, -t46 + t247, -t11 * t166 + t12 * t69 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166 * t69 - t177, t46 + t247, -t139 ^ 2 - t245, -t12 * t139 + t2 + t244;];
tauc_reg  = t13;
