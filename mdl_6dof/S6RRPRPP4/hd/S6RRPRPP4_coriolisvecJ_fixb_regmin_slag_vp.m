% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:40
% EndTime: 2019-03-09 10:01:48
% DurationCPUTime: 2.65s
% Computational Cost: add. (4827->371), mult. (10851->497), div. (0->0), fcn. (6574->6), ass. (0->200)
t269 = pkin(3) + pkin(7);
t173 = sin(qJ(2));
t232 = t173 * qJD(1);
t150 = qJD(4) + t232;
t174 = cos(qJ(4));
t172 = sin(qJ(4));
t233 = t172 * qJD(2);
t175 = cos(qJ(2));
t239 = qJD(1) * t175;
t118 = t174 * t239 + t233;
t230 = t174 * qJD(2);
t120 = -t172 * t239 + t230;
t170 = sin(pkin(9));
t171 = cos(pkin(9));
t69 = t171 * t118 + t170 * t120;
t279 = t150 * t69;
t219 = t174 * t232;
t250 = t170 * t172;
t235 = qJD(4) * t174;
t236 = qJD(4) * t172;
t273 = -t170 * t236 + t171 * t235;
t258 = -t171 * t219 + t232 * t250 - t273;
t194 = t170 * t174 + t171 * t172;
t257 = -t170 * t235 - t171 * t236 - t194 * t232;
t157 = pkin(7) * t232;
t278 = qJD(3) + t157;
t195 = -t170 * t118 + t171 * t120;
t277 = t195 ^ 2;
t226 = qJD(1) * qJD(2);
t276 = -0.2e1 * t226;
t229 = t174 * qJD(5);
t176 = -pkin(2) - pkin(8);
t243 = qJ(5) - t176;
t183 = t243 * t236 - t229;
t158 = pkin(7) * t239;
t126 = pkin(3) * t239 + t158;
t162 = pkin(2) * t232;
t196 = pkin(8) * t173 - qJ(3) * t175;
t99 = qJD(1) * t196 + t162;
t206 = t174 * t126 - t172 * t99;
t249 = t172 * t173;
t52 = (pkin(4) * t175 - qJ(5) * t249) * qJD(1) + t206;
t256 = t172 * t126 + t174 * t99;
t58 = qJ(5) * t219 + t256;
t205 = t174 * t243;
t93 = -qJD(4) * t205 - t172 * qJD(5);
t264 = (t183 - t52) * t171 + (t58 - t93) * t170;
t167 = qJD(2) * qJ(3);
t106 = t167 + t126;
t74 = t118 * pkin(4) + qJD(5) + t106;
t30 = t69 * pkin(5) - qJ(6) * t195 + t74;
t275 = t30 * t195;
t208 = -t173 * qJ(3) - pkin(1);
t114 = t176 * t175 + t208;
t137 = t269 * t173;
t241 = t174 * t114 + t172 * t137;
t234 = qJD(4) * t175;
t217 = t172 * t234;
t220 = pkin(4) * t174 + pkin(3);
t238 = qJD(2) * t173;
t274 = (-pkin(7) - t220) * t238 - pkin(4) * t217;
t272 = pkin(4) * t235 + t220 * t232 + t278;
t209 = t175 * t226;
t148 = pkin(7) * t209;
t112 = pkin(3) * t209 + t148;
t210 = t173 * t226;
t149 = pkin(2) * t210;
t231 = t173 * qJD(3);
t182 = qJD(2) * t196 - t231;
t77 = qJD(1) * t182 + t149;
t207 = t174 * t112 - t172 * t77;
t90 = t114 * qJD(1);
t228 = pkin(3) * t232 + t278;
t92 = t176 * qJD(2) + t228;
t56 = t172 * t92 + t174 * t90;
t180 = -qJD(4) * t56 + t207;
t80 = -qJD(4) * t118 + t172 * t210;
t14 = pkin(4) * t209 - t80 * qJ(5) - t120 * qJD(5) + t180;
t225 = -t172 * t112 - t174 * t77 - t92 * t235;
t188 = -t90 * t236 - t225;
t81 = t120 * qJD(4) - t174 * t210;
t17 = -t81 * qJ(5) - t118 * qJD(5) + t188;
t4 = t170 * t14 + t171 * t17;
t224 = qJ(6) * t209 + t4;
t1 = t150 * qJD(6) + t224;
t46 = -t118 * qJ(5) + t56;
t262 = t170 * t46;
t55 = -t172 * t90 + t174 * t92;
t45 = -t120 * qJ(5) + t55;
t41 = t150 * pkin(4) + t45;
t18 = t171 * t41 - t262;
t12 = -t150 * pkin(5) + qJD(6) - t18;
t43 = t171 * t46;
t19 = t170 * t41 + t43;
t13 = t150 * qJ(6) + t19;
t193 = -t171 * t174 + t250;
t3 = t171 * t14 - t170 * t17;
t2 = -pkin(5) * t209 - t3;
t271 = t1 * t194 - t257 * t12 - t258 * t13 + t193 * t2;
t270 = t257 * t18 - t258 * t19 - t193 * t3 + t194 * t4;
t21 = t170 * t45 + t43;
t268 = t21 * t195;
t28 = t170 * t52 + t171 * t58;
t23 = qJ(6) * t239 + t28;
t60 = t170 * t183 + t171 * t93;
t267 = -t23 + t60;
t266 = -pkin(5) * t239 + t264;
t237 = qJD(2) * t175;
t127 = t269 * t237;
t111 = t174 * t127;
t202 = qJ(5) * t175 - t114;
t161 = pkin(2) * t238;
t83 = t161 + t182;
t26 = pkin(4) * t237 + t111 + t202 * t235 + (-qJ(5) * t238 - qJD(4) * t137 + qJD(5) * t175 - t83) * t172;
t184 = -t114 * t236 + t172 * t127 + t137 * t235 + t174 * t83;
t218 = t173 * t230;
t31 = -t175 * t229 + (t217 + t218) * qJ(5) + t184;
t8 = t170 * t26 + t171 * t31;
t265 = t258 * pkin(5) + t257 * qJ(6) - t193 * qJD(6) - t272;
t122 = t174 * t137;
t62 = t173 * pkin(4) + t172 * t202 + t122;
t248 = t174 * t175;
t66 = -qJ(5) * t248 + t241;
t35 = t170 * t62 + t171 * t66;
t263 = qJD(2) * pkin(2);
t261 = t80 * t174;
t125 = t269 * t238;
t166 = qJD(2) * qJD(3);
t96 = -qJD(1) * t125 + t166;
t260 = t96 * t172;
t259 = t96 * t174;
t255 = t118 * t150;
t254 = t120 * t150;
t253 = t120 * t175;
t252 = t150 * t173;
t251 = t150 * t176;
t178 = qJD(1) ^ 2;
t247 = t175 * t178;
t177 = qJD(2) ^ 2;
t246 = t177 * t173;
t245 = t177 * t175;
t244 = t172 * pkin(4) + qJ(3);
t22 = t171 * t45 - t262;
t242 = qJD(6) - t22;
t138 = t269 * t175;
t168 = t173 ^ 2;
t169 = t175 ^ 2;
t240 = t168 - t169;
t131 = -t175 * pkin(2) + t208;
t107 = qJD(1) * t131;
t223 = t174 * t252;
t222 = t173 * t247;
t221 = pkin(4) * t248 + t138;
t216 = t150 * t235;
t215 = t174 * t234;
t212 = t173 * t233;
t48 = t170 * t80 + t171 * t81;
t204 = pkin(1) * t276;
t203 = qJD(3) - t263;
t49 = -t170 * t81 + t171 * t80;
t129 = t243 * t172;
t75 = -t170 * t129 + t171 * t205;
t76 = -t171 * t129 - t170 * t205;
t200 = -t76 * t48 + t75 * t49 - t60 * t69;
t197 = -t69 ^ 2 - t277;
t7 = -t170 * t31 + t171 * t26;
t34 = -t170 * t66 + t171 * t62;
t192 = -qJD(1) * t169 + t252;
t191 = -0.2e1 * qJD(2) * t107;
t190 = t150 * t172;
t185 = -qJ(3) * t237 - t231;
t101 = t161 + t185;
t88 = qJD(1) * t185 + t149;
t189 = pkin(7) * t177 + qJD(1) * t101 + t88;
t187 = t106 * t173 + t176 * t237;
t95 = t194 * t175;
t61 = t81 * pkin(4) + t96;
t181 = t193 * t49 - t194 * t48 - t195 * t257 + t258 * t69;
t128 = pkin(7) * t210 - t166;
t130 = t157 + t203;
t136 = -t158 - t167;
t179 = -t128 * t175 + (t130 * t175 + (t136 + t158) * t173) * qJD(2);
t9 = t48 * pkin(5) - t49 * qJ(6) - qJD(6) * t195 + t61;
t156 = -t171 * pkin(4) - pkin(5);
t152 = t170 * pkin(4) + qJ(6);
t140 = t174 * t209;
t123 = -qJ(3) * t239 + t162;
t94 = t193 * t175;
t91 = t107 * t232;
t67 = pkin(5) * t194 + qJ(6) * t193 + t244;
t64 = -t170 * t218 - t171 * t212 + t175 * t273;
t63 = qJD(4) * t95 - t193 * t238;
t50 = -t94 * pkin(5) + t95 * qJ(6) + t221;
t37 = t120 * pkin(4) + pkin(5) * t195 + qJ(6) * t69;
t33 = -t173 * pkin(5) - t34;
t32 = t173 * qJ(6) + t35;
t20 = -t63 * pkin(5) + t64 * qJ(6) + t95 * qJD(6) + t274;
t6 = -pkin(5) * t237 - t7;
t5 = qJ(6) * t237 + t173 * qJD(6) + t8;
t10 = [0, 0, 0, 0.2e1 * t173 * t209, t240 * t276, t245, -t246, 0, -pkin(7) * t245 + t173 * t204, pkin(7) * t246 + t175 * t204, t179, t173 * t191 + t175 * t189, -t173 * t189 + t175 * t191, pkin(7) * t179 + t107 * t101 + t88 * t131, -t80 * t172 * t175 + (t212 - t215) * t120 (-t118 * t172 + t120 * t174) * t238 + (t172 * t81 - t261 + (t118 * t174 + t120 * t172) * qJD(4)) * t175, -t150 * t215 + t80 * t173 + (t172 * t192 + t253) * qJD(2), t150 * t217 - t81 * t173 + (-t118 * t175 + t174 * t192) * qJD(2) (t150 + t232) * t237 (-t172 * t83 + t111) * t150 - t125 * t118 + t138 * t81 + (-t106 * t230 + t207) * t173 + (-t150 * t241 - t56 * t173) * qJD(4) + (-t106 * t236 + t259 + ((-t172 * t114 + t122) * qJD(1) + t55) * qJD(2)) * t175, -t184 * t150 - t125 * t120 + t138 * t80 + ((qJD(2) * t106 + qJD(4) * t90) * t172 + t225) * t173 + (-t106 * t235 - t260 + (-t241 * qJD(1) - t56) * qJD(2)) * t175, t18 * t64 + t19 * t63 - t195 * t7 + t3 * t95 - t34 * t49 - t35 * t48 + t4 * t94 - t8 * t69, t18 * t7 + t19 * t8 + t221 * t61 + t274 * t74 + t3 * t34 + t4 * t35, -t6 * t150 - t2 * t173 + t20 * t69 - t30 * t63 + t50 * t48 - t9 * t94 + (-qJD(1) * t33 - t12) * t237, t1 * t94 - t12 * t64 + t13 * t63 + t195 * t6 - t2 * t95 - t32 * t48 + t33 * t49 - t5 * t69, t1 * t173 + t5 * t150 - t20 * t195 + t30 * t64 - t50 * t49 + t9 * t95 + (qJD(1) * t32 + t13) * t237, t1 * t32 + t12 * t6 + t13 * t5 + t2 * t33 + t30 * t20 + t9 * t50; 0, 0, 0, -t222, t240 * t178, 0, 0, 0, t178 * pkin(1) * t173, pkin(1) * t247 ((-t136 - t167) * t173 + (-t130 + t203) * t175) * qJD(1), -t123 * t239 + t91, 0.2e1 * t166 + (t107 * t175 + t123 * t173) * qJD(1), -t128 * qJ(3) - t136 * qJD(3) - t107 * t123 + (-t136 * t173 + (-t130 - t263) * t175) * qJD(1) * pkin(7), -t120 * t190 + t261 (-t81 - t254) * t174 + (-t80 + t255) * t172, -t150 * t236 + t140 + (-t150 * t249 - t253) * qJD(1), -t216 + (-t223 + (t118 - t233) * t175) * qJD(1), -t150 * t239, qJ(3) * t81 + t260 - t206 * t150 + t228 * t118 + (t106 * t174 - t172 * t251) * qJD(4) + (t174 * t187 - t55 * t175) * qJD(1), qJ(3) * t80 + t259 + t256 * t150 + t228 * t120 + (-t106 * t172 - t174 * t251) * qJD(4) + (-t172 * t187 + t56 * t175) * qJD(1), -t195 * t264 + t28 * t69 + t200 - t270, t4 * t76 - t3 * t75 + t61 * t244 + t272 * t74 + (t60 - t28) * t19 + t264 * t18, t9 * t194 + t67 * t48 - t265 * t69 - t258 * t30 + t266 * t150 + (-qJD(2) * t75 + t12) * t239, -t195 * t266 + t23 * t69 + t200 - t271, t9 * t193 - t67 * t49 + t265 * t195 - t257 * t30 + t267 * t150 + (qJD(2) * t76 - t13) * t239, t1 * t76 - t266 * t12 + t267 * t13 + t2 * t75 - t265 * t30 + t9 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t168 * t178 - t177, t136 * qJD(2) + t148 + t91, 0, 0, 0, 0, 0, -qJD(2) * t118 - t150 * t190 + t140, -t216 - qJD(2) * t120 + (-t175 * t233 - t223) * qJD(1), t181, -t74 * qJD(2) + t270, t257 * t150 + (-t193 * t239 - t69) * qJD(2), t181, -t258 * t150 + (t194 * t239 + t195) * qJD(2), -t30 * qJD(2) + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120 * t118, -t118 ^ 2 + t120 ^ 2, t80 + t255, -t81 + t254, t209, -t106 * t120 + t56 * t150 + t180, t106 * t118 + t55 * t150 - t188, t19 * t195 - t268 + (-t170 * t48 - t171 * t49) * pkin(4) + (-t18 + t22) * t69, t18 * t21 - t19 * t22 + (-t120 * t74 + t170 * t4 + t171 * t3) * pkin(4), t21 * t150 - t275 - t37 * t69 + (pkin(5) - t156) * t209 + t3, t13 * t195 - t152 * t48 + t156 * t49 - t268 + (t12 - t242) * t69, t152 * t209 - t30 * t69 + t37 * t195 + (0.2e1 * qJD(6) - t22) * t150 + t224, t1 * t152 - t12 * t21 + t13 * t242 + t2 * t156 - t30 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t18 * t195 + t19 * t69 + t61, t150 * t195 + t48, t197, -t49 + t279, -t12 * t195 + t13 * t69 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195 * t69 - t209, t49 + t279, -t150 ^ 2 - t277, -t13 * t150 + t2 + t275;];
tauc_reg  = t10;
