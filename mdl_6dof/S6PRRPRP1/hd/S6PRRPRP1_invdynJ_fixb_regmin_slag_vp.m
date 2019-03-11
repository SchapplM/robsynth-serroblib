% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [6x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:37
% EndTime: 2019-03-08 21:26:47
% DurationCPUTime: 4.15s
% Computational Cost: add. (4216->418), mult. (9954->578), div. (0->0), fcn. (7821->14), ass. (0->224)
t170 = sin(qJ(3));
t171 = sin(qJ(2));
t163 = sin(pkin(6));
t242 = qJD(1) * t163;
t224 = t171 * t242;
t271 = qJD(3) * pkin(3);
t296 = t170 * t271 - t224;
t173 = cos(qJ(3));
t280 = qJ(4) + pkin(8);
t214 = qJD(3) * t280;
t110 = t173 * qJD(4) - t170 * t214;
t111 = -t170 * qJD(4) - t173 * t214;
t161 = sin(pkin(11));
t164 = cos(pkin(11));
t54 = t110 * t164 + t111 * t161;
t250 = t164 * t173;
t125 = t161 * t170 - t250;
t174 = cos(qJ(2));
t223 = t174 * t242;
t88 = t125 * t223;
t276 = t54 + t88;
t126 = t161 * t173 + t164 * t170;
t117 = t126 * qJD(3);
t120 = t125 * qJD(3);
t295 = pkin(4) * t117 + pkin(9) * t120 + t296;
t169 = sin(qJ(5));
t172 = cos(qJ(5));
t294 = -t169 * t88 + t295 * t172;
t237 = qJD(5) * t172;
t152 = t173 * pkin(3) + pkin(2);
t292 = -pkin(9) * t126 - t152;
t69 = pkin(4) * t125 + t292;
t293 = t295 * t169 + t276 * t172 + t69 * t237;
t162 = sin(pkin(10));
t165 = cos(pkin(10));
t166 = cos(pkin(6));
t247 = t166 * t174;
t112 = t162 * t171 - t165 * t247;
t114 = t162 * t247 + t165 * t171;
t210 = g(1) * t114 + g(2) * t112;
t251 = t163 * t174;
t188 = -g(3) * t251 + t210;
t211 = t280 * qJD(2) + t224;
t241 = qJD(1) * t166;
t86 = t170 * t241 + t173 * t211;
t268 = t164 * t86;
t85 = -t170 * t211 + t173 * t241;
t80 = t85 + t271;
t38 = t161 * t80 + t268;
t33 = qJD(3) * pkin(9) + t38;
t108 = -t152 * qJD(2) + qJD(4) - t223;
t240 = qJD(2) * t170;
t116 = qJD(2) * t250 - t161 * t240;
t118 = t126 * qJD(2);
t45 = -pkin(4) * t116 - pkin(9) * t118 + t108;
t19 = t169 * t45 + t172 * t33;
t215 = qJDD(1) * t251;
t235 = qJD(1) * qJD(2);
t219 = t171 * t235;
t206 = t163 * t219 - t215;
t234 = qJD(2) * qJD(3);
t218 = t170 * t234;
t189 = pkin(3) * t218 + qJDD(4) + t206;
t137 = t161 * t218;
t217 = t173 * t234;
t197 = t164 * t217 - t137;
t231 = t173 * qJDD(2);
t232 = t170 * qJDD(2);
t205 = -t161 * t232 + t164 * t231;
t71 = -qJD(3) * t118 + t205;
t22 = -t71 * pkin(4) - t197 * pkin(9) + qJDD(2) * t292 + t189;
t21 = t172 * t22;
t233 = t166 * qJDD(1);
t141 = t173 * t233;
t101 = qJDD(2) * pkin(8) + (qJDD(1) * t171 + t174 * t235) * t163;
t183 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t241 + t101;
t201 = t211 * qJD(3);
t29 = qJDD(3) * pkin(3) - t170 * t183 - t173 * t201 + t141;
t30 = (-t201 + t233) * t170 + t183 * t173;
t11 = t161 * t29 + t164 * t30;
t9 = qJDD(3) * pkin(9) + t11;
t185 = -qJD(5) * t19 - t169 * t9 + t21;
t236 = t172 * qJD(3);
t238 = qJD(5) * t169;
t289 = qJDD(2) * t126 + t197;
t34 = -qJD(5) * t236 - t169 * qJDD(3) + t118 * t238 - t172 * t289;
t67 = qJD(2) * t117 + qJDD(5) - t205;
t94 = qJD(3) * t169 + t118 * t172;
t1 = pkin(5) * t67 + qJ(6) * t34 - qJD(6) * t94 + t185;
t109 = qJD(5) - t116;
t92 = t118 * t169 - t236;
t13 = -qJ(6) * t92 + t19;
t291 = t109 * t13 + t1;
t175 = qJD(3) ^ 2;
t290 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t175 + t163 * (-g(3) * t174 + t219) - t206 + t210;
t288 = t94 ^ 2;
t18 = -t169 * t33 + t172 * t45;
t12 = -qJ(6) * t94 + t18;
t6 = pkin(5) * t109 + t12;
t287 = t12 - t6;
t202 = qJ(6) * t120 - qJD(6) * t126;
t133 = t280 * t170;
t134 = t280 * t173;
t90 = -t133 * t161 + t134 * t164;
t81 = t172 * t90;
t286 = pkin(5) * t117 - t169 * t54 + t202 * t172 + (-t81 + (qJ(6) * t126 - t69) * t169) * qJD(5) + t294;
t220 = t126 * t237;
t285 = -qJ(6) * t220 + (-qJD(5) * t90 + t202) * t169 + t293;
t284 = pkin(3) * t164;
t283 = pkin(5) * t169;
t281 = g(3) * t163;
t178 = -t172 * qJDD(3) + t169 * t289;
t35 = qJD(5) * t94 + t178;
t279 = -t169 * t35 - t92 * t237;
t76 = t161 * t86;
t42 = t164 * t85 - t76;
t55 = pkin(3) * t240 + pkin(4) * t118 - pkin(9) * t116;
t278 = t169 * t55 + t172 * t42;
t277 = t161 * t110 - t164 * t111 - t126 * t223;
t275 = t169 * t69 + t81;
t148 = pkin(3) * t161 + pkin(9);
t245 = qJ(6) + t148;
t213 = qJD(5) * t245;
t260 = t116 * t169;
t274 = qJ(6) * t260 + t172 * qJD(6) - t169 * t213 - t278;
t51 = t172 * t55;
t273 = -pkin(5) * t118 - t51 + (qJ(6) * t116 - t213) * t172 + (-qJD(6) + t42) * t169;
t272 = qJD(2) * pkin(2);
t270 = t118 * t92;
t269 = t118 * t94;
t267 = t169 * t67;
t266 = t169 * t94;
t265 = t172 * t92;
t264 = t34 * t169;
t249 = t166 * t171;
t113 = t162 * t174 + t165 * t249;
t263 = -t112 * t152 + t113 * t280;
t115 = -t162 * t249 + t165 * t174;
t262 = -t114 * t152 + t115 * t280;
t259 = t120 * t172;
t258 = t126 * t169;
t257 = t126 * t172;
t158 = qJ(3) + pkin(11);
t154 = cos(t158);
t256 = t154 * t169;
t255 = t162 * t163;
t254 = t163 * t165;
t253 = t163 * t171;
t252 = t163 * t173;
t248 = t166 * t173;
t246 = t169 * t174;
t244 = qJDD(1) - g(3);
t159 = t170 ^ 2;
t243 = -t173 ^ 2 + t159;
t239 = qJD(2) * t171;
t228 = t165 * t252;
t227 = t163 * t246;
t226 = t170 * t253;
t225 = t172 * t251;
t151 = pkin(5) * t172 + pkin(4);
t222 = t163 * t239;
t221 = qJD(2) * t251;
t216 = t174 * t234;
t40 = t161 * t85 + t268;
t10 = -t161 * t30 + t164 * t29;
t37 = t164 * t80 - t76;
t89 = t164 * t133 + t134 * t161;
t212 = t109 * t172;
t209 = g(1) * t115 + g(2) * t113;
t208 = g(1) * t162 - g(2) * t165;
t194 = -t169 * t22 - t172 * t9 - t45 * t237 + t33 * t238;
t2 = -qJ(6) * t35 - qJD(6) * t92 - t194;
t207 = -t109 * t6 + t2;
t204 = t265 + t266;
t153 = sin(t158);
t167 = -qJ(6) - pkin(9);
t203 = t151 * t154 - t153 * t167;
t176 = qJD(2) ^ 2;
t200 = qJDD(2) * t174 - t171 * t176;
t199 = t172 * t67 + (-t238 + t260) * t109;
t32 = -qJD(3) * pkin(4) - t37;
t8 = -qJDD(3) * pkin(4) - t10;
t196 = -t172 * t34 - t94 * t238;
t121 = -t226 + t248;
t122 = t166 * t170 + t171 * t252;
t61 = t121 * t161 + t122 * t164;
t46 = -t169 * t61 - t225;
t195 = -t172 * t61 + t227;
t192 = -t120 * t169 + t220;
t191 = t109 * t32 - t148 * t67;
t103 = t153 * t253 - t166 * t154;
t72 = t113 * t153 + t154 * t254;
t74 = t115 * t153 - t154 * t255;
t190 = g(1) * t74 + g(2) * t72 + g(3) * t103;
t187 = -g(3) * t253 - t209;
t131 = -t223 - t272;
t186 = -qJD(2) * t131 - t101 + t209;
t3 = pkin(5) * t35 + qJDD(6) + t8;
t182 = qJD(5) * t109 * t148 - t190 + t8;
t181 = -pkin(8) * qJDD(3) + (t131 + t223 - t272) * qJD(3);
t70 = -qJDD(2) * t152 + t189;
t104 = t153 * t166 + t154 * t253;
t73 = t113 * t154 - t153 * t254;
t75 = t115 * t154 + t153 * t255;
t177 = -g(1) * (t114 * t172 - t169 * t75) - g(2) * (t112 * t172 - t169 * t73) - g(3) * (-t104 * t169 - t225);
t149 = -pkin(4) - t284;
t144 = pkin(3) * t248;
t132 = t162 * pkin(3) * t252;
t128 = t152 * t251;
t124 = t245 * t172;
t123 = t245 * t169;
t91 = t92 ^ 2;
t84 = -qJD(3) * t122 - t170 * t221;
t83 = qJD(3) * t121 + t173 * t221;
t64 = t172 * t69;
t60 = -t164 * t121 + t122 * t161;
t41 = t161 * t84 + t164 * t83;
t39 = t161 * t83 - t164 * t84;
t25 = -qJ(6) * t258 + t275;
t24 = pkin(5) * t92 + qJD(6) + t32;
t23 = pkin(5) * t125 - qJ(6) * t257 - t169 * t90 + t64;
t17 = qJD(5) * t195 - t169 * t41 + t172 * t222;
t16 = qJD(5) * t46 + t169 * t222 + t172 * t41;
t4 = [t244, 0, t200 * t163 (-qJDD(2) * t171 - t174 * t176) * t163, 0, 0, 0, 0, 0, qJD(3) * t84 + qJDD(3) * t121 + (-t170 * t216 + t173 * t200) * t163, -qJD(3) * t83 - qJDD(3) * t122 + (-t170 * t200 - t173 * t216) * t163, t41 * t116 + t39 * t118 + t289 * t60 + t61 * t71, -t10 * t60 + t11 * t61 - t37 * t39 + t38 * t41 - g(3) + (t108 * t239 - t174 * t70) * t163, 0, 0, 0, 0, 0, t109 * t17 + t35 * t60 + t39 * t92 + t46 * t67, -t109 * t16 + t195 * t67 - t34 * t60 + t39 * t94, -t16 * t92 - t17 * t94 + t195 * t35 + t34 * t46, t1 * t46 + t13 * t16 + t17 * t6 - t195 * t2 + t24 * t39 + t3 * t60 - g(3); 0, qJDD(2), t188 + t215, -t244 * t253 + t209, qJDD(2) * t159 + 0.2e1 * t170 * t217, 0.2e1 * t170 * t231 - 0.2e1 * t243 * t234, qJDD(3) * t170 + t173 * t175, qJDD(3) * t173 - t170 * t175, 0, t181 * t170 + t173 * t290, -t290 * t170 + t181 * t173, -t10 * t126 - t11 * t125 + t276 * t116 - t38 * t117 + t277 * t118 + t37 * t120 + t289 * t89 + t90 * t71 + t187, t11 * t90 - t10 * t89 - t70 * t152 - g(1) * t262 - g(2) * t263 - g(3) * (t253 * t280 + t128) + t276 * t38 - t277 * t37 + t296 * t108, t126 * t196 - t94 * t259, t204 * t120 + (t264 - t172 * t35 + (t169 * t92 - t172 * t94) * qJD(5)) * t126, t67 * t257 + t117 * t94 - t125 * t34 + (-t126 * t238 - t259) * t109, -t109 * t192 - t117 * t92 - t125 * t35 - t67 * t258, t109 * t117 + t125 * t67, t18 * t117 + t21 * t125 + t89 * t35 + t64 * t67 + t277 * t92 + t294 * t109 + (t188 * t154 + (-t109 * t90 - t125 * t33 + t126 * t32) * qJD(5)) * t172 + ((-qJD(5) * t69 - t54) * t109 - t90 * t67 + (-qJD(5) * t45 - t9) * t125 + t8 * t126 - t32 * t120 + t187) * t169, -t275 * t67 + t194 * t125 - t19 * t117 - t89 * t34 - t32 * t259 - g(1) * (t114 * t256 + t115 * t172) - g(2) * (t112 * t256 + t113 * t172) + t277 * t94 - (-t154 * t246 + t171 * t172) * t281 + (t8 * t172 - t32 * t238) * t126 + (t90 * t238 - t293) * t109, t23 * t34 - t25 * t35 - t286 * t94 - t285 * t92 - (-t13 * t169 - t172 * t6) * t120 + t188 * t153 + (-t1 * t172 - t169 * t2 + (-t13 * t172 + t169 * t6) * qJD(5)) * t126, t2 * t25 + t1 * t23 + t3 * (pkin(5) * t258 + t89) - g(1) * (-t114 * t203 + t115 * t283 + t262) - g(2) * (-t112 * t203 + t113 * t283 + t263) - g(3) * t128 + t286 * t6 + (pkin(5) * t192 + t277) * t24 + t285 * t13 - (t203 * t174 + (t280 + t283) * t171) * t281; 0, 0, 0, 0, -t170 * t176 * t173, t243 * t176, t232, t231, qJDD(3), -g(3) * t121 + t170 * t186 - t208 * t252 + t141, g(3) * t122 + (t163 * t208 - t233) * t170 + t186 * t173 (t38 - t40) * t118 + (-t42 + t37) * t116 + (t161 * t71 + (-t161 * t231 + t137 + (-t217 - t232) * t164) * t164) * pkin(3), -g(1) * t132 - g(3) * t144 + t37 * t40 - t38 * t42 + (g(2) * t228 + t10 * t164 + t11 * t161 + (-qJD(2) * t108 - t187) * t170) * pkin(3), t212 * t94 - t264, t116 * t204 + t196 + t279, t109 * t212 + t267 - t269, t199 + t270, -t109 * t118, -t51 * t109 - t18 * t118 + t149 * t35 - t40 * t92 + (t42 * t109 + t191) * t169 - t182 * t172, t278 * t109 + t19 * t118 - t149 * t34 + t182 * t169 + t191 * t172 - t40 * t94, -g(1) * t75 - g(2) * t73 - g(3) * t104 - t123 * t34 - t124 * t35 - t169 * t291 + t207 * t172 - t273 * t94 - t274 * t92, t2 * t124 - t1 * t123 + t3 * (-t151 - t284) - g(1) * (-pkin(3) * t115 * t170 - t151 * t74 - t167 * t75 + t132) - g(2) * (-t151 * t72 - t167 * t73 + (-t113 * t170 - t228) * pkin(3)) - g(3) * (-pkin(3) * t226 - t103 * t151 - t104 * t167 + t144) + t273 * t6 + (t109 * t283 - t40) * t24 + t274 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 ^ 2 - t118 ^ 2, -t116 * t38 + t118 * t37 - t188 + t70, 0, 0, 0, 0, 0, t199 - t270, -t109 ^ 2 * t172 - t267 - t269 (t265 - t266) * t116 - t196 + t279, -t118 * t24 + t207 * t169 + t172 * t291 - t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t92, -t91 + t288, t109 * t92 - t34, -t178 + (-qJD(5) + t109) * t94, t67, t19 * t109 - t32 * t94 + t177 + t185, t18 * t109 + t32 * t92 - g(1) * (-t114 * t169 - t172 * t75) - g(2) * (-t112 * t169 - t172 * t73) - g(3) * (-t104 * t172 + t227) + t194, pkin(5) * t34 + t287 * t92, -t287 * t13 + (-t24 * t94 + t1 + t177) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91 - t288, t13 * t92 + t6 * t94 - t190 + t3;];
tau_reg  = t4;
