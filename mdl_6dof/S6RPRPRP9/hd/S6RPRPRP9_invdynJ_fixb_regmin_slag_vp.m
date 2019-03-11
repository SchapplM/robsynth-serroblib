% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:31
% EndTime: 2019-03-09 03:29:41
% DurationCPUTime: 4.13s
% Computational Cost: add. (5035->489), mult. (9853->612), div. (0->0), fcn. (6582->10), ass. (0->228)
t173 = cos(qJ(3));
t251 = qJD(1) * qJD(3);
t230 = t173 * t251;
t170 = sin(qJ(3));
t246 = t170 * qJDD(1);
t192 = t230 + t246;
t120 = qJDD(5) + t192;
t163 = pkin(9) + qJ(5);
t155 = sin(t163);
t300 = g(3) * t170;
t174 = cos(qJ(1));
t162 = g(2) * t174;
t171 = sin(qJ(1));
t301 = g(1) * t171;
t316 = t162 - t301;
t306 = -t173 * t316 - t300;
t166 = sin(pkin(9));
t297 = pkin(8) + qJ(4);
t133 = t297 * t166;
t167 = cos(pkin(9));
t134 = t297 * t167;
t169 = sin(qJ(5));
t172 = cos(qJ(5));
t78 = -t133 * t169 + t134 * t172;
t319 = -t78 * t120 + t306 * t155;
t259 = qJD(1) * t173;
t232 = t166 * t259;
t257 = qJD(3) * t167;
t117 = -t232 + t257;
t239 = t167 * t259;
t258 = qJD(3) * t166;
t118 = t239 + t258;
t63 = -t172 * t117 + t118 * t169;
t318 = t63 ^ 2;
t260 = qJD(1) * t170;
t148 = qJD(5) + t260;
t317 = t148 * t63;
t240 = t166 * t260;
t271 = t172 * t167;
t252 = qJD(5) * t172;
t253 = qJD(5) * t169;
t309 = -t166 * t253 + t167 * t252;
t285 = t169 * t240 - t260 * t271 - t309;
t122 = t166 * t172 + t167 * t169;
t106 = t122 * qJD(1);
t108 = t122 * qJD(5);
t284 = t170 * t106 + t108;
t201 = t117 * t169 + t118 * t172;
t304 = t201 ^ 2;
t244 = pkin(8) * t167 * t170;
t211 = pkin(3) * t173 + qJ(4) * t170;
t124 = t211 * qJD(1);
t175 = -pkin(1) - pkin(7);
t142 = t175 * qJD(1) + qJD(2);
t278 = t166 * t173;
t74 = t167 * t124 - t142 * t278;
t42 = (pkin(4) * t173 + t244) * qJD(1) + t74;
t277 = t167 * t173;
t75 = t166 * t124 + t142 * t277;
t54 = pkin(8) * t240 + t75;
t315 = -qJD(4) * t122 - qJD(5) * t78 + t169 * t54 - t172 * t42;
t121 = t166 * t169 - t271;
t200 = -t133 * t172 - t134 * t169;
t314 = qJD(4) * t121 - qJD(5) * t200 + t169 * t42 + t172 * t54;
t313 = t148 * t201;
t312 = t316 * t167;
t272 = t171 * t173;
t225 = -g(1) * t272 + t300;
t140 = t175 * qJDD(1) + qJDD(2);
t227 = -t140 - t162;
t310 = t227 * t173 - t225;
t245 = t173 * qJDD(1);
t218 = t167 * qJDD(3) - t166 * t245;
t231 = t170 * t251;
t189 = t166 * t231 + t218;
t219 = t166 * qJDD(3) - t167 * t231;
t191 = t167 * t245 + t219;
t21 = -t117 * t252 + t118 * t253 - t169 * t189 - t172 * t191;
t256 = qJD(3) * t170;
t99 = t121 * t173;
t286 = qJD(3) * t99 + t170 * t108 + t106;
t98 = t121 * t170;
t308 = t120 * t98 + t148 * t286 + t173 * t21 + t201 * t256;
t22 = qJD(5) * t201 + t169 * t191 - t172 * t189;
t210 = pkin(3) * t170 - qJ(4) * t173;
t127 = qJ(2) + t210;
t115 = t167 * t127;
t229 = -t166 * t175 + pkin(4);
t61 = -pkin(8) * t277 + t170 * t229 + t115;
t273 = t170 * t175;
t83 = t166 * t127 + t167 * t273;
t73 = -pkin(8) * t278 + t83;
t203 = t169 * t61 + t172 * t73;
t100 = qJD(3) * t211 - t173 * qJD(4) + qJD(2);
t85 = t167 * t100;
t39 = t85 + (t173 * t229 + t244) * qJD(3);
t238 = t166 * t256;
t254 = qJD(3) * t175;
t235 = t173 * t254;
t70 = t166 * t100 + t167 * t235;
t51 = pkin(8) * t238 + t70;
t305 = -qJD(5) * t203 - t169 * t51 + t172 * t39;
t303 = pkin(4) * t166;
t302 = pkin(5) * t120;
t299 = g(3) * t173;
t298 = t201 * t63;
t126 = t170 * t142;
t88 = -pkin(4) * t240 + t126;
t296 = pkin(5) * t284 + qJ(6) * t285 - qJD(6) * t122 - t88;
t295 = -qJ(6) * t259 - t314;
t294 = pkin(5) * t259 - t315;
t52 = qJD(1) * t100 + qJDD(1) * t127;
t280 = t142 * t173;
t71 = qJDD(3) * qJ(4) + t170 * t140 + (qJD(4) + t280) * qJD(3);
t24 = t166 * t52 + t167 * t71;
t111 = t127 * qJD(1);
t112 = qJD(3) * qJ(4) + t126;
t56 = t166 * t111 + t167 * t112;
t55 = t167 * t111 - t112 * t166;
t31 = pkin(4) * t260 - pkin(8) * t118 + t55;
t35 = pkin(8) * t117 + t56;
t11 = t169 * t31 + t172 * t35;
t292 = t11 * t148;
t291 = t167 * t55;
t202 = -qJDD(3) * pkin(3) + t142 * t256 + qJDD(4);
t270 = t173 * t140;
t76 = t202 - t270;
t290 = t173 * t76;
t289 = t76 * t166;
t97 = t122 * t173;
t287 = t121 * qJD(1) - qJD(3) * t97 + qJD(5) * t98;
t283 = pkin(1) * qJDD(1);
t177 = qJD(1) ^ 2;
t282 = qJ(2) * t177;
t281 = qJ(6) * t120;
t156 = cos(t163);
t279 = t156 * t171;
t275 = t170 * t171;
t274 = t170 * t174;
t269 = t173 * t297;
t268 = t173 * t174;
t267 = t174 * t156;
t176 = qJD(3) ^ 2;
t266 = t175 * t176;
t10 = -t169 * t35 + t172 * t31;
t265 = qJD(6) - t10;
t264 = g(1) * t268 + g(2) * t272;
t263 = t174 * pkin(1) + t171 * qJ(2);
t164 = t170 ^ 2;
t165 = t173 ^ 2;
t262 = t164 - t165;
t261 = -t176 - t177;
t255 = qJD(3) * t173;
t250 = qJDD(1) * qJ(2);
t249 = qJDD(1) * t166;
t248 = qJDD(1) * t167;
t247 = qJDD(3) * t170;
t243 = g(1) * t279;
t242 = 0.2e1 * qJD(1) * qJD(2);
t241 = t174 * pkin(7) + t263;
t154 = pkin(4) * t167 + pkin(3);
t237 = t169 * t256;
t236 = t172 * t256;
t226 = -g(2) * t274 + t299;
t23 = -t166 * t71 + t167 * t52;
t14 = pkin(4) * t192 - pkin(8) * t191 + t23;
t17 = pkin(8) * t189 + t24;
t224 = -t172 * t14 + t169 * t17 + t35 * t252 + t31 * t253;
t223 = qJD(1) * t83 + t56;
t222 = qJD(3) * pkin(3) - qJD(4);
t82 = -t166 * t273 + t115;
t221 = -qJDD(1) * t82 - t23;
t116 = pkin(4) * t278 - t173 * t175;
t220 = qJDD(2) - t283;
t89 = t155 * t275 - t267;
t91 = t155 * t274 + t279;
t217 = g(1) * t91 + g(2) * t89;
t90 = t155 * t174 + t156 * t275;
t92 = -t155 * t171 + t170 * t267;
t216 = -g(1) * t92 - g(2) * t90;
t215 = g(2) * t173 * t267 + t120 * t200 + t156 * t300;
t214 = g(1) * t174 + g(2) * t171;
t209 = t282 + t301;
t208 = -t23 * t166 + t24 * t167;
t207 = -t166 * t55 + t167 * t56;
t204 = -t169 * t73 + t172 * t61;
t104 = -pkin(4) * t238 + t170 * t254;
t199 = t242 + 0.2e1 * t250;
t198 = pkin(5) * t156 + qJ(6) * t155 + t154;
t197 = t316 * t166;
t103 = -t222 - t280;
t196 = t169 * t14 + t172 * t17 + t31 * t252 - t253 * t35;
t195 = t169 * t39 + t172 * t51 + t61 * t252 - t253 * t73;
t194 = -g(1) * t275 - t226;
t193 = 0.2e1 * qJ(2) * t251 + qJDD(3) * t175;
t72 = -pkin(4) * t117 + t103;
t188 = g(1) * t89 - g(2) * t91 + t155 * t299 - t224;
t187 = t21 + t317;
t96 = t122 * t170;
t185 = -t96 * t120 + t148 * t287 - t173 * t22 + t63 * t256;
t184 = t199 - t214;
t18 = pkin(5) * t63 - qJ(6) * t201 + t72;
t182 = t18 * t201 + qJDD(6) - t188;
t181 = -g(1) * t90 + g(2) * t92 - t156 * t299 + t196;
t180 = -pkin(4) * t189 + t202;
t179 = t22 + t313;
t178 = t22 * pkin(5) + t21 * qJ(6) - qJD(6) * t201 + t180;
t159 = t174 * qJ(2);
t157 = qJDD(3) * t173;
t69 = -t166 * t235 + t85;
t59 = pkin(5) * t121 - qJ(6) * t122 - t154;
t50 = -t166 * t236 - t167 * t237 + t309 * t173;
t48 = t108 * t173 - t166 * t237 + t167 * t236;
t34 = pkin(5) * t97 + qJ(6) * t99 + t116;
t33 = t180 - t270;
t27 = pkin(5) * t201 + qJ(6) * t63;
t26 = -pkin(5) * t170 - t204;
t25 = qJ(6) * t170 + t203;
t9 = pkin(5) * t50 + qJ(6) * t48 + qJD(6) * t99 + t104;
t8 = qJ(6) * t148 + t11;
t7 = -pkin(5) * t148 + t265;
t6 = -t21 + t317;
t5 = -pkin(5) * t255 - t305;
t4 = qJ(6) * t255 + qJD(6) * t170 + t195;
t3 = t178 - t270;
t2 = qJDD(6) + t224 - t302;
t1 = qJD(6) * t148 + t196 + t281;
t12 = [qJDD(1), -t316, t214, qJDD(2) + t316 - 0.2e1 * t283, t184, -t220 * pkin(1) - g(1) * (-t171 * pkin(1) + t159) - g(2) * t263 + (t242 + t250) * qJ(2), qJDD(1) * t165 - 0.2e1 * t170 * t230, -0.2e1 * t170 * t245 + 0.2e1 * t251 * t262, -t170 * t176 + t157, -t173 * t176 - t247, 0, t193 * t173 + (t184 - t266) * t170, -t193 * t170 + (t199 - t266) * t173 - t264, -t197 + (t175 * t218 + t289 + (qJD(1) * t82 + t55) * qJD(3)) * t173 + (t69 * qJD(1) - t214 * t167 + (-t117 * t175 + (t175 * t259 - t103) * t166) * qJD(3) - t221) * t170, -t312 + (-t175 * t219 + (-t175 * t245 + t76) * t167 - t223 * qJD(3)) * t173 + (-t70 * qJD(1) - t83 * qJDD(1) - t24 + t214 * t166 + (-t103 * t167 + t118 * t175) * qJD(3)) * t170, t70 * t117 + t83 * t218 - t69 * t118 - t82 * t219 + (-t24 * t166 + t167 * t221) * t173 + (t166 * t223 + t291) * t256 + t264, t24 * t83 + t56 * t70 + t23 * t82 + t55 * t69 - g(1) * (pkin(3) * t274 - qJ(4) * t268 + t159) - g(2) * t241 + (t103 * t256 - t290) * t175 + (-g(1) * t175 - g(2) * t210) * t171, -t201 * t48 + t21 * t99, -t201 * t50 + t21 * t97 + t22 * t99 + t48 * t63, -t120 * t99 - t148 * t48 - t170 * t21 + t201 * t255, -t120 * t97 - t148 * t50 - t170 * t22 - t255 * t63, t120 * t170 + t148 * t255, t10 * t255 + t104 * t63 + t116 * t22 + t204 * t120 + t305 * t148 - t224 * t170 + t33 * t97 + t72 * t50 + t216, t104 * t201 - t11 * t255 - t116 * t21 - t120 * t203 - t195 * t148 - t196 * t170 - t33 * t99 - t72 * t48 + t217, -t120 * t26 - t148 * t5 - t170 * t2 + t18 * t50 + t22 * t34 - t255 * t7 + t3 * t97 + t63 * t9 + t216, -t1 * t97 - t2 * t99 + t201 * t5 - t21 * t26 - t22 * t25 - t4 * t63 - t48 * t7 - t50 * t8 + t264, t1 * t170 + t120 * t25 + t148 * t4 + t18 * t48 - t201 * t9 + t21 * t34 + t255 * t8 + t3 * t99 - t217, t1 * t25 + t8 * t4 + t3 * t34 + t18 * t9 + t2 * t26 + t7 * t5 - g(1) * (pkin(5) * t92 + qJ(6) * t91 + t154 * t274 - t268 * t297 + t159) - g(2) * (pkin(5) * t90 + qJ(6) * t89 + t174 * t303 + t241) + (-g(1) * (-t303 + t175) - g(2) * (t154 * t170 - t269)) * t171; 0, 0, 0, qJDD(1), -t177, t162 - t209 + t220, 0, 0, 0, 0, 0, t170 * t261 + t157, t173 * t261 - t247, -t164 * t249 + t173 * t218 + (-t167 * t177 + (-t117 - t232) * qJD(3)) * t170, -t164 * t248 - t173 * t191 + (t166 * t177 + (t118 - 0.2e1 * t239) * qJD(3)) * t170 (t167 * t218 + ((t231 + t245) * t167 + t219) * t166) * t170 + (qJD(1) * t167 + t166 * t255) * t118 + (-qJD(1) * t166 + t167 * t255) * t117, -t290 + t208 * t170 + (-t166 * t56 - t291) * qJD(1) + (t103 * t170 + t173 * t207) * qJD(3) + t316, 0, 0, 0, 0, 0, t185, t308, t185, -t201 * t287 - t21 * t96 + t22 * t98 + t286 * t63, -t308, -t1 * t98 - t173 * t3 + t18 * t256 + t2 * t96 - t286 * t8 - t287 * t7 + t316; 0, 0, 0, 0, 0, 0, t173 * t177 * t170, -t262 * t177, t245, -t246, qJDD(3) (-t227 - t282) * t173 + t225 (-t140 + t209) * t170 + t226, pkin(3) * t218 - t76 * t167 + (t312 + (-qJ(4) * t258 - t55) * qJD(1)) * t173 + (-qJ(4) * t249 + g(3) * t167 + t142 * t117 + (-t74 + (t103 + t222) * t166) * qJD(1)) * t170, -pkin(3) * t219 + t289 + (-pkin(3) * t248 - t197 + (-qJ(4) * t257 + t56) * qJD(1)) * t173 + (-qJ(4) * t248 - g(3) * t166 - t142 * t118 + (t75 + (-qJD(4) + t103) * t167) * qJD(1)) * t170, -t75 * t117 + t74 * t118 + (qJD(4) * t117 - t260 * t55 + t24) * t167 + (qJD(4) * t118 - t260 * t56 - t23) * t166 + (t166 * t191 + t167 * t189) * qJ(4) + t194, -t103 * t126 - t55 * t74 - t56 * t75 + t207 * qJD(4) + (-t76 - t306) * pkin(3) + (t170 * t316 + t208 - t299) * qJ(4), -t21 * t122 - t201 * t285, t121 * t21 - t122 * t22 - t201 * t284 + t285 * t63, t122 * t120 - t148 * t285 - t201 * t259, -t121 * t120 - t148 * t284 + t259 * t63, -t148 * t259, t33 * t121 - t154 * t22 - t88 * t63 + t284 * t72 + (-qJD(1) * t10 - t243) * t173 + t315 * t148 + t215, t11 * t259 + t33 * t122 + t314 * t148 + t154 * t21 - t88 * t201 - t285 * t72 + t319, t121 * t3 + t22 * t59 + t296 * t63 + t284 * t18 + (qJD(1) * t7 - t243) * t173 - t294 * t148 + t215, -t1 * t121 + t122 * t2 + t200 * t21 + t201 * t294 - t22 * t78 - t284 * t8 - t285 * t7 - t295 * t63 + t194, -t122 * t3 + t295 * t148 + t285 * t18 - t201 * t296 + t21 * t59 - t8 * t259 - t319, -g(3) * t269 + t1 * t78 + t296 * t18 + t198 * t300 - t2 * t200 + t294 * t7 + t295 * t8 + t3 * t59 + t316 * (t170 * t297 + t173 * t198); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t118 - t258) * t260 - t218, t117 * t260 + t191, -t117 ^ 2 - t118 ^ 2, -t117 * t56 + t118 * t55 + t202 + t310, 0, 0, 0, 0, 0, t179, -t187, t179, -t304 - t318, t187, -t201 * t7 + t63 * t8 + t178 + t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, t304 - t318, t6, -t22 + t313, t120, -t201 * t72 + t188 + t292, t10 * t148 + t63 * t72 - t181, -t27 * t63 - t182 + t292 + 0.2e1 * t302, pkin(5) * t21 - qJ(6) * t22 + (-t11 + t8) * t201 + (t7 - t265) * t63, 0.2e1 * t281 - t18 * t63 + t27 * t201 + (0.2e1 * qJD(6) - t10) * t148 + t181, t1 * qJ(6) - t2 * pkin(5) - t18 * t27 - t7 * t11 - g(1) * (-pkin(5) * t89 + qJ(6) * t90) - g(2) * (pkin(5) * t91 - qJ(6) * t92) + t265 * t8 - (-pkin(5) * t155 + qJ(6) * t156) * t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120 + t298, t6, -t148 ^ 2 - t304, -t148 * t8 + t182 - t302;];
tau_reg  = t12;
