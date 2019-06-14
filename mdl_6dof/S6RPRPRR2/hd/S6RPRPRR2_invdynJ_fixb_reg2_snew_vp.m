% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:24:37
% EndTime: 2019-05-05 18:24:51
% DurationCPUTime: 6.48s
% Computational Cost: add. (35582->503), mult. (78008->730), div. (0->0), fcn. (54673->12), ass. (0->289)
t247 = sin(pkin(11));
t249 = cos(pkin(11));
t258 = cos(qJ(3));
t298 = qJD(1) * t258;
t255 = sin(qJ(3));
t299 = qJD(1) * t255;
t214 = t247 * t299 - t249 * t298;
t216 = t247 * t298 + t249 * t299;
t192 = t216 * t214;
t332 = qJDD(3) - t192;
t341 = t247 * t332;
t340 = t249 * t332;
t253 = sin(qJ(6));
t254 = sin(qJ(5));
t257 = cos(qJ(5));
t198 = -t257 * qJD(3) + t216 * t254;
t200 = qJD(3) * t254 + t216 * t257;
t256 = cos(qJ(6));
t167 = t256 * t198 + t200 * t253;
t169 = -t198 * t253 + t200 * t256;
t131 = t169 * t167;
t289 = qJD(1) * qJD(3);
t279 = t258 * t289;
t288 = t255 * qJDD(1);
t222 = t279 + t288;
t239 = t258 * qJDD(1);
t280 = t255 * t289;
t223 = t239 - t280;
t272 = t222 * t247 - t249 * t223;
t191 = qJDD(5) + t272;
t190 = qJDD(6) + t191;
t335 = -t131 + t190;
t339 = t253 * t335;
t173 = t200 * t198;
t333 = -t173 + t191;
t338 = t254 * t333;
t337 = t256 * t335;
t336 = t257 * t333;
t296 = qJD(3) * t216;
t174 = t272 + t296;
t194 = t222 * t249 + t223 * t247;
t269 = -qJDD(3) * t254 - t194 * t257;
t161 = -qJD(5) * t198 - t269;
t274 = -t257 * qJDD(3) + t194 * t254;
t267 = qJD(5) * t200 + t274;
t106 = -t167 * qJD(6) + t256 * t161 - t253 * t267;
t211 = qJD(5) + t214;
t207 = qJD(6) + t211;
t155 = t207 * t167;
t334 = -t155 + t106;
t183 = t211 * t198;
t136 = t161 + t183;
t227 = qJD(3) * pkin(3) - qJ(4) * t299;
t244 = t258 ^ 2;
t248 = sin(pkin(10));
t325 = sin(qJ(1));
t326 = cos(qJ(1));
t265 = t325 * g(1) - t326 * g(2);
t263 = qJDD(1) * pkin(1) + t265;
t266 = t326 * g(1) + t325 * g(2);
t330 = qJD(1) ^ 2;
t220 = -t330 * pkin(1) - t266;
t250 = cos(pkin(10));
t304 = t250 * t220;
t262 = qJDD(1) * pkin(7) + t248 * t263 + t304;
t301 = -g(3) + qJDD(2);
t278 = t255 * t301;
t149 = t258 * t262 + t278 - qJD(3) * t227 + t223 * qJ(4) + (-t258 * pkin(2) - t244 * pkin(3)) * t330;
t277 = t258 * t301;
t260 = -t255 * t262 + t277 - t222 * qJ(4) + qJDD(3) * pkin(3) + (qJD(3) * t258 * qJ(4) + (t258 * pkin(3) + pkin(2)) * t299) * qJD(1);
t109 = -0.2e1 * qJD(4) * t214 + t249 * t149 + t247 * t260;
t275 = t161 * t253 + t256 * t267;
t85 = (qJD(6) - t207) * t169 + t275;
t132 = (qJD(5) - t211) * t200 + t274;
t165 = t167 ^ 2;
t166 = t169 ^ 2;
t331 = t198 ^ 2;
t197 = t200 ^ 2;
t206 = t207 ^ 2;
t210 = t211 ^ 2;
t212 = t214 ^ 2;
t213 = t216 ^ 2;
t329 = 0.2e1 * qJD(4);
t217 = t250 * t263;
t273 = -t248 * t220 + t217;
t184 = -qJDD(1) * pkin(2) - t330 * pkin(7) - t273;
t242 = t244 * t330;
t157 = -t223 * pkin(3) - qJ(4) * t242 + t227 * t299 + qJDD(4) + t184;
t297 = qJD(3) * t214;
t271 = -t194 + t297;
t112 = t174 * pkin(4) + t271 * pkin(8) + t157;
t185 = pkin(4) * t214 - pkin(8) * t216;
t259 = qJD(3) ^ 2;
t96 = -pkin(4) * t259 + qJDD(3) * pkin(8) - t185 * t214 + t109;
t66 = -t257 * t112 + t254 * t96;
t47 = pkin(5) * t333 - pkin(9) * t136 - t66;
t180 = pkin(5) * t211 - pkin(9) * t200;
t67 = t254 * t112 + t257 * t96;
t49 = -t331 * pkin(5) - pkin(9) * t267 - t211 * t180 + t67;
t26 = t253 * t49 - t256 * t47;
t27 = t253 * t47 + t256 * t49;
t13 = t253 * t27 - t256 * t26;
t328 = pkin(5) * t13;
t88 = t155 + t106;
t57 = -t253 * t85 - t256 * t88;
t327 = pkin(5) * t57;
t324 = pkin(4) * t247;
t323 = t13 * t254;
t322 = t13 * t257;
t276 = t149 * t247 - t249 * t260;
t95 = -qJDD(3) * pkin(4) - t259 * pkin(8) + (t329 + t185) * t216 + t276;
t68 = pkin(5) * t267 - t331 * pkin(9) + t180 * t200 + t95;
t321 = t253 * t68;
t320 = t254 * t95;
t108 = t216 * t329 + t276;
t71 = -t108 * t249 + t109 * t247;
t319 = t255 * t71;
t318 = t256 * t68;
t317 = t257 * t95;
t122 = t131 + t190;
t316 = t122 * t253;
t315 = t122 * t256;
t147 = t173 + t191;
t314 = t147 * t254;
t313 = t147 * t257;
t312 = t157 * t247;
t311 = t157 * t249;
t188 = qJDD(3) + t192;
t310 = t188 * t247;
t309 = t188 * t249;
t308 = t207 * t253;
t307 = t207 * t256;
t306 = t211 * t254;
t305 = t211 * t257;
t235 = t258 * t330 * t255;
t228 = qJDD(3) + t235;
t303 = t255 * t228;
t229 = qJDD(3) - t235;
t302 = t258 * t229;
t295 = qJD(3) * t247;
t294 = qJD(3) * t249;
t291 = qJD(5) + t211;
t287 = t247 * t131;
t286 = t249 * t131;
t285 = t247 * t173;
t284 = t249 * t173;
t283 = -pkin(1) * t250 - pkin(2);
t282 = pkin(1) * t248 + pkin(7);
t281 = -pkin(4) * t249 - pkin(3);
t14 = t253 * t26 + t256 * t27;
t39 = t254 * t66 + t257 * t67;
t72 = t108 * t247 + t249 * t109;
t261 = -t330 * pkin(2) + t262;
t178 = t255 * t261 - t277;
t179 = t258 * t261 + t278;
t140 = t255 * t178 + t258 * t179;
t38 = t254 * t67 - t257 * t66;
t126 = -t206 - t165;
t76 = t126 * t253 + t337;
t268 = pkin(5) * t76 - t26;
t175 = -t272 + t296;
t141 = -t166 - t206;
t92 = t141 * t256 - t316;
t264 = pkin(5) * t92 - t27;
t243 = t255 ^ 2;
t240 = t243 * t330;
t234 = -t242 - t259;
t233 = -t240 - t259;
t226 = t240 + t242;
t225 = (t243 + t244) * qJDD(1);
t224 = t239 - 0.2e1 * t280;
t221 = 0.2e1 * t279 + t288;
t205 = -t213 - t259;
t204 = -t213 + t259;
t203 = t212 - t259;
t202 = -t233 * t255 - t302;
t201 = t234 * t258 - t303;
t186 = -t259 - t212;
t182 = -t197 + t210;
t181 = -t210 + t331;
t177 = t194 + t297;
t172 = -t212 - t213;
t171 = t197 - t331;
t164 = -t197 - t210;
t163 = -t205 * t247 - t309;
t162 = t205 * t249 - t310;
t159 = -t210 - t331;
t154 = t197 + t331;
t153 = -t166 + t206;
t152 = t165 - t206;
t151 = t186 * t249 - t341;
t150 = t186 * t247 + t340;
t142 = (-t198 * t257 + t200 * t254) * t211;
t139 = t175 * t249 + t177 * t247;
t138 = t175 * t247 - t177 * t249;
t137 = t291 * t198 + t269;
t135 = t161 - t183;
t133 = -t291 * t200 - t274;
t130 = t166 - t165;
t129 = t161 * t257 - t200 * t306;
t128 = t198 * t305 + t254 * t267;
t127 = -t162 * t255 + t163 * t258;
t125 = t181 * t257 - t314;
t124 = -t182 * t254 + t336;
t120 = -t164 * t254 - t313;
t119 = t164 * t257 - t314;
t118 = (-t167 * t256 + t169 * t253) * t207;
t117 = (-t167 * t253 - t169 * t256) * t207;
t116 = t159 * t257 - t338;
t115 = t159 * t254 + t336;
t114 = -t150 * t255 + t151 * t258;
t113 = -t165 - t166;
t105 = -qJD(6) * t169 - t275;
t104 = -t138 * t255 + t139 * t258;
t103 = -t132 * t257 + t136 * t254;
t102 = t133 * t257 - t135 * t254;
t101 = -t132 * t254 - t136 * t257;
t100 = t152 * t256 - t316;
t99 = -t153 * t253 + t337;
t98 = t152 * t253 + t315;
t97 = t153 * t256 + t339;
t93 = -t141 * t253 - t315;
t91 = t120 * t249 - t137 * t247;
t90 = t120 * t247 + t137 * t249;
t84 = (qJD(6) + t207) * t169 + t275;
t83 = t116 * t249 - t133 * t247;
t82 = t116 * t247 + t133 * t249;
t81 = t106 * t256 - t169 * t308;
t80 = t106 * t253 + t169 * t307;
t79 = -t105 * t253 + t167 * t307;
t78 = t105 * t256 + t167 * t308;
t77 = t126 * t256 - t339;
t75 = t103 * t249 - t154 * t247;
t74 = t103 * t247 + t154 * t249;
t73 = -t117 * t254 + t118 * t257;
t70 = -pkin(8) * t119 + t317;
t69 = -pkin(8) * t115 + t320;
t64 = t100 * t257 - t254 * t98;
t63 = -t254 * t97 + t257 * t99;
t62 = -t254 * t92 + t257 * t93;
t61 = t254 * t93 + t257 * t92;
t60 = -t255 * t90 + t258 * t91;
t59 = t253 * t88 - t256 * t85;
t58 = -t253 * t334 - t256 * t84;
t56 = -t253 * t84 + t256 * t334;
t55 = -t255 * t82 + t258 * t83;
t54 = -t254 * t80 + t257 * t81;
t53 = -t254 * t78 + t257 * t79;
t52 = -pkin(4) * t119 + t67;
t51 = -t254 * t76 + t257 * t77;
t50 = t254 * t77 + t257 * t76;
t48 = -pkin(4) * t115 + t66;
t44 = -pkin(9) * t92 + t318;
t43 = -pkin(9) * t76 + t321;
t42 = t258 * t72 - t319;
t41 = t247 * t334 + t249 * t62;
t40 = t247 * t62 - t249 * t334;
t37 = t247 * t84 + t249 * t51;
t36 = t247 * t51 - t249 * t84;
t35 = -pkin(5) * t334 + pkin(9) * t93 + t321;
t34 = -pkin(5) * t84 + pkin(9) * t77 - t318;
t33 = -pkin(8) * t101 - t38;
t32 = t247 * t95 + t249 * t39;
t31 = t247 * t39 - t249 * t95;
t30 = -t254 * t57 + t257 * t59;
t29 = -t254 * t56 + t257 * t58;
t28 = t254 * t59 + t257 * t57;
t24 = t113 * t247 + t249 * t30;
t23 = -t113 * t249 + t247 * t30;
t22 = -t255 * t40 + t258 * t41;
t21 = -pkin(4) * t28 - t327;
t20 = -t255 * t36 + t258 * t37;
t19 = -pkin(4) * t61 - t264;
t18 = -pkin(4) * t50 - t268;
t17 = -pkin(8) * t61 - t254 * t35 + t257 * t44;
t15 = -pkin(8) * t50 - t254 * t34 + t257 * t43;
t12 = -t255 * t23 + t258 * t24;
t11 = -pkin(5) * t68 + pkin(9) * t14;
t10 = -pkin(9) * t57 - t13;
t9 = -pkin(5) * t113 + pkin(9) * t59 + t14;
t8 = t14 * t257 - t323;
t7 = t14 * t254 + t322;
t6 = t247 * t68 + t249 * t8;
t5 = t247 * t8 - t249 * t68;
t4 = -pkin(8) * t28 + t10 * t257 - t254 * t9;
t3 = -pkin(4) * t7 - t328;
t2 = -pkin(8) * t7 - pkin(9) * t322 - t11 * t254;
t1 = -t255 * t5 + t258 * t6;
t16 = [0, 0, 0, 0, 0, qJDD(1), t265, t266, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t250 - t330 * t248) + t273, -t304 - t248 * t265 + (-0.2e1 * qJDD(1) * t248 - t330 * t250) * pkin(1), 0, pkin(1) * (t248 ^ 2 * t263 + t250 * t217), (t222 + t279) * t255, t221 * t258 + t224 * t255, t303 + t258 * (-t240 + t259), (t223 - t280) * t258, t255 * (t242 - t259) + t302, 0, -t258 * t184 + pkin(2) * t224 + pkin(7) * t201 + pkin(1) * (t201 * t248 + t224 * t250), t255 * t184 - pkin(2) * t221 + pkin(7) * t202 + pkin(1) * (t202 * t248 - t221 * t250), pkin(2) * t226 + pkin(7) * t225 + pkin(1) * (t225 * t248 + t226 * t250) + t140, -pkin(2) * t184 + pkin(7) * t140 + pkin(1) * (t140 * t248 - t184 * t250), t255 * (t194 * t249 - t216 * t295) + t258 * (t194 * t247 + t216 * t294), t255 * (-t174 * t249 + t247 * t271) + t258 * (-t174 * t247 - t249 * t271), t255 * (-t204 * t247 + t340) + t258 * (t204 * t249 + t341), t255 * (t214 * t294 + t247 * t272) + t258 * (t214 * t295 - t249 * t272), t255 * (t203 * t249 - t310) + t258 * (t203 * t247 + t309), (t255 * (-t214 * t249 + t216 * t247) + t258 * (-t214 * t247 - t216 * t249)) * qJD(3), t255 * (-qJ(4) * t150 + t312) + t258 * (-pkin(3) * t174 + qJ(4) * t151 - t311) - pkin(2) * t174 + pkin(7) * t114 + pkin(1) * (t114 * t248 - t174 * t250), t255 * (-qJ(4) * t162 + t311) + t258 * (pkin(3) * t271 + qJ(4) * t163 + t312) + pkin(2) * t271 + pkin(7) * t127 + pkin(1) * (t127 * t248 + t250 * t271), t255 * (-qJ(4) * t138 - t71) + t258 * (-pkin(3) * t172 + qJ(4) * t139 + t72) - pkin(2) * t172 + pkin(7) * t104 + pkin(1) * (t104 * t248 - t172 * t250), -qJ(4) * t319 + t258 * (-pkin(3) * t157 + qJ(4) * t72) - pkin(2) * t157 + pkin(7) * t42 + pkin(1) * (-t157 * t250 + t248 * t42), t255 * (t129 * t249 + t285) + t258 * (t129 * t247 - t284), t255 * (t102 * t249 + t171 * t247) + t258 * (t102 * t247 - t171 * t249), t255 * (t124 * t249 + t136 * t247) + t258 * (t124 * t247 - t136 * t249), t255 * (t128 * t249 - t285) + t258 * (t128 * t247 + t284), t255 * (t125 * t249 - t132 * t247) + t258 * (t125 * t247 + t132 * t249), t255 * (t142 * t249 + t191 * t247) + t258 * (t142 * t247 - t191 * t249), t255 * (-qJ(4) * t82 - t247 * t48 + t249 * t69) + t258 * (-pkin(3) * t115 + qJ(4) * t83 + t247 * t69 + t249 * t48) - pkin(2) * t115 + pkin(7) * t55 + pkin(1) * (-t115 * t250 + t248 * t55), t255 * (-qJ(4) * t90 - t247 * t52 + t249 * t70) + t258 * (-pkin(3) * t119 + qJ(4) * t91 + t247 * t70 + t249 * t52) - pkin(2) * t119 + pkin(7) * t60 + pkin(1) * (-t119 * t250 + t248 * t60), t255 * (-qJ(4) * t74 + t249 * t33) + t258 * (qJ(4) * t75 + t247 * t33) + t282 * (-t255 * t74 + t258 * t75) + (t255 * t324 + t258 * t281 + t283) * t101, (t255 * (-pkin(8) * t249 + t324) + t258 * (-pkin(8) * t247 + t281) + t283) * t38 + (t282 + qJ(4)) * (-t255 * t31 + t258 * t32), t255 * (t249 * t54 + t287) + t258 * (t247 * t54 - t286), t255 * (t130 * t247 + t249 * t29) + t258 * (-t130 * t249 + t247 * t29), t255 * (t247 * t88 + t249 * t63) + t258 * (t247 * t63 - t249 * t88), t255 * (t249 * t53 - t287) + t258 * (t247 * t53 + t286), t255 * (-t247 * t85 + t249 * t64) + t258 * (t247 * t64 + t249 * t85), t255 * (t190 * t247 + t249 * t73) + t258 * (-t190 * t249 + t247 * t73), t255 * (-qJ(4) * t36 + t15 * t249 - t18 * t247) + t258 * (-pkin(3) * t50 + qJ(4) * t37 + t15 * t247 + t18 * t249) - pkin(2) * t50 + pkin(7) * t20 + pkin(1) * (t20 * t248 - t250 * t50), t255 * (-qJ(4) * t40 + t17 * t249 - t19 * t247) + t258 * (-pkin(3) * t61 + qJ(4) * t41 + t17 * t247 + t19 * t249) - pkin(2) * t61 + pkin(7) * t22 + pkin(1) * (t22 * t248 - t250 * t61), t255 * (-qJ(4) * t23 - t21 * t247 + t249 * t4) + t258 * (-pkin(3) * t28 + qJ(4) * t24 + t21 * t249 + t247 * t4) - pkin(2) * t28 + pkin(7) * t12 + pkin(1) * (t12 * t248 - t250 * t28), t255 * (-qJ(4) * t5 + t2 * t249 - t247 * t3) + t258 * (-pkin(3) * t7 + qJ(4) * t6 + t2 * t247 + t249 * t3) - pkin(2) * t7 + pkin(7) * t1 + pkin(1) * (t1 * t248 - t250 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, 0, 0, 0, 0, 0, 0, t228 * t258 + t234 * t255, -t229 * t255 + t233 * t258, 0, -t178 * t258 + t179 * t255, 0, 0, 0, 0, 0, 0, t150 * t258 + t151 * t255, t162 * t258 + t163 * t255, t138 * t258 + t139 * t255, t255 * t72 + t258 * t71, 0, 0, 0, 0, 0, 0, t255 * t83 + t258 * t82, t255 * t91 + t258 * t90, t255 * t75 + t258 * t74, t255 * t32 + t258 * t31, 0, 0, 0, 0, 0, 0, t255 * t37 + t258 * t36, t255 * t41 + t258 * t40, t23 * t258 + t24 * t255, t255 * t6 + t258 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t240 - t242, t288, t235, t239, qJDD(3), -t178, -t179, 0, 0, t192, t213 - t212, t177, -t192, t175, qJDD(3), pkin(3) * t150 - t108, pkin(3) * t162 - t109, pkin(3) * t138, pkin(3) * t71, t161 * t254 + t200 * t305, t133 * t254 + t135 * t257, t182 * t257 + t338, t198 * t306 - t257 * t267, t181 * t254 + t313, (-t198 * t254 - t200 * t257) * t211, pkin(3) * t82 + pkin(4) * t133 + pkin(8) * t116 - t317, pkin(3) * t90 + pkin(4) * t137 + pkin(8) * t120 + t320, pkin(3) * t74 + pkin(4) * t154 + pkin(8) * t103 + t39, pkin(3) * t31 - pkin(4) * t95 + pkin(8) * t39, t254 * t81 + t257 * t80, t254 * t58 + t257 * t56, t254 * t99 + t257 * t97, t254 * t79 + t257 * t78, t100 * t254 + t257 * t98, t117 * t257 + t118 * t254, pkin(3) * t36 - pkin(4) * t84 + pkin(8) * t51 + t254 * t43 + t257 * t34, pkin(3) * t40 - pkin(4) * t334 + pkin(8) * t62 + t254 * t44 + t257 * t35, pkin(3) * t23 - pkin(4) * t113 + pkin(8) * t30 + t10 * t254 + t257 * t9, pkin(3) * t5 - pkin(4) * t68 + pkin(8) * t8 - pkin(9) * t323 + t11 * t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, -t271, t172, t157, 0, 0, 0, 0, 0, 0, t115, t119, t101, t38, 0, 0, 0, 0, 0, 0, t50, t61, t28, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t171, t136, -t173, -t132, t191, -t66, -t67, 0, 0, t131, t130, t88, -t131, -t85, t190, t268, t264, t327, t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t130, t88, -t131, -t85, t190, -t26, -t27, 0, 0;];
tauJ_reg  = t16;
