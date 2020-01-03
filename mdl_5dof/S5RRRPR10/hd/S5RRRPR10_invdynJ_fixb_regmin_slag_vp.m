% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:54
% EndTime: 2019-12-31 21:30:09
% DurationCPUTime: 6.41s
% Computational Cost: add. (6545->462), mult. (16627->662), div. (0->0), fcn. (13288->14), ass. (0->247)
t220 = cos(qJ(2));
t296 = qJD(1) * qJD(2);
t274 = t220 * t296;
t216 = sin(qJ(2));
t294 = qJDD(1) * t216;
t370 = t274 + t294;
t210 = sin(pkin(5));
t304 = qJD(1) * t220;
t280 = t210 * t304;
t369 = qJD(3) - t280;
t215 = sin(qJ(3));
t219 = cos(qJ(3));
t259 = t215 * t280;
t213 = -qJ(4) - pkin(8);
t271 = qJD(3) * t213;
t306 = qJD(1) * t210;
t281 = t216 * t306;
t212 = cos(pkin(5));
t305 = qJD(1) * t212;
t291 = pkin(1) * t305;
t145 = -pkin(7) * t281 + t220 * t291;
t254 = pkin(2) * t216 - pkin(8) * t220;
t146 = t254 * t306;
t310 = t219 * t145 + t215 * t146;
t368 = qJ(4) * t259 + qJD(4) * t219 + t215 * t271 - t310;
t134 = t219 * t146;
t315 = t219 * t220;
t367 = t219 * t271 - t134 - (pkin(3) * t216 - qJ(4) * t315) * t306 + (-qJD(4) + t145) * t215;
t194 = qJD(2) + t305;
t260 = t215 * t281;
t128 = -t219 * t194 + t260;
t130 = t194 * t215 + t219 * t281;
t209 = sin(pkin(10));
t211 = cos(pkin(10));
t267 = -t211 * t128 - t130 * t209;
t356 = qJD(5) - t267;
t214 = sin(qJ(5));
t218 = cos(qJ(5));
t244 = -t128 * t209 + t211 * t130;
t61 = t214 * t244 - t218 * t369;
t366 = t356 * t61;
t245 = -t214 * t369 - t218 * t244;
t365 = t245 * t356;
t161 = t209 * t215 - t211 * t219;
t113 = t161 * t280;
t153 = t161 * qJD(3);
t364 = t113 - t153;
t162 = t209 * t219 + t211 * t215;
t309 = t369 * t162;
t363 = t370 * t210;
t269 = t218 * t356;
t295 = qJDD(1) * t212;
t193 = qJDD(2) + t295;
t300 = qJD(3) * t219;
t74 = -qJD(3) * t260 + t215 * t193 + t194 * t300 + t363 * t219;
t302 = qJD(2) * t220;
t277 = t215 * t302;
t301 = qJD(3) * t215;
t75 = -t219 * t193 + t194 * t301 + (qJD(1) * (t216 * t300 + t277) + t215 * t294) * t210;
t47 = -t209 * t74 - t211 * t75;
t46 = qJDD(5) - t47;
t334 = t214 * t46;
t362 = -t269 * t356 - t334;
t217 = sin(qJ(1));
t317 = t217 * t220;
t221 = cos(qJ(1));
t318 = t216 * t221;
t157 = t212 * t318 + t317;
t206 = qJ(3) + pkin(10);
t203 = sin(t206);
t204 = cos(t206);
t322 = t210 * t221;
t103 = -t157 * t204 + t203 * t322;
t313 = t220 * t221;
t319 = t216 * t217;
t156 = -t212 * t313 + t319;
t361 = t103 * t214 + t156 * t218;
t360 = t103 * t218 - t156 * t214;
t205 = t210 ^ 2;
t292 = 0.2e1 * t205;
t340 = -t368 * t209 + t367 * t211;
t338 = t367 * t209 + t368 * t211;
t323 = t210 * t220;
t349 = pkin(1) * t216;
t308 = pkin(7) * t323 + t212 * t349;
t141 = pkin(8) * t212 + t308;
t243 = -pkin(2) * t220 - pkin(8) * t216 - pkin(1);
t142 = t243 * t210;
t311 = t219 * t141 + t215 * t142;
t148 = pkin(7) * t280 + t216 * t291;
t358 = -t148 + (-t259 + t301) * pkin(3);
t357 = g(1) * t221 + g(2) * t217;
t355 = pkin(3) * t75 + qJDD(4);
t159 = -t212 * t319 + t313;
t324 = t210 * t219;
t109 = -t159 * t215 + t217 * t324;
t326 = t210 * t216;
t154 = -t212 * t219 + t215 * t326;
t314 = t219 * t221;
t329 = t157 * t215;
t354 = g(3) * t154 - g(2) * (-t210 * t314 - t329) - g(1) * t109;
t200 = pkin(3) * t209 + pkin(9);
t293 = qJDD(1) * t220;
t192 = t210 * t293;
t275 = t216 * t296;
t258 = t210 * t275;
t144 = qJDD(3) - t192 + t258;
t118 = pkin(8) * t194 + t148;
t122 = qJD(1) * t142;
t73 = t118 * t219 + t122 * t215;
t290 = pkin(1) * qJD(2) * t212;
t262 = qJD(1) * t290;
t288 = pkin(1) * t295;
t283 = -pkin(7) * t192 - t216 * t288 - t220 * t262;
t230 = -pkin(7) * t258 - t283;
t90 = pkin(8) * t193 + t230;
t241 = t254 * qJD(2);
t92 = (qJD(1) * t241 + qJDD(1) * t243) * t210;
t227 = -qJD(3) * t73 - t215 * t90 + t219 * t92;
t16 = pkin(3) * t144 - qJ(4) * t74 - qJD(4) * t130 + t227;
t240 = t118 * t301 - t122 * t300 - t215 * t92 - t219 * t90;
t18 = -qJ(4) * t75 - qJD(4) * t128 - t240;
t5 = t16 * t211 - t18 * t209;
t3 = -pkin(4) * t144 - t5;
t325 = t210 * t217;
t353 = (pkin(3) * t130 + pkin(4) * t244 - pkin(9) * t267 + qJD(5) * t200) * t356 + g(1) * (-t159 * t203 + t204 * t325) + g(2) * (-t157 * t203 - t204 * t322) + g(3) * (-t203 * t326 + t204 * t212) + t3;
t48 = -t209 * t75 + t211 * t74;
t20 = -qJD(5) * t245 - t218 * t144 + t214 * t48;
t261 = t363 * pkin(7) + t216 * t262 - t220 * t288;
t348 = pkin(2) * t193;
t91 = t261 - t348;
t50 = t91 + t355;
t10 = -pkin(4) * t47 - pkin(9) * t48 + t50;
t72 = -t118 * t215 + t219 * t122;
t57 = -qJ(4) * t130 + t72;
t53 = pkin(3) * t369 + t57;
t58 = -qJ(4) * t128 + t73;
t55 = t211 * t58;
t25 = t209 * t53 + t55;
t22 = pkin(9) * t369 + t25;
t117 = -pkin(2) * t194 - t145;
t82 = pkin(3) * t128 + qJD(4) + t117;
t28 = -pkin(4) * t267 - pkin(9) * t244 + t82;
t248 = t214 * t22 - t218 * t28;
t6 = t209 * t16 + t211 * t18;
t4 = pkin(9) * t144 + t6;
t1 = -t248 * qJD(5) + t214 * t10 + t218 * t4;
t202 = pkin(3) * t219 + pkin(2);
t100 = pkin(4) * t161 - pkin(9) * t162 - t202;
t179 = t213 * t219;
t276 = t213 * t215;
t116 = -t211 * t179 + t209 * t276;
t158 = t212 * t317 + t318;
t253 = g(1) * t158 + g(2) * t156;
t351 = t204 * t253 - (-t309 * pkin(4) + t364 * pkin(9) + qJD(5) * t116 - t358) * t356 + t100 * t46;
t343 = g(3) * t210;
t342 = t61 * t244;
t341 = t245 * t244;
t278 = t210 * t302;
t108 = -qJD(3) * t154 + t219 * t278;
t155 = t212 * t215 + t216 * t324;
t147 = t210 * t241;
t195 = pkin(7) * t326;
t321 = t212 * t220;
t149 = (pkin(1) * t321 - t195) * qJD(2);
t226 = -t311 * qJD(3) + t219 * t147 - t149 * t215;
t303 = qJD(2) * t216;
t279 = t210 * t303;
t32 = pkin(3) * t279 - qJ(4) * t108 - qJD(4) * t155 + t226;
t107 = qJD(3) * t155 + t210 * t277;
t239 = -t141 * t301 + t142 * t300 + t215 * t147 + t219 * t149;
t36 = -qJ(4) * t107 - qJD(4) * t154 + t239;
t14 = t209 * t32 + t211 * t36;
t266 = -t141 * t215 + t219 * t142;
t60 = -pkin(3) * t323 - qJ(4) * t155 + t266;
t68 = -qJ(4) * t154 + t311;
t35 = t209 * t60 + t211 * t68;
t339 = pkin(4) * t281 - t340;
t298 = qJD(5) * t218;
t299 = qJD(5) * t214;
t19 = t214 * t144 + t218 * t48 - t244 * t299 + t298 * t369;
t336 = t19 * t214;
t335 = t209 * t58;
t333 = t128 * t369;
t332 = t130 * t369;
t328 = t162 * t218;
t327 = t205 * qJD(1) ^ 2;
t320 = t214 * t220;
t316 = t218 * t220;
t150 = pkin(7) * t278 + t216 * t290;
t207 = t216 ^ 2;
t307 = -t220 ^ 2 + t207;
t297 = qJD(2) - t194;
t287 = t220 * t327;
t286 = t210 * t320;
t285 = t210 * t316;
t265 = t157 * t219 - t215 * t322;
t264 = t194 + t305;
t263 = t193 + t295;
t256 = pkin(3) * t107 + t150;
t252 = -g(1) * t156 + g(2) * t158;
t251 = g(1) * t159 + g(2) * t157;
t13 = -t209 * t36 + t211 * t32;
t24 = t211 * t53 - t335;
t34 = -t209 * t68 + t211 * t60;
t9 = t214 * t28 + t218 * t22;
t31 = -pkin(9) * t323 + t35;
t140 = t195 + (-pkin(1) * t220 - pkin(2)) * t212;
t231 = pkin(3) * t154 + t140;
t98 = t211 * t154 + t155 * t209;
t99 = -t154 * t209 + t155 * t211;
t49 = pkin(4) * t98 - pkin(9) * t99 + t231;
t247 = t214 * t49 + t218 * t31;
t246 = -t214 * t31 + t218 * t49;
t242 = t218 * t46 + (t214 * t267 - t299) * t356;
t79 = t214 * t99 + t285;
t96 = -t113 * t214 - t218 * t281;
t238 = -t214 * t153 + t162 * t298 - t96;
t97 = -t113 * t218 + t214 * t281;
t237 = -t153 * t218 - t162 * t299 - t97;
t232 = g(3) * t323 - t253;
t21 = -pkin(4) * t369 - t24;
t27 = t211 * t57 - t335;
t229 = -t200 * t46 + (t21 + t27) * t356;
t228 = -pkin(8) * t144 + t117 * t369;
t2 = -qJD(5) * t9 + t218 * t10 - t214 * t4;
t225 = -t232 - t261;
t224 = -pkin(8) * qJD(3) * t369 - t232 - t91;
t223 = -t116 * t46 + t3 * t162 + (pkin(9) * t281 - qJD(5) * t100 - t338) * t356 - t251;
t201 = -pkin(3) * t211 - pkin(4);
t139 = t203 * t212 + t204 * t326;
t115 = -t179 * t209 - t211 * t276;
t110 = t159 * t219 + t215 * t325;
t105 = t159 * t204 + t203 * t325;
t80 = t218 * t99 - t286;
t77 = t105 * t218 + t158 * t214;
t76 = -t105 * t214 + t158 * t218;
t67 = -t107 * t209 + t108 * t211;
t66 = t211 * t107 + t108 * t209;
t38 = -qJD(5) * t286 + t214 * t67 - t218 * t279 + t99 * t298;
t37 = -qJD(5) * t79 + t214 * t279 + t218 * t67;
t30 = pkin(4) * t323 - t34;
t26 = t209 * t57 + t55;
t23 = pkin(4) * t66 - pkin(9) * t67 + t256;
t12 = pkin(9) * t279 + t14;
t11 = -pkin(4) * t279 - t13;
t7 = [qJDD(1), g(1) * t217 - g(2) * t221, t357, (qJDD(1) * t207 + 0.2e1 * t216 * t274) * t205, (t216 * t293 - t307 * t296) * t292, (t216 * t263 + t264 * t302) * t210, (t220 * t263 - t264 * t303) * t210, t193 * t212, -t150 * t194 - t195 * t193 - t261 * t212 + g(1) * t157 - g(2) * t159 + (t193 * t321 + (-t275 + t293) * t292) * pkin(1), -t370 * pkin(1) * t292 - t149 * t194 - t308 * t193 - t230 * t212 + t252, t108 * t130 + t155 * t74, -t107 * t130 - t108 * t128 - t154 * t74 - t155 * t75, t108 * t369 + t144 * t155 + (t130 * t303 - t220 * t74) * t210, -t107 * t369 - t144 * t154 + (-t128 * t303 + t220 * t75) * t210, (-t144 * t220 + t303 * t369) * t210, t226 * t369 + t266 * t144 + t150 * t128 + t140 * t75 + t91 * t154 + t117 * t107 + g(1) * t265 - g(2) * t110 + (-t220 * t227 + t303 * t72) * t210, -t239 * t369 - t311 * t144 + t150 * t130 + t140 * t74 + t91 * t155 + t117 * t108 - g(1) * t329 - g(2) * t109 + (-g(1) * t314 - t220 * t240 - t303 * t73) * t210, -t13 * t244 + t14 * t267 - t24 * t67 - t25 * t66 - t34 * t48 + t35 * t47 - t5 * t99 - t6 * t98 - t252, t6 * t35 + t25 * t14 + t5 * t34 + t24 * t13 + t50 * t231 + t82 * t256 - g(1) * (-pkin(1) * t217 + t156 * t213 - t157 * t202) - g(2) * (pkin(1) * t221 - t158 * t213 + t159 * t202) - t357 * t210 * (pkin(3) * t215 + pkin(7)), t19 * t80 - t245 * t37, -t19 * t79 - t20 * t80 + t245 * t38 - t37 * t61, t19 * t98 - t245 * t66 + t356 * t37 + t46 * t80, -t20 * t98 - t356 * t38 - t46 * t79 - t61 * t66, t356 * t66 + t46 * t98, (-qJD(5) * t247 - t12 * t214 + t218 * t23) * t356 + t246 * t46 + t2 * t98 - t248 * t66 + t11 * t61 + t30 * t20 + t3 * t79 + t21 * t38 - g(1) * t360 - g(2) * t77, -(qJD(5) * t246 + t12 * t218 + t214 * t23) * t356 - t247 * t46 - t1 * t98 - t9 * t66 - t11 * t245 + t30 * t19 + t3 * t80 + t21 * t37 + g(1) * t361 - g(2) * t76; 0, 0, 0, -t216 * t287, t307 * t327, (t297 * t304 + t294) * t210, -t297 * t281 + t192, t193, t148 * t194 + t327 * t349 + t225, pkin(1) * t287 + t145 * t194 + (pkin(7) * t296 + g(3)) * t326 + t251 + t283, t215 * t74 + t219 * t332, (t74 - t333) * t219 + (-t332 - t75) * t215, t369 * t300 + t144 * t215 + (-t130 * t216 - t315 * t369) * t306, -t369 * t301 + t144 * t219 + (t215 * t220 * t369 + t128 * t216) * t306, -t369 * t281, -t72 * t281 - pkin(2) * t75 - t148 * t128 - t134 * t369 + (t145 * t369 + t228) * t215 + t224 * t219, -pkin(2) * t74 - t148 * t130 - t215 * t224 + t219 * t228 + t281 * t73 + t310 * t369, -g(3) * t326 + t115 * t48 + t116 * t47 - t161 * t6 - t162 * t5 - t364 * t24 - t340 * t244 - t309 * t25 + t338 * t267 - t251, t6 * t116 - t5 * t115 - t50 * t202 - g(1) * (-t158 * t202 - t159 * t213) - g(2) * (-t156 * t202 - t157 * t213) + t358 * t82 + t338 * t25 + t340 * t24 - (t202 * t220 - t213 * t216) * t343, t19 * t328 - t237 * t245, t61 * t97 - t245 * t96 - (t214 * t245 - t218 * t61) * t153 + (-t336 - t20 * t218 + (t214 * t61 + t218 * t245) * qJD(5)) * t162, t161 * t19 + t237 * t356 - t245 * t309 + t328 * t46, -t161 * t20 - t162 * t334 - t238 * t356 - t309 * t61, t161 * t46 + t309 * t356, t115 * t20 + t2 * t161 - t309 * t248 + t339 * t61 + t351 * t218 + t223 * t214 - (t204 * t316 + t214 * t216) * t343 + t238 * t21, -t1 * t161 + t115 * t19 - t309 * t9 - t339 * t245 - t351 * t214 + t223 * t218 - (-t204 * t320 + t216 * t218) * t343 + t237 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t128, -t128 ^ 2 + t130 ^ 2, t74 + t333, -t75 + t332, t144, -t117 * t130 + t369 * t73 + t227 + t354, g(1) * t110 + g(2) * t265 + g(3) * t155 + t117 * t128 + t369 * t72 + t240, (t209 * t47 - t211 * t48) * pkin(3) + (t24 - t27) * t267 + (t25 - t26) * t244, t24 * t26 - t25 * t27 + (-t82 * t130 + t6 * t209 + t5 * t211 + t354) * pkin(3), -t245 * t269 + t336, (t19 - t366) * t218 + (-t20 + t365) * t214, t341 - t362, t242 + t342, -t356 * t244, t201 * t20 + t229 * t214 - t353 * t218 + t244 * t248 - t26 * t61, t201 * t19 + t353 * t214 + t229 * t218 + t9 * t244 + t245 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244 ^ 2 - t267 ^ 2, t24 * t244 - t25 * t267 - t225 - t348 + t355, 0, 0, 0, 0, 0, t242 - t342, t341 + t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245 * t61, t245 ^ 2 - t61 ^ 2, t19 + t366, -t20 - t365, t46, t9 * t356 + t21 * t245 - g(1) * t76 - g(2) * t361 - g(3) * (-t139 * t214 - t285) + t2, -t248 * t356 + t21 * t61 + g(1) * t77 - g(2) * t360 - g(3) * (-t139 * t218 + t286) - t1;];
tau_reg = t7;
