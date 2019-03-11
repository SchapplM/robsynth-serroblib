% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:35
% EndTime: 2019-03-09 05:02:46
% DurationCPUTime: 4.95s
% Computational Cost: add. (5958->449), mult. (12894->613), div. (0->0), fcn. (9163->16), ass. (0->243)
t230 = cos(qJ(3));
t310 = qJD(1) * t230;
t198 = -qJD(4) + t310;
t193 = -qJD(6) + t198;
t225 = sin(qJ(4));
t229 = cos(qJ(4));
t298 = t229 * qJD(3);
t226 = sin(qJ(3));
t311 = qJD(1) * t226;
t173 = t225 * t311 - t298;
t306 = qJD(3) * t225;
t175 = t229 * t311 + t306;
t219 = sin(pkin(11));
t221 = cos(pkin(11));
t113 = t173 * t219 - t175 * t221;
t224 = sin(qJ(6));
t228 = cos(qJ(6));
t261 = -t173 * t221 - t175 * t219;
t320 = t228 * t261;
t59 = t113 * t224 + t320;
t334 = t193 * t59;
t299 = qJD(6) * t224;
t297 = qJD(1) * qJD(3);
t279 = t230 * t297;
t294 = t226 * qJDD(1);
t302 = qJD(4) * t226;
t359 = -qJD(1) * t302 + qJDD(3);
t108 = qJD(4) * t298 + (t279 + t294) * t229 + t359 * t225;
t109 = t225 * (qJD(3) * (qJD(4) + t310) + t294) - t359 * t229;
t49 = -t108 * t219 - t109 * t221;
t50 = t108 * t221 - t109 * t219;
t8 = qJD(6) * t320 + t113 * t299 + t224 * t49 + t228 * t50;
t374 = t8 + t334;
t358 = -t228 * t113 + t224 * t261;
t373 = t358 * t59;
t216 = qJ(1) + pkin(10);
t208 = sin(t216);
t209 = cos(t216);
t267 = g(1) * t209 + g(2) * t208;
t249 = t267 * t226;
t239 = -g(3) * t230 + t249;
t220 = sin(pkin(10));
t202 = pkin(1) * t220 + pkin(7);
t187 = t202 * qJD(1);
t304 = qJD(3) * t230;
t185 = t202 * qJDD(1);
t361 = -qJD(2) * qJD(3) - t185;
t289 = -t187 * t304 + t361 * t226;
t257 = -qJDD(3) * pkin(3) - t289;
t295 = qJDD(2) * t230;
t91 = t257 - t295;
t372 = qJD(4) * pkin(8) * t198 + t239 - t91;
t371 = t358 ^ 2 - t59 ^ 2;
t213 = qJ(4) + pkin(11) + qJ(6);
t200 = sin(t213);
t201 = cos(t213);
t325 = t208 * t230;
t120 = t200 * t209 - t201 * t325;
t324 = t209 * t230;
t122 = t200 * t208 + t201 * t324;
t212 = t230 * qJDD(1);
t164 = t226 * t297 + qJDD(4) - t212;
t222 = cos(pkin(10));
t204 = -pkin(1) * t222 - pkin(2);
t161 = -pkin(3) * t230 - pkin(8) * t226 + t204;
t268 = pkin(3) * t226 - pkin(8) * t230;
t178 = t268 * qJD(3);
t110 = qJD(1) * t178 + qJDD(1) * t161;
t101 = t229 * t110;
t147 = t226 * qJD(2) + t230 * t187;
t130 = qJD(3) * pkin(8) + t147;
t133 = t161 * qJD(1);
t82 = t130 * t229 + t133 * t225;
t146 = qJD(2) * t230 - t226 * t187;
t360 = qJD(3) * t146;
t90 = qJDD(3) * pkin(8) + qJDD(2) * t226 + t185 * t230 + t360;
t14 = pkin(4) * t164 - qJ(5) * t108 - qJD(4) * t82 - qJD(5) * t175 - t225 * t90 + t101;
t301 = qJD(4) * t229;
t291 = t225 * t110 + t133 * t301 + t229 * t90;
t303 = qJD(4) * t225;
t247 = -t130 * t303 + t291;
t16 = -qJ(5) * t109 - qJD(5) * t173 + t247;
t4 = t221 * t14 - t16 * t219;
t2 = pkin(5) * t164 - pkin(9) * t50 + t4;
t69 = -qJ(5) * t173 + t82;
t333 = t221 * t69;
t81 = -t130 * t225 + t229 * t133;
t68 = -qJ(5) * t175 + t81;
t61 = -pkin(4) * t198 + t68;
t27 = t219 * t61 + t333;
t357 = pkin(9) * t261;
t21 = t27 + t357;
t20 = t21 * t299;
t343 = g(3) * t226;
t129 = -qJD(3) * pkin(3) - t146;
t104 = pkin(4) * t173 + qJD(5) + t129;
t62 = -pkin(5) * t261 + t104;
t370 = g(1) * t122 - g(2) * t120 - t224 * t2 + t201 * t343 - t62 * t59 + t20;
t335 = t193 * t358;
t9 = qJD(6) * t358 + t224 * t50 - t228 * t49;
t368 = -t9 - t335;
t341 = qJ(5) + pkin(8);
t275 = qJD(4) * t341;
t284 = t225 * t310;
t300 = qJD(5) * t229;
t177 = t268 * qJD(1);
t317 = t229 * t146 + t225 * t177;
t367 = qJ(5) * t284 - t225 * t275 + t300 - t317;
t158 = t229 * t177;
t319 = t229 * t230;
t348 = pkin(4) * t226;
t255 = -qJ(5) * t319 + t348;
t366 = -qJD(1) * t255 - t229 * t275 - t158 + (-qJD(5) + t146) * t225;
t119 = t200 * t325 + t201 * t209;
t121 = -t200 * t324 + t201 * t208;
t5 = t219 * t14 + t221 * t16;
t3 = pkin(9) * t49 + t5;
t285 = t228 * t2 - t224 * t3;
t365 = -g(1) * t121 + g(2) * t119 + t200 * t343 - t62 * t358 + t285;
t364 = pkin(9) * t113;
t166 = t219 * t229 + t221 * t225;
t246 = t166 * t230;
t363 = qJD(1) * t246 - t166 * qJD(4);
t258 = t219 * t225 - t221 * t229;
t362 = t198 * t258;
t337 = -t367 * t219 + t366 * t221;
t336 = t366 * t219 + t367 * t221;
t318 = qJDD(2) - g(3);
t354 = t318 * t230;
t353 = -t147 + (-t284 + t303) * pkin(4);
t322 = t225 * t230;
t138 = t208 * t322 + t209 * t229;
t140 = t208 * t229 - t209 * t322;
t352 = -g(1) * t140 + g(2) * t138;
t281 = t230 * t298;
t244 = -t225 * t302 + t281;
t321 = t226 * t229;
t351 = t164 * t321 - t198 * t244;
t349 = pkin(4) * t219;
t160 = qJDD(6) + t164;
t144 = t166 * t226;
t145 = t258 * t226;
t263 = -t228 * t144 + t145 * t224;
t89 = -t144 * t224 - t145 * t228;
t92 = -qJD(3) * t246 + t258 * t302;
t282 = t225 * t304;
t93 = t166 * t302 + t219 * t282 - t221 * t281;
t29 = qJD(6) * t89 - t224 * t93 - t228 * t92;
t340 = t160 * t263 + t29 * t193;
t176 = t202 * t319;
t305 = qJD(3) * t226;
t326 = t202 * t225;
t315 = t229 * t178 + t305 * t326;
t40 = -t226 * t300 + t255 * qJD(3) + (-t176 + (qJ(5) * t226 - t161) * t225) * qJD(4) + t315;
t307 = qJD(3) * t202;
t316 = t161 * t301 + t225 * t178;
t44 = (-qJ(5) * qJD(4) - t307) * t321 + (-qJD(5) * t226 + (-qJ(5) * qJD(3) - qJD(4) * t202) * t230) * t225 + t316;
t19 = t219 * t40 + t221 * t44;
t262 = -t166 * t224 - t228 * t258;
t339 = qJD(6) * t262 + t363 * t224 + t362 * t228;
t112 = t166 * t228 - t224 * t258;
t338 = qJD(6) * t112 + t362 * t224 - t363 * t228;
t64 = t219 * t69;
t31 = t221 * t68 - t64;
t314 = t225 * t161 + t176;
t323 = t225 * t226;
t103 = -qJ(5) * t323 + t314;
t149 = t229 * t161;
t94 = -qJ(5) * t321 + t149 + (-pkin(4) - t326) * t230;
t43 = t221 * t103 + t219 * t94;
t26 = t221 * t61 - t64;
t17 = -pkin(5) * t198 + t26 + t364;
t332 = t228 * t17;
t331 = -t363 * pkin(5) + t353;
t330 = t108 * t225;
t329 = t173 * t198;
t328 = t175 * t198;
t327 = t198 * t229;
t189 = t341 * t225;
t190 = t341 * t229;
t118 = -t219 * t189 + t221 * t190;
t313 = pkin(4) * t323 + t226 * t202;
t217 = t226 ^ 2;
t312 = -t230 ^ 2 + t217;
t188 = qJD(1) * t204;
t308 = qJD(3) * t173;
t287 = pkin(4) * t282 + t202 * t304 + t301 * t348;
t207 = pkin(4) * t229 + pkin(3);
t286 = pkin(4) * t225 + pkin(7);
t283 = t198 * t306;
t280 = -t230 * t8 + t305 * t358;
t277 = qJD(6) * t17 + t3;
t18 = -t219 * t44 + t221 * t40;
t30 = -t219 * t68 - t333;
t42 = -t103 * t219 + t221 * t94;
t274 = -qJD(4) * t133 - t90;
t273 = -t108 * t230 + t175 * t305;
t272 = t198 * t202 + t130;
t117 = -t221 * t189 - t190 * t219;
t98 = -pkin(9) * t258 + t118;
t270 = pkin(5) * t311 + t362 * pkin(9) + qJD(6) * t98 - t337;
t97 = -pkin(9) * t166 + t117;
t269 = t363 * pkin(9) + qJD(6) * t97 + t336;
t266 = g(1) * t208 - g(2) * t209;
t227 = sin(qJ(1));
t231 = cos(qJ(1));
t265 = g(1) * t227 - g(2) * t231;
t28 = qJD(6) * t263 + t224 * t92 - t228 * t93;
t264 = -t160 * t89 + t193 * t28;
t7 = t224 * t17 + t228 * t21;
t260 = t207 * t230 + t226 * t341;
t256 = t265 * pkin(1);
t203 = pkin(4) * t221 + pkin(5);
t254 = t203 * t224 + t228 * t349;
t253 = t203 * t228 - t224 * t349;
t251 = t230 * t9 + t305 * t59;
t250 = pkin(2) + t260;
t248 = -t164 * t225 + t198 * t301;
t245 = -qJD(1) * t188 + t267;
t243 = -pkin(8) * t164 - t129 * t198;
t240 = 0.2e1 * t188 * qJD(3) - qJDD(3) * t202;
t237 = pkin(4) * t109 + qJDD(5) + t257;
t232 = qJD(3) ^ 2;
t235 = -0.2e1 * qJDD(1) * t204 - t202 * t232 + t266;
t48 = t237 - t295;
t233 = qJD(1) ^ 2;
t181 = qJDD(3) * t230 - t226 * t232;
t180 = qJDD(3) * t226 + t230 * t232;
t141 = t208 * t225 + t209 * t319;
t139 = -t208 * t319 + t209 * t225;
t134 = pkin(5) * t258 - t207;
t105 = pkin(5) * t144 + t313;
t76 = pkin(4) * t175 - pkin(5) * t113;
t63 = -pkin(5) * t92 + t287;
t35 = -pkin(9) * t144 + t43;
t34 = -pkin(5) * t230 + pkin(9) * t145 + t42;
t24 = -pkin(5) * t49 + t48;
t23 = t31 + t364;
t22 = t30 - t357;
t11 = pkin(9) * t92 + t19;
t10 = pkin(5) * t305 + pkin(9) * t93 + t18;
t6 = -t21 * t224 + t332;
t1 = [qJDD(1), t265, g(1) * t231 + g(2) * t227 (t220 ^ 2 + t222 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t256, qJDD(1) * t217 + 0.2e1 * t226 * t279, 0.2e1 * t226 * t212 - 0.2e1 * t312 * t297, t180, t181, 0, t226 * t240 + t230 * t235, -t226 * t235 + t230 * t240, t108 * t321 + t175 * t244 (-t173 * t229 - t175 * t225) * t304 + (-t330 - t109 * t229 + (t173 * t225 - t175 * t229) * qJD(4)) * t226, t273 + t351 (t109 + t283) * t230 + (t248 - t308) * t226, -t164 * t230 - t198 * t305 -(-t161 * t303 + t315) * t198 + t149 * t164 - g(1) * t139 - g(2) * t141 + (t173 * t307 - t101 + t272 * t301 + (qJD(3) * t129 - t164 * t202 - t274) * t225) * t230 + (qJD(3) * t81 + t109 * t202 + t129 * t301 + t91 * t225) * t226, t316 * t198 - t314 * t164 - g(1) * t138 - g(2) * t140 + (-t272 * t303 + (t129 * t229 + t175 * t202) * qJD(3) + t291) * t230 + (-t129 * t303 + t202 * t108 + t91 * t229 + (-t202 * t327 - t82) * qJD(3)) * t226, t113 * t18 - t144 * t5 + t145 * t4 + t19 * t261 + t226 * t266 + t26 * t93 + t27 * t92 - t42 * t50 + t43 * t49, t5 * t43 + t27 * t19 + t4 * t42 + t26 * t18 + t48 * t313 + t104 * t287 + t256 + (-g(1) * t286 - g(2) * t250) * t209 + (g(1) * t250 - g(2) * t286) * t208, t28 * t358 + t8 * t89, t263 * t8 + t28 * t59 - t29 * t358 - t89 * t9, -t264 + t280, t251 + t340, -t160 * t230 - t193 * t305 -(t10 * t228 - t11 * t224) * t193 + (-t224 * t35 + t228 * t34) * t160 - t285 * t230 + t6 * t305 - t63 * t59 + t105 * t9 - t24 * t263 + t62 * t29 - g(1) * t120 - g(2) * t122 + (-(-t224 * t34 - t228 * t35) * t193 + t7 * t230) * qJD(6), -t7 * t305 - g(1) * t119 - g(2) * t121 + t105 * t8 - t20 * t230 + t24 * t89 + t62 * t28 + t63 * t358 + ((-qJD(6) * t35 + t10) * t193 - t34 * t160 + t2 * t230) * t224 + ((qJD(6) * t34 + t11) * t193 - t35 * t160 + t277 * t230) * t228; 0, 0, 0, t318, 0, 0, 0, 0, 0, t181, -t180, 0, 0, 0, 0, 0 (-t109 + t283) * t230 + (t248 + t308) * t226, t273 - t351, t113 * t92 + t144 * t50 - t145 * t49 - t261 * t93, t104 * t305 - t144 * t4 - t145 * t5 - t230 * t48 + t26 * t92 - t27 * t93 - g(3), 0, 0, 0, 0, 0, -t251 + t340, t264 + t280; 0, 0, 0, 0, -t226 * t233 * t230, t312 * t233, t294, t212, qJDD(3), qJD(3) * t147 + t245 * t226 + t289 + t354, t360 + (qJD(3) * t187 - t318) * t226 + (t245 + t361) * t230, -t175 * t327 + t330 (t108 + t329) * t229 + (-t109 + t328) * t225 (-t175 * t226 + t198 * t319) * qJD(1) - t248, t198 * t303 + t164 * t229 + (t173 * t226 - t198 * t322) * qJD(1), t198 * t311, -t81 * t311 - pkin(3) * t109 - t147 * t173 + t158 * t198 + (-t146 * t198 + t243) * t225 + t372 * t229, -pkin(3) * t108 - t147 * t175 - t317 * t198 - t372 * t225 + t243 * t229 + t82 * t311, t337 * t113 - t117 * t50 + t118 * t49 - t166 * t4 - t267 * t230 - t258 * t5 - t362 * t26 + t336 * t261 + t363 * t27 - t343, t5 * t118 + t4 * t117 - t48 * t207 - g(3) * t260 + t336 * t27 + t337 * t26 + t353 * t104 + t267 * (t207 * t226 - t230 * t341) t112 * t8 + t339 * t358, -t112 * t9 + t262 * t8 - t338 * t358 + t339 * t59, t112 * t160 - t193 * t339 - t311 * t358, t160 * t262 + t193 * t338 - t311 * t59, t193 * t311 (-t224 * t98 + t228 * t97) * t160 + t134 * t9 - t24 * t262 - t6 * t311 + t338 * t62 - t331 * t59 + (t224 * t269 + t228 * t270) * t193 + t239 * t201 -(t224 * t97 + t228 * t98) * t160 + t134 * t8 + t24 * t112 + t7 * t311 + t339 * t62 + t331 * t358 + (-t224 * t270 + t228 * t269) * t193 - t239 * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175 * t173, -t173 ^ 2 + t175 ^ 2, t108 - t329, -t109 - t328, t164, -t130 * t301 - t129 * t175 - t198 * t82 + t101 + (t274 + t343) * t225 + t352, g(1) * t141 - g(2) * t139 + g(3) * t321 + t129 * t173 - t198 * t81 - t247 (t219 * t49 - t221 * t50) * pkin(4) + (-t31 + t26) * t261 + (-t27 - t30) * t113, -t26 * t30 - t27 * t31 + (g(3) * t323 - t104 * t175 + t5 * t219 + t4 * t221 + t352) * pkin(4), -t373, t371, t374, t368, t160, t253 * t160 + (t22 * t228 - t224 * t23) * t193 + t76 * t59 + (t193 * t254 - t7) * qJD(6) + t365, -t254 * t160 - t228 * t3 - (t22 * t224 + t228 * t23) * t193 - t76 * t358 + (t193 * t253 - t332) * qJD(6) + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113 ^ 2 - t261 ^ 2, -t113 * t26 - t261 * t27 + t237 - t249 - t354, 0, 0, 0, 0, 0, t9 - t335, t8 - t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t373, t371, t374, t368, t160 (-qJD(6) - t193) * t7 + t365, -t193 * t6 - t228 * t277 + t370;];
tau_reg  = t1;
