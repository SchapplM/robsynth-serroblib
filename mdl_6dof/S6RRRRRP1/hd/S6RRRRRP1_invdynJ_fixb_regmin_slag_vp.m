% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:18
% EndTime: 2019-03-10 00:58:34
% DurationCPUTime: 6.86s
% Computational Cost: add. (12534->504), mult. (29245->633), div. (0->0), fcn. (22028->14), ass. (0->283)
t252 = cos(qJ(5));
t348 = qJD(5) * t252;
t249 = sin(qJ(3));
t250 = sin(qJ(2));
t355 = qJD(1) * t250;
t330 = t249 * t355;
t254 = cos(qJ(3));
t255 = cos(qJ(2));
t354 = qJD(1) * t255;
t332 = t254 * t354;
t164 = -t330 + t332;
t165 = -t249 * t354 - t254 * t355;
t248 = sin(qJ(4));
t253 = cos(qJ(4));
t426 = t253 * t164 + t165 * t248;
t435 = t426 * t252;
t440 = t348 - t435;
t247 = sin(qJ(5));
t288 = t164 * t248 - t253 * t165;
t241 = qJD(2) + qJD(3);
t324 = qJD(4) + t241;
t293 = t252 * t324;
t113 = t247 * t288 - t293;
t115 = t247 * t324 + t252 * t288;
t343 = t255 * qJDD(1);
t345 = qJD(1) * qJD(2);
t327 = t255 * t345;
t344 = t250 * qJDD(1);
t432 = t327 + t344;
t111 = qJD(3) * t332 - t241 * t330 + t249 * t343 + t254 * t432;
t177 = t249 * t255 + t250 * t254;
t141 = t241 * t177;
t297 = t249 * t344 - t254 * t343;
t112 = qJD(1) * t141 + t297;
t291 = t253 * t111 - t248 * t112;
t265 = qJD(4) * t426 + t291;
t240 = qJDD(2) + qJDD(3);
t305 = qJDD(4) + t240;
t349 = qJD(5) * t247;
t42 = -qJD(5) * t293 - t247 * t305 - t252 * t265 + t288 * t349;
t260 = t247 * t265 - t252 * t305;
t43 = qJD(5) * t115 + t260;
t346 = qJD(5) - t426;
t439 = t247 * t346;
t6 = -t113 * t440 - t115 * t439 - t247 * t43 - t42 * t252;
t319 = t248 * t111 + t253 * t112;
t56 = qJD(4) * t288 + t319;
t55 = qJDD(5) + t56;
t15 = t113 * t288 + t252 * t55 - t346 * t439;
t245 = qJ(2) + qJ(3);
t238 = qJ(4) + t245;
t224 = sin(t238);
t251 = sin(qJ(1));
t256 = cos(qJ(1));
t301 = g(1) * t256 + g(2) * t251;
t438 = t301 * t224;
t412 = pkin(7) + pkin(8);
t200 = t412 * t255;
t187 = qJD(1) * t200;
t170 = t254 * t187;
t199 = t412 * t250;
t185 = qJD(1) * t199;
t316 = t185 * t249 - t170;
t406 = pkin(9) * t164;
t118 = t316 - t406;
t160 = t165 * pkin(9);
t166 = t249 * t187;
t360 = -t254 * t185 - t166;
t119 = t160 + t360;
t229 = pkin(2) * t254 + pkin(3);
t350 = qJD(4) * t253;
t351 = qJD(4) * t248;
t368 = t248 * t249;
t428 = t118 * t248 + t119 * t253 - t229 * t350 - (-t249 * t351 + (t253 * t254 - t368) * qJD(3)) * pkin(2);
t391 = qJD(2) * pkin(2);
t172 = -t185 + t391;
t287 = -t172 * t249 - t170;
t109 = -t287 + t406;
t105 = t248 * t109;
t317 = t254 * t172 - t166;
t108 = t160 + t317;
t67 = t108 * t253 - t105;
t436 = -pkin(3) * t350 + t67;
t40 = t42 * t247;
t18 = t115 * t440 - t40;
t16 = -t115 * t288 + t247 * t55 + t346 * t440;
t100 = pkin(3) * t241 + t108;
t62 = t253 * t100 - t105;
t59 = -pkin(4) * t324 - t62;
t390 = t426 * t59;
t376 = t288 * t426;
t236 = cos(t245);
t225 = cos(t238);
t407 = pkin(5) * t252;
t227 = pkin(4) + t407;
t246 = -qJ(6) - pkin(10);
t314 = -t224 * t246 + t225 * t227;
t298 = pkin(3) * t236 + t314;
t367 = t249 * t253;
t433 = (t249 * t350 + (t248 * t254 + t367) * qJD(3)) * pkin(2) + t253 * t118;
t50 = t288 ^ 2 - t426 ^ 2;
t92 = pkin(4) * t288 - pkin(10) * t426;
t239 = t255 * pkin(2);
t399 = pkin(1) + t239;
t198 = t399 * qJD(1);
t144 = -t164 * pkin(3) - t198;
t217 = g(3) * t224;
t335 = t225 * t301 + t217;
t142 = qJDD(2) * pkin(2) - t412 * t432;
t328 = t250 * t345;
t143 = t412 * (-t328 + t343);
t269 = qJD(3) * t287 + t254 * t142 - t249 * t143;
t48 = pkin(3) * t240 - pkin(9) * t111 + t269;
t353 = qJD(3) * t249;
t414 = (qJD(3) * t172 + t143) * t254 + t249 * t142 - t187 * t353;
t54 = -t112 * pkin(9) + t414;
t415 = -(qJD(4) * t100 + t54) * t253 + t109 * t351 - t248 * t48;
t263 = -t144 * t426 + t335 + t415;
t237 = t252 * qJ(6);
t299 = pkin(5) * t288 - t426 * t237;
t44 = -t241 * t426 + t291;
t385 = -t119 * t248 + t229 * t351 + t433;
t380 = t346 * t288;
t300 = g(1) * t251 - g(2) * t256;
t427 = t300 * t224;
t359 = -t249 * t199 + t254 * t200;
t234 = t252 * qJD(6);
t378 = t426 * t247;
t425 = qJ(6) * t378 + t234;
t235 = sin(t245);
t286 = t224 * t227 + t225 * t246;
t424 = pkin(3) * t235 + t286;
t231 = pkin(2) * t355;
t411 = pkin(3) * t165;
t80 = -t411 + t92;
t75 = t231 + t80;
t423 = t247 * t75 + t252 * t428;
t422 = t247 * t80 + t252 * t436;
t364 = t252 * t256;
t370 = t247 * t251;
t151 = t225 * t370 + t364;
t366 = t251 * t252;
t369 = t247 * t256;
t153 = -t225 * t369 + t366;
t421 = -g(1) * t153 + g(2) * t151;
t401 = g(3) * t225;
t420 = t401 - t438;
t106 = t253 * t109;
t63 = t248 * t100 + t106;
t60 = pkin(10) * t324 + t63;
t72 = -pkin(4) * t426 - pkin(10) * t288 + t144;
t36 = -t247 * t60 + t252 * t72;
t57 = t59 * t349;
t418 = t252 * t438 - t36 * t288 + t57;
t321 = t100 * t351 + t109 * t350 + t248 * t54 - t253 * t48;
t11 = -pkin(4) * t305 + t321;
t37 = t247 * t72 + t252 * t60;
t400 = g(3) * t247;
t417 = t11 * t247 + t225 * t400 + t37 * t288 + t59 * t348;
t271 = -t144 * t288 - t321 - t420;
t45 = t241 * t288 - t319;
t413 = t115 ^ 2;
t409 = pkin(3) * t253;
t408 = pkin(5) * t247;
t23 = -qJ(6) * t115 + t36;
t20 = pkin(5) * t346 + t23;
t398 = -t23 + t20;
t397 = t247 * t92 + t252 * t62;
t315 = -t254 * t199 - t200 * t249;
t122 = -pkin(9) * t177 + t315;
t176 = t249 * t250 - t254 * t255;
t123 = -pkin(9) * t176 + t359;
t89 = t122 * t248 + t123 * t253;
t86 = t252 * t89;
t136 = t253 * t176 + t177 * t248;
t137 = -t176 * t248 + t177 * t253;
t148 = pkin(3) * t176 - t399;
t87 = pkin(4) * t136 - pkin(10) * t137 + t148;
t394 = t247 * t87 + t86;
t358 = pkin(2) * t367 + t248 * t229;
t157 = pkin(10) + t358;
t362 = -qJ(6) - t157;
t312 = qJD(5) * t362;
t393 = t247 * t312 - t423 + t425;
t71 = t252 * t75;
t392 = t252 * t312 - t299 - t71 + (-qJD(6) + t428) * t247;
t389 = t20 * t252;
t140 = t241 * t176;
t76 = -qJD(4) * t136 - t253 * t140 - t248 * t141;
t388 = t252 * t76;
t226 = pkin(3) * t248 + pkin(10);
t361 = -qJ(6) - t226;
t311 = qJD(5) * t361;
t384 = t247 * t311 - t422 + t425;
t79 = t252 * t80;
t383 = t252 * t311 - t299 - t79 + (-qJD(6) + t436) * t247;
t325 = qJD(5) * t246;
t382 = t234 - t397 + (qJ(6) * t426 + t325) * t247;
t91 = t252 * t92;
t381 = t252 * t325 - t91 + (-qJD(6) + t62) * t247 - t299;
t375 = t137 * t247;
t374 = t137 * t252;
t373 = t165 * t164;
t357 = (t349 - t378) * pkin(5);
t243 = t250 ^ 2;
t356 = -t255 ^ 2 + t243;
t352 = qJD(3) * t254;
t290 = t122 * t253 - t123 * t248;
t333 = qJD(2) * t412;
t186 = t250 * t333;
t188 = t255 * t333;
t276 = -t254 * t186 - t249 * t188 - t199 * t352 - t200 * t353;
t84 = -pkin(9) * t141 + t276;
t268 = -qJD(3) * t359 + t249 * t186 - t254 * t188;
t85 = t140 * pkin(9) + t268;
t27 = qJD(4) * t290 + t248 * t85 + t253 * t84;
t233 = t250 * t391;
t131 = pkin(3) * t141 + t233;
t77 = qJD(4) * t137 - t248 * t140 + t253 * t141;
t34 = pkin(4) * t77 - pkin(10) * t76 + t131;
t341 = t247 * t34 + t252 * t27 + t87 * t348;
t338 = qJD(5) * pkin(10) * t346;
t331 = t137 * t348;
t329 = -t11 - t401;
t326 = pkin(9) + t412 + t408;
t10 = pkin(10) * t305 - t415;
t323 = -qJD(5) * t72 - t10;
t313 = -pkin(2) * t368 + t229 * t253;
t161 = pkin(2) * t328 - qJDD(1) * t399;
t97 = t112 * pkin(3) + t161;
t14 = t56 * pkin(4) - pkin(10) * t265 + t97;
t285 = t252 * t10 + t247 * t14 + t72 * t348 - t349 * t60;
t3 = -qJ(6) * t43 - qJD(6) * t113 + t285;
t308 = t3 * t252 - t335;
t66 = t108 * t248 + t106;
t304 = pkin(3) * t351 - t66;
t156 = -pkin(4) - t313;
t13 = t252 * t14;
t303 = -t348 * t60 + t13;
t302 = -pkin(10) * t55 - t390;
t296 = -t157 * t55 - t390;
t295 = -t226 * t55 - t390;
t24 = -qJ(6) * t113 + t37;
t294 = -t24 * t247 - t389;
t292 = -qJ(6) * t76 - qJD(6) * t137;
t283 = -0.2e1 * pkin(1) * t345 - pkin(7) * qJDD(2);
t282 = t247 * t76 + t331;
t281 = -t137 * t349 + t388;
t280 = -t399 - t298;
t257 = qJD(2) ^ 2;
t273 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t257 + t300;
t258 = qJD(1) ^ 2;
t272 = pkin(1) * t258 - pkin(7) * qJDD(1) + t301;
t270 = t252 * t329 + t418;
t28 = qJD(4) * t89 + t248 * t84 - t253 * t85;
t264 = -t247 * t438 + t417;
t7 = t43 * pkin(5) + qJDD(6) + t11;
t1 = pkin(5) * t55 + qJ(6) * t42 - qJD(5) * t37 - qJD(6) * t115 - t247 * t10 + t13;
t262 = qJD(5) * t294 - t1 * t247 + t20 * t435 + t24 * t378 + t308;
t261 = g(3) * t235 + t198 * t164 + t236 * t301 - t414;
t259 = -g(3) * t236 - t198 * t165 + t235 * t301 + t269;
t228 = -pkin(4) - t409;
t197 = pkin(10) * t252 + t237;
t196 = t246 * t247;
t174 = t226 * t252 + t237;
t173 = t361 * t247;
t154 = t225 * t364 + t370;
t152 = -t225 * t366 + t369;
t147 = t231 - t411;
t146 = t157 * t252 + t237;
t145 = t362 * t247;
t117 = -t164 ^ 2 + t165 ^ 2;
t110 = t113 ^ 2;
t96 = -t297 + (-qJD(1) * t177 - t165) * t241;
t95 = -t164 * t241 + t111;
t83 = t252 * t87;
t46 = t113 * pkin(5) + qJD(6) + t59;
t38 = -qJ(6) * t375 + t394;
t33 = t252 * t34;
t31 = pkin(5) * t136 - t137 * t237 - t247 * t89 + t83;
t5 = -qJ(6) * t331 + (-qJD(5) * t89 + t292) * t247 + t341;
t4 = pkin(5) * t77 - t247 * t27 + t33 + t292 * t252 + (-t86 + (qJ(6) * t137 - t87) * t247) * qJD(5);
t2 = [qJDD(1), t300, t301, qJDD(1) * t243 + 0.2e1 * t250 * t327, 0.2e1 * t250 * t343 - 0.2e1 * t345 * t356, qJDD(2) * t250 + t255 * t257, qJDD(2) * t255 - t250 * t257, 0, t250 * t283 + t255 * t273, -t250 * t273 + t255 * t283, t111 * t177 + t140 * t165, -t111 * t176 - t112 * t177 - t140 * t164 + t141 * t165, -t140 * t241 + t177 * t240, -t141 * t241 - t176 * t240, 0, -t112 * t399 - t198 * t141 + t161 * t176 - t164 * t233 + t236 * t300 + t240 * t315 + t241 * t268, -t111 * t399 + t198 * t140 + t161 * t177 - t165 * t233 - t235 * t300 - t240 * t359 - t241 * t276, t137 * t265 + t288 * t76, -t136 * t265 - t137 * t56 - t288 * t77 + t426 * t76, t137 * t305 + t324 * t76, -t136 * t305 - t324 * t77, 0, -t131 * t426 + t97 * t136 + t144 * t77 + t148 * t56 + t225 * t300 - t28 * t324 + t290 * t305, t131 * t288 + t97 * t137 + t144 * t76 + t148 * t265 - t27 * t324 - t305 * t89 - t427, t115 * t281 - t374 * t42 (-t113 * t252 - t115 * t247) * t76 + (t40 - t252 * t43 + (t113 * t247 - t115 * t252) * qJD(5)) * t137, t115 * t77 - t42 * t136 + t281 * t346 + t374 * t55, -t113 * t77 - t43 * t136 - t282 * t346 - t375 * t55, t136 * t55 + t346 * t77 (-t348 * t89 + t33) * t346 + t83 * t55 + t303 * t136 + t36 * t77 + t28 * t113 - t290 * t43 + t59 * t331 - g(1) * t152 - g(2) * t154 + ((-qJD(5) * t87 - t27) * t346 - t89 * t55 + t323 * t136 + t11 * t137 + t59 * t76) * t247 -(-t349 * t89 + t341) * t346 - t394 * t55 - t285 * t136 - t37 * t77 + t28 * t115 + t290 * t42 + t59 * t388 - g(1) * t151 - g(2) * t153 + (t11 * t252 - t57) * t137, -t113 * t5 - t115 * t4 + t31 * t42 - t38 * t43 + t294 * t76 + t427 + (-t1 * t252 - t247 * t3 + (t20 * t247 - t24 * t252) * qJD(5)) * t137, t3 * t38 + t24 * t5 + t1 * t31 + t20 * t4 + t7 * (pkin(5) * t375 - t290) + t46 * (pkin(5) * t282 + t28) + (-g(1) * t326 + g(2) * t280) * t256 + (-g(1) * t280 - g(2) * t326) * t251; 0, 0, 0, -t250 * t258 * t255, t356 * t258, t344, t343, qJDD(2), -g(3) * t255 + t250 * t272, g(3) * t250 + t255 * t272, t373, t117, t95, t96, t240, -t316 * t241 + (t164 * t355 + t240 * t254 - t241 * t353) * pkin(2) + t259, t360 * t241 + (t165 * t355 - t240 * t249 - t241 * t352) * pkin(2) + t261, -t376, t50, t44, t45, t305, t147 * t426 + t305 * t313 - t324 * t385 + t271, -t147 * t288 - t305 * t358 + t324 * t428 + t263, t18, t6, t16, t15, -t380, t156 * t43 + t296 * t247 + t385 * t113 + (-t157 * t348 + t247 * t428 - t71) * t346 + t270, -t156 * t42 + t296 * t252 + t385 * t115 + (t157 * t349 + t423) * t346 + t264, -t113 * t393 - t115 * t392 + t145 * t42 - t146 * t43 + t262, t3 * t146 + t1 * t145 + t7 * (t156 - t407) - g(3) * (t239 + t298) + t393 * t24 + t392 * t20 + ((qJD(4) * t229 - t119) * t248 + t357 + t433) * t46 + t301 * (pkin(2) * t250 + t424); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t373, t117, t95, t96, t240, -t241 * t287 + t259, t241 * t317 + t261, -t376, t50, t44, t45, t305, t66 * t324 + (-t165 * t426 + t253 * t305 - t324 * t351) * pkin(3) + t271, t67 * t324 + (t165 * t288 - t248 * t305 - t324 * t350) * pkin(3) + t263, t18, t6, t16, t15, -t380, t228 * t43 + t295 * t247 + t304 * t113 + (-t226 * t348 + t247 * t436 - t79) * t346 + t270, -t228 * t42 + t295 * t252 + t304 * t115 + (t226 * t349 + t422) * t346 + t264, -t113 * t384 - t115 * t383 + t173 * t42 - t174 * t43 + t262, t3 * t174 + t1 * t173 + t7 * (-t227 - t409) - g(3) * t298 + (-t106 + (pkin(3) * qJD(4) - t108) * t248 + t357) * t46 + t384 * t24 + t383 * t20 + t301 * t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t376, t50, t44, t45, t305, t324 * t63 + t271, t324 * t62 + t263, t18, t6, t16, t15, -t380, -pkin(4) * t43 - t63 * t113 - t91 * t346 + (t346 * t62 + t302) * t247 + (t329 - t338) * t252 + t418, pkin(4) * t42 + t397 * t346 - t63 * t115 + t302 * t252 + (-t438 + t338) * t247 + t417, t196 * t42 - t197 * t43 - t346 * t389 - t381 * t115 - t382 * t113 + (-t24 * t346 - t1) * t247 + t308, t3 * t197 + t1 * t196 - t7 * t227 - g(3) * t314 + (t346 * t408 - t63) * t46 + t382 * t24 + t381 * t20 + t301 * t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t113, -t110 + t413, t113 * t346 - t42, -t260 + (-qJD(5) + t346) * t115, t55, -t115 * t59 + t346 * t37 + (t323 + t217) * t247 + t303 + t421, g(1) * t154 - g(2) * t152 + t113 * t59 + t217 * t252 + t346 * t36 - t285, pkin(5) * t42 - t113 * t398, t398 * t24 + (-t46 * t115 + t224 * t400 + t1 + t421) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 - t413, t24 * t113 + t20 * t115 + t420 + t7;];
tau_reg  = t2;
