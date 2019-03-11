% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:18
% EndTime: 2019-03-09 21:59:40
% DurationCPUTime: 9.60s
% Computational Cost: add. (15685->567), mult. (36819->735), div. (0->0), fcn. (28256->18), ass. (0->312)
t308 = cos(qJ(3));
t303 = sin(qJ(3));
t304 = sin(qJ(2));
t404 = qJD(1) * t304;
t384 = t303 * t404;
t309 = cos(qJ(2));
t403 = qJD(1) * t309;
t204 = t308 * t403 - t384;
t205 = -t303 * t403 - t308 * t404;
t302 = sin(qJ(4));
t307 = cos(qJ(4));
t169 = t307 * t204 + t205 * t302;
t165 = qJD(6) - t169;
t299 = sin(pkin(11));
t300 = cos(pkin(11));
t301 = sin(qJ(6));
t306 = cos(qJ(6));
t218 = t299 * t301 - t306 * t300;
t482 = t165 * t218;
t219 = t299 * t306 + t300 * t301;
t202 = t219 * qJD(6);
t481 = -t219 * t169 + t202;
t476 = t169 * t299;
t160 = pkin(5) * t476;
t391 = pkin(10) * t476;
t222 = t303 * t304 - t308 * t309;
t291 = t309 * pkin(2);
t443 = pkin(1) + t291;
t192 = pkin(3) * t222 - t443;
t245 = qJD(1) * t443;
t480 = qJD(6) - t165;
t479 = qJDD(1) * t443;
t392 = qJD(3) + qJD(4);
t286 = qJD(2) + t392;
t353 = t204 * t302 - t307 * t205;
t155 = t286 * t299 + t300 * t353;
t184 = -pkin(3) * t204 - t245;
t107 = -pkin(4) * t169 - qJ(5) * t353 + t184;
t198 = t205 * pkin(9);
t457 = pkin(7) + pkin(8);
t247 = t457 * t309;
t231 = qJD(1) * t247;
t206 = t303 * t231;
t246 = t457 * t304;
t229 = qJD(1) * t246;
t441 = qJD(2) * pkin(2);
t214 = -t229 + t441;
t374 = t308 * t214 - t206;
t150 = t198 + t374;
t294 = qJD(2) + qJD(3);
t140 = pkin(3) * t294 + t150;
t210 = t308 * t231;
t352 = -t214 * t303 - t210;
t449 = pkin(9) * t204;
t151 = -t352 + t449;
t146 = t307 * t151;
t95 = t302 * t140 + t146;
t92 = qJ(5) * t286 + t95;
t52 = t300 * t107 - t299 * t92;
t33 = -pkin(5) * t169 - pkin(10) * t155 + t52;
t153 = -t300 * t286 + t299 * t353;
t53 = t299 * t107 + t300 * t92;
t38 = -pkin(10) * t153 + t53;
t12 = t301 * t33 + t306 * t38;
t399 = qJD(4) * t307;
t400 = qJD(4) * t302;
t393 = t309 * qJDD(1);
t396 = qJD(1) * qJD(2);
t382 = t309 * t396;
t394 = t304 * qJDD(1);
t340 = -t382 - t394;
t395 = qJD(1) * qJD(3);
t472 = -t309 * t395 + t340;
t156 = -t294 * t384 + t303 * t393 - t308 * t472;
t292 = qJDD(2) + qJDD(3);
t180 = qJDD(2) * pkin(2) + t340 * t457;
t383 = t304 * t396;
t339 = -t383 + t393;
t183 = t457 * t339;
t323 = qJD(3) * t352 + t308 * t180 - t303 * t183;
t77 = pkin(3) * t292 - pkin(9) * t156 + t323;
t402 = qJD(3) * t303;
t200 = t231 * t402;
t368 = -qJD(3) * t214 - t183;
t81 = -t200 + (pkin(9) * t472 + t180) * t303 + ((-t304 * t395 + t339) * pkin(9) - t368) * t308;
t369 = -t140 * t400 - t151 * t399 - t302 * t81 + t307 * t77;
t283 = qJDD(4) + t292;
t464 = -pkin(4) * t283 + qJDD(5);
t28 = -t369 + t464;
t223 = t303 * t309 + t304 * t308;
t179 = t294 * t223;
t315 = qJD(1) * t179;
t313 = -t222 * qJDD(1) - t315;
t83 = t307 * t156 + t204 * t399 + t205 * t400 + t302 * t313;
t73 = -t300 * t283 + t299 * t83;
t15 = pkin(5) * t73 + t28;
t298 = qJ(2) + qJ(3);
t290 = qJ(4) + t298;
t277 = cos(t290);
t293 = pkin(11) + qJ(6);
t284 = sin(t293);
t276 = sin(t290);
t305 = sin(qJ(1));
t310 = cos(qJ(1));
t363 = g(1) * t310 + g(2) * t305;
t349 = t276 * t363;
t445 = g(3) * t284;
t144 = t302 * t151;
t94 = t140 * t307 - t144;
t91 = -pkin(4) * t286 + qJD(5) - t94;
t67 = pkin(5) * t153 + t91;
t316 = t12 * t353 + t15 * t219 + t277 * t445 - t284 * t349 - t482 * t67;
t285 = cos(t293);
t358 = t301 * t38 - t306 * t33;
t444 = g(3) * t285;
t421 = t276 * t310;
t422 = t276 * t305;
t473 = g(1) * t421 + g(2) * t422;
t321 = t15 * t218 - t277 * t444 + t285 * t473 + t358 * t353 + t481 * t67;
t355 = t153 * t301 - t155 * t306;
t84 = qJD(4) * t353 + t156 * t302 - t307 * t313;
t82 = qJDD(6) + t84;
t18 = -t165 * t482 + t219 * t82 + t353 * t355;
t99 = t306 * t153 + t155 * t301;
t19 = -t165 * t481 - t218 * t82 + t353 * t99;
t397 = qJD(6) * t306;
t398 = qJD(6) * t301;
t74 = t283 * t299 + t300 * t83;
t26 = -t153 * t397 - t155 * t398 - t301 * t73 + t306 * t74;
t8 = t26 * t219 + t355 * t482;
t27 = -qJD(6) * t355 + t301 * t74 + t306 * t73;
t1 = -t26 * t218 - t219 * t27 + t481 * t355 + t482 * t99;
t478 = t165 * t99;
t477 = t165 * t355;
t430 = t353 * t169;
t288 = cos(t298);
t407 = t277 * pkin(4) + t276 * qJ(5);
t388 = pkin(3) * t288 + t407;
t79 = -t169 ^ 2 + t353 ^ 2;
t471 = -g(3) * t277 + t473;
t61 = -t169 * t286 + t83;
t130 = pkin(4) * t353 - qJ(5) * t169;
t418 = t277 * t310;
t419 = t277 * t305;
t389 = g(1) * t418 + g(2) * t419 + g(3) * t276;
t458 = -(qJD(4) * t140 + t81) * t307 + t151 * t400 - t302 * t77;
t318 = -t184 * t169 + t389 + t458;
t468 = pkin(5) * t353;
t258 = pkin(3) * t399 + qJD(5);
t103 = t150 * t307 - t144;
t456 = pkin(3) * t205;
t115 = t130 - t456;
t56 = t300 * t103 + t299 * t115;
t467 = -t258 * t300 + t56;
t268 = t302 * t303 * pkin(2);
t279 = pkin(2) * t308 + pkin(3);
t401 = qJD(3) * t308;
t459 = t307 * pkin(2) * t401 - t268 * t392 + t279 * t399;
t172 = qJD(5) + t459;
t373 = t229 * t303 - t210;
t158 = t373 - t449;
t410 = -t308 * t229 - t206;
t159 = t198 + t410;
t109 = t158 * t302 + t159 * t307;
t281 = pkin(2) * t404;
t110 = t115 + t281;
t58 = t300 * t109 + t299 * t110;
t466 = -t172 * t300 + t58;
t434 = t165 * t353;
t362 = g(1) * t305 - g(2) * t310;
t465 = t362 * t276;
t101 = t150 * t302 + t146;
t365 = pkin(3) * t400 - t101;
t413 = t303 * t307;
t411 = t158 * t307 - t159 * t302 + t279 * t400 + (t303 * t399 + (t302 * t308 + t413) * qJD(3)) * pkin(2);
t409 = -t303 * t246 + t308 * t247;
t287 = sin(t298);
t364 = -pkin(3) * t287 - pkin(4) * t276;
t420 = t277 * t299;
t320 = t353 * t53 + g(3) * t420 + (t28 - t349) * t299;
t327 = -t353 * t52 + (-t28 + t471) * t300;
t332 = t369 + t471;
t325 = -t184 * t353 + t332;
t62 = t286 * t353 - t84;
t453 = pkin(3) * t307;
t450 = pkin(5) * t300;
t289 = t300 * pkin(10);
t24 = qJ(5) * t283 + qJD(5) * t286 - t458;
t267 = pkin(2) * t383;
t135 = pkin(3) * t315 + qJDD(1) * t192 + t267;
t31 = t84 * pkin(4) - t83 * qJ(5) - qJD(5) * t353 + t135;
t7 = t300 * t24 + t299 * t31;
t5 = t7 * t300;
t176 = t307 * t222 + t223 * t302;
t178 = t294 * t222;
t111 = -qJD(4) * t176 - t178 * t307 - t179 * t302;
t177 = -t222 * t302 + t223 * t307;
t112 = qJD(4) * t177 - t178 * t302 + t307 * t179;
t282 = t304 * t441;
t171 = pkin(3) * t179 + t282;
t43 = pkin(4) * t112 - qJ(5) * t111 - qJD(5) * t177 + t171;
t385 = qJD(2) * t457;
t230 = t304 * t385;
t232 = t309 * t385;
t338 = -t308 * t230 - t303 * t232 - t246 * t401 - t247 * t402;
t118 = -pkin(9) * t179 + t338;
t322 = -qJD(3) * t409 + t230 * t303 - t308 * t232;
t119 = pkin(9) * t178 + t322;
t372 = -t308 * t246 - t247 * t303;
t162 = -pkin(9) * t223 + t372;
t163 = -pkin(9) * t222 + t409;
t354 = t162 * t307 - t163 * t302;
t46 = qJD(4) * t354 + t118 * t307 + t119 * t302;
t17 = t299 * t43 + t300 * t46;
t440 = t169 * t91;
t437 = -t160 + t411;
t60 = t299 * t130 + t300 * t94;
t436 = qJ(5) * t300;
t435 = t111 * t299;
t431 = t169 * t300;
t428 = t177 * t299;
t427 = t177 * t300;
t406 = pkin(2) * t413 + t302 * t279;
t196 = qJ(5) + t406;
t426 = t196 * t300;
t425 = t205 * t204;
t271 = pkin(3) * t302 + qJ(5);
t423 = t271 * t300;
t417 = t284 * t305;
t416 = t284 * t310;
t415 = t285 * t305;
t414 = t285 * t310;
t412 = -qJD(5) + t91;
t125 = pkin(4) * t176 - qJ(5) * t177 + t192;
t127 = t162 * t302 + t163 * t307;
t66 = t299 * t125 + t300 * t127;
t296 = t304 ^ 2;
t405 = -t309 ^ 2 + t296;
t272 = -pkin(4) - t450;
t6 = -t24 * t299 + t300 * t31;
t3 = pkin(5) * t84 - pkin(10) * t74 + t6;
t4 = -pkin(10) * t73 + t7;
t386 = t306 * t3 - t301 * t4;
t380 = -pkin(2) * t304 + t364;
t16 = -t299 * t46 + t300 * t43;
t59 = t300 * t130 - t299 * t94;
t55 = -t103 * t299 + t300 * t115;
t57 = -t109 * t299 + t300 * t110;
t65 = t300 * t125 - t127 * t299;
t371 = t279 * t307 - t268;
t370 = t5 - t389;
t367 = -t160 + t365;
t366 = -pkin(10) * t431 + t468;
t197 = -pkin(4) - t371;
t361 = -t299 * t6 + t5;
t360 = t3 * t301 + t306 * t4;
t359 = -t299 * t52 + t300 * t53;
t50 = pkin(5) * t176 - pkin(10) * t427 + t65;
t54 = -pkin(10) * t428 + t66;
t357 = -t301 * t54 + t306 * t50;
t356 = t301 * t50 + t306 * t54;
t351 = t52 * t431 + t476 * t53 + t370;
t350 = t443 + t388;
t348 = t362 * t277;
t347 = -0.2e1 * pkin(1) * t396 - pkin(7) * qJDD(2);
t181 = (-pkin(10) - t196) * t299;
t346 = -qJD(6) * t181 - t391 + t466;
t182 = t289 + t426;
t345 = qJD(6) * t182 + t172 * t299 + t366 + t57;
t212 = (-pkin(10) - t271) * t299;
t344 = -qJD(6) * t212 - t391 + t467;
t213 = t289 + t423;
t343 = qJD(6) * t213 + t258 * t299 + t366 + t55;
t239 = (-pkin(10) - qJ(5)) * t299;
t342 = -qJD(5) * t300 - qJD(6) * t239 - t391 + t60;
t240 = t289 + t436;
t341 = qJD(5) * t299 + qJD(6) * t240 - t169 * t289 + t468 + t59;
t334 = t169 * t172 - t196 * t84 - t440;
t333 = t169 * t258 - t271 * t84 - t440;
t311 = qJD(2) ^ 2;
t331 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t311 + t362;
t312 = qJD(1) ^ 2;
t330 = pkin(1) * t312 - pkin(7) * qJDD(1) + t363;
t329 = t111 * t91 + t177 * t28 - t363;
t47 = qJD(4) * t127 + t118 * t302 - t119 * t307;
t317 = g(3) * t287 - t303 * t180 + t245 * t204 + t288 * t363 + t308 * t368 + t200;
t314 = -g(3) * t288 - t245 * t205 + t287 * t363 + t323;
t295 = -pkin(9) - t457;
t278 = -pkin(4) - t453;
t237 = qJ(5) * t418;
t236 = qJ(5) * t419;
t235 = t272 - t453;
t199 = t267 - t479;
t190 = t197 - t450;
t189 = t277 * t414 + t417;
t188 = -t277 * t416 + t415;
t187 = -t277 * t415 + t416;
t186 = t277 * t417 + t414;
t185 = t281 - t456;
t157 = -t204 ^ 2 + t205 ^ 2;
t134 = -t205 * t294 + t313;
t133 = -t204 * t294 + t156;
t132 = t218 * t177;
t131 = t219 * t177;
t90 = pkin(5) * t428 - t354;
t76 = t95 + t160;
t35 = t111 * t219 + t397 * t427 - t398 * t428;
t34 = -t111 * t218 - t177 * t202;
t32 = pkin(5) * t435 + t47;
t10 = -pkin(10) * t435 + t17;
t9 = pkin(5) * t112 - t111 * t289 + t16;
t2 = [qJDD(1), t362, t363, qJDD(1) * t296 + 0.2e1 * t304 * t382, 0.2e1 * t304 * t393 - 0.2e1 * t396 * t405, qJDD(2) * t304 + t309 * t311, qJDD(2) * t309 - t304 * t311, 0, t304 * t347 + t309 * t331, -t304 * t331 + t309 * t347, t156 * t223 + t178 * t205, -t156 * t222 - t178 * t204 + t205 * t179 + t223 * t313, -t178 * t294 + t223 * t292, -t179 * t294 - t222 * t292, 0, -t204 * t282 + t288 * t362 + t292 * t372 + t294 * t322 + (t199 - t479) * t222 - 0.2e1 * t245 * t179, -t156 * t443 + t245 * t178 + t199 * t223 - t205 * t282 - t287 * t362 - t292 * t409 - t294 * t338, t111 * t353 + t177 * t83, t111 * t169 - t112 * t353 - t176 * t83 - t177 * t84, t111 * t286 + t177 * t283, -t112 * t286 - t176 * t283, 0, t112 * t184 + t135 * t176 - t169 * t171 + t192 * t84 + t283 * t354 - t286 * t47 + t348, t111 * t184 - t127 * t283 + t135 * t177 + t171 * t353 + t192 * t83 - t286 * t46 - t465, t52 * t112 + t47 * t153 - t16 * t169 + t6 * t176 + t299 * t329 + t300 * t348 - t354 * t73 + t65 * t84, -t53 * t112 + t47 * t155 + t169 * t17 - t7 * t176 + t300 * t329 - t354 * t74 - t362 * t420 - t66 * t84, -t153 * t17 - t155 * t16 - t65 * t74 - t66 * t73 + t465 + (-t299 * t7 - t300 * t6) * t177 + (-t299 * t53 - t300 * t52) * t111, -t28 * t354 + t52 * t16 + t53 * t17 + t91 * t47 + t6 * t65 + t7 * t66 + (g(1) * t295 - g(2) * t350) * t310 + (g(1) * t350 + g(2) * t295) * t305, -t132 * t26 - t34 * t355, -t131 * t26 + t132 * t27 - t34 * t99 + t35 * t355, -t112 * t355 - t132 * t82 + t165 * t34 + t176 * t26, -t112 * t99 - t131 * t82 - t165 * t35 - t176 * t27, t112 * t165 + t176 * t82 (-t10 * t301 + t306 * t9) * t165 + t357 * t82 + t386 * t176 - t358 * t112 + t32 * t99 + t90 * t27 + t15 * t131 + t67 * t35 - g(1) * t187 - g(2) * t189 + (-t12 * t176 - t165 * t356) * qJD(6) -(t10 * t306 + t301 * t9) * t165 - t356 * t82 - t360 * t176 - t12 * t112 - t32 * t355 + t90 * t26 - t15 * t132 + t67 * t34 - g(1) * t186 - g(2) * t188 + (-t165 * t357 + t176 * t358) * qJD(6); 0, 0, 0, -t304 * t312 * t309, t405 * t312, t394, t393, qJDD(2), -g(3) * t309 + t304 * t330, g(3) * t304 + t309 * t330, t425, t157, t133, t134, t292, -t373 * t294 + (t204 * t404 + t292 * t308 - t294 * t402) * pkin(2) + t314, t410 * t294 + (t205 * t404 - t292 * t303 - t294 * t401) * pkin(2) + t317, -t430, t79, t61, t62, t283, t169 * t185 + t283 * t371 - t286 * t411 + t325, -t406 * t283 - t185 * t353 + (t109 - t459) * t286 + t318, t153 * t411 + t169 * t57 + t197 * t73 + t299 * t334 + t327, t155 * t411 - t169 * t58 + t197 * t74 + t300 * t334 + t320, -t73 * t426 + t155 * t57 + t466 * t153 + (t155 * t172 + t196 * t74 - t6) * t299 + t351, t28 * t197 - t53 * t58 - t52 * t57 - g(1) * (t310 * t380 + t237) - g(2) * (t305 * t380 + t236) - g(3) * (t291 + t388) + t411 * t91 + t361 * t196 + t359 * t172, t8, t1, t18, t19, -t434 (t181 * t306 - t182 * t301) * t82 + t190 * t27 + t437 * t99 + (t301 * t346 - t306 * t345) * t165 + t321 -(t181 * t301 + t182 * t306) * t82 + t190 * t26 + (t301 * t345 + t306 * t346) * t165 - t437 * t355 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425, t157, t133, t134, t292, -t294 * t352 + t314, t294 * t374 + t317, -t430, t79, t61, t62, t283, t101 * t286 + (-t169 * t205 + t283 * t307 - t286 * t400) * pkin(3) + t325, t103 * t286 + (t205 * t353 - t283 * t302 - t286 * t399) * pkin(3) + t318, t153 * t365 + t169 * t55 + t278 * t73 + t299 * t333 + t327, t155 * t365 - t169 * t56 + t278 * t74 + t300 * t333 + t320, -t73 * t423 + t155 * t55 + t467 * t153 + (t155 * t258 + t271 * t74 - t6) * t299 + t351, t28 * t278 - t53 * t56 - t52 * t55 - g(1) * (t310 * t364 + t237) - g(2) * (t305 * t364 + t236) - g(3) * t388 + t365 * t91 + t361 * t271 + t359 * t258, t8, t1, t18, t19, -t434 (t212 * t306 - t213 * t301) * t82 + t235 * t27 + t367 * t99 + (t301 * t344 - t306 * t343) * t165 + t321 -(t212 * t301 + t213 * t306) * t82 + t235 * t26 + (t301 * t343 + t306 * t344) * t165 - t367 * t355 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t430, t79, t61, t62, t283, t286 * t95 + t325, t286 * t94 + t318, -qJ(5) * t299 * t84 - pkin(4) * t73 - t153 * t95 - (t299 * t412 - t59) * t169 + t327, -t84 * t436 - pkin(4) * t74 - t155 * t95 - (t300 * t412 + t60) * t169 + t320, t153 * t60 + t155 * t59 + (-qJ(5) * t73 - qJD(5) * t153 + t169 * t52) * t300 + (qJ(5) * t74 + qJD(5) * t155 + t169 * t53 - t6) * t299 + t370, -t28 * pkin(4) - t53 * t60 - t52 * t59 - t91 * t95 - g(1) * (-pkin(4) * t421 + t237) - g(2) * (-pkin(4) * t422 + t236) - g(3) * t407 + t359 * qJD(5) + t361 * qJ(5), t8, t1, t18, t19, -t434 (t239 * t306 - t240 * t301) * t82 + t272 * t27 - t76 * t99 + (t301 * t342 - t306 * t341) * t165 + t321 -(t239 * t301 + t240 * t306) * t82 + t272 * t26 + t76 * t355 + (t301 * t341 + t306 * t342) * t165 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155 * t169 + t73, t153 * t169 + t74, -t153 ^ 2 - t155 ^ 2, t153 * t53 + t155 * t52 - t332 + t464, 0, 0, 0, 0, 0, t27 - t477, t26 - t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355 * t99, t355 ^ 2 - t99 ^ 2, t26 + t478, -t27 - t477, t82, -g(1) * t188 + g(2) * t186 - t12 * t480 + t276 * t445 + t355 * t67 + t386, g(1) * t189 - g(2) * t187 + t276 * t444 + t358 * t480 + t67 * t99 - t360;];
tau_reg  = t2;
