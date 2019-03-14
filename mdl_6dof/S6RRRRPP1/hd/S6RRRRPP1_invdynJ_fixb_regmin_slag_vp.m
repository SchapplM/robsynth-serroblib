% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:52
% EndTime: 2019-03-09 20:47:09
% DurationCPUTime: 7.38s
% Computational Cost: add. (14022->578), mult. (31091->715), div. (0->0), fcn. (22516->14), ass. (0->303)
t283 = qJD(2) + qJD(3);
t292 = sin(qJ(3));
t296 = cos(qJ(2));
t438 = cos(qJ(3));
t362 = qJD(1) * t438;
t293 = sin(qJ(2));
t383 = qJD(1) * t293;
t457 = -t292 * t383 + t296 * t362;
t458 = t283 * t457;
t395 = t292 * t296;
t222 = t293 * t438 + t395;
t162 = t283 * t222;
t356 = qJDD(1) * t438;
t378 = t293 * qJDD(1);
t121 = qJD(1) * t162 + t292 * t378 - t296 * t356;
t119 = qJDD(4) + t121;
t288 = sin(pkin(10));
t289 = cos(pkin(10));
t291 = sin(qJ(4));
t295 = cos(qJ(4));
t445 = -t288 * t291 + t289 * t295;
t140 = t445 * t222;
t219 = t288 * t295 + t289 * t291;
t200 = t219 * qJD(4);
t454 = -t219 * t457 + t200;
t380 = qJD(4) * t295;
t381 = qJD(4) * t291;
t453 = -t288 * t381 + t289 * t380 - t445 * t457;
t287 = qJ(2) + qJ(3);
t280 = cos(t287);
t268 = g(3) * t280;
t282 = qJDD(2) + qJDD(3);
t379 = qJD(1) * qJD(2);
t358 = t296 * t379;
t439 = pkin(8) + pkin(7);
t166 = qJDD(2) * pkin(2) + t439 * (-t358 - t378);
t359 = t293 * t379;
t377 = t296 * qJDD(1);
t170 = t439 * (-t359 + t377);
t236 = t439 * t293;
t224 = qJD(1) * t236;
t430 = qJD(2) * pkin(2);
t209 = -t224 + t430;
t237 = t439 * t296;
t226 = qJD(1) * t237;
t361 = qJD(3) * t438;
t382 = qJD(3) * t292;
t348 = -t438 * t166 + t292 * t170 + t209 * t382 + t226 * t361;
t84 = -pkin(3) * t282 + t348;
t456 = t84 + t268;
t204 = -qJD(1) * t395 - t293 * t362;
t171 = -t204 * t291 - t295 * t283;
t335 = t204 * t295 - t283 * t291;
t110 = t289 * t171 - t288 * t335;
t195 = qJD(4) - t457;
t455 = t110 * t195;
t411 = t457 * t291;
t448 = (t381 - t411) * pkin(4);
t279 = sin(t287);
t297 = cos(qJ(1));
t404 = t279 * t297;
t294 = sin(qJ(1));
t405 = t279 * t294;
t452 = g(1) * t404 + g(2) * t405;
t323 = -t292 * t293 + t296 * t438;
t161 = t283 * t323;
t364 = t222 * t380;
t451 = t161 * t291 + t364;
t336 = -t171 * t288 - t289 * t335;
t450 = t336 ^ 2;
t278 = t295 * qJD(5);
t349 = pkin(2) * t361;
t339 = t295 * t349;
t270 = pkin(2) * t292 + pkin(9);
t389 = -qJ(5) - t270;
t351 = qJD(4) * t389;
t160 = t291 * t351 + t278 + t339;
t303 = (-t349 - qJD(5)) * t291 + t295 * t351;
t150 = -pkin(3) * t204 - pkin(9) * t457;
t134 = pkin(2) * t383 + t150;
t130 = t295 * t134;
t207 = t292 * t226;
t159 = -t224 * t438 - t207;
t281 = t295 * qJ(5);
t343 = -t204 * pkin(4) - t281 * t457;
t77 = -t159 * t291 + t130 + t343;
t371 = qJ(5) * t411;
t386 = t291 * t134 + t295 * t159;
t89 = -t371 + t386;
t423 = (-t303 + t77) * t289 + (t160 - t89) * t288;
t290 = -qJ(5) - pkin(9);
t357 = qJD(4) * t290;
t193 = t291 * t357 + t278;
t317 = -t291 * qJD(5) + t295 * t357;
t142 = t295 * t150;
t154 = t209 * t438 - t207;
t80 = -t154 * t291 + t142 + t343;
t387 = t291 * t150 + t295 * t154;
t91 = -t371 + t387;
t418 = (-t317 + t80) * t289 + (t193 - t91) * t288;
t449 = pkin(5) * t454 - qJ(6) * t453 - qJD(6) * t219 + t448;
t208 = t438 * t226;
t158 = -t292 * t224 + t208;
t346 = pkin(2) * t382 - t158;
t284 = qJ(4) + pkin(10);
t277 = cos(t284);
t406 = t277 * t280;
t276 = sin(t284);
t407 = t276 * t280;
t447 = pkin(5) * t406 + qJ(6) * t407;
t446 = -t438 * t236 - t292 * t237;
t345 = g(1) * t297 + g(2) * t294;
t133 = t289 * t193 + t288 * t317;
t234 = pkin(9) * t295 + t281;
t360 = t290 * t291;
t168 = t288 * t234 - t289 * t360;
t169 = t289 * t234 + t288 * t360;
t120 = t292 * t377 + t293 * t356 + t458;
t93 = t295 * t120 + t204 * t381 + t291 * t282 + t283 * t380;
t94 = -qJD(4) * t335 + t291 * t120 - t295 * t282;
t56 = t288 * t93 + t289 * t94;
t57 = -t288 * t94 + t289 * t93;
t444 = -t133 * t110 + t168 * t57 - t169 * t56;
t106 = t289 * t160 + t288 * t303;
t210 = t270 * t295 + t281;
t355 = t389 * t291;
t146 = t288 * t210 - t289 * t355;
t147 = t289 * t210 + t288 * t355;
t443 = -t106 * t110 + t146 * t57 - t147 * t56;
t442 = t268 - t452;
t402 = t280 * t294;
t176 = t276 * t402 + t277 * t297;
t390 = t297 * t276;
t394 = t294 * t277;
t178 = t280 * t390 - t394;
t267 = g(3) * t279;
t436 = pkin(2) * t296;
t273 = pkin(1) + t436;
t197 = pkin(2) * t359 - qJDD(1) * t273;
t72 = t121 * pkin(3) - t120 * pkin(9) + t197;
t68 = t295 * t72;
t305 = t292 * t166 + t170 * t438 + t209 * t361 - t226 * t382;
t83 = t282 * pkin(9) + t305;
t235 = t273 * qJD(1);
t131 = -pkin(3) * t457 + t204 * pkin(9) - t235;
t155 = t292 * t209 + t208;
t137 = pkin(9) * t283 + t155;
t96 = t131 * t291 + t137 * t295;
t18 = t119 * pkin(4) - t93 * qJ(5) - qJD(4) * t96 + qJD(5) * t335 - t291 * t83 + t68;
t321 = t131 * t380 - t137 * t381 + t291 * t72 + t295 * t83;
t21 = -qJ(5) * t94 - qJD(5) * t171 + t321;
t7 = t289 * t18 - t288 * t21;
t369 = -qJDD(6) + t7;
t136 = -t283 * pkin(3) - t154;
t104 = t171 * pkin(4) + qJD(5) + t136;
t53 = t110 * pkin(5) - qJ(6) * t336 + t104;
t441 = g(1) * t178 + g(2) * t176 + t267 * t276 - t53 * t336 + t369;
t440 = t195 ^ 2;
t437 = pkin(2) * t293;
t435 = pkin(5) * t119;
t434 = pkin(5) * t204;
t431 = t295 * pkin(4);
t8 = t288 * t18 + t289 * t21;
t374 = t293 * t430;
t103 = pkin(3) * t162 - pkin(9) * t161 + t374;
t100 = t295 * t103;
t367 = qJD(2) * t439;
t225 = t293 * t367;
t227 = t296 * t367;
t114 = t446 * qJD(3) - t438 * t225 - t292 * t227;
t153 = -pkin(3) * t323 - pkin(9) * t222 - t273;
t175 = -t292 * t236 + t237 * t438;
t167 = t295 * t175;
t333 = -qJ(5) * t161 - qJD(5) * t222;
t25 = t162 * pkin(4) - t291 * t114 + t100 + t333 * t295 + (-t167 + (qJ(5) * t222 - t153) * t291) * qJD(4);
t372 = t291 * t103 + t295 * t114 + t153 * t380;
t33 = -qJ(5) * t364 + (-qJD(4) * t175 + t333) * t291 + t372;
t15 = t288 * t25 + t289 * t33;
t95 = t295 * t131 - t137 * t291;
t74 = qJ(5) * t335 + t95;
t65 = pkin(4) * t195 + t74;
t75 = -qJ(5) * t171 + t96;
t70 = t289 * t75;
t35 = t288 * t65 + t70;
t44 = t288 * t77 + t289 * t89;
t46 = t288 * t80 + t289 * t91;
t144 = t295 * t153;
t90 = -pkin(4) * t323 - t175 * t291 - t222 * t281 + t144;
t385 = t291 * t153 + t167;
t409 = t222 * t291;
t98 = -qJ(5) * t409 + t385;
t55 = t288 * t90 + t289 * t98;
t427 = t288 * t75;
t36 = t288 * t74 + t70;
t426 = t36 * t336;
t425 = t93 * t291;
t424 = -t434 + t423;
t188 = t204 * qJ(6);
t38 = -t188 + t44;
t422 = t106 - t38;
t421 = t346 + t449;
t420 = -t155 + t449;
t419 = -t434 + t418;
t40 = -t188 + t46;
t417 = t133 - t40;
t416 = t136 * t457;
t414 = t171 * t195;
t413 = t335 * t195;
t412 = t195 * t204;
t410 = t204 * t457;
t408 = t222 * t295;
t271 = pkin(3) + t431;
t232 = t280 * t271;
t403 = t280 * t290;
t401 = t280 * t297;
t398 = t291 * t119;
t397 = t291 * t294;
t396 = t291 * t297;
t393 = t294 * t295;
t392 = t295 * t119;
t391 = t295 * t297;
t37 = t289 * t74 - t427;
t388 = qJD(6) - t37;
t285 = t293 ^ 2;
t384 = -t296 ^ 2 + t285;
t376 = t119 * qJ(6) + t8;
t375 = t438 * pkin(2);
t373 = qJD(4) * pkin(9) * t195;
t370 = t280 * t396;
t368 = g(1) * t401 + g(2) * t402 + t267;
t128 = t136 * t380;
t354 = -qJD(4) * t131 - t83;
t352 = -t279 * t290 + t232;
t350 = t195 * t295;
t272 = -t375 - pkin(3);
t347 = -g(1) * t405 + g(2) * t404;
t344 = g(1) * t294 - g(2) * t297;
t342 = -t137 * t380 + t68;
t341 = -pkin(9) * t119 - t416;
t14 = t25 * t289 - t288 * t33;
t34 = t289 * t65 - t427;
t54 = -t288 * t98 + t289 * t90;
t338 = -t110 ^ 2 - t450;
t337 = -t119 * t270 - t416;
t334 = t271 * t279 + t403;
t332 = -t96 * t204 + t291 * t456 + t128;
t331 = t136 * t381 + t95 * t204 + t295 * t452;
t330 = pkin(4) * t409 - t446;
t329 = t352 + t436;
t327 = pkin(5) * t277 + qJ(6) * t276 + t271;
t326 = t345 * t279;
t325 = t345 * t280;
t324 = -0.2e1 * pkin(1) * t379 - pkin(7) * qJDD(2);
t189 = t280 * t397 + t391;
t322 = t161 * t295 - t222 * t381;
t115 = -t292 * t225 + t227 * t438 - t236 * t382 + t237 * t361;
t149 = -pkin(5) * t445 - t219 * qJ(6) - t271;
t316 = pkin(4) * t451 + t115;
t315 = pkin(4) * t397 + t271 * t401 + t297 * t273 - t290 * t404 + t294 * t439;
t299 = qJD(2) ^ 2;
t314 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t299 + t344;
t300 = qJD(1) ^ 2;
t313 = pkin(1) * t300 - pkin(7) * qJDD(1) + t345;
t312 = -t235 * t204 - t348 - t442;
t52 = pkin(4) * t94 + qJDD(5) + t84;
t311 = t297 * t439 + t290 * t405 + pkin(4) * t396 + (-t273 - t232) * t294;
t2 = qJD(6) * t195 + t376;
t28 = -pkin(5) * t195 + qJD(6) - t34;
t29 = qJ(6) * t195 + t35;
t4 = -t369 - t435;
t309 = t2 * t445 + t4 * t219 + t28 * t453 - t29 * t454 - t368;
t308 = -t7 * t219 - t34 * t453 - t35 * t454 + t445 * t8 - t368;
t11 = pkin(5) * t56 - qJ(6) * t57 - qJD(6) * t336 + t52;
t307 = -g(3) * t406 - t11 * t445 - t28 * t204 + t277 * t452 + t454 * t53;
t306 = -g(3) * t407 - t11 * t219 + t29 * t204 + t276 * t452 - t453 * t53;
t302 = t235 * t457 - t305 + t368;
t266 = -pkin(4) * t289 - pkin(5);
t264 = pkin(4) * t288 + qJ(6);
t261 = pkin(4) * t393;
t192 = t280 * t391 + t397;
t191 = -t370 + t393;
t190 = -t280 * t393 + t396;
t179 = t276 * t294 + t277 * t401;
t177 = t280 * t394 - t390;
t139 = t219 * t222;
t135 = -t375 + t149;
t122 = t204 ^ 2 - t457 ^ 2;
t102 = -t204 * t283 - t121;
t101 = t120 - t458;
t79 = -t445 * t161 + t200 * t222;
t78 = -qJD(4) * t140 - t161 * t219;
t69 = t139 * pkin(5) - t140 * qJ(6) + t330;
t61 = -pkin(4) * t335 + pkin(5) * t336 + qJ(6) * t110;
t60 = t195 * t350 - t204 * t335 + t398;
t59 = -t171 * t204 - t291 * t440 + t392;
t58 = -t335 * t350 + t425;
t51 = pkin(5) * t323 - t54;
t50 = -qJ(6) * t323 + t55;
t23 = (t93 - t414) * t295 + (-t94 + t413) * t291;
t22 = -t78 * pkin(5) + t79 * qJ(6) - t140 * qJD(6) + t316;
t13 = -pkin(5) * t162 - t14;
t12 = qJ(6) * t162 - qJD(6) * t323 + t15;
t1 = [qJDD(1), t344, t345, qJDD(1) * t285 + 0.2e1 * t293 * t358, 0.2e1 * t293 * t377 - 0.2e1 * t379 * t384, qJDD(2) * t293 + t296 * t299, qJDD(2) * t296 - t293 * t299, 0, t293 * t324 + t296 * t314, -t293 * t314 + t296 * t324, t120 * t222 - t161 * t204, t120 * t323 - t121 * t222 + t161 * t457 + t162 * t204, t161 * t283 + t222 * t282, -t162 * t283 + t282 * t323, 0, -t115 * t283 - t273 * t121 - t235 * t162 - t197 * t323 + t280 * t344 + t282 * t446 - t374 * t457, -t114 * t283 - t120 * t273 - t161 * t235 - t175 * t282 + t197 * t222 - t204 * t374 + t347, -t322 * t335 + t408 * t93 (-t171 * t295 + t291 * t335) * t161 + (-t425 - t295 * t94 + (t171 * t291 + t295 * t335) * qJD(4)) * t222, -t162 * t335 + t195 * t322 + t222 * t392 - t323 * t93, -t171 * t162 - t195 * t451 - t222 * t398 + t94 * t323, -t119 * t323 + t162 * t195 (-t175 * t380 + t100) * t195 + t144 * t119 - t342 * t323 + t95 * t162 + t115 * t171 - t446 * t94 + t222 * t128 - g(1) * t190 - g(2) * t192 + ((-qJD(4) * t153 - t114) * t195 - t175 * t119 - t354 * t323 + t84 * t222 + t136 * t161) * t291 -(-t175 * t381 + t372) * t195 - t385 * t119 + t321 * t323 - t96 * t162 - t115 * t335 - t446 * t93 + t84 * t408 - g(1) * t189 - g(2) * t191 + t322 * t136, -t110 * t15 - t139 * t8 - t14 * t336 - t140 * t7 + t34 * t79 + t35 * t78 - t54 * t57 - t55 * t56 - t347, -g(1) * t311 - g(2) * t315 + t104 * t316 + t34 * t14 + t35 * t15 + t330 * t52 + t7 * t54 + t8 * t55, g(1) * t177 - g(2) * t179 + t11 * t139 + t110 * t22 - t119 * t51 - t13 * t195 - t162 * t28 + t323 * t4 - t53 * t78 + t56 * t69, -t110 * t12 + t13 * t336 - t139 * t2 + t140 * t4 - t28 * t79 + t29 * t78 - t50 * t56 + t51 * t57 - t347, g(1) * t176 - g(2) * t178 - t11 * t140 + t119 * t50 + t12 * t195 + t162 * t29 - t2 * t323 - t22 * t336 + t53 * t79 - t57 * t69, t2 * t50 + t29 * t12 + t11 * t69 + t53 * t22 + t4 * t51 + t28 * t13 - g(1) * (-t177 * pkin(5) - t176 * qJ(6) + t311) - g(2) * (pkin(5) * t179 + qJ(6) * t178 + t315); 0, 0, 0, -t293 * t300 * t296, t384 * t300, t378, t377, qJDD(2), -g(3) * t296 + t293 * t313, g(3) * t293 + t296 * t313, t410, t122, t101, t102, t282, t158 * t283 + (t282 * t438 - t283 * t382 + t383 * t457) * pkin(2) + t312, t159 * t283 + (t204 * t383 - t282 * t292 - t283 * t361) * pkin(2) + t302, t58, t23, t60, t59, t412, t272 * t94 - t456 * t295 + t337 * t291 + t346 * t171 + (-t270 * t380 - t130 + (-t349 + t159) * t291) * t195 + t331, t272 * t93 + t337 * t295 - t291 * t326 - t346 * t335 + (t270 * t381 - t339 + t386) * t195 + t332, t44 * t110 + t336 * t423 + t308 + t443, t8 * t147 - t7 * t146 + t52 * (t272 - t431) - g(3) * t329 + (t106 - t44) * t35 - t423 * t34 + (t448 + t346) * t104 + t345 * (t334 + t437) t110 * t421 - t146 * t119 + t135 * t56 - t195 * t424 + t307, t38 * t110 + t336 * t424 + t309 + t443, t147 * t119 - t135 * t57 + t195 * t422 - t336 * t421 + t306, t2 * t147 + t11 * t135 + t4 * t146 - g(3) * (t329 + t447) + t421 * t53 + t422 * t29 + t424 * t28 + t345 * (t279 * t327 + t403 + t437); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t410, t122, t101, t102, t282, t155 * t283 + t312, t154 * t283 + t302, t58, t23, t60, t59, t412, -pkin(3) * t94 - t142 * t195 - t155 * t171 + (t154 * t195 + t341) * t291 + (-t456 - t373) * t295 + t331, -pkin(3) * t93 + t387 * t195 + t155 * t335 + t341 * t295 + (-t326 + t373) * t291 + t332, t46 * t110 + t336 * t418 + t308 + t444, t8 * t169 - t7 * t168 - t52 * t271 - g(3) * t352 + (t133 - t46) * t35 - t418 * t34 + (-t155 + t448) * t104 + t345 * t334, t110 * t420 - t168 * t119 + t149 * t56 - t195 * t419 + t307, t40 * t110 + t336 * t419 + t309 + t444, t169 * t119 - t149 * t57 + t195 * t417 - t336 * t420 + t306, t2 * t169 + t11 * t149 + t4 * t168 - g(3) * (t232 + t447) + t420 * t53 + t417 * t29 + t290 * t325 + t419 * t28 + (g(3) * t290 + t345 * t327) * t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t335 * t171, -t171 ^ 2 + t335 ^ 2, t93 + t414, -t94 - t413, t119, -g(1) * t191 + g(2) * t189 + t136 * t335 + t96 * t195 + (t354 + t267) * t291 + t342, g(1) * t192 - g(2) * t190 + t136 * t171 + t195 * t95 + t267 * t295 - t321, t35 * t336 - t426 + (-t288 * t56 - t289 * t57) * pkin(4) + (t37 - t34) * t110, -g(1) * t261 + t34 * t36 - t35 * t37 + (g(2) * t391 + t104 * t335 + t8 * t288 + t7 * t289 + (t325 + t267) * t291) * pkin(4), -t61 * t110 + t36 * t195 + (pkin(5) - t266) * t119 + t441, -t264 * t56 + t266 * t57 + t29 * t336 - t426 + (t28 - t388) * t110, -t277 * t267 - g(1) * t179 - g(2) * t177 - t53 * t110 + t61 * t336 + t264 * t119 + (0.2e1 * qJD(6) - t37) * t195 + t376, t2 * t264 + t4 * t266 - t53 * t61 - t28 * t36 - g(1) * (-pkin(4) * t370 - pkin(5) * t178 + qJ(6) * t179 + t261) - g(2) * (-pkin(4) * t189 - t176 * pkin(5) + t177 * qJ(6)) + t388 * t29 - (-pkin(4) * t291 - pkin(5) * t276 + qJ(6) * t277) * t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t110 * t35 + t336 * t34 + t442 + t52, t195 * t336 + t56, t338, -t57 + t455, t110 * t29 - t28 * t336 + t11 + t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t336 - t119, t57 + t455, -t440 - t450, -t195 * t29 - t435 - t441;];
tau_reg  = t1;