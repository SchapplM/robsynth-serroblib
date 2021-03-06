% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:36:05
% EndTime: 2019-12-31 22:36:37
% DurationCPUTime: 15.25s
% Computational Cost: add. (18105->740), mult. (44840->1008), div. (0->0), fcn. (35440->14), ass. (0->332)
t294 = cos(qJ(2));
t290 = sin(qJ(2));
t285 = sin(pkin(5));
t412 = qJD(1) * t285;
t374 = t290 * t412;
t286 = cos(pkin(5));
t411 = qJD(1) * t286;
t391 = pkin(1) * t411;
t209 = -pkin(7) * t374 + t294 * t391;
t341 = pkin(2) * t290 - pkin(8) * t294;
t210 = t341 * t412;
t289 = sin(qJ(3));
t293 = cos(qJ(3));
t136 = t293 * t209 + t289 * t210;
t410 = qJD(1) * t294;
t373 = t285 * t410;
t352 = t289 * t373;
t295 = -pkin(9) - pkin(8);
t380 = qJD(3) * t295;
t512 = pkin(9) * t352 + t289 * t380 - t136;
t292 = cos(qJ(5));
t402 = qJD(5) * t292;
t266 = qJD(2) + t411;
t351 = t289 * t374;
t188 = t266 * t293 - t351;
t189 = t266 * t289 + t293 * t374;
t288 = sin(qJ(4));
t476 = cos(qJ(4));
t122 = -t476 * t188 + t189 * t288;
t503 = t122 * t292;
t511 = t402 + t503;
t375 = t476 * t293;
t349 = qJD(3) * t375;
t369 = qJD(4) * t476;
t395 = qJD(3) + qJD(4);
t426 = t288 * t289;
t158 = -t293 * t369 + t395 * t426 - t349;
t318 = t375 - t426;
t170 = t318 * t373;
t419 = t158 + t170;
t425 = t288 * t293;
t227 = t476 * t289 + t425;
t418 = (-t373 + t395) * t227;
t399 = qJD(1) * qJD(2);
t366 = t294 * t399;
t397 = qJDD(1) * t290;
t510 = t366 + t397;
t135 = -t209 * t289 + t293 * t210;
t420 = t293 * t294;
t102 = (pkin(3) * t290 - pkin(9) * t420) * t412 + t135;
t252 = t295 * t289;
t253 = t295 * t293;
t319 = t476 * t252 + t288 * t253;
t456 = t319 * qJD(4) - t288 * t102 + t380 * t425 + t512 * t476;
t429 = t285 * t294;
t475 = pkin(1) * t290;
t415 = pkin(7) * t429 + t286 * t475;
t212 = t415 * qJD(1);
t407 = qJD(3) * t289;
t343 = -t212 + (-t352 + t407) * pkin(3);
t287 = sin(qJ(5));
t367 = t290 * t399;
t348 = t285 * t367;
t396 = qJDD(1) * t294;
t263 = t285 * t396;
t394 = t263 - qJDD(3);
t203 = t348 - t394;
t303 = qJDD(4) + t203;
t248 = -qJD(3) + t373;
t324 = -qJD(4) + t248;
t308 = t292 * t324;
t320 = t288 * t188 + t476 * t189;
t403 = qJD(5) * t287;
t398 = qJDD(1) * t286;
t350 = qJDD(2) + t398;
t406 = qJD(3) * t293;
t502 = t510 * t285;
t111 = qJD(3) * t351 - t266 * t406 - t289 * t350 - t293 * t502;
t408 = qJD(2) * t294;
t371 = t289 * t408;
t112 = t266 * t407 + (qJD(1) * (t290 * t406 + t371) + t289 * t397) * t285 - t293 * t350;
t405 = qJD(4) * t288;
t53 = t476 * t111 + t288 * t112 - t188 * t369 + t189 * t405;
t29 = qJD(5) * t308 - t287 * t303 + t292 * t53 + t320 * t403;
t24 = t29 * t287;
t96 = -t287 * t324 + t292 * t320;
t509 = t511 * t96 - t24;
t359 = t288 * t111 - t476 * t112;
t54 = qJD(4) * t320 - t359;
t50 = qJDD(5) + t54;
t500 = qJD(5) + t122;
t508 = t287 * t50 - t320 * t96 + t511 * t500;
t432 = t285 * t290;
t215 = -t286 * t293 + t289 * t432;
t507 = pkin(3) * t215;
t477 = cos(qJ(1));
t376 = t477 * t294;
t291 = sin(qJ(1));
t424 = t290 * t291;
t221 = -t286 * t424 + t376;
t284 = qJ(3) + qJ(4);
t278 = sin(t284);
t279 = cos(t284);
t431 = t285 * t291;
t153 = t221 * t278 - t279 * t431;
t198 = -t278 * t432 + t279 * t286;
t377 = t477 * t290;
t422 = t291 * t294;
t219 = t286 * t377 + t422;
t379 = t285 * t477;
t358 = -t219 * t278 - t279 * t379;
t314 = g(1) * t153 - g(2) * t358 - g(3) * t198;
t168 = pkin(8) * t266 + t212;
t327 = -pkin(2) * t294 - pkin(8) * t290 - pkin(1);
t202 = t327 * t285;
t180 = qJD(1) * t202;
t109 = t168 * t293 + t180 * t289;
t355 = qJD(2) * t391;
t387 = pkin(1) * t398;
t381 = pkin(7) * t263 + t290 * t387 + t294 * t355;
t141 = -pkin(7) * t348 + t381;
t130 = pkin(8) * t350 + t141;
t321 = t341 * qJD(2);
t134 = (qJD(1) * t321 + qJDD(1) * t327) * t285;
t58 = -qJD(3) * t109 - t289 * t130 + t293 * t134;
t34 = pkin(3) * t203 + pkin(9) * t111 + t58;
t315 = -t293 * t130 - t289 * t134 + t168 * t407 - t180 * t406;
t38 = -pkin(9) * t112 - t315;
t362 = t288 * t38 - t476 * t34;
t87 = pkin(9) * t188 + t109;
t388 = t476 * t87;
t108 = -t168 * t289 + t293 * t180;
t86 = -pkin(9) * t189 + t108;
t80 = -pkin(3) * t248 + t86;
t46 = t288 * t80 + t388;
t10 = -qJD(4) * t46 - t362;
t8 = -pkin(4) * t303 - t10;
t306 = t314 - t8;
t506 = t418 * pkin(4) + t419 * pkin(10) + t343;
t457 = t288 * t87;
t45 = t476 * t80 - t457;
t43 = pkin(4) * t324 - t45;
t505 = t122 * t43;
t504 = -pkin(10) * t374 + t456;
t446 = t122 * t320;
t430 = t285 * t293;
t163 = -t221 * t289 + t291 * t430;
t501 = t108 * t248 - t315;
t267 = pkin(7) * t432;
t474 = pkin(1) * t294;
t223 = t286 * t474 - t267;
t213 = qJD(2) * t223;
t499 = -t122 ^ 2 + t320 ^ 2;
t78 = pkin(4) * t320 + pkin(10) * t122;
t498 = -t122 * t324 - t53;
t167 = -pkin(2) * t266 - t209;
t126 = -pkin(3) * t188 + t167;
t150 = t219 * t279 - t278 * t379;
t154 = t221 * t279 + t278 * t431;
t199 = t278 * t286 + t279 * t432;
t313 = -g(1) * t154 - g(2) * t150 - g(3) * t199;
t9 = t288 * t34 + t80 * t369 + t476 * t38 - t87 * t405;
t497 = t126 * t122 - t313 - t9;
t452 = qJD(5) * t96;
t30 = -t287 * t53 - t292 * t303 + t452;
t94 = t287 * t320 + t308;
t331 = t287 * t96 + t292 * t94;
t496 = -t122 * t331 - t287 * t30 - t29 * t292 - t94 * t402 - t96 * t403;
t218 = -t286 * t376 + t424;
t495 = t150 * t287 - t218 * t292;
t494 = t150 * t292 + t218 * t287;
t281 = t285 ^ 2;
t492 = 0.2e1 * t281;
t353 = pkin(7) * t502 + t290 * t355 - t294 * t387;
t131 = -pkin(2) * t350 + t353;
t77 = t112 * pkin(3) + t131;
t12 = t54 * pkin(4) + t53 * pkin(10) + t77;
t44 = -pkin(10) * t324 + t46;
t59 = pkin(4) * t122 - pkin(10) * t320 + t126;
t332 = t287 * t44 - t292 * t59;
t7 = pkin(10) * t303 + t9;
t2 = -t332 * qJD(5) + t287 * t12 + t292 * t7;
t1 = t2 * t292;
t18 = t287 * t59 + t292 * t44;
t3 = -qJD(5) * t18 + t292 * t12 - t287 * t7;
t490 = -t3 * t287 + t1;
t489 = -t18 * t500 - t3;
t173 = t288 * t252 - t476 * t253;
t454 = t173 * qJD(4) + t476 * t102 + t512 * t288 - t295 * t349;
t55 = t288 * t86 + t388;
t345 = pkin(3) * t405 - t55;
t488 = t109 * t248 - t58;
t487 = t500 * t320;
t220 = t286 * t422 + t377;
t305 = -g(1) * t220 - g(2) * t218 + g(3) * t429;
t486 = t305 * t278;
t201 = pkin(8) * t286 + t415;
t133 = t293 * t201 + t289 * t202;
t214 = t415 * qJD(2);
t56 = t476 * t86 - t457;
t473 = pkin(3) * t189;
t68 = t473 + t78;
t20 = t287 * t68 + t292 * t56;
t472 = pkin(3) * t288;
t275 = pkin(10) + t472;
t356 = pkin(3) * t369;
t485 = -t275 * t403 + t292 * t356 - t20;
t484 = t320 * t332 + t43 * t403;
t48 = t292 * t50;
t483 = t320 * t94 + t48;
t482 = t18 * t320 - t306 * t287 + t43 * t402;
t480 = -t126 * t320 + t314 - t362;
t479 = -t320 * t248 + t359;
t101 = -pkin(9) * t215 + t133;
t132 = -t201 * t289 + t293 * t202;
t216 = t286 * t289 + t290 * t430;
t92 = -pkin(3) * t429 - pkin(9) * t216 + t132;
t462 = t476 * t101 + t288 * t92;
t161 = -qJD(3) * t215 + t408 * t430;
t409 = qJD(2) * t290;
t372 = t285 * t409;
t211 = t285 * t321;
t82 = -t133 * qJD(3) + t293 * t211 - t213 * t289;
t67 = pkin(3) * t372 - pkin(9) * t161 + t82;
t160 = qJD(3) * t216 + t285 * t371;
t81 = -t201 * t407 + t202 * t406 + t289 * t211 + t293 * t213;
t72 = -pkin(9) * t160 + t81;
t16 = -qJD(4) * t462 - t288 * t72 + t476 * t67;
t468 = g(3) * t285;
t466 = t96 * t94;
t277 = pkin(3) * t293 + pkin(2);
t146 = -pkin(4) * t318 - pkin(10) * t227 - t277;
t98 = t146 * t292 - t173 * t287;
t465 = qJD(5) * t98 + t506 * t287 + t504 * t292;
t99 = t146 * t287 + t173 * t292;
t464 = -qJD(5) * t99 - t504 * t287 + t506 * t292;
t461 = t332 * t500;
t460 = t332 * t287;
t458 = t287 * t94;
t26 = t30 * t292;
t455 = pkin(4) * t374 + t454;
t449 = t500 * t287;
t445 = t188 * t248;
t444 = t189 * t188;
t443 = t189 * t248;
t440 = t219 * t289;
t438 = t227 * t287;
t437 = t227 * t292;
t436 = t248 * t289;
t435 = t279 * t287;
t434 = t279 * t292;
t433 = t281 * qJD(1) ^ 2;
t427 = t287 * t294;
t423 = t290 * t295;
t421 = t292 * t294;
t417 = -t218 * t277 - t219 * t295;
t416 = -t220 * t277 - t221 * t295;
t414 = t477 * pkin(1) + pkin(7) * t431;
t282 = t290 ^ 2;
t283 = t294 ^ 2;
t413 = t282 - t283;
t404 = qJD(5) * t500;
t401 = qJD(2) - t266;
t392 = t476 * pkin(3);
t386 = t294 * t433;
t384 = t289 * t431;
t382 = t285 * t421;
t255 = t285 * t427;
t378 = t293 * t477;
t368 = pkin(1) * t492;
t363 = -pkin(1) * t291 + pkin(7) * t379;
t360 = qJDD(4) - t394;
t256 = t289 * t379;
t357 = t219 * t293 - t256;
t354 = t290 * t386;
t347 = t290 * t366;
t344 = t163 * pkin(3);
t340 = pkin(4) * t279 + pkin(10) * t278;
t339 = g(1) * t358 + g(2) * t153;
t338 = g(1) * t218 - g(2) * t220;
t337 = g(1) * t221 + g(2) * t219;
t334 = -t275 * t50 + t505;
t333 = t18 * t287 - t292 * t332;
t61 = -pkin(10) * t429 + t462;
t139 = t476 * t215 + t216 * t288;
t140 = -t288 * t215 + t476 * t216;
t200 = t267 + (-pkin(2) - t474) * t286;
t143 = t200 + t507;
t73 = pkin(4) * t139 - pkin(10) * t140 + t143;
t28 = t287 * t73 + t292 * t61;
t27 = -t287 * t61 + t292 * t73;
t329 = pkin(3) * t384 - t220 * t295 + t221 * t277 + t414;
t326 = g(1) * t477 + g(2) * t291;
t62 = -t288 * t101 + t476 * t92;
t117 = t140 * t287 + t382;
t15 = -t101 * t405 + t288 * t67 + t92 * t369 + t476 * t72;
t317 = pkin(3) * t256 + t218 * t295 - t219 * t277 + t363;
t312 = -t285 * t378 - t440;
t137 = t170 * t287 - t292 * t374;
t311 = -t287 * t158 + t227 * t402 - t137;
t138 = t170 * t292 + t287 * t374;
t310 = -t292 * t158 - t227 * t403 - t138;
t307 = t1 + t313;
t304 = -g(3) * t432 - t337;
t302 = t312 * pkin(3);
t301 = -t131 - t305;
t300 = -pkin(8) * t203 - t167 * t248;
t127 = pkin(3) * t160 + t214;
t297 = pkin(8) * qJD(3) * t248 + t301;
t276 = -t392 - pkin(4);
t229 = t277 * t429;
t190 = t198 * pkin(4);
t164 = t221 * t293 + t384;
t148 = t153 * pkin(4);
t147 = t358 * pkin(4);
t118 = t140 * t292 - t255;
t115 = t154 * t292 + t220 * t287;
t114 = -t154 * t287 + t220 * t292;
t76 = t140 * qJD(4) + t476 * t160 + t288 * t161;
t75 = t288 * t160 - t476 * t161 + t215 * t369 + t216 * t405;
t60 = pkin(4) * t429 - t62;
t52 = -qJD(5) * t255 + t140 * t402 - t287 * t75 - t292 * t372;
t51 = qJD(5) * t117 - t287 * t372 + t292 * t75;
t31 = pkin(4) * t76 + pkin(10) * t75 + t127;
t22 = t287 * t78 + t292 * t45;
t21 = -t287 * t45 + t292 * t78;
t19 = -t287 * t56 + t292 * t68;
t14 = -pkin(4) * t372 - t16;
t13 = pkin(10) * t372 + t15;
t5 = -qJD(5) * t28 - t13 * t287 + t292 * t31;
t4 = qJD(5) * t27 + t13 * t292 + t287 * t31;
t6 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t291 - g(2) * t477, t326, 0, 0, (qJDD(1) * t282 + 0.2e1 * t347) * t281, (t290 * t396 - t399 * t413) * t492, (t266 * t408 + t290 * qJDD(2) + (t366 + 0.2e1 * t397) * t286) * t285, (qJDD(1) * t283 - 0.2e1 * t347) * t281, (-t266 * t409 + t294 * qJDD(2) + (-t367 + 0.2e1 * t396) * t286) * t285, t350 * t286, -t214 * t266 + t223 * t350 - t353 * t286 + g(1) * t219 - g(2) * t221 + (-t367 + t396) * t368, -t141 * t286 - t213 * t266 - t415 * t350 - t510 * t368 - t338, ((-qJD(2) * t209 + t415 * qJDD(1) + t141) * t294 + (-qJD(2) * t212 - t223 * qJDD(1) + t353) * t290 - t326) * t285, t281 * qJDD(1) * pkin(1) ^ 2 - g(1) * t363 - g(2) * t414 + t141 * t415 - t209 * t214 + t212 * t213 - t223 * t353, -t111 * t216 + t161 * t189, t111 * t215 - t112 * t216 - t160 * t189 + t161 * t188, -t161 * t248 + t203 * t216 + (t111 * t294 + t189 * t409) * t285, t112 * t215 - t160 * t188, t160 * t248 - t203 * t215 + (t112 * t294 + t188 * t409) * t285, (-t203 * t294 - t248 * t409) * t285, -t82 * t248 + t132 * t203 - t214 * t188 + t200 * t112 + t131 * t215 + t167 * t160 + g(1) * t357 - g(2) * t164 + (t108 * t409 - t294 * t58) * t285, -g(1) * t440 - g(2) * t163 - t200 * t111 + t131 * t216 - t133 * t203 + t167 * t161 + t214 * t189 + t81 * t248 + (-g(1) * t378 - t109 * t409 - t294 * t315) * t285, -t108 * t161 - t109 * t160 + t111 * t132 - t112 * t133 + t188 * t81 - t189 * t82 + t215 * t315 - t216 * t58 + t338, -t315 * t133 + t109 * t81 + t58 * t132 + t108 * t82 + t131 * t200 + t167 * t214 - g(1) * (-pkin(2) * t219 - pkin(8) * t218 + t363) - g(2) * (pkin(2) * t221 + pkin(8) * t220 + t414), -t140 * t53 - t320 * t75, t122 * t75 + t139 * t53 - t140 * t54 - t320 * t76, -t75 * t395 + t140 * t360 + (t320 * t409 + t53 * t294 + (t140 * t409 + t294 * t75) * qJD(1)) * t285, t122 * t76 + t139 * t54, -t76 * t395 - t139 * t360 + (-t122 * t409 + t54 * t294 + (-t139 * t409 + t294 * t76) * qJD(1)) * t285, (-t360 * t294 + (-0.2e1 * t373 + t395) * t409) * t285, t16 * t395 + t62 * t360 + t127 * t122 + t143 * t54 + t77 * t139 + t126 * t76 + g(1) * t150 - g(2) * t154 + (t45 * t409 - t10 * t294 + (-t16 * t294 + t409 * t62) * qJD(1)) * t285, -t15 * t395 - t462 * t360 + t127 * t320 - t143 * t53 + t77 * t140 - t126 * t75 + (-t46 * t409 + t9 * t294 + (t15 * t294 - t409 * t462) * qJD(1)) * t285 + t339, -t10 * t140 - t122 * t15 - t139 * t9 - t16 * t320 + t45 * t75 - t46 * t76 - t462 * t54 + t53 * t62 + t338, -g(1) * t317 - g(2) * t329 + t10 * t62 + t126 * t127 + t77 * t143 + t46 * t15 + t45 * t16 + t462 * t9, -t118 * t29 - t51 * t96, t117 * t29 - t118 * t30 + t51 * t94 - t52 * t96, t118 * t50 - t139 * t29 - t500 * t51 + t76 * t96, t117 * t30 + t52 * t94, -t117 * t50 - t139 * t30 - t500 * t52 - t76 * t94, t139 * t50 + t500 * t76, g(1) * t494 - g(2) * t115 + t8 * t117 + t3 * t139 + t14 * t94 + t27 * t50 + t60 * t30 - t332 * t76 + t43 * t52 + t5 * t500, -g(1) * t495 - g(2) * t114 + t8 * t118 - t2 * t139 + t14 * t96 - t18 * t76 - t28 * t50 - t60 * t29 - t4 * t500 - t43 * t51, -t117 * t2 - t118 * t3 - t18 * t52 + t27 * t29 - t28 * t30 - t332 * t51 - t4 * t94 - t5 * t96 - t339, t2 * t28 + t18 * t4 + t3 * t27 - t332 * t5 + t8 * t60 + t43 * t14 - g(1) * (-pkin(4) * t150 + pkin(10) * t358 + t317) - g(2) * (pkin(4) * t154 + pkin(10) * t153 + t329); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t354, t413 * t433, (t401 * t410 + t397) * t285, t354, -t374 * t401 + t263, t350, t212 * t266 + t433 * t475 - t305 - t353, pkin(1) * t386 + t209 * t266 + (pkin(7) * t399 + g(3)) * t432 + t337 - t381, 0, 0, -t111 * t289 - t293 * t443, (-t111 - t445) * t293 + (-t112 + t443) * t289, -t248 * t406 + t203 * t289 + (-t189 * t290 + t248 * t420) * t412, -t112 * t293 + t188 * t436, t248 * t407 + t203 * t293 + (-t188 * t290 - t294 * t436) * t412, t248 * t374, -pkin(2) * t112 - t108 * t374 + t135 * t248 + t188 * t212 + t289 * t300 + t293 * t297, pkin(2) * t111 + t109 * t374 - t136 * t248 - t189 * t212 - t289 * t297 + t293 * t300, t135 * t189 - t136 * t188 + ((qJD(3) * t189 - t112) * pkin(8) + t501) * t293 + ((-qJD(3) * t188 - t111) * pkin(8) + t488) * t289 + t304, -t108 * t135 - t109 * t136 - t167 * t212 + t301 * pkin(2) + (-t58 * t289 - t315 * t293 + (-t108 * t293 - t109 * t289) * qJD(3) + t304) * pkin(8), -t227 * t53 - t320 * t419, t122 * t419 - t227 * t54 - t318 * t53 - t320 * t418, t227 * t303 - t320 * t374 + t419 * t324, t122 * t418 - t318 * t54, t122 * t374 + t303 * t318 + t418 * t324, t324 * t374, t122 * t343 + t126 * t418 - t277 * t54 - t279 * t305 + t303 * t319 - t318 * t77 + t454 * t324 - t374 * t45, -t126 * t419 - t173 * t303 + t77 * t227 + t277 * t53 + t320 * t343 + t456 * t324 + t374 * t46 + t486, -t10 * t227 - t122 * t456 - t173 * t54 + t318 * t9 + t319 * t53 + t320 * t454 - t418 * t46 + t419 * t45 + t304, t9 * t173 + t10 * t319 - t77 * t277 - g(1) * t416 - g(2) * t417 - g(3) * (-t285 * t423 + t229) + t456 * t46 - t454 * t45 + t343 * t126, -t29 * t437 + t310 * t96, t137 * t96 + t138 * t94 + t331 * t158 + (t24 - t26 + (-t292 * t96 + t458) * qJD(5)) * t227, t29 * t318 + t310 * t500 + t418 * t96 + t437 * t50, t30 * t438 + t311 * t94, t30 * t318 - t311 * t500 - t418 * t94 - t438 * t50, -t318 * t50 + t418 * t500, t98 * t50 - t3 * t318 - t319 * t30 + t8 * t438 - g(1) * (-t220 * t434 + t221 * t287) - g(2) * (-t218 * t434 + t219 * t287) + t455 * t94 - (t279 * t421 + t287 * t290) * t468 - t418 * t332 + t464 * t500 + t311 * t43, -t99 * t50 + t2 * t318 + t319 * t29 + t8 * t437 - g(1) * (t220 * t435 + t221 * t292) - g(2) * (t218 * t435 + t219 * t292) + t455 * t96 - (-t279 * t427 + t290 * t292) * t468 - t418 * t18 - t465 * t500 + t310 * t43, t137 * t18 - t138 * t332 + t29 * t98 - t30 * t99 - t464 * t96 - t465 * t94 + t333 * t158 - t486 + (-t2 * t287 - t292 * t3 + (-t18 * t292 - t460) * qJD(5)) * t227, t2 * t99 + t3 * t98 - t8 * t319 - g(1) * (-t220 * t340 + t416) - g(2) * (-t218 * t340 + t417) - g(3) * t229 + t455 * t43 - (t294 * t340 - t423) * t468 + t465 * t18 - t464 * t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444, -t188 ^ 2 + t189 ^ 2, -t111 + t445, t444, -t112 - t443, t203, -g(1) * t163 - g(2) * t312 + g(3) * t215 - t167 * t189 - t488, g(1) * t164 + g(2) * t357 + g(3) * t216 - t167 * t188 - t501, 0, 0, t446, t499, t498, -t446, t479, t303, -t55 * t248 + (t55 - t46) * qJD(4) + (-t189 * t122 + t476 * t303 + t324 * t405) * pkin(3) + t480, -t56 * t324 + (-t189 * t320 - t288 * t303 + t324 * t369) * pkin(3) + t497, t46 * t320 + t56 * t122 - t45 * t122 - t55 * t320 + (t476 * t53 - t288 * t54 + (-t476 * t122 + t288 * t320) * qJD(4)) * pkin(3), -g(1) * t344 - g(2) * t302 + g(3) * t507 + t10 * t392 - t126 * t473 + t472 * t9 + (t356 - t56) * t46 - t345 * t45, t509, t496, t508, t458 * t500 - t26, -t449 * t500 + t483, -t487, -t19 * t500 + t276 * t30 + t345 * t94 + (-t356 * t500 + t334) * t287 + (-t275 * t404 + t306) * t292 + t484, -t276 * t29 + t334 * t292 + t345 * t96 - t485 * t500 + t482, t19 * t96 + t20 * t94 + (-t94 * t356 + t122 * t332 - t275 * t30 + (t275 * t96 + t332) * qJD(5)) * t292 + (t96 * t356 - t122 * t18 - t275 * t29 - t3 + (t275 * t94 - t18) * qJD(5)) * t287 + t307, t356 * t460 + t8 * t276 + t332 * t19 - g(1) * (pkin(10) * t154 - t148 + t344) - g(2) * (t150 * pkin(10) + t147 + t302) - g(3) * (pkin(10) * t199 + t190 - t507) + t345 * t43 + (t332 * t402 + t490) * t275 + t485 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t446, t499, t498, -t446, t479, t303, -t248 * t46 + t480, -t324 * t45 + t497, 0, 0, t509, t496, t508, t449 * t94 - t26, -t287 * t500 ^ 2 + t483, -t487, -pkin(4) * t30 - t500 * t21 - t46 * t94 + (-pkin(10) * t50 + t505) * t287 + (-pkin(10) * t404 + t306) * t292 + t484, t43 * t503 + pkin(4) * t29 + t500 * t22 - t46 * t96 + (t403 * t500 - t48) * pkin(10) + t482, t21 * t96 + t22 * t94 + (t461 + (-t30 + t452) * pkin(10)) * t292 + ((qJD(5) * t94 - t29) * pkin(10) + t489) * t287 + t307, -t8 * pkin(4) + g(1) * t148 - g(2) * t147 - g(3) * t190 + t332 * t21 - t18 * t22 - t43 * t46 + (-qJD(5) * t333 + t313 + t490) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t466, -t94 ^ 2 + t96 ^ 2, t500 * t94 - t29, -t466, t500 * t96 - t30, t50, -t43 * t96 - g(1) * t114 + g(2) * t495 - g(3) * (-t199 * t287 - t382) - t489, -t461 + t43 * t94 + g(1) * t115 + g(2) * t494 - g(3) * (-t199 * t292 + t255) - t2, 0, 0;];
tau_reg = t6;
