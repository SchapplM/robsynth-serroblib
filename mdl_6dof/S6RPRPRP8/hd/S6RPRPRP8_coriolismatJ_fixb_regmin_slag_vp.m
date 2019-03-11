% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRP8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:07
% EndTime: 2019-03-09 03:26:21
% DurationCPUTime: 6.47s
% Computational Cost: add. (7176->429), mult. (12579->548), div. (0->0), fcn. (13616->6), ass. (0->338)
t295 = cos(qJ(5));
t294 = sin(qJ(3));
t476 = sin(pkin(9));
t276 = t476 * t294;
t296 = cos(qJ(3));
t477 = cos(pkin(9));
t356 = t477 * t296;
t265 = t356 - t276;
t263 = t265 ^ 2;
t266 = -t294 * t477 - t296 * t476;
t264 = t266 ^ 2;
t517 = t263 + t264;
t524 = t517 * t295;
t529 = t524 * qJD(1);
t293 = sin(qJ(5));
t525 = t517 * t293;
t528 = t525 * qJD(1);
t527 = t525 * qJD(4);
t427 = qJD(2) * t266;
t526 = qJD(4) * t524 + t293 * t427;
t258 = t265 * qJ(6);
t496 = t296 * pkin(3);
t172 = t265 * pkin(4) - t266 * pkin(8) + t496;
t137 = t293 * t172;
t297 = -pkin(1) - pkin(7);
t439 = -qJ(4) + t297;
t268 = t439 * t294;
t269 = t439 * t296;
t200 = t268 * t476 - t477 * t269;
t522 = t200 * t295;
t438 = t137 / 0.2e1 - t522 / 0.2e1;
t523 = -t258 - t438;
t516 = t477 * t268 + t476 * t269;
t468 = t516 * t293;
t467 = t516 * t295;
t345 = t264 / 0.2e1 + t263 / 0.2e1;
t390 = -t263 / 0.2e1;
t521 = -t264 / 0.2e1 + t390;
t161 = t293 * t266;
t411 = t161 * qJD(1);
t422 = qJD(5) * t293;
t127 = t411 - t422;
t165 = t293 * t265;
t399 = t266 * qJD(5);
t376 = t295 * t399;
t340 = -t165 * qJD(3) + t376;
t327 = 0.1e1 / 0.2e1 + t345;
t312 = t327 * t295;
t514 = t312 * qJD(1);
t520 = t514 - t340;
t519 = 0.2e1 * t293;
t495 = qJD(3) * pkin(3);
t518 = (t265 * t476 + t266 * t477) * t495;
t166 = t295 * t265;
t378 = t293 * t399;
t339 = t166 * qJD(3) + t378;
t515 = qJD(2) * t312;
t400 = t266 * qJD(1);
t249 = t295 * t400;
t513 = -qJD(5) * t312 + t249;
t454 = t293 * qJ(6);
t497 = t295 * pkin(5);
t338 = t454 + t497;
t393 = t295 * qJD(6);
t512 = qJD(5) * t338 - t393;
t498 = t266 * pkin(5);
t281 = t294 * pkin(3) + qJ(2);
t308 = -t266 * pkin(4) - t265 * pkin(8) + t281;
t88 = -t295 * t308 + t468;
t73 = t88 + t498;
t383 = t88 / 0.2e1 - t73 / 0.2e1;
t460 = t266 * qJ(6);
t89 = t293 * t308 + t467;
t72 = t89 - t460;
t384 = t72 / 0.2e1 - t89 / 0.2e1;
t511 = t293 * t384 + t295 * t383;
t250 = pkin(5) * t165;
t344 = -qJ(6) * t166 + t250;
t94 = t200 + t344;
t510 = -t94 / 0.2e1;
t509 = -t172 / 0.2e1;
t508 = -t250 / 0.2e1;
t507 = -t265 / 0.2e1;
t506 = t265 / 0.2e1;
t452 = t295 * qJ(6);
t500 = pkin(5) * t293;
t337 = -t452 + t500;
t505 = t337 / 0.2e1;
t291 = t293 ^ 2;
t504 = t291 / 0.2e1;
t292 = t295 ^ 2;
t503 = t292 / 0.2e1;
t289 = -t293 / 0.2e1;
t288 = t293 / 0.2e1;
t502 = -t295 / 0.2e1;
t501 = t295 / 0.2e1;
t499 = t265 * pkin(5);
t157 = t338 * t265;
t300 = t157 * t507 + t266 * t511;
t322 = -t454 / 0.2e1 - t497 / 0.2e1;
t7 = t300 + t322;
t494 = t7 * qJD(1);
t493 = t72 * t266;
t492 = t72 * t293;
t491 = t72 * t295;
t490 = t73 * t293;
t489 = t73 * t295;
t435 = -t522 + t137;
t76 = t258 + t435;
t488 = t76 * t295;
t451 = t295 * t172;
t453 = t200 * t293;
t77 = -t451 - t453 - t499;
t487 = t77 * t293;
t8 = t94 * t157 - t72 * t88 + t73 * t89;
t486 = t8 * qJD(1);
t485 = t88 * t266;
t484 = t89 * t266;
t163 = t295 * t266;
t9 = -t73 * t163 - t77 * t166 + (t265 * t76 + t493) * t293;
t483 = t9 * qJD(1);
t482 = t94 * t293;
t481 = t94 * t295;
t95 = t266 * t337 + t516;
t480 = t95 * t293;
t479 = t95 * t295;
t280 = -pkin(3) * t477 - pkin(4);
t259 = -t338 + t280;
t359 = t166 / 0.2e1;
t478 = t259 * t359 + t482 / 0.2e1;
t10 = -t89 * t166 + (t491 + (t73 - t88) * t293) * t265;
t475 = t10 * qJD(1);
t319 = t460 / 0.2e1 - t384;
t323 = -t498 / 0.2e1 + t383;
t11 = -t293 * t319 + t295 * t323;
t474 = t11 * qJD(1);
t12 = t293 * t323 + t295 * t319;
t473 = t12 * qJD(1);
t15 = (-t76 - t481) * t266 + (t72 - t479) * t265;
t472 = t15 * qJD(1);
t16 = (t77 + t482) * t266 + (-t73 + t480) * t265;
t471 = t16 * qJD(1);
t19 = t94 * t265 + (t490 + t491) * t266;
t470 = t19 * qJD(1);
t22 = -t451 * t266 + (-t88 + t468) * t265;
t466 = t22 * qJD(1);
t23 = (t435 + t522) * t266 + (-t89 + t467) * t265;
t465 = t23 * qJD(1);
t357 = t503 + t504;
t365 = t259 * t506;
t279 = pkin(3) * t476 + pkin(8);
t456 = t279 * t266;
t304 = t357 * t456 + t365;
t320 = t289 * t76 + t501 * t77;
t25 = t304 + t320;
t464 = t25 * qJD(1);
t26 = t484 + (t157 * t293 + t481) * t265;
t463 = t26 * qJD(1);
t462 = t265 * t266;
t461 = t265 * t279;
t459 = t266 * t259;
t27 = t485 + (-t157 * t295 + t482) * t265;
t458 = t27 * qJD(1);
t457 = t337 * t293;
t455 = t280 * t265;
t30 = -t166 * t94 - t493;
t450 = t30 * qJD(1);
t31 = -t489 + t492;
t449 = t31 * qJD(1);
t39 = -t165 * t200 - t485;
t446 = t39 * qJD(1);
t40 = t166 * t200 + t484;
t445 = t40 * qJD(1);
t311 = -t264 * t357 + t390;
t59 = t311 - t357;
t444 = t59 * qJD(1);
t71 = t200 * t265 + t266 * t516;
t443 = t71 * qJD(1);
t442 = t88 * qJD(5);
t348 = t521 * t293;
t96 = t289 + t348;
t441 = t96 * qJD(1);
t436 = t451 / 0.2e1 + t453 / 0.2e1;
t433 = t265 * t503 + t291 * t507;
t431 = t291 + t292;
t274 = t292 - t291;
t430 = qJD(1) * qJ(2);
t429 = qJD(1) * t265;
t428 = qJD(1) * t295;
t426 = qJD(3) * t293;
t425 = qJD(3) * t295;
t424 = qJD(4) * t293;
t423 = qJD(4) * t295;
t421 = qJD(5) * t295;
t420 = t327 * qJD(1);
t347 = t345 * t293;
t112 = t288 + t347;
t419 = t112 * qJD(1);
t115 = t295 * t345 + t502;
t417 = t115 * qJD(2);
t385 = t263 - t264;
t119 = t385 * t293;
t416 = t119 * qJD(1);
t121 = t385 * t295;
t414 = t121 * qJD(1);
t305 = t476 * t266 / 0.2e1 + t477 * t507;
t130 = (-t296 / 0.2e1 + t305) * pkin(3);
t412 = t130 * qJD(1);
t148 = t163 * qJD(1);
t410 = t165 * qJD(1);
t408 = t166 * qJD(1);
t174 = t431 * t265;
t406 = t174 * qJD(1);
t252 = t291 * t266;
t253 = t292 * t266;
t175 = -t252 - t253;
t405 = t175 * qJD(1);
t402 = t517 * qJD(1);
t260 = t276 / 0.2e1 - t356 / 0.2e1;
t401 = t260 * qJD(1);
t257 = t266 * qJD(6);
t273 = t294 ^ 2 - t296 ^ 2;
t398 = t273 * qJD(1);
t397 = t281 * qJD(1);
t396 = t293 * qJD(6);
t395 = t294 * qJD(1);
t394 = t294 * qJD(3);
t392 = t296 * qJD(1);
t391 = t296 * qJD(3);
t254 = t499 / 0.2e1;
t362 = -t456 / 0.2e1;
t382 = -t259 * t166 / 0.2e1 + t295 * t362 - t482 / 0.2e1;
t381 = t254 + t436;
t380 = qJ(2) * t395;
t379 = qJ(2) * t392;
t377 = t279 * t422;
t375 = t279 * t421;
t374 = t265 * t400;
t373 = qJD(3) * t462;
t372 = t293 * t421;
t371 = t293 * t429;
t370 = t265 * t396;
t247 = t293 * t400;
t275 = t293 * t425;
t368 = t293 * t393;
t367 = t265 * t428;
t248 = t265 * t425;
t366 = t294 * t392;
t364 = t337 * t506;
t363 = t457 / 0.2e1;
t361 = t456 / 0.2e1;
t360 = t452 / 0.2e1;
t117 = -t260 + t433;
t354 = t117 * qJD(1) + t275;
t158 = (t504 - t292 / 0.2e1) * t265;
t353 = t158 * qJD(1) - t275;
t225 = t263 * t293 * t428;
t352 = t158 * qJD(3) + t225;
t351 = -qJD(5) + t400;
t350 = t293 * t367;
t349 = t293 * t248;
t346 = -t499 / 0.2e1 - t453 / 0.2e1;
t343 = -t157 / 0.2e1 + t361;
t342 = -t165 * qJD(5) + t266 * t425;
t321 = t491 / 0.2e1 + t490 / 0.2e1;
t4 = (-t488 / 0.2e1 + t510 - t487 / 0.2e1) * t266 + (-t95 / 0.2e1 + t321) * t265;
t5 = t72 * t76 + t73 * t77 + t94 * t95;
t336 = t5 * qJD(1) + t4 * qJD(2);
t335 = t487 + t488;
t74 = (0.1e1 - t431) * t462;
t334 = t4 * qJD(1) + t74 * qJD(2);
t41 = t281 * t496;
t333 = t41 * qJD(1);
t332 = t459 - t461;
t331 = t266 * t280 - t461;
t17 = (t363 - pkin(5) / 0.2e1) * t265 + (t509 + t343) * t295 + t346 + t478;
t183 = -t259 * t293 + t295 * t337;
t330 = -t17 * qJD(1) + t183 * qJD(3);
t182 = t259 * t295 + t457;
t299 = (-t364 + t510) * t295 + (t365 + t343) * t293;
t21 = t299 + t523;
t329 = -t21 * qJD(1) + t182 * qJD(3);
t328 = t265 * t351;
t326 = -t161 * qJD(4) + t515;
t325 = qJD(5) * t337 - t396;
t324 = -t76 * qJ(6) / 0.2e1 + t77 * pkin(5) / 0.2e1;
t29 = t381 + t382;
t318 = -t29 * qJD(1) + t259 * t426;
t221 = t295 * t361;
t37 = t221 + (t455 / 0.2e1 + t509) * t295;
t317 = -t37 * qJD(1) - t280 * t426;
t301 = (t362 - t455 / 0.2e1) * t293 + t522 / 0.2e1;
t35 = t301 + t438;
t316 = -t35 * qJD(1) - t280 * t425;
t315 = t295 * t328;
t314 = t260 * qJD(5) + t374;
t313 = t166 * qJD(5) + t266 * t426;
t173 = t274 * t263;
t310 = t173 * qJD(1) + 0.2e1 * t349;
t309 = -t274 * qJD(3) + 0.2e1 * t350;
t104 = t508 + (t360 + t505) * t265;
t298 = -t511 * t279 + t157 * t259 / 0.2e1 + t94 * t505;
t2 = t298 + t324;
t306 = -qJD(3) * t259 * t337 - t2 * qJD(1) + t104 * qJD(2);
t184 = t292 * t263 + t264;
t302 = t184 * qJD(1) + t349 - t399;
t251 = t260 * qJD(3);
t246 = t265 * t426;
t226 = t265 * t368;
t224 = t351 * qJ(6);
t213 = t291 * qJD(3) + t350;
t150 = t163 * qJD(5);
t143 = t161 * qJD(5);
t138 = t158 * qJD(5);
t129 = t496 / 0.2e1 + t305 * pkin(3);
t128 = -t148 + t421;
t118 = t260 + t433;
t113 = t288 + t348;
t106 = 0.1e1 / 0.2e1 - t345;
t105 = qJ(6) * t359 - t364 + t508;
t99 = t295 * t521 + t501;
t98 = t289 + t347;
t93 = t163 * qJD(3) - t247 * t265;
t87 = t89 * qJD(5);
t75 = t293 * t328;
t58 = t311 + t357;
t38 = t200 * t288 + t280 * t359 + t221 + t436;
t36 = t301 - t438;
t28 = -t451 / 0.2e1 + t346 + t382;
t24 = t304 - t320;
t20 = t299 - t523;
t18 = t265 * t363 + t295 * t343 + t254 + t381 + t478;
t14 = t88 * t502 - t492 / 0.2e1 + t89 * t288 + t489 / 0.2e1 + t322 * t266;
t13 = t88 * t289 + t89 * t502 + (t360 - t500 / 0.2e1) * t266 + t321;
t6 = t300 - t322;
t3 = t4 * qJD(3);
t1 = t298 - t324;
t32 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t294 * t391, t273 * qJD(3), 0, 0, 0, qJ(2) * t391 + qJD(2) * t294, -qJ(2) * t394 + qJD(2) * t296, qJD(4) * t517, t281 * qJD(2) + t41 * qJD(3) + t71 * qJD(4), -t263 * t372 + t292 * t373, -t173 * qJD(5) - 0.2e1 * t266 * t349, t121 * qJD(3) + t265 * t378, -t119 * qJD(3) + t265 * t376, -t373, t22 * qJD(3) + t40 * qJD(5) - t295 * t427 + t527, t23 * qJD(3) + t39 * qJD(5) + t526, t16 * qJD(3) + t527 + t26 * qJD(5) + (-t263 * t396 - t427) * t295, -t174 * qJD(2) - t9 * qJD(3) - t10 * qJD(5) + t266 * t370, t15 * qJD(3) + t27 * qJD(5) + t184 * qJD(6) - t526, qJD(2) * t31 + qJD(3) * t5 + qJD(4) * t19 + qJD(5) * t8 + qJD(6) * t30; 0, 0, 0, 0, qJD(1), t430, 0, 0, 0, 0, 0, t395, t392, 0, t106 * qJD(4) + t397, 0, 0, 0, 0, 0, -t115 * qJD(5) - t249, t98 * qJD(5) + t247, t99 * qJD(5) - t249, -t406, t113 * qJD(5) - t247, t58 * qJD(4) + t6 * qJD(5) + t115 * qJD(6) + t3 + t449; 0, 0, 0, 0, 0, 0, -t366, t398, -t394, -t391, 0, -t297 * t394 + t379, -t297 * t391 - t380, -t518 (-t200 * t476 - t477 * t516) * t495 + t129 * qJD(4) + t333, -t138 + (t292 * t429 + t275) * t266 (-t252 + t253) * qJD(3) + (-qJD(5) - t400) * t166 * t519, t246 + t414, t248 - t416, -t314, t466 + (t293 * t331 - t467) * qJD(3) + t38 * qJD(5), t465 + (t295 * t331 + t468) * qJD(3) + t36 * qJD(5), t471 + (t293 * t332 - t479) * qJD(3) + t18 * qJD(5) + t118 * qJD(6), qJD(3) * t335 + t14 * qJD(5) - t483, t472 + (-t295 * t332 - t480) * qJD(3) + t20 * qJD(5) + t226 (t95 * t259 + t279 * t335) * qJD(3) + t24 * qJD(4) + t1 * qJD(5) + t28 * qJD(6) + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t402, t106 * qJD(2) + t129 * qJD(3) + t443, 0, 0, 0, 0, 0, t528, t529, t528, 0, -t529, t58 * qJD(2) + t24 * qJD(3) + t13 * qJD(5) + t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t352, -t310, t75, t315, -t251, t38 * qJD(3) - t417 + t445 - t87, t98 * qJD(2) + t36 * qJD(3) + t442 + t446, t99 * qJD(2) + t18 * qJD(3) + t463 - t87, t14 * qJD(3) + qJD(5) * t344 - t370 - t475, t113 * qJD(2) + t20 * qJD(3) - t257 - t442 + t458, t486 + t6 * qJD(2) + t1 * qJD(3) + t13 * qJD(4) + (-pkin(5) * t89 - t88 * qJ(6)) * qJD(5) + t72 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * qJD(3) - t225, t75, t302, t28 * qJD(3) + t72 * qJD(5) + t417 + t450; 0, 0, 0, 0, -qJD(1), -t430, 0, 0, 0, 0, 0, -t395, -t392, 0, -qJD(4) * t327 - t397, 0, 0, 0, 0, 0, t513, -t96 * qJD(5) - t247, t513, t406, -t112 * qJD(5) + t247, t59 * qJD(4) + t7 * qJD(5) + qJD(6) * t312 + t3 - t449; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t394, -t391, 0, t518, 0, 0, 0, 0, 0, t342, -t313, t342, t174 * qJD(3), t313 (t174 * t279 - t459) * qJD(3) + t105 * qJD(5) + t165 * qJD(6) + t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t420, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t520, -t441 - t339, -t520, 0, t339 - t419, t105 * qJD(3) + t266 * t512 + t494; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t520; 0, 0, 0, 0, 0, 0, t366, -t398, 0, 0, 0, -t379, t380, 0, t130 * qJD(4) - t333, -t292 * t374 - t138, t315 * t519, -t150 - t414, t143 + t416, t314, t37 * qJD(5) - t265 * t423 - t466, t165 * qJD(4) + t35 * qJD(5) - t465, -t166 * qJD(4) + t17 * qJD(5) + t117 * qJD(6) - t471, -t175 * qJD(4) - t11 * qJD(5) - t163 * qJD(6) + t483, t21 * qJD(5) - t265 * t424 + t226 - t472, qJD(4) * t25 + qJD(5) * t2 + qJD(6) * t29 - t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * qJD(5) - t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, t274 * qJD(5), 0, 0, 0, t280 * t422, t280 * t421, -t183 * qJD(5) + t368, 0, -t182 * qJD(5) + t291 * qJD(6), t325 * t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t412, 0, 0, 0, 0, 0, -t367, t410, -t408, -t405, -t371, t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, -t309, t128, t127, t401, -t317 - t375, -t316 + t377, -t330 - t375, -t512 - t474, -t329 - t377, -t279 * t512 - t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, t128, t213, -t318 + t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t402, qJD(2) * t327 - t130 * qJD(3) - t443, 0, 0, 0, 0, 0, t143 + t248 - t528, t340 - t529, t339 - t528, t175 * qJD(3), -t150 + t246 + t529, -t59 * qJD(2) - t25 * qJD(3) - t12 * qJD(5) - t161 * qJD(6) - t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t420, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t412, 0, 0, 0, 0, 0, t367, -t410, t408, t405, t371, -t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t249 - t421, t247 - t422, 0, t128, -t325 - t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t352, t310, t93, -t161 * qJD(3) - t249 * t265, -t251, -t37 * qJD(3) + t326 - t445, t96 * qJD(2) - t35 * qJD(3) - t266 * t423 - t446, -t17 * qJD(3) - t266 * t424 - t463 + t515, qJD(3) * t11 + t475, t112 * qJD(2) - t21 * qJD(3) + t163 * qJD(4) - t257 - t458, -qJ(6) * t257 - t7 * qJD(2) - t2 * qJD(3) + t12 * qJD(4) - t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t514, t441, t514, 0, t419, t104 * qJD(3) - t494; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, t309, t148, -t411, -t401, t317, t316, t330, t474, t329, t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t411, -t249, -t247, 0, t148, t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t351, -t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 * qJD(3) + t225, t93, -t302, qJ(6) * t399 - t29 * qJD(3) - t326 - t450; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t514; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t354, t148, -t213, t318; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t351, t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t32;
