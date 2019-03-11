% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:23
% EndTime: 2019-03-08 21:49:30
% DurationCPUTime: 3.87s
% Computational Cost: add. (2490->401), mult. (5789->535), div. (0->0), fcn. (3731->8), ass. (0->178)
t471 = pkin(4) + pkin(8);
t363 = cos(qJ(2));
t356 = sin(pkin(6));
t443 = qJD(2) * t356;
t407 = qJD(1) * t443;
t394 = t363 * t407;
t357 = cos(pkin(6));
t440 = qJD(3) * t357;
t480 = qJD(1) * t440 + t394;
t360 = sin(qJ(2));
t446 = qJD(1) * t356;
t415 = t360 * t446;
t327 = qJD(2) * pkin(8) + t415;
t359 = sin(qJ(3));
t362 = cos(qJ(3));
t445 = qJD(1) * t357;
t299 = t359 * t327 - t362 * t445;
t479 = qJD(4) + t299;
t442 = qJD(2) * t359;
t346 = qJD(5) + t442;
t478 = t346 ^ 2;
t354 = t359 ^ 2;
t355 = t362 ^ 2;
t477 = (t354 - t355) * MDP(6);
t438 = qJD(3) * t359;
t348 = pkin(3) * t438;
t388 = pkin(9) * t359 - qJ(4) * t362;
t435 = qJD(4) * t359;
t368 = qJD(3) * t388 - t435;
t302 = t348 + t368;
t364 = -pkin(3) - pkin(9);
t402 = -qJ(4) * t359 - pkin(2);
t315 = t362 * t364 + t402;
t436 = qJD(3) * t362;
t324 = t471 * t436;
t333 = t471 * t359;
t358 = sin(qJ(5));
t361 = cos(qJ(5));
t432 = qJD(5) * t361;
t433 = qJD(5) * t358;
t458 = t358 * t363;
t476 = (t359 * t458 + t360 * t361) * t446 - t361 * t302 + t315 * t433 - t358 * t324 - t333 * t432;
t475 = t361 * t315 + t358 * t333;
t426 = qJD(2) * qJD(3);
t405 = t359 * t426;
t439 = qJD(3) * t358;
t441 = qJD(2) * t362;
t317 = t361 * t441 + t439;
t434 = qJD(5) * t317;
t294 = -t358 * t405 + t434;
t300 = t362 * t327 + t359 * t445;
t353 = qJD(3) * qJ(4);
t292 = -t353 - t300;
t291 = -qJD(3) * pkin(3) + t479;
t429 = pkin(4) * t442 + t479;
t423 = MDP(21) + MDP(23);
t422 = -MDP(22) + MDP(25);
t376 = -qJ(4) * t436 - t435;
t448 = pkin(3) * t405 + t360 * t407;
t281 = qJD(2) * t376 + t448;
t307 = t348 + t376;
t365 = qJD(3) ^ 2;
t473 = qJD(2) * (-t307 + t415) - pkin(8) * t365 - t281;
t455 = t361 * t363;
t472 = -qJD(5) * t475 + (t358 * t360 - t359 * t455) * t446 - t358 * t302 + t361 * t324;
t470 = qJD(2) * pkin(2);
t352 = qJD(3) * qJD(4);
t416 = -t327 * t438 + t362 * t480;
t267 = -t352 - t416;
t260 = -pkin(4) * t405 - t267;
t408 = t358 * t441;
t295 = qJD(3) * t432 - qJD(5) * t408 - t361 * t405;
t437 = qJD(3) * t361;
t319 = -t408 + t437;
t247 = pkin(5) * t295 + qJ(6) * t294 - qJD(6) * t319 + t260;
t469 = t247 * t358;
t468 = t247 * t361;
t467 = t260 * t358;
t466 = t260 * t361;
t465 = t294 * t361;
t409 = t363 * t446;
t331 = -pkin(3) * t362 + t402;
t444 = qJD(2) * t331;
t301 = -t409 + t444;
t464 = t301 * t360;
t463 = t317 * t346;
t462 = t319 * t362;
t461 = t346 * t359;
t460 = t346 * t364;
t459 = t356 * t360;
t457 = t359 * t361;
t456 = t359 * t365;
t454 = t362 * t365;
t453 = qJ(6) * t436 + qJD(6) * t359 - t476;
t452 = -pkin(5) * t436 - t472;
t390 = pkin(5) * t361 + qJ(6) * t358;
t379 = -pkin(4) - t390;
t451 = -qJD(5) * t390 + qJD(6) * t361 + t379 * t442 - t479;
t288 = pkin(4) * t441 + t300;
t349 = pkin(3) * t442;
t305 = qJD(2) * t388 + t349;
t450 = t358 * t288 + t361 * t305;
t334 = t471 * t362;
t431 = qJD(6) * t346;
t279 = t353 + t288;
t430 = t279 * qJD(5);
t275 = qJD(3) * t364 + t429;
t289 = qJD(2) * t315 - t409;
t254 = t275 * t361 - t289 * t358;
t428 = qJD(6) - t254;
t425 = -MDP(10) + MDP(13);
t424 = MDP(11) - MDP(14);
t421 = t358 * t460;
t420 = t361 * t460;
t419 = t359 * t459;
t366 = qJD(2) ^ 2;
t418 = t359 * t362 * t366;
t273 = t327 * t436 + t359 * t480;
t404 = t362 * t426;
t265 = pkin(4) * t404 + t273;
t274 = qJD(2) * t368 + t448;
t417 = -t358 * t265 - t361 * t274 - t275 * t432;
t414 = t360 * t443;
t413 = t363 * t443;
t412 = t364 * t436;
t411 = t346 * t433;
t410 = t362 * t432;
t403 = MDP(20) * t441;
t399 = pkin(5) * t404;
t398 = t361 * t265 - t358 * t274 - t275 * t433 - t289 * t432;
t397 = t317 * t409;
t396 = t319 * t409;
t395 = qJ(6) * t404;
t339 = t361 * t404;
t246 = -t398 - t399;
t255 = t275 * t358 + t289 * t361;
t250 = qJ(6) * t346 + t255;
t391 = -t250 * t442 + t246;
t389 = -pkin(5) * t358 + qJ(6) * t361;
t249 = -pkin(5) * t346 + t428;
t387 = t249 * t358 + t250 * t361;
t386 = t288 * t361 - t305 * t358;
t383 = -t315 * t358 + t333 * t361;
t382 = -qJD(2) * t355 + t461;
t381 = t346 * t319;
t378 = qJD(3) * t299 + t416;
t377 = qJD(3) * t300 - t273;
t310 = -t357 * t362 + t419;
t286 = t310 * t358 - t356 * t455;
t285 = t310 * t361 + t356 * t458;
t311 = t357 * t359 + t362 * t459;
t375 = t255 * t346 + t398;
t374 = t289 * t433 + t417;
t328 = -t409 - t470;
t371 = qJD(3) * (t328 + t409 - t470);
t370 = qJD(3) * (-t301 - t409 - t444);
t369 = t254 * t346 + t374;
t367 = -t267 * t362 + t273 * t359 + (t291 * t362 + t292 * t359) * qJD(3);
t330 = qJ(4) - t389;
t329 = t364 * t339;
t323 = t471 * t438;
t322 = -qJ(4) * t441 + t349;
t298 = t362 * t390 + t334;
t290 = t301 * t442;
t284 = qJD(3) * t311 + t359 * t413;
t283 = -qJD(3) * t419 + (t413 + t440) * t362;
t282 = pkin(5) * t319 + qJ(6) * t317;
t277 = -pkin(5) * t359 - t383;
t276 = qJ(6) * t359 + t475;
t264 = t463 - t294;
t259 = (qJD(5) * t389 + qJD(6) * t358) * t362 + (-pkin(8) + t379) * t438;
t258 = -pkin(5) * t441 - t386;
t257 = qJ(6) * t441 + t450;
t256 = pkin(5) * t317 - qJ(6) * t319 + t279;
t253 = qJD(5) * t285 + t284 * t358 + t361 * t414;
t252 = qJD(5) * t286 - t284 * t361 + t358 * t414;
t245 = -t374 + t395 + t431;
t1 = [(-t267 * t311 + t273 * t310 - t283 * t292 + t284 * t291) * MDP(15) + (t252 * t319 - t253 * t317 + t285 * t294 - t286 * t295) * MDP(24) + (t245 * t286 - t246 * t285 + t247 * t311 + t249 * t252 + t250 * t253 + t256 * t283) * MDP(26) + (t283 * t362 + t284 * t359) * MDP(12) * qJD(2) + (t425 * t284 - t424 * t283 + (-t311 * t359 * MDP(12) + (t310 * MDP(12) + t285 * t423 + t286 * t422) * t362) * qJD(2)) * qJD(3) + (-t281 * t363 * MDP(15) + (MDP(15) * t464 + (t359 * t425 - t362 * t424) * t363 * qJD(3)) * qJD(2) + (-t363 * MDP(4) + (t359 * t424 + t362 * t425 - MDP(3)) * t360) * t366) * t356 + t423 * (-t252 * t346 + t283 * t317 + t311 * t295) + t422 * (t253 * t346 - t283 * t319 + t294 * t311); 0.2e1 * t359 * MDP(5) * t404 - 0.2e1 * t426 * t477 + MDP(7) * t454 - MDP(8) * t456 + (-pkin(8) * t454 + t359 * t371) * MDP(10) + (pkin(8) * t456 + t362 * t371) * MDP(11) + ((-t354 - t355) * t394 + t367) * MDP(12) + (t359 * t370 - t473 * t362) * MDP(13) + (t473 * t359 + t362 * t370) * MDP(14) + (t281 * t331 + t301 * t307 + (-t464 + (-t291 * t359 + t292 * t362) * t363) * t446 + t367 * pkin(8)) * MDP(15) + (t294 * t358 * t362 + (t358 * t438 - t410) * t319) * MDP(16) + ((-t317 * t358 + t319 * t361) * t438 + (t465 + t295 * t358 + (t317 * t361 + t319 * t358) * qJD(5)) * t362) * MDP(17) + (-t346 * t410 - t294 * t359 + (t358 * t382 + t462) * qJD(3)) * MDP(18) + (t362 * t411 - t295 * t359 + (-t317 * t362 + t361 * t382) * qJD(3)) * MDP(19) + (t346 + t442) * MDP(20) * t436 + (t334 * t295 - t323 * t317 + (-t279 * t437 + t398) * t359 + t472 * t346 + (-t397 - t358 * t430 + t466 + (qJD(2) * t383 + t254) * qJD(3)) * t362) * MDP(21) + (-t334 * t294 - t323 * t319 + ((qJD(3) * t279 + qJD(5) * t289) * t358 + t417) * t359 + t476 * t346 + (-t396 - t361 * t430 - t467 + (-qJD(2) * t475 - t255) * qJD(3)) * t362) * MDP(22) + (t259 * t317 + t295 * t298 + (-t256 * t437 - t246) * t359 - t452 * t346 + (-t397 - t256 * t433 + t468 + (-qJD(2) * t277 - t249) * qJD(3)) * t362) * MDP(23) + (-t276 * t295 - t277 * t294 + t452 * t319 - t453 * t317 + t387 * t438 + (-t245 * t361 - t246 * t358 + (-t249 * t361 + t250 * t358) * qJD(5)) * t362) * MDP(24) + (-t259 * t319 + t294 * t298 + (-t256 * t439 + t245) * t359 + t453 * t346 + (t396 + t256 * t432 + t469 + (qJD(2) * t276 + t250) * qJD(3)) * t362) * MDP(25) + (t245 * t276 + t246 * t277 + t247 * t298 + (-t362 * t409 + t259) * t256 + t453 * t250 + t452 * t249) * MDP(26); -MDP(5) * t418 + t366 * t477 + (-t328 * t442 + t377) * MDP(10) + (-t328 * t441 - t378) * MDP(11) + (-t322 * t441 + t290 - t377) * MDP(13) + (0.2e1 * t352 + (t301 * t362 + t322 * t359) * qJD(2) + t378) * MDP(14) + (-pkin(3) * t273 - qJ(4) * t267 - t291 * t300 - t292 * t479 - t301 * t322) * MDP(15) + (-t358 * t381 - t465) * MDP(16) + ((-t295 - t381) * t361 + (t294 + t463) * t358) * MDP(17) + (-t411 + t339 + (-t358 * t461 - t462) * qJD(2)) * MDP(18) + (-t346 * t432 + (-t346 * t457 + (t317 - t439) * t362) * qJD(2)) * MDP(19) - t346 * t403 + (t329 + qJ(4) * t295 + t467 - t386 * t346 + t429 * t317 + (t279 * t361 - t421) * qJD(5) + (-t254 * t362 + t279 * t457) * qJD(2)) * MDP(21) + (-qJ(4) * t294 + t466 + t450 * t346 + t429 * t319 + (-t279 * t358 - t420) * qJD(5) + (t255 * t362 + (-t279 * t359 - t412) * t358) * qJD(2)) * MDP(22) + (t469 + t258 * t346 + t295 * t330 + t329 - t451 * t317 + (t256 * t361 - t421) * qJD(5) + (t249 * t362 + t256 * t457) * qJD(2)) * MDP(23) + (t257 * t317 - t258 * t319 + (t294 * t364 + (-t317 * t364 - t250) * qJD(5) + t391) * t361 + (-t249 * t442 - t295 * t364 - t245 + (t319 * t364 - t249) * qJD(5)) * t358) * MDP(24) + (-t468 - t257 * t346 + t294 * t330 + t451 * t319 + (t256 * t358 + t420) * qJD(5) + (-t250 * t362 + (t256 * t359 + t412) * t358) * qJD(2)) * MDP(25) + (t247 * t330 - t249 * t258 - t250 * t257 - t451 * t256 + (qJD(5) * t387 + t245 * t358 - t246 * t361) * t364) * MDP(26); MDP(13) * t418 + (-t354 * t366 - t365) * MDP(14) + (t290 + t273) * MDP(15) + t423 * t339 + (t292 * MDP(15) - t256 * MDP(26) - t317 * t423 + t319 * t422) * qJD(3) + ((-t317 * t442 + t294 - t434) * MDP(24) + (qJD(5) * t250 - t391) * MDP(26) + t422 * t478) * t361 + ((qJD(5) * t319 - t295) * MDP(24) + (qJD(5) * t249 + t245) * MDP(26) + ((MDP(24) * t319 + MDP(26) * t249) * t359 + t422 * t436) * qJD(2) - t423 * t478) * t358; t264 * MDP(18) - t295 * MDP(19) + qJD(3) * t403 + t375 * MDP(21) + t369 * MDP(22) + (t375 + 0.2e1 * t399) * MDP(23) + (pkin(5) * t294 - qJ(6) * t295) * MDP(24) + (-t369 + 0.2e1 * t395 + 0.2e1 * t431) * MDP(25) + (-pkin(5) * t246 + qJ(6) * t245 - t249 * t255 + t250 * t428 - t256 * t282) * MDP(26) + (t346 * MDP(19) - t279 * MDP(21) - t256 * MDP(23) + (t250 - t255) * MDP(24) + t282 * MDP(25) + MDP(17) * t319) * t319 + (t319 * MDP(16) + t279 * MDP(22) - t282 * MDP(23) + (t249 - t428) * MDP(24) - t256 * MDP(25) - MDP(17) * t317) * t317; (t317 * t319 - t404) * MDP(23) + t264 * MDP(24) + (-t319 ^ 2 - t478) * MDP(25) + (-t250 * t346 + t256 * t319 + t246) * MDP(26);];
tauc  = t1;
