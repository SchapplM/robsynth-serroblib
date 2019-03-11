% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:07
% EndTime: 2019-03-09 05:21:13
% DurationCPUTime: 3.24s
% Computational Cost: add. (4915->344), mult. (10778->468), div. (0->0), fcn. (7651->8), ass. (0->171)
t392 = sin(qJ(4));
t393 = sin(qJ(3));
t395 = cos(qJ(4));
t396 = cos(qJ(3));
t359 = t392 * t396 + t393 * t395;
t354 = t359 * qJD(1);
t385 = qJD(3) + qJD(4);
t464 = t395 * t396;
t419 = t392 * t393 - t464;
t489 = qJD(1) * t419;
t403 = t385 * t489;
t397 = -pkin(1) - pkin(7);
t368 = qJD(1) * t397 + qJD(2);
t441 = pkin(8) * qJD(1) - t368;
t453 = qJD(3) * t393;
t335 = t441 * t453;
t455 = qJD(1) * t393;
t343 = -pkin(8) * t455 + t368 * t393;
t451 = qJD(4) * t392;
t439 = t392 * t335 - t343 * t451;
t454 = qJD(1) * t396;
t344 = -pkin(8) * t454 + t396 * t368;
t334 = qJD(3) * pkin(3) + t344;
t452 = qJD(3) * t396;
t336 = t441 * t452;
t490 = t395 * (qJD(4) * t334 - t336);
t267 = qJ(5) * t403 - t354 * qJD(5) + t439 + t490;
t390 = sin(pkin(10));
t326 = t385 * t359;
t319 = t326 * qJD(1);
t353 = t392 * t455 - t395 * t454;
t421 = t395 * t335 + t392 * t336;
t332 = t395 * t343;
t422 = -t334 * t392 - t332;
t406 = qJD(4) * t422 + t421;
t401 = qJ(5) * t319 + qJD(5) * t353 + t406;
t482 = cos(pkin(10));
t254 = t267 * t390 - t401 * t482;
t322 = t359 * t390 + t419 * t482;
t501 = t254 * t322;
t412 = -t353 * t482 - t390 * t354;
t391 = sin(qJ(6));
t450 = qJD(6) * t391;
t292 = -t319 * t482 + t390 * t403;
t394 = cos(qJ(6));
t449 = qJD(6) * t394;
t461 = t394 * t292 + t385 * t449;
t274 = -t412 * t450 + t461;
t500 = t274 * t322;
t472 = t412 * t391;
t308 = -t394 * t385 + t472;
t437 = t353 * t390 - t482 * t354;
t488 = qJD(6) - t437;
t499 = t308 * t488;
t310 = t385 * t391 + t394 * t412;
t498 = t310 * t488;
t497 = qJ(2) * MDP(6) + MDP(5);
t436 = t488 * t394;
t291 = -t319 * t390 - t403 * t482;
t465 = t391 * t291;
t496 = -t488 * t436 - t465;
t255 = t267 * t482 + t390 * t401;
t347 = t353 * qJ(5);
t331 = t392 * t343;
t440 = t395 * t334 - t331;
t302 = t347 + t440;
t298 = pkin(4) * t385 + t302;
t480 = qJ(5) * t354;
t303 = -t422 - t480;
t467 = t390 * t303;
t272 = t298 * t482 - t467;
t299 = t482 * t303;
t273 = t390 * t298 + t299;
t327 = t385 * t464 - t392 * t453 - t393 * t451;
t295 = t482 * t326 + t327 * t390;
t411 = t359 * t482 - t390 * t419;
t413 = -t390 * t326 + t327 * t482;
t495 = t255 * t411 - t272 * t295 + t273 * t413;
t494 = MDP(7) * t393;
t493 = MDP(8) * (t393 ^ 2 - t396 ^ 2);
t492 = t291 * t411;
t386 = qJD(1) * qJD(2);
t447 = qJD(1) * qJD(3);
t444 = t396 * t447;
t362 = pkin(3) * t444 + t386;
t400 = -pkin(4) * t403 + t362;
t484 = pkin(8) - t397;
t357 = t484 * t453;
t364 = t484 * t396;
t358 = qJD(3) * t364;
t363 = t484 * t393;
t469 = t364 * t395;
t410 = -qJD(4) * t469 + t392 * t357 - t395 * t358 + t363 * t451;
t281 = -qJ(5) * t327 - qJD(5) * t359 + t410;
t420 = t363 * t395 + t364 * t392;
t407 = qJD(4) * t420 + t395 * t357 + t358 * t392;
t404 = qJ(5) * t326 + qJD(5) * t419 + t407;
t262 = t281 * t482 + t390 * t404;
t270 = -t385 * pkin(5) - t272;
t365 = pkin(3) * t455 + qJD(1) * qJ(2);
t328 = pkin(4) * t354 + qJD(5) + t365;
t282 = -pkin(5) * t437 - pkin(9) * t412 + t328;
t378 = t393 * pkin(3) + qJ(2);
t430 = pkin(4) * t359 + t378;
t287 = pkin(5) * t411 + pkin(9) * t322 + t430;
t312 = -qJ(5) * t359 - t420;
t409 = qJ(5) * t419 + t363 * t392 - t469;
t286 = t312 * t482 + t390 * t409;
t427 = -t286 * t291 - t501;
t486 = -t270 * t295 - (qJD(6) * t287 + t262) * t488 - t411 * (qJD(6) * t282 + t255) + t427;
t485 = pkin(4) * t353;
t483 = pkin(3) * qJD(4);
t479 = t270 * t437;
t478 = t270 * t322;
t477 = t274 * t391;
t476 = t287 * t291;
t475 = t292 * t391;
t474 = t308 * t412;
t473 = t310 * t412;
t471 = t326 * t385;
t470 = t327 * t385;
t468 = t365 * t353;
t466 = t390 * t392;
t288 = t394 * t291;
t399 = qJD(1) ^ 2;
t463 = t396 * t399;
t398 = qJD(3) ^ 2;
t462 = t397 * t398;
t460 = t395 * t344 - t331;
t304 = t347 + t460;
t438 = -t344 * t392 - t332;
t417 = t438 + t480;
t442 = t482 * t392;
t459 = -t304 * t390 + t482 * t417 + (t390 * t395 + t442) * t483;
t458 = -t482 * t304 - t390 * t417 + (t395 * t482 - t466) * t483;
t379 = pkin(3) * t395 + pkin(4);
t346 = pkin(3) * t442 + t390 * t379;
t369 = pkin(3) * t452 + qJD(2);
t446 = 0.2e1 * qJD(1);
t382 = pkin(3) * t454;
t443 = -pkin(3) * t385 - t334;
t284 = pkin(5) * t412 - pkin(9) * t437 - t485;
t340 = pkin(9) + t346;
t432 = qJD(6) * t340 + t284 + t382;
t376 = pkin(4) * t390 + pkin(9);
t431 = qJD(6) * t376 + t284;
t429 = pkin(4) * t327 + t369;
t271 = pkin(9) * t385 + t273;
t260 = t271 * t394 + t282 * t391;
t428 = t254 * t391 + t260 * t412 + t270 * t449;
t426 = -t291 * t340 - t479;
t425 = -t291 * t376 - t479;
t424 = t271 * t391 - t282 * t394;
t423 = t272 * t437 + t273 * t412;
t418 = t288 + (t391 * t437 - t450) * t488;
t416 = -t254 * t394 + t270 * t450 + t412 * t424;
t415 = t365 * t354 - t439;
t414 = -t295 * t394 + t322 * t450;
t345 = -pkin(3) * t466 + t379 * t482;
t275 = qJD(6) * t310 + t475;
t405 = -t353 * t354 * MDP(14) - t488 * t412 * MDP(27) + ((t274 - t499) * t394 + (-t275 - t498) * t391) * MDP(24) + (t418 + t474) * MDP(26) + (-t473 - t496) * MDP(25) + (t310 * t436 + t477) * MDP(23) + (t354 * t385 - t319) * MDP(16) + (-t353 * t385 + t403) * MDP(17) + (t353 ^ 2 - t354 ^ 2) * MDP(15);
t377 = -pkin(4) * t482 - pkin(5);
t339 = -pkin(5) - t345;
t285 = t312 * t390 - t409 * t482;
t277 = t302 * t482 - t467;
t276 = t302 * t390 + t299;
t265 = pkin(5) * t413 + pkin(9) * t295 + t429;
t264 = t291 * pkin(5) - t292 * pkin(9) + t400;
t263 = t394 * t264;
t261 = t281 * t390 - t404 * t482;
t1 = [-0.2e1 * t444 * t494 + 0.2e1 * t447 * t493 + (-t393 * t462 + (qJ(2) * t452 + qJD(2) * t393) * t446) * MDP(12) + (-t396 * t462 + (-qJ(2) * t453 + qJD(2) * t396) * t446) * MDP(13) + (t319 * t419 + t326 * t353) * MDP(14) + (t319 * t359 + t326 * t354 + t353 * t327 - t403 * t419) * MDP(15) - MDP(16) * t471 - MDP(17) * t470 + (t365 * t327 + t369 * t354 + t362 * t359 + (-t378 * t489 + t407) * t385) * MDP(19) + (-t378 * t319 - t365 * t326 - t369 * t353 - t362 * t419 - t385 * t410) * MDP(20) + (t261 * t412 + t262 * t437 + t285 * t292 + t427 - t495) * MDP(21) + (t254 * t285 + t255 * t286 - t272 * t261 + t273 * t262 + t328 * t429 + t400 * t430) * MDP(22) + (t310 * t414 - t394 * t500) * MDP(23) + (-(-t308 * t394 - t310 * t391) * t295 - (-t477 - t275 * t394 + (t308 * t391 - t310 * t394) * qJD(6)) * t322) * MDP(24) + (t274 * t411 - t288 * t322 + t310 * t413 + t414 * t488) * MDP(25) + (t322 * t465 - t275 * t411 - t413 * t308 + (t295 * t391 + t322 * t449) * t488) * MDP(26) + (t413 * t488 + t492) * MDP(27) + (-t424 * t413 + t261 * t308 + t263 * t411 + t285 * t275 + (t265 * t488 + t476 + (-t271 * t411 - t286 * t488 - t478) * qJD(6)) * t394 + t486 * t391) * MDP(28) + (-t260 * t413 + t261 * t310 + t285 * t274 + (-(-qJD(6) * t286 + t265) * t488 - t476 - (-qJD(6) * t271 + t264) * t411 + qJD(6) * t478) * t391 + t486 * t394) * MDP(29) + 0.2e1 * t497 * t386 + (-t396 * MDP(10) - t393 * MDP(9)) * t398; (-qJD(1) * t354 - t471) * MDP(19) + (qJD(1) * t353 - t470) * MDP(20) + (t292 * t322 + t295 * t412 + t413 * t437 - t492) * MDP(21) + (-qJD(1) * t328 + t495 + t501) * MDP(22) + (t275 * t322 + t295 * t308 - t411 * t465) * MDP(28) + (-t288 * t411 + t295 * t310 + t500) * MDP(29) - t497 * t399 + ((-qJD(1) * t394 - t391 * t413 - t411 * t449) * MDP(28) + (qJD(1) * t391 - t394 * t413 + t411 * t450) * MDP(29)) * t488 + (MDP(12) * t393 + MDP(13) * t396) * (-t398 - t399); (-t354 * t382 + t468 - t438 * t385 + (t392 * t443 - t332) * qJD(4) + t421) * MDP(19) + (t353 * t382 + t460 * t385 + (qJD(4) * t443 + t336) * t395 + t415) * MDP(20) + (-t291 * t346 - t292 * t345 + t412 * t459 + t437 * t458 + t423) * MDP(21) - t399 * t493 + t463 * t494 + (t255 * t346 - t254 * t345 - t328 * (t382 - t485) + t458 * t273 - t459 * t272) * MDP(22) + (t339 * t275 + t426 * t391 + t459 * t308 + (-t391 * t458 - t394 * t432) * t488 + t416) * MDP(28) + (t339 * t274 + t426 * t394 + t459 * t310 + (t391 * t432 - t394 * t458) * t488 + t428) * MDP(29) + t405 + (MDP(13) * t393 * t399 - MDP(12) * t463) * qJ(2); (t385 * t440 + t415 - t490) * MDP(20) + (t377 * t274 - t276 * t310 + t425 * t394 + (t394 * t277 + t391 * t431) * t488 + t428) * MDP(29) + (-t276 * t412 - t277 * t437 + (-t291 * t390 - t292 * t482) * pkin(4) + t423) * MDP(21) + (-t385 * t422 + t406 + t468) * MDP(19) + (t272 * t276 - t273 * t277 + (-t254 * t482 + t255 * t390 + t328 * t353) * pkin(4)) * MDP(22) + (t377 * t275 - t276 * t308 + t425 * t391 + (t391 * t277 - t394 * t431) * t488 + t416) * MDP(28) + t405; (-t412 ^ 2 - t437 ^ 2) * MDP(21) + (t418 - t474) * MDP(28) + (-t473 + t496) * MDP(29) + (t272 * t412 - t273 * t437 + t400) * MDP(22); t310 * t308 * MDP(23) + (-t308 ^ 2 + t310 ^ 2) * MDP(24) + (t461 + t499) * MDP(25) + (-t475 + t498) * MDP(26) + t291 * MDP(27) + (-t255 * t391 + t260 * t488 - t270 * t310 + t263) * MDP(28) + (-t255 * t394 - t264 * t391 + t270 * t308 - t424 * t488) * MDP(29) + (-MDP(25) * t472 - MDP(26) * t310 - MDP(28) * t260 + MDP(29) * t424) * qJD(6);];
tauc  = t1;
