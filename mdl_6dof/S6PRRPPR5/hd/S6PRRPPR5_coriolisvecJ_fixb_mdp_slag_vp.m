% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:12
% EndTime: 2019-03-08 21:21:21
% DurationCPUTime: 4.22s
% Computational Cost: add. (2147->393), mult. (5329->551), div. (0->0), fcn. (3715->10), ass. (0->190)
t510 = pkin(4) + pkin(8);
t408 = sin(qJ(3));
t482 = qJD(2) * t408;
t391 = qJD(6) + t482;
t404 = cos(pkin(11));
t411 = cos(qJ(3));
t481 = qJD(2) * t411;
t458 = t404 * t481;
t402 = sin(pkin(11));
t480 = qJD(3) * t402;
t358 = t458 + t480;
t456 = t402 * t481;
t479 = qJD(3) * t404;
t360 = -t456 + t479;
t407 = sin(qJ(6));
t410 = cos(qJ(6));
t432 = t358 * t407 - t360 * t410;
t519 = t391 * t432;
t412 = cos(qJ(2));
t403 = sin(pkin(6));
t483 = qJD(2) * t403;
t455 = qJD(1) * t483;
t445 = t412 * t455;
t405 = cos(pkin(6));
t478 = qJD(3) * t405;
t518 = qJD(1) * t478 + t445;
t409 = sin(qJ(2));
t486 = qJD(1) * t403;
t463 = t409 * t486;
t376 = qJD(2) * pkin(8) + t463;
t485 = qJD(1) * t405;
t331 = t408 * t376 - t411 * t485;
t517 = qJD(4) + t331;
t400 = t408 ^ 2;
t401 = t411 ^ 2;
t516 = MDP(6) * (t400 - t401);
t515 = t408 * MDP(5);
t477 = qJD(3) * t408;
t394 = pkin(3) * t477;
t440 = -qJ(4) * t411 + qJ(5) * t408;
t475 = qJD(4) * t408;
t415 = qJD(3) * t440 - qJD(5) * t411 - t475;
t324 = t394 + t415;
t476 = qJD(3) * t411;
t368 = t510 * t476;
t499 = t408 * t412;
t493 = -t324 * t402 + t404 * t368 - (-t402 * t409 + t404 * t499) * t486;
t492 = t404 * t324 + t402 * t368 - (t402 * t499 + t404 * t409) * t486;
t332 = t411 * t376 + t408 * t485;
t399 = qJD(3) * qJ(4);
t327 = -t399 - t332;
t326 = -qJD(3) * pkin(3) + t517;
t362 = t402 * t410 + t404 * t407;
t347 = t362 * qJD(6);
t471 = pkin(4) * t482 + t517;
t513 = qJD(5) + t399;
t393 = pkin(4) * t481;
t323 = t393 + t332;
t309 = t323 + t513;
t406 = -pkin(3) - qJ(5);
t512 = t408 * (-t309 + t513) - t406 * t476;
t425 = -qJ(4) * t476 - t475;
t470 = qJD(2) * qJD(3);
t453 = t408 * t470;
t488 = pkin(3) * t453 + t409 * t455;
t316 = qJD(2) * t425 + t488;
t344 = t394 + t425;
t413 = qJD(3) ^ 2;
t511 = qJD(2) * (-t344 + t463) - pkin(8) * t413 - t316;
t509 = -pkin(9) + t406;
t508 = qJD(2) * pkin(2);
t398 = qJD(3) * qJD(4);
t465 = -t376 * t477 + t518 * t411;
t299 = -t398 - t465;
t288 = -pkin(4) * t453 - t299;
t507 = t288 * t402;
t506 = t288 * t404;
t504 = t360 * t407;
t306 = t410 * t358 + t504;
t505 = t306 * t391;
t503 = t402 * t407;
t502 = t403 * t409;
t501 = t403 * t412;
t500 = t404 * t411;
t498 = t408 * t413;
t497 = t411 * t413;
t414 = qJD(2) ^ 2;
t496 = t411 * t414;
t429 = -pkin(9) * t402 * t408 + pkin(5) * t411;
t420 = t429 * qJD(3);
t495 = t420 + t493;
t446 = pkin(9) * t404 * t477;
t494 = -t446 - t492;
t303 = t376 * t476 + t518 * t408;
t289 = (-qJD(5) + t393) * qJD(3) + t303;
t295 = qJD(2) * t415 + t488;
t271 = t402 * t289 + t404 * t295;
t305 = qJD(3) * t406 + t471;
t450 = -qJ(4) * t408 - pkin(2);
t355 = t406 * t411 + t450;
t462 = t412 * t486;
t321 = qJD(2) * t355 - t462;
t275 = t402 * t305 + t404 * t321;
t431 = -t404 * t410 + t503;
t452 = t411 * t470;
t491 = -t347 * t391 - t431 * t452;
t395 = pkin(3) * t482;
t342 = qJD(2) * t440 + t395;
t284 = t402 * t323 + t404 * t342;
t381 = t510 * t408;
t314 = t404 * t355 + t402 * t381;
t422 = t362 * t408;
t336 = qJD(2) * t422;
t490 = -t347 - t336;
t459 = t404 * t482;
t473 = qJD(6) * t410;
t474 = qJD(6) * t407;
t489 = -t402 * t474 + t404 * t473 + t410 * t459 - t482 * t503;
t382 = t510 * t411;
t378 = -pkin(3) * t411 + t450;
t484 = qJD(2) * t378;
t464 = -pkin(5) * t404 - pkin(4);
t428 = t464 * t482;
t472 = -t428 + t517;
t469 = -MDP(10) + MDP(13);
t468 = MDP(11) - MDP(14);
t467 = t408 * t502;
t443 = t410 * t453;
t444 = t407 * t453;
t466 = -t358 * t473 + t402 * t443 + t404 * t444;
t461 = t409 * t483;
t460 = t412 * t483;
t451 = MDP(24) * t481;
t270 = t404 * t289 - t295 * t402;
t268 = qJD(2) * t420 + t270;
t269 = qJD(2) * t446 + t271;
t447 = t410 * t268 - t269 * t407;
t274 = t404 * t305 - t321 * t402;
t283 = t404 * t323 - t342 * t402;
t439 = t268 * t407 + t269 * t410;
t438 = t270 * t404 + t271 * t402;
t272 = pkin(5) * t482 - pkin(9) * t360 + t274;
t273 = -pkin(9) * t358 + t275;
t265 = t272 * t410 - t273 * t407;
t266 = t272 * t407 + t273 * t410;
t437 = -t274 * t402 + t275 * t404;
t365 = t404 * t381;
t296 = pkin(5) * t408 + t365 + (pkin(9) * t411 - t355) * t402;
t302 = -pkin(9) * t500 + t314;
t436 = t296 * t410 - t302 * t407;
t435 = t296 * t407 + t302 * t410;
t351 = -t405 * t411 + t467;
t317 = t351 * t404 + t402 * t501;
t318 = t351 * t402 - t404 * t501;
t434 = t317 * t410 - t318 * t407;
t433 = t317 * t407 + t318 * t410;
t427 = qJD(3) * t331 + t465;
t426 = qJD(3) * t332 - t303;
t352 = t405 * t408 + t411 * t502;
t369 = t509 * t402;
t424 = qJD(2) * t429 + qJD(5) * t404 + qJD(6) * t369 + t283;
t370 = t509 * t404;
t423 = pkin(9) * t459 + qJD(5) * t402 - qJD(6) * t370 + t284;
t339 = t431 * t411;
t279 = -t360 * t474 + t466;
t421 = t402 * t444 - t404 * t443;
t419 = t327 * MDP(15) - t309 * MDP(19) - t306 * MDP(25);
t377 = -t462 - t508;
t418 = qJD(3) * (t377 + t462 - t508);
t333 = -t462 + t484;
t417 = qJD(3) * (-t333 - t462 - t484);
t416 = -t299 * t411 + t303 * t408 + (t326 * t411 + t327 * t408) * qJD(3);
t280 = -qJD(6) * t432 + t421;
t392 = pkin(5) * t402 + qJ(4);
t367 = t510 * t477;
t366 = -qJ(4) * t481 + t395;
t346 = pkin(5) * t500 + t382;
t340 = t362 * t411;
t337 = (-pkin(8) + t464) * t477;
t325 = t333 * t482;
t320 = qJD(3) * t352 + t408 * t460;
t319 = -qJD(3) * t467 + (t460 + t478) * t411;
t313 = -t355 * t402 + t365;
t298 = t347 * t411 - t431 * t477;
t297 = qJD(3) * t422 + qJD(6) * t339;
t294 = t320 * t402 + t404 * t461;
t293 = t320 * t404 - t402 * t461;
t287 = pkin(5) * t358 + t309;
t281 = qJD(3) * t428 - t299;
t1 = [(-t299 * t352 + t303 * t351 + t320 * t326) * MDP(15) + (-t293 * t360 - t294 * t358) * MDP(18) + (t270 * t317 + t271 * t318 + t274 * t293 + t275 * t294 + t288 * t352) * MDP(19) + ((-qJD(6) * t433 + t293 * t410 - t294 * t407) * t391 + t352 * t280) * MDP(25) + (-(qJD(6) * t434 + t293 * t407 + t294 * t410) * t391 + t352 * t279) * MDP(26) + t469 * t320 * qJD(3) + (-t316 * t412 * MDP(15) + (-MDP(4) * t412 + (t408 * t468 + t411 * t469 - MDP(3)) * t409) * t414) * t403 + (t358 * MDP(16) + t360 * MDP(17) - MDP(26) * t432 - qJD(3) * t468 - t419) * t319 + (t333 * MDP(15) * t502 + t319 * t411 * MDP(12) + (MDP(12) * t320 + MDP(16) * t293 - MDP(17) * t294) * t408 + ((t351 * MDP(12) + t317 * MDP(16) - t318 * MDP(17) + MDP(25) * t434 - MDP(26) * t433) * t411 + ((-t317 * t402 + t318 * t404) * MDP(18) + (-t404 * MDP(16) + t402 * MDP(17) - MDP(12)) * t352) * t408 + (t408 * t469 - t411 * t468) * t501) * qJD(3)) * qJD(2); 0.2e1 * t452 * t515 - 0.2e1 * t470 * t516 + MDP(7) * t497 - MDP(8) * t498 + (-pkin(8) * t497 + t408 * t418) * MDP(10) + (pkin(8) * t498 + t411 * t418) * MDP(11) + ((-t400 - t401) * t445 + t416) * MDP(12) + (t408 * t417 - t411 * t511) * MDP(13) + (t408 * t511 + t411 * t417) * MDP(14) + (t316 * t378 + t333 * t344 + (-t333 * t409 + (-t326 * t408 + t327 * t411) * t412) * t486 + t416 * pkin(8)) * MDP(15) + (-t358 * t367 + (-t358 * t462 + t506 + (qJD(2) * t313 + t274) * qJD(3)) * t411 + (-t309 * t479 + t270 + (-t382 * t479 + t493) * qJD(2)) * t408) * MDP(16) + (-t360 * t367 + (-t360 * t462 - t507 + (-qJD(2) * t314 - t275) * qJD(3)) * t411 + (t309 * t480 - t271 + (t382 * t480 - t492) * qJD(2)) * t408) * MDP(17) + ((t270 * t402 - t271 * t404) * t411 - t493 * t360 - t492 * t358 + ((-t313 * t402 + t314 * t404) * qJD(2) + t437) * t477) * MDP(18) + (t270 * t313 + t271 * t314 + t288 * t382 + (-t411 * t462 - t367) * t309 + t492 * t275 + t493 * t274) * MDP(19) + (-t279 * t340 - t297 * t432) * MDP(20) + (t279 * t339 + t280 * t340 - t297 * t306 - t298 * t432) * MDP(21) + (t279 * t408 + t297 * t391 + (-qJD(2) * t340 - t432) * t476) * MDP(22) + (-t280 * t408 + t298 * t391 + (qJD(2) * t339 - t306) * t476) * MDP(23) + (t391 + t482) * MDP(24) * t476 + (t447 * t408 + t337 * t306 + t346 * t280 - t281 * t339 - t287 * t298 + (t407 * t494 + t410 * t495) * t391 + (-t266 * t408 - t391 * t435) * qJD(6) + (-t306 * t462 + (qJD(2) * t436 + t265) * qJD(3)) * t411) * MDP(25) + (-t439 * t408 - t337 * t432 + t346 * t279 - t281 * t340 + t287 * t297 + (-t407 * t495 + t410 * t494) * t391 + (-t265 * t408 - t391 * t436) * qJD(6) + (t432 * t462 + (-qJD(2) * t435 - t266) * qJD(3)) * t411) * MDP(26); -t496 * t515 + t414 * t516 + (-t377 * t482 + t426) * MDP(10) + (-t377 * t481 - t427) * MDP(11) + (-t366 * t481 + t325 - t426) * MDP(13) + (0.2e1 * t398 + (t333 * t411 + t366 * t408) * qJD(2) + t427) * MDP(14) + (-pkin(3) * t303 - qJ(4) * t299 - t326 * t332 - t327 * t517 - t333 * t366) * MDP(15) + (t507 + t471 * t358 + (-t274 * t411 - t283 * t408 - t404 * t512) * qJD(2)) * MDP(16) + (t506 + t471 * t360 + (t275 * t411 + t284 * t408 + t402 * t512) * qJD(2)) * MDP(17) + (t283 * t360 + t284 * t358 + (qJD(5) * t360 - t275 * t482 - t270) * t404 + (qJD(5) * t358 + t274 * t482 - t271) * t402) * MDP(18) + (qJ(4) * t288 - t274 * t283 - t275 * t284 + t438 * t406 + t471 * t309 + (-t274 * t404 - t275 * t402) * qJD(5)) * MDP(19) + (-t279 * t431 - t432 * t490) * MDP(20) + (-t279 * t362 + t280 * t431 - t306 * t490 + t432 * t489) * MDP(21) + (-t336 * t391 + t432 * t481 + t491) * MDP(22) + (-t489 * t391 + (-qJD(3) * t362 + t306) * t481) * MDP(23) - t391 * t451 + (t392 * t280 + t281 * t362 + (t407 * t423 - t410 * t424) * t391 + t472 * t306 + t489 * t287 + ((-t369 * t407 + t370 * t410) * qJD(3) - t265) * t481) * MDP(25) + (t392 * t279 - t281 * t431 + (t407 * t424 + t410 * t423) * t391 - t472 * t432 + t490 * t287 + (-(t369 * t410 + t370 * t407) * qJD(3) + t266) * t481) * MDP(26); -t413 * MDP(14) + (t325 + t303) * MDP(15) + t438 * MDP(19) + t491 * MDP(25) + (-MDP(16) * t402 - MDP(17) * t404 - MDP(14)) * t414 * t400 + (-t336 * MDP(25) - MDP(26) * t489) * t391 + ((-t358 + t458) * MDP(16) + (-t360 - t456) * MDP(17) + (-t362 * t481 + t432) * MDP(26) + t419) * qJD(3) + (MDP(13) * t496 + ((-t358 * t404 + t360 * t402) * MDP(18) + t437 * MDP(19)) * qJD(2)) * t408; (-t358 ^ 2 - t360 ^ 2) * MDP(18) + (t274 * t360 + t275 * t358 + t288) * MDP(19) + (t280 - t519) * MDP(25) + (t279 - t505) * MDP(26) + ((t360 - t479) * MDP(16) + (-t358 + t480) * MDP(17)) * t482; -t432 * t306 * MDP(20) + (-t306 ^ 2 + t432 ^ 2) * MDP(21) + (t466 + t505) * MDP(22) + (-t421 - t519) * MDP(23) + qJD(3) * t451 + (t266 * t391 + t287 * t432 + t447) * MDP(25) + (t265 * t391 + t287 * t306 - t439) * MDP(26) + (-MDP(22) * t504 + MDP(23) * t432 - MDP(25) * t266 - MDP(26) * t265) * qJD(6);];
tauc  = t1;
