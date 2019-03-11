% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:27
% EndTime: 2019-03-08 23:59:37
% DurationCPUTime: 4.81s
% Computational Cost: add. (4463->373), mult. (10921->517), div. (0->0), fcn. (8166->10), ass. (0->183)
t421 = sin(qJ(3));
t422 = sin(qJ(2));
t417 = sin(pkin(6));
t490 = qJD(1) * t417;
t474 = t422 * t490;
t517 = qJD(3) * pkin(3);
t525 = -t421 * t517 + t474;
t420 = sin(qJ(4));
t424 = cos(qJ(3));
t521 = cos(qJ(4));
t441 = -t420 * t421 + t424 * t521;
t481 = qJD(3) + qJD(4);
t360 = t481 * t441;
t504 = t420 * t424;
t391 = t421 * t521 + t504;
t361 = t481 * t391;
t536 = pkin(4) * t361 - pkin(10) * t360 - t525;
t522 = -pkin(9) - pkin(8);
t476 = qJD(3) * t522;
t393 = t421 * t476;
t394 = t424 * t476;
t399 = t522 * t421;
t400 = t522 * t424;
t524 = t521 * t399 + t420 * t400;
t321 = qJD(4) * t524 + t521 * t393 + t420 * t394;
t425 = cos(qJ(2));
t473 = t425 * t490;
t369 = t441 * t473;
t535 = t321 - t369;
t463 = -qJD(2) * t522 + t474;
t418 = cos(pkin(6));
t489 = qJD(1) * t418;
t366 = t421 * t489 + t424 * t463;
t355 = t420 * t366;
t365 = -t463 * t421 + t424 * t489;
t319 = t365 * t521 - t355;
t467 = t521 * qJD(4);
t533 = -pkin(3) * t467 + t319;
t419 = sin(qJ(5));
t423 = cos(qJ(5));
t534 = t369 * t419 + t423 * t536;
t412 = -pkin(3) * t424 - pkin(2);
t352 = -pkin(4) * t441 - pkin(10) * t391 + t412;
t483 = qJD(5) * t423;
t532 = t352 * t483 + t419 * t536 + t535 * t423;
t531 = MDP(5) * t424;
t530 = MDP(6) * (t421 ^ 2 - t424 ^ 2);
t375 = t420 * t399 - t400 * t521;
t364 = t423 * t375;
t445 = -qJ(6) * t360 - qJD(6) * t391;
t529 = pkin(5) * t361 - t321 * t419 + t445 * t423 + (-t364 + (qJ(6) * t391 - t352) * t419) * qJD(5) + t534;
t469 = t391 * t483;
t528 = -qJ(6) * t469 + (-qJD(5) * t375 + t445) * t419 + t532;
t497 = qJD(4) * t375 - t391 * t473 + t420 * t393 - t394 * t521;
t468 = qJD(2) * t521;
t487 = qJD(2) * t421;
t385 = t420 * t487 - t424 * t468;
t512 = t385 * t419;
t527 = -qJ(6) * t512 + t423 * qJD(6);
t387 = -qJD(2) * t504 - t421 * t468;
t351 = -pkin(4) * t387 + pkin(10) * t385;
t339 = pkin(3) * t487 + t351;
t526 = t419 * t339 + t423 * t533;
t429 = t360 * qJD(2);
t435 = t423 * t387 - t419 * t481;
t314 = -qJD(5) * t435 + t419 * t429;
t523 = t435 ^ 2;
t520 = t423 * pkin(5);
t519 = -qJ(6) - pkin(10);
t518 = qJD(2) * pkin(2);
t459 = t423 * t481;
t484 = qJD(5) * t419;
t313 = -qJD(5) * t459 - t387 * t484 - t423 * t429;
t516 = t313 * t419;
t349 = t361 * qJD(2);
t515 = t349 * t423;
t370 = -t387 * t419 - t459;
t381 = qJD(5) + t385;
t514 = t370 * t381;
t513 = t435 * t381;
t511 = t385 * t423;
t510 = t391 * t419;
t509 = t391 * t423;
t508 = t417 * t422;
t507 = t417 * t425;
t427 = qJD(2) ^ 2;
t506 = t417 * t427;
t505 = t419 * t349;
t426 = qJD(3) ^ 2;
t503 = t421 * t426;
t414 = t423 * qJ(6);
t502 = t424 * t426;
t410 = pkin(3) * t420 + pkin(10);
t501 = -qJ(6) - t410;
t356 = t521 * t366;
t357 = t365 + t517;
t317 = t420 * t357 + t356;
t310 = pkin(10) * t481 + t317;
t377 = qJD(2) * t412 - t473;
t329 = pkin(4) * t385 + pkin(10) * t387 + t377;
t294 = -t310 * t419 + t423 * t329;
t286 = qJ(6) * t435 + t294;
t279 = pkin(5) * t381 + t286;
t500 = t279 - t286;
t316 = t357 * t521 - t355;
t499 = t423 * t316 + t419 * t351;
t461 = qJD(5) * t501;
t496 = t419 * t461 - t526 + t527;
t336 = t423 * t339;
t452 = -t387 * pkin(5) + t385 * t414;
t495 = t423 * t461 - t336 - t452 + (-qJD(6) + t533) * t419;
t494 = t419 * t352 + t364;
t464 = qJD(5) * t519;
t493 = t419 * t464 - t499 + t527;
t462 = -t316 * t419 + t423 * t351;
t492 = -t419 * qJD(6) + t423 * t464 - t452 - t462;
t482 = qJD(2) * qJD(3);
t465 = t421 * t482;
t488 = qJD(2) * t417;
t466 = qJD(1) * t488;
t382 = pkin(3) * t465 + t422 * t466;
t486 = qJD(2) * t422;
t485 = qJD(4) * t420;
t478 = t422 * t506;
t471 = t417 * t486;
t470 = t425 * t488;
t309 = -pkin(4) * t481 - t316;
t306 = t309 * t483;
t460 = t381 * t423;
t455 = t425 * t466;
t333 = qJD(3) * t365 + t424 * t455;
t334 = -qJD(3) * t366 - t421 * t455;
t285 = t420 * t333 - t521 * t334 + t357 * t485 + t366 * t467;
t457 = t421 * t470;
t456 = t424 * t470;
t411 = -pkin(3) * t521 - pkin(4);
t318 = t420 * t365 + t356;
t454 = pkin(3) * t485 - t318;
t453 = (t484 + t512) * pkin(5);
t295 = t310 * t423 + t329 * t419;
t451 = t285 * t419 - t295 * t387 + t306;
t287 = -qJ(6) * t370 + t295;
t447 = -t279 * t423 - t287 * t419;
t446 = t309 * t385 - t349 * t410;
t383 = t418 * t424 - t421 * t508;
t384 = t418 * t421 + t424 * t508;
t341 = t420 * t383 + t384 * t521;
t327 = -t341 * t419 - t423 * t507;
t444 = -t341 * t423 + t419 * t507;
t443 = -t285 * t423 + t294 * t387 + t309 * t484;
t442 = t383 * t521 - t420 * t384;
t440 = t360 * t419 + t469;
t439 = t360 * t423 - t391 * t484;
t276 = pkin(5) * t314 + t285;
t438 = qJD(2) * t518;
t437 = t377 * t387 - t285;
t284 = t333 * t521 + t420 * t334 + t357 * t467 - t366 * t485;
t303 = t349 * pkin(4) - pkin(10) * t429 + t382;
t436 = t423 * t284 + t419 * t303 - t310 * t484 + t329 * t483;
t434 = -0.2e1 * qJD(3) * t518;
t302 = t423 * t303;
t433 = -qJD(5) * t295 - t284 * t419 + t302;
t270 = pkin(5) * t349 + qJ(6) * t313 + qJD(6) * t435 + t433;
t272 = -qJ(6) * t314 - qJD(6) * t370 + t436;
t432 = qJD(5) * t447 - t270 * t419 + t272 * t423 - t279 * t511 - t287 * t512;
t431 = ((-t313 - t514) * t423 + (-t314 + t513) * t419) * MDP(20) + (-t435 * t460 - t516) * MDP(19) + (-t381 ^ 2 * t419 - t370 * t387 + t515) * MDP(22) + (t381 * t460 - t387 * t435 + t505) * MDP(21) + t429 * MDP(14) + (-t385 ^ 2 + t387 ^ 2) * MDP(13) + (-MDP(12) * t385 + t381 * MDP(23)) * t387 + (t385 * MDP(14) + (-qJD(2) * t391 - t387) * MDP(15)) * t481;
t430 = t377 * t385 - t284;
t398 = pkin(10) * t423 + t414;
t397 = t519 * t419;
t389 = t410 * t423 + t414;
t388 = t501 * t419;
t367 = t370 ^ 2;
t363 = -qJD(3) * t384 - t457;
t362 = qJD(3) * t383 + t456;
t347 = t423 * t352;
t304 = -qJ(6) * t510 + t494;
t300 = -pkin(5) * t441 - t375 * t419 - t391 * t414 + t347;
t299 = t370 * pkin(5) + qJD(6) + t309;
t298 = qJD(4) * t341 + t420 * t362 - t363 * t521;
t297 = qJD(4) * t442 + t362 * t521 + t420 * t363;
t281 = qJD(5) * t444 - t297 * t419 + t423 * t471;
t280 = qJD(5) * t327 + t297 * t423 + t419 * t471;
t1 = [-MDP(3) * t478 - t425 * MDP(4) * t506 + (-t424 * t478 + (t363 - t457) * qJD(3)) * MDP(10) + (t421 * t478 + (-t362 - t456) * qJD(3)) * MDP(11) + (-t298 * t481 + (-t349 * t425 + t385 * t486) * t417) * MDP(17) + (-t297 * t481 + (-t360 * t425 - t422 * t387) * t488) * MDP(18) + (t281 * t381 + t298 * t370 - t314 * t442 + t327 * t349) * MDP(24) + (-t280 * t381 - t298 * t435 + t313 * t442 + t349 * t444) * MDP(25) + (-t280 * t370 + t281 * t435 + t313 * t327 + t314 * t444) * MDP(26) + (t270 * t327 - t272 * t444 - t276 * t442 + t279 * t281 + t280 * t287 + t298 * t299) * MDP(27); 0.2e1 * t465 * t531 - 0.2e1 * t482 * t530 + MDP(7) * t502 - MDP(8) * t503 + (-pkin(8) * t502 + t421 * t434) * MDP(10) + (pkin(8) * t503 + t424 * t434) * MDP(11) + (-t387 * t360 + t391 * t429) * MDP(12) + (-t391 * t349 - t360 * t385 + t387 * t361 + t429 * t441) * MDP(13) + (t412 * t349 + t377 * t361 - t382 * t441 - t385 * t525) * MDP(17) + (t377 * t360 + t382 * t391 + t525 * t387 + t412 * t429) * MDP(18) + (-t313 * t509 - t435 * t439) * MDP(19) + ((-t370 * t423 + t419 * t435) * t360 + (t516 - t314 * t423 + (t370 * t419 + t423 * t435) * qJD(5)) * t391) * MDP(20) + (t313 * t441 + t349 * t509 - t361 * t435 + t381 * t439) * MDP(21) + (t314 * t441 - t361 * t370 - t381 * t440 - t391 * t505) * MDP(22) + (-t349 * t441 + t361 * t381) * MDP(23) + (t347 * t349 - (-t310 * t483 + t302) * t441 + t294 * t361 - t524 * t314 + t391 * t306 + (-t375 * t483 + t534) * t381 + t497 * t370 + ((-qJD(5) * t352 - t321) * t381 - t375 * t349 - (-qJD(5) * t329 - t284) * t441 + t285 * t391 + t309 * t360) * t419) * MDP(24) + (-t494 * t349 + t436 * t441 - t295 * t361 + t524 * t313 + t285 * t509 + (t375 * t484 - t532) * t381 - t497 * t435 + t439 * t309) * MDP(25) + (t300 * t313 - t304 * t314 + t529 * t435 - t528 * t370 + t447 * t360 + (-t270 * t423 - t272 * t419 + (t279 * t419 - t287 * t423) * qJD(5)) * t391) * MDP(26) + (t272 * t304 + t270 * t300 + t276 * (pkin(5) * t510 - t524) + (pkin(5) * t440 + t497) * t299 + t528 * t287 + t529 * t279) * MDP(27) + (t360 * MDP(14) - t361 * MDP(15) - t497 * MDP(17) - MDP(18) * t535) * t481; (t319 * t481 + (t387 * t487 - t467 * t481) * pkin(3) + t430) * MDP(18) + t427 * t530 + t424 * t438 * MDP(11) + t431 + (t313 * t388 - t314 * t389 - t370 * t496 + t435 * t495 + t432) * MDP(26) + (t272 * t389 + t270 * t388 + t276 * (t411 - t520) + (-t356 + (pkin(3) * qJD(4) - t365) * t420 + t453) * t299 + t496 * t287 + t495 * t279) * MDP(27) + (t318 * t481 + (-t385 * t487 - t481 * t485) * pkin(3) + t437) * MDP(17) + (t411 * t314 + t446 * t419 + t454 * t370 + (-t410 * t483 + t419 * t533 - t336) * t381 + t443) * MDP(24) + (-t411 * t313 + t446 * t423 - t454 * t435 + (t410 * t484 + t526) * t381 + t451) * MDP(25) + (MDP(10) * t438 - t427 * t531) * t421; (t316 * t481 + t430) * MDP(18) + t431 + (t313 * t397 - t314 * t398 - t370 * t493 + t435 * t492 + t432) * MDP(26) + (t272 * t398 + t270 * t397 + t276 * (-pkin(4) - t520) + (t453 - t317) * t299 + t493 * t287 + t492 * t279) * MDP(27) + (-pkin(4) * t314 - t462 * t381 - t317 * t370 + t309 * t512 + (-t381 * t483 - t505) * pkin(10) + t443) * MDP(24) + (pkin(4) * t313 + t499 * t381 + t317 * t435 + t309 * t511 + (t381 * t484 - t515) * pkin(10) + t451) * MDP(25) + (t317 * t481 + t437) * MDP(17); -t435 * t370 * MDP(19) + (-t367 + t523) * MDP(20) + (-t313 + t514) * MDP(21) + (-t314 - t513) * MDP(22) + t349 * MDP(23) + (t295 * t381 + t309 * t435 + t433) * MDP(24) + (t294 * t381 + t309 * t370 - t436) * MDP(25) + (pkin(5) * t313 - t370 * t500) * MDP(26) + (t500 * t287 + (t299 * t435 + t270) * pkin(5)) * MDP(27); (-t367 - t523) * MDP(26) + (-t279 * t435 + t287 * t370 + t276) * MDP(27);];
tauc  = t1;
