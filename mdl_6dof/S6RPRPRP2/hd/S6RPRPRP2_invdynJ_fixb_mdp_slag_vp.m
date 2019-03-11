% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:22
% EndTime: 2019-03-09 03:06:30
% DurationCPUTime: 6.16s
% Computational Cost: add. (5178->470), mult. (10935->580), div. (0->0), fcn. (7590->14), ass. (0->200)
t556 = -MDP(20) + MDP(23);
t425 = sin(pkin(10));
t430 = sin(qJ(3));
t433 = cos(qJ(3));
t531 = cos(pkin(10));
t387 = t425 * t433 + t430 * t531;
t377 = t387 * qJD(3);
t484 = t531 * t433;
t494 = qJDD(1) * t430;
t468 = qJDD(1) * t484 - t425 * t494;
t350 = qJD(1) * t377 + qJDD(5) - t468;
t432 = cos(qJ(5));
t344 = t432 * t350;
t400 = qJD(1) * t484;
t501 = qJD(1) * t430;
t376 = -t425 * t501 + t400;
t373 = qJD(5) - t376;
t429 = sin(qJ(5));
t499 = qJD(5) * t429;
t459 = -t425 * t430 + t484;
t380 = t459 * qJD(3);
t521 = t380 * t432;
t460 = t387 * t499 - t521;
t458 = t344 * t387 - t373 * t460;
t426 = sin(pkin(9));
t405 = pkin(1) * t426 + pkin(7);
t514 = qJ(4) + t405;
t427 = cos(pkin(9));
t407 = -t427 * pkin(1) - pkin(2);
t553 = qJDD(1) * t407;
t421 = qJ(3) + pkin(10);
t410 = sin(t421);
t422 = qJ(1) + pkin(9);
t411 = sin(t422);
t413 = cos(t422);
t474 = g(1) * t413 + g(2) * t411;
t552 = t410 * t474;
t378 = t387 * qJD(1);
t361 = qJD(3) * t429 + t378 * t432;
t496 = qJD(1) * qJD(3);
t487 = t430 * t496;
t545 = qJD(3) * t400 + qJDD(1) * t387 - t425 * t487;
t438 = -qJDD(3) * t432 + t429 * t545;
t310 = qJD(5) * t361 + t438;
t359 = -qJD(3) * t432 + t378 * t429;
t551 = t310 * t459 - t359 * t377;
t495 = qJD(3) * qJD(5);
t309 = -t429 * qJDD(3) + t378 * t499 + (-t545 - t495) * t432;
t509 = t309 * t459 + t361 * t377;
t419 = t433 * pkin(3);
t550 = t407 - t419;
t342 = -pkin(4) * t459 - pkin(8) * t387 + t550;
t382 = t514 * t430;
t383 = t514 * t433;
t347 = -t382 * t425 + t383 * t531;
t508 = t342 * t429 + t347 * t432;
t486 = -g(1) * t411 + g(2) * t413;
t492 = pkin(3) * t487 + qJDD(4);
t549 = qJDD(1) * t550 + t492;
t412 = cos(t421);
t540 = pkin(8) * t410;
t548 = pkin(4) * t412 + t540;
t547 = MDP(19) + MDP(21);
t479 = t514 * qJD(1);
t364 = qJD(2) * t430 + t433 * t479;
t354 = t425 * t364;
t363 = qJD(2) * t433 - t430 * t479;
t532 = qJD(3) * pkin(3);
t357 = t363 + t532;
t319 = t357 * t531 - t354;
t313 = -qJD(3) * pkin(4) - t319;
t296 = pkin(5) * t359 - qJ(6) * t361 + t313;
t543 = pkin(3) * t425;
t404 = pkin(8) + t543;
t519 = t404 * t350;
t546 = t296 * t373 - t519;
t544 = t361 ^ 2;
t541 = pkin(5) * t350;
t536 = g(3) * t410;
t535 = g(3) * t412;
t534 = g(3) * t433;
t530 = qJ(6) * t350;
t485 = t531 * t364;
t320 = t357 * t425 + t485;
t314 = qJD(3) * pkin(8) + t320;
t375 = qJD(1) * t550 + qJD(4);
t331 = -pkin(4) * t376 - pkin(8) * t378 + t375;
t295 = t314 * t432 + t331 * t429;
t529 = t295 * t373;
t528 = t309 * t429;
t527 = t359 * t376;
t526 = t359 * t429;
t525 = t361 * t359;
t480 = t361 * t373;
t524 = t361 * t432;
t523 = t373 * t429;
t522 = t380 * t429;
t520 = t387 * t432;
t518 = t411 * t429;
t517 = t411 * t432;
t516 = t413 * t429;
t515 = t413 * t432;
t343 = t429 * t350;
t513 = qJDD(2) - g(3);
t512 = -t310 * t520 - t359 * t521;
t498 = qJD(5) * t432;
t511 = -t310 * t429 - t359 * t498;
t416 = t433 * qJDD(2);
t391 = t405 * qJDD(1);
t449 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t391;
t464 = t479 * qJD(3);
t317 = qJDD(3) * pkin(3) - t430 * t449 - t433 * t464 + t416;
t321 = (qJDD(2) - t464) * t430 + t449 * t433;
t291 = t317 * t425 + t321 * t531;
t323 = t363 * t531 - t354;
t337 = pkin(3) * t501 + pkin(4) * t378 - pkin(8) * t376;
t510 = t323 * t432 + t337 * t429;
t507 = t373 * t498 + t343;
t506 = t376 * t523 + t344;
t322 = t363 * t425 + t485;
t469 = pkin(5) * t429 - qJ(6) * t432;
t505 = -qJD(6) * t429 + t373 * t469 - t322;
t504 = (g(1) * t515 + g(2) * t517) * t410;
t409 = t419 + pkin(2);
t434 = cos(qJ(1));
t503 = pkin(1) * t434 + t409 * t413;
t423 = t430 ^ 2;
t502 = -t433 ^ 2 + t423;
t394 = qJD(1) * t407;
t500 = qJD(5) * t404;
t294 = -t314 * t429 + t331 * t432;
t497 = qJD(6) - t294;
t493 = qJDD(1) * t433;
t491 = t430 * t532;
t490 = t361 * t522;
t489 = t531 * pkin(3);
t289 = qJDD(3) * pkin(8) + t291;
t351 = -qJD(3) * t378 + t468;
t303 = -pkin(3) * t493 - t351 * pkin(4) - pkin(8) * t545 + t492 + t553;
t457 = t289 * t432 + t303 * t429 - t314 * t499 + t331 * t498;
t278 = qJD(6) * t373 + t457 + t530;
t284 = -pkin(5) * t373 + t497;
t483 = -t284 * t376 + t278;
t478 = -t289 * t429 + t303 * t432 - t314 * t498 - t331 * t499;
t279 = qJDD(6) - t478 - t541;
t285 = qJ(6) * t373 + t295;
t482 = t285 * t376 + t279;
t481 = qJD(3) * t514;
t365 = qJD(4) * t433 - t430 * t481;
t366 = -qJD(4) * t430 - t433 * t481;
t329 = t365 * t425 - t366 * t531;
t346 = t382 * t531 + t383 * t425;
t368 = t412 * t518 + t515;
t370 = t412 * t516 - t517;
t476 = -g(1) * t368 + g(2) * t370;
t369 = t412 * t517 - t516;
t371 = t412 * t515 + t518;
t475 = g(1) * t369 - g(2) * t371;
t431 = sin(qJ(1));
t472 = g(1) * t431 - g(2) * t434;
t428 = -qJ(4) - pkin(7);
t471 = -pkin(1) * t431 - t413 * t428;
t470 = pkin(5) * t432 + qJ(6) * t429;
t290 = t317 * t531 - t321 * t425;
t467 = t284 * t432 - t285 * t429;
t466 = t284 * t429 + t285 * t432;
t463 = pkin(4) + t470;
t462 = t373 * t500 + t535;
t461 = t387 * t498 + t522;
t330 = t365 * t531 + t366 * t425;
t338 = pkin(4) * t377 - pkin(8) * t380 + t491;
t456 = t330 * t432 + t338 * t429 + t342 * t498 - t347 * t499;
t455 = t313 * t373 - t519;
t288 = -qJDD(3) * pkin(4) - t290;
t280 = pkin(5) * t310 + qJ(6) * t309 - qJD(6) * t361 + t288;
t454 = -t280 - t462;
t451 = -qJD(1) * t394 - t391 + t474;
t450 = 0.2e1 * qJD(3) * t394 - qJDD(3) * t405;
t448 = g(1) * t370 + g(2) * t368 + t429 * t536 + t478;
t435 = qJD(3) ^ 2;
t447 = -t405 * t435 - t486 - 0.2e1 * t553;
t446 = -t343 * t387 - t373 * t461;
t444 = qJD(5) * t467 + t278 * t432 + t279 * t429;
t443 = t296 * t361 + qJDD(6) - t448;
t441 = -g(1) * t371 - g(2) * t369 - t432 * t536 + t457;
t406 = -t489 - pkin(4);
t390 = qJDD(3) * t433 - t430 * t435;
t389 = qJDD(3) * t430 + t433 * t435;
t381 = -t489 - t463;
t327 = pkin(5) * t361 + qJ(6) * t359;
t308 = t387 * t469 + t346;
t301 = pkin(5) * t459 - t342 * t432 + t347 * t429;
t298 = -qJ(6) * t459 + t508;
t297 = t359 * t373 - t309;
t293 = -pkin(5) * t378 + t323 * t429 - t337 * t432;
t292 = qJ(6) * t378 + t510;
t283 = t469 * t380 + (qJD(5) * t470 - qJD(6) * t432) * t387 + t329;
t282 = -pkin(5) * t377 + qJD(5) * t508 + t330 * t429 - t338 * t432;
t281 = qJ(6) * t377 - qJD(6) * t459 + t456;
t1 = [qJDD(1) * MDP(1) + t472 * MDP(2) + (g(1) * t434 + g(2) * t431) * MDP(3) + (t472 + (t426 ^ 2 + t427 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t423 + 0.2e1 * t433 * t487) * MDP(5) + 0.2e1 * (t430 * t493 - t496 * t502) * MDP(6) + t389 * MDP(7) + t390 * MDP(8) + (t430 * t450 + t433 * t447) * MDP(10) + (-t430 * t447 + t433 * t450) * MDP(11) + (-t290 * t387 + t291 * t459 - t319 * t380 - t320 * t377 + t329 * t378 + t330 * t376 + t346 * t545 + t347 * t351 - t474) * MDP(12) + (t291 * t347 + t320 * t330 - t290 * t346 - t319 * t329 + t549 * t550 + t375 * t491 - g(1) * (-t409 * t411 + t471) - g(2) * (-t411 * t428 + t503)) * MDP(13) + (-t309 * t520 - t361 * t460) * MDP(14) + (-t490 + (t528 + (-t524 + t526) * qJD(5)) * t387 + t512) * MDP(15) + (t458 + t509) * MDP(16) + (t446 + t551) * MDP(17) + (-t350 * t459 + t373 * t377) * MDP(18) + (-t478 * t459 + t294 * t377 + t329 * t359 + t346 * t310 + ((-qJD(5) * t347 + t338) * t373 + t342 * t350 + t313 * qJD(5) * t387) * t432 + ((-qJD(5) * t342 - t330) * t373 - t347 * t350 + t288 * t387 + t313 * t380) * t429 + t475) * MDP(19) + (t288 * t520 - t295 * t377 - t346 * t309 - t313 * t460 + t329 * t361 - t350 * t508 - t373 * t456 + t457 * t459 + t476) * MDP(20) + (t280 * t387 * t429 + t279 * t459 - t282 * t373 + t283 * t359 - t284 * t377 + t296 * t461 - t301 * t350 + t308 * t310 + t475) * MDP(21) + (-t281 * t359 + t282 * t361 - t298 * t310 - t301 * t309 - t486 * t410 + t467 * t380 + (-qJD(5) * t466 - t278 * t429 + t279 * t432) * t387) * MDP(22) + (-t278 * t459 - t280 * t520 + t281 * t373 - t283 * t361 + t285 * t377 + t296 * t460 + t298 * t350 + t308 * t309 - t476) * MDP(23) + (t278 * t298 + t285 * t281 + t280 * t308 + t296 * t283 + t279 * t301 + t284 * t282 - g(1) * (-pkin(5) * t369 - qJ(6) * t368 + t471) - g(2) * (pkin(5) * t371 + qJ(6) * t370 + t413 * t548 + t503) + (-g(1) * (-t409 - t548) + g(2) * t428) * t411) * MDP(24); t513 * MDP(4) + t390 * MDP(10) - t389 * MDP(11) + (t387 * t351 + t380 * t376 + t377 * t378 - t459 * t545) * MDP(12) + (t290 * t459 + t291 * t387 - t319 * t377 + t320 * t380 - g(3)) * MDP(13) + (t490 + (-t528 + (t524 + t526) * qJD(5)) * t387 + t512) * MDP(22) + (-t280 * t459 + t296 * t377 + t380 * t466 + t387 * t444 - g(3)) * MDP(24) + t547 * (t446 - t551) + t556 * (t458 - t509); MDP(7) * t494 + MDP(8) * t493 + qJDD(3) * MDP(9) + (t430 * t451 + t416 - t534) * MDP(10) + (-t430 * t513 + t433 * t451) * MDP(11) + (t351 * t543 - t545 * t489 - (-t320 + t322) * t378 + (t319 - t323) * t376) * MDP(12) + (t319 * t322 - t320 * t323 + (t531 * t290 - t534 + t291 * t425 + (-qJD(1) * t375 + t474) * t430) * pkin(3)) * MDP(13) + (t432 * t480 - t528) * MDP(14) + ((-t309 + t527) * t432 - t361 * t523 + t511) * MDP(15) + (-t373 * t376 * t432 - t361 * t378 + t507) * MDP(16) + (t359 * t378 - t373 * t499 + t506) * MDP(17) - t373 * t378 * MDP(18) + (-t294 * t378 + t406 * t310 - t322 * t359 + (-t535 - t288 + (-t337 - t500) * t373) * t432 + (t323 * t373 + t455) * t429 + t504) * MDP(19) + (-t406 * t309 + t510 * t373 + t295 * t378 - t322 * t361 + t455 * t432 + (t288 + t462 - t552) * t429) * MDP(20) + (t284 * t378 + t293 * t373 + t310 * t381 + t505 * t359 + t429 * t546 + t454 * t432 + t504) * MDP(21) + (-t536 + t292 * t359 - t293 * t361 - t474 * t412 + (-t310 * t404 + (t361 * t404 + t284) * qJD(5) + t483) * t432 + (-t309 * t404 + (t359 * t404 - t285) * qJD(5) + t482) * t429) * MDP(22) + (-t285 * t378 - t292 * t373 + t309 * t381 - t505 * t361 - t546 * t432 + (t454 + t552) * t429) * MDP(23) + (t280 * t381 - t285 * t292 - t284 * t293 - g(3) * (t419 + t540) - t463 * t535 + t505 * t296 + t444 * t404 + t474 * (pkin(3) * t430 - pkin(8) * t412 + t410 * t463)) * MDP(24) + (-MDP(5) * t430 * t433 + MDP(6) * t502) * qJD(1) ^ 2; -t376 ^ 2 * MDP(12) + t506 * MDP(19) + t511 * MDP(22) + t507 * MDP(23) + t486 * MDP(24) + (-MDP(12) * t378 - t296 * MDP(24) - t547 * t359 + t361 * t556) * t378 + (t350 * MDP(21) + (t309 + t527) * MDP(22) + (qJD(5) * t285 - t482) * MDP(24) + (-MDP(20) * t373 - MDP(23) * t376) * t373) * t432 + (-t350 * MDP(20) + (qJD(5) * t284 + t483) * MDP(24) + MDP(22) * t480 + (-MDP(19) * qJD(5) - MDP(21) * t373) * t373) * t429 + (t319 * t378 - t320 * t376 + t486 + t549) * MDP(13); MDP(14) * t525 + (-t359 ^ 2 + t544) * MDP(15) + t297 * MDP(16) + (-t378 * t498 - t429 * t495 - t438 + t480) * MDP(17) + t350 * MDP(18) + (-t313 * t361 + t448 + t529) * MDP(19) + (t294 * t373 + t313 * t359 - t441) * MDP(20) + (-t327 * t359 - t443 + t529 + 0.2e1 * t541) * MDP(21) + (pkin(5) * t309 - qJ(6) * t310 + (t285 - t295) * t361 + (t284 - t497) * t359) * MDP(22) + (0.2e1 * t530 - t296 * t359 + t327 * t361 + (0.2e1 * qJD(6) - t294) * t373 + t441) * MDP(23) + (t278 * qJ(6) - t279 * pkin(5) - t296 * t327 - t284 * t295 - g(1) * (-pkin(5) * t370 + qJ(6) * t371) - g(2) * (-pkin(5) * t368 + qJ(6) * t369) + t469 * t536 + t497 * t285) * MDP(24); (-qJDD(5) + t351 + t525) * MDP(21) + t297 * MDP(22) + (-t373 ^ 2 - t544) * MDP(23) + (-t285 * t373 + t443 - t541) * MDP(24);];
tau  = t1;
