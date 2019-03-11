% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:40
% EndTime: 2019-03-09 06:32:53
% DurationCPUTime: 5.89s
% Computational Cost: add. (5311->480), mult. (11222->625), div. (0->0), fcn. (7152->6), ass. (0->197)
t434 = sin(qJ(4));
t437 = cos(qJ(4));
t504 = qJD(3) * t437;
t438 = cos(qJ(3));
t509 = qJD(1) * t438;
t391 = -t434 * t509 + t504;
t481 = t437 * t509;
t506 = qJD(3) * t434;
t392 = t481 + t506;
t433 = sin(qJ(5));
t436 = cos(qJ(5));
t344 = -t436 * t391 + t392 * t433;
t457 = t391 * t433 + t436 * t392;
t571 = t344 * t457;
t435 = sin(qJ(3));
t463 = pkin(3) * t438 + pkin(8) * t435;
t396 = t463 * qJD(1);
t382 = t437 * t396;
t551 = -pkin(9) - pkin(8);
t483 = qJD(4) * t551;
t439 = -pkin(1) - pkin(7);
t416 = qJD(1) * t439 + qJD(2);
t536 = t434 * t438;
t487 = t416 * t536;
t533 = t435 * t437;
t491 = pkin(9) * t533;
t570 = -t487 + t382 + (pkin(4) * t438 + t491) * qJD(1) - t437 * t483;
t510 = qJD(1) * t435;
t482 = t434 * t510;
t529 = t437 * t438;
t516 = t434 * t396 + t416 * t529;
t569 = pkin(9) * t482 - t434 * t483 + t516;
t499 = qJD(4) * t438;
t477 = t434 * t499;
t479 = t435 * t504;
t447 = -t477 - t479;
t494 = qJD(3) * qJD(4);
t444 = qJD(1) * t447 + t437 * t494;
t505 = qJD(3) * t435;
t480 = t434 * t505;
t474 = qJD(1) * t499;
t515 = t434 * t494 + t437 * t474;
t448 = qJD(1) * t480 - t515;
t497 = qJD(5) * t436;
t498 = qJD(5) * t433;
t305 = -t391 * t497 + t392 * t498 - t433 * t448 - t436 * t444;
t424 = qJD(4) + t510;
t414 = qJD(5) + t424;
t289 = t344 * t414 - t305;
t306 = qJD(5) * t457 + t433 * t444 - t436 * t448;
t495 = qJD(1) * qJD(3);
t475 = t438 * t495;
t552 = t457 ^ 2;
t568 = MDP(25) * t475 + (t414 * t457 - t306) * MDP(24) + MDP(21) * t571 + (-t344 ^ 2 + t552) * MDP(22) + t289 * MDP(23);
t550 = qJD(3) * pkin(3);
t384 = -t416 * t438 - t550;
t353 = -pkin(4) * t391 + t384;
t302 = pkin(5) * t344 - qJ(6) * t457 + t353;
t567 = t302 * t344;
t566 = t344 * t353;
t493 = qJD(4) + qJD(5);
t500 = qJD(4) * t437;
t531 = t436 * t437;
t538 = t433 * t434;
t518 = t433 * t482 - t436 * t500 - t437 * t497 + t493 * t538 - t510 * t531;
t395 = t433 * t437 + t434 * t436;
t351 = t493 * t395;
t378 = t395 * qJD(1);
t517 = t435 * t378 + t351;
t476 = t437 * t499;
t565 = t476 - t480;
t314 = pkin(5) * t457 + qJ(6) * t344;
t492 = 0.2e1 * qJD(1);
t432 = t438 ^ 2;
t513 = t435 ^ 2 - t432;
t561 = MDP(8) * t513;
t548 = t302 * t457;
t560 = t457 * t353;
t409 = t551 * t434;
t410 = t551 * t437;
t456 = t409 * t436 + t410 * t433;
t559 = -qJD(5) * t456 + t433 * t570 + t436 * t569;
t358 = t409 * t433 - t410 * t436;
t558 = -qJD(5) * t358 + t433 * t569 - t436 * t570;
t462 = pkin(3) * t435 - pkin(8) * t438;
t400 = qJ(2) + t462;
t388 = t437 * t400;
t535 = t434 * t439;
t473 = pkin(4) - t535;
t342 = -pkin(9) * t529 + t435 * t473 + t388;
t532 = t435 * t439;
t413 = t437 * t532;
t514 = t434 * t400 + t413;
t352 = -pkin(9) * t536 + t514;
t557 = t433 * t342 + t436 * t352;
t399 = t435 * t416;
t502 = qJD(4) * t434;
t465 = -t399 + (t482 + t502) * pkin(4);
t389 = qJD(3) * t463 + qJD(2);
t374 = t437 * t389;
t312 = t374 + (-t413 + (pkin(9) * t438 - t400) * t434) * qJD(4) + (t438 * t473 + t491) * qJD(3);
t503 = qJD(3) * t438;
t478 = t437 * t503;
t484 = t434 * t389 + t400 * t500 + t439 * t478;
t501 = qJD(4) * t435;
t315 = -pkin(9) * t565 - t501 * t535 + t484;
t556 = -qJD(5) * t557 + t312 * t436 - t315 * t433;
t555 = t435 * t493;
t441 = qJD(1) ^ 2;
t549 = qJ(2) * t441;
t383 = qJD(3) * pkin(8) + t399;
t530 = t437 * t383;
t380 = t400 * qJD(1);
t537 = t434 * t380;
t338 = t530 + t537;
t330 = pkin(9) * t391 + t338;
t546 = t330 * t433;
t545 = t330 * t436;
t544 = t384 * t434;
t543 = t384 * t435;
t542 = t391 * t438;
t541 = t392 * t424;
t540 = t424 * t434;
t539 = t424 * t437;
t534 = t435 * t424;
t528 = t438 * t441;
t440 = qJD(3) ^ 2;
t527 = t439 * t440;
t337 = t437 * t380 - t383 * t434;
t329 = -pkin(9) * t392 + t337;
t297 = t329 * t436 - t546;
t526 = -pkin(4) * t497 - qJD(6) + t297;
t525 = qJ(6) * t509 + t559;
t524 = -pkin(5) * t509 + t558;
t523 = pkin(5) * t517 + qJ(6) * t518 - qJD(6) * t395 + t465;
t394 = -t531 + t538;
t371 = t394 * t438;
t522 = -qJD(3) * t371 - t395 * t555 - t378;
t369 = t395 * t438;
t521 = qJD(3) * t369 + (-qJD(1) - t555) * t394;
t508 = qJD(3) * t456;
t507 = qJD(3) * t358;
t323 = pkin(4) * t424 + t329;
t292 = t323 * t436 - t546;
t496 = qJD(6) - t292;
t488 = t434 * t534;
t486 = t424 * t536;
t367 = t389 * qJD(1);
t485 = t434 * t367 + t380 * t500 + t416 * t478;
t430 = -pkin(4) * t437 - pkin(3);
t472 = t424 * t439 + t383;
t390 = pkin(4) * t536 - t438 * t439;
t471 = -t391 + t504;
t470 = -t392 + t506;
t469 = -qJD(4) + t510;
t468 = pkin(5) * t475;
t361 = t437 * t367;
t298 = t361 + (-t530 + (pkin(9) * t509 - t380) * t434) * qJD(4) + ((pkin(4) * qJD(1) - t416 * t434) * t438 + t469 * pkin(9) * t437) * qJD(3);
t450 = -t383 * t502 + t485;
t301 = pkin(9) * t448 + t450;
t467 = -t433 * t298 - t436 * t301 - t323 * t497 + t330 * t498;
t466 = -t436 * t298 + t433 * t301 + t323 * t498 + t330 * t497;
t411 = t435 * t475;
t296 = t329 * t433 + t545;
t464 = pkin(4) * t498 - t296;
t293 = t323 * t433 + t545;
t459 = t342 * t436 - t352 * t433;
t455 = t435 * t470;
t454 = t471 * t435;
t401 = t414 * qJD(6);
t419 = qJ(6) * t475;
t283 = t419 + t401 - t467;
t359 = pkin(4) * t565 + t439 * t505;
t452 = t292 * t414 + t467;
t451 = t293 * t414 - t466;
t449 = t433 * t312 + t436 * t315 + t342 * t497 - t352 * t498;
t284 = t466 - t468;
t446 = qJD(3) * t462 + t543;
t334 = -pkin(4) * t448 + t416 * t505;
t429 = -pkin(4) * t436 - pkin(5);
t426 = pkin(4) * t433 + qJ(6);
t370 = t394 * t435;
t368 = t395 * t435;
t341 = pkin(5) * t394 - qJ(6) * t395 + t430;
t328 = pkin(5) * t369 + qJ(6) * t371 + t390;
t322 = -t498 * t536 + (t493 * t529 - t480) * t436 + t447 * t433;
t320 = t351 * t438 - t433 * t480 + t436 * t479;
t311 = -pkin(5) * t435 - t459;
t310 = qJ(6) * t435 + t557;
t308 = pkin(4) * t392 + t314;
t291 = qJ(6) * t414 + t293;
t290 = -pkin(5) * t414 + t496;
t288 = pkin(5) * t322 + qJ(6) * t320 + qJD(6) * t371 + t359;
t287 = -pkin(5) * t503 - t556;
t286 = qJ(6) * t503 + qJD(6) * t435 + t449;
t285 = t306 * pkin(5) + t305 * qJ(6) - qJD(6) * t457 + t334;
t1 = [(-t424 * t476 - t515 * t435 + (t542 + (qJD(1) * t513 + t534) * t434) * qJD(3)) * MDP(17) + (-t484 * t424 - t485 * t435 + ((t439 * t509 - t384) * t438 + t472 * t435) * t502 + (t392 * t532 + (-t514 * qJD(1) - t338) * t438 + (-t543 + (t439 * t469 + t399) * t438) * t437) * qJD(3)) * MDP(20) + ((-t400 * t502 + t374) * t424 + (-t439 * t515 + t384 * t500 + (qJD(1) * t388 - t424 * t535 + t337) * qJD(3)) * t438 + (t361 + (-t391 * t439 - t544) * qJD(3) + (-t437 * t472 - t537) * qJD(4)) * t435) * MDP(19) + (-t515 * t529 + (-0.2e1 * t391 * t434 - t392 * t437) * t499 + (-t437 * t391 + (t392 + 0.2e1 * t481) * t434) * t505) * MDP(15) + (t392 * t447 + t444 * t529) * MDP(14) + (-t438 * t527 + (-qJ(2) * t505 + qJD(2) * t438) * t492) * MDP(13) + (-t435 * t527 + (qJ(2) * t503 + qJD(2) * t435) * t492) * MDP(12) + (-t449 * t414 + t467 * t435 + t359 * t457 - t390 * t305 - t334 * t371 - t353 * t320 + (-qJD(1) * t557 - t293) * t503) * MDP(27) + ((-t424 - t510) * t477 + (t392 * t438 + (t432 * qJD(1) + (-t424 - t469) * t435) * t437) * qJD(3)) * MDP(16) + (-t305 * t435 - t320 * t414 + (-qJD(1) * t371 + t457) * t503) * MDP(23) + (-t306 * t435 - t322 * t414 + (-qJD(1) * t369 - t344) * t503) * MDP(24) + (t556 * t414 - t466 * t435 + t359 * t344 + t390 * t306 + t334 * t369 + t353 * t322 + (qJD(1) * t459 + t292) * t503) * MDP(26) + (t283 * t435 + t285 * t371 + t286 * t414 - t288 * t457 + t302 * t320 + t305 * t328 + (qJD(1) * t310 + t291) * t503) * MDP(30) + (-t284 * t435 + t285 * t369 - t287 * t414 + t288 * t344 + t302 * t322 + t306 * t328 + (-qJD(1) * t311 - t290) * t503) * MDP(28) + (t424 * t503 + t411) * MDP(18) + (t414 * t503 + t411) * MDP(25) + (-t283 * t369 - t284 * t371 - t286 * t344 + t287 * t457 - t290 * t320 - t291 * t322 - t305 * t311 - t306 * t310) * MDP(29) + (t305 * t369 + t306 * t371 + t320 * t344 - t322 * t457) * MDP(22) + (t305 * t371 - t320 * t457) * MDP(21) + (t283 * t310 + t284 * t311 + t285 * t328 + t286 * t291 + t287 * t290 + t288 * t302) * MDP(31) - 0.2e1 * MDP(7) * t411 + 0.2e1 * t495 * t561 + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t492 + (-MDP(10) * t438 - MDP(9) * t435) * t440; -t441 * MDP(5) - MDP(6) * t549 + (-t438 * t515 + (-qJD(1) - t501) * t539 + (-t435 * t391 - t486) * qJD(3)) * MDP(19) + (qJD(1) * t540 + (t392 * t435 - t424 * t529) * qJD(3) + (t488 - t542) * qJD(4)) * MDP(20) + (-t305 * t368 + t306 * t370 - t344 * t522 + t457 * t521) * MDP(29) + (-t283 * t370 + t284 * t368 - t285 * t438 + t290 * t521 + t291 * t522 + t302 * t505) * MDP(31) + (t435 * MDP(12) + t438 * MDP(13)) * (-t440 - t441) + (MDP(26) + MDP(28)) * (t344 * t505 + (-t368 * t495 - t306) * t438 - t521 * t414) + (-MDP(27) + MDP(30)) * (qJD(3) * (-t370 * t509 - t435 * t457) + t414 * t522 - t305 * t438); -t441 * t561 - qJ(2) * MDP(12) * t528 + (-t434 ^ 2 * t474 + (-t469 * t506 + t541) * t437) * MDP(14) + ((qJD(1) * t455 - t392 * qJD(4) - t515) * t434 + ((t391 + t504) * qJD(4) + (-t454 - t477) * qJD(1)) * t437) * MDP(15) + (t424 * t500 + (t424 * t533 + t438 * t470) * qJD(1)) * MDP(16) + (-t424 * t502 + (t438 * t471 - t488) * qJD(1)) * MDP(17) + (-pkin(3) * t515 - t382 * t424 + (-t454 + t486) * t416 + (-pkin(8) * t539 + t544) * qJD(4) + (-t337 * t438 + t434 * t446) * qJD(1)) * MDP(19) + (t516 * t424 + t416 * t455 + (pkin(8) * t540 + (t384 - t550) * t437) * qJD(4) + ((pkin(3) * t502 + t338) * t438 + t446 * t437) * qJD(1)) * MDP(20) + (-t305 * t395 - t457 * t518) * MDP(21) + (t305 * t394 - t306 * t395 + t344 * t518 - t457 * t517) * MDP(22) + (t430 * t306 + t334 * t394 + t465 * t344 + t517 * t353) * MDP(26) + (-t430 * t305 + t334 * t395 - t518 * t353 + t457 * t465) * MDP(27) + (t285 * t394 + t517 * t302 + t306 * t341 + t523 * t344) * MDP(28) + (-t283 * t394 + t284 * t395 - t290 * t518 - t291 * t517 + t305 * t456 - t306 * t358 + t344 * t525 - t457 * t524) * MDP(29) + (-t285 * t395 + t518 * t302 + t305 * t341 - t457 * t523) * MDP(30) + (t283 * t358 - t284 * t456 + t285 * t341 - t290 * t524 - t291 * t525 + t302 * t523) * MDP(31) + (-t518 * MDP(23) - t517 * MDP(24) + MDP(26) * t558 + MDP(27) * t559 + t524 * MDP(28) - t525 * MDP(30)) * t414 + (-t424 * MDP(18) + (qJD(3) * t395 - t457) * MDP(23) + (-qJD(3) * t394 + t344) * MDP(24) - t414 * MDP(25) + (-t292 + t508) * MDP(26) + (t293 - t507) * MDP(27) + (t290 + t508) * MDP(28) + (-t291 + t507) * MDP(30)) * t509 + (MDP(13) * t549 + MDP(7) * t528) * t435; -t392 * t391 * MDP(14) + (-t391 ^ 2 + t392 ^ 2) * MDP(15) + (-t391 * t424 + t444) * MDP(16) + (t448 + t541) * MDP(17) + MDP(18) * t475 + (-qJD(3) * t487 - t384 * t392 + t361 + (-qJD(4) + t424) * t338) * MDP(19) + (t337 * t424 - t384 * t391 - t450) * MDP(20) + (t296 * t414 - t560 + (-t344 * t392 - t414 * t498 + t436 * t475) * pkin(4) - t466) * MDP(26) + (t297 * t414 + t566 + (-t392 * t457 - t414 * t497 - t433 * t475) * pkin(4) + t467) * MDP(27) + (-t548 - t308 * t344 - t464 * t414 + (pkin(5) - t429) * t475 - t466) * MDP(28) + (-t305 * t429 - t306 * t426 + (t291 + t464) * t457 + (t290 + t526) * t344) * MDP(29) + (t308 * t457 - t414 * t526 + t426 * t475 + t283 - t567) * MDP(30) + (t283 * t426 + t284 * t429 + t290 * t464 - t291 * t526 - t302 * t308) * MDP(31) + t568; (t451 - t560) * MDP(26) + (t452 + t566) * MDP(27) + (-t314 * t344 + t451 + 0.2e1 * t468 - t548) * MDP(28) + (pkin(5) * t305 - qJ(6) * t306 + (t291 - t293) * t457 + (t290 - t496) * t344) * MDP(29) + (t314 * t457 + 0.2e1 * t401 + 0.2e1 * t419 - t452 - t567) * MDP(30) + (-pkin(5) * t284 + qJ(6) * t283 - t290 * t293 + t291 * t496 - t302 * t314) * MDP(31) + t568; (-t475 + t571) * MDP(28) + t289 * MDP(29) + (-t414 ^ 2 - t552) * MDP(30) + (-t291 * t414 + t284 + t548) * MDP(31);];
tauc  = t1;
