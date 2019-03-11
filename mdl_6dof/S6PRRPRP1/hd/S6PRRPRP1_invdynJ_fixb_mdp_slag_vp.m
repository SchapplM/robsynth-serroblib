% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:47
% EndTime: 2019-03-08 21:26:56
% DurationCPUTime: 7.11s
% Computational Cost: add. (4273->476), mult. (10020->642), div. (0->0), fcn. (7821->14), ass. (0->224)
t466 = sin(qJ(3));
t467 = sin(qJ(2));
t459 = sin(pkin(6));
t536 = qJD(1) * t459;
t516 = t467 * t536;
t576 = qJD(3) * pkin(3);
t594 = t466 * t576 - t516;
t469 = cos(qJ(3));
t578 = qJ(4) + pkin(8);
t507 = qJD(3) * t578;
t406 = qJD(4) * t469 - t466 * t507;
t407 = -qJD(4) * t466 - t469 * t507;
t457 = sin(pkin(11));
t460 = cos(pkin(11));
t350 = t406 * t460 + t407 * t457;
t558 = t460 * t469;
t421 = t457 * t466 - t558;
t470 = cos(qJ(2));
t515 = t470 * t536;
t384 = t421 * t515;
t544 = t350 + t384;
t422 = t457 * t469 + t460 * t466;
t413 = t422 * qJD(3);
t416 = t421 * qJD(3);
t593 = pkin(4) * t413 + pkin(9) * t416 + t594;
t414 = t422 * qJD(2);
t465 = sin(qJ(5));
t468 = cos(qJ(5));
t528 = t468 * qJD(3);
t530 = qJD(5) * t465;
t526 = qJD(2) * qJD(3);
t510 = t466 * t526;
t433 = t457 * t510;
t509 = t469 * t526;
t494 = t460 * t509 - t433;
t584 = qJDD(2) * t422 + t494;
t330 = -qJD(5) * t528 - t465 * qJDD(3) + t414 * t530 - t468 * t584;
t524 = qJDD(2) * t469;
t525 = qJDD(2) * t466;
t499 = -t457 * t525 + t460 * t524;
t363 = qJD(2) * t413 + qJDD(5) - t499;
t390 = qJD(3) * t465 + t414 * t468;
t462 = cos(pkin(6));
t523 = t462 * qJDD(1);
t437 = t469 * t523;
t527 = qJD(1) * qJD(2);
t397 = qJDD(2) * pkin(8) + (qJDD(1) * t467 + t470 * t527) * t459;
t535 = qJD(1) * t462;
t480 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t535 + t397;
t504 = qJD(2) * t578 + t516;
t496 = t504 * qJD(3);
t325 = qJDD(3) * pkin(3) - t466 * t480 - t469 * t496 + t437;
t326 = (-t496 + t523) * t466 + t480 * t469;
t307 = t457 * t325 + t460 * t326;
t305 = qJDD(3) * pkin(9) + t307;
t381 = -t466 * t504 + t469 * t535;
t376 = t381 + t576;
t382 = t466 * t535 + t469 * t504;
t559 = t460 * t382;
t334 = t457 * t376 + t559;
t329 = qJD(3) * pkin(9) + t334;
t448 = t469 * pkin(3) + pkin(2);
t404 = -qJD(2) * t448 + qJD(4) - t515;
t532 = qJD(2) * t466;
t412 = qJD(2) * t558 - t457 * t532;
t341 = -pkin(4) * t412 - pkin(9) * t414 + t404;
t315 = t329 * t468 + t341 * t465;
t367 = -qJD(3) * t414 + t499;
t560 = t459 * t470;
t508 = qJDD(1) * t560;
t511 = t467 * t527;
t500 = t459 * t511 - t508;
t485 = pkin(3) * t510 + qJDD(4) + t500;
t587 = -pkin(9) * t422 - t448;
t318 = -t367 * pkin(4) - t494 * pkin(9) + qJDD(2) * t587 + t485;
t317 = t468 * t318;
t479 = -qJD(5) * t315 - t465 * t305 + t317;
t297 = pkin(5) * t363 + qJ(6) * t330 - qJD(6) * t390 + t479;
t388 = t414 * t465 - t528;
t309 = -qJ(6) * t388 + t315;
t405 = qJD(5) - t412;
t592 = t309 * t405 + t297;
t331 = qJD(5) * t390 - t468 * qJDD(3) + t465 * t584;
t529 = qJD(5) * t468;
t488 = -t468 * t305 - t465 * t318 + t329 * t530 - t341 * t529;
t298 = -qJ(6) * t331 - qJD(6) * t388 - t488;
t314 = -t329 * t465 + t468 * t341;
t308 = -qJ(6) * t390 + t314;
t302 = pkin(5) * t405 + t308;
t591 = -t302 * t405 + t298;
t588 = -t384 * t465 + t593 * t468;
t458 = sin(pkin(10));
t461 = cos(pkin(10));
t555 = t462 * t470;
t408 = t458 * t467 - t461 * t555;
t410 = t458 * t555 + t461 * t467;
t503 = g(1) * t410 + g(2) * t408;
t484 = -g(3) * t560 + t503;
t365 = pkin(4) * t421 + t587;
t586 = t365 * t529 + t593 * t465 + t544 * t468;
t471 = qJD(3) ^ 2;
t585 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t471 + t459 * (-g(3) * t470 + t511) - t500 + t503;
t583 = t390 ^ 2;
t582 = pkin(3) * t460;
t581 = pkin(5) * t465;
t579 = g(3) * t459;
t577 = qJD(2) * pkin(2);
t574 = t330 * t465;
t573 = t388 * t412;
t572 = t388 * t414;
t571 = t390 * t405;
t570 = t390 * t414;
t569 = t405 * t465;
t568 = t412 * t465;
t567 = t422 * t465;
t566 = t422 * t468;
t454 = qJ(3) + pkin(11);
t450 = cos(t454);
t565 = t450 * t465;
t372 = t457 * t382;
t564 = t458 * t459;
t563 = t459 * t461;
t562 = t459 * t467;
t561 = t459 * t469;
t557 = t462 * t467;
t556 = t462 * t469;
t554 = t465 * t363;
t553 = t465 * t470;
t358 = t468 * t363;
t429 = t578 * t466;
t430 = t578 * t469;
t386 = -t429 * t457 + t430 * t460;
t377 = t468 * t386;
t444 = pkin(3) * t457 + pkin(9);
t552 = qJ(6) + t444;
t551 = qJDD(1) - g(3);
t497 = qJ(6) * t416 - qJD(6) * t422;
t550 = pkin(5) * t413 - t350 * t465 + t497 * t468 + (-t377 + (qJ(6) * t422 - t365) * t465) * qJD(5) + t588;
t512 = t422 * t529;
t549 = -qJ(6) * t512 + (-qJD(5) * t386 + t497) * t465 + t586;
t548 = -t308 + t302;
t547 = -t465 * t331 - t388 * t529;
t338 = t381 * t460 - t372;
t351 = pkin(3) * t532 + pkin(4) * t414 - pkin(9) * t412;
t546 = t468 * t338 + t465 * t351;
t545 = t406 * t457 - t460 * t407 - t422 * t515;
t543 = t405 * t568 + t358;
t542 = t465 * t365 + t377;
t409 = t458 * t470 + t461 * t557;
t541 = -t408 * t448 + t409 * t578;
t411 = -t458 * t557 + t461 * t470;
t540 = -t410 * t448 + t411 * t578;
t506 = qJD(5) * t552;
t539 = qJ(6) * t568 + qJD(6) * t468 - t465 * t506 - t546;
t347 = t468 * t351;
t538 = -pkin(5) * t414 - t347 + (qJ(6) * t412 - t506) * t468 + (-qJD(6) + t338) * t465;
t455 = t466 ^ 2;
t537 = -t469 ^ 2 + t455;
t534 = qJD(2) * t404;
t533 = qJD(2) * t459;
t531 = qJD(5) * t405;
t521 = t461 * t561;
t520 = t459 * t553;
t519 = t466 * t562;
t518 = t468 * t560;
t447 = pkin(5) * t468 + pkin(4);
t514 = t467 * t533;
t513 = t470 * t533;
t306 = t325 * t460 - t457 * t326;
t333 = t376 * t460 - t372;
t336 = t381 * t457 + t559;
t385 = t460 * t429 + t430 * t457;
t505 = t405 * t468;
t502 = g(1) * t411 + g(2) * t409;
t501 = g(1) * t458 - g(2) * t461;
t449 = sin(t454);
t463 = -qJ(6) - pkin(9);
t498 = t447 * t450 - t449 * t463;
t328 = -qJD(3) * pkin(4) - t333;
t417 = -t519 + t556;
t418 = t462 * t466 + t467 * t561;
t357 = t417 * t457 + t418 * t460;
t342 = -t357 * t465 - t518;
t492 = -t357 * t468 + t520;
t304 = -qJDD(3) * pkin(4) - t306;
t491 = -t416 * t465 + t512;
t490 = -t416 * t468 - t422 * t530;
t489 = -t509 - t525;
t487 = t328 * t405 - t444 * t363;
t368 = t409 * t449 + t450 * t563;
t370 = t411 * t449 - t450 * t564;
t399 = t449 * t562 - t462 * t450;
t486 = g(1) * t370 + g(2) * t368 + g(3) * t399;
t483 = -g(3) * t562 - t502;
t427 = -t515 - t577;
t482 = -qJD(2) * t427 - t397 + t502;
t299 = pkin(5) * t331 + qJDD(6) + t304;
t478 = t444 * t531 + t304 - t486;
t477 = -pkin(8) * qJDD(3) + (t427 + t515 - t577) * qJD(3);
t366 = -qJDD(2) * t448 + t485;
t369 = t409 * t450 - t449 * t563;
t371 = t411 * t450 + t449 * t564;
t400 = t449 * t462 + t450 * t562;
t473 = -g(1) * (-t371 * t465 + t410 * t468) - g(2) * (-t369 * t465 + t408 * t468) - g(3) * (-t400 * t465 - t518);
t472 = qJD(2) ^ 2;
t445 = -pkin(4) - t582;
t440 = pkin(3) * t556;
t428 = t458 * pkin(3) * t561;
t424 = t448 * t560;
t420 = t552 * t468;
t419 = t552 * t465;
t387 = t388 ^ 2;
t380 = -qJD(3) * t418 - t466 * t513;
t379 = qJD(3) * t417 + t469 * t513;
t360 = t468 * t365;
t356 = -t460 * t417 + t418 * t457;
t337 = t379 * t460 + t380 * t457;
t335 = t379 * t457 - t460 * t380;
t321 = -qJ(6) * t567 + t542;
t320 = pkin(5) * t388 + qJD(6) + t328;
t319 = pkin(5) * t421 - qJ(6) * t566 - t386 * t465 + t360;
t313 = qJD(5) * t492 - t337 * t465 + t468 * t514;
t312 = qJD(5) * t342 + t337 * t468 + t465 * t514;
t1 = [t551 * MDP(1) + (qJD(3) * t380 + qJDD(3) * t417) * MDP(10) + (-qJD(3) * t379 - qJDD(3) * t418) * MDP(11) + (t335 * t414 + t337 * t412 + t356 * t584 + t357 * t367) * MDP(12) + (-t306 * t356 + t307 * t357 - t333 * t335 + t334 * t337 - g(3)) * MDP(13) + (t313 * t405 + t331 * t356 + t335 * t388 + t342 * t363) * MDP(19) + (-t312 * t405 - t330 * t356 + t335 * t390 + t363 * t492) * MDP(20) + (-t312 * t388 - t313 * t390 + t330 * t342 + t331 * t492) * MDP(21) + (t297 * t342 - t298 * t492 + t299 * t356 + t302 * t313 + t309 * t312 + t320 * t335 - g(3)) * MDP(22) + ((MDP(13) * t534 - qJDD(2) * MDP(4) + (-t469 * MDP(10) + t466 * MDP(11) - MDP(3)) * t472) * t467 + (qJDD(2) * MDP(3) - t472 * MDP(4) + (-t510 + t524) * MDP(10) + t489 * MDP(11) - t366 * MDP(13)) * t470) * t459; qJDD(2) * MDP(2) + (t484 + t508) * MDP(3) + (-t551 * t562 + t502) * MDP(4) + (qJDD(2) * t455 + 0.2e1 * t466 * t509) * MDP(5) + 0.2e1 * (t466 * t524 - t526 * t537) * MDP(6) + (qJDD(3) * t466 + t469 * t471) * MDP(7) + (qJDD(3) * t469 - t466 * t471) * MDP(8) + (t477 * t466 + t469 * t585) * MDP(10) + (-t466 * t585 + t477 * t469) * MDP(11) + (-t306 * t422 - t307 * t421 + t333 * t416 - t334 * t413 + t386 * t367 + t385 * t584 + t412 * t544 + t414 * t545 + t483) * MDP(12) + (t307 * t386 - t306 * t385 - t366 * t448 - g(1) * t540 - g(2) * t541 - g(3) * (t562 * t578 + t424) + t594 * t404 + t544 * t334 - t545 * t333) * MDP(13) + (-t330 * t566 + t390 * t490) * MDP(14) + (-(-t388 * t468 - t390 * t465) * t416 + (t574 - t331 * t468 + (t388 * t465 - t390 * t468) * qJD(5)) * t422) * MDP(15) + (-t330 * t421 + t358 * t422 + t390 * t413 + t405 * t490) * MDP(16) + (-t331 * t421 - t388 * t413 - t405 * t491 - t422 * t554) * MDP(17) + (t363 * t421 + t405 * t413) * MDP(18) + (t314 * t413 + t317 * t421 + t385 * t331 + t360 * t363 + t588 * t405 + t545 * t388 + (t484 * t450 + (t328 * t422 - t329 * t421 - t386 * t405) * qJD(5)) * t468 + ((-qJD(5) * t365 - t350) * t405 - t386 * t363 + (-qJD(5) * t341 - t305) * t421 + t304 * t422 - t328 * t416 + t483) * t465) * MDP(19) + (-t542 * t363 + t488 * t421 - t315 * t413 - t385 * t330 + t304 * t566 - g(1) * (t410 * t565 + t411 * t468) - g(2) * (t408 * t565 + t409 * t468) - (-t450 * t553 + t467 * t468) * t579 + (t386 * t530 - t586) * t405 + t545 * t390 + t490 * t328) * MDP(20) + (t319 * t330 - t321 * t331 - (-t302 * t468 - t309 * t465) * t416 - t550 * t390 - t549 * t388 + t484 * t449 + (-t297 * t468 - t298 * t465 + (t302 * t465 - t309 * t468) * qJD(5)) * t422) * MDP(21) + (t298 * t321 + t297 * t319 + t299 * (pkin(5) * t567 + t385) - g(1) * (-t410 * t498 + t411 * t581 + t540) - g(2) * (-t408 * t498 + t409 * t581 + t541) - g(3) * t424 + (pkin(5) * t491 + t545) * t320 + t549 * t309 + t550 * t302 - (t498 * t470 + (t578 + t581) * t467) * t579) * MDP(22); MDP(7) * t525 + MDP(8) * t524 + qJDD(3) * MDP(9) + (-g(3) * t417 + t466 * t482 - t501 * t561 + t437) * MDP(10) + (g(3) * t418 + (t459 * t501 - t523) * t466 + t482 * t469) * MDP(11) + ((t334 - t336) * t414 + (-t338 + t333) * t412 + (t457 * t367 + (-t457 * t524 + t460 * t489 + t433) * t460) * pkin(3)) * MDP(12) + (-g(1) * t428 - g(3) * t440 + t333 * t336 - t334 * t338 + (g(2) * t521 + t306 * t460 + t307 * t457 + (-t483 - t534) * t466) * pkin(3)) * MDP(13) + (t390 * t505 - t574) * MDP(14) + ((-t330 + t573) * t468 - t390 * t569 + t547) * MDP(15) + (t405 * t505 + t554 - t570) * MDP(16) + (-t405 * t530 + t543 + t572) * MDP(17) - t405 * t414 * MDP(18) + (-t314 * t414 + t445 * t331 - t336 * t388 - t347 * t405 + (t338 * t405 + t487) * t465 - t478 * t468) * MDP(19) + (t315 * t414 - t445 * t330 - t336 * t390 + t405 * t546 + t465 * t478 + t468 * t487) * MDP(20) + (-g(1) * t371 - g(2) * t369 - g(3) * t400 - t330 * t419 - t331 * t420 - t539 * t388 - t538 * t390 - t465 * t592 + t468 * t591) * MDP(21) + (t298 * t420 - t297 * t419 + t299 * (-t447 - t582) - g(1) * (-pkin(3) * t411 * t466 - t370 * t447 - t371 * t463 + t428) - g(2) * (-t368 * t447 - t369 * t463 + (-t409 * t466 - t521) * pkin(3)) - g(3) * (-pkin(3) * t519 - t399 * t447 - t400 * t463 + t440) + (pkin(5) * t569 - t336) * t320 + t539 * t309 + t538 * t302) * MDP(22) + (-MDP(5) * t466 * t469 + MDP(6) * t537) * t472; (-t412 ^ 2 - t414 ^ 2) * MDP(12) + (t333 * t414 - t334 * t412 + t366 - t484) * MDP(13) + (t543 - t572) * MDP(19) - MDP(20) * t570 + t547 * MDP(21) + (-t320 * t414 - t484) * MDP(22) + (-MDP(19) * t531 - t363 * MDP(20) + MDP(21) * t571 + MDP(22) * t591) * t465 + ((t330 + t573) * MDP(21) + t592 * MDP(22) - t405 ^ 2 * MDP(20)) * t468; t390 * t388 * MDP(14) + (-t387 + t583) * MDP(15) + (t388 * t405 - t330) * MDP(16) + (-t331 + t571) * MDP(17) + t363 * MDP(18) + (t315 * t405 - t328 * t390 + t473 + t479) * MDP(19) + (t314 * t405 + t328 * t388 - g(1) * (-t371 * t468 - t410 * t465) - g(2) * (-t369 * t468 - t408 * t465) - g(3) * (-t400 * t468 + t520) + t488) * MDP(20) + (pkin(5) * t330 - t388 * t548) * MDP(21) + (t548 * t309 + (-t320 * t390 + t297 + t473) * pkin(5)) * MDP(22); (-t387 - t583) * MDP(21) + (t302 * t390 + t309 * t388 + t299 - t486) * MDP(22);];
tau  = t1;
