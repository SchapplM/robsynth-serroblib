% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:59
% EndTime: 2019-03-09 02:27:07
% DurationCPUTime: 5.90s
% Computational Cost: add. (3263->499), mult. (6095->660), div. (0->0), fcn. (4110->12), ass. (0->223)
t469 = cos(qJ(4));
t470 = -pkin(1) - pkin(2);
t427 = qJD(1) * t470 + qJD(2);
t461 = sin(pkin(10));
t462 = cos(pkin(10));
t554 = qJ(2) * qJD(1);
t386 = t461 * t427 + t462 * t554;
t376 = -qJD(1) * pkin(7) + t386;
t466 = sin(qJ(4));
t546 = qJD(4) * t469;
t539 = qJD(1) * qJD(2);
t432 = t462 * t539;
t424 = qJDD(1) * t470 + qJDD(2);
t534 = qJDD(1) * t462;
t560 = qJ(2) * t534 + t461 * t424;
t368 = t432 + t560;
t363 = -qJDD(1) * pkin(7) + t368;
t606 = -qJD(3) * qJD(4) - t363;
t525 = -t376 * t546 + t466 * t606;
t320 = -qJDD(4) * pkin(4) - qJDD(3) * t469 - t525;
t552 = qJD(1) * t469;
t428 = qJD(5) + t552;
t592 = sin(qJ(1));
t593 = cos(qJ(1));
t402 = t461 * t593 - t462 * t592;
t588 = g(2) * t402;
t401 = -t461 * t592 - t462 * t593;
t591 = g(1) * t401;
t498 = t588 + t591;
t477 = -g(3) * t469 + t466 * t498;
t608 = qJD(5) * pkin(8) * t428 + t320 + t477;
t465 = sin(qJ(5));
t468 = cos(qJ(5));
t547 = qJD(4) * t468;
t553 = qJD(1) * t466;
t403 = t465 * t553 + t547;
t467 = cos(qJ(6));
t541 = t465 * qJD(4);
t404 = t468 * t553 - t541;
t464 = sin(qJ(6));
t578 = t404 * t464;
t355 = t467 * t403 + t578;
t426 = qJD(6) + t428;
t607 = t355 * t426;
t358 = qJD(3) * t469 - t466 * t376;
t605 = qJD(4) * t358;
t359 = t466 * qJD(3) + t469 * t376;
t350 = qJD(4) * pkin(8) + t359;
t385 = t427 * t462 - t461 * t554;
t375 = qJD(1) * pkin(3) - t385;
t501 = pkin(4) * t469 + pkin(8) * t466;
t351 = qJD(1) * t501 + t375;
t322 = t468 * t350 + t465 * t351;
t317 = pkin(9) * t403 + t322;
t543 = qJD(6) * t464;
t315 = t317 * t543;
t349 = -qJD(4) * pkin(4) - t358;
t332 = -pkin(5) * t403 + t349;
t459 = qJ(5) + qJ(6);
t450 = sin(t459);
t451 = cos(t459);
t575 = t451 * t469;
t335 = t401 * t450 + t402 * t575;
t337 = -t401 * t575 + t402 * t450;
t587 = g(3) * t466;
t604 = g(1) * t337 - g(2) * t335 - t332 * t355 - t451 * t587 + t315;
t447 = t469 * qJDD(1);
t538 = qJD(1) * qJD(4);
t513 = t466 * t538;
t599 = -t513 + t447;
t400 = -qJDD(5) - t599;
t399 = -qJDD(6) + t400;
t490 = t403 * t464 - t467 * t404;
t603 = (-t355 ^ 2 + t490 ^ 2) * MDP(25) - t399 * MDP(28) - t355 * MDP(24) * t490;
t601 = t426 * t490;
t545 = qJD(5) * t465;
t514 = t466 * t545;
t598 = t468 * t546 - t514;
t597 = qJD(5) + qJD(6);
t576 = t450 * t469;
t334 = -t401 * t451 + t402 * t576;
t336 = t401 * t576 + t402 * t451;
t319 = qJDD(4) * pkin(8) + qJDD(3) * t466 + t363 * t469 + t605;
t430 = t461 * t539;
t535 = qJDD(1) * t461;
t508 = -qJ(2) * t535 + t424 * t462;
t367 = -t430 + t508;
t362 = qJDD(1) * pkin(3) - t367;
t512 = t469 * t538;
t533 = qJDD(1) * t466;
t485 = t512 + t533;
t329 = pkin(4) * t599 + t485 * pkin(8) + t362;
t328 = t468 * t329;
t536 = qJD(4) * qJD(5);
t343 = qJD(1) * t514 + t465 * qJDD(4) + (-t485 + t536) * t468;
t305 = -pkin(5) * t400 - pkin(9) * t343 - qJD(5) * t322 - t465 * t319 + t328;
t517 = t469 * t541;
t544 = qJD(5) * t468;
t480 = t544 * t466 + t517;
t344 = qJD(1) * t480 + qJDD(4) * t468 + (t533 - t536) * t465;
t527 = -t468 * t319 - t465 * t329 - t351 * t544;
t484 = t350 * t545 + t527;
t306 = pkin(9) * t344 - t484;
t511 = t467 * t305 - t464 * t306;
t596 = -g(1) * t336 - g(2) * t334 - t332 * t490 - t450 * t587 + t511;
t414 = t462 * qJ(2) + t461 * t470;
t406 = -pkin(7) + t414;
t590 = g(1) * t402;
t595 = qJD(5) * (t406 * t428 + t350) + t590;
t510 = t343 * t464 - t467 * t344;
t311 = qJD(6) * t490 + t510;
t594 = pkin(8) + pkin(9);
t589 = g(2) * t401;
t585 = pkin(1) * qJDD(1);
t321 = -t350 * t465 + t468 * t351;
t316 = pkin(9) * t404 + t321;
t314 = pkin(5) * t428 + t316;
t584 = t314 * t467;
t583 = t317 * t467;
t582 = t343 * t465;
t581 = t401 * t465;
t580 = t403 * t428;
t579 = t404 * t428;
t577 = t428 * t468;
t574 = t461 * t466;
t573 = t465 * t466;
t572 = t465 * t469;
t571 = t466 * t468;
t570 = t467 * t468;
t569 = t468 * t469;
t568 = qJDD(3) + g(3);
t500 = -pkin(4) * t466 + pkin(8) * t469;
t410 = t500 * qJD(1);
t567 = t468 * t358 + t465 * t410;
t407 = t464 * t465 - t570;
t486 = t407 * t469;
t566 = -qJD(1) * t486 - t407 * t597;
t408 = t464 * t468 + t465 * t467;
t361 = t597 * t408;
t565 = t408 * t552 + t361;
t395 = -t461 * t572 - t462 * t468;
t516 = t466 * t547;
t564 = -qJD(5) * t395 + t461 * t516 + (t461 * t465 + t462 * t569) * qJD(1);
t396 = t461 * t569 - t462 * t465;
t518 = t466 * t541;
t563 = -qJD(5) * t396 + t461 * t518 - (t461 * t468 - t462 * t572) * qJD(1);
t391 = qJD(2) * t461 + qJD(4) * t500;
t562 = t468 * t391 + t406 * t518;
t413 = -t461 * qJ(2) + t462 * t470;
t405 = pkin(3) - t413;
t377 = t405 + t501;
t387 = t406 * t569;
t561 = t465 * t377 + t387;
t559 = t593 * pkin(1) + t592 * qJ(2);
t558 = g(1) * t592 - g(2) * t593;
t457 = t466 ^ 2;
t557 = -t469 ^ 2 + t457;
t471 = qJD(4) ^ 2;
t472 = qJD(1) ^ 2;
t556 = t471 + t472;
t551 = qJD(2) * t462;
t549 = qJD(4) * t403;
t548 = qJD(4) * t466;
t542 = qJD(6) * t467;
t540 = qJ(2) * qJDD(1);
t532 = qJDD(4) * t466;
t531 = qJDD(4) * t469;
t529 = 0.2e1 * t539;
t528 = t464 * t573;
t526 = t467 * t343 + t464 * t344 + t403 * t542;
t521 = t469 * t551;
t524 = t377 * t544 + t465 * t391 + t468 * t521;
t523 = qJD(5) * t594;
t522 = t465 * t552;
t520 = t428 * t541;
t519 = t428 * t547;
t507 = -qJD(5) * t351 - t319;
t506 = qJD(6) * t314 + t306;
t505 = 0.2e1 * t512;
t504 = qJDD(2) - t585;
t503 = -pkin(1) * t592 + t593 * qJ(2);
t502 = -t359 + (t522 + t545) * pkin(5);
t499 = t589 - t590;
t398 = t468 * t410;
t422 = t594 * t468;
t489 = -pkin(5) * t466 + pkin(9) * t569;
t497 = qJD(1) * t489 + qJD(6) * t422 - t358 * t465 + t468 * t523 + t398;
t421 = t594 * t465;
t496 = pkin(9) * t522 + qJD(6) * t421 + t465 * t523 + t567;
t495 = -qJD(6) * t395 + t564;
t494 = qJD(6) * t396 - t563;
t308 = t314 * t464 + t583;
t324 = qJD(4) * t486 + t361 * t466;
t389 = t466 * t570 - t528;
t493 = t324 * t426 + t389 * t399;
t325 = -qJD(6) * t528 + (t571 * t597 + t517) * t467 + t598 * t464;
t388 = t408 * t466;
t492 = -t325 * t426 + t388 * t399;
t491 = t385 * t461 - t386 * t462;
t488 = -t465 * t400 + t428 * t544;
t487 = t400 * t468 + t428 * t545;
t310 = t404 * t543 + t526;
t483 = g(1) * t593 + g(2) * t592;
t482 = qJD(1) * t375 - t498;
t481 = -t461 * t546 + t462 * t553;
t479 = pkin(8) * t400 + t349 * t428;
t474 = -qJDD(4) * t406 + (-qJD(1) * t405 - t375 - t551) * qJD(4);
t473 = qJDD(1) * t405 - t406 * t471 + t362 + t430 + t499;
t443 = -pkin(5) * t468 - pkin(4);
t417 = -t466 * t471 + t531;
t416 = -t469 * t471 - t532;
t374 = (-pkin(5) * t465 + t406) * t466;
t373 = t468 * t377;
t346 = -t401 * t569 + t402 * t465;
t345 = t401 * t572 + t402 * t468;
t333 = -pkin(5) * t480 + t406 * t546 + t551 * t466;
t331 = pkin(9) * t573 + t561;
t330 = pkin(9) * t571 + t373 + (-t406 * t465 + pkin(5)) * t469;
t313 = -pkin(5) * t344 + t320;
t312 = (-t469 * t545 - t516) * t406 + t480 * pkin(9) + t524;
t309 = -t465 * t521 + t489 * qJD(4) + (-t387 + (-pkin(9) * t466 - t377) * t465) * qJD(5) + t562;
t307 = -t317 * t464 + t584;
t1 = [(t310 * t388 + t311 * t389 + t324 * t355 + t325 * t490) * MDP(25) + ((t309 * t467 - t312 * t464) * t426 - (t330 * t467 - t331 * t464) * t399 + t511 * t469 - t307 * t548 - t333 * t355 + t374 * t311 - t313 * t388 - t332 * t325 - g(1) * t335 - g(2) * t337 + ((-t330 * t464 - t331 * t467) * t426 - t308 * t469) * qJD(6)) * MDP(29) + (-t311 * t469 - t355 * t548 - t492) * MDP(27) + (-t310 * t389 + t324 * t490) * MDP(24) + (t310 * t469 - t490 * t548 + t493) * MDP(26) + (t315 * t469 + t308 * t548 + t333 * t490 + t374 * t310 - t313 * t389 + t332 * t324 + g(1) * t334 - g(2) * t336 + (-(-qJD(6) * t331 + t309) * t426 + t330 * t399 - t305 * t469) * t464 + (-(qJD(6) * t330 + t312) * t426 + t331 * t399 - t506 * t469) * t467) * MDP(30) + t416 * MDP(12) - t417 * MDP(13) + (-t343 * t571 + t404 * t598) * MDP(17) + (-t524 * t428 + t561 * t400 - t468 * t591 - g(2) * t345 + ((-t349 * t468 - t404 * t406) * qJD(4) + t595 * t465 + t527) * t469 + (-t404 * t551 + t349 * t545 - t320 * t468 + t406 * t343 + (t406 * t577 + t322) * qJD(4)) * t466) * MDP(23) + ((-t377 * t545 + t562) * t428 - t373 * t400 - g(1) * t581 - g(2) * t346 + (-t406 * t549 + t328 - t595 * t468 + (-t349 * qJD(4) + t406 * t400 - t428 * t551 + t507) * t465) * t469 + (-qJD(4) * t321 - t320 * t465 - t344 * t406 - t349 * t544 - t403 * t551) * t466) * MDP(22) + (-t483 + t529 + 0.2e1 * t540) * MDP(5) + ((t343 - t519) * t469 + (qJD(4) * t404 + t487) * t466) * MDP(19) + qJDD(1) * MDP(1) + (-qJDD(1) * t413 + 0.2e1 * t430 + t499 - t508) * MDP(7) + (qJDD(1) * t457 + t466 * t505) * MDP(10) + 0.2e1 * (t447 * t466 - t538 * t557) * MDP(11) + t558 * MDP(2) + (-t504 * pkin(1) - g(1) * t503 - g(2) * t559 + (t529 + t540) * qJ(2)) * MDP(6) + (qJDD(1) * t414 + 0.2e1 * t432 + t498 + t560) * MDP(8) + t483 * MDP(3) + (t466 * t474 + t469 * t473) * MDP(15) + (-t466 * t473 + t469 * t474) * MDP(16) + (-t400 * t469 - t428 * t548) * MDP(21) + (-t399 * t469 - t426 * t548) * MDP(28) + ((t344 + t520) * t469 + (t488 - t549) * t466) * MDP(20) + ((-t403 * t468 - t404 * t465) * t546 + (t582 - t344 * t468 + (t403 * t465 - t404 * t468) * qJD(5)) * t466) * MDP(18) + (-qJDD(2) + t558 + 0.2e1 * t585) * MDP(4) + (t368 * t414 + t367 * t413 - g(1) * (-pkin(2) * t592 + t503) - g(2) * (pkin(2) * t593 + t559) - t491 * qJD(2)) * MDP(9); -qJDD(1) * MDP(4) - t472 * MDP(5) + (-qJ(2) * t472 + t504 - t558) * MDP(6) + (-t461 * t472 - t534) * MDP(7) + (-t462 * t472 + t535) * MDP(8) + (qJD(1) * t491 + t367 * t462 + t368 * t461 - t558) * MDP(9) + ((0.2e1 * t513 - t447) * t462 + (-t469 * t556 - t532) * t461) * MDP(15) + ((t505 + t533) * t462 + (t466 * t556 - t531) * t461) * MDP(16) + (-t344 * t574 - t395 * t400 + t403 * t481 + t428 * t563) * MDP(22) + (t343 * t574 + t396 * t400 + t404 * t481 + t428 * t564) * MDP(23) + (-(t395 * t467 - t396 * t464) * t399 + t311 * t574 + (t464 * t495 - t467 * t494) * t426 + t481 * t355) * MDP(29) + ((t395 * t464 + t396 * t467) * t399 + t310 * t574 + (t464 * t494 + t467 * t495) * t426 - t481 * t490) * MDP(30); t568 * MDP(9) + t417 * MDP(15) + t416 * MDP(16) + t492 * MDP(29) + t493 * MDP(30) + ((t344 - t520) * MDP(22) + (-t343 - t519) * MDP(23) - t311 * MDP(29) - t310 * MDP(30)) * t469 + (-t488 * MDP(22) + t487 * MDP(23) + (-t403 * MDP(22) - t404 * MDP(23) - MDP(29) * t355 + MDP(30) * t490) * qJD(4)) * t466; -MDP(12) * t533 - MDP(13) * t447 + qJDD(4) * MDP(14) + (qJD(4) * t359 + t466 * t482 + t469 * t568 + t525) * MDP(15) + (t605 + (qJD(4) * t376 - t568) * t466 + (t482 + t606) * t469) * MDP(16) + (-t404 * t577 + t582) * MDP(17) + ((t343 + t580) * t468 + (t344 + t579) * t465) * MDP(18) + ((-t404 * t466 + t428 * t569) * qJD(1) + t488) * MDP(19) + ((t403 * t466 - t428 * t572) * qJD(1) - t487) * MDP(20) + (pkin(4) * t344 + t359 * t403 - t398 * t428 + (t358 * t428 + t479) * t465 - t608 * t468) * MDP(22) + (-pkin(4) * t343 + t359 * t404 + t567 * t428 + t608 * t465 + t479 * t468) * MDP(23) + (t310 * t408 + t490 * t566) * MDP(24) + (-t310 * t407 - t311 * t408 + t355 * t566 - t490 * t565) * MDP(25) + (-t399 * t408 + t426 * t566) * MDP(26) + (t399 * t407 - t426 * t565) * MDP(27) + (-(-t421 * t467 - t422 * t464) * t399 + t443 * t311 + t313 * t407 + (t464 * t496 - t467 * t497) * t426 - t502 * t355 + t565 * t332 - t477 * t451) * MDP(29) + ((-t421 * t464 + t422 * t467) * t399 + t443 * t310 + t313 * t408 + (t464 * t497 + t467 * t496) * t426 + t502 * t490 + t566 * t332 + t477 * t450) * MDP(30) + (t428 * MDP(21) + t321 * MDP(22) - MDP(23) * t322 + MDP(26) * t490 + MDP(27) * t355 + t426 * MDP(28) + t307 * MDP(29) - t308 * MDP(30)) * t553 + (-MDP(10) * t466 * t469 + MDP(11) * t557) * t472; t404 * t403 * MDP(17) + (-t403 ^ 2 + t404 ^ 2) * MDP(18) + (t343 - t580) * MDP(19) + (t344 - t579) * MDP(20) - t400 * MDP(21) + (-g(1) * t345 + t322 * t428 + t349 * t404 + t328 + (-qJD(5) * t350 + t589) * t468 + (-t469 * t588 + t507 - t587) * t465) * MDP(22) + (t321 * t428 - t349 * t403 + g(1) * t346 - g(2) * (t402 * t569 + t581) - g(3) * t571 + t484) * MDP(23) + (t310 - t607) * MDP(26) + (-t311 + t601) * MDP(27) + (-(-t316 * t464 - t583) * t426 - t308 * qJD(6) + (-t355 * t404 - t399 * t467 - t426 * t543) * pkin(5) + t596) * MDP(29) + ((-t317 * t426 - t305) * t464 + (t316 * t426 - t506) * t467 + (t399 * t464 + t404 * t490 - t426 * t542) * pkin(5) + t604) * MDP(30) + t603; (t526 - t607) * MDP(26) + (-t510 + t601) * MDP(27) + (t308 * t426 + t596) * MDP(29) + (-t464 * t305 - t467 * t306 + t307 * t426 + t604) * MDP(30) + (MDP(26) * t578 - MDP(27) * t490 - MDP(29) * t308 - MDP(30) * t584) * qJD(6) + t603;];
tau  = t1;
