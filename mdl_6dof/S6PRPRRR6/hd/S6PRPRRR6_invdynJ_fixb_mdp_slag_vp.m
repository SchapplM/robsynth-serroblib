% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:59
% EndTime: 2019-03-08 20:49:07
% DurationCPUTime: 6.28s
% Computational Cost: add. (2905->517), mult. (6073->709), div. (0->0), fcn. (4657->14), ass. (0->228)
t463 = sin(qJ(5));
t467 = cos(qJ(5));
t548 = t467 * qJD(4);
t468 = cos(qJ(4));
t563 = qJD(2) * t468;
t417 = t463 * t563 - t548;
t466 = cos(qJ(6));
t560 = qJD(4) * t463;
t419 = t467 * t563 + t560;
t462 = sin(qJ(6));
t595 = t419 * t462;
t354 = t466 * t417 + t595;
t464 = sin(qJ(4));
t564 = qJD(2) * t464;
t446 = qJD(5) + t564;
t440 = qJD(6) + t446;
t623 = t354 * t440;
t497 = t417 * t462 - t466 * t419;
t622 = t440 * t497;
t542 = qJDD(2) * t468;
t547 = qJD(2) * qJD(4);
t621 = -t464 * t547 + t542;
t470 = -pkin(2) - pkin(8);
t469 = cos(qJ(2));
t460 = sin(pkin(6));
t569 = qJD(1) * t460;
t532 = t469 * t569;
t501 = qJD(3) - t532;
t406 = qJD(2) * t470 + t501;
t461 = cos(pkin(6));
t568 = qJD(1) * t468;
t443 = t461 * t568;
t376 = t464 * t406 + t443;
t362 = qJD(4) * pkin(9) + t376;
t426 = pkin(4) * t464 - pkin(9) * t468 + qJ(3);
t465 = sin(qJ(2));
t533 = t465 * t569;
t383 = qJD(2) * t426 + t533;
t333 = t362 * t467 + t383 * t463;
t324 = -pkin(10) * t417 + t333;
t551 = qJD(6) * t462;
t322 = t324 * t551;
t586 = t461 * t464;
t615 = -qJD(1) * t586 + t406 * t468;
t361 = -qJD(4) * pkin(4) - t615;
t345 = pkin(5) * t417 + t361;
t607 = cos(pkin(11));
t515 = t607 * t465;
t459 = sin(pkin(11));
t589 = t459 * t469;
t400 = t461 * t589 + t515;
t591 = t459 * t460;
t367 = t400 * t464 + t468 * t591;
t514 = t607 * t469;
t590 = t459 * t465;
t398 = -t461 * t514 + t590;
t516 = t460 * t607;
t369 = -t398 * t464 + t468 * t516;
t399 = t461 * t515 + t589;
t401 = -t461 * t590 + t514;
t587 = t460 * t469;
t405 = t461 * t468 - t464 * t587;
t458 = qJ(5) + qJ(6);
t454 = sin(t458);
t455 = cos(t458);
t588 = t460 * t465;
t620 = t345 * t354 - g(1) * (-t367 * t455 - t401 * t454) - g(2) * (t369 * t455 - t399 * t454) - g(3) * (-t405 * t455 - t454 * t588) + t322;
t552 = qJD(5) * t468;
t485 = -t463 * t552 - t464 * t548;
t349 = qJD(2) * t485 + qJD(5) * t548 + t463 * qJDD(4) + t467 * t542;
t521 = t468 * t547;
t543 = qJDD(2) * t464;
t489 = t521 + t543;
t413 = qJDD(5) + t489;
t565 = qJD(2) * t460;
t530 = t465 * t565;
t436 = qJD(1) * t530;
t546 = qJDD(1) * t460;
t519 = t469 * t546;
t495 = qJDD(3) + t436 - t519;
t378 = qJDD(2) * t470 + t495;
t545 = qJDD(1) * t461;
t518 = t468 * t545;
t329 = qJDD(4) * pkin(9) + qJD(4) * t615 + t378 * t464 + t518;
t506 = pkin(4) * t468 + pkin(9) * t464;
t414 = qJD(4) * t506 + qJD(3);
t520 = t465 * t546;
t344 = t520 + t426 * qJDD(2) + (t414 + t532) * qJD(2);
t343 = t467 * t344;
t475 = -qJD(5) * t333 - t463 * t329 + t343;
t314 = pkin(5) * t413 - pkin(10) * t349 + t475;
t350 = qJD(5) * t419 - t467 * qJDD(4) + t463 * t621;
t553 = qJD(5) * t467;
t538 = -t467 * t329 - t463 * t344 - t383 * t553;
t555 = qJD(5) * t463;
t488 = t362 * t555 + t538;
t315 = -pkin(10) * t350 - t488;
t513 = t466 * t314 - t462 * t315;
t619 = t345 * t497 - g(1) * (-t367 * t454 + t401 * t455) - g(2) * (t369 * t454 + t399 * t455) - g(3) * (-t405 * t454 + t455 * t588) + t513;
t410 = qJDD(6) + t413;
t618 = t410 * MDP(26) + (-t354 ^ 2 + t497 ^ 2) * MDP(23) - t354 * MDP(22) * t497;
t421 = t462 * t467 + t463 * t466;
t390 = t421 * t468;
t583 = t464 * t465;
t617 = -(-t463 * t583 + t467 * t469) * t569 + t467 * t414;
t525 = t468 * t548;
t581 = t465 * t467;
t616 = (t463 * t469 + t464 * t581) * t569 - t463 * t414 - t426 * t553 - t470 * t525;
t482 = g(1) * t400 + g(2) * t398 - g(3) * t587;
t614 = -qJD(6) * t467 - t553;
t613 = qJD(5) + qJD(6);
t504 = g(1) * t401 + g(2) * t399;
t579 = qJDD(1) - g(3);
t612 = -t579 * t588 + t504;
t567 = qJD(2) * qJ(3);
t425 = t533 + t567;
t611 = qJD(4) * (-t425 + t533 - t567) - qJDD(4) * t470;
t481 = g(3) * t588 + t504;
t610 = qJD(5) * (t446 * t470 + t362) + t481;
t512 = t349 * t462 + t466 * t350;
t319 = -qJD(6) * t497 + t512;
t609 = pkin(9) + pkin(10);
t608 = pkin(10) * t468;
t606 = qJDD(2) * pkin(2);
t332 = -t362 * t463 + t467 * t383;
t323 = -pkin(10) * t419 + t332;
t321 = pkin(5) * t446 + t323;
t605 = t321 * t466;
t604 = t324 * t466;
t603 = t349 * t463;
t420 = t462 * t463 - t466 * t467;
t601 = t410 * t420;
t600 = t410 * t421;
t599 = t413 * t463;
t598 = t413 * t467;
t597 = t417 * t446;
t596 = t419 * t446;
t594 = t419 * t467;
t593 = t454 * t464;
t592 = t455 * t464;
t585 = t462 * t314;
t584 = t463 * t470;
t582 = t464 * t467;
t580 = t467 * t468;
t439 = t470 * t582;
t517 = pkin(5) - t584;
t539 = pkin(10) * t582;
t578 = (-t439 + (-t426 + t608) * t463) * qJD(5) + (t468 * t517 + t539) * qJD(4) + t617;
t559 = qJD(4) * t464;
t527 = t463 * t559;
t484 = -t467 * t552 + t527;
t554 = qJD(5) * t464;
t577 = -pkin(10) * t484 + t554 * t584 + t616;
t422 = t506 * qJD(2);
t576 = t463 * t422 + t467 * t615;
t490 = t420 * t464;
t575 = -qJD(2) * t490 - t420 * t613;
t487 = t421 * qJD(2);
t574 = t421 * t613 + t464 * t487;
t573 = t463 * t426 + t439;
t457 = t468 ^ 2;
t572 = t464 ^ 2 - t457;
t471 = qJD(4) ^ 2;
t472 = qJD(2) ^ 2;
t571 = -t471 - t472;
t566 = qJD(2) * t425;
t562 = qJD(4) * t417;
t561 = qJD(4) * t419;
t558 = qJD(4) * t468;
t557 = qJD(4) * t470;
t556 = qJD(5) * t446;
t550 = qJD(6) * t466;
t544 = qJDD(2) * qJ(3);
t541 = qJDD(4) * t464;
t537 = t466 * t349 - t462 * t350 - t417 * t550;
t536 = -qJD(4) * t443 - t406 * t559 - t464 * t545;
t534 = qJD(5) * t609;
t531 = t465 * t568;
t529 = t469 * t565;
t528 = t463 * t564;
t523 = t463 * t551;
t510 = -t378 + t566;
t509 = qJD(6) * t321 + t315;
t508 = qJD(2) + t554;
t507 = -t376 + (t528 + t555) * pkin(5);
t408 = t467 * t422;
t435 = t609 * t467;
t503 = qJD(6) * t435 - t615 * t463 + t408 + (pkin(5) * t468 + t539) * qJD(2) + t467 * t534;
t434 = t609 * t463;
t502 = pkin(10) * t528 + qJD(6) * t434 + t463 * t534 + t576;
t317 = t321 * t462 + t604;
t412 = t467 * t426;
t352 = -pkin(10) * t580 + t464 * t517 + t412;
t365 = -t463 * t608 + t573;
t500 = t352 * t462 + t365 * t466;
t372 = -t405 * t463 + t460 * t581;
t373 = t405 * t467 + t463 * t588;
t499 = t372 * t466 - t373 * t462;
t498 = t372 * t462 + t373 * t466;
t404 = t468 * t587 + t586;
t493 = t446 * t553 + t599;
t492 = -t446 * t555 + t598;
t491 = t440 * t420;
t318 = -t419 * t551 + t537;
t486 = g(1) * (t400 * t468 - t464 * t591) + g(2) * (t398 * t468 + t464 * t516) - g(3) * t404;
t330 = -qJDD(4) * pkin(4) - t378 * t468 - t536;
t479 = -pkin(9) * t413 + t361 * t446;
t478 = t519 + t482;
t476 = qJDD(3) - t478;
t474 = pkin(9) * t556 + t330 + t486;
t379 = t520 + t544 + (qJD(3) + t532) * qJD(2);
t473 = qJD(2) * t501 - t470 * t471 + t379 - t481 + t544;
t452 = qJDD(4) * t468;
t449 = -pkin(5) * t467 - pkin(4);
t416 = -qJD(2) * pkin(2) + t501;
t415 = (pkin(5) * t463 - t470) * t468;
t391 = t420 * t468;
t384 = t495 - t606;
t380 = -pkin(5) * t484 + t464 * t557;
t371 = qJD(4) * t405 - t468 * t530;
t370 = -qJD(4) * t404 + t464 * t530;
t335 = -t468 * t523 + (t580 * t613 - t527) * t466 + t485 * t462;
t334 = qJD(4) * t490 - t390 * t613;
t326 = qJD(5) * t372 + t370 * t467 + t463 * t529;
t325 = -qJD(5) * t373 - t370 * t463 + t467 * t529;
t320 = pkin(5) * t350 + t330;
t316 = -t324 * t462 + t605;
t1 = [t579 * MDP(1) + (qJDD(1) * t461 ^ 2 - g(3)) * MDP(7) + (-qJD(4) * t371 - qJDD(4) * t404) * MDP(13) + (-qJD(4) * t370 - qJDD(4) * t405) * MDP(14) + (t325 * t446 + t350 * t404 + t371 * t417 + t372 * t413) * MDP(20) + (-t326 * t446 + t349 * t404 + t371 * t419 - t373 * t413) * MDP(21) + ((-qJD(6) * t498 + t325 * t466 - t326 * t462) * t440 + t499 * t410 + t371 * t354 + t404 * t319) * MDP(27) + (-(qJD(6) * t499 + t325 * t462 + t326 * t466) * t440 - t498 * t410 - t371 * t497 + t404 * t318) * MDP(28) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t465 + t469 * t472) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t469 + t465 * t472) + ((-t384 + t566) * MDP(7) + (MDP(13) * t464 + MDP(14) * t468) * t472) * t469 + ((qJD(2) * t416 + t379) * MDP(7) + t489 * MDP(13) + t621 * MDP(14)) * t465) * t460; (qJDD(2) * t457 - 0.2e1 * t464 * t521) * MDP(8) + (t379 * qJ(3) + t425 * qJD(3) - t384 * pkin(2) - g(1) * (-pkin(2) * t400 + qJ(3) * t401) - g(2) * (-pkin(2) * t398 + qJ(3) * t399) + (-g(3) * (pkin(2) * t469 + qJ(3) * t465) + (-t416 * t465 - t425 * t469) * qJD(1)) * t460) * MDP(7) + (t476 - 0.2e1 * t606) * MDP(5) + ((t417 * t467 + t419 * t463) * t559 + (-t603 - t350 * t467 + (t417 * t463 - t594) * qJD(5)) * t468) * MDP(16) + ((t352 * t466 - t365 * t462) * t410 + t513 * t464 + t316 * t558 + t380 * t354 + t415 * t319 + t320 * t390 + t345 * t335 - g(1) * (-t400 * t454 + t401 * t592) - g(2) * (-t398 * t454 + t399 * t592) + (t462 * t577 + t466 * t578) * t440 + (-t317 * t464 - t440 * t500) * qJD(6) + (t354 * t531 - g(3) * (t454 * t469 + t455 * t583)) * t460) * MDP(27) + (t464 * t611 + t473 * t468) * MDP(14) + (t473 * t464 - t468 * t611) * MDP(13) + (t349 * t580 + t419 * t485) * MDP(15) + 0.2e1 * (-t464 * t542 + t547 * t572) * MDP(9) + ((t446 * t560 - t350) * t464 + (-t493 - t562) * t468) * MDP(18) + (-t319 * t464 - t335 * t440 - t354 * t558 - t390 * t410) * MDP(25) + (t413 * t464 + t446 * t558) * MDP(19) + (t410 * t464 + t440 * t558) * MDP(26) + ((-t446 * t548 + t349) * t464 + (t492 + t561) * t468) * MDP(17) + t478 * MDP(3) + (-t500 * t410 - (t509 * t466 - t322 + t585) * t464 - t317 * t558 - t380 * t497 + t415 * t318 - t320 * t391 + t345 * t334 - g(1) * (-t400 * t455 - t401 * t593) - g(2) * (-t398 * t455 - t399 * t593) + ((-qJD(6) * t352 + t577) * t466 + (qJD(6) * t365 - t578) * t462) * t440 + (-t497 * t531 - g(3) * (-t454 * t583 + t455 * t469)) * t460) * MDP(28) + (t318 * t464 + t334 * t440 - t391 * t410 - t497 * t558) * MDP(24) + (-t318 * t390 + t319 * t391 - t334 * t354 + t335 * t497) * MDP(23) + (-t318 * t391 - t334 * t497) * MDP(22) + (-t468 * t471 - t541) * MDP(11) + (-t464 * t471 + t452) * MDP(10) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t544 - t612) * MDP(6) + t612 * MDP(4) + qJDD(2) * MDP(2) + (-t573 * t413 + t616 * t446 + t482 * t467 + ((-t361 * t467 + t419 * t470) * qJD(4) + t610 * t463 + t538) * t464 + (-qJD(4) * t333 + t330 * t467 - t349 * t470 - t361 * t555 + t419 * t533) * t468) * MDP(21) + (t412 * t413 + t617 * t446 + (-t426 * t556 + t482) * t463 + (t417 * t557 + t343 + (-qJD(4) * t361 - qJD(5) * t383 - t413 * t470 - t329) * t463 - t610 * t467) * t464 + (t417 * t533 + t361 * t553 + t330 * t463 - t470 * t350 + (-t446 * t584 + t332) * qJD(4)) * t468) * MDP(20); qJDD(2) * MDP(5) - t472 * MDP(6) + (t436 + t476 - t566 - t606) * MDP(7) + (t464 * t571 + t452) * MDP(13) + (t468 * t571 - t541) * MDP(14) + (-t350 * t468 + (t562 - t599) * t464 + (-t463 * t558 - t467 * t508) * t446) * MDP(20) + (-t349 * t468 + (t561 - t598) * t464 + (t463 * t508 - t525) * t446) * MDP(21) + (qJD(2) * t491 + (-qJD(4) * t421 * t440 - t319) * t468 + ((t462 * t555 + t466 * t614 + t523) * t440 - t600 + qJD(4) * t354) * t464) * MDP(27) + (t440 * t487 + (qJD(4) * t491 - t318) * t468 + (-(t462 * t614 - t463 * t550 - t466 * t555) * t440 + t601 - qJD(4) * t497) * t464) * MDP(28); MDP(10) * t542 - MDP(11) * t543 + qJDD(4) * MDP(12) + (qJD(4) * t376 - t468 * t510 - t486 + t536) * MDP(13) + (g(1) * t367 - g(2) * t369 + g(3) * t405 + t510 * t464 - t518) * MDP(14) + (t446 * t594 + t603) * MDP(15) + ((t349 - t597) * t467 + (-t350 - t596) * t463) * MDP(16) + ((-t419 * t468 + t446 * t582) * qJD(2) + t493) * MDP(17) + ((-t446 * t463 * t464 + t417 * t468) * qJD(2) + t492) * MDP(18) + (-pkin(4) * t350 - t376 * t417 - t408 * t446 + (t446 * t615 + t479) * t463 - t474 * t467) * MDP(20) + (-pkin(4) * t349 - t376 * t419 + t446 * t576 + t463 * t474 + t467 * t479) * MDP(21) + (t318 * t421 - t497 * t575) * MDP(22) + (-t318 * t420 - t319 * t421 - t354 * t575 + t497 * t574) * MDP(23) + (t440 * t575 + t600) * MDP(24) + (-t440 * t574 - t601) * MDP(25) + ((-t434 * t466 - t435 * t462) * t410 + t449 * t319 + t320 * t420 + (t462 * t502 - t466 * t503) * t440 + t507 * t354 + t574 * t345 - t486 * t455) * MDP(27) + (-(-t434 * t462 + t435 * t466) * t410 + t449 * t318 + t320 * t421 + (t462 * t503 + t466 * t502) * t440 - t507 * t497 + t575 * t345 + t486 * t454) * MDP(28) + (-t446 * MDP(19) - MDP(20) * t332 + t333 * MDP(21) + MDP(24) * t497 + t354 * MDP(25) - t440 * MDP(26) - t316 * MDP(27) + t317 * MDP(28)) * t563 + (MDP(8) * t464 * t468 - MDP(9) * t572) * t472; t419 * t417 * MDP(15) + (-t417 ^ 2 + t419 ^ 2) * MDP(16) + (t349 + t597) * MDP(17) + (-t350 + t596) * MDP(18) + t413 * MDP(19) + (t333 * t446 - t361 * t419 - g(1) * (-t367 * t463 + t401 * t467) - g(2) * (t369 * t463 + t399 * t467) - g(3) * t372 + t475) * MDP(20) + (t332 * t446 + t361 * t417 - g(1) * (-t367 * t467 - t401 * t463) - g(2) * (t369 * t467 - t399 * t463) + g(3) * t373 + t488) * MDP(21) + (t318 + t623) * MDP(24) + (-t319 - t622) * MDP(25) + (-(-t323 * t462 - t604) * t440 - t317 * qJD(6) + (-t354 * t419 + t410 * t466 - t440 * t551) * pkin(5) + t619) * MDP(27) + ((-t324 * t440 - t314) * t462 + (t323 * t440 - t509) * t466 + (-t410 * t462 + t419 * t497 - t440 * t550) * pkin(5) + t620) * MDP(28) + t618; (t537 + t623) * MDP(24) + (-t512 - t622) * MDP(25) + (t317 * t440 + t619) * MDP(27) + (-t466 * t315 + t316 * t440 - t585 + t620) * MDP(28) + (-MDP(24) * t595 + MDP(25) * t497 - MDP(27) * t317 - MDP(28) * t605) * qJD(6) + t618;];
tau  = t1;
