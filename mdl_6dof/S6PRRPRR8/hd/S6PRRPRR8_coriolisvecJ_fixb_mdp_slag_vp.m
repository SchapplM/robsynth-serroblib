% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:51
% EndTime: 2019-03-08 22:41:02
% DurationCPUTime: 6.75s
% Computational Cost: add. (3766->502), mult. (10189->716), div. (0->0), fcn. (8056->12), ass. (0->223)
t470 = sin(pkin(7));
t476 = sin(qJ(3));
t567 = qJD(3) * t476;
t543 = t470 * t567;
t455 = pkin(3) * t543;
t480 = cos(qJ(3));
t515 = pkin(10) * t476 - qJ(4) * t480;
t565 = qJD(4) * t476;
t484 = (qJD(3) * t515 - t565) * t470;
t477 = sin(qJ(2));
t471 = sin(pkin(6));
t573 = qJD(1) * t471;
t547 = t477 * t573;
t527 = t470 * t547;
t623 = -t455 - t484 + t527;
t472 = cos(pkin(7));
t586 = t477 * t480;
t481 = cos(qJ(2));
t587 = t476 * t481;
t499 = t472 * t586 + t587;
t393 = t499 * t573;
t592 = t472 * t476;
t464 = pkin(2) * t592;
t593 = t470 * t480;
t609 = pkin(4) + pkin(9);
t622 = (t593 * t609 + t464) * qJD(3) - t393;
t589 = t476 * t477;
t551 = t471 * t589;
t528 = t472 * t551;
t566 = qJD(3) * t480;
t541 = t472 * t566;
t546 = t481 * t573;
t621 = -pkin(2) * t541 - qJD(1) * t528 + t480 * t546;
t620 = -t472 * qJD(4) + t621;
t570 = qJD(2) * t471;
t539 = qJD(1) * t570;
t473 = cos(pkin(6));
t572 = qJD(1) * t473;
t548 = t470 * t572;
t619 = qJD(3) * t548 + t481 * t539;
t571 = qJD(2) * t470;
t422 = pkin(9) * t571 + t547;
t608 = qJD(2) * pkin(2);
t438 = t546 + t608;
t361 = t476 * t422 + (-t438 * t472 - t548) * t480;
t558 = qJD(4) + t361;
t569 = qJD(2) * t472;
t459 = qJD(3) + t569;
t475 = sin(qJ(5));
t479 = cos(qJ(5));
t544 = t480 * t571;
t400 = t459 * t475 + t479 * t544;
t399 = qJD(6) + t400;
t468 = t476 ^ 2;
t618 = MDP(6) * (-t480 ^ 2 + t468);
t568 = qJD(2) * t476;
t545 = t470 * t568;
t449 = qJD(5) + t545;
t478 = cos(qJ(6));
t525 = t475 * t544;
t402 = t459 * t479 - t525;
t474 = sin(qJ(6));
t599 = t402 * t474;
t368 = -t478 * t449 + t599;
t617 = t368 * t449;
t616 = t476 * t480;
t594 = t470 * t476;
t461 = pkin(9) * t594;
t549 = -pkin(2) * t480 - pkin(3);
t378 = pkin(4) * t594 + t461 + (-pkin(10) + t549) * t472;
t482 = -pkin(3) - pkin(10);
t607 = qJ(4) * t476;
t501 = t480 * t482 - t607;
t389 = (-pkin(2) + t501) * t470;
t562 = qJD(5) * t479;
t564 = qJD(5) * t475;
t615 = -t378 * t562 + t389 * t564 - t622 * t475 + t479 * t623;
t581 = -t543 * t609 - t620;
t614 = t475 * t378 + t479 * t389;
t613 = pkin(9) * t593 + t464;
t559 = pkin(4) * t545 + t558;
t557 = qJD(2) * qJD(3);
t538 = t470 * t557;
t451 = t480 * t538;
t342 = t459 * t482 + t559;
t458 = t472 * t572;
t358 = t458 + (qJD(2) * t501 - t438) * t470;
t325 = t342 * t475 + t358 * t479;
t412 = t438 * t592;
t524 = t477 * t539;
t503 = t480 * t524;
t339 = qJD(3) * t412 + t422 * t566 + t472 * t503 + t476 * t619;
t335 = pkin(4) * t451 + t339;
t521 = t476 * t538;
t576 = pkin(3) * t521 + t470 * t524;
t357 = qJD(2) * t484 + t576;
t485 = -qJD(5) * t325 + t479 * t335 - t475 * t357;
t313 = -pkin(5) * t451 - t485;
t612 = t399 * (pkin(5) * t402 + pkin(11) * t399) + t313;
t578 = -qJD(3) * t613 + t393;
t611 = -t339 * t472 + t578 * t459;
t610 = -qJD(5) * t614 + t475 * t623 + t622 * t479;
t372 = -qJD(5) * t400 + t475 * t521;
t560 = qJD(6) * t478;
t550 = t478 * t372 + t449 * t560 + t474 * t451;
t561 = qJD(6) * t474;
t331 = -t402 * t561 + t550;
t606 = t331 * t474;
t604 = t368 * t399;
t370 = t402 * t478 + t449 * t474;
t603 = t370 * t399;
t602 = t399 * t482;
t601 = t400 * t449;
t600 = t402 * t449;
t598 = t449 * t475;
t597 = t449 * t479;
t596 = t449 * t482;
t467 = t470 ^ 2;
t483 = qJD(2) ^ 2;
t595 = t467 * t483;
t591 = t472 * t481;
t373 = -qJD(5) * t525 + t459 * t562 - t479 * t521;
t590 = t474 * t373;
t588 = t476 * t478;
t585 = t478 * t373;
t584 = t479 * t331;
t542 = t470 * t566;
t583 = -pkin(5) * t542 - t610;
t520 = pkin(5) * t479 + pkin(11) * t475;
t582 = (-pkin(4) - t520) * t545 - qJD(5) * t520 - t558;
t362 = t480 * t422 + t476 * t548 + t412;
t354 = pkin(4) * t544 + t362;
t457 = pkin(3) * t545;
t390 = t515 * t571 + t457;
t579 = t475 * t354 + t479 * t390;
t577 = pkin(9) * t543 + t620;
t563 = qJD(5) * t478;
t556 = MDP(10) - MDP(13);
t555 = MDP(11) - MDP(14);
t554 = t467 * t608;
t553 = t476 * t597;
t552 = t475 * t593;
t403 = -t472 * qJ(4) - t613;
t540 = t449 * t562;
t491 = t475 * t335 + t342 * t562 + t479 * t357 - t358 * t564;
t312 = pkin(11) * t451 + t491;
t504 = t476 * t524;
t514 = -t422 * t567 + t438 * t541 - t472 * t504 + t480 * t619;
t336 = -t459 * qJD(4) - t514;
t329 = -pkin(4) * t521 - t336;
t317 = pkin(5) * t373 - pkin(11) * t372 + t329;
t536 = -t312 * t474 + t478 * t317;
t534 = t372 * t474 - t478 * t451;
t533 = t399 * t478;
t441 = pkin(5) * t475 - pkin(11) * t479 + qJ(4);
t531 = pkin(11) * t544 - qJD(6) * t441 + t579;
t529 = t595 * t616;
t388 = pkin(4) * t593 - t403;
t526 = t470 * t477 * t570;
t396 = (t474 * t480 + t475 * t588) * t571;
t519 = -t475 * t563 - t396;
t417 = t472 * t475 + t479 * t593;
t418 = t472 * t479 - t552;
t346 = pkin(5) * t417 - pkin(11) * t418 + t388;
t518 = -pkin(11) * t542 - qJD(6) * t346 + t615;
t344 = pkin(11) * t594 + t614;
t384 = -qJD(5) * t417 + t475 * t543;
t385 = -qJD(5) * t552 + t472 * t562 - t479 * t543;
t517 = -pkin(5) * t385 + pkin(11) * t384 + qJD(6) * t344 - t581;
t516 = -pkin(3) * t480 - t607;
t452 = t459 * qJ(4);
t345 = t452 + t354;
t513 = t312 * t478 + t317 * t474;
t321 = pkin(11) * t449 + t325;
t326 = pkin(5) * t400 - pkin(11) * t402 + t345;
t315 = t321 * t478 + t326 * t474;
t512 = t321 * t474 - t326 * t478;
t324 = t342 * t479 - t358 * t475;
t511 = t354 * t479 - t390 * t475;
t379 = -t473 * t593 + (-t480 * t591 + t589) * t471;
t416 = -t470 * t471 * t481 + t472 * t473;
t356 = t379 * t475 + t416 * t479;
t498 = t472 * t587 + t586;
t380 = t471 * t498 + t473 * t594;
t510 = t356 * t478 + t380 * t474;
t509 = -t356 * t474 + t380 * t478;
t507 = t378 * t479 - t389 * t475;
t505 = t379 * t479 - t416 * t475;
t502 = -t459 * t544 + t451;
t500 = -t418 * t474 + t470 * t588;
t387 = t418 * t478 + t474 * t594;
t497 = t345 * t476 + t482 * t566;
t496 = -t399 * t560 - t590;
t495 = t399 * t561 - t585;
t492 = t449 * t370;
t488 = t362 * t459 - t339;
t487 = (-qJ(4) * t566 - t565) * t470;
t320 = -pkin(5) * t449 - t324;
t486 = -pkin(11) * t373 + (t320 + t324) * t399;
t435 = t479 * t451;
t410 = -qJ(4) * t544 + t457;
t405 = t472 * t549 + t461;
t404 = (-pkin(2) + t516) * t470;
t397 = -t438 * t470 + t458;
t395 = t474 * t475 * t545 - t478 * t544;
t392 = t455 + t487;
t366 = t458 + (qJD(2) * t516 - t438) * t470;
t363 = qJD(2) * t487 + t576;
t359 = t366 * t545;
t351 = -t452 - t362;
t350 = -pkin(3) * t459 + t558;
t348 = -qJD(2) * t528 - qJD(3) * t551 + (t481 * t570 + (t470 * t473 + t471 * t591) * qJD(3)) * t480;
t347 = t473 * t543 + (qJD(2) * t499 + qJD(3) * t498) * t471;
t343 = -pkin(5) * t594 - t507;
t341 = qJD(6) * t387 + t384 * t474 - t478 * t542;
t340 = qJD(6) * t500 + t384 * t478 + t474 * t542;
t332 = qJD(6) * t370 + t534;
t327 = -pkin(5) * t544 - t511;
t323 = qJD(5) * t505 + t347 * t475 + t479 * t526;
t322 = qJD(5) * t356 - t347 * t479 + t475 * t526;
t311 = -qJD(6) * t315 + t536;
t310 = -qJD(6) * t512 + t513;
t1 = [(-t336 * t380 + t339 * t379 + t347 * t350 - t348 * t351 + t363 * t416) * MDP(15) + (-t322 * t449 + t348 * t400 + t373 * t380) * MDP(21) + (-t323 * t449 + t348 * t402 + t372 * t380) * MDP(22) + ((-qJD(6) * t510 - t323 * t474 + t348 * t478) * t399 + t509 * t373 + t322 * t368 - t505 * t332) * MDP(28) + (-(qJD(6) * t509 + t323 * t478 + t348 * t474) * t399 - t510 * t373 + t322 * t370 - t505 * t331) * MDP(29) + (-t347 * t556 - t348 * t555) * t459 + (-MDP(4) * t481 + (-MDP(3) + (t476 * t555 - t480 * t556) * t467) * t477) * t483 * t471 + ((t347 * t476 + t348 * t480) * MDP(12) + t366 * t471 * t477 * MDP(15) + (-t380 * t476 * MDP(12) + (MDP(12) * t379 + MDP(21) * t505 - MDP(22) * t356) * t480 + (t476 * t556 + t480 * t555) * t416) * qJD(3)) * t571; ((t397 * t470 - t554) * t567 + t611) * MDP(10) + (-t514 * t472 + t621 * t459 + (-t480 * t554 + (pkin(9) * t459 * t476 + t397 * t480) * t470) * qJD(3)) * MDP(11) + (-t336 * t480 + t339 * t476 + (t350 * t480 + t351 * t476) * qJD(3) + ((qJD(3) * t405 - t577) * t480 + (qJD(3) * t403 - t578) * t476) * qJD(2)) * t470 * MDP(12) + (-t467 * t503 + (-t366 * t567 + t363 * t480 + (t392 * t480 - t404 * t567) * qJD(2)) * t470 - t611) * MDP(13) + (t467 * t504 - t336 * t472 - t577 * t459 + (-t366 * t566 - t363 * t476 + (-t392 * t476 - t404 * t566) * qJD(2)) * t470) * MDP(14) + (t336 * t403 + t339 * t405 + t363 * t404 + (t392 - t527) * t366 + t577 * t351 - t578 * t350) * MDP(15) + (t372 * t418 + t384 * t402) * MDP(16) + (-t372 * t417 - t373 * t418 - t384 * t400 - t385 * t402) * MDP(17) + (t384 * t449 + (t372 * t476 + (qJD(2) * t418 + t402) * t566) * t470) * MDP(18) + (-t385 * t449 + (-t373 * t476 + (-qJD(2) * t417 - t400) * t566) * t470) * MDP(19) + (t449 * t470 + t467 * t568) * MDP(20) * t566 + (t329 * t417 + t345 * t385 + t388 * t373 + t610 * t449 + t581 * t400 + (t485 * t476 + (qJD(2) * t507 + t324) * t566) * t470) * MDP(21) + (t329 * t418 + t345 * t384 + t388 * t372 + t615 * t449 + t581 * t402 + (-t491 * t476 + (-qJD(2) * t614 - t325) * t566) * t470) * MDP(22) + (t331 * t387 + t340 * t370) * MDP(23) + (t331 * t500 - t332 * t387 - t340 * t368 - t341 * t370) * MDP(24) + (t331 * t417 + t340 * t399 + t370 * t385 + t373 * t387) * MDP(25) + (-t332 * t417 - t341 * t399 - t368 * t385 + t373 * t500) * MDP(26) + (t373 * t417 + t385 * t399) * MDP(27) + ((-t344 * t474 + t346 * t478) * t373 + t311 * t417 - t512 * t385 + t343 * t332 - t313 * t500 + t320 * t341 + (t474 * t518 - t478 * t517) * t399 + t583 * t368) * MDP(28) + (-(t344 * t478 + t346 * t474) * t373 - t310 * t417 - t315 * t385 + t343 * t331 + t313 * t387 + t320 * t340 + (t474 * t517 + t478 * t518) * t399 + t583 * t370) * MDP(29) + (MDP(7) * t542 - MDP(8) * t543) * (t459 + t569) + 0.2e1 * (MDP(5) * t616 - t618) * t467 * t557; t595 * t618 + t502 * MDP(7) + (-qJD(3) + t459) * MDP(8) * t545 + (-t397 * t545 + t488) * MDP(10) + (-t361 * t459 - t397 * t544 - t514) * MDP(11) + ((-qJ(4) * qJD(3) - t351 - t362) * t476 + (-pkin(3) * qJD(3) - t350 + t558) * t480) * MDP(12) * t571 + (-t410 * t544 + t359 - t488) * MDP(13) + (t558 * t459 + (t366 * t480 + t410 * t476) * t571 - t336) * MDP(14) + (-pkin(3) * t339 - qJ(4) * t336 - t350 * t362 - t351 * t558 - t366 * t410) * MDP(15) + (t372 * t479 - t402 * t598) * MDP(16) + ((-t373 - t600) * t479 + (-t372 + t601) * t475) * MDP(17) + (-t449 * t564 + t435 + (-t402 * t480 - t476 * t598) * t571) * MDP(18) + (-t540 + (-t553 + (-qJD(3) * t475 + t400) * t480) * t571) * MDP(19) + (qJ(4) * t373 + t329 * t475 - t511 * t449 + t559 * t400 + (t345 * t479 - t475 * t596) * qJD(5) + (-t324 * t480 + t479 * t497) * t571) * MDP(21) + (qJ(4) * t372 + t329 * t479 + t579 * t449 + t559 * t402 + (-t345 * t475 - t479 * t596) * qJD(5) + (t325 * t480 - t475 * t497) * t571) * MDP(22) + (t478 * t584 + (-t479 * t561 + t519) * t370) * MDP(23) + (t368 * t396 + t370 * t395 + (t368 * t478 + t370 * t474) * t564 + (-t606 - t332 * t478 + (t368 * t474 - t370 * t478) * qJD(6)) * t479) * MDP(24) + (t331 * t475 + t519 * t399 + (t492 - t495) * t479) * MDP(25) + (-t332 * t475 + (t474 * t564 + t395) * t399 + (t496 - t617) * t479) * MDP(26) + (t373 * t475 + t399 * t597) * MDP(27) + (t441 * t585 - t320 * t395 - t327 * t368 + (t474 * t531 - t478 * t582) * t399 + (-t320 * t474 * qJD(5) + t311 + (qJD(5) * t368 + t496) * t482) * t475 + (-t512 * t545 + t320 * t560 + t313 * t474 - t482 * t332 + (-t474 * t602 - t512) * qJD(5)) * t479) * MDP(28) + (-t441 * t590 - t320 * t396 - t327 * t370 + (t474 * t582 + t478 * t531) * t399 + (-t320 * t563 - t310 + (qJD(5) * t370 + t495) * t482) * t475 + (-t315 * t545 - t320 * t561 + t313 * t478 - t482 * t331 + (-t478 * t602 - t315) * qJD(5)) * t479) * MDP(29) - MDP(5) * t529 - t449 * MDP(20) * t544; t502 * MDP(12) + MDP(13) * t529 + (-t459 ^ 2 - t468 * t595) * MDP(14) + (t351 * t459 + t339 + t359) * MDP(15) + (-t400 * t459 - t449 * t598 + t435) * MDP(21) + (-t540 - t402 * t459 + (-t475 * t566 - t553) * t571) * MDP(22) + (-t479 * t332 + (-t478 * t459 - t474 * t597) * t399 + (t496 + t617) * t475) * MDP(28) + (-t584 + (t474 * t459 - t478 * t597) * t399 + (t492 + t495) * t475) * MDP(29); -t400 ^ 2 * MDP(17) + (t372 + t601) * MDP(18) + (-t373 + t600) * MDP(19) + MDP(20) * t451 + (t325 * t449 + t485) * MDP(21) + (t324 * t449 + t345 * t400 - t491) * MDP(22) + (t370 * t533 + t606) * MDP(23) + ((t331 - t604) * t478 + (-t332 - t603) * t474) * MDP(24) + (t399 * t533 + t590) * MDP(25) + (-t399 ^ 2 * t474 + t585) * MDP(26) + (-pkin(5) * t332 - t325 * t368 + t486 * t474 - t478 * t612) * MDP(28) + (-pkin(5) * t331 - t325 * t370 + t474 * t612 + t486 * t478) * MDP(29) + (MDP(16) * t400 + t402 * MDP(17) - t345 * MDP(21) - t370 * MDP(25) + t368 * MDP(26) - t399 * MDP(27) + MDP(28) * t512 + MDP(29) * t315) * t402; t370 * t368 * MDP(23) + (-t368 ^ 2 + t370 ^ 2) * MDP(24) + (t550 + t604) * MDP(25) + (-t534 + t603) * MDP(26) + t373 * MDP(27) + (t315 * t399 - t320 * t370 + t536) * MDP(28) + (t320 * t368 - t399 * t512 - t513) * MDP(29) + (-MDP(25) * t599 - MDP(26) * t370 - MDP(28) * t315 + MDP(29) * t512) * qJD(6);];
tauc  = t1;
