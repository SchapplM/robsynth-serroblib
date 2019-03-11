% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:35
% EndTime: 2019-03-09 04:25:49
% DurationCPUTime: 7.93s
% Computational Cost: add. (8741->532), mult. (29010->727), div. (0->0), fcn. (24963->12), ass. (0->225)
t516 = sin(pkin(12));
t644 = cos(pkin(7));
t646 = cos(qJ(3));
t576 = t644 * t646;
t518 = sin(pkin(6));
t519 = cos(pkin(12));
t630 = t518 * t519;
t546 = t576 * t630;
t540 = qJD(1) * t546;
t517 = sin(pkin(7));
t645 = cos(pkin(6));
t596 = t645 * t517;
t560 = t646 * t596;
t522 = sin(qJ(3));
t620 = qJD(1) * t518;
t601 = t522 * t620;
t453 = -qJD(1) * t560 + t516 * t601 - t540;
t572 = t645 * t644;
t602 = t519 * t620;
t479 = -qJD(1) * t572 + t517 * t602 - qJD(3);
t521 = sin(qJ(5));
t524 = cos(qJ(5));
t414 = -t524 * t453 - t479 * t521;
t597 = t522 * t644;
t466 = t518 * (t516 * t646 + t519 * t597) + t522 * t596;
t654 = qJD(1) * t466;
t664 = qJD(5) + t654;
t667 = t414 * t664;
t632 = t516 * t518;
t666 = t522 * t632 - t560;
t532 = t518 * (t516 * t576 + t519 * t522);
t474 = qJD(1) * t532;
t618 = qJD(3) * t522;
t600 = t517 * t618;
t665 = t474 - t600;
t484 = t517 * t630 - t572;
t598 = t518 * t644;
t531 = (t519 * t598 + t596) * pkin(9);
t605 = pkin(1) * t645;
t622 = qJ(2) * t630 + t516 * t605;
t463 = t531 + t622;
t511 = t519 * t605;
t529 = t645 * pkin(2) + (-pkin(9) * t644 - qJ(2)) * t632;
t467 = t511 + t529;
t633 = t516 * t517;
t476 = (-pkin(2) * t519 - pkin(9) * t633 - pkin(1)) * t518;
t584 = t646 * t630;
t501 = qJD(2) * t584;
t558 = qJD(3) * t576;
t579 = t516 * t598;
t559 = qJD(2) * t579;
t599 = qJD(3) * t646;
t580 = t517 * t599;
t527 = -t522 * (qJD(3) * t463 + t559) + t467 * t558 + t476 * t580 + t501;
t383 = t484 * qJD(4) - t527;
t663 = t666 * qJD(3);
t662 = MDP(4) * t516 + MDP(5) * t519;
t413 = qJD(6) + t414;
t523 = cos(qJ(6));
t416 = t453 * t521 - t479 * t524;
t520 = sin(qJ(6));
t636 = t416 * t520;
t401 = -t523 * t664 + t636;
t659 = t401 * t664;
t589 = t524 * t664;
t658 = (t516 ^ 2 + t519 ^ 2) * MDP(6) * t518 ^ 2;
t604 = t517 * t646;
t528 = t522 * t463 - t467 * t576 - t476 * t604;
t647 = pkin(3) + pkin(10);
t371 = t466 * pkin(4) + t484 * t647 + t528;
t465 = -t546 + t666;
t417 = -t467 * t517 + t644 * t476;
t553 = -qJ(4) * t466 + t417;
t376 = t465 * t647 + t553;
t657 = t521 * t371 + t524 * t376;
t538 = -t521 * t644 - t524 * t604;
t603 = t516 * t620;
t583 = t517 * t603;
t656 = -qJD(5) * t538 + t521 * t665 + t524 * t583;
t486 = -t521 * t604 + t524 * t644;
t655 = qJD(5) * t486 - t521 * t583 + t524 * t665;
t475 = (-t522 * t579 + t584) * qJD(1);
t557 = t580 - t475;
t593 = qJD(1) * t645;
t581 = pkin(1) * t593;
t482 = qJ(2) * t602 + t516 * t581;
t440 = qJD(1) * t531 + t482;
t507 = t519 * t581;
t450 = qJD(1) * t529 + t507;
t472 = qJD(1) * t476 + qJD(2);
t393 = t522 * t440 - t450 * t576 - t472 * t604;
t611 = -qJD(4) - t393;
t438 = qJD(1) * t663 - qJD(3) * t540;
t610 = pkin(4) * t654 - t611;
t359 = t479 * t647 + t610;
t411 = -t450 * t517 + t644 * t472;
t554 = -qJ(4) * t654 + t411;
t364 = t453 * t647 + t554;
t344 = t359 * t521 + t364 * t524;
t429 = t450 * t597;
t544 = qJD(1) * t559;
t382 = t519 * qJD(2) * t601 + qJD(3) * t429 + t440 * t599 + t472 * t600 + t646 * t544;
t362 = -pkin(4) * t438 + t382;
t456 = t466 * qJD(3);
t439 = qJD(1) * t456;
t619 = qJD(2) * t518;
t499 = t619 * t633;
t494 = qJD(1) * t499;
t548 = qJ(4) * t438 - qJD(4) * t654 + t494;
t370 = t439 * t647 + t548;
t591 = -t524 * t362 + t370 * t521;
t649 = -t344 * qJD(5) - t591;
t334 = pkin(5) * t438 - t649;
t653 = t413 * (pkin(5) * t416 + pkin(11) * t413) + t334;
t631 = t517 * t522;
t394 = t646 * t440 + t472 * t631 + t429;
t379 = -pkin(4) * t453 + t394;
t478 = t479 * qJ(4);
t365 = t379 - t478;
t652 = -t365 * t664 - t438 * t647;
t648 = t654 ^ 2;
t642 = qJ(4) * t453;
t615 = qJD(5) * t524;
t617 = qJD(5) * t521;
t399 = t521 * t439 + t453 * t615 + t479 * t617;
t612 = qJD(6) * t523;
t607 = t523 * t399 - t520 * t438 + t612 * t664;
t613 = qJD(6) * t520;
t352 = -t416 * t613 + t607;
t641 = t352 * t520;
t384 = pkin(3) * t453 + t554;
t640 = t384 * t654;
t639 = t401 * t413;
t403 = t416 * t523 + t520 * t664;
t638 = t403 * t413;
t637 = t413 * t647;
t635 = t453 * t654;
t634 = t654 * t521;
t428 = t524 * t439;
t400 = qJD(5) * t416 - t428;
t629 = t520 * t400;
t628 = t523 * t400;
t627 = t524 * t352;
t575 = pkin(5) * t524 + pkin(11) * t521;
t625 = (-pkin(4) - t575) * t654 - qJD(5) * t575 + t611;
t398 = t647 * t654 + t642;
t623 = t521 * t379 + t524 * t398;
t616 = qJD(5) * t523;
t614 = qJD(5) * t647;
t448 = t646 * t463;
t595 = t644 * t467;
t543 = t359 * t615 + t521 * t362 - t364 * t617 + t524 * t370;
t333 = -pkin(11) * t438 + t543;
t571 = qJD(1) * t501 - t440 * t618 + t450 * t558 + t472 * t580 - t522 * t544;
t369 = qJD(4) * t479 - t571;
t357 = -pkin(4) * t439 - t369;
t340 = pkin(5) * t400 - pkin(11) * t399 + t357;
t592 = -t333 * t520 + t523 * t340;
t590 = t399 * t520 + t523 * t438;
t588 = t664 * t416;
t587 = t664 * t403;
t586 = t413 * t523;
t503 = pkin(5) * t521 - pkin(11) * t524 + qJ(4);
t585 = -pkin(11) * t453 - qJD(6) * t503 + t623;
t391 = t484 * qJ(4) - t476 * t631 - t522 * t595 - t448;
t409 = -t453 * t520 + t523 * t634;
t573 = -t521 * t616 - t409;
t570 = t333 * t523 + t340 * t520;
t342 = pkin(11) * t664 + t344;
t350 = pkin(5) * t414 - pkin(11) * t416 + t365;
t336 = t342 * t523 + t350 * t520;
t569 = t342 * t520 - t350 * t523;
t347 = pkin(11) * t466 + t657;
t381 = -pkin(4) * t465 - t391;
t419 = t465 * t521 - t484 * t524;
t562 = t465 * t524 + t484 * t521;
t351 = -pkin(5) * t562 - pkin(11) * t419 + t381;
t568 = t347 * t523 + t351 * t520;
t567 = -t347 * t520 + t351 * t523;
t343 = t359 * t524 - t364 * t521;
t566 = t371 * t524 - t376 * t521;
t388 = qJD(2) * t532 + (t448 + (t476 * t517 + t595) * t522) * qJD(3);
t455 = -qJD(3) * t546 + t663;
t373 = -t455 * pkin(4) + t388;
t547 = qJ(4) * t455 - qJD(4) * t466 + t499;
t380 = t456 * t647 + t547;
t564 = t373 * t524 - t380 * t521;
t563 = t382 * t484 + t388 * t479;
t407 = t419 * t523 + t466 * t520;
t406 = t419 * t520 - t466 * t523;
t561 = (-qJ(2) * t603 + t507) * t516 - t482 * t519;
t556 = -t486 * t520 + t523 * t631;
t555 = t486 * t523 + t520 * t631;
t551 = -t521 * t664 ^ 2 - t524 * t438;
t550 = -t413 * t612 - t629;
t549 = t413 * t613 - t628;
t542 = t371 * t615 + t521 * t373 - t376 * t617 + t524 * t380;
t541 = -t394 * t479 - t382;
t539 = t521 * t438 - t589 * t664;
t341 = -pkin(5) * t664 - t343;
t535 = -pkin(11) * t400 + (t341 + t343) * t413;
t530 = -t453 * t479 - t438;
t363 = -t456 * pkin(4) - t383;
t412 = t438 * t466;
t410 = pkin(3) * t654 + t642;
t408 = t523 * t453 + t520 * t634;
t405 = qJD(5) * t419 - t456 * t524;
t404 = qJD(5) * t562 + t456 * t521;
t395 = pkin(3) * t456 + t547;
t392 = t484 * pkin(3) + t528;
t390 = pkin(3) * t465 + t553;
t387 = pkin(3) * t439 + t548;
t386 = t478 - t394;
t385 = pkin(3) * t479 - t611;
t356 = -qJD(6) * t406 + t404 * t523 - t455 * t520;
t355 = qJD(6) * t407 + t404 * t520 + t455 * t523;
t353 = qJD(6) * t403 + t590;
t348 = pkin(5) * t453 - t379 * t524 + t398 * t521;
t346 = -pkin(5) * t466 - t566;
t345 = t405 * pkin(5) - t404 * pkin(11) + t363;
t338 = pkin(5) * t455 + qJD(5) * t657 - t564;
t337 = -pkin(11) * t455 + t542;
t332 = -qJD(6) * t336 + t592;
t331 = -qJD(6) * t569 + t570;
t1 = [(((t519 * t622 + (qJ(2) * t632 - t511) * t516) * qJD(1) - t561) * MDP(7) - 0.2e1 * t662 * t593) * t619 + (-t400 * t562 + t405 * t413) * MDP(30) + (t399 * t562 - t400 * t419 - t404 * t414 - t405 * t416) * MDP(20) + ((-qJD(6) * t568 - t337 * t520 + t345 * t523) * t413 + t567 * t400 - t332 * t562 - t569 * t405 + t338 * t401 + t346 * t353 + t334 * t406 + t341 * t355) * MDP(31) + (-(qJD(6) * t567 + t337 * t523 + t345 * t520) * t413 - t568 * t400 + t331 * t562 - t336 * t405 + t338 * t403 + t346 * t352 + t334 * t407 + t341 * t356) * MDP(32) + (t353 * t562 - t355 * t413 - t400 * t406 - t401 * t405) * MDP(29) + (-t352 * t562 + t356 * t413 + t400 * t407 + t403 * t405) * MDP(28) + (t564 * t664 - t566 * t438 - t591 * t466 - t343 * t455 + t363 * t414 + t381 * t400 - t357 * t562 + t365 * t405 + (-t344 * t466 - t657 * t664) * qJD(5)) * MDP(24) + (t344 * t455 + t357 * t419 + t363 * t416 + t365 * t404 + t381 * t399 + t438 * t657 - t466 * t543 - t542 * t664) * MDP(25) + (-t400 * t466 - t405 * t664 + t414 * t455 - t438 * t562) * MDP(22) + (t399 * t466 + t404 * t664 - t416 * t455 - t419 * t438) * MDP(21) + (-t455 * t664 - t412) * MDP(23) + (-t384 * t456 - t387 * t465 - t390 * t439 - t395 * t453 - t563) * MDP(16) + (t439 * t484 + t456 * t479) * MDP(11) + (t438 * t484 + t455 * t479) * MDP(10) + 0.2e1 * qJD(2) * qJD(1) * t658 + (t399 * t419 + t404 * t416) * MDP(19) + (-t352 * t406 - t353 * t407 - t355 * t403 - t356 * t401) * MDP(27) + (t352 * t407 + t356 * t403) * MDP(26) + (t369 * t391 + t382 * t392 + t383 * t386 + t384 * t395 + t385 * t388 + t387 * t390) * MDP(18) + (-t411 * t455 - t417 * t438 + t466 * t494 + t479 * t527 + t484 * t571 + t499 * t654) * MDP(14) + (-t455 * t654 - t412) * MDP(8) + (t369 * t484 + t383 * t479 + t384 * t455 - t387 * t466 + t390 * t438 - t395 * t654) * MDP(17) + (t369 * t465 + t382 * t466 + t383 * t453 - t385 * t455 + t386 * t456 + t388 * t654 + t391 * t439 - t392 * t438) * MDP(15) + (t438 * t465 - t439 * t466 + t453 * t455 - t456 * t654) * MDP(9) + (t411 * t456 + t417 * t439 + (qJD(1) * t465 + t453) * t499 + t563) * MDP(13); t561 * MDP(7) * t620 + (t475 * t453 - t474 * t654 + (t646 * t438 - t439 * t522 + (-t453 * t646 + t522 * t654) * qJD(3)) * t517) * MDP(15) + (t387 * t644 - t385 * t474 + t386 * t475 + (-t384 * t603 - t646 * t382 - t369 * t522 + (t385 * t522 - t386 * t646) * qJD(3)) * t517) * MDP(18) + (t400 * t631 + t557 * t414 - t438 * t538 - t655 * t664) * MDP(24) + (t399 * t631 + t557 * t416 + t486 * t438 + t656 * t664) * MDP(25) + (-t538 * t353 + t556 * t400 + (-t555 * qJD(6) + t520 * t656 + t523 * t557) * t413 + t655 * t401) * MDP(31) + (-t538 * t352 - t555 * t400 + (-t556 * qJD(6) - t520 * t557 + t523 * t656) * t413 + t655 * t403) * MDP(32) + (MDP(13) - MDP(16)) * (t439 * t644 - t453 * t583 - t479 * t665) + (MDP(14) - MDP(17)) * (-t438 * t644 + t479 * t557 - t583 * t654) + (t518 * t645 * t662 - t658) * qJD(1) ^ 2; (-t453 ^ 2 + t648) * MDP(9) + t530 * MDP(10) - (qJD(3) + t479) * t654 * MDP(11) + (-t411 * t654 + t541) * MDP(13) + (t393 * t479 + t411 * t453 - t571) * MDP(14) + (pkin(3) * t438 - qJ(4) * t439 + (-t386 - t394) * t654 + (t385 + t611) * t453) * MDP(15) + (t410 * t453 - t541 + t640) * MDP(16) + (-t384 * t453 + t410 * t654 + t479 * t611 - t369) * MDP(17) + (-pkin(3) * t382 - qJ(4) * t369 - t384 * t410 - t385 * t394 + t386 * t611) * MDP(18) + (t399 * t524 - t521 * t588) * MDP(19) + ((-t400 - t588) * t524 + (-t399 + t667) * t521) * MDP(20) + (t416 * t453 + t551) * MDP(21) + (-t414 * t453 + t539) * MDP(22) + (qJ(4) * t400 + t343 * t453 + t610 * t414 + (t357 + (t398 + t614) * t664) * t521 + (-t379 * t664 - t652) * t524) * MDP(24) + (qJ(4) * t399 - t344 * t453 + t357 * t524 + (t524 * t614 + t623) * t664 + t610 * t416 + t652 * t521) * MDP(25) + (t523 * t627 + (-t524 * t613 + t573) * t403) * MDP(26) + (t401 * t409 + t403 * t408 + (t401 * t523 + t403 * t520) * t617 + (-t641 - t353 * t523 + (t401 * t520 - t403 * t523) * qJD(6)) * t524) * MDP(27) + (t352 * t521 + t573 * t413 + (t587 - t549) * t524) * MDP(28) + (-t353 * t521 + (t520 * t617 + t408) * t413 + (t550 - t659) * t524) * MDP(29) + (t400 * t521 + t413 * t589) * MDP(30) + (t503 * t628 - t341 * t408 - t348 * t401 + (t520 * t585 - t523 * t625) * t413 + (-t341 * t520 * qJD(5) + t332 - (qJD(5) * t401 + t550) * t647) * t521 + (t341 * t612 + t334 * t520 - t569 * t654 + t647 * t353 + (t520 * t637 - t569) * qJD(5)) * t524) * MDP(31) + (-t503 * t629 - t341 * t409 - t348 * t403 + (t520 * t625 + t523 * t585) * t413 + (-t341 * t616 - t331 - (qJD(5) * t403 + t549) * t647) * t521 + (-t341 * t613 + t334 * t523 - t336 * t654 + t647 * t352 + (t523 * t637 - t336) * qJD(5)) * t524) * MDP(32) + t664 * t453 * MDP(23) + MDP(8) * t635; t530 * MDP(15) - MDP(16) * t635 + (-t479 ^ 2 - t648) * MDP(17) + (-t386 * t479 + t382 + t640) * MDP(18) + (t414 * t479 + t551) * MDP(24) + (t416 * t479 + t539) * MDP(25) + (-t524 * t353 + (t523 * t479 - t520 * t589) * t413 + (t550 + t659) * t521) * MDP(31) + (-t627 + (-t520 * t479 - t523 * t589) * t413 + (t587 + t549) * t521) * MDP(32); -t414 ^ 2 * MDP(20) + (t399 + t667) * MDP(21) + t428 * MDP(22) - t438 * MDP(23) + (t344 * t664 + t649) * MDP(24) + (t343 * t664 + t365 * t414 - t543) * MDP(25) + (t403 * t586 + t641) * MDP(26) + ((t352 - t639) * t523 + (-t353 - t638) * t520) * MDP(27) + (t413 * t586 + t629) * MDP(28) + (-t413 ^ 2 * t520 + t628) * MDP(29) + (-pkin(5) * t353 - t344 * t401 + t535 * t520 - t523 * t653) * MDP(31) + (-pkin(5) * t352 - t344 * t403 + t520 * t653 + t535 * t523) * MDP(32) + (t414 * MDP(19) + (-qJD(5) + t664) * MDP(22) - t365 * MDP(24) - t403 * MDP(28) + t401 * MDP(29) - t413 * MDP(30) + t569 * MDP(31) + t336 * MDP(32) + t416 * MDP(20)) * t416; t403 * t401 * MDP(26) + (-t401 ^ 2 + t403 ^ 2) * MDP(27) + (t607 + t639) * MDP(28) + (-t590 + t638) * MDP(29) + t400 * MDP(30) + (t336 * t413 - t341 * t403 + t592) * MDP(31) + (t341 * t401 - t413 * t569 - t570) * MDP(32) + (-MDP(28) * t636 - MDP(29) * t403 - MDP(31) * t336 + MDP(32) * t569) * qJD(6);];
tauc  = t1;
