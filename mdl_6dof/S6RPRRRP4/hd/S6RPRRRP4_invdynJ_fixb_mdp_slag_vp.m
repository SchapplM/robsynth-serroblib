% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:59
% EndTime: 2019-03-09 06:09:10
% DurationCPUTime: 7.94s
% Computational Cost: add. (8817->498), mult. (21385->607), div. (0->0), fcn. (16944->14), ass. (0->232)
t554 = cos(qJ(5));
t630 = qJD(5) * t554;
t548 = cos(pkin(10));
t556 = cos(qJ(3));
t649 = t556 * t548;
t547 = sin(pkin(10));
t552 = sin(qJ(3));
t654 = t547 * t552;
t498 = -t649 + t654;
t488 = t498 * qJD(1);
t499 = t547 * t556 + t548 * t552;
t489 = t499 * qJD(1);
t551 = sin(qJ(4));
t555 = cos(qJ(4));
t458 = t555 * t488 + t489 * t551;
t710 = t458 * t554;
t720 = t630 + t710;
t550 = sin(qJ(5));
t705 = qJD(5) + t458;
t714 = t705 * t550;
t719 = pkin(5) * t714;
t579 = -t488 * t551 + t555 * t489;
t622 = qJD(3) + qJD(4);
t450 = t550 * t622 + t554 * t579;
t718 = t450 * t714;
t546 = pkin(10) + qJ(3);
t539 = qJ(4) + t546;
t529 = sin(t539);
t553 = sin(qJ(1));
t557 = cos(qJ(1));
t587 = g(1) * t557 + g(2) * t553;
t717 = t587 * t529;
t611 = qJD(1) * t654;
t623 = qJDD(1) * t556;
t624 = qJDD(1) * t552;
t634 = qJD(3) * t556;
t613 = t547 * t623 + (qJD(1) * t634 + t624) * t548;
t465 = -qJD(3) * t611 + t613;
t491 = t499 * qJD(3);
t583 = t547 * t624 - t548 * t623;
t466 = qJD(1) * t491 + t583;
t632 = qJD(4) * t555;
t633 = qJD(4) * t551;
t569 = -t555 * t465 + t551 * t466 + t488 * t632 + t489 * t633;
t621 = qJDD(3) + qJDD(4);
t563 = -t550 * t569 - t554 * t621;
t391 = qJD(5) * t450 + t563;
t593 = t554 * t622;
t448 = t550 * t579 - t593;
t670 = pkin(7) + qJ(2);
t513 = t670 * t547;
t500 = qJD(1) * t513;
t514 = t670 * t548;
t501 = qJD(1) * t514;
t578 = t500 * t552 - t501 * t556;
t625 = qJD(1) * qJD(2);
t681 = qJDD(1) * t670 + t625;
t475 = t681 * t547;
t476 = t681 * t548;
t602 = -t556 * t475 - t552 * t476;
t402 = qJDD(3) * pkin(3) - pkin(8) * t465 + qJD(3) * t578 + t602;
t580 = -t552 * t475 + t556 * t476;
t696 = -t556 * t500 - t501 * t552;
t405 = -pkin(8) * t466 + qJD(3) * t696 + t580;
t446 = -pkin(8) * t489 + t696;
t444 = qJD(3) * pkin(3) + t446;
t447 = -pkin(8) * t488 - t578;
t682 = -t555 * (qJD(4) * t444 + t405) - t551 * t402 + t447 * t633;
t368 = pkin(9) * t621 - t682;
t603 = t551 * t465 + t555 * t466;
t690 = t579 * qJD(4);
t406 = t603 + t690;
t612 = pkin(2) * t548 + pkin(1);
t506 = -qJDD(1) * t612 + qJDD(2);
t452 = t466 * pkin(3) + t506;
t377 = t406 * pkin(4) + pkin(9) * t569 + t452;
t443 = t555 * t447;
t414 = t551 * t444 + t443;
t408 = pkin(9) * t622 + t414;
t507 = -qJD(1) * t612 + qJD(2);
t471 = pkin(3) * t488 + t507;
t415 = pkin(4) * t458 - pkin(9) * t579 + t471;
t631 = qJD(5) * t550;
t568 = t554 * t368 + t550 * t377 - t408 * t631 + t415 * t630;
t362 = -qJ(6) * t391 - qJD(6) * t448 + t568;
t523 = g(3) * t529;
t530 = cos(t539);
t614 = t530 * t587 + t523;
t373 = t554 * t377;
t382 = t408 * t554 + t415 * t550;
t390 = -qJD(5) * t593 - t550 * t621 + t554 * t569 + t579 * t631;
t404 = qJDD(5) + t406;
t360 = pkin(5) * t404 + qJ(6) * t390 - qJD(5) * t382 - qJD(6) * t450 - t368 * t550 + t373;
t375 = -qJ(6) * t448 + t382;
t687 = t375 * t705 + t360;
t716 = t362 * t554 - t550 * t687 - t614;
t442 = t551 * t447;
t417 = t446 * t555 - t442;
t709 = -pkin(3) * t632 + t417;
t389 = t554 * t390;
t646 = -t550 * t391 - t448 * t630;
t713 = -t448 * t710 - t389 + t646;
t388 = t390 * t550;
t397 = t550 * t404;
t628 = t579 * qJD(3);
t663 = t450 * t579;
t698 = t458 * t622;
t712 = t621 * MDP(19) + (t628 - t603) * MDP(18) - t458 ^ 2 * MDP(16) + (MDP(15) * t458 + t579 * MDP(16) - t705 * MDP(26)) * t579 + (-t569 + t698) * MDP(17) + (t450 * t720 - t388) * MDP(22) + (t705 * t720 + t397 - t663) * MDP(24);
t413 = t555 * t444 - t442;
t407 = -pkin(4) * t622 - t413;
t711 = t407 * t458;
t538 = cos(t546);
t534 = pkin(5) * t554 + pkin(4);
t549 = -qJ(6) - pkin(9);
t600 = -t529 * t549 + t530 * t534;
t707 = pkin(3) * t538 + t600;
t667 = qJDD(1) * pkin(1);
t695 = g(1) * t553 - g(2) * t557;
t575 = -qJDD(2) + t667 + t695;
t431 = pkin(4) * t579 + pkin(9) * t458;
t704 = t471 * t458 + t614 + t682;
t541 = t554 * qJ(6);
t703 = -pkin(5) * t579 - t458 * t541;
t664 = t448 * t579;
t697 = t695 * t529;
t636 = -t552 * t513 + t556 * t514;
t420 = pkin(3) * t489 + t431;
t694 = t550 * t420 + t554 * t709;
t692 = -MDP(4) * t548 + MDP(5) * t547;
t650 = t554 * t557;
t653 = t550 * t553;
t480 = t530 * t653 + t650;
t651 = t553 * t554;
t652 = t550 * t557;
t482 = -t530 * t652 + t651;
t691 = -g(1) * t482 + g(2) * t480;
t689 = qJ(2) * qJDD(1);
t673 = g(3) * t530;
t688 = -t673 + t717;
t381 = -t408 * t550 + t554 * t415;
t686 = -t381 * t579 + t407 * t631 + t554 * t717;
t590 = -t555 * t402 + t551 * t405 + t444 * t633 + t447 * t632;
t369 = -pkin(4) * t621 + t590;
t400 = t407 * t630;
t672 = g(3) * t550;
t685 = t369 * t550 + t382 * t579 + t530 * t672 + t400;
t684 = -t471 * t579 - t590 + t688;
t680 = t450 ^ 2;
t679 = pkin(3) * t555;
t671 = t491 * pkin(3);
t468 = -t498 * t551 + t499 * t555;
t666 = t369 * t468;
t374 = -qJ(6) * t450 + t381;
t370 = pkin(5) * t705 + t374;
t665 = t370 * t554;
t662 = t450 * t550;
t659 = t458 * t550;
t657 = t468 * t550;
t398 = t554 * t404;
t601 = -t556 * t513 - t514 * t552;
t454 = -pkin(8) * t499 + t601;
t455 = -pkin(8) * t498 + t636;
t426 = t454 * t551 + t455 * t555;
t422 = t554 * t426;
t533 = pkin(3) * t551 + pkin(9);
t648 = -qJ(6) - t533;
t647 = -t374 + t370;
t644 = t554 * t413 + t550 * t431;
t467 = t555 * t498 + t499 * t551;
t473 = pkin(3) * t498 - t612;
t427 = pkin(4) * t467 - pkin(9) * t468 + t473;
t642 = t550 * t427 + t422;
t540 = t554 * qJD(6);
t599 = qJD(5) * t648;
t640 = -qJ(6) * t659 + t550 * t599 + t540 - t694;
t419 = t554 * t420;
t639 = t554 * t599 - t419 + (-qJD(6) + t709) * t550 + t703;
t606 = qJD(5) * t549;
t638 = t540 - t644 + (-qJ(6) * t458 + t606) * t550;
t429 = t554 * t431;
t637 = t554 * t606 - t429 + (-qJD(6) + t413) * t550 + t703;
t635 = t547 ^ 2 + t548 ^ 2;
t619 = qJD(5) * pkin(9) * t705;
t564 = -t513 * t634 + qJD(2) * t649 + (-qJD(2) * t547 - qJD(3) * t514) * t552;
t434 = -pkin(8) * t491 + t564;
t490 = t498 * qJD(3);
t560 = -t499 * qJD(2) - qJD(3) * t636;
t435 = pkin(8) * t490 + t560;
t581 = t454 * t555 - t455 * t551;
t385 = qJD(4) * t581 + t434 * t555 + t435 * t551;
t432 = -qJD(4) * t467 - t490 * t555 - t491 * t551;
t433 = qJD(4) * t468 - t490 * t551 + t555 * t491;
t394 = pkin(4) * t433 - pkin(9) * t432 + t671;
t617 = t554 * t385 + t550 * t394 + t427 * t630;
t610 = t468 * t630;
t609 = pkin(5) * t550 + pkin(8) + t670;
t607 = -t369 - t673;
t605 = t635 * qJD(1) ^ 2;
t597 = t705 * t554;
t594 = -qJD(5) * t415 - t368;
t589 = 0.2e1 * t635;
t416 = t446 * t551 + t443;
t588 = pkin(3) * t633 - t416;
t585 = -t408 * t630 + t373;
t584 = -pkin(9) * t404 + t711;
t582 = -t404 * t533 + t711;
t577 = t529 * t534 + t530 * t549;
t576 = -qJ(6) * t432 - qJD(6) * t468;
t574 = t398 + (-t631 - t659) * t705;
t572 = -t612 - t707;
t571 = t432 * t550 + t610;
t570 = t432 * t554 - t468 * t631;
t562 = t589 * t625 - t587;
t365 = t391 * pkin(5) + qJDD(6) + t369;
t386 = qJD(4) * t426 + t434 * t551 - t435 * t555;
t537 = sin(t546);
t535 = -pkin(4) - t679;
t516 = pkin(9) * t554 + t541;
t515 = t549 * t550;
t494 = t533 * t554 + t541;
t493 = t648 * t550;
t483 = t530 * t650 + t653;
t481 = -t530 * t651 + t652;
t445 = t448 ^ 2;
t424 = t554 * t427;
t395 = t448 * pkin(5) + qJD(6) + t407;
t393 = t554 * t394;
t384 = -qJ(6) * t657 + t642;
t378 = pkin(5) * t467 - t426 * t550 - t468 * t541 + t424;
t364 = -qJ(6) * t610 + (-qJD(5) * t426 + t576) * t550 + t617;
t363 = pkin(5) * t433 - t385 * t550 + t393 + t576 * t554 + (-t422 + (qJ(6) * t468 - t427) * t550) * qJD(5);
t1 = [(-t386 * t622 + t473 * t406 + t471 * t433 + t452 * t467 + t458 * t671 + t530 * t695 + t581 * t621) * MDP(20) + (qJD(3) * t560 + qJDD(3) * t601 - t466 * t612 + t507 * t491 + t506 * t498 + t538 * t695) * MDP(13) + (-qJD(3) * t564 - qJDD(3) * t636 - t465 * t612 - t507 * t490 + t506 * t499 - t537 * t695) * MDP(14) + t695 * MDP(2) + (-t385 * t622 - t426 * t621 + t471 * t432 + t452 * t468 - t473 * t569 + t579 * t671 - t697) * MDP(21) + (-t363 * t450 - t364 * t448 + t378 * t390 - t384 * t391 + t697 + (-t375 * t550 - t665) * t432 + (-t360 * t554 - t362 * t550 + (t370 * t550 - t375 * t554) * qJD(5)) * t468) * MDP(29) - t692 * (t575 + t667) + (pkin(1) * t575 + (t635 * t689 + t562) * qJ(2)) * MDP(7) + (t589 * t689 + t562) * MDP(6) + (-(-t426 * t631 + t617) * t705 - t642 * t404 - t568 * t467 - t382 * t433 + t386 * t450 + t581 * t390 + t554 * t666 - g(1) * t480 - g(2) * t482 + t570 * t407) * MDP(28) + ((-t426 * t630 + t393) * t705 + t424 * t404 + t585 * t467 + t381 * t433 + t386 * t448 - t581 * t391 + t468 * t400 - g(1) * t481 - g(2) * t483 + ((-qJD(5) * t427 - t385) * t705 - t426 * t404 + t594 * t467 + t666 + t407 * t432) * t550) * MDP(27) + (-t391 * t467 - t397 * t468 - t433 * t448 - t571 * t705) * MDP(25) + (-t390 * t467 + t398 * t468 + t433 * t450 + t570 * t705) * MDP(24) + (t404 * t467 + t433 * t705) * MDP(26) + (t465 * t499 - t489 * t490) * MDP(8) + (-qJD(3) * t490 + qJDD(3) * t499) * MDP(10) + (-t468 * t406 - t432 * t458 - t433 * t579 + t467 * t569) * MDP(16) + (t432 * t579 - t468 * t569) * MDP(15) + ((-t448 * t554 - t662) * t432 + (t388 - t391 * t554 + (t448 * t550 - t450 * t554) * qJD(5)) * t468) * MDP(23) + (t362 * t384 + t375 * t364 + t360 * t378 + t370 * t363 + t365 * (pkin(5) * t657 - t581) + t395 * (pkin(5) * t571 + t386) + (-g(1) * t609 + g(2) * t572) * t557 + (-g(1) * t572 - g(2) * t609) * t553) * MDP(30) + (-t389 * t468 + t450 * t570) * MDP(22) + (-t433 * t622 - t467 * t621) * MDP(18) + (t432 * t622 + t468 * t621) * MDP(17) + (-qJD(3) * t491 - qJDD(3) * t498) * MDP(11) + (-t465 * t498 - t466 * t499 + t488 * t490 - t489 * t491) * MDP(9) + t587 * MDP(3) + qJDD(1) * MDP(1); -MDP(6) * t605 + (-qJ(2) * t605 - t575) * MDP(7) + (0.2e1 * qJD(3) * t489 + t583) * MDP(13) + ((-t488 - t611) * qJD(3) + t613) * MDP(14) + (t603 + t628 + 0.2e1 * t690) * MDP(20) + (-t569 - t698) * MDP(21) + (t574 - t664) * MDP(27) + (-t597 * t705 - t397 - t663) * MDP(28) + ((-t448 * t458 + t390) * t554 + t718 + t646) * MDP(29) + (-t395 * t579 + t687 * t554 + (-t370 * t705 + t362) * t550 - t695) * MDP(30) + t692 * qJDD(1); (t417 * t622 + (-t489 * t579 - t551 * t621 - t622 * t632) * pkin(3) + t704) * MDP(21) - t583 * MDP(11) + (t362 * t494 + t360 * t493 + t365 * (-t534 - t679) - g(3) * t707 + (-t443 + (pkin(3) * qJD(4) - t446) * t551 + t719) * t395 + t640 * t375 + t639 * t370 + t587 * (pkin(3) * t537 + t577)) * MDP(30) + (-t535 * t390 + t582 * t554 - t550 * t717 + t588 * t450 + (t533 * t631 + t694) * t705 + t685) * MDP(28) + (t535 * t391 + t607 * t554 + t582 * t550 + t588 * t448 + (-t533 * t630 + t550 * t709 - t419) * t705 + t686) * MDP(27) + (-t370 * t597 + t390 * t493 - t391 * t494 - t640 * t448 - t639 * t450 + t716) * MDP(29) + (t416 * t622 + (-t489 * t458 + t555 * t621 - t622 * t633) * pkin(3) + t684) * MDP(20) + (-g(3) * t538 - t507 * t489 + t537 * t587 + t602) * MDP(13) + (-t662 * t705 + t713) * MDP(23) + ((t488 - t611) * qJD(3) + t613) * MDP(10) + (g(3) * t537 + t507 * t488 + t538 * t587 - t580) * MDP(14) + (t574 + t664) * MDP(25) + (-t488 ^ 2 + t489 ^ 2) * MDP(9) + qJDD(3) * MDP(12) + t488 * t489 * MDP(8) + t712; (t414 * t622 + t684) * MDP(20) + (t413 * t622 + t704) * MDP(21) + (t713 - t718) * MDP(23) + (-t705 * t714 + t398 + t664) * MDP(25) + (-pkin(4) * t391 - t414 * t448 - t429 * t705 + (t413 * t705 + t584) * t550 + (t607 - t619) * t554 + t686) * MDP(27) + (pkin(4) * t390 + t644 * t705 - t414 * t450 + t584 * t554 + (-t717 + t619) * t550 + t685) * MDP(28) + (t390 * t515 - t391 * t516 - t638 * t448 - t637 * t450 - t705 * t665 + t716) * MDP(29) + (t362 * t516 + t360 * t515 - t365 * t534 - g(3) * t600 + (-t414 + t719) * t395 + t638 * t375 + t637 * t370 + t587 * t577) * MDP(30) + t712; t450 * t448 * MDP(22) + (-t445 + t680) * MDP(23) + (t448 * t705 - t390) * MDP(24) + (-t563 + (-qJD(5) + t705) * t450) * MDP(25) + t404 * MDP(26) + (t382 * t705 - t407 * t450 + (t594 + t523) * t550 + t585 + t691) * MDP(27) + (g(1) * t483 - g(2) * t481 + t381 * t705 + t407 * t448 + t523 * t554 - t568) * MDP(28) + (pkin(5) * t390 - t448 * t647) * MDP(29) + (t647 * t375 + (-t395 * t450 + t529 * t672 + t360 + t691) * pkin(5)) * MDP(30); (-t445 - t680) * MDP(29) + (t370 * t450 + t375 * t448 + t365 - t688) * MDP(30);];
tau  = t1;
