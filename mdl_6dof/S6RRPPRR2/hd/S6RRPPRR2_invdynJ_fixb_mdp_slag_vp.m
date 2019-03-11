% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:55
% EndTime: 2019-03-09 08:53:11
% DurationCPUTime: 11.90s
% Computational Cost: add. (8954->588), mult. (21242->766), div. (0->0), fcn. (16901->18), ass. (0->242)
t616 = cos(qJ(2));
t702 = cos(pkin(10));
t655 = t702 * t616;
t584 = qJD(1) * t655;
t607 = sin(pkin(10));
t612 = sin(qJ(2));
t673 = qJD(1) * t612;
t554 = t607 * t673 - t584;
t549 = qJD(5) + t554;
t570 = t607 * t616 + t612 * t702;
t557 = t570 * qJD(1);
t606 = sin(pkin(11));
t608 = cos(pkin(11));
t527 = qJD(2) * t606 + t557 * t608;
t528 = qJD(2) * t608 - t557 * t606;
t611 = sin(qJ(5));
t615 = cos(qJ(5));
t464 = t527 * t611 - t528 * t615;
t614 = cos(qJ(6));
t610 = sin(qJ(6));
t720 = -t527 * t615 - t528 * t611;
t694 = t720 * t610;
t410 = t464 * t614 - t694;
t543 = qJD(6) + t549;
t724 = t410 * t543;
t728 = t464 * t549;
t636 = t464 * t610 + t614 * t720;
t723 = t543 * t636;
t569 = t606 * t611 - t608 * t615;
t677 = t549 * t569;
t571 = t606 * t615 + t608 * t611;
t561 = t571 * qJD(5);
t676 = t554 * t571 + t561;
t727 = t549 * t720;
t491 = pkin(2) * t673 + pkin(3) * t557 + qJ(4) * t554;
t609 = -qJ(3) - pkin(7);
t580 = t609 * t616;
t576 = qJD(1) * t580;
t562 = t607 * t576;
t579 = t609 * t612;
t575 = qJD(1) * t579;
t518 = t575 * t702 + t562;
t446 = t491 * t606 + t518 * t608;
t689 = t554 * t606;
t428 = pkin(8) * t689 + t446;
t725 = -qJD(4) * t608 + t428;
t657 = qJD(2) * t609;
t551 = -qJD(3) * t612 + t616 * t657;
t505 = qJDD(2) * pkin(2) + qJD(1) * t551 + qJDD(1) * t579;
t550 = qJD(3) * t616 + t612 * t657;
t516 = qJD(1) * t550 - qJDD(1) * t580;
t448 = t505 * t702 - t516 * t607;
t447 = -qJDD(2) * pkin(3) + qJDD(4) - t448;
t603 = qJ(2) + pkin(10);
t597 = sin(t603);
t599 = cos(t603);
t613 = sin(qJ(1));
t617 = cos(qJ(1));
t646 = g(1) * t617 + g(2) * t613;
t626 = -g(3) * t599 + t597 * t646;
t717 = t447 - t626;
t515 = -t569 * t610 + t571 * t614;
t680 = qJD(6) * t515 - t610 * t677 + t614 * t676;
t712 = pkin(2) * t616;
t593 = pkin(1) + t712;
t578 = -qJD(1) * t593 + qJD(3);
t480 = pkin(3) * t554 - qJ(4) * t557 + t578;
t703 = qJD(2) * pkin(2);
t567 = t575 + t703;
t656 = t702 * t576;
t510 = t567 * t607 - t656;
t501 = qJD(2) * qJ(4) + t510;
t433 = t480 * t608 - t501 * t606;
t404 = pkin(4) * t554 - pkin(8) * t527 + t433;
t434 = t480 * t606 + t501 * t608;
t417 = pkin(8) * t528 + t434;
t386 = t404 * t611 + t417 * t615;
t376 = -pkin(9) * t464 + t386;
t669 = qJD(6) * t610;
t374 = t376 * t669;
t509 = t567 * t702 + t562;
t494 = -qJD(2) * pkin(3) + qJD(4) - t509;
t457 = -pkin(4) * t528 + t494;
t408 = pkin(5) * t464 + t457;
t602 = pkin(11) + qJ(5);
t600 = qJ(6) + t602;
t590 = sin(t600);
t591 = cos(t600);
t685 = t599 * t613;
t530 = t590 * t617 - t591 * t685;
t684 = t599 * t617;
t532 = t590 * t613 + t591 * t684;
t707 = g(3) * t597;
t722 = g(1) * t532 - g(2) * t530 + t408 * t410 + t591 * t707 + t374;
t529 = t590 * t685 + t591 * t617;
t531 = -t590 * t684 + t591 * t613;
t666 = qJD(1) * qJD(2);
t660 = t612 * t666;
t512 = qJD(2) * t584 + qJDD(1) * t570 - t607 * t660;
t487 = -qJDD(2) * t608 + t512 * t606;
t488 = qJDD(2) * t606 + t512 * t608;
t670 = qJD(5) * t615;
t671 = qJD(5) * t611;
t396 = -t487 * t611 + t488 * t615 - t527 * t671 + t528 * t670;
t556 = t570 * qJD(2);
t665 = qJDD(1) * t612;
t511 = qJD(1) * t556 - qJDD(1) * t655 + t607 * t665;
t506 = qJDD(5) + t511;
t663 = pkin(2) * t660 + qJDD(3);
t627 = -qJDD(1) * t593 + t663;
t420 = pkin(3) * t511 - qJ(4) * t512 - qJD(4) * t557 + t627;
t449 = t505 * t607 + t516 * t702;
t444 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t449;
t391 = t420 * t608 - t444 * t606;
t381 = pkin(4) * t511 - pkin(8) * t488 + t391;
t392 = t420 * t606 + t444 * t608;
t389 = -pkin(8) * t487 + t392;
t653 = t381 * t615 - t389 * t611;
t620 = -qJD(5) * t386 + t653;
t365 = pkin(5) * t506 - pkin(9) * t396 + t620;
t397 = -qJD(5) * t720 + t487 * t615 + t488 * t611;
t630 = t381 * t611 + t389 * t615 + t404 * t670 - t417 * t671;
t366 = -pkin(9) * t397 + t630;
t654 = t365 * t614 - t610 * t366;
t721 = -g(1) * t531 + g(2) * t529 + t408 * t636 + t590 * t707 + t654;
t502 = qJDD(6) + t506;
t719 = t502 * MDP(28) + (-t410 ^ 2 + t636 ^ 2) * MDP(25) - t410 * MDP(24) * t636;
t716 = g(1) * t613 - g(2) * t617;
t718 = t599 * t716;
t632 = -t607 * t612 + t655;
t508 = -pkin(3) * t632 - qJ(4) * t570 - t593;
t523 = t579 * t607 - t580 * t702;
t455 = t508 * t608 - t523 * t606;
t686 = t570 * t608;
t427 = -pkin(4) * t632 - pkin(8) * t686 + t455;
t456 = t508 * t606 + t523 * t608;
t687 = t570 * t606;
t437 = -pkin(8) * t687 + t456;
t682 = t427 * t611 + t437 * t615;
t587 = pkin(2) * t607 + qJ(4);
t704 = pkin(8) + t587;
t565 = t704 * t606;
t566 = t704 * t608;
t678 = -t565 * t611 + t566 * t615;
t514 = t569 * t614 + t571 * t610;
t681 = -qJD(6) * t514 - t610 * t676 - t614 * t677;
t715 = -t502 * t515 - t543 * t681;
t714 = -t506 * t571 + t549 * t677;
t652 = t610 * t396 + t397 * t614;
t372 = -qJD(6) * t636 + t652;
t445 = t491 * t608 - t518 * t606;
t711 = pkin(8) * t608;
t416 = pkin(4) * t557 + t554 * t711 + t445;
t634 = qJD(4) * t606 + qJD(5) * t566;
t713 = t565 * t670 + t725 * t615 + (t416 + t634) * t611;
t705 = g(3) * t616;
t701 = qJDD(1) * pkin(1);
t385 = t404 * t615 - t417 * t611;
t375 = pkin(9) * t720 + t385;
t373 = pkin(5) * t549 + t375;
t700 = t373 * t614;
t699 = t376 * t614;
t698 = t410 * t557;
t697 = t636 * t557;
t696 = t464 * t557;
t695 = t720 * t557;
t691 = t511 * t606;
t690 = t511 * t608;
t559 = t632 * qJD(2);
t688 = t559 * t606;
t662 = t612 * t703;
t469 = pkin(3) * t556 - qJ(4) * t559 - qJD(4) * t570 + t662;
t493 = t550 * t702 + t551 * t607;
t423 = t469 * t606 + t493 * t608;
t604 = t612 ^ 2;
t675 = -t616 ^ 2 + t604;
t668 = qJD(6) * t614;
t667 = -qJD(4) + t494;
t664 = qJDD(1) * t616;
t661 = t396 * t614 - t397 * t610 - t464 * t668;
t517 = t575 * t607 - t656;
t470 = -pkin(4) * t689 + t517;
t659 = pkin(5) * t676 - t470;
t422 = t469 * t608 - t493 * t606;
t401 = pkin(4) * t556 - t559 * t711 + t422;
t406 = -pkin(8) * t688 + t423;
t651 = t401 * t615 - t406 * t611;
t650 = t427 * t615 - t437 * t611;
t492 = t550 * t607 - t551 * t702;
t648 = -t565 * t615 - t566 * t611;
t522 = -t579 * t702 - t580 * t607;
t647 = qJD(6) * t373 + t366;
t592 = -pkin(2) * t702 - pkin(3);
t644 = -t502 * t514 - t543 * t680;
t643 = -t569 * t506 - t549 * t676;
t459 = pkin(4) * t688 + t492;
t486 = pkin(4) * t687 + t522;
t415 = t615 * t416;
t472 = -pkin(9) * t569 + t678;
t642 = pkin(5) * t557 - pkin(9) * t677 + qJD(4) * t571 + qJD(5) * t678 + qJD(6) * t472 - t428 * t611 + t415;
t471 = -pkin(9) * t571 + t648;
t641 = pkin(9) * t676 - qJD(6) * t471 + t713;
t640 = pkin(3) * t599 + qJ(4) * t597;
t368 = t373 * t610 + t699;
t639 = -t391 * t608 - t392 * t606;
t637 = -t433 * t606 + t434 * t608;
t495 = t571 * t570;
t496 = t569 * t570;
t441 = t495 * t614 - t496 * t610;
t442 = -t495 * t610 - t496 * t614;
t633 = -0.2e1 * pkin(1) * t666 - pkin(7) * qJDD(2);
t577 = -pkin(4) * t608 + t592;
t629 = t401 * t611 + t406 * t615 + t427 * t670 - t437 * t671;
t371 = t669 * t720 + t661;
t618 = qJD(2) ^ 2;
t623 = -pkin(7) * t618 + 0.2e1 * t701 + t716;
t619 = qJD(1) ^ 2;
t622 = pkin(1) * t619 - pkin(7) * qJDD(1) + t646;
t621 = t447 * t570 + t494 * t559 - t646;
t407 = pkin(4) * t487 + t447;
t598 = cos(t602);
t596 = sin(t602);
t582 = t617 * t593;
t552 = t554 ^ 2;
t536 = t596 * t613 + t598 * t684;
t535 = -t596 * t684 + t598 * t613;
t534 = t596 * t617 - t598 * t685;
t533 = t596 * t685 + t598 * t617;
t524 = pkin(5) * t569 + t577;
t443 = pkin(5) * t495 + t486;
t440 = t559 * t571 + t670 * t686 - t671 * t687;
t439 = -t559 * t569 - t561 * t570;
t398 = pkin(5) * t440 + t459;
t390 = -pkin(9) * t495 + t682;
t387 = -pkin(5) * t632 + pkin(9) * t496 + t650;
t384 = qJD(6) * t442 + t439 * t610 + t440 * t614;
t383 = -qJD(6) * t441 + t439 * t614 - t440 * t610;
t380 = pkin(5) * t397 + t407;
t370 = -pkin(9) * t440 + t629;
t369 = pkin(5) * t556 - pkin(9) * t439 - qJD(5) * t682 + t651;
t367 = -t376 * t610 + t700;
t1 = [(-t371 * t441 - t372 * t442 - t383 * t410 + t384 * t636) * MDP(25) + (t371 * t442 - t383 * t636) * MDP(24) + (-g(1) * t529 - g(2) * t531 - t368 * t556 + t443 * t371 - t374 * t632 + t380 * t442 + t408 * t383 - t398 * t636 + (-(-qJD(6) * t390 + t369) * t543 - t387 * t502 + t365 * t632) * t610 + (-(qJD(6) * t387 + t370) * t543 - t390 * t502 + t647 * t632) * t614) * MDP(30) + (-t371 * t632 + t383 * t543 + t442 * t502 - t556 * t636) * MDP(26) + (-t391 * t632 + t422 * t554 + t433 * t556 + t455 * t511 + t522 * t487 - t492 * t528 + t606 * t621 + t608 * t718) * MDP(13) + (t392 * t632 - t423 * t554 - t434 * t556 - t456 * t511 + t522 * t488 + t492 * t527 - t606 * t718 + t608 * t621) * MDP(14) + ((t369 * t614 - t370 * t610) * t543 + (t387 * t614 - t390 * t610) * t502 - t654 * t632 + t367 * t556 + t398 * t410 + t443 * t372 + t380 * t441 + t408 * t384 - g(1) * t530 - g(2) * t532 + ((-t387 * t610 - t390 * t614) * t543 + t368 * t632) * qJD(6)) * MDP(29) + (t397 * t632 - t440 * t549 - t464 * t556 - t495 * t506) * MDP(20) + (-t506 * t632 + t549 * t556) * MDP(21) + (t372 * t632 - t384 * t543 - t410 * t556 - t441 * t502) * MDP(27) + (-t502 * t632 + t543 * t556) * MDP(28) + (-t448 * t570 + t449 * t632 + t492 * t557 - t493 * t554 - t509 * t559 - t510 * t556 - t511 * t523 + t512 * t522 - t646) * MDP(11) + (-t422 * t527 + t423 * t528 - t455 * t488 - t456 * t487 + t716 * t597 + t639 * t570 + (-t433 * t608 - t434 * t606) * t559) * MDP(15) + t716 * MDP(2) + 0.2e1 * (t612 * t664 - t666 * t675) * MDP(5) + (t449 * t523 + t510 * t493 - t448 * t522 - t509 * t492 - t627 * t593 + t578 * t662 - g(1) * (-t593 * t613 - t609 * t617) - g(2) * (-t609 * t613 + t582)) * MDP(12) + (qJDD(1) * t604 + 0.2e1 * t616 * t660) * MDP(4) + t646 * MDP(3) + (-g(2) * t582 + t391 * t455 + t392 * t456 + t433 * t422 + t434 * t423 + t447 * t522 + t494 * t492 + (g(1) * t609 - g(2) * t640) * t617 + (-g(1) * (-t593 - t640) + g(2) * t609) * t613) * MDP(16) + (t612 * t633 + t616 * t623) * MDP(9) + (-t612 * t623 + t616 * t633) * MDP(10) + (-t396 * t495 + t397 * t496 - t439 * t464 + t440 * t720) * MDP(18) + (-t396 * t496 - t439 * t720) * MDP(17) + (-g(1) * t533 - g(2) * t535 - t386 * t556 + t486 * t396 - t407 * t496 + t457 * t439 - t459 * t720 - t506 * t682 - t549 * t629 + t630 * t632) * MDP(23) + (-t396 * t632 + t439 * t549 - t496 * t506 - t556 * t720) * MDP(19) + (t651 * t549 + t650 * t506 - t653 * t632 + t385 * t556 + t459 * t464 + t486 * t397 + t407 * t495 + t457 * t440 - g(1) * t534 - g(2) * t536 + (t386 * t632 - t549 * t682) * qJD(5)) * MDP(22) + (qJDD(2) * t612 + t616 * t618) * MDP(6) + (qJDD(2) * t616 - t612 * t618) * MDP(7) + qJDD(1) * MDP(1); (t509 * t517 - t510 * t518 + (t702 * t448 - t705 + t449 * t607 + (-qJD(1) * t578 + t646) * t612) * pkin(2)) * MDP(12) + (t612 * t622 - t705) * MDP(9) + ((-t509 + t518) * t554 + (-t511 * t607 - t512 * t702) * pkin(2)) * MDP(11) + ((t471 * t614 - t472 * t610) * t502 + t524 * t372 + t380 * t514 + (t610 * t641 - t614 * t642) * t543 + t659 * t410 + t680 * t408 + t626 * t591) * MDP(29) + (t371 * t515 - t636 * t681) * MDP(24) + (-t371 * t514 - t372 * t515 - t410 * t681 + t636 * t680) * MDP(25) + (-t707 + t445 * t527 - t446 * t528 - t646 * t599 + (qJD(4) * t528 - t433 * t554 - t487 * t587 + t392) * t608 + (qJD(4) * t527 - t434 * t554 + t488 * t587 - t391) * t606) * MDP(15) + (t386 * MDP(23) + t434 * MDP(14) + (t510 - t517) * MDP(11) - t367 * MDP(29) + t368 * MDP(30) - t433 * MDP(13) - t385 * MDP(22) - t543 * MDP(28) - t549 * MDP(21)) * t557 + (t695 - t714) * MDP(19) + (t697 - t715) * MDP(26) + (-t587 * t691 + t487 * t592 + t517 * t528 + (t606 * t667 - t445) * t554 - t717 * t608) * MDP(13) + (t648 * t506 + t577 * t397 + t407 * t569 - t470 * t464 + (-t415 - t634 * t615 + (qJD(5) * t565 + t725) * t611) * t549 + t676 * t457 + t626 * t598) * MDP(22) + MDP(7) * t664 + MDP(6) * t665 + (g(3) * t612 + t616 * t622) * MDP(10) + (t396 * t571 + t677 * t720) * MDP(17) + (-t396 * t569 - t397 * t571 + t464 * t677 + t676 * t720) * MDP(18) + (t577 * t396 + t407 * t571 - t677 * t457 + t470 * t720 - t678 * t506 + t549 * t713 - t596 * t626) * MDP(23) + (t447 * t592 - t434 * t446 - t433 * t445 - t494 * t517 - g(3) * (t640 + t712) + (-t391 * t606 + t392 * t608) * t587 + t637 * qJD(4) + t646 * (pkin(2) * t612 + pkin(3) * t597 - qJ(4) * t599)) * MDP(16) + (-MDP(4) * t612 * t616 + MDP(5) * t675) * t619 + (-t587 * t690 + t488 * t592 - t517 * t527 + (t608 * t667 + t446) * t554 + t717 * t606) * MDP(14) + (-(t471 * t610 + t472 * t614) * t502 + t524 * t371 + t380 * t515 + (t610 * t642 + t614 * t641) * t543 - t659 * t636 + t681 * t408 - t626 * t590) * MDP(30) + (t643 + t696) * MDP(20) + (t644 + t698) * MDP(27) + qJDD(2) * MDP(8); (-t557 ^ 2 - t552) * MDP(11) + (-pkin(2) * t664 + t509 * t557 + t663 - t701 - t716) * MDP(12) + (t528 * t557 + t690) * MDP(13) + (-t527 * t557 - t552 * t608 - t691) * MDP(14) + (-t487 * t606 - t488 * t608) * MDP(15) + (-t494 * t557 - t639 - t716) * MDP(16) + (t643 - t696) * MDP(22) + (t695 + t714) * MDP(23) + (t644 - t698) * MDP(29) + (t697 + t715) * MDP(30) + (t510 * MDP(12) + (t527 * t606 + t528 * t608) * MDP(15) + t637 * MDP(16) - MDP(13) * t689) * t554; (t527 * t554 + t487) * MDP(13) + (t528 * t554 + t488) * MDP(14) + (-t527 ^ 2 - t528 ^ 2) * MDP(15) + (t433 * t527 - t434 * t528 + t717) * MDP(16) + (t397 - t727) * MDP(22) + (t396 - t728) * MDP(23) + (t372 - t723) * MDP(29) + (t371 - t724) * MDP(30); -t720 * t464 * MDP(17) + (-t464 ^ 2 + t720 ^ 2) * MDP(18) + (t396 + t728) * MDP(19) + (-t397 - t727) * MDP(20) + t506 * MDP(21) + (-g(1) * t535 + g(2) * t533 + t386 * t549 + t457 * t720 + t596 * t707 + t620) * MDP(22) + (g(1) * t536 - g(2) * t534 + t385 * t549 + t457 * t464 + t598 * t707 - t630) * MDP(23) + (t371 + t724) * MDP(26) + (-t372 - t723) * MDP(27) + (-(-t375 * t610 - t699) * t543 - t368 * qJD(6) + (t410 * t720 + t502 * t614 - t543 * t669) * pkin(5) + t721) * MDP(29) + ((-t376 * t543 - t365) * t610 + (t375 * t543 - t647) * t614 + (-t502 * t610 - t543 * t668 - t636 * t720) * pkin(5) + t722) * MDP(30) + t719; (t661 + t724) * MDP(26) + (-t652 - t723) * MDP(27) + (t368 * t543 + t721) * MDP(29) + (-t610 * t365 - t614 * t366 + t367 * t543 + t722) * MDP(30) + (MDP(26) * t694 + MDP(27) * t636 - MDP(29) * t368 - MDP(30) * t700) * qJD(6) + t719;];
tau  = t1;
