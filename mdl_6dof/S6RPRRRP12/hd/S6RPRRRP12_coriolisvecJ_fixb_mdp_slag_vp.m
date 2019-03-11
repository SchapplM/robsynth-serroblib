% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:51:02
% EndTime: 2019-03-09 06:51:23
% DurationCPUTime: 13.64s
% Computational Cost: add. (19191->620), mult. (62951->830), div. (0->0), fcn. (53612->12), ass. (0->243)
t568 = sin(pkin(12));
t573 = sin(qJ(3));
t569 = sin(pkin(6));
t671 = qJD(1) * t569;
t655 = t573 * t671;
t570 = cos(pkin(12));
t711 = sin(pkin(7));
t713 = cos(pkin(6));
t630 = t713 * t711;
t716 = cos(qJ(3));
t605 = t716 * t630;
t712 = cos(pkin(7));
t638 = t712 * t716;
t610 = t638 * t671;
t675 = -qJD(1) * t605 - t570 * t610;
t512 = t568 * t655 + t675;
t604 = qJD(4) + t512;
t648 = qJD(1) * t713;
t640 = pkin(1) * t648;
t559 = t570 * t640;
t694 = t568 * t569;
t588 = t713 * pkin(2) + (-pkin(9) * t712 - qJ(2)) * t694;
t511 = qJD(1) * t588 + t559;
t535 = (-pkin(9) * t568 * t711 - pkin(2) * t570 - pkin(1)) * t569;
t526 = qJD(1) * t535 + qJD(2);
t486 = -t511 * t711 + t712 * t526;
t652 = t573 * t712;
t522 = t569 * (t568 * t716 + t570 * t652) + t573 * t630;
t516 = qJD(1) * t522;
t438 = t512 * pkin(3) - t516 * pkin(10) + t486;
t649 = t711 * t526;
t656 = qJ(2) * t671;
t539 = t568 * t640 + t570 * t656;
t692 = t569 * t570;
t587 = (t692 * t712 + t630) * pkin(9);
t507 = qJD(1) * t587 + t539;
t501 = t716 * t507;
t650 = t712 * t511;
t729 = t573 * t650 + t501;
t462 = t573 * t649 + t729;
t631 = t713 * t712;
t647 = qJD(1) * t711;
t636 = t569 * t647;
t552 = t570 * t636;
t661 = qJD(3) - t552;
t721 = qJD(1) * t631 + t661;
t442 = pkin(10) * t721 + t462;
t572 = sin(qJ(4));
t575 = cos(qJ(4));
t408 = t575 * t438 - t572 * t442;
t409 = t572 * t438 + t575 * t442;
t488 = t572 * t516 - t575 * t721;
t490 = t575 * t516 + t572 * t721;
t745 = t486 * MDP(13) + t490 * MDP(17) - t488 * MDP(18) + t408 * MDP(20) - t409 * MDP(21) - MDP(9) * t516;
t601 = t568 * t652 - t570 * t716;
t594 = t569 * t601;
t637 = t711 * t716;
t744 = qJD(1) * t594 + qJD(3) * t637;
t743 = t462 - t604 * (pkin(4) * t572 - pkin(11) * t575);
t515 = t522 * qJD(3);
t742 = qJD(1) * t515;
t668 = qJD(4) * t572;
t741 = pkin(10) * t668;
t651 = t573 * t711;
t542 = t572 * t651 - t575 * t712;
t620 = t568 * t636;
t739 = qJD(4) * t542 + t572 * t620 - t575 * t744;
t543 = t572 * t712 + t575 * t651;
t677 = qJD(4) * t543 + t572 * t744 + t575 * t620;
t599 = t568 * t638 + t570 * t573;
t593 = t569 * t599;
t635 = qJD(3) * t651;
t738 = qJD(1) * t593 - t635;
t596 = qJD(4) * t604;
t737 = t572 * t742 + t575 * t596;
t506 = t512 * qJD(3);
t689 = t575 * t506;
t581 = -qJD(4) * t488 - t689;
t736 = qJD(5) * t604 + t581;
t697 = t512 * t572;
t735 = t697 + t668;
t734 = MDP(4) * t568 + MDP(5) * t570;
t487 = qJD(5) + t488;
t731 = (t568 ^ 2 + t570 ^ 2) * MDP(6) * t569 ^ 2;
t621 = t570 * t638;
t693 = t568 * t573;
t658 = t569 * t693;
t521 = -t569 * t621 - t605 + t658;
t657 = pkin(1) * t713;
t562 = t570 * t657;
t523 = t562 + t588;
t491 = -t523 * t711 + t712 * t535;
t453 = t521 * pkin(3) - t522 * pkin(10) + t491;
t653 = t569 * t711;
t541 = t570 * t653 - t631;
t674 = qJ(2) * t692 + t568 * t657;
t519 = t587 + t674;
t583 = t716 * t519 + (t523 * t712 + t535 * t711) * t573;
t460 = -t541 * pkin(10) + t583;
t681 = t572 * t453 + t575 * t460;
t412 = pkin(11) * t521 + t681;
t727 = -t573 * t519 + t523 * t638 + t535 * t637;
t459 = t541 * pkin(3) - t727;
t493 = t522 * t572 + t541 * t575;
t494 = t522 * t575 - t541 * t572;
t424 = t493 * pkin(4) - t494 * pkin(11) + t459;
t571 = sin(qJ(5));
t574 = cos(qJ(5));
t730 = t574 * t412 + t571 * t424;
t728 = -t573 * t507 + t511 * t638 + t526 * t637;
t481 = pkin(3) * t516 + pkin(10) * t512;
t680 = t572 * t481 + t575 * t728;
t422 = pkin(11) * t516 + t680;
t555 = -pkin(4) * t575 - pkin(11) * t572 - pkin(3);
t664 = qJD(5) * t574;
t726 = t574 * t422 - t555 * t664 + t571 * t743;
t691 = t572 * t506;
t464 = qJD(4) * t490 - t691;
t725 = t464 * MDP(26);
t724 = t487 * MDP(26);
t403 = -pkin(4) * t604 - t408;
t465 = t490 * t571 - t574 * t604;
t467 = t574 * t490 + t571 * t604;
t387 = t465 * pkin(5) - t467 * qJ(6) + t403;
t714 = pkin(11) * t464;
t723 = t387 * t487 - t714;
t722 = t573 * (t650 + t649) + t501;
t590 = qJD(2) * t594;
t447 = qJD(3) * t727 - t590;
t514 = (t605 + (t621 - t693) * t569) * qJD(3);
t619 = t568 * qJD(2) * t653;
t477 = t515 * pkin(3) - t514 * pkin(10) + t619;
t666 = qJD(4) * t575;
t607 = t575 * t447 + t453 * t666 - t460 * t668 + t572 * t477;
t389 = pkin(11) * t515 + t607;
t589 = qJD(2) * t593;
t448 = qJD(3) * t583 + t589;
t471 = qJD(4) * t494 + t514 * t572;
t472 = -qJD(4) * t493 + t514 * t575;
t406 = t471 * pkin(4) - t472 * pkin(11) + t448;
t720 = -qJD(5) * t730 - t389 * t571 + t406 * t574;
t718 = t467 ^ 2;
t717 = t487 ^ 2;
t715 = pkin(5) * t464;
t709 = qJ(6) * t464;
t432 = -qJD(1) * t590 + qJD(3) * t728;
t609 = qJD(1) * t619;
t473 = pkin(3) * t742 + t506 * pkin(10) + t609;
t642 = t572 * t432 + t438 * t668 + t442 * t666 - t575 * t473;
t384 = -pkin(4) * t742 + t642;
t665 = qJD(5) * t571;
t417 = t490 * t665 - t571 * t742 - t574 * t736;
t418 = t490 * t664 + t571 * t736 - t574 * t742;
t375 = t418 * pkin(5) + t417 * qJ(6) - t467 * qJD(6) + t384;
t708 = t375 * t571;
t707 = t375 * t574;
t404 = pkin(11) * t604 + t409;
t441 = -pkin(3) * t721 - t728;
t416 = t488 * pkin(4) - t490 * pkin(11) + t441;
t380 = t404 * t574 + t416 * t571;
t378 = qJ(6) * t487 + t380;
t706 = t378 * t487;
t705 = t380 * t487;
t704 = t384 * t571;
t703 = t417 * t571;
t702 = t464 * t571;
t701 = t464 * t574;
t700 = t465 * t487;
t699 = t467 * t465;
t698 = t467 * t487;
t696 = t512 * t575;
t695 = t555 * t574;
t690 = t574 * t575;
t663 = qJD(5) * t575;
t667 = qJD(4) * t574;
t688 = qJD(6) * t575 - (-t571 * t663 - t572 * t667) * pkin(10) + t726 - t735 * qJ(6);
t687 = -t555 * t665 + (t422 + t741) * t571 + (-pkin(10) * t663 - t743) * t574 + t735 * pkin(5);
t628 = pkin(5) * t571 - qJ(6) * t574;
t686 = qJD(6) * t571 - t487 * t628 + t409;
t645 = t575 * t481 - t572 * t728;
t421 = -pkin(4) * t516 - t645;
t478 = -t574 * t516 - t571 * t696;
t479 = -t512 * t690 + t516 * t571;
t618 = pkin(10) + t628;
t629 = pkin(5) * t574 + qJ(6) * t571;
t685 = pkin(5) * t478 - qJ(6) * t479 + t421 - (qJD(5) * t629 - qJD(6) * t574) * t572 - t618 * t666;
t451 = pkin(4) * t490 + pkin(11) * t488;
t684 = t574 * t408 + t571 * t451;
t528 = t574 * t543 - t571 * t637;
t679 = -qJD(5) * t528 + t571 * t739 - t574 * t738;
t602 = -t571 * t543 - t574 * t637;
t678 = -qJD(5) * t602 + t571 * t738 + t574 * t739;
t673 = pkin(10) * t690 + t571 * t555;
t670 = qJD(2) * t569;
t669 = qJD(4) * t571;
t662 = t441 * qJD(4);
t379 = -t404 * t571 + t416 * t574;
t660 = qJD(6) - t379;
t659 = pkin(11) * t665;
t654 = t487 * t665;
t646 = t453 * t575 - t572 * t460;
t644 = t487 * t574;
t608 = -t575 * t432 - t438 * t666 + t442 * t668 - t572 * t473;
t383 = pkin(11) * t742 - t608;
t401 = t464 * pkin(4) - pkin(11) * t581 + t526 * t635 + t729 * qJD(3) + (t568 * t610 + t570 * t655) * qJD(2);
t643 = t571 * t383 - t574 * t401 + t404 * t664 + t416 * t665;
t633 = t571 * t666 - t478;
t632 = t574 * t666 - t479;
t377 = -pkin(5) * t487 + t660;
t627 = t377 * t574 - t378 * t571;
t624 = -t412 * t571 + t424 * t574;
t475 = t494 * t574 + t521 * t571;
t474 = t494 * t571 - t521 * t574;
t622 = (-t568 * t656 + t559) * t568 - t539 * t570;
t617 = -t572 * t447 - t453 * t668 - t460 * t666 + t477 * t575;
t411 = -pkin(4) * t521 - t646;
t613 = -t487 * t664 - t702;
t612 = t403 * t487 - t714;
t611 = t387 * t467 + t643;
t371 = t574 * t383 + t571 * t401 - t404 * t665 + t416 * t664;
t606 = t574 * t389 + t571 * t406 - t412 * t665 + t424 * t664;
t390 = -pkin(4) * t515 - t617;
t597 = t512 * t604;
t554 = -pkin(4) - t629;
t537 = t618 * t572;
t531 = -t695 + (pkin(10) * t571 + pkin(5)) * t575;
t530 = -qJ(6) * t575 + t673;
t433 = qJD(1) * t589 + t722 * qJD(3);
t427 = pkin(5) * t467 + qJ(6) * t465;
t426 = -qJD(5) * t474 + t472 * t574 + t515 * t571;
t425 = qJD(5) * t475 + t472 * t571 - t515 * t574;
t396 = -t417 + t700;
t395 = pkin(5) * t474 - qJ(6) * t475 + t411;
t392 = -pkin(5) * t490 + t408 * t571 - t451 * t574;
t391 = qJ(6) * t490 + t684;
t386 = -pkin(5) * t493 - t624;
t385 = qJ(6) * t493 + t730;
t376 = pkin(5) * t425 - qJ(6) * t426 - qJD(6) * t475 + t390;
t374 = -pkin(5) * t471 - t720;
t373 = qJ(6) * t471 + qJD(6) * t493 + t606;
t370 = t643 - t715;
t369 = qJD(6) * t487 + t371 + t709;
t1 = [(t379 * t471 + t384 * t474 + t390 * t465 + t403 * t425 + t411 * t418 + t624 * t464 + t487 * t720 - t493 * t643) * MDP(27) + (-t371 * t493 - t380 * t471 + t384 * t475 + t390 * t467 + t403 * t426 - t411 * t417 - t464 * t730 - t487 * t606) * MDP(28) + (-t506 * t522 + t514 * t516) * MDP(8) + (t506 * t541 + t514 * t721) * MDP(10) + (t432 * t541 - t447 * t721 + t486 * t514 - t491 * t506) * MDP(14) + 0.2e1 * qJD(2) * qJD(1) * t731 + (t472 * t604 + t494 * t742 + t521 * t581) * MDP(17) + ((t516 * t711 + t522 * t647) * t568 * MDP(14) + ((t570 * t674 + (qJ(2) * t694 - t562) * t568) * qJD(1) - t622) * MDP(7) - 0.2e1 * t734 * t648) * t670 + (-t418 * t493 - t425 * t487 - t464 * t474 - t465 * t471) * MDP(25) + (-t417 * t493 + t426 * t487 + t464 * t475 + t467 * t471) * MDP(24) + (t369 * t493 + t373 * t487 - t375 * t475 - t376 * t467 + t378 * t471 + t385 * t464 - t387 * t426 + t395 * t417) * MDP(31) + (t464 * t493 + t471 * t487) * MDP(26) + (-t370 * t493 - t374 * t487 + t375 * t474 + t376 * t465 - t377 * t471 - t386 * t464 + t387 * t425 + t395 * t418) * MDP(29) + (t433 * t541 - t448 * t721 + t491 * t742 + t512 * t619 + t521 * t609) * MDP(13) + (t433 * t493 + t441 * t471 + t448 * t488 + t459 * t464 - t521 * t642 + t604 * t617 + t646 * t742) * MDP(20) + (t433 * t494 + t441 * t472 + t448 * t490 + t459 * t581 + t521 * t608 - t604 * t607 - t681 * t742) * MDP(21) + (t506 * t521 - t514 * t512 - t522 * t742) * MDP(9) + ((qJD(4) + t675 + (t521 + t658) * qJD(1)) * MDP(19) + (-t661 + (t541 - t631) * qJD(1)) * MDP(11) + t745) * t515 + (-t369 * t474 + t370 * t475 - t373 * t465 + t374 * t467 + t377 * t426 - t378 * t425 - t385 * t418 - t386 * t417) * MDP(30) + (t417 * t474 - t418 * t475 - t425 * t467 - t426 * t465) * MDP(23) + (-t417 * t475 + t426 * t467) * MDP(22) + (-t464 * t521 - t471 * t604 - t493 * t742) * MDP(18) + (t490 * t472 + t494 * t581) * MDP(15) + (-t494 * t464 - t490 * t471 - t472 * t488 - t493 * t581) * MDP(16) + (t369 * t385 + t370 * t386 + t373 * t378 + t374 * t377 + t375 * t395 + t376 * t387) * MDP(32); t622 * MDP(7) * t671 + (-t512 * t620 + t712 * t742) * MDP(13) + (-t506 * t712 - t516 * t620) * MDP(14) + (-t464 * t637 - t488 * t738 - t542 * t742) * MDP(20) + (-t490 * t738 - t543 * t742 - t581 * t637) * MDP(21) + (t417 * t602 - t418 * t528 + t465 * t678 - t467 * t679) * MDP(30) + (t369 * t528 - t370 * t602 + t375 * t542 - t377 * t679 - t378 * t678 + t387 * t677) * MDP(32) + (-t677 * MDP(20) + MDP(21) * t739) * t604 - (-t738 * MDP(13) + MDP(14) * t744) * t721 + (MDP(27) + MDP(29)) * (t542 * t418 + t464 * t602 + t465 * t677 + t487 * t679) + (MDP(28) - MDP(31)) * (-t417 * t542 - t464 * t528 + t467 * t677 + t487 * t678) + (t569 * t713 * t734 - t731) * qJD(1) ^ 2; -t512 ^ 2 * MDP(9) + (t512 * t721 - t506) * MDP(10) - MDP(11) * t742 + (-t462 * t552 + (t462 * t631 - t599 * t670) * qJD(1) + (t462 - t722) * qJD(3)) * MDP(13) + (-t728 * t552 + t486 * t512 + (t601 * t670 + t631 * t728) * qJD(1)) * MDP(14) + (-t735 * t490 + (-t666 - t696) * t488) * MDP(16) + t737 * MDP(17) + (-pkin(3) * t464 - pkin(10) * t737 + t441 * t697 - t462 * t488 - t604 * t645) * MDP(20) + (-pkin(3) * t581 + t441 * t696 - t462 * t490 + (t680 + t741) * t604) * MDP(21) + t632 * t467 * MDP(22) + (t465 * t479 + t467 * t478 + (-t465 * t574 - t467 * t571) * t666) * MDP(23) + (-t403 * t478 - t421 * t465 + t464 * t695) * MDP(27) + (-t403 * t479 - t421 * t467 - t673 * t464) * MDP(28) + (t633 * t387 + t418 * t537 - t464 * t531 - t685 * t465) * MDP(29) + (-t377 * t479 + t378 * t478 - t417 * t531 - t418 * t530 + t688 * t465 - t687 * t467 + t627 * t666) * MDP(30) + (-t632 * t387 + t417 * t537 + t464 * t530 + t685 * t467) * MDP(31) + (t369 * t530 + t370 * t531 + t375 * t537 - t377 * t687 - t378 * t688 - t387 * t685) * MDP(32) + (t632 * MDP(24) - t633 * MDP(25) + (-t743 * t574 + (-qJD(5) * t555 + t422) * t571) * MDP(27) + t726 * MDP(28) + t687 * MDP(29) - t688 * MDP(31)) * t487 + (-t464 * MDP(16) + (-t596 - t597) * MDP(18) + t662 * MDP(20) + t433 * MDP(21) + (-t417 * t574 - t467 * t665) * MDP(22) + (t703 - t418 * t574 + (t465 * t571 - t467 * t574) * qJD(5)) * MDP(23) + (t467 * t604 - t654 + t701) * MDP(24) + (-t465 * t604 + t613) * MDP(25) + t604 * t724 + (t403 * t664 + t704 + t604 * t379 + (t487 * t669 + t418) * pkin(10)) * MDP(27) + (-t403 * t665 + t384 * t574 - t604 * t380 + (t487 * t667 - t417) * pkin(10)) * MDP(28) + (-t377 * t604 + t387 * t664 + t708) * MDP(29) + (-t369 * t571 + t370 * t574 + (-t377 * t571 - t378 * t574) * qJD(5)) * MDP(30) + (t378 * t604 + t387 * t665 - t707) * MDP(31)) * t572 + (((qJD(4) * t721 - t506) * t572 + t604 * t490) * MDP(15) + t581 * MDP(16) + t597 * MDP(17) + t742 * MDP(18) - t433 * MDP(20) + (-pkin(10) * t742 + t662) * MDP(21) + t417 * MDP(24) + t418 * MDP(25) - t725 + (t403 * t669 + t643 + (qJD(4) * t465 + t613) * pkin(10)) * MDP(27) + (t403 * t667 + t371 + (qJD(4) * t467 + t654) * pkin(10)) * MDP(28) + t370 * MDP(29) - t369 * MDP(31)) * t575 + (-qJD(4) * t572 ^ 2 * MDP(15) + MDP(11) * t721 - t604 * MDP(19) + t512 * MDP(8) - t745) * t516; -t488 ^ 2 * MDP(16) + (t512 * t488 - t689) * MDP(17) + t691 * MDP(18) + MDP(19) * t742 + (t409 * t604 - t642) * MDP(20) + (t408 * t604 + t441 * t488 + t608) * MDP(21) + (t467 * t644 - t703) * MDP(22) + ((-t417 - t700) * t574 + (-t418 - t698) * t571) * MDP(23) + (t487 * t644 + t702) * MDP(24) + (-t571 * t717 + t701) * MDP(25) + (-pkin(4) * t418 - t409 * t465 + (-t384 + (-pkin(11) * qJD(5) - t451) * t487) * t574 + (t408 * t487 + t612) * t571) * MDP(27) + (pkin(4) * t417 + t704 - t409 * t467 + (t659 + t684) * t487 + t612 * t574) * MDP(28) + (-t707 + t418 * t554 + (-pkin(11) * t664 + t392) * t487 - t686 * t465 + t723 * t571) * MDP(29) + (t391 * t465 - t392 * t467 + (t369 + t487 * t377 + (qJD(5) * t467 - t418) * pkin(11)) * t574 + (t370 - t706 + (qJD(5) * t465 - t417) * pkin(11)) * t571) * MDP(30) + (-t708 + t417 * t554 + (-t391 - t659) * t487 + t686 * t467 - t723 * t574) * MDP(31) + (t375 * t554 - t377 * t392 - t378 * t391 - t686 * t387 + (qJD(5) * t627 + t369 * t574 + t370 * t571) * pkin(11)) * MDP(32) + (MDP(15) * t488 + MDP(16) * t490 + MDP(18) * t512 - t441 * MDP(20) - t467 * MDP(24) + t465 * MDP(25) - MDP(27) * t379 + MDP(28) * t380 + t377 * MDP(29) - t378 * MDP(31) - t724) * t490; MDP(22) * t699 + (-t465 ^ 2 + t718) * MDP(23) + t396 * MDP(24) + (-t418 + t698) * MDP(25) + t725 + (-t403 * t467 - t643 + t705) * MDP(27) + (t379 * t487 + t403 * t465 - t371) * MDP(28) + (-t427 * t465 - t611 + t705 + 0.2e1 * t715) * MDP(29) + (pkin(5) * t417 - qJ(6) * t418 + (t378 - t380) * t467 + (t377 - t660) * t465) * MDP(30) + (0.2e1 * t709 - t387 * t465 + t427 * t467 + (0.2e1 * qJD(6) - t379) * t487 + t371) * MDP(31) + (-pkin(5) * t370 + qJ(6) * t369 - t377 * t380 + t378 * t660 - t387 * t427) * MDP(32); (-t464 + t699) * MDP(29) + t396 * MDP(30) + (-t717 - t718) * MDP(31) + (t611 - t706 - t715) * MDP(32);];
tauc  = t1;
