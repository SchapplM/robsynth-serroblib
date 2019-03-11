% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:51
% EndTime: 2019-03-09 05:32:07
% DurationCPUTime: 8.88s
% Computational Cost: add. (17004->528), mult. (56228->748), div. (0->0), fcn. (48730->14), ass. (0->235)
t609 = sin(pkin(7));
t608 = sin(pkin(12));
t610 = sin(pkin(6));
t612 = cos(pkin(12));
t619 = cos(qJ(3));
t613 = cos(pkin(7));
t616 = sin(qJ(3));
t712 = t613 * t616;
t628 = (t608 * t712 - t612 * t619) * t610;
t696 = qJD(3) * t619;
t753 = qJD(1) * t628 + t609 * t696;
t730 = cos(pkin(6));
t676 = t730 * t609;
t655 = t619 * t676;
t711 = t613 * t619;
t688 = t612 * t711;
t752 = t610 * t688 + t655;
t675 = qJD(1) * t730;
t658 = pkin(1) * t675;
t699 = qJD(1) * t610;
t683 = t612 * t699;
t572 = qJ(2) * t683 + t608 * t658;
t714 = t610 * t612;
t626 = (t613 * t714 + t676) * pkin(9);
t532 = qJD(1) * t626 + t572;
t595 = t612 * t658;
t718 = t608 * t610;
t623 = t730 * pkin(2) + (-pkin(9) * t613 - qJ(2)) * t718;
t539 = qJD(1) * t623 + t595;
t566 = (-pkin(9) * t608 * t609 - pkin(2) * t612 - pkin(1)) * t610;
t558 = qJD(1) * t566 + qJD(2);
t737 = t616 * t532 - (t539 * t613 + t558 * t609) * t619;
t751 = qJD(3) * t737;
t684 = t608 * t699;
t741 = qJD(1) * t752 - t616 * t684;
t538 = qJD(4) - t741;
t553 = t610 * (t608 * t619 + t612 * t712) + t616 * t676;
t545 = qJD(1) * t553;
t498 = pkin(3) * t545 - pkin(10) * t741;
t618 = cos(qJ(4));
t497 = t618 * t498;
t615 = sin(qJ(4));
t731 = -qJ(5) - pkin(10);
t677 = qJD(4) * t731;
t750 = -pkin(4) * t545 - t497 + (qJ(5) * t741 + t677) * t618 + (-qJD(5) - t737) * t615;
t708 = t615 * t498 - t618 * t737;
t720 = t741 * t615;
t749 = -qJ(5) * t720 - qJD(5) * t618 - t615 * t677 + t708;
t617 = cos(qJ(6));
t614 = sin(qJ(6));
t567 = t609 * t683 - t613 * t675 - qJD(3);
t504 = t545 * t615 + t618 * t567;
t506 = t545 * t618 - t567 * t615;
t607 = sin(pkin(13));
t611 = cos(pkin(13));
t646 = -t504 * t607 + t611 * t506;
t724 = t646 * t614;
t445 = -t617 * t538 + t724;
t667 = -t611 * t504 - t506 * t607;
t735 = qJD(6) - t667;
t748 = t445 * t735;
t447 = t538 * t614 + t617 * t646;
t747 = t447 * t735;
t583 = t607 * t615 - t611 * t618;
t746 = t538 * t583;
t716 = t609 * t616;
t578 = t613 * t615 + t618 * t716;
t661 = t609 * t684;
t745 = qJD(4) * t578 + t615 * t753 + t618 * t661;
t577 = t613 * t618 - t615 * t716;
t744 = -qJD(4) * t577 + t615 * t661 - t618 * t753;
t584 = t607 * t618 + t611 * t615;
t702 = t538 * t584;
t698 = qJD(2) * t610;
t682 = t608 * t698;
t713 = t612 * t616;
t743 = t682 * t711 + t698 * t713;
t665 = t735 * t617;
t527 = t741 * qJD(3);
t694 = qJD(4) * t618;
t695 = qJD(4) * t615;
t481 = t618 * t527 - t545 * t695 - t567 * t694;
t482 = qJD(4) * t506 + t527 * t615;
t440 = t481 * t607 + t611 * t482;
t710 = t614 * t440;
t742 = -t735 * t665 - t710;
t738 = (t608 ^ 2 + t612 ^ 2) * MDP(6) * t610 ^ 2;
t685 = pkin(1) * t730;
t599 = t612 * t685;
t554 = t599 + t623;
t507 = -t554 * t609 + t613 * t566;
t717 = t608 * t616;
t552 = t610 * t717 - t752;
t469 = pkin(3) * t552 - pkin(10) * t553 + t507;
t701 = qJ(2) * t714 + t608 * t685;
t549 = t626 + t701;
t574 = t609 * t714 - t613 * t730;
t644 = t554 * t613 + t566 * t609;
t475 = -pkin(10) * t574 + t549 * t619 + t616 * t644;
t709 = t615 * t469 + t618 * t475;
t705 = t749 * t607 + t611 * t750;
t703 = t607 * t750 - t749 * t611;
t477 = t619 * t532 + t539 * t712 + t558 * t716;
t736 = -t477 + (t695 - t720) * pkin(4);
t544 = t553 * qJD(3);
t528 = qJD(1) * t544;
t499 = -t539 * t609 + t613 * t558;
t455 = -pkin(3) * t741 - pkin(10) * t545 + t499;
t457 = -pkin(10) * t567 + t477;
t431 = t455 * t615 + t457 * t618;
t625 = qJD(2) * t628;
t451 = -qJD(1) * t625 - t751;
t659 = t609 * t682;
t488 = pkin(3) * t528 - pkin(10) * t527 + qJD(1) * t659;
t672 = -t451 * t615 + t618 * t488;
t622 = -qJD(4) * t431 + t672;
t392 = pkin(4) * t528 - qJ(5) * t481 - qJD(5) * t506 + t622;
t634 = t618 * t451 + t455 * t694 - t457 * t695 + t615 * t488;
t396 = -qJ(5) * t482 - qJD(5) * t504 + t634;
t379 = t392 * t611 - t396 * t607;
t377 = -pkin(5) * t528 - t379;
t602 = pkin(4) * t607 + pkin(11);
t734 = t735 * (pkin(4) * t506 + pkin(5) * t646 - pkin(11) * t667 + qJD(6) * t602) + t377;
t733 = -t616 * t549 + t619 * t644;
t732 = MDP(4) * t608 + MDP(5) * t612;
t441 = t481 * t611 - t482 * t607;
t692 = qJD(6) * t617;
t687 = t617 * t441 + t614 * t528 + t538 * t692;
t693 = qJD(6) * t614;
t413 = -t646 * t693 + t687;
t728 = t413 * t614;
t423 = -qJ(5) * t504 + t431;
t727 = t423 * t607;
t726 = t445 * t646;
t725 = t447 * t646;
t723 = t504 * t538;
t722 = t506 * t538;
t719 = t584 * t617;
t715 = t609 * t619;
t420 = t611 * t423;
t438 = t617 * t440;
t380 = t607 * t392 + t611 * t396;
t509 = t553 * t615 + t574 * t618;
t543 = (t655 + (t688 - t717) * t610) * qJD(3);
t487 = -qJD(4) * t509 + t543 * t618;
t510 = t553 * t618 - t574 * t615;
t464 = qJD(3) * t733 - t625;
t493 = pkin(3) * t544 - pkin(10) * t543 + t659;
t671 = -t464 * t615 + t618 * t493;
t402 = pkin(4) * t544 - qJ(5) * t487 - qJD(4) * t709 - qJD(5) * t510 + t671;
t486 = qJD(4) * t510 + t543 * t615;
t633 = t618 * t464 + t469 * t694 - t475 * t695 + t615 * t493;
t406 = -qJ(5) * t486 - qJD(5) * t509 + t633;
t384 = t607 * t402 + t611 * t406;
t430 = t618 * t455 - t457 * t615;
t422 = -qJ(5) * t506 + t430;
t419 = pkin(4) * t538 + t422;
t391 = t607 * t419 + t420;
t670 = t618 * t469 - t475 * t615;
t425 = pkin(4) * t552 - qJ(5) * t510 + t670;
t429 = -qJ(5) * t509 + t709;
t404 = t607 * t425 + t611 * t429;
t707 = -t607 * t745 - t611 * t744;
t706 = -t607 * t744 + t611 * t745;
t704 = pkin(5) * t545 - t705;
t697 = qJD(3) * t616;
t690 = qJD(1) * qJD(2);
t686 = -pkin(4) * t618 - pkin(3);
t681 = t609 * t697;
t679 = t613 * t697;
t678 = t615 * t731;
t378 = pkin(11) * t528 + t380;
t452 = qJD(1) * t743 + t532 * t696 + t539 * t679 + t558 * t681;
t435 = pkin(4) * t482 + t452;
t398 = pkin(5) * t440 - pkin(11) * t441 + t435;
t674 = -t378 * t614 + t617 * t398;
t673 = t441 * t614 - t617 * t528;
t669 = -t617 * t545 + t614 * t746;
t668 = t545 * t614 + t617 * t746;
t666 = t538 * t618;
t654 = -(t608 * t711 + t713) * t699 + t681;
t546 = pkin(5) * t583 - pkin(11) * t584 + t686;
t653 = pkin(11) * t545 - qJD(6) * t546 - t703;
t591 = t731 * t618;
t560 = -t611 * t591 + t607 * t678;
t652 = -pkin(5) * t702 - pkin(11) * t746 + qJD(6) * t560 - t736;
t465 = t549 * t696 + t554 * t679 + t566 * t681 + t743;
t651 = t378 * t617 + t398 * t614;
t388 = pkin(11) * t538 + t391;
t456 = pkin(3) * t567 + t737;
t442 = pkin(4) * t504 + qJD(5) + t456;
t412 = -pkin(5) * t667 - pkin(11) * t646 + t442;
t386 = t388 * t617 + t412 * t614;
t650 = t388 * t614 - t412 * t617;
t400 = pkin(11) * t552 + t404;
t478 = t611 * t509 + t510 * t607;
t479 = -t509 * t607 + t510 * t611;
t474 = pkin(3) * t574 - t733;
t621 = pkin(4) * t509 + t474;
t415 = pkin(5) * t478 - pkin(11) * t479 + t621;
t649 = t400 * t617 + t415 * t614;
t648 = -t400 * t614 + t415 * t617;
t383 = t402 * t611 - t406 * t607;
t390 = t419 * t611 - t727;
t403 = t425 * t611 - t429 * t607;
t449 = t479 * t617 + t552 * t614;
t448 = t479 * t614 - t552 * t617;
t643 = (-qJ(2) * t684 + t595) * t608 - t572 * t612;
t642 = t438 + (t614 * t667 - t693) * t735;
t641 = qJD(6) * t715 - t707;
t637 = -pkin(10) * t528 + t456 * t538;
t632 = t584 * t692 - t669;
t631 = -t584 * t693 - t668;
t630 = pkin(4) * t486 + t465;
t530 = t577 * t607 + t578 * t611;
t629 = -qJD(6) * t530 + t654;
t387 = -pkin(5) * t538 - t390;
t395 = t422 * t611 - t727;
t624 = -t602 * t440 + (t387 + t395) * t735;
t603 = -pkin(4) * t611 - pkin(5);
t559 = -t591 * t607 - t611 * t678;
t529 = -t611 * t577 + t578 * t607;
t444 = -t486 * t607 + t487 * t611;
t443 = t611 * t486 + t487 * t607;
t418 = qJD(6) * t449 + t444 * t614 - t544 * t617;
t417 = -qJD(6) * t448 + t444 * t617 + t544 * t614;
t414 = qJD(6) * t447 + t673;
t407 = pkin(5) * t443 - pkin(11) * t444 + t630;
t399 = -pkin(5) * t552 - t403;
t394 = t422 * t607 + t420;
t382 = pkin(11) * t544 + t384;
t381 = -pkin(5) * t544 - t383;
t376 = -qJD(6) * t386 + t674;
t375 = -qJD(6) * t650 + t651;
t1 = [(t379 * t403 + t380 * t404 + t390 * t383 + t391 * t384 + t435 * t621 + t442 * t630) * MDP(23) + (t451 * t574 + t464 * t567 + t499 * t543 + t507 * t527 + 0.2e1 * t545 * t659) * MDP(14) + (t671 * t538 + t670 * t528 + t672 * t552 + t430 * t544 + t465 * t504 + t474 * t482 + t452 * t509 + t456 * t486 + (-t431 * t552 - t538 * t709) * qJD(4)) * MDP(20) + (-t431 * t544 + t452 * t510 + t456 * t487 + t465 * t506 + t474 * t481 - t528 * t709 - t538 * t633 - t552 * t634) * MDP(21) + (-t527 * t574 - t543 * t567) * MDP(10) + (t528 * t574 + t544 * t567) * MDP(11) + (t527 * t553 + t543 * t545) * MDP(8) + (-t527 * t552 - t528 * t553 + t543 * t741 - t544 * t545) * MDP(9) + (-t482 * t552 - t486 * t538 - t504 * t544 - t509 * t528) * MDP(18) + (t481 * t552 + t487 * t538 + t506 * t544 + t510 * t528) * MDP(17) + (t528 * t552 + t538 * t544) * MDP(19) + (-t481 * t509 - t482 * t510 - t486 * t506 - t487 * t504) * MDP(16) + (t481 * t510 + t487 * t506) * MDP(15) + (-(qJD(6) * t648 + t382 * t617 + t407 * t614) * t735 - t649 * t440 - t375 * t478 - t386 * t443 + t381 * t447 + t399 * t413 + t377 * t449 + t387 * t417) * MDP(30) + 0.2e1 * t690 * t738 + (t413 * t449 + t417 * t447) * MDP(24) + (-t413 * t448 - t414 * t449 - t417 * t445 - t418 * t447) * MDP(25) + (-t379 * t479 - t380 * t478 - t383 * t646 + t384 * t667 - t390 * t444 - t391 * t443 - t403 * t441 - t404 * t440) * MDP(22) + ((-qJD(6) * t649 - t382 * t614 + t407 * t617) * t735 + t648 * t440 + t376 * t478 - t650 * t443 + t381 * t445 + t399 * t414 + t377 * t448 + t387 * t418) * MDP(29) + (t440 * t478 + t443 * t735) * MDP(28) + (-t414 * t478 - t418 * t735 - t440 * t448 - t443 * t445) * MDP(27) + (t413 * t478 + t417 * t735 + t440 * t449 + t443 * t447) * MDP(26) + (t452 * t574 + t465 * t567 + t499 * t544 + t507 * t528 + (qJD(1) * t552 - t741) * t659) * MDP(13) + (((t612 * t701 + (qJ(2) * t718 - t599) * t608) * qJD(1) - t643) * MDP(7) - 0.2e1 * t732 * t675) * t698; t643 * MDP(7) * t699 + (t528 * t613 + t567 * t654 + t661 * t741) * MDP(13) + (t527 * t613 - t545 * t661 + t567 * t753) * MDP(14) + (-t482 * t715 + t654 * t504 + t528 * t577 - t538 * t745) * MDP(20) + (-t481 * t715 + t654 * t506 - t528 * t578 + t538 * t744) * MDP(21) + (-t440 * t530 + t441 * t529 + t646 * t706 + t667 * t707) * MDP(22) + (-t379 * t529 + t380 * t530 - t390 * t706 + t391 * t707 - t435 * t715 + t442 * t654) * MDP(23) + ((-t530 * t614 - t617 * t715) * t440 + t529 * t414 + (t614 * t641 + t617 * t629) * t735 + t706 * t445) * MDP(29) + (-(t530 * t617 - t614 * t715) * t440 + t529 * t413 + (-t614 * t629 + t617 * t641) * t735 + t706 * t447) * MDP(30) + (t610 * t730 * t732 - t738) * qJD(1) ^ 2; -t741 ^ 2 * MDP(9) + (t567 * t741 + t527) * MDP(10) - MDP(11) * t528 + (-t477 * t567 - t452) * MDP(13) + (-t499 * t741 + t567 * t737 + t628 * t690 + t751) * MDP(14) + (t481 * t615 + t506 * t666) * MDP(15) + ((t481 - t723) * t618 + (-t482 - t722) * t615) * MDP(16) + (t615 * t528 + t538 * t666) * MDP(17) + (-t538 ^ 2 * t615 + t528 * t618) * MDP(18) + (-pkin(3) * t482 - t452 * t618 - t477 * t504 + (-pkin(10) * t694 - t497) * t538 + (-t538 * t737 + t637) * t615) * MDP(20) + (-pkin(3) * t481 + t452 * t615 - t477 * t506 + (pkin(10) * t695 + t708) * t538 + t637 * t618) * MDP(21) + (-t379 * t584 - t380 * t583 + t390 * t746 - t702 * t391 - t560 * t440 + t441 * t559 - t705 * t646 + t703 * t667) * MDP(22) + (-t379 * t559 + t380 * t560 + t705 * t390 + t703 * t391 + t435 * t686 + t442 * t736) * MDP(23) + (t413 * t719 + t447 * t631) * MDP(24) + (t669 * t447 + t668 * t445 + (-t728 - t414 * t617 + (t445 * t614 - t447 * t617) * qJD(6)) * t584) * MDP(25) + (t413 * t583 + t438 * t584 + t447 * t702 + t631 * t735) * MDP(26) + (-t414 * t583 - t445 * t702 - t584 * t710 - t632 * t735) * MDP(27) + (t440 * t583 + t702 * t735) * MDP(28) + ((t546 * t617 - t560 * t614) * t440 + t376 * t583 + t559 * t414 + t377 * t614 * t584 + (t614 * t653 - t617 * t652) * t735 + t704 * t445 - t702 * t650 + t632 * t387) * MDP(29) + (-(t546 * t614 + t560 * t617) * t440 - t375 * t583 + t559 * t413 + t377 * t719 + (t614 * t652 + t617 * t653) * t735 + t704 * t447 - t702 * t386 + t631 * t387) * MDP(30) + (-MDP(11) * t567 - MDP(13) * t499 - MDP(17) * t506 + MDP(18) * t504 - MDP(19) * t538 - MDP(20) * t430 + MDP(21) * t431 - MDP(8) * t741 + MDP(9) * t545) * t545; t506 * t504 * MDP(15) + (-t504 ^ 2 + t506 ^ 2) * MDP(16) + (t481 + t723) * MDP(17) + (-t482 + t722) * MDP(18) + t528 * MDP(19) + (t431 * t538 - t456 * t506 + t622) * MDP(20) + (t430 * t538 + t456 * t504 - t634) * MDP(21) + ((-t440 * t607 - t441 * t611) * pkin(4) + (t390 - t395) * t667 + (t391 - t394) * t646) * MDP(22) + (t390 * t394 - t391 * t395 + (t379 * t611 + t380 * t607 - t442 * t506) * pkin(4)) * MDP(23) + (t447 * t665 + t728) * MDP(24) + ((t413 - t748) * t617 + (-t414 - t747) * t614) * MDP(25) + (-t725 - t742) * MDP(26) + (t642 + t726) * MDP(27) - t735 * t646 * MDP(28) + (-t394 * t445 + t603 * t414 + t624 * t614 - t617 * t734 + t646 * t650) * MDP(29) + (t386 * t646 - t394 * t447 + t603 * t413 + t614 * t734 + t624 * t617) * MDP(30); (-t646 ^ 2 - t667 ^ 2) * MDP(22) + (t390 * t646 - t391 * t667 + t435) * MDP(23) + (t642 - t726) * MDP(29) + (-t725 + t742) * MDP(30); t447 * t445 * MDP(24) + (-t445 ^ 2 + t447 ^ 2) * MDP(25) + (t687 + t748) * MDP(26) + (-t673 + t747) * MDP(27) + t440 * MDP(28) + (t386 * t735 - t387 * t447 + t674) * MDP(29) + (t387 * t445 - t650 * t735 - t651) * MDP(30) + (-MDP(26) * t724 - MDP(27) * t447 - MDP(29) * t386 + MDP(30) * t650) * qJD(6);];
tauc  = t1;
