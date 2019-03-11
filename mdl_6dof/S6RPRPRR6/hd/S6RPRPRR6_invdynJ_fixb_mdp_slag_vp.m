% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:53:11
% EndTime: 2019-03-09 03:53:25
% DurationCPUTime: 11.39s
% Computational Cost: add. (8564->567), mult. (20596->725), div. (0->0), fcn. (16880->18), ass. (0->233)
t617 = cos(pkin(10));
t729 = cos(qJ(3));
t676 = t729 * t617;
t594 = qJD(1) * t676;
t615 = sin(pkin(10));
t621 = sin(qJ(3));
t702 = t621 * t615;
t675 = qJD(1) * t702;
t559 = -t594 + t675;
t553 = qJD(5) + t559;
t576 = t729 * t615 + t621 * t617;
t561 = t576 * qJD(1);
t614 = sin(pkin(11));
t616 = cos(pkin(11));
t529 = -t616 * qJD(3) + t561 * t614;
t531 = qJD(3) * t614 + t561 * t616;
t620 = sin(qJ(5));
t624 = cos(qJ(5));
t474 = t624 * t529 + t531 * t620;
t623 = cos(qJ(6));
t619 = sin(qJ(6));
t741 = t529 * t620 - t531 * t624;
t713 = t741 * t619;
t424 = t623 * t474 - t713;
t549 = qJD(6) + t553;
t745 = t424 * t549;
t749 = t474 * t553;
t650 = t474 * t619 + t623 * t741;
t744 = t549 * t650;
t573 = t614 * t620 - t624 * t616;
t696 = t553 * t573;
t575 = t614 * t624 + t616 * t620;
t563 = t575 * qJD(5);
t695 = t575 * t559 + t563;
t748 = t553 * t741;
t514 = pkin(3) * t561 + qJ(4) * t559;
t724 = pkin(7) + qJ(2);
t584 = t724 * t615;
t577 = qJD(1) * t584;
t586 = t724 * t617;
t578 = qJD(1) * t586;
t643 = t729 * t577 + t621 * t578;
t461 = t614 * t514 - t616 * t643;
t708 = t559 * t614;
t443 = pkin(8) * t708 + t461;
t746 = -qJD(4) * t616 + t443;
t683 = qJD(1) * qJD(2);
t731 = qJDD(1) * t724 + t683;
t546 = t731 * t615;
t547 = t731 * t617;
t674 = qJD(3) * t729;
t690 = qJD(3) * t621;
t640 = t546 * t729 + t621 * t547 - t577 * t690 + t578 * t674;
t458 = -qJDD(3) * pkin(3) + qJDD(4) + t640;
t613 = pkin(10) + qJ(3);
t605 = sin(t613);
t607 = cos(t613);
t622 = sin(qJ(1));
t625 = cos(qJ(1));
t660 = g(1) * t625 + g(2) * t622;
t636 = -g(3) * t607 + t605 * t660;
t630 = t458 - t636;
t521 = -t573 * t619 + t575 * t623;
t697 = qJD(6) * t521 - t696 * t619 + t695 * t623;
t720 = qJDD(1) * pkin(1);
t738 = g(1) * t622 - g(2) * t625;
t647 = -qJDD(2) + t720 + t738;
t565 = t576 * qJD(3);
t671 = qJDD(1) * t729;
t682 = qJDD(1) * t621;
t654 = t615 * t682 - t617 * t671;
t518 = qJD(1) * t565 + t654;
t599 = pkin(2) * t617 + pkin(1);
t582 = -qJD(1) * t599 + qJD(2);
t491 = pkin(3) * t559 - qJ(4) * t561 + t582;
t523 = -t621 * t577 + t729 * t578;
t513 = qJD(3) * qJ(4) + t523;
t449 = t616 * t491 - t513 * t614;
t418 = pkin(4) * t559 - pkin(8) * t531 + t449;
t450 = t614 * t491 + t616 * t513;
t430 = -pkin(8) * t529 + t450;
t400 = t418 * t620 + t430 * t624;
t390 = -pkin(9) * t474 + t400;
t686 = qJD(6) * t619;
t388 = t390 * t686;
t511 = -qJD(3) * pkin(3) + qJD(4) + t643;
t468 = t529 * pkin(4) + t511;
t422 = t474 * pkin(5) + t468;
t612 = pkin(11) + qJ(5);
t608 = qJ(6) + t612;
t597 = sin(t608);
t598 = cos(t608);
t704 = t607 * t622;
t534 = t597 * t625 - t598 * t704;
t703 = t607 * t625;
t536 = t597 * t622 + t598 * t703;
t726 = g(3) * t605;
t743 = g(1) * t536 - g(2) * t534 + t422 * t424 + t598 * t726 + t388;
t533 = t597 * t704 + t598 * t625;
t535 = -t597 * t703 + t598 * t622;
t680 = qJD(3) * t594 + t615 * t671 + t617 * t682;
t517 = -qJD(3) * t675 + t680;
t497 = -t616 * qJDD(3) + t517 * t614;
t498 = qJDD(3) * t614 + t517 * t616;
t687 = qJD(5) * t624;
t688 = qJD(5) * t620;
t410 = -t620 * t497 + t624 * t498 - t529 * t687 - t531 * t688;
t515 = qJDD(5) + t518;
t581 = -qJDD(1) * t599 + qJDD(2);
t440 = pkin(3) * t518 - qJ(4) * t517 - qJD(4) * t561 + t581;
t644 = -t621 * t546 + t547 * t729;
t453 = qJDD(3) * qJ(4) + (qJD(4) - t643) * qJD(3) + t644;
t405 = t616 * t440 - t453 * t614;
t395 = pkin(4) * t518 - pkin(8) * t498 + t405;
t406 = t614 * t440 + t616 * t453;
t403 = -pkin(8) * t497 + t406;
t668 = t624 * t395 - t403 * t620;
t629 = -qJD(5) * t400 + t668;
t379 = pkin(5) * t515 - pkin(9) * t410 + t629;
t411 = -qJD(5) * t741 + t624 * t497 + t498 * t620;
t639 = t620 * t395 + t624 * t403 + t418 * t687 - t430 * t688;
t380 = -pkin(9) * t411 + t639;
t669 = t623 * t379 - t619 * t380;
t742 = -g(1) * t535 + g(2) * t533 + t422 * t650 + t597 * t726 + t669;
t510 = qJDD(6) + t515;
t740 = t510 * MDP(30) + (-t424 ^ 2 + t650 ^ 2) * MDP(27) - t424 * MDP(26) * t650;
t739 = t738 * t605;
t645 = t738 * t607;
t642 = t676 - t702;
t516 = -pkin(3) * t642 - qJ(4) * t576 - t599;
t528 = -t621 * t584 + t586 * t729;
t465 = t616 * t516 - t528 * t614;
t705 = t576 * t616;
t439 = -pkin(4) * t642 - pkin(8) * t705 + t465;
t466 = t614 * t516 + t616 * t528;
t706 = t576 * t614;
t452 = -pkin(8) * t706 + t466;
t699 = t620 * t439 + t624 * t452;
t723 = pkin(8) + qJ(4);
t583 = t723 * t614;
t585 = t723 * t616;
t693 = -t620 * t583 + t624 * t585;
t737 = -t729 * t584 - t621 * t586;
t736 = MDP(4) * t617 - MDP(5) * t615;
t734 = qJ(2) * qJDD(1);
t520 = t623 * t573 + t575 * t619;
t698 = -qJD(6) * t520 - t619 * t695 - t623 * t696;
t733 = -t510 * t521 - t549 * t698;
t732 = -t515 * t575 + t553 * t696;
t634 = -t660 * t607 - t726;
t667 = t619 * t410 + t623 * t411;
t386 = -qJD(6) * t650 + t667;
t460 = t616 * t514 + t614 * t643;
t728 = pkin(8) * t616;
t431 = pkin(4) * t561 + t559 * t728 + t460;
t648 = qJD(4) * t614 + qJD(5) * t585;
t730 = t583 * t687 + t746 * t624 + (t431 + t648) * t620;
t554 = t559 ^ 2;
t399 = t624 * t418 - t430 * t620;
t389 = pkin(9) * t741 + t399;
t387 = pkin(5) * t553 + t389;
t719 = t387 * t623;
t718 = t390 * t623;
t717 = t424 * t561;
t716 = t650 * t561;
t715 = t474 * t561;
t714 = t741 * t561;
t710 = t518 * t614;
t709 = t518 * t616;
t564 = t642 * qJD(3);
t707 = t564 * t614;
t482 = pkin(3) * t565 - qJ(4) * t564 - qJD(4) * t576;
t492 = t642 * qJD(2) + qJD(3) * t737;
t438 = t614 * t482 + t616 * t492;
t692 = t615 ^ 2 + t617 ^ 2;
t685 = qJD(6) * t623;
t684 = -qJD(4) + t511;
t681 = t623 * t410 - t619 * t411 - t474 * t685;
t481 = -pkin(4) * t708 + t523;
t600 = -pkin(4) * t616 - pkin(3);
t673 = pkin(5) * t695 - t481;
t670 = t692 * qJD(1) ^ 2;
t437 = t616 * t482 - t492 * t614;
t415 = pkin(4) * t565 - t564 * t728 + t437;
t421 = -pkin(8) * t707 + t438;
t666 = t624 * t415 - t421 * t620;
t665 = t624 * t439 - t452 * t620;
t663 = -t624 * t583 - t585 * t620;
t662 = qJD(6) * t387 + t380;
t493 = qJD(2) * t576 - t584 * t690 + t586 * t674;
t661 = 0.2e1 * t692;
t658 = -t520 * t510 - t549 * t697;
t657 = -t573 * t515 - t553 * t695;
t429 = t624 * t431;
t501 = -pkin(9) * t573 + t693;
t656 = pkin(5) * t561 - pkin(9) * t696 + t575 * qJD(4) + qJD(5) * t693 + qJD(6) * t501 - t443 * t620 + t429;
t500 = -pkin(9) * t575 + t663;
t655 = pkin(9) * t695 - qJD(6) * t500 + t730;
t467 = pkin(4) * t707 + t493;
t382 = t387 * t619 + t718;
t653 = -t405 * t616 - t406 * t614;
t651 = -t449 * t614 + t450 * t616;
t505 = t575 * t576;
t506 = t573 * t576;
t456 = t623 * t505 - t506 * t619;
t457 = -t505 * t619 - t506 * t623;
t494 = pkin(4) * t706 - t737;
t646 = pkin(3) * t607 + qJ(4) * t605 + t599;
t638 = t620 * t415 + t624 * t421 + t439 * t687 - t452 * t688;
t385 = t686 * t741 + t681;
t632 = t458 * t576 + t511 * t564 - t660;
t628 = t661 * t683 - t660;
t419 = t497 * pkin(4) + t458;
t606 = cos(t612);
t604 = sin(t612);
t542 = pkin(5) * t573 + t600;
t541 = t604 * t622 + t606 * t703;
t540 = -t604 * t703 + t606 * t622;
t539 = t604 * t625 - t606 * t704;
t538 = t604 * t704 + t606 * t625;
t459 = t505 * pkin(5) + t494;
t455 = t564 * t575 + t687 * t705 - t688 * t706;
t454 = -t563 * t576 - t564 * t573;
t412 = pkin(5) * t455 + t467;
t404 = -pkin(9) * t505 + t699;
t401 = -pkin(5) * t642 + pkin(9) * t506 + t665;
t397 = qJD(6) * t457 + t454 * t619 + t623 * t455;
t396 = -qJD(6) * t456 + t454 * t623 - t455 * t619;
t391 = t411 * pkin(5) + t419;
t384 = -pkin(9) * t455 + t638;
t383 = pkin(5) * t565 - pkin(9) * t454 - qJD(5) * t699 + t666;
t381 = -t390 * t619 + t719;
t1 = [(-t385 * t456 - t386 * t457 - t396 * t424 + t397 * t650) * MDP(27) + (t385 * t457 - t396 * t650) * MDP(26) + (-t385 * t642 + t396 * t549 + t457 * t510 - t565 * t650) * MDP(28) + (-g(1) * t533 - g(2) * t535 - t382 * t565 + t459 * t385 - t388 * t642 + t391 * t457 + t422 * t396 - t412 * t650 + (-(-qJD(6) * t404 + t383) * t549 - t401 * t510 + t379 * t642) * t619 + (-(qJD(6) * t401 + t384) * t549 - t404 * t510 + t662 * t642) * t623) * MDP(32) + ((t383 * t623 - t384 * t619) * t549 + (t401 * t623 - t404 * t619) * t510 - t669 * t642 + t381 * t565 + t412 * t424 + t459 * t386 + t391 * t456 + t422 * t397 - g(1) * t534 - g(2) * t536 + ((-t401 * t619 - t404 * t623) * t549 + t382 * t642) * qJD(6)) * MDP(31) + (-t515 * t642 + t553 * t565) * MDP(23) + (t411 * t642 - t455 * t553 - t474 * t565 - t505 * t515) * MDP(22) + (-t510 * t642 + t549 * t565) * MDP(30) + (t386 * t642 - t397 * t549 - t424 * t565 - t456 * t510) * MDP(29) + (-qJD(3) * t565 + qJDD(3) * t642) * MDP(11) + (t517 * t642 - t518 * t576 - t559 * t564 - t561 * t565) * MDP(9) + t738 * MDP(2) + (-t405 * t642 + t437 * t559 + t449 * t565 + t465 * t518 + t493 * t529 - t497 * t737 + t614 * t632 + t616 * t645) * MDP(15) + (-qJD(3) * t493 + qJDD(3) * t737 - t518 * t599 + t565 * t582 - t581 * t642 + t645) * MDP(13) + (t406 * t642 - t438 * t559 - t450 * t565 - t466 * t518 + t493 * t531 - t498 * t737 - t614 * t645 + t616 * t632) * MDP(16) + (t405 * t465 + t406 * t466 + t449 * t437 + t450 * t438 - t458 * t737 + t511 * t493 + (-g(1) * t724 - g(2) * t646) * t625 + (g(1) * t646 - g(2) * t724) * t622) * MDP(18) + (-t410 * t505 + t411 * t506 - t454 * t474 + t455 * t741) * MDP(20) + (-t410 * t506 - t454 * t741) * MDP(19) + (-t410 * t642 + t454 * t553 - t506 * t515 - t565 * t741) * MDP(21) + (-g(1) * t538 - g(2) * t540 - t400 * t565 + t494 * t410 - t419 * t506 + t468 * t454 - t467 * t741 - t515 * t699 - t553 * t638 + t639 * t642) * MDP(25) + (t666 * t553 + t665 * t515 - t668 * t642 + t399 * t565 + t467 * t474 + t494 * t411 + t419 * t505 + t468 * t455 - g(1) * t539 - g(2) * t541 + (t400 * t642 - t553 * t699) * qJD(5)) * MDP(24) + qJDD(1) * MDP(1) + (qJD(3) * t564 + qJDD(3) * t576) * MDP(10) + (t517 * t576 + t561 * t564) * MDP(8) + t660 * MDP(3) + (pkin(1) * t647 + (t692 * t734 + t628) * qJ(2)) * MDP(7) + (t661 * t734 + t628) * MDP(6) + t736 * (t647 + t720) + (-qJD(3) * t492 - qJDD(3) * t528 - t517 * t599 + t564 * t582 + t576 * t581 - t739) * MDP(14) + (-t437 * t531 - t438 * t529 - t465 * t498 - t466 * t497 + t739 + t653 * t576 + (-t449 * t616 - t450 * t614) * t564) * MDP(17); -MDP(6) * t670 + (-qJ(2) * t670 - t647) * MDP(7) + (0.2e1 * qJD(3) * t561 + t654) * MDP(13) + ((-t559 - t675) * qJD(3) + t680) * MDP(14) + (-t529 * t561 - t554 * t614 + t709) * MDP(15) + (-t531 * t561 - t554 * t616 - t710) * MDP(16) + (-t497 * t614 - t498 * t616 + (-t529 * t616 + t531 * t614) * t559) * MDP(17) + (-t511 * t561 + t559 * t651 - t653 - t738) * MDP(18) + (t657 - t715) * MDP(24) + (t714 + t732) * MDP(25) + (t658 - t717) * MDP(31) + (t716 + t733) * MDP(32) - t736 * qJDD(1); (-qJ(4) * t709 - pkin(3) * t498 - t523 * t531 + (t616 * t684 + t461) * t559 + t630 * t614) * MDP(16) + (-qJ(4) * t710 - pkin(3) * t497 - t523 * t529 + (t614 * t684 - t460) * t559 - t630 * t616) * MDP(15) + (t523 * qJD(3) + t636 - t640) * MDP(13) + (t714 - t732) * MDP(21) + (t657 + t715) * MDP(22) + (t716 - t733) * MDP(28) + (t658 + t717) * MDP(29) + (-t449 * t460 - t450 * t461 - t511 * t523 + t651 * qJD(4) - t630 * pkin(3) + (-t405 * t614 + t406 * t616 + t634) * qJ(4)) * MDP(18) + (t663 * t515 + t600 * t411 + t419 * t573 - t481 * t474 + (-t429 - t648 * t624 + (qJD(5) * t583 + t746) * t620) * t553 + t695 * t468 + t636 * t606) * MDP(24) + (t410 * t575 + t696 * t741) * MDP(19) + (-t410 * t573 - t411 * t575 + t474 * t696 + t695 * t741) * MDP(20) + ((t559 - t675) * qJD(3) + t680) * MDP(10) + qJDD(3) * MDP(12) - t554 * MDP(9) + ((t500 * t623 - t501 * t619) * t510 + t542 * t386 + t391 * t520 + (t619 * t655 - t623 * t656) * t549 + t673 * t424 + t697 * t422 + t636 * t598) * MDP(31) + (-(t500 * t619 + t501 * t623) * t510 + t542 * t385 + t391 * t521 + (t619 * t656 + t623 * t655) * t549 - t673 * t650 + t698 * t422 - t636 * t597) * MDP(32) + (t385 * t521 - t650 * t698) * MDP(26) + (-t385 * t520 - t386 * t521 - t424 * t698 + t650 * t697) * MDP(27) + (t600 * t410 + t419 * t575 - t696 * t468 + t481 * t741 - t693 * t515 + t553 * t730 - t604 * t636) * MDP(25) - t518 * MDP(11) + (t582 * t559 - t634 - t644) * MDP(14) + (t460 * t531 + t461 * t529 + (-qJ(4) * t497 - qJD(4) * t529 - t449 * t559 + t406) * t616 + (qJ(4) * t498 + qJD(4) * t531 - t450 * t559 - t405) * t614 + t634) * MDP(17) + (MDP(11) * qJD(3) - MDP(13) * t582 - MDP(15) * t449 + MDP(16) * t450 - MDP(23) * t553 - MDP(24) * t399 + MDP(25) * t400 - MDP(30) * t549 - MDP(31) * t381 + MDP(32) * t382 + MDP(8) * t559 + MDP(9) * t561) * t561; (t531 * t559 + t497) * MDP(15) + (-t529 * t559 + t498) * MDP(16) + (-t529 ^ 2 - t531 ^ 2) * MDP(17) + (t449 * t531 + t450 * t529 + t630) * MDP(18) + (t411 - t748) * MDP(24) + (t410 - t749) * MDP(25) + (t386 - t744) * MDP(31) + (t385 - t745) * MDP(32); -t741 * t474 * MDP(19) + (-t474 ^ 2 + t741 ^ 2) * MDP(20) + (t410 + t749) * MDP(21) + (-t411 - t748) * MDP(22) + t515 * MDP(23) + (-g(1) * t540 + g(2) * t538 + t400 * t553 + t468 * t741 + t604 * t726 + t629) * MDP(24) + (g(1) * t541 - g(2) * t539 + t399 * t553 + t468 * t474 + t606 * t726 - t639) * MDP(25) + (t385 + t745) * MDP(28) + (-t386 - t744) * MDP(29) + (-(-t389 * t619 - t718) * t549 - t382 * qJD(6) + (t424 * t741 + t510 * t623 - t549 * t686) * pkin(5) + t742) * MDP(31) + ((-t390 * t549 - t379) * t619 + (t389 * t549 - t662) * t623 + (-t510 * t619 - t549 * t685 - t650 * t741) * pkin(5) + t743) * MDP(32) + t740; (t681 + t745) * MDP(28) + (-t667 - t744) * MDP(29) + (t382 * t549 + t742) * MDP(31) + (-t619 * t379 - t623 * t380 + t381 * t549 + t743) * MDP(32) + (MDP(28) * t713 + MDP(29) * t650 - MDP(31) * t382 - MDP(32) * t719) * qJD(6) + t740;];
tau  = t1;
