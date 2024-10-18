% Calculate Coriolis joint torque vector for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR15_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR15_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:32
% EndTime: 2024-09-27 22:26:39
% DurationCPUTime: 5.12s
% Computational Cost: add. (10450->440), mult. (35102->632), div. (0->0), fcn. (28421->12), ass. (0->245)
t572 = sin(qJ(3));
t573 = sin(qJ(2));
t567 = sin(pkin(5));
t671 = qJD(1) * t567;
t713 = cos(qJ(3));
t627 = t713 * t671;
t575 = cos(qJ(2));
t649 = t575 * t671;
t532 = -t572 * t649 - t573 * t627;
t571 = sin(qJ(4));
t712 = cos(qJ(4));
t650 = t573 * t671;
t741 = -t572 * t650 + t575 * t627;
t491 = t532 * t571 + t712 * t741;
t566 = sin(pkin(6));
t574 = cos(qJ(5));
t602 = -t532 * t712 + t571 * t741;
t664 = qJD(5) * t574;
t568 = cos(pkin(6));
t570 = sin(qJ(5));
t685 = t568 * t570;
t742 = -t491 * t574 + t566 * t664 + t602 * t685;
t606 = t572 * t575 + t573 * t713;
t721 = qJD(2) + qJD(3);
t583 = t721 * t606;
t579 = qJD(1) * t583;
t578 = t567 * t579;
t569 = cos(pkin(5));
t710 = pkin(1) * t575;
t660 = t569 * t710;
t555 = qJD(1) * t660;
t714 = pkin(8) + pkin(9);
t656 = t567 * t714;
t631 = t573 * t656;
t520 = -qJD(1) * t631 + t555;
t670 = qJD(1) * t569;
t557 = qJD(2) + t670;
t503 = pkin(2) * t557 + t520;
t548 = qJD(2) * t555;
t620 = qJD(2) * t631;
t511 = -qJD(1) * t620 + t548;
t711 = pkin(1) * t573;
t661 = t569 * t711;
t523 = (-t575 * t656 - t661) * qJD(2);
t512 = qJD(1) * t523;
t687 = t567 * t575;
t528 = t687 * t714 + t661;
t521 = t528 * qJD(1);
t644 = t713 * qJD(3);
t668 = qJD(3) * t572;
t629 = t503 * t644 + t713 * t511 + t572 * t512 - t521 * t668;
t423 = -pkin(10) * t578 + t629;
t493 = t741 * t721;
t516 = t713 * t521;
t608 = -t572 * t503 - t516;
t635 = -t572 * t511 + t713 * t512;
t586 = qJD(3) * t608 + t635;
t424 = -t493 * pkin(10) + t586;
t527 = t532 * pkin(10);
t513 = t572 * t521;
t636 = t713 * t503 - t513;
t461 = t527 + t636;
t552 = qJD(3) + t557;
t456 = pkin(3) * t552 + t461;
t707 = t741 * pkin(10);
t462 = -t608 + t707;
t645 = qJD(4) * t712;
t667 = qJD(4) * t571;
t630 = -t712 * t423 - t571 * t424 - t456 * t645 + t462 * t667;
t581 = t567 * t583;
t577 = t712 * t581;
t637 = qJD(1) * t577 + t571 * t493;
t434 = t602 * qJD(4) + t637;
t700 = t434 * t568;
t381 = -pkin(11) * t700 - t630;
t433 = t712 * t493 + t532 * t667 - t571 * t578 + t645 * t741;
t460 = t712 * t462;
t604 = -t571 * t456 - t460;
t642 = -t571 * t423 + t712 * t424;
t588 = qJD(4) * t604 + t642;
t708 = pkin(11) * t568;
t382 = -t433 * t708 + t588;
t609 = qJD(4) + t552;
t593 = t566 * t609;
t686 = t568 * t491;
t589 = -t593 - t686;
t407 = -pkin(11) * t589 - t604;
t482 = t602 * t708;
t458 = t571 * t462;
t640 = t712 * t456 - t458;
t409 = t640 - t482;
t408 = pkin(4) * t609 + t409;
t541 = (-pkin(2) * t575 - pkin(1)) * t567;
t539 = qJD(1) * t541;
t500 = -pkin(3) * t741 + t539;
t726 = t602 * t566;
t435 = -pkin(4) * t491 - pkin(11) * t726 + t500;
t618 = t408 * t568 + t435 * t566;
t383 = -t407 * t570 + t574 * t618;
t669 = qJD(2) * t567;
t648 = t573 * t669;
t625 = qJD(1) * t648;
t547 = pkin(2) * t625;
t476 = pkin(3) * t578 + t547;
t561 = t566 * pkin(11);
t401 = t434 * pkin(4) - t433 * t561 + t476;
t370 = (t382 * t568 + t401 * t566) * t570 + t381 * t574 + qJD(5) * t383;
t374 = -t382 * t566 + t401 * t568;
t396 = -t408 * t566 + t435 * t568;
t689 = t566 * t570;
t594 = -t370 * t568 + t374 * t689 + t742 * t396;
t430 = t434 * t685;
t665 = qJD(5) * t570;
t684 = t568 * t574;
t737 = t491 * t684 + t574 * t593;
t392 = t737 * qJD(5) + t574 * t433 - t602 * t665 - t430;
t681 = t574 * t602;
t719 = t570 * t589 - t681;
t393 = -t719 * qJD(5) + t570 * t433 + t434 * t684;
t444 = t570 * t602 - t737;
t688 = t566 * t574;
t740 = t392 * t688 - t393 * t689 - t742 * t444;
t695 = t491 * t566;
t474 = -t568 * t609 - qJD(5) + t695;
t562 = t566 ^ 2;
t672 = MDP(29) * t566;
t699 = t434 * t570;
t739 = (t474 * t602 + t700) * t672 + (t602 * t552 - t637) * MDP(21) + (-t491 ^ 2 + t602 ^ 2) * MDP(19) + (-t491 * t609 + t433) * MDP(20) - t491 * t602 * MDP(18) + (t392 * t568 - t742 * t474 + t562 * t699 + t719 * t726) * MDP(27) + (t392 * t689 - t742 * t719) * MDP(25);
t738 = t491 * t570;
t639 = -t461 * t571 - t460;
t659 = t491 * t708;
t411 = t639 - t659;
t709 = pkin(3) * t532;
t730 = pkin(4) * t602;
t442 = -t491 * t561 - t709 + t730;
t736 = -t442 * t568 + (t667 * pkin(3) + t411) * t566;
t735 = t434 * t562 * t574 - t393 * t568 + t444 * t726;
t598 = -t500 * t491 + t630;
t729 = MDP(5) * (t573 ^ 2 - t575 ^ 2);
t384 = t407 * t574 + t570 * t618;
t728 = t384 * t602;
t727 = t573 * MDP(4);
t634 = -t520 * t572 - t516;
t467 = t634 - t707;
t674 = t713 * t520 - t513;
t468 = t527 + t674;
t560 = pkin(2) * t713 + pkin(3);
t683 = t571 * t572;
t724 = -t560 * t645 - (-t572 * t667 + (t712 * t713 - t683) * qJD(3)) * pkin(2) + t571 * t467 + t712 * t468;
t723 = -t572 * t573 + t713 * t575;
t720 = -t500 * t602 + t642;
t519 = (pkin(2) + t710) * t569 - t631;
t535 = t606 * t567;
t469 = t569 * pkin(3) - t535 * pkin(10) + t519 * t713 - t572 * t528;
t534 = t723 * t567;
t607 = -t572 * t519 - t528 * t713;
t473 = t534 * pkin(10) - t607;
t603 = -t571 * t469 - t473 * t712;
t601 = t534 * t712 - t571 * t535;
t615 = t566 * t569 + t568 * t601;
t415 = pkin(11) * t615 - t603;
t496 = t571 * t534 + t535 * t712;
t416 = t569 * pkin(4) + t469 * t712 - t571 * t473 - t496 * t708;
t507 = -pkin(3) * t534 + t541;
t443 = -pkin(4) * t601 - t496 * t561 + t507;
t616 = t416 * t568 + t443 * t566;
t718 = t415 * t574 + t570 * t616;
t716 = t567 * t721;
t706 = pkin(3) * qJD(4);
t701 = t434 * t566;
t497 = t723 * t716;
t439 = qJD(4) * t496 + t571 * t497 + t577;
t698 = t439 * t566;
t697 = t474 * t570;
t693 = t496 * t570;
t692 = t539 * t532;
t559 = pkin(3) * t712 + pkin(4);
t691 = t559 * t566;
t563 = t567 ^ 2;
t576 = qJD(1) ^ 2;
t690 = t563 * t576;
t371 = -t384 * qJD(5) - t381 * t570 + t382 * t684 + t401 * t688;
t647 = t566 * t665;
t680 = t371 * t568 + t396 * t647;
t638 = t712 * t467 - t468 * t571;
t413 = t638 - t659;
t651 = t712 * t572;
t502 = -t560 * t667 + (-t572 * t645 + (-t571 * t713 - t651) * qJD(3)) * pkin(2);
t678 = t413 - t502;
t677 = t482 - t724;
t676 = t712 * t461 - t458;
t666 = qJD(5) * t474;
t658 = pkin(11) * t689;
t657 = t563 * t711;
t655 = t575 * t690;
t643 = qJD(1) * qJD(2) * t563;
t554 = pkin(2) * t648;
t553 = pkin(2) * t650;
t626 = t575 * t643;
t447 = t568 * t681 + t738;
t622 = -t447 + t647;
t448 = t602 * t684 + t738;
t621 = -t448 + t647;
t617 = t411 * t568 + t442 * t566;
t612 = -t383 * t726 - t396 * t448 + t680;
t549 = pkin(3) * t571 + t561;
t611 = -t549 * t570 + t559 * t684;
t610 = t549 * t574 + t559 * t685;
t600 = -pkin(8) * t687 - t661;
t599 = -pkin(8) * t625 + t548;
t597 = -t539 * t741 - t629;
t556 = qJD(2) * t660;
t522 = t556 - t620;
t596 = t519 * t644 + t713 * t522 + t572 * t523 - t528 * t668;
t428 = -pkin(10) * t581 + t596;
t585 = qJD(3) * t607 - t572 * t522 + t713 * t523;
t429 = -t497 * pkin(10) + t585;
t595 = t712 * t428 + t571 * t429 + t469 * t645 - t473 * t667;
t592 = t600 * t557;
t590 = -t415 * t570 + t574 * t616;
t452 = t496 * t574 + t570 * t615;
t587 = qJD(4) * t603 - t571 * t428 + t712 * t429;
t584 = -t456 * t667 - t462 * t645 + t720;
t582 = t532 * t741 * MDP(11) + (t621 * t719 + t740) * MDP(26) + (t474 * t621 + t735) * MDP(28) + (-t552 * t741 + t493) * MDP(13) + (-t532 * t552 - t578) * MDP(14) + (t532 ^ 2 - t741 ^ 2) * MDP(12) + t739;
t483 = pkin(3) * t581 + t554;
t537 = -pkin(2) * t683 + t560 * t712 + pkin(4);
t533 = pkin(2) * t651 + t571 * t560 + t561;
t508 = t553 - t709;
t481 = t566 * t601 - t568 * t569;
t453 = -pkin(11) * t695 + t730;
t451 = -t569 * t688 - t601 * t684 + t693;
t438 = qJD(4) * t601 + t497 * t712 - t571 * t581;
t437 = t442 + t553;
t412 = -t482 + t676;
t410 = -pkin(11) * t686 + t604;
t405 = t439 * pkin(4) - t438 * t561 + t483;
t404 = -t416 * t566 + t443 * t568;
t403 = -t410 * t566 + t453 * t568;
t402 = -t413 * t566 + t437 * t568;
t398 = qJD(5) * t452 + t438 * t570 + t439 * t684;
t397 = -t439 * t685 + t438 * t574 + (t574 * t615 - t693) * qJD(5);
t386 = -t438 * t708 + t587;
t385 = -t439 * t708 + t595;
t377 = -t386 * t566 + t405 * t568;
t1 = [-0.2e1 * t643 * t729 + (t507 * t433 + t500 * t438 + t476 * t496 + t483 * t602 + t569 * t630 - t595 * t609) * MDP(24) + (-t596 * t552 - t629 * t569 + t541 * t493 + t539 * t497 + (qJD(1) * t535 - t532) * t554) * MDP(17) + (0.2e1 * t539 * t606 * t716 - t534 * t547 + t552 * t585 - t554 * t741 + t569 * t586) * MDP(16) + (-0.2e1 * pkin(1) * t626 - (-pkin(8) * t648 + t556) * t557 - t599 * t569) * MDP(10) + (t433 * t601 - t434 * t496 + t438 * t491 - t439 * t602) * MDP(19) + (t433 * t496 + t438 * t602) * MDP(18) + (t392 * t452 - t397 * t719) * MDP(25) + (-t392 * t451 - t393 * t452 - t397 * t444 + t398 * t719) * MDP(26) + (t493 * t569 + t497 * t552) * MDP(13) + (t493 * t535 - t497 * t532) * MDP(11) + (-t552 * t583 - t569 * t579) * t567 * MDP(14) + (-t434 * t481 - t439 * t474) * t672 + (t393 * t481 + t398 * t474 + (-t434 * t451 - t439 * t444) * t566) * MDP(28) + (-(-t385 * t570 + t386 * t684 + t405 * t688) * t474 + t590 * t701 - t371 * t481 + t383 * t698 + t377 * t444 + t404 * t393 + t374 * t451 + t396 * t398 + t718 * t666) * MDP(30) + ((t385 * t574 + t386 * t685 + t405 * t689) * t474 - t718 * t701 + t370 * t481 - t384 * t698 - t377 * t719 + t404 * t392 + t374 * t452 + t396 * t397 + t590 * t666) * MDP(31) + (-t392 * t481 - t397 * t474 + (t434 * t452 - t439 * t719) * t566) * MDP(27) + (t592 + (t569 * t600 - 0.2e1 * t657) * qJD(1)) * qJD(2) * MDP(9) + 0.2e1 * t626 * t727 + (t433 * t569 + t438 * t609) * MDP(20) + (-t434 * t569 - t439 * t609) * MDP(21) + (t507 * t434 + t500 * t439 - t476 * t601 - t483 * t491 + t569 * t588 + t587 * t609) * MDP(23) + (t493 * t534 + t497 * t741 + (t532 * t583 - t535 * t579) * t567) * MDP(12) + (t575 * MDP(6) * t669 - MDP(7) * t648) * (t557 + t670); t690 * t729 + (t402 * t719 + (-t533 * t665 + t677 * t574 + (t537 * t664 - t570 * t678) * t568) * t474 + (-(t533 * t574 + t537 * t685) * t434 + t502 * t719 - t537 * t392 - t437 * t697 + t728) * t566 + t594) * MDP(31) + (pkin(1) * t655 + (-pkin(8) * t650 + t555) * t557 - t599) * MDP(10) + (t491 * t508 + t584) * MDP(23) + (-t508 * t602 + t598) * MDP(24) + (t576 * t657 + (qJD(2) * t600 - t592) * qJD(1)) * MDP(9) + (-t402 * t444 + (t533 * t664 + t677 * t570 + (t537 * t665 + t574 * t678) * t568) * t474 + (-t533 * t699 - t537 * t393 - t502 * t444 + (t437 * t474 + t537 * t700 - t374) * t574) * t566 + t612) * MDP(30) + (t674 * t552 + (t532 * t650 - t552 * t644) * pkin(2) + t597) * MDP(17) + (-t634 * t552 + t741 * t553 + t692 + (-t516 + (-pkin(2) * t552 - t503) * t572) * qJD(3) + t635) * MDP(16) - t655 * t727 + t582 + (MDP(6) * t649 - MDP(7) * t650) * (qJD(2) - t557) + ((t502 - t638) * MDP(23) + t724 * MDP(24)) * t609; (t552 * t636 + t597) * MDP(17) + (-t374 * t688 - t393 * t691 + t736 * t444 + t611 * t701 + t612) * MDP(30) + (t676 * t609 + (t532 * t602 - t609 * t645) * pkin(3) + t598) * MDP(24) + (-t639 * t609 + (-t491 * t532 - t609 * t667) * pkin(3) + t584) * MDP(23) + (t384 * t726 - t392 * t691 - t610 * t701 - t719 * t736 + t594) * MDP(31) + t582 + (-t552 * t608 + t586 + t692) * MDP(16) + ((t610 * qJD(5) - (-t570 * t712 - t571 * t684) * t706 - t412 * t570 + t617 * t574) * MDP(30) + (t611 * qJD(5) + (-t571 * t685 + t574 * t712) * t706 - t412 * t574 - t617 * t570) * MDP(31)) * t474; (-t552 * t604 + t720) * MDP(23) + (t609 * t640 + t598) * MDP(24) + (t622 * t719 + t740) * MDP(26) + (t474 * t622 + t735) * MDP(28) + (-t396 * t447 - t403 * t444 + (-t409 * t570 + (pkin(4) * t665 + t410 * t574) * t568) * t474 + (-t434 * t658 - pkin(4) * t393 - t383 * t602 + (pkin(4) * t700 - t374 + (pkin(11) * qJD(5) + t453) * t474) * t574) * t566 + t680) * MDP(30) + (-(t409 * t574 + t410 * t685) * t474 + t403 * t719 + (-pkin(11) * t434 * t688 - t453 * t697 + t728 + (-t392 - t430) * pkin(4)) * t566 + (pkin(4) * t684 - t658) * t666 + t594) * MDP(31) + t739; -t719 * t444 * MDP(25) + (-t444 ^ 2 + t719 ^ 2) * MDP(26) + (-t444 * t474 + t392) * MDP(27) + (t474 * t719 - t393) * MDP(28) + t434 * t672 + (-t384 * t474 + t396 * t719 + t371) * MDP(30) + (-t383 * t474 + t396 * t444 - t370) * MDP(31);];
tauc = t1;
