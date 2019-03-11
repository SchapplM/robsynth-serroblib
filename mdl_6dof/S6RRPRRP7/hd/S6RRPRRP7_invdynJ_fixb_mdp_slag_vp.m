% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:24
% EndTime: 2019-03-09 12:20:36
% DurationCPUTime: 10.40s
% Computational Cost: add. (7121->624), mult. (14770->741), div. (0->0), fcn. (9917->8), ass. (0->254)
t565 = sin(qJ(4));
t569 = cos(qJ(4));
t570 = cos(qJ(2));
t673 = qJD(1) * t570;
t566 = sin(qJ(2));
t674 = qJD(1) * t566;
t491 = -t565 * t673 + t569 * t674;
t564 = sin(qJ(5));
t568 = cos(qJ(5));
t689 = t569 * t570;
t495 = t565 * t566 + t689;
t600 = t495 * qJD(4);
t450 = qJD(2) * t495 - t600;
t614 = t570 * t565 - t566 * t569;
t596 = t614 * qJDD(1);
t576 = t450 * qJD(1) - t596;
t657 = qJD(2) - qJD(4);
t635 = t568 * t657;
t656 = qJDD(2) - qJDD(4);
t665 = qJD(5) * t564;
t405 = qJD(5) * t635 + t491 * t665 + t564 * t656 - t568 * t576;
t460 = t568 * t491 - t564 * t657;
t667 = qJD(5) * t460;
t406 = t564 * t576 + t568 * t656 + t667;
t660 = qJD(1) * qJD(2);
t649 = t566 * t660;
t668 = qJD(4) * t569;
t669 = qJD(4) * t565;
t670 = qJD(2) * t570;
t761 = t565 * t670 + t566 * t668 - t570 * t669;
t426 = qJD(1) * t761 + qJDD(1) * t495 - t569 * t649;
t458 = t491 * t564 + t635;
t488 = t495 * qJD(1);
t423 = qJDD(5) + t426;
t690 = t568 * t423;
t696 = t564 * t423;
t708 = t405 * t564;
t750 = qJD(5) + t488;
t741 = t750 * t460;
t743 = t750 ^ 2;
t763 = t458 * t750;
t764 = -(t491 * t657 + t426) * MDP(18) - (t564 * (t406 + t741) + (t405 + t763) * t568) * MDP(23) + (t568 * t741 - t708) * MDP(22) + (-t460 * t491 + t568 * t743 + t696) * MDP(24) - (-t458 * t491 + t564 * t743 - t690) * MDP(25) - (t488 ^ 2 - t491 ^ 2) * MDP(16) - t656 * MDP(19) + (MDP(15) * t488 - MDP(26) * t750) * t491;
t493 = -qJD(1) * pkin(1) - pkin(2) * t673 - qJ(3) * t674;
t468 = pkin(3) * t673 - t493;
t421 = pkin(4) * t488 - pkin(9) * t491 + t468;
t731 = pkin(2) + pkin(3);
t652 = t731 * qJD(2);
t543 = pkin(7) * t674;
t500 = pkin(8) * t674 - t543;
t759 = qJD(3) - t500;
t472 = -t652 + t759;
t544 = pkin(7) * t673;
t502 = -pkin(8) * t673 + t544;
t559 = qJD(2) * qJ(3);
t492 = t502 + t559;
t435 = t565 * t472 + t569 * t492;
t430 = -pkin(9) * t657 + t435;
t399 = t421 * t564 + t430 * t568;
t393 = qJ(6) * t750 + t399;
t762 = t750 * t393;
t567 = sin(qJ(1));
t479 = t614 * t567;
t571 = cos(qJ(1));
t686 = t570 * t571;
t693 = t566 * t571;
t481 = t565 * t686 - t569 * t693;
t598 = g(1) * t481 + g(2) * t479 + g(3) * t495;
t624 = pkin(5) * t568 + qJ(6) * t564;
t611 = pkin(4) + t624;
t756 = t598 * t611;
t648 = t570 * t660;
t659 = qJDD(1) * t566;
t753 = t648 + t659;
t677 = t569 * qJ(3) - t565 * t731;
t499 = -pkin(9) + t677;
t697 = t499 * t423;
t645 = -t565 * qJ(3) - t569 * t731;
t473 = qJD(3) * t569 + qJD(4) * t645;
t700 = t473 * t750;
t434 = t569 * t472 - t565 * t492;
t429 = pkin(4) * t657 - t434;
t402 = t458 * pkin(5) - t460 * qJ(6) + t429;
t742 = t750 * t402;
t752 = t742 + t697 + t700;
t727 = pkin(9) * t423;
t751 = t742 - t727;
t676 = t570 * pkin(2) + t566 * qJ(3);
t745 = -pkin(1) - t676;
t744 = t429 * t750;
t447 = t500 * t569 + t502 * t565;
t740 = t447 - t473;
t552 = t570 * pkin(3);
t654 = t552 + t676;
t613 = pkin(9) * t614 + t654;
t433 = pkin(4) * t495 + pkin(1) + t613;
t730 = pkin(7) - pkin(8);
t508 = t730 * t566;
t509 = t730 * t570;
t463 = t508 * t565 + t509 * t569;
t682 = t564 * t433 + t568 * t463;
t681 = qJD(4) * t677 + t502 * t569 + t565 * t759;
t626 = g(1) * t571 + g(2) * t567;
t717 = pkin(7) * qJDD(2);
t736 = qJD(2) * (qJD(1) * t745 + t493) - t717;
t735 = t657 ^ 2;
t623 = pkin(5) * t564 - qJ(6) * t568;
t734 = -qJD(6) * t564 + t623 * t750;
t732 = t460 ^ 2;
t729 = pkin(5) * t423;
t728 = pkin(5) * t491;
t725 = g(1) * t567;
t721 = g(2) * t571;
t719 = g(3) * t614;
t718 = pkin(9) * qJD(5);
t716 = qJ(6) * t423;
t715 = qJ(6) * t491;
t562 = qJDD(1) * pkin(1);
t714 = qJDD(2) * pkin(2);
t398 = t421 * t568 - t430 * t564;
t661 = qJD(6) - t398;
t392 = -pkin(5) * t750 + t661;
t713 = t392 * t491;
t712 = t393 * t491;
t711 = t398 * t491;
t710 = t399 * t750;
t709 = t399 * t491;
t706 = t406 * t568;
t705 = t458 * t460;
t704 = t458 * t564;
t703 = t458 * t568;
t702 = t460 * t564;
t701 = t460 * t568;
t699 = t614 * t564;
t698 = t614 * t568;
t695 = t566 * t567;
t574 = qJD(1) ^ 2;
t692 = t566 * t574;
t691 = t567 * t570;
t442 = pkin(4) * t491 + pkin(9) * t488;
t685 = t568 * t434 + t564 * t442;
t684 = -t734 + t681;
t533 = qJ(3) * t673;
t478 = -t674 * t731 + t533;
t424 = -t442 + t478;
t683 = t564 * t424 + t568 * t447;
t680 = -t435 + t734;
t547 = t566 * qJD(3);
t678 = qJ(3) * t670 + t547;
t560 = t566 ^ 2;
t561 = t570 ^ 2;
t675 = t560 - t561;
t672 = qJD(2) * t566;
t671 = qJD(2) * t569;
t666 = qJD(5) * t499;
t664 = qJD(5) * t568;
t658 = qJDD(1) * t570;
t655 = t570 * t692;
t653 = -g(1) * t693 - g(2) * t695 + g(3) * t570;
t538 = pkin(7) * t659;
t647 = pkin(7) * t648 + qJDD(3) + t538;
t646 = -qJD(2) * pkin(2) + qJD(3);
t480 = t495 * t567;
t452 = t480 * t564 - t571 * t568;
t622 = pkin(2) * t658 + qJ(3) * t753 + qJD(1) * t547 + t562;
t609 = pkin(3) * t658 + t622;
t390 = t426 * pkin(4) + pkin(9) * t596 + (pkin(9) * t600 + (-pkin(9) * t689 + (-pkin(9) * t565 - t731) * t566) * qJD(2)) * qJD(1) + t609;
t444 = -pkin(8) * t753 - t731 * qJDD(2) + t647;
t539 = pkin(7) * t658;
t557 = qJDD(2) * qJ(3);
t558 = qJD(2) * qJD(3);
t470 = -pkin(7) * t649 + t539 + t557 + t558;
t445 = (t649 - t658) * pkin(8) + t470;
t603 = t565 * t444 + t569 * t445 + t472 * t668 - t492 * t669;
t396 = -pkin(9) * t656 + t603;
t634 = -t568 * t390 + t564 * t396 + t421 * t665 + t430 * t664;
t633 = t571 * pkin(1) + pkin(2) * t686 + t567 * pkin(7) + qJ(3) * t693;
t632 = -t538 - t653;
t631 = t566 * t652;
t573 = qJD(2) ^ 2;
t630 = pkin(7) * t573 + t721;
t482 = t495 * t571;
t456 = t482 * t564 + t567 * t568;
t629 = -g(1) * t452 + g(2) * t456;
t453 = t480 * t568 + t564 * t571;
t457 = t482 * t568 - t567 * t564;
t628 = g(1) * t453 - g(2) * t457;
t627 = g(1) * t479 - g(2) * t481;
t621 = t392 * t568 - t393 * t564;
t620 = t392 * t564 + t393 * t568;
t504 = t543 + t646;
t507 = t544 + t559;
t616 = t504 * t570 - t507 * t566;
t615 = t508 * t569 - t509 * t565;
t610 = -t569 * t444 + t565 * t445 + t472 * t669 + t492 * t668;
t477 = t647 - t714;
t608 = -0.2e1 * pkin(1) * t660 - t717;
t606 = t450 * t564 - t614 * t664;
t605 = -t450 * t568 - t614 * t665;
t604 = -t727 + t744;
t466 = -t631 + t678;
t602 = t564 * t390 + t568 * t396 + t421 * t664 - t430 * t665;
t449 = -t566 * t671 + t761;
t410 = pkin(4) * t449 - pkin(9) * t450 + t466;
t501 = t730 * t672;
t503 = qJD(2) * t509;
t417 = qJD(4) * t615 - t501 * t569 + t503 * t565;
t601 = t564 * t410 + t568 * t417 + t433 * t664 - t463 * t665;
t599 = -t697 - t744;
t597 = g(1) * t482 + g(2) * t480 - t719;
t595 = -t630 + 0.2e1 * t562;
t397 = pkin(4) * t656 + t610;
t592 = -t397 + t598;
t591 = g(1) * t456 + g(2) * t452 - g(3) * t699 - t634;
t590 = -t718 * t750 + t598;
t589 = -t666 * t750 - t598;
t440 = pkin(2) * t649 - t622;
t484 = pkin(2) * t672 - t678;
t588 = -qJD(1) * t484 - qJDD(1) * t745 - t440 - t630;
t384 = t406 * pkin(5) + t405 * qJ(6) - t460 * qJD(6) + t397;
t587 = -t384 + t590;
t586 = t384 + t589;
t382 = qJD(6) * t750 + t602 + t716;
t383 = qJDD(6) + t634 - t729;
t585 = qJD(5) * t621 + t382 * t568 + t383 * t564;
t584 = qJD(2) * t616 + t470 * t570 + t477 * t566;
t418 = qJD(4) * t463 - t501 * t565 - t503 * t569;
t583 = t402 * t460 + qJDD(6) - t591;
t582 = -g(1) * t457 - g(2) * t453 + g(3) * t698 + t602;
t581 = -t468 * t491 + t598 - t610;
t580 = t468 * t488 + t597 - t603;
t554 = t571 * pkin(7);
t527 = g(1) * t691;
t521 = qJ(3) * t686;
t519 = qJ(3) * t691;
t498 = pkin(4) - t645;
t497 = pkin(2) * t674 - t533;
t494 = pkin(1) + t654;
t490 = t564 * t674 + t568 * t671;
t487 = t564 * t671 - t568 * t674;
t469 = t611 - t645;
t427 = -qJD(1) * t631 + t609;
t419 = pkin(5) * t460 + qJ(6) * t458;
t414 = -t614 * t623 - t615;
t408 = -pkin(5) * t495 - t433 * t568 + t463 * t564;
t407 = qJ(6) * t495 + t682;
t404 = t434 * t564 - t442 * t568 - t728;
t403 = t685 + t715;
t401 = -t424 * t568 + t447 * t564 + t728;
t400 = t683 - t715;
t391 = -t405 + t763;
t387 = t623 * t450 - (qJD(5) * t624 - qJD(6) * t568) * t614 + t418;
t386 = -pkin(5) * t449 + qJD(5) * t682 - t410 * t568 + t417 * t564;
t385 = qJ(6) * t449 + qJD(6) * t495 + t601;
t1 = [(-t385 * t458 + t386 * t460 - t405 * t408 - t406 * t407 + t621 * t450 - (-qJD(5) * t620 - t382 * t564 + t383 * t568) * t614 + t627) * MDP(30) + ((-t702 - t703) * t450 - (t708 - t706 + (-t701 + t704) * qJD(5)) * t614) * MDP(23) + (t491 * t450 - t576 * t614) * MDP(15) + (t417 * t657 - t427 * t614 + t468 * t450 + t463 * t656 + t466 * t491 + t494 * t576 - t627) * MDP(21) + (t426 * t614 - t491 * t449 - t450 * t488 - t495 * t576) * MDP(16) + (t382 * t407 + t393 * t385 + t384 * t414 + t402 * t387 + t383 * t408 + t392 * t386 - g(1) * (-pkin(4) * t480 - pkin(5) * t453 - pkin(8) * t571 - pkin(9) * t479 - qJ(6) * t452 + t554) - g(2) * (pkin(3) * t686 + pkin(4) * t482 + pkin(5) * t457 + pkin(9) * t481 + qJ(6) * t456 + t633) + (-g(1) * (t745 - t552) + g(2) * pkin(8)) * t567) * MDP(32) + (g(1) * t480 - g(2) * t482 + t418 * t657 + t494 * t426 + t427 * t495 + t468 * t449 + t466 * t488 - t615 * t656) * MDP(20) + (-t450 * t657 + t614 * t656) * MDP(17) + (t566 * t736 + t588 * t570 + t527) * MDP(11) + (-t736 * t570 + (t588 + t725) * t566) * MDP(13) + (pkin(7) * t584 - g(1) * t554 - g(2) * t633 + t493 * t484 + (t440 - t725) * t745) * MDP(14) + (qJDD(1) * t560 + 0.2e1 * t566 * t648) * MDP(4) + (-t406 * t495 - t449 * t458 - t606 * t750 + t614 * t696) * MDP(25) + (-t405 * t495 + t449 * t460 - t605 * t750 - t614 * t690) * MDP(24) + (-t397 * t698 - t399 * t449 + t405 * t615 + t418 * t460 - t423 * t682 - t429 * t605 - t495 * t602 - t601 * t750 + t629) * MDP(28) + (-t634 * t495 + t398 * t449 + t418 * t458 - t615 * t406 + ((-qJD(5) * t463 + t410) * t750 + t433 * t423 - t429 * qJD(5) * t614) * t568 + ((-qJD(5) * t433 - t417) * t750 - t463 * t423 - t397 * t614 + t429 * t450) * t564 + t628) * MDP(27) + (t423 * t495 + t449 * t750) * MDP(26) + (t382 * t495 + t384 * t698 + t385 * t750 - t387 * t460 + t393 * t449 + t402 * t605 + t405 * t414 + t407 * t423 - t629) * MDP(31) + (-t383 * t495 - t384 * t699 - t386 * t750 + t387 * t458 - t392 * t449 + t402 * t606 + t406 * t414 - t408 * t423 + t628) * MDP(29) + (t608 * t570 + (-t595 - t725) * t566) * MDP(10) + 0.2e1 * (t566 * t658 - t660 * t675) * MDP(5) + (-t721 + t725) * MDP(2) + (t449 * t657 + t495 * t656) * MDP(18) + ((t560 + t561) * qJDD(1) * pkin(7) + t584 - t626) * MDP(12) + t626 * MDP(3) + (t566 * t608 + t570 * t595 + t527) * MDP(9) + qJDD(1) * MDP(1) + (qJDD(2) * t566 + t570 * t573) * MDP(6) + (qJDD(2) * t570 - t566 * t573) * MDP(7) + (t405 * t698 - t460 * t605) * MDP(22); -t764 + t596 * MDP(17) + (t384 * t469 - g(1) * (-pkin(9) * t482 + t521) - g(2) * (-pkin(9) * t480 + t519) - g(3) * t613 + t684 * t402 + (t473 * t568 - t400) * t393 + (t473 * t564 - t401) * t392 + t585 * t499 + t626 * t566 * t731 - t756) * MDP(32) + (0.2e1 * t714 - qJDD(3) + (-t493 * t566 + t497 * t570) * qJD(1) + t632) * MDP(11) + (t470 * qJ(3) + t507 * qJD(3) - t477 * pkin(2) - t493 * t497 - g(1) * (-pkin(2) * t693 + t521) - g(2) * (-pkin(2) * t695 + t519) - g(3) * t676 - t616 * qJD(1) * pkin(7)) * MDP(14) + (pkin(1) * t692 + t632) * MDP(9) - MDP(4) * t655 + (-t400 * t750 + t405 * t469 - t684 * t460 + t586 * t564 + t568 * t752 + t712) * MDP(31) + (t401 * t750 + t406 * t469 + t684 * t458 - t564 * t752 + t586 * t568 - t713) * MDP(29) + (-t498 * t405 + t683 * t750 - t709 + t681 * t460 + (t599 - t700) * t568 + (-t397 - t589) * t564) * MDP(28) + (t711 + t498 * t406 + t681 * t458 + (t740 * t750 + t599) * t564 + ((-t424 - t666) * t750 - t592) * t568) * MDP(27) + ((-t566 * pkin(2) + qJ(3) * t570) * qJDD(1) + ((t507 - t559) * t566 + (-t504 + t646) * t570) * qJD(1)) * MDP(12) + (t400 * t458 - t401 * t460 + (-t392 * t488 - t406 * t499 - t458 * t473 - t382 + (t460 * t499 - t392) * qJD(5)) * t568 + (t393 * t488 - t405 * t499 + t460 * t473 - t383 + (t458 * t499 + t393) * qJD(5)) * t564 + t597) * MDP(30) + (t539 + 0.2e1 * t557 + 0.2e1 * t558 + (qJD(1) * t497 - g(3)) * t566 + (qJD(1) * t493 - t626) * t570) * MDP(13) + (g(3) * t566 - t539 + (pkin(1) * t574 + t626) * t570) * MDP(10) + qJDD(2) * MDP(8) + t675 * MDP(5) * t574 + MDP(7) * t658 + MDP(6) * t659 + (-t478 * t488 - t645 * t656 + t657 * t681 - t581) * MDP(20) + (-t478 * t491 + t656 * t677 - t657 * t740 - t580) * MDP(21); (-qJDD(2) - t655) * MDP(11) + MDP(12) * t659 + (-t560 * t574 - t573) * MDP(13) + (-qJD(2) * t507 + t493 * t674 + t477 + t653) * MDP(14) + (-t488 * t674 - t565 * t735 - t569 * t656) * MDP(20) + (-t491 * t674 + t565 * t656 - t569 * t735) * MDP(21) + (t458 * t490 - t460 * t487 + (t702 - t703) * t668 + (-t708 - t706 + (t701 + t704) * qJD(5)) * t565) * MDP(30) + (-t392 * t487 - t393 * t490 + (qJD(4) * t620 - t384) * t569 + (-t402 * t657 + t585) * t565 + t653) * MDP(32) + (MDP(27) + MDP(29)) * (-t406 * t569 + t458 * t669 + (-qJD(2) * t458 - t696) * t565 + (-t564 * t668 - t565 * t664 + t487) * t750) + (-MDP(28) + MDP(31)) * (t565 * (t460 * t657 - t665 * t750 + t690) + (t568 * t668 - t490) * t750 - t405 * t569); (-t488 * t657 + t576) * MDP(17) + (-t435 * t657 + t581) * MDP(20) + (-t434 * t657 + t580) * MDP(21) + (-pkin(4) * t406 - t711 - t435 * t458 + (t434 * t750 + t604) * t564 + ((-t442 - t718) * t750 + t592) * t568) * MDP(27) + (pkin(4) * t405 + t685 * t750 + t709 - t435 * t460 + t604 * t568 + (t397 - t590) * t564) * MDP(28) + (t404 * t750 - t406 * t611 + t680 * t458 + t564 * t751 + t587 * t568 + t713) * MDP(29) + (t403 * t458 - t404 * t460 + (t382 + t750 * t392 + (-t406 + t667) * pkin(9)) * t568 + (t383 - t762 + (qJD(5) * t458 - t405) * pkin(9)) * t564 - t597) * MDP(30) + (-t403 * t750 - t405 * t611 - t680 * t460 + t587 * t564 - t568 * t751 - t712) * MDP(31) + (-t384 * t611 - t392 * t404 - t393 * t403 + t680 * t402 + (t585 - t597) * pkin(9) + t756) * MDP(32) + t764; MDP(22) * t705 + (-t458 ^ 2 + t732) * MDP(23) + t391 * MDP(24) + (-t406 + t741) * MDP(25) + t423 * MDP(26) + (-t429 * t460 + t591 + t710) * MDP(27) + (t398 * t750 + t429 * t458 - t582) * MDP(28) + (-t419 * t458 - t583 + t710 + 0.2e1 * t729) * MDP(29) + (pkin(5) * t405 - qJ(6) * t406 + (t393 - t399) * t460 + (t392 - t661) * t458) * MDP(30) + (0.2e1 * t716 - t402 * t458 + t419 * t460 + (0.2e1 * qJD(6) - t398) * t750 + t582) * MDP(31) + (t382 * qJ(6) - t383 * pkin(5) - t402 * t419 - t392 * t399 - g(1) * (-pkin(5) * t456 + qJ(6) * t457) - g(2) * (-pkin(5) * t452 + qJ(6) * t453) - t623 * t719 + t661 * t393) * MDP(32); (t705 - t423) * MDP(29) + t391 * MDP(30) + (-t732 - t743) * MDP(31) + (t583 - t729 - t762) * MDP(32);];
tau  = t1;
