% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:28
% EndTime: 2019-03-09 17:01:45
% DurationCPUTime: 10.98s
% Computational Cost: add. (11199->538), mult. (29429->735), div. (0->0), fcn. (23150->10), ass. (0->233)
t605 = sin(pkin(6));
t613 = cos(qJ(2));
t679 = qJD(1) * t613;
t593 = t605 * t679;
t632 = t593 - qJD(3);
t609 = sin(qJ(3));
t610 = sin(qJ(2));
t612 = cos(qJ(3));
t680 = qJD(1) * t610;
t661 = t605 * t680;
t607 = cos(pkin(6));
t681 = qJD(1) * t607;
t668 = pkin(1) * t681;
t552 = -pkin(8) * t661 + t613 * t668;
t628 = (pkin(2) * t610 - pkin(9) * t613) * t605;
t553 = qJD(1) * t628;
t647 = -t552 * t609 + t612 * t553;
t713 = -qJ(4) - pkin(9);
t655 = qJD(3) * t713;
t682 = qJD(1) * t605;
t702 = t612 * t613;
t734 = (pkin(3) * t610 - qJ(4) * t702) * t682 + t647 + qJD(4) * t609 - t612 * t655;
t638 = t609 * t593;
t685 = t612 * t552 + t609 * t553;
t733 = qJ(4) * t638 + qJD(4) * t612 + t609 * t655 - t685;
t604 = sin(pkin(11));
t606 = cos(pkin(11));
t568 = t604 * t612 + t606 * t609;
t684 = t632 * t568;
t567 = t604 * t609 - t606 * t612;
t517 = t567 * t593;
t561 = t567 * qJD(3);
t732 = t561 - t517;
t690 = -t604 * t734 + t733 * t606;
t590 = t610 * t668;
t555 = pkin(8) * t593 + t590;
t676 = qJD(3) * t609;
t731 = -t555 + (-t638 + t676) * pkin(3);
t729 = pkin(10) * t661 - t690;
t730 = -t684 * pkin(4) + t732 * pkin(10) + t731;
t691 = t733 * t604 + t606 * t734;
t669 = qJD(1) * qJD(2);
t656 = t605 * t669;
t636 = t613 * t656;
t546 = pkin(8) * t636 + qJD(2) * t590;
t640 = qJD(3) * t632;
t728 = -pkin(9) * t640 + t546;
t594 = qJD(2) + t681;
t635 = qJD(3) * t661;
t675 = qJD(3) * t612;
t509 = t594 * t675 - t609 * t635 + t612 * t636;
t537 = t594 * t609 + t612 * t661;
t522 = pkin(9) * t594 + t555;
t550 = (-pkin(2) * t613 - pkin(9) * t610 - pkin(1)) * t605;
t529 = qJD(1) * t550;
t485 = t522 * t612 + t529 * t609;
t554 = qJD(2) * t628;
t544 = qJD(1) * t554;
t706 = t605 * t610;
t595 = pkin(8) * t706;
t714 = pkin(1) * t613;
t556 = (t607 * t714 - t595) * qJD(2);
t545 = qJD(1) * t556;
t649 = t612 * t544 - t609 * t545;
t615 = -qJD(3) * t485 + t649;
t637 = t610 * t656;
t413 = pkin(3) * t637 - qJ(4) * t509 - qJD(4) * t537 + t615;
t535 = -t612 * t594 + t609 * t661;
t627 = t522 * t676 - t529 * t675 - t609 * t544 - t612 * t545;
t664 = t594 * t676 + t609 * t636 + t612 * t635;
t419 = -qJ(4) * t664 - t535 * qJD(4) - t627;
t388 = t604 * t413 + t606 * t419;
t385 = pkin(10) * t637 + t388;
t484 = -t609 * t522 + t612 * t529;
t463 = -t537 * qJ(4) + t484;
t453 = -pkin(3) * t632 + t463;
t464 = -qJ(4) * t535 + t485;
t704 = t606 * t464;
t412 = t604 * t453 + t704;
t406 = -pkin(10) * t632 + t412;
t521 = -pkin(2) * t594 - t552;
t490 = pkin(3) * t535 + qJD(4) + t521;
t630 = -t535 * t604 + t606 * t537;
t650 = -t606 * t535 - t537 * t604;
t425 = -pkin(4) * t650 - pkin(10) * t630 + t490;
t608 = sin(qJ(5));
t611 = cos(qJ(5));
t390 = t406 * t611 + t425 * t608;
t465 = t509 * t604 + t606 * t664;
t481 = pkin(3) * t664 + t546;
t621 = t606 * t509 - t604 * t664;
t403 = t465 * pkin(4) - pkin(10) * t621 + t481;
t376 = -qJD(5) * t390 - t385 * t608 + t611 * t403;
t645 = t611 * t632;
t674 = qJD(5) * t608;
t421 = qJD(5) * t645 - t608 * t637 - t611 * t621 + t630 * t674;
t473 = -t608 * t632 + t611 * t630;
t371 = pkin(5) * t465 + qJ(6) * t421 - qJD(6) * t473 + t376;
t471 = t608 * t630 + t645;
t380 = -qJ(6) * t471 + t390;
t717 = qJD(5) - t650;
t727 = t380 * t717 + t371;
t673 = qJD(5) * t611;
t375 = t611 * t385 + t608 * t403 - t406 * t674 + t425 * t673;
t618 = t608 * t621 - t611 * t637;
t422 = qJD(5) * t473 + t618;
t372 = -qJ(6) * t422 - qJD(6) * t471 + t375;
t389 = -t406 * t608 + t611 * t425;
t379 = -qJ(6) * t473 + t389;
t378 = pkin(5) * t717 + t379;
t726 = -t378 * t717 + t372;
t601 = t605 ^ 2;
t725 = -0.2e1 * t601 * t669;
t723 = MDP(5) * (t610 ^ 2 - t613 ^ 2);
t722 = t535 * t632;
t721 = t537 * t632;
t720 = t608 * t717;
t563 = t607 * t609 + t612 * t706;
t705 = t605 * t613;
t715 = pkin(1) * t610;
t549 = pkin(8) * t705 + (pkin(9) + t715) * t607;
t648 = -t549 * t609 + t612 * t550;
t468 = -pkin(3) * t705 - qJ(4) * t563 + t648;
t562 = -t607 * t612 + t609 * t706;
t688 = t612 * t549 + t609 * t550;
t479 = -qJ(4) * t562 + t688;
t434 = t604 * t468 + t606 * t479;
t429 = -pkin(10) * t705 + t434;
t504 = t606 * t562 + t563 * t604;
t505 = -t562 * t604 + t563 * t606;
t548 = t595 + (-pkin(2) - t714) * t607;
t619 = pkin(3) * t562 + t548;
t449 = pkin(4) * t504 - pkin(10) * t505 + t619;
t695 = t611 * t429 + t608 * t449;
t719 = t730 * t611;
t692 = pkin(4) * t661 + t691;
t662 = -pkin(3) * t612 - pkin(2);
t510 = pkin(4) * t567 - pkin(10) * t568 + t662;
t718 = t510 * t673 + t730 * t608 - t729 * t611;
t716 = t473 ^ 2;
t712 = t421 * t608;
t711 = t471 * t650;
t710 = t650 * t608;
t709 = t568 * t608;
t708 = t568 * t611;
t614 = qJD(1) ^ 2;
t707 = t601 * t614;
t458 = t604 * t464;
t703 = t608 * t465;
t461 = t611 * t465;
t586 = t713 * t609;
t587 = t713 * t612;
t520 = t586 * t604 - t587 * t606;
t512 = t611 * t520;
t599 = pkin(3) * t604 + pkin(10);
t701 = qJ(6) + t599;
t700 = t378 - t379;
t502 = -t517 * t611 + t608 * t661;
t629 = qJ(6) * t561 - qJD(6) * t568;
t699 = qJ(6) * t502 - t512 * qJD(5) + t629 * t611 + t719 + ((qJ(6) * t568 - t510) * qJD(5) + t729) * t608 - t684 * pkin(5);
t501 = -t517 * t608 - t611 * t661;
t658 = t568 * t673;
t698 = (-qJD(5) * t520 + t629) * t608 + t718 + (t501 - t658) * qJ(6);
t418 = t463 * t606 - t458;
t443 = pkin(3) * t537 + pkin(4) * t630 - pkin(10) * t650;
t697 = t611 * t418 + t608 * t443;
t696 = -t608 * t422 - t471 * t673;
t659 = qJD(2) * t705;
t514 = -qJD(3) * t562 + t612 * t659;
t616 = -qJD(3) * t688 + t612 * t554 - t556 * t609;
t678 = qJD(2) * t610;
t660 = t605 * t678;
t430 = pkin(3) * t660 - qJ(4) * t514 - qJD(4) * t563 + t616;
t513 = qJD(3) * t563 + t609 * t659;
t626 = -t549 * t676 + t550 * t675 + t609 * t554 + t612 * t556;
t435 = -qJ(4) * t513 - qJD(4) * t562 + t626;
t397 = t604 * t430 + t606 * t435;
t693 = t710 * t717 + t461;
t689 = t608 * t510 + t512;
t646 = qJD(5) * t701;
t687 = qJ(6) * t710 + qJD(6) * t611 - t608 * t646 - t697;
t442 = t611 * t443;
t686 = -pkin(5) * t630 - t442 + (qJ(6) * t650 - t646) * t611 + (-qJD(6) + t418) * t608;
t557 = t607 * pkin(1) * t678 + pkin(8) * t659;
t677 = qJD(2) * t612;
t672 = t521 * qJD(3);
t666 = t608 * t705;
t600 = -pkin(3) * t606 - pkin(4);
t387 = t413 * t606 - t604 * t419;
t654 = -t429 * t608 + t611 * t449;
t396 = t430 * t606 - t604 * t435;
t411 = t606 * t453 - t458;
t417 = t463 * t604 + t704;
t433 = t468 * t606 - t604 * t479;
t653 = t561 * t608 + t501;
t652 = t561 * t611 + t502;
t651 = t611 * t510 - t520 * t608;
t519 = -t606 * t586 - t587 * t604;
t644 = t613 * t632;
t643 = t717 * t611;
t642 = t632 * t605;
t639 = t601 * t610 * t613 * MDP(4);
t634 = t513 * pkin(3) + t557;
t633 = pkin(1) * t725;
t428 = pkin(4) * t705 - t433;
t631 = t513 * t604 - t514 * t606;
t487 = t505 * t608 + t611 * t705;
t394 = pkin(10) * t660 + t397;
t478 = t606 * t513 + t514 * t604;
t409 = t478 * pkin(4) + pkin(10) * t631 + t634;
t625 = t611 * t394 + t608 * t409 - t429 * t674 + t449 * t673;
t405 = pkin(4) * t632 - t411;
t624 = t405 * t717 - t599 * t465;
t623 = -t653 + t658;
t622 = -t568 * t674 - t652;
t393 = -pkin(4) * t660 - t396;
t384 = -pkin(4) * t637 - t387;
t617 = -qJD(5) * t695 - t394 * t608 + t611 * t409;
t377 = pkin(5) * t422 + t384;
t565 = t701 * t611;
t564 = t701 * t608;
t488 = t505 * t611 - t666;
t470 = t471 ^ 2;
t456 = -qJ(6) * t709 + t689;
t451 = pkin(5) * t567 - qJ(6) * t708 + t651;
t438 = -qJD(5) * t666 + t505 * t673 - t608 * t631 - t611 * t660;
t437 = qJD(5) * t487 - t608 * t660 + t611 * t631;
t399 = t471 * pkin(5) + qJD(6) + t405;
t391 = -qJ(6) * t487 + t695;
t382 = pkin(5) * t504 - qJ(6) * t488 + t654;
t374 = -qJ(6) * t438 - qJD(6) * t487 + t625;
t373 = pkin(5) * t478 + qJ(6) * t437 - qJD(6) * t488 + t617;
t1 = [0.2e1 * t639 * t669 + (t376 * t504 + t384 * t487 + t389 * t478 + t393 * t471 + t405 * t438 + t428 * t422 + t465 * t654 + t617 * t717) * MDP(25) + (-t375 * t504 + t384 * t488 - t390 * t478 + t393 * t473 - t405 * t437 - t428 * t421 - t465 * t695 - t625 * t717) * MDP(26) + (-t371 * t488 - t372 * t487 - t373 * t473 - t374 * t471 + t378 * t437 - t380 * t438 + t382 * t421 - t391 * t422) * MDP(27) + (t372 * t391 + t380 * t374 + t371 * t382 + t378 * t373 + t377 * (pkin(5) * t487 + t428) + t399 * (pkin(5) * t438 + t393)) * MDP(28) + t723 * t725 + (-t546 * t607 - t557 * t594 + t610 * t633) * MDP(9) + (-t545 * t607 - t556 * t594 + t613 * t633) * MDP(10) + (t509 * t563 + t514 * t537) * MDP(11) + (-t509 * t562 - t537 * t513 - t514 * t535 - t563 * t664) * MDP(12) + (-t514 * t632 + (-t509 * t613 + (qJD(1) * t563 + t537) * t678) * t605) * MDP(13) + (t513 * t632 + (t664 * t613 + (-qJD(1) * t562 - t535) * t678) * t605) * MDP(14) + (-t601 * t679 - t642) * MDP(15) * t678 + (-t616 * t632 + t557 * t535 + t548 * t664 + t546 * t562 + t521 * t513 + (-t615 * t613 + (qJD(1) * t648 + t484) * t678) * t605) * MDP(16) + (t626 * t632 + t557 * t537 + t548 * t509 + t546 * t563 + t521 * t514 + (-t627 * t613 + (-qJD(1) * t688 - t485) * t678) * t605) * MDP(17) + (-t387 * t505 - t388 * t504 - t396 * t630 + t397 * t650 + t411 * t631 - t412 * t478 - t433 * t621 - t434 * t465) * MDP(18) + (t387 * t433 + t388 * t434 + t411 * t396 + t412 * t397 + t481 * t619 + t490 * t634) * MDP(19) + (-t421 * t488 - t437 * t473) * MDP(20) + (t421 * t487 - t422 * t488 + t437 * t471 - t438 * t473) * MDP(21) + (-t421 * t504 - t437 * t717 + t465 * t488 + t473 * t478) * MDP(22) + (-t422 * t504 - t438 * t717 - t465 * t487 - t471 * t478) * MDP(23) + (t465 * t504 + t478 * t717) * MDP(24) + (MDP(6) * t659 - MDP(7) * t660) * (t594 + t681); t707 * t723 + (t555 * t594 + t707 * t715 - t546) * MDP(9) + (pkin(8) * t637 + t552 * t594 + (-t607 * t669 + t707) * t714) * MDP(10) + (t509 * t609 - t612 * t721) * MDP(11) + ((t509 + t722) * t612 + (-t664 + t721) * t609) * MDP(12) + (-t612 * t640 + (t612 * t644 + (qJD(2) * t609 - t537) * t610) * t682) * MDP(13) + (t609 * t640 + (-t609 * t644 + (t535 + t677) * t610) * t682) * MDP(14) + (-pkin(2) * t664 + t609 * t672 + t647 * t632 - t555 * t535 - t728 * t612 + (-t484 * t610 + (-pkin(9) * t678 - t521 * t613) * t609) * t682) * MDP(16) + (-pkin(2) * t509 + t612 * t672 - t685 * t632 - t555 * t537 + t728 * t609 + (-t521 * t702 + (-pkin(9) * t677 + t485) * t610) * t682) * MDP(17) + (-t387 * t568 - t388 * t567 + t732 * t411 + t684 * t412 - t520 * t465 + t519 * t621 + t691 * t630 + t690 * t650) * MDP(18) + (-t387 * t519 + t388 * t520 - t691 * t411 + t690 * t412 + t481 * t662 + t490 * t731) * MDP(19) + (-t421 * t708 + t473 * t622) * MDP(20) + (t653 * t473 + t652 * t471 + (t712 - t422 * t611 + (t471 * t608 - t473 * t611) * qJD(5)) * t568) * MDP(21) + (-t421 * t567 + t461 * t568 - t473 * t684 + t622 * t717) * MDP(22) + (-t422 * t567 + t471 * t684 - t568 * t703 - t623 * t717) * MDP(23) + (t465 * t567 - t684 * t717) * MDP(24) + (t651 * t465 + t376 * t567 + t519 * t422 + t384 * t709 + (-t520 * t673 + (-qJD(5) * t510 + t729) * t608 + t719) * t717 + t692 * t471 - t684 * t389 + t623 * t405) * MDP(25) + (-t689 * t465 - t375 * t567 - t519 * t421 + t384 * t708 + (t520 * t674 - t718) * t717 + t692 * t473 + t684 * t390 + t622 * t405) * MDP(26) + (t421 * t451 - t422 * t456 - t699 * t473 - t698 * t471 + t653 * t380 + t652 * t378 + (-t371 * t611 - t372 * t608 + (t378 * t608 - t380 * t611) * qJD(5)) * t568) * MDP(27) + (t372 * t456 + t371 * t451 + t377 * (pkin(5) * t709 + t519) + (pkin(5) * t623 + t692) * t399 + t698 * t380 + t699 * t378) * MDP(28) - t614 * t639 + MDP(15) * t642 * t680 + (MDP(6) * t593 - MDP(7) * t661) * (qJD(2) - t594); t537 * t535 * MDP(11) + (-t535 ^ 2 + t537 ^ 2) * MDP(12) + (t509 - t722) * MDP(13) + (-t664 - t721) * MDP(14) + MDP(15) * t637 + (-t485 * t593 - t521 * t537 + t649) * MDP(16) + (-t484 * t632 + t521 * t535 + t627) * MDP(17) + ((-t604 * t465 - t606 * t621) * pkin(3) + (t411 - t418) * t650 + (t412 - t417) * t630) * MDP(18) + (t411 * t417 - t412 * t418 + (t387 * t606 + t388 * t604 - t490 * t537) * pkin(3)) * MDP(19) + (t473 * t643 - t712) * MDP(20) + ((-t421 + t711) * t611 - t473 * t720 + t696) * MDP(21) + (-t473 * t630 + t643 * t717 + t703) * MDP(22) + (t471 * t630 - t674 * t717 + t693) * MDP(23) - t717 * t630 * MDP(24) + (-t384 * t611 - t389 * t630 - t417 * t471 + t600 * t422 + (-t599 * t673 - t442) * t717 + (t418 * t717 + t624) * t608) * MDP(25) + (t384 * t608 + t390 * t630 - t417 * t473 - t600 * t421 + (t599 * t674 + t697) * t717 + t624 * t611) * MDP(26) + (-t421 * t564 - t422 * t565 - t687 * t471 - t686 * t473 - t608 * t727 + t611 * t726) * MDP(27) + (t372 * t565 - t371 * t564 + t377 * (-pkin(5) * t611 + t600) + (pkin(5) * t720 - t417) * t399 + t687 * t380 + t686 * t378) * MDP(28); -t650 ^ 2 * MDP(18) + (-t412 * t650 + t481) * MDP(19) + t693 * MDP(25) + t696 * MDP(27) + (-MDP(18) * t630 + MDP(19) * t411 - MDP(25) * t471 - MDP(26) * t473 - MDP(28) * t399) * t630 + (-t465 * MDP(26) + t726 * MDP(28) + (-MDP(25) * qJD(5) + MDP(27) * t473) * t717) * t608 + ((t421 + t711) * MDP(27) + t727 * MDP(28) - t717 ^ 2 * MDP(26)) * t611; t473 * t471 * MDP(20) + (-t470 + t716) * MDP(21) + (t471 * t717 - t421) * MDP(22) + (-t618 + (-qJD(5) + t717) * t473) * MDP(23) + t465 * MDP(24) + (t390 * t717 - t405 * t473 + t376) * MDP(25) + (t389 * t717 + t405 * t471 - t375) * MDP(26) + (pkin(5) * t421 - t471 * t700) * MDP(27) + (t700 * t380 + (-t399 * t473 + t371) * pkin(5)) * MDP(28); (-t470 - t716) * MDP(27) + (t378 * t473 + t380 * t471 + t377) * MDP(28);];
tauc  = t1;
