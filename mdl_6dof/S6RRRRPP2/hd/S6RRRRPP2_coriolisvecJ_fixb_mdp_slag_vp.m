% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:51
% EndTime: 2019-03-09 20:52:04
% DurationCPUTime: 7.15s
% Computational Cost: add. (7974->550), mult. (19001->651), div. (0->0), fcn. (13092->6), ass. (0->234)
t565 = sin(qJ(4));
t648 = qJD(4) * t565;
t566 = sin(qJ(3));
t569 = cos(qJ(2));
t692 = cos(qJ(3));
t626 = qJD(1) * t692;
t567 = sin(qJ(2));
t650 = qJD(1) * t567;
t518 = t566 * t650 - t569 * t626;
t672 = t518 * t565;
t715 = t648 + t672;
t668 = t566 * t569;
t535 = t567 * t692 + t668;
t641 = qJD(2) + qJD(3);
t483 = t641 * t535;
t714 = t483 * qJD(1);
t515 = qJD(4) + t518;
t713 = qJ(5) * t714 + t515 * qJD(5);
t693 = -pkin(8) - pkin(7);
t545 = t693 * t567;
t537 = qJD(1) * t545;
t546 = t693 * t569;
t539 = qJD(1) * t546;
t630 = t692 * t539;
t477 = t537 * t566 - t630;
t649 = qJD(3) * t566;
t614 = pkin(2) * t649 - t477;
t520 = -qJD(1) * t668 - t567 * t626;
t568 = cos(qJ(4));
t600 = -t566 * t567 + t569 * t692;
t589 = t600 * qJD(3);
t482 = qJD(2) * t600 + t589;
t575 = t482 * qJD(1);
t618 = qJD(4) * t641;
t605 = t520 * t648 + (t575 + t618) * t568;
t490 = -t520 * t565 - t568 * t641;
t677 = t490 * t515;
t712 = t605 - t677;
t647 = qJD(4) * t568;
t671 = t518 * t568;
t711 = t647 + t671;
t642 = qJD(1) * qJD(2);
t709 = -0.2e1 * t642;
t678 = t714 * t568;
t708 = pkin(9) * (t515 * t648 - t678);
t707 = MDP(5) * (t567 ^ 2 - t569 ^ 2);
t521 = t566 * t539;
t688 = qJD(2) * pkin(2);
t525 = t537 + t688;
t474 = t692 * t525 + t521;
t457 = -pkin(3) * t641 - t474;
t593 = t568 * t520 - t565 * t641;
t586 = qJ(5) * t593 + t457;
t694 = pkin(4) + pkin(5);
t393 = -t490 * t694 + qJD(6) - t586;
t706 = (qJD(6) + t393) * t593;
t418 = t714 * pkin(3) + (-pkin(9) * t589 + (t567 * pkin(2) - pkin(9) * t600) * qJD(2)) * qJD(1);
t633 = qJD(2) * t693;
t615 = qJD(1) * t633;
t526 = t567 * t615;
t527 = t569 * t615;
t625 = qJD(3) * t692;
t426 = t525 * t625 + t526 * t692 + t566 * t527 + t539 * t649;
t557 = -pkin(2) * t569 - pkin(1);
t544 = t557 * qJD(1);
t453 = pkin(3) * t518 + pkin(9) * t520 + t544;
t475 = t566 * t525 - t630;
t458 = pkin(9) * t641 + t475;
t595 = -t565 * t418 - t568 * t426 - t453 * t647 + t458 * t648;
t371 = -t595 + t713;
t468 = t714 * pkin(4);
t616 = -t568 * t418 + t565 * t426 + t453 * t648 + t458 * t647;
t375 = -t468 + t616;
t705 = t371 * t568 + t375 * t565;
t652 = -pkin(4) * t672 + qJ(5) * t671;
t704 = t652 - t614;
t574 = t565 * t575;
t435 = -qJD(4) * t593 + t574;
t703 = t435 * qJ(6) + t490 * qJD(6);
t702 = 0.2e1 * t713;
t701 = t483 * qJ(5) - qJD(5) * t600;
t552 = qJ(6) * t648;
t700 = qJ(6) * t672 + t552;
t560 = t565 * qJ(5);
t698 = t568 * pkin(4) + t560;
t697 = t692 * t545 + t566 * t546;
t513 = pkin(4) * t648 - qJ(5) * t647 - t565 * qJD(5);
t696 = pkin(5) * t715 + t513;
t413 = t568 * t453 - t565 * t458;
t644 = qJD(5) - t413;
t495 = t566 * t545 - t546 * t692;
t538 = t567 * t633;
t540 = t569 * t633;
t440 = t495 * qJD(3) + t566 * t538 - t692 * t540;
t695 = t490 ^ 2;
t489 = t593 ^ 2;
t510 = t515 ^ 2;
t691 = pkin(4) * t520;
t690 = pkin(9) - qJ(6);
t689 = pkin(2) * qJD(3);
t687 = qJ(5) * t435;
t686 = qJ(5) * t490;
t685 = qJ(5) * t568;
t684 = qJ(6) * t518;
t414 = t565 * t453 + t568 * t458;
t683 = t414 * t515;
t682 = t605 * t565;
t681 = t435 * t568;
t555 = pkin(2) * t566 + pkin(9);
t680 = t714 * t555;
t679 = t714 * t565;
t676 = t490 * t565;
t675 = t593 * t490;
t674 = t593 * t515;
t673 = t593 * t568;
t670 = t535 * t565;
t669 = t535 * t568;
t571 = qJD(2) ^ 2;
t667 = t567 * t571;
t666 = t569 * t571;
t572 = qJD(1) ^ 2;
t665 = t569 * t572;
t664 = -qJ(6) + t555;
t470 = -pkin(3) * t520 + pkin(9) * t518;
t456 = pkin(2) * t650 + t470;
t478 = t537 * t692 + t521;
t471 = t565 * t478;
t533 = t664 * t568;
t617 = pkin(2) * t625;
t606 = t617 - qJD(6);
t636 = t694 * t520;
t663 = t471 + (-t456 + t684) * t568 + t636 - qJD(4) * t533 - t565 * t606;
t466 = t565 * t474;
t543 = t690 * t568;
t662 = t466 + (-t470 + t684) * t568 + t636 - qJD(4) * t543 + t565 * qJD(6);
t511 = t520 * qJ(5);
t654 = t565 * t456 + t568 * t478;
t406 = -t511 + t654;
t627 = t555 * t648;
t661 = t568 * t606 - t406 - t627 + t700;
t655 = t565 * t470 + t568 * t474;
t409 = -t511 + t655;
t660 = -pkin(9) * t648 - t568 * qJD(6) - t409 + t700;
t425 = t475 + t652;
t659 = -t425 + t696;
t658 = t696 - t704;
t657 = t425 - t513;
t656 = -t513 + t704;
t473 = -pkin(3) * t600 - pkin(9) * t535 + t557;
t653 = t565 * t473 + t568 * t495;
t646 = qJD(5) * t568;
t394 = -qJ(6) * t593 + t413;
t645 = qJD(5) - t394;
t639 = pkin(3) + t698;
t638 = t567 * t688;
t432 = pkin(3) * t483 - pkin(9) * t482 + t638;
t439 = qJD(3) * t697 + t692 * t538 + t566 * t540;
t635 = t565 * t432 + t568 * t439 + t473 * t647;
t634 = t565 * t439 + t473 * t648 + t495 * t647;
t419 = -qJ(5) * t600 + t653;
t632 = t565 * t692;
t631 = t568 * t692;
t628 = t535 * t647;
t624 = t567 * t642;
t623 = pkin(1) * t709;
t622 = t470 * t568 - t466;
t484 = t565 * t495;
t621 = t473 * t568 - t484;
t619 = t515 * t568;
t427 = t525 * t649 + t566 * t526 - t692 * t527 - t539 * t625;
t556 = -pkin(2) * t692 - pkin(3);
t395 = qJ(6) * t490 + t414;
t613 = -t414 * t520 + t427 * t565 + t457 * t647;
t612 = pkin(4) * t565 - t685;
t611 = t565 * t617;
t402 = -pkin(4) * t515 + t644;
t505 = t515 * qJ(5);
t403 = t505 + t414;
t610 = t402 * t568 - t403 * t565;
t609 = t457 * t518 - t680;
t608 = -qJ(6) * t482 - qJD(6) * t535;
t607 = t432 * t568 - t634;
t604 = -t565 * t694 + t685;
t529 = t556 - t698;
t377 = t435 * pkin(4) - qJ(5) * t605 + qJD(5) * t593 + t427;
t408 = t490 * pkin(4) + t586;
t603 = -t377 * t565 + t403 * t520 - t408 * t671;
t602 = -t377 * t568 - t402 * t520 + t408 * t648;
t601 = t413 * t520 - t427 * t568 + t457 * t648;
t598 = t482 * t565 + t628;
t597 = -t482 * t568 + t535 * t648;
t596 = t544 * t520 - t427;
t594 = -t495 * t648 + t635;
t592 = (-t515 * t647 - t679) * pkin(9);
t591 = -qJ(6) * t605 + t375;
t590 = -t408 * t593 + t375;
t372 = -pkin(5) * t435 - t377;
t388 = t395 + t505;
t588 = t372 * t565 + t388 * t520 + t393 * t711;
t382 = -t515 * t694 + t645;
t587 = t372 * t568 - t382 * t520 - t393 * t715;
t401 = t605 + t677;
t585 = t568 * t617 - t627;
t584 = t413 * t515 + t595;
t583 = -pkin(5) * t714 + t591;
t582 = t402 * t711 - t403 * t715 + t705;
t581 = qJD(4) * t610 + t705;
t580 = t682 - t681 + (-t673 + t676) * qJD(4);
t364 = qJD(6) * t593 + t583;
t365 = t371 + t703;
t579 = t388 * t648 + (-qJD(4) * t382 - t365) * t568 + (t388 * t518 - t364) * t565 - t382 * t671;
t578 = (t712 * t568 + (-t435 + t674) * t565) * MDP(19) + (-t593 * t619 + t682) * MDP(18) + (-t490 * t520 - t510 * t565 + t678) * MDP(21) + (t515 * t619 - t520 * t593 + t679) * MDP(20) + (t518 * t641 + t575) * MDP(13) + (-t520 * t641 - t714) * MDP(14) + (-t518 ^ 2 + t520 ^ 2) * MDP(12) + (-MDP(11) * t518 + t515 * MDP(22)) * t520;
t577 = t544 * t518 - t426;
t561 = t568 * pkin(5);
t542 = t690 * t565;
t532 = t664 * t565;
t530 = t561 + t639;
t504 = t561 - t529;
t443 = -pkin(4) * t593 + t686;
t438 = t535 * t612 - t697;
t421 = t535 * t604 + t697;
t420 = pkin(4) * t600 - t621;
t417 = t593 * t694 - t686;
t411 = -t622 + t691;
t407 = -t456 * t568 + t471 + t691;
t405 = qJ(6) * t670 + t419;
t400 = t484 + (-qJ(6) * t535 - t473) * t568 + t694 * t600;
t380 = t612 * t482 + (qJD(4) * t698 - t646) * t535 + t440;
t379 = t604 * t482 + (t646 + (-t568 * t694 - t560) * qJD(4)) * t535 - t440;
t378 = -pkin(4) * t483 - t607;
t376 = t594 + t701;
t367 = qJ(6) * t628 + (-qJD(4) * t495 - t608) * t565 + t635 + t701;
t366 = t535 * t552 - t694 * t483 + (-t432 + t608) * t568 + t634;
t1 = [(-pkin(7) * t666 + t567 * t623) * MDP(9) - MDP(7) * t667 + (-t520 * t482 + t535 * t575) * MDP(11) + (pkin(2) * t535 * t624 + t544 * t482 - t520 * t638 + t557 * t575) * MDP(17) + (t482 * MDP(13) - t483 * MDP(14) - t440 * MDP(16) - MDP(17) * t439) * t641 + (t366 * t593 + t367 * t490 - t400 * t605 + t405 * t435 + (-t382 * t568 + t388 * t565) * t482 + (-t364 * t568 + t365 * t565 + (t382 * t565 + t388 * t568) * qJD(4)) * t535) * MDP(31) + (-t376 * t490 - t378 * t593 - t419 * t435 + t420 * t605 + t610 * t482 + (-t371 * t565 + t375 * t568 + (-t402 * t565 - t403 * t568) * qJD(4)) * t535) * MDP(26) + (t593 * t597 + t605 * t669) * MDP(18) + ((-t490 * t568 + t565 * t593) * t482 + (-t682 - t681 + (t673 + t676) * qJD(4)) * t535) * MDP(19) + (-t414 * t483 + t427 * t669 - t440 * t593 - t457 * t597 - t515 * t594 - t595 * t600 - t605 * t697 - t653 * t714) * MDP(24) + (t413 * t483 + t427 * t670 - t435 * t697 + t440 * t490 + t457 * t598 + t515 * t607 + t600 * t616 + t621 * t714) * MDP(23) + (-t483 * t593 - t515 * t597 - t600 * t605 + t669 * t714) * MDP(20) + (-t365 * t600 + t367 * t515 + t372 * t669 - t379 * t593 + t388 * t483 - t393 * t597 + t405 * t714 + t421 * t605) * MDP(30) + (-t371 * t600 + t376 * t515 - t377 * t669 + t380 * t593 + t403 * t483 + t408 * t597 + t419 * t714 - t438 * t605) * MDP(27) + (t435 * t600 - t483 * t490 - t515 * t598 - t670 * t714) * MDP(21) + (t364 * t600 - t366 * t515 - t372 * t670 - t379 * t490 - t382 * t483 - t393 * t598 - t400 * t714 - t421 * t435) * MDP(29) + (t375 * t600 + t377 * t670 - t378 * t515 + t380 * t490 - t402 * t483 + t408 * t598 - t420 * t714 + t435 * t438) * MDP(25) + (t483 * t515 - t600 * t714) * MDP(22) + (-t482 * t518 + t520 * t483 - t535 * t714 + t575 * t600) * MDP(12) + (t557 * t714 + t544 * t483 + (-qJD(1) * t600 + t518) * t638) * MDP(16) + (pkin(7) * t667 + t569 * t623) * MDP(10) + 0.2e1 * t569 * MDP(4) * t624 + MDP(6) * t666 + (t371 * t419 + t375 * t420 + t376 * t403 + t377 * t438 + t378 * t402 + t380 * t408) * MDP(28) + (t364 * t400 + t365 * t405 + t366 * t382 + t367 * t388 + t372 * t421 + t379 * t393) * MDP(32) + t707 * t709; (t556 * t435 + t609 * t565 + t614 * t490 + (-t611 + t471 + (-qJD(4) * t555 - t456) * t568) * t515 + t601) * MDP(23) + (t504 * t605 + t515 * t661 + t533 * t714 + t593 * t658 + t588) * MDP(30) + (t556 * t605 + t609 * t568 - t614 * t593 + (-t585 + t654) * t515 + t613) * MDP(24) + (t377 * t529 - t402 * t407 - t403 * t406 - t656 * t408 + (t402 * t632 + t403 * t631) * t689 + t581 * t555) * MDP(28) + t572 * t707 + (t478 * t641 + (t520 * t650 - t625 * t641) * pkin(2) + t577) * MDP(17) + (t529 * t435 + (t408 * t518 - t680) * t565 - t656 * t490 + (-t555 * t647 + t407 - t611) * t515 + t602) * MDP(25) + (t364 * t532 + t365 * t533 + t372 * t504 - t382 * t663 + t388 * t661 - t393 * t658) * MDP(32) + t578 + (t477 * t641 + (-t518 * t650 - t641 * t649) * pkin(2) + t596) * MDP(16) + (t435 * t533 + t490 * t661 - t532 * t605 - t593 * t663 + t579) * MDP(31) + (-t435 * t504 + t490 * t658 + t515 * t663 - t532 * t714 + t587) * MDP(29) + (t406 * t490 + t407 * t593 + (-t490 * t631 - t593 * t632) * t689 + t580 * t555 + t582) * MDP(26) + (-t529 * t605 + (-qJD(4) * t408 + t680) * t568 - t656 * t593 + (-t406 + t585) * t515 + t603) * MDP(27) - t567 * MDP(4) * t665 + (MDP(9) * t567 * t572 + MDP(10) * t665) * pkin(1); (-pkin(3) * t435 + t457 * t672 - t475 * t490 - t515 * t622 + t592 + t601) * MDP(23) + (t515 * t660 + t530 * t605 + t543 * t714 + t593 * t659 + t588) * MDP(30) + (t364 * t542 + t365 * t543 + t372 * t530 - t382 * t662 + t388 * t660 - t393 * t659) * MDP(32) + (t474 * t641 + t577) * MDP(17) + (t408 * t672 + t411 * t515 - t435 * t639 - t490 * t657 + t592 + t602) * MDP(25) + (pkin(9) * t581 - t377 * t639 - t402 * t411 - t403 * t409 - t408 * t657) * MDP(28) + t578 + (t435 * t543 + t490 * t660 - t542 * t605 - t593 * t662 + t579) * MDP(31) + (-t435 * t530 + t490 * t659 + t515 * t662 - t542 * t714 + t587) * MDP(29) + (pkin(9) * t580 + t409 * t490 + t411 * t593 + t582) * MDP(26) + (-t408 * t647 - t409 * t515 - t593 * t657 + t605 * t639 + t603 - t708) * MDP(27) + (-pkin(3) * t605 + t457 * t671 + t475 * t593 + t515 * t655 + t613 + t708) * MDP(24) + (t475 * t641 + t596) * MDP(16); -MDP(18) * t675 + (t489 - t695) * MDP(19) + t401 * MDP(20) + (-t435 - t674) * MDP(21) + t714 * MDP(22) + (t457 * t593 - t616 + t683) * MDP(23) + (t457 * t490 + t584) * MDP(24) + (-t443 * t490 + t468 - t590 + t683) * MDP(25) + (-pkin(4) * t605 - t687 - (t403 - t414) * t593 + (t402 - t644) * t490) * MDP(26) + (-t408 * t490 - t443 * t593 - t584 + t702) * MDP(27) + (-pkin(4) * t375 + qJ(5) * t371 - t402 * t414 + t403 * t644 - t408 * t443) * MDP(28) + (t395 * t515 + t417 * t490 - t706 + (pkin(5) + t694) * t714 - t591) * MDP(29) + (t393 * t490 - t394 * t515 + t417 * t593 - t595 + t702 + t703) * MDP(30) + (t687 + t605 * t694 - (-t388 + t395) * t593 + (-t382 + t645) * t490) * MDP(31) + (qJ(5) * t365 - t364 * t694 - t382 * t395 + t388 * t645 - t393 * t417) * MDP(32); (-t403 * t515 + t590) * MDP(28) + (-t388 * t515 + t583 + t706) * MDP(32) + (MDP(25) + MDP(29)) * (-t714 - t675) + (MDP(26) - MDP(31)) * t401 + (MDP(27) + MDP(30)) * (-t510 - t489); t712 * MDP(30) + (-t489 - t695) * MDP(31) + (-t382 * t593 - t388 * t490 + t372) * MDP(32) + (t520 * t647 - t565 * t618 - t574 + t674) * MDP(29);];
tauc  = t1;
