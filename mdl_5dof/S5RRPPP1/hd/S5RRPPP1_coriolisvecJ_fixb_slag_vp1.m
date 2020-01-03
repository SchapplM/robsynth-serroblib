% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:24:25
% DurationCPUTime: 49.60s
% Computational Cost: add. (16023->801), mult. (47955->1021), div. (0->0), fcn. (50241->8), ass. (0->361)
t802 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t801 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t800 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t799 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t798 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t797 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t488 = sin(pkin(5));
t489 = sin(qJ(2));
t491 = cos(qJ(2));
t710 = cos(pkin(8));
t711 = cos(pkin(5));
t561 = t711 * t710;
t709 = sin(pkin(8));
t402 = t489 * t561 + t491 * t709;
t492 = cos(qJ(1));
t508 = t402 * t492;
t490 = sin(qJ(1));
t595 = t490 * t710;
t337 = -t488 * t595 + t508;
t592 = t492 * t710;
t594 = t490 * t709;
t509 = t488 * t594 + t491 * t592;
t560 = t711 * t709;
t540 = t489 * t560;
t338 = -t492 * t540 + t509;
t687 = t489 * t492;
t628 = t488 * t687;
t405 = t490 * t711 + t628;
t335 = t402 * t490 + t488 * t592;
t403 = t491 * t710 - t540;
t591 = t492 * t709;
t569 = t488 * t591;
t336 = t403 * t490 - t569;
t688 = t489 * t490;
t404 = t488 * t688 - t492 * t711;
t751 = -t800 * t335 + t802 * t336 - t799 * t404;
t774 = t798 * t335 - t799 * t336 + t801 * t404;
t811 = -t797 * t335 + t800 * t336 - t798 * t404;
t822 = t337 * t811 - t338 * t751 - t405 * t774;
t784 = -t335 * t811 + t336 * t751 + t404 * t774;
t750 = -t800 * t337 + t802 * t338 - t799 * t405;
t748 = t797 * t337 - t800 * t338 + t798 * t405;
t746 = t798 * t337 - t799 * t338 + t801 * t405;
t539 = t491 * t561;
t400 = t489 * t709 - t539;
t401 = t489 * t710 + t491 * t560;
t689 = t488 * t491;
t738 = t798 * t400 - t799 * t401 - t801 * t689;
t737 = t797 * t400 - t800 * t401 - t798 * t689;
t735 = -t800 * t400 + t802 * t401 + t799 * t689;
t690 = t488 * t489;
t423 = pkin(2) * t491 + qJ(3) * t690;
t527 = -pkin(1) - t423;
t790 = t490 * t527;
t482 = Icges(3,4) * t491;
t550 = -Icges(3,2) * t489 + t482;
t379 = Icges(3,6) * t490 + t492 * t550;
t707 = Icges(3,4) * t489;
t450 = Icges(3,1) * t491 - t707;
t381 = Icges(3,5) * t490 + t450 * t492;
t686 = t490 * t491;
t358 = t381 * t686;
t446 = Icges(3,5) * t491 - Icges(3,6) * t489;
t377 = Icges(3,3) * t490 + t446 * t492;
t583 = t377 * t492 - t358;
t185 = -t379 * t688 - t583;
t783 = t748 * t335 + t750 * t336 + t746 * t404;
t821 = t185 + t783;
t685 = t491 * t492;
t659 = t490 * t377 + t381 * t685;
t187 = -t379 * t687 + t659;
t782 = t748 * t337 + t750 * t338 + t746 * t405;
t820 = t187 + t782;
t818 = t822 * t492;
t817 = t784 * t492;
t447 = Icges(3,2) * t491 + t707;
t449 = Icges(3,1) * t489 + t482;
t546 = t447 * t489 - t449 * t491;
t445 = Icges(3,5) * t489 + Icges(3,6) * t491;
t691 = t445 * t492;
t816 = -t737 * t335 - t735 * t336 - t738 * t404 + t490 * t546 + t691;
t692 = t445 * t490;
t815 = t737 * t337 + t735 * t338 + t738 * t405 - t492 * t546 + t692;
t792 = rSges(6,1) + pkin(4);
t320 = qJD(5) * t336;
t791 = rSges(6,3) + qJ(5);
t679 = t337 * rSges(6,2) + t791 * t338 + t405 * t792;
t640 = qJD(2) * t490;
t670 = rSges(6,2) * t400 + t791 * t401 - t689 * t792;
t785 = t640 * t670;
t208 = t338 * pkin(3) + qJ(4) * t337;
t315 = pkin(3) * t401 + qJ(4) * t400;
t636 = qJD(4) * t335;
t475 = pkin(2) * t685;
t589 = t711 * qJ(3);
t362 = qJ(3) * t628 + t490 * t589 + t475;
t726 = pkin(2) * t489;
t538 = qJ(3) * t689 - t726;
t647 = t492 * pkin(1) + t490 * pkin(7);
t786 = qJD(1) * t647;
t803 = -qJD(1) * t362 - qJD(3) * t404 - t538 * t640 - t786;
t788 = -qJD(1) * t208 + t315 * t640 - t636 + t803;
t40 = t679 * qJD(1) + t320 - t785 - t788;
t646 = qJD(1) * t335;
t755 = t400 * qJD(2);
t220 = t492 * t755 + t646;
t373 = t401 * t492;
t521 = qJD(1) * t540;
t564 = qJD(1) * t595;
t221 = -qJD(1) * t569 + qJD(2) * t373 - t490 * t521 + t491 * t564;
t639 = qJD(2) * t492;
t609 = t491 * t639;
t441 = t488 * t609;
t581 = -qJD(1) * t404 + t441;
t780 = -t220 * t797 + t221 * t800 + t581 * t798;
t222 = qJD(1) * t508 - t488 * t564 - t490 * t755;
t537 = qJD(2) * t560;
t587 = qJD(2) * t710;
t223 = qJD(1) * t509 - t492 * t521 - t537 * t686 - t587 * t688;
t630 = t488 * t686;
t756 = qJD(1) * t405 + qJD(2) * t630;
t779 = -t222 * t797 + t223 * t800 - t756 * t798;
t778 = t800 * t220 - t221 * t802 - t799 * t581;
t777 = t800 * t222 - t223 * t802 + t799 * t756;
t776 = -t220 * t798 + t221 * t799 + t581 * t801;
t812 = t222 * t798 - t223 * t799 + t756 * t801;
t384 = t402 * qJD(2);
t385 = -t489 * t537 + t491 * t587;
t641 = qJD(2) * t489;
t611 = t488 * t641;
t810 = t384 * t798 - t385 * t799 + t611 * t801;
t809 = t384 * t797 - t385 * t800 + t611 * t798;
t808 = -t800 * t384 + t385 * t802 - t799 * t611;
t524 = qJD(2) * t445;
t472 = Icges(3,4) * t688;
t706 = Icges(3,5) * t492;
t380 = Icges(3,1) * t686 - t472 - t706;
t704 = Icges(3,6) * t492;
t378 = Icges(3,4) * t686 - Icges(3,2) * t688 - t704;
t694 = t378 * t489;
t548 = -t380 * t491 + t694;
t806 = (-t490 * t524 + (t377 + t548) * qJD(1)) * t492;
t469 = t492 * t589;
t361 = t423 * t490 - t469;
t391 = qJD(3) * t405;
t486 = t492 * pkin(7);
t459 = pkin(1) * t490 - t486;
t432 = qJD(1) * t459;
t642 = qJD(1) * t492;
t478 = pkin(7) * t642;
t637 = qJD(3) * t489;
t467 = t488 * t637;
t586 = qJD(3) * t711;
t571 = qJ(3) * t441 + qJD(1) * t469 + t492 * t467 + t490 * t586;
t610 = t489 * t639;
t805 = -pkin(2) * t610 + qJD(1) * t361 - t391 + t432 + t478 + t571;
t451 = rSges(3,1) * t489 + rSges(3,2) * t491;
t804 = t451 * t640 - t786;
t796 = t815 * qJD(1);
t703 = Icges(3,3) * t492;
t376 = Icges(3,5) * t686 - Icges(3,6) * t688 - t703;
t660 = -t490 * t376 - t380 * t685;
t186 = -t378 * t687 - t660;
t795 = (-t186 * t492 + t820 * t490 + t818) * qJD(2);
t184 = -t376 * t492 - t490 * t548;
t794 = (-t184 * t492 + t821 * t490 - t817) * qJD(2);
t793 = t816 * qJD(1);
t324 = t335 * qJ(4);
t206 = pkin(3) * t336 + t324;
t323 = qJD(4) * t337;
t87 = -t221 * pkin(3) - qJ(4) * t220 + t323;
t789 = qJD(1) * t206 - t323 + t805 + t87;
t787 = -rSges(6,2) * t335 - t404 * t792;
t773 = -t793 + t794;
t772 = t795 + t796;
t525 = qJD(2) * t447;
t526 = qJD(2) * t449;
t771 = t548 * qJD(2) - (qJD(1) * t379 - t490 * t525) * t491 - (qJD(1) * t381 - t490 * t526) * t489 + (t491 * t812 - t641 * t774) * t488 + t777 * t401 + t779 * t400 - t751 * t385 + t811 * t384;
t693 = t379 * t489;
t547 = -t381 * t491 + t693;
t770 = -qJD(2) * t547 + (-t492 * t525 + (-t490 * t550 + t704) * qJD(1)) * t491 + (-t492 * t526 + (-t450 * t490 + t706) * qJD(1)) * t489 + (-t491 * t776 + t641 * t746) * t488 + t778 * t401 + t780 * t400 + t750 * t385 + t748 * t384;
t420 = t550 * qJD(2);
t421 = t450 * qJD(2);
t496 = qJD(1) * t445 - t420 * t489 + t421 * t491 + (-t447 * t491 - t449 * t489) * qJD(2);
t730 = t546 * qJD(1) + t446 * qJD(2);
t769 = -t220 * t737 - t221 * t735 + t337 * t809 + t338 * t808 + t405 * t810 + t490 * t730 + t496 * t492 + t581 * t738;
t768 = t222 * t737 + t223 * t735 + t335 * t809 + t336 * t808 + t404 * t810 + t496 * t490 - t492 * t730 + t738 * t756;
t767 = t186 - t822;
t766 = t378 * t491 + t380 * t489 - t400 * t811 + t401 * t751 - t689 * t774;
t765 = t379 * t491 + t381 * t489 + t400 * t748 + t401 * t750 - t689 * t746;
t763 = t469 + t486 + t790;
t650 = t492 * t586 + t640 * t726;
t762 = t527 * t642 + t650 + ((-t589 - pkin(7)) * qJD(1) + (-qJ(3) * qJD(2) * t491 - t637) * t488) * t490;
t761 = 0.2e1 * qJD(2);
t681 = t336 * t791 - t787;
t760 = -(t335 * t490 + t337 * t492) * qJD(2) + t384;
t759 = t400 * t639 - t220 + t646;
t758 = -qJD(1) * t337 + t400 * t640 + t222;
t754 = -t222 * rSges(6,2) - t756 * t792 - t320;
t177 = t338 * rSges(4,1) - t337 * rSges(4,2) + t405 * rSges(4,3);
t752 = -qJD(1) * t177 + t803;
t643 = qJD(1) * t490;
t370 = t489 * t594 - t490 * t539;
t371 = t401 * t490;
t745 = -t370 * t800 + t371 * t802 + t630 * t799;
t372 = t489 * t591 - t492 * t539;
t629 = t488 * t685;
t744 = t372 * t800 - t373 * t802 - t629 * t799;
t743 = t370 * t797 - t371 * t800 - t630 * t798;
t742 = -t372 * t797 + t373 * t800 + t629 * t798;
t741 = -t370 * t798 + t371 * t799 + t630 * t801;
t740 = -t372 * t798 + t373 * t799 + t629 * t801;
t739 = t402 * t798 - t403 * t799 + t690 * t801;
t736 = t402 * t797 - t403 * t800 + t690 * t798;
t734 = -t402 * t800 + t403 * t802 - t690 * t799;
t183 = t405 * rSges(5,1) - t338 * rSges(5,2) + t337 * rSges(5,3);
t733 = -qJD(1) * t183 + t788;
t732 = -t492 * t524 + (-t446 * t490 + t547 + t703) * qJD(1);
t720 = rSges(3,1) * t491;
t300 = -rSges(5,1) * t689 - rSges(5,2) * t401 + rSges(5,3) * t400;
t179 = rSges(5,1) * t404 - rSges(5,2) * t336 + rSges(5,3) * t335;
t666 = -t315 + t538;
t621 = -t300 + t666;
t565 = t621 * t492;
t658 = -t361 - t459;
t624 = -t206 + t658;
t664 = t323 + t391;
t44 = qJD(2) * t565 + (-t179 + t624) * qJD(1) + t664;
t719 = t44 * t300;
t483 = t490 * rSges(3,3);
t301 = rSges(4,1) * t401 - rSges(4,2) * t400 - rSges(4,3) * t689;
t175 = rSges(4,1) * t336 - rSges(4,2) * t335 + rSges(4,3) * t404;
t669 = -t301 + t538;
t585 = t492 * t669;
t78 = t391 + qJD(2) * t585 + (-t175 + t658) * qJD(1);
t712 = t78 * t301;
t648 = rSges(3,2) * t688 + t492 * rSges(3,3);
t382 = rSges(3,1) * t686 - t648;
t612 = t451 * t639;
t271 = -t612 + (-t382 - t459) * qJD(1);
t698 = t271 * t490;
t697 = t271 * t492;
t383 = rSges(3,1) * t685 - rSges(3,2) * t687 + t483;
t272 = qJD(1) * t383 - t804;
t413 = t451 * t492;
t696 = t272 * t413;
t321 = qJD(5) * t338;
t684 = -t220 * rSges(6,2) - t791 * t221 + t581 * t792 + t321;
t683 = t223 * t791 - t754;
t682 = -t177 - t362;
t205 = qJ(3) * t756 + qJD(1) * t475 + t490 * t467 - t650;
t678 = -t205 - t786;
t677 = -t208 - t362;
t607 = qJD(3) * t689;
t365 = qJD(2) * t423 - t607;
t387 = qJD(4) * t400;
t676 = -pkin(3) * t385 - qJ(4) * t384 - t365 - t387;
t675 = t222 * qJ(4) + t636;
t674 = -rSges(6,2) * t370 - t791 * t371 + t630 * t792;
t673 = -rSges(6,2) * t372 - t791 * t373 + t629 * t792;
t672 = -rSges(4,1) * t385 + rSges(4,2) * t384 - rSges(4,3) * t611 - t365;
t273 = -pkin(3) * t371 - qJ(4) * t370;
t392 = t538 * t490;
t671 = -t273 - t392;
t668 = -rSges(4,1) * t403 + rSges(4,2) * t402 - rSges(4,3) * t690 - t423;
t399 = t538 * t643;
t667 = t315 * t643 - t399;
t665 = -pkin(3) * t403 - qJ(4) * t402 - t423;
t633 = qJD(1) * qJD(2);
t605 = t490 * t633;
t663 = qJD(3) * t581 - t538 * t605;
t661 = t490 * t361 + t492 * t362;
t443 = t492 * t607;
t657 = -qJD(4) * t372 + t443;
t393 = t538 * t492;
t656 = qJD(1) * t393 + t490 * t607;
t655 = -Icges(3,2) * t686 + t380 - t472;
t654 = -t447 * t492 + t381;
t653 = -t447 + t450;
t652 = t449 + t550;
t616 = t489 * t643;
t651 = rSges(3,2) * t616 + rSges(3,3) * t642;
t644 = qJD(1) * t446;
t635 = qJD(5) * t401;
t632 = -pkin(3) - t791;
t88 = t223 * pkin(3) + t675;
t631 = -t88 + t678;
t627 = -t183 + t677;
t519 = -t491 * t643 - t610;
t204 = -qJ(3) * t488 * t616 + pkin(2) * t519 + t571;
t626 = t492 * t204 + t490 * t205 + t361 * t642;
t414 = qJD(1) * (-pkin(1) * t643 + t478);
t625 = qJD(1) * t204 + qJD(3) * t756 + t414;
t623 = -rSges(5,1) * t611 + rSges(5,2) * t385 - rSges(5,3) * t384 + t676;
t107 = -t221 * rSges(4,1) + t220 * rSges(4,2) + rSges(4,3) * t581;
t622 = -rSges(5,1) * t690 + rSges(5,2) * t403 - rSges(5,3) * t402 + t665;
t110 = rSges(5,1) * t581 + t221 * rSges(5,2) - t220 * rSges(5,3);
t620 = t392 * t640 + t393 * t639 + t467;
t615 = t300 * t640;
t614 = t301 * t640;
t606 = -pkin(1) - t720;
t604 = t492 * t633;
t601 = -t640 / 0.2e1;
t598 = t639 / 0.2e1;
t79 = -t614 - t752;
t590 = t79 * t669;
t582 = -t376 + t693;
t580 = t677 - t679;
t579 = qJD(2) * t467 + t204 * t639 + t205 * t640 + t361 * t604;
t577 = t490 * t206 + t492 * t208 + t661;
t576 = -rSges(6,2) * t384 - t791 * t385 - t611 * t792 - t635 + t676;
t274 = -pkin(3) * t373 - qJ(4) * t372;
t575 = qJD(1) * t274 - qJD(4) * t370 + t656;
t574 = -rSges(6,2) * t402 - t791 * t403 - t690 * t792 + t665;
t573 = t666 - t670;
t572 = -qJD(4) * t220 + t315 * t605 + t663;
t45 = -t615 - t733;
t570 = t45 * t621;
t559 = qJD(1) * t87 + qJD(4) * t222 + t625;
t557 = -rSges(3,2) * t489 + t720;
t552 = t362 + t647;
t549 = -t272 * t490 - t697;
t545 = t40 * t573;
t544 = t573 * t492;
t543 = t206 * t642 + t490 * t88 + t492 * t87 + t626;
t542 = t361 * t640 + t362 * t639 - t607;
t536 = qJD(4) * t402 + t273 * t640 + t274 * t639 + t620;
t412 = t451 * t490;
t269 = (t382 * t490 + t383 * t492) * qJD(2);
t108 = t223 * rSges(4,1) - t222 * rSges(4,2) + rSges(4,3) * t756;
t112 = rSges(5,1) * t756 - t223 * rSges(5,2) + t222 * rSges(5,3);
t518 = qJD(4) * t384 + t206 * t604 + t87 * t639 + t88 * t640 + t579;
t516 = t378 * t492 - t379 * t490;
t515 = t208 + t552;
t514 = t206 * t640 + t208 * t639 + t387 + t542;
t510 = (-t489 * t652 + t491 * t653) * qJD(1);
t494 = t516 * t491 + (-t490 * t654 + t492 * t655) * t489;
t424 = t557 * qJD(2);
t314 = -qJD(2) * t412 + (t492 * t557 + t483) * qJD(1);
t313 = rSges(3,1) * t519 - rSges(3,2) * t609 + t651;
t304 = (t404 * t490 + t405 * t492) * qJD(2);
t247 = -rSges(4,1) * t373 + rSges(4,2) * t372 + rSges(4,3) * t629;
t246 = -rSges(4,1) * t371 + rSges(4,2) * t370 + rSges(4,3) * t630;
t245 = rSges(5,1) * t629 + rSges(5,2) * t373 - rSges(5,3) * t372;
t243 = rSges(5,1) * t630 + rSges(5,2) * t371 - rSges(5,3) * t370;
t147 = -t424 * t639 + (-t314 + t804) * qJD(1);
t146 = -t424 * t640 + t414 + (t313 - t612) * qJD(1);
t68 = (t175 * t490 + t177 * t492) * qJD(2) + t542;
t43 = (t179 * t490 + t183 * t492) * qJD(2) + t514;
t42 = t672 * t639 + (-t108 + t614 + t678) * qJD(1) + t663;
t41 = qJD(1) * t107 + (qJD(1) * t585 + t490 * t672) * qJD(2) + t625;
t39 = t321 + qJD(2) * t544 + (t624 - t681) * qJD(1) + t664;
t32 = t635 + (t490 * t681 + t492 * t679) * qJD(2) + t514;
t31 = (t107 * t492 + t108 * t490 + (t175 * t492 + t490 * t682) * qJD(1)) * qJD(2) + t579;
t30 = t623 * t639 + (-t112 + t615 + t631) * qJD(1) + t572;
t29 = qJD(1) * t110 + (qJD(1) * t565 + t490 * t623) * qJD(2) + t559;
t4 = -qJD(5) * t221 + t576 * t639 + (t631 - t683 + t785) * qJD(1) + t572;
t3 = qJD(5) * t223 + t684 * qJD(1) + (qJD(1) * t544 + t490 * t576) * qJD(2) + t559;
t2 = (t110 * t492 + t112 * t490 + (t179 * t492 + t490 * t627) * qJD(1)) * qJD(2) + t518;
t1 = qJD(5) * t385 + (t684 * t492 + t683 * t490 + (t490 * t580 + t492 * t681) * qJD(1)) * qJD(2) + t518;
t5 = [(-qJD(2) * t546 + t420 * t491 + t421 * t489 + (-t491 * t810 + t738 * t641) * t488 + t808 * t401 + t809 * t400 + t735 * t385 + t737 * t384) * qJD(1) + (t3 * (t515 + t679) + (t632 * t336 - t324 + t763 + t787) * t4 - t545 * t639 + (qJD(1) * t681 + t527 * t643 - t321 + t684 + t789) * t40 + (t632 * t223 + t40 - t675 + t754 + t762) * t39) * m(6) + (-(t490 * t719 + t570 * t492) * qJD(2) + t29 * (t515 + t183) + (t110 + (t179 + t790) * qJD(1) + t789) * t45 + (-t179 - t206 + t763) * t30 + (-t733 - t112 - t88 + t762) * t44) * m(5) + (-(t490 * t712 + t590 * t492) * qJD(2) + t41 * (t552 + t177) + (t107 + (t175 + t790) * qJD(1) + t805) * t79 + (-t175 + t763) * t42 + (-t752 - t108 + t762) * t78) * m(4) + (-(-qJD(1) * t382 - t271 - t432 - t612) * t272 + t147 * (t490 * t606 + t486 + t648) + t146 * (t383 + t647) + t272 * (t478 + t651) + (t451 * t698 - t696) * qJD(2) + ((-pkin(1) - t557) * t697 + (t271 * (-rSges(3,3) - pkin(7)) + t272 * t606) * t490) * qJD(1)) * m(3) + (((t185 - t358 + (t377 + t694) * t492 + t660) * t492 + (t659 + t782) * t490 + t818) * qJD(2) + t796) * t598 + (t769 + t770) * t640 / 0.2e1 - (t768 - t771 + t772) * t639 / 0.2e1 + (((t492 * t582 + t187 - t659) * t492 + (t490 * t582 + t583 + t767 - t783 + t822) * t490 + t817) * qJD(2) + t773 + t793) * t601 + ((t766 - t816) * t490 + (t765 + t815) * t492) * t633 / 0.2e1; (t39 * t667 + t1 * t577 + t32 * t543 + (t4 * t573 + t39 * t576 + t1 * t679 + t32 * t684 + (t32 * t681 + t545) * qJD(1)) * t492 + (t3 * t573 + t40 * t576 + t1 * t681 + t32 * t683 + (t32 * t580 + t39 * t670) * qJD(1)) * t490 - t39 * (-qJD(5) * t373 + t657) - t40 * (-qJD(5) * t371 + t575) - t32 * (qJD(5) * t403 + t536) - (t39 * (t671 - t674) + t40 * t673) * qJD(1) - ((t32 * t673 + t39 * t574) * t492 + (t32 * t674 + t40 * t574) * t490) * qJD(2)) * m(6) + (t44 * t667 + t2 * t577 + t43 * t543 + (t30 * t621 + t44 * t623 + t2 * t183 + t43 * t110 + (t43 * t179 + t570) * qJD(1)) * t492 + (t29 * t621 + t45 * t623 + t2 * t179 + t43 * t112 + (t43 * t627 + t719) * qJD(1)) * t490 - t44 * t657 - t45 * t575 - t43 * t536 - (t44 * (-t243 + t671) + t45 * t245) * qJD(1) - ((t43 * t245 + t44 * t622) * t492 + (t43 * t243 + t45 * t622) * t490) * qJD(2)) * m(5) + (-t78 * t399 + t31 * t661 + t68 * t626 + (t42 * t669 + t78 * t672 + t31 * t177 + t68 * t107 + (t68 * t175 + t590) * qJD(1)) * t492 + (t41 * t669 + t79 * t672 + t31 * t175 + t68 * t108 + (t68 * t682 + t712) * qJD(1)) * t490 - t78 * (t443 + (-t246 - t392) * qJD(1)) - t79 * (qJD(1) * t247 + t656) - t68 * t620 - ((t68 * t247 + t668 * t78) * t492 + (t68 * t246 + t668 * t79) * t490) * qJD(2)) * m(4) + (-(t271 * t412 - t696) * qJD(1) - (t269 * (-t412 * t490 - t413 * t492) + t549 * t557) * qJD(2) + 0.2e1 * t269 * (t313 * t492 + t314 * t490 + (t382 * t492 - t383 * t490) * qJD(1)) + t549 * t424 + (-t146 * t490 - t147 * t492 + (-t272 * t492 + t698) * qJD(1)) * t451) * m(3) - ((t516 * t489 + (t743 * t400 + t745 * t401 + t402 * t811 - t751 * t403 - t655 * t491) * t492 + (t742 * t400 + t744 * t401 + t748 * t402 + t750 * t403 + t654 * t491) * t490 + ((-t489 * t774 + t491 * t741) * t492 + (t489 * t746 - t491 * t740) * t490) * t488) * qJD(2) + (t489 * t653 + t491 * t652 + (t489 * t738 - t491 * t739) * t488 + t735 * t403 + t737 * t402 + t734 * t401 + t736 * t400) * qJD(1)) * qJD(1) / 0.2e1 + (t771 * t492 + t770 * t490 + (t766 * t490 + t765 * t492) * qJD(1)) * qJD(1) / 0.2e1 + ((-t640 * t691 + t644) * t490 + t510 * t492 + ((t337 * t742 + t338 * t744 - t372 * t748 - t373 * t750 + t405 * t740 + t629 * t746) * t490 + (t337 * t743 + t338 * t745 - t372 * t811 + t373 * t751 - t405 * t741 + t490 * t692 - t629 * t774 + t494) * t492) * qJD(2) + (t337 * t736 + t338 * t734 - t372 * t737 - t373 * t735 + t405 * t739 + t629 * t738) * qJD(1)) * t601 + ((-t639 * t692 - t644) * t492 + t510 * t490 + ((t335 * t743 + t336 * t745 - t370 * t811 + t371 * t751 - t404 * t741 - t630 * t774) * t492 + (t335 * t742 + t336 * t744 - t370 * t748 - t371 * t750 + t404 * t740 + t492 * t691 + t630 * t746 + t494) * t490) * qJD(2) + (t335 * t736 + t336 * t734 - t370 * t737 - t371 * t735 + t404 * t739 + t630 * t738) * qJD(1)) * t598 + (t769 * qJD(1) + ((-t220 * t811 + t751 * t221 + t779 * t337 + t777 * t338 - t405 * t812 - t774 * t581) * t492 + (-t748 * t220 - t750 * t221 + t780 * t337 + t778 * t338 + t776 * t405 + t732 * t490 + t746 * t581 - t806) * t490 + (t767 * t490 + t820 * t492) * qJD(1)) * t761) * t490 / 0.2e1 - (t768 * qJD(1) + ((t222 * t811 - t751 * t223 + t779 * t335 + t777 * t336 - t404 * t812 - t774 * t756 + t806) * t492 + (t748 * t222 + t750 * t223 + t780 * t335 + t778 * t336 + t776 * t404 - t492 * t732 + t746 * t756) * t490 + (t821 * t492 + (t184 + t784) * t490) * qJD(1)) * t761) * t492 / 0.2e1 + (t773 + t794) * t643 / 0.2e1 + (t772 + t795) * t642 / 0.2e1; (t3 * t404 + t4 * t405 + (-t1 * t491 + t32 * t641) * t488 - t304 * t32) * m(6) + (t29 * t404 + t30 * t405 + (-t2 * t491 + t43 * t641) * t488 - t304 * t43) * m(5) + (t404 * t41 + t405 * t42 + (-t31 * t491 + t641 * t68) * t488 - t304 * t68) * m(4); (t1 * t400 + t3 * t335 + t32 * t760 + t337 * t4 + t39 * t759 + t40 * t758) * m(6) + (t2 * t400 + t29 * t335 + t30 * t337 + t43 * t760 + t44 * t759 + t45 * t758) * m(5); (t1 * t401 - t221 * t39 + t223 * t40 + t3 * t336 + t32 * t385 + t338 * t4 - (-t336 * t39 + t338 * t40) * qJD(1) - (t32 * (t336 * t490 + t338 * t492) + (-t39 * t492 - t40 * t490) * t401) * qJD(2)) * m(6);];
tauc = t5(:);
