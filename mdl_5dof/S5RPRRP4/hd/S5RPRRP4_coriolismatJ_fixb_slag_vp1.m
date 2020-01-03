% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:17
% EndTime: 2020-01-03 11:50:01
% DurationCPUTime: 25.56s
% Computational Cost: add. (56417->621), mult. (76013->838), div. (0->0), fcn. (82179->8), ass. (0->374)
t855 = Icges(5,4) + Icges(6,4);
t534 = cos(pkin(8));
t532 = qJ(3) + qJ(4);
t523 = sin(t532);
t538 = cos(qJ(1));
t658 = t538 * t523;
t594 = t534 * t658;
t524 = cos(t532);
t536 = sin(qJ(1));
t661 = t536 * t524;
t476 = t594 - t661;
t664 = t534 * t538;
t477 = t523 * t536 + t524 * t664;
t395 = t476 * rSges(6,1) + t477 * rSges(6,2);
t454 = t476 * pkin(4);
t305 = -t395 - t454;
t665 = t534 * t536;
t474 = -t523 * t665 - t524 * t538;
t475 = t534 * t661 - t658;
t394 = t474 * rSges(5,1) - t475 * rSges(5,2);
t396 = rSges(5,1) * t476 + rSges(5,2) * t477;
t393 = t474 * rSges(6,1) - t475 * rSges(6,2);
t453 = t474 * pkin(4);
t621 = -t393 - t453;
t752 = m(6) / 0.2e1;
t753 = m(5) / 0.2e1;
t649 = (-t305 * t538 + t536 * t621) * t752 + (-t394 * t536 + t396 * t538) * t753;
t533 = sin(pkin(8));
t539 = -pkin(7) - pkin(6);
t530 = -qJ(5) + t539;
t603 = t530 - t539;
t537 = cos(qJ(3));
t528 = t537 * pkin(3);
t503 = pkin(4) * t524 + t528;
t502 = pkin(2) + t503;
t522 = t528 + pkin(2);
t606 = t502 - t522;
t619 = (t603 - rSges(6,3)) * t534 + (rSges(6,1) * t524 - rSges(6,2) * t523 + t606) * t533;
t535 = sin(qJ(3));
t657 = t538 * t535;
t518 = pkin(3) * t657;
t671 = t533 * t539;
t605 = t536 * t671 + t518;
t673 = t533 * t536;
t704 = pkin(4) * t523;
t571 = pkin(3) * t535 + t704;
t824 = rSges(6,1) * t475 + rSges(6,2) * t474 - t530 * t673 - t538 * t571;
t641 = rSges(6,3) * t673 + t606 * t665 + t605 + t824;
t163 = t534 * t641 + t619 * t673;
t569 = rSges(5,1) * t475 + rSges(5,2) * t474;
t352 = rSges(5,3) * t673 + t569;
t443 = -rSges(5,3) * t534 + (rSges(5,1) * t524 - rSges(5,2) * t523) * t533;
t267 = t352 * t534 + t443 * t673;
t672 = t533 * t538;
t354 = t477 * rSges(5,1) - t476 * rSges(5,2) + rSges(5,3) * t672;
t331 = t354 * t534;
t835 = -t443 * t672 - t331;
t663 = t535 * t536;
t516 = pkin(3) * t663;
t607 = t522 * t664 + t516;
t834 = t477 * rSges(6,1) - t476 * rSges(6,2) + rSges(6,3) * t672 + t502 * t664 + t536 * t571;
t825 = t603 * t672 + t607 - t834;
t645 = t825 * t534;
t841 = -t619 * t672 + t645;
t656 = (t163 * t538 + t841 * t536) * t752 + (t267 * t538 + t835 * t536) * t753;
t18 = t656 - t649;
t854 = qJD(1) * t18;
t848 = -Icges(5,1) - Icges(6,1);
t847 = Icges(5,5) + Icges(6,5);
t846 = Icges(5,2) + Icges(6,2);
t845 = Icges(5,6) + Icges(6,6);
t853 = t855 * t524;
t852 = t855 * t523;
t701 = pkin(6) + t539;
t705 = pkin(2) * t534;
t398 = (t533 * t701 + t705) * t538 - t607;
t371 = t534 * t398;
t702 = -pkin(2) + t522;
t431 = t533 * t702 + t534 * t701;
t588 = t431 + t619;
t562 = t588 * t672;
t133 = t371 - t562 + t645;
t134 = (t398 + t825) * t534 - t562;
t655 = t133 - t134;
t851 = t655 * t752;
t843 = t845 * t534 + (t846 * t523 - t853) * t533;
t850 = t847 * t534 + (t848 * t524 + t852) * t533;
t786 = (t848 * t523 - t853) * t533;
t618 = t431 + t443;
t573 = t618 * t672;
t203 = -t331 + t371 - t573;
t204 = (-t354 + t398) * t534 - t573;
t652 = t203 - t204;
t849 = t652 * t753;
t830 = -t850 + (-t846 * t524 - t852) * t533;
t829 = -t786 - t843;
t693 = Icges(4,4) * t537;
t465 = -Icges(4,6) * t534 + (-Icges(4,2) * t535 + t693) * t533;
t694 = Icges(4,4) * t535;
t466 = -Icges(4,5) * t534 + (Icges(4,1) * t537 - t694) * t533;
t492 = (-Icges(4,2) * t537 - t694) * t533;
t493 = (-Icges(4,1) * t535 - t693) * t533;
t491 = (-Icges(4,5) * t535 - Icges(4,6) * t537) * t533;
t666 = t534 * t491;
t844 = t666 / 0.2e1 - (-(t466 / 0.2e1 + t492 / 0.2e1) * t535 + (t493 / 0.2e1 - t465 / 0.2e1) * t537) * t533;
t580 = t672 / 0.2e1;
t581 = -t672 / 0.2e1;
t582 = t673 / 0.2e1;
t334 = Icges(6,5) * t477 - Icges(6,6) * t476 + Icges(6,3) * t672;
t337 = Icges(5,5) * t477 - Icges(5,6) * t476 + Icges(5,3) * t672;
t678 = t337 * t536;
t336 = Icges(5,5) * t475 + Icges(5,6) * t474 + Icges(5,3) * t673;
t679 = t336 * t538;
t680 = t334 * t536;
t333 = Icges(6,5) * t475 + Icges(6,6) * t474 + Icges(6,3) * t673;
t681 = t333 * t538;
t761 = t538 ^ 2;
t798 = t533 * t761;
t452 = Icges(5,4) * t477;
t344 = Icges(5,2) * t476 - Icges(5,6) * t672 - t452;
t451 = Icges(5,4) * t476;
t349 = Icges(5,1) * t477 + Icges(5,5) * t672 - t451;
t818 = -t474 * t344 + t349 * t475;
t449 = Icges(6,4) * t477;
t341 = Icges(6,2) * t476 - Icges(6,6) * t672 - t449;
t448 = Icges(6,4) * t476;
t346 = Icges(6,1) * t477 + Icges(6,5) * t672 - t448;
t820 = -t474 * t341 + t346 * t475;
t836 = t337 * t672;
t837 = t336 * t672;
t838 = t334 * t672;
t839 = t333 * t672;
t172 = -t334 * t673 - t820;
t174 = -t337 * t673 - t818;
t840 = -t172 - t174;
t777 = ((t334 + t337) * t798 + (-t836 - t838) * t538 + (t820 + t818 + (t680 + t681 + t678 + t679) * t533 - t837 - t839 - t840) * t536) * t533;
t468 = (-Icges(6,5) * t523 - Icges(6,6) * t524) * t533;
t469 = (-Icges(5,5) * t523 - Icges(5,6) * t524) * t533;
t831 = t468 + t469;
t794 = -t830 * t476 - t829 * t477 + t831 * t672;
t795 = -t830 * t474 + t829 * t475 - t831 * t673;
t782 = t848 * t476 + t341 + t344 - t449 - t452;
t828 = t846 * t477 - t346 - t349 + t448 + t451;
t832 = t847 * t476 + t845 * t477;
t796 = t832 * t534 + (t828 * t523 + t782 * t524) * t533;
t689 = Icges(6,4) * t475;
t339 = Icges(6,2) * t474 + Icges(6,6) * t673 + t689;
t692 = Icges(5,4) * t475;
t342 = Icges(5,2) * t474 + Icges(5,6) * t673 + t692;
t781 = t848 * t474 + t339 + t342 + t689 + t692;
t447 = Icges(6,4) * t474;
t345 = Icges(6,1) * t475 + Icges(6,5) * t673 + t447;
t450 = Icges(5,4) * t474;
t348 = Icges(5,1) * t475 + Icges(5,5) * t673 + t450;
t783 = -t846 * t475 + t345 + t348 + t447 + t450;
t833 = t847 * t474 - t845 * t475;
t797 = -t833 * t534 + (-t783 * t523 - t781 * t524) * t533;
t171 = t333 * t673 + t474 * t339 + t475 * t345;
t173 = t336 * t673 + t474 * t342 + t475 * t348;
t823 = (t850 * t475 + t843 * t474 + ((Icges(5,3) + Icges(6,3)) * t534 + (t845 * t523 - t847 * t524) * t533) * t673) * t534;
t811 = t823 + (t840 * t538 + (t171 + t173) * t536) * t533;
t812 = ((t171 + t838) * t536 + ((-t680 + t681) * t533 - t172 - t839) * t538) * t533 + ((t173 + t836) * t536 + ((-t678 + t679) * t533 - t174 - t837) * t538) * t533 + t823;
t541 = t812 * t580 - t777 * t673 / 0.2e1 + (-t795 + t797) * t582 + (-t794 - t796 + t811) * t581;
t667 = t534 * t469;
t668 = t534 * t468;
t650 = t667 + t668 + (t523 * t830 + t829 * t524) * t533;
t779 = t650 * t534;
t842 = t541 + t779;
t531 = t533 ^ 2;
t593 = t534 * t657;
t660 = t536 * t537;
t496 = t593 - t660;
t659 = t537 * t538;
t497 = t534 * t659 + t663;
t385 = Icges(4,5) * t497 - Icges(4,6) * t496 + Icges(4,3) * t672;
t826 = t385 * t672;
t400 = t497 * rSges(4,1) - t496 * rSges(4,2) + rSges(4,3) * t672;
t467 = -rSges(4,3) * t534 + (rSges(4,1) * t537 - rSges(4,2) * t535) * t533;
t821 = -t400 * t534 - t467 * t672;
t485 = Icges(4,4) * t497;
t389 = Icges(4,2) * t496 - Icges(4,6) * t672 - t485;
t484 = Icges(4,4) * t496;
t391 = Icges(4,1) * t497 + Icges(4,5) * t672 - t484;
t495 = t534 * t660 - t657;
t552 = t534 * t663 + t659;
t816 = t389 * t552 + t391 * t495;
t814 = -t832 * t672 + t833 * t673;
t754 = m(4) / 0.2e1;
t384 = Icges(4,5) * t495 - Icges(4,6) * t552 + Icges(4,3) * t673;
t791 = t384 * t672;
t131 = t133 * t536;
t703 = pkin(6) * t533;
t397 = (t534 * t702 - t703) * t536 - t605;
t132 = (t397 + t641) * t534 + t588 * t673;
t201 = t203 * t536;
t202 = (t352 + t397) * t534 + t618 * t673;
t570 = rSges(4,1) * t495 - rSges(4,2) * t552;
t399 = rSges(4,3) * t673 + t570;
t280 = t399 * t534 + t467 * t673;
t598 = (t202 * t538 + t201) * t753 + (t280 * t538 + t821 * t536) * t754 + (t132 * t538 + t131) * t752;
t601 = (-t134 * t536 + t131) * t752 + (-t204 * t536 + t201) * t753;
t785 = t598 - t601;
t478 = (-rSges(6,1) * t523 - rSges(6,2) * t524) * t533;
t675 = t523 * t533;
t576 = t533 * (pkin(4) * t675 - t478);
t572 = t536 * t576;
t486 = t552 * pkin(3);
t559 = t534 * t571;
t575 = -t503 * t538 - t536 * t559;
t363 = t575 + t486;
t628 = -t363 - t393;
t216 = t534 * t628 + t572;
t517 = pkin(3) * t660;
t662 = t536 * t503;
t364 = pkin(4) * t594 + t517 - t662;
t369 = t534 * t395;
t217 = t534 * t364 + t538 * t576 + t369;
t526 = t538 * qJ(2);
t253 = -t526 + (rSges(6,3) * t533 + t502 * t534 + pkin(1)) * t536 + t824;
t277 = t575 + t393;
t278 = -t538 * t559 - t395 + t662;
t487 = pkin(3) * t593 - t517;
t309 = -t396 - t487;
t620 = -t394 + t486;
t265 = -t526 + (rSges(5,3) * t533 + t522 * t534 + pkin(1)) * t536 + t569 - t605;
t479 = (-rSges(5,1) * t523 - rSges(5,2) * t524) * t533;
t595 = t479 * t673;
t283 = -t394 * t534 - t595;
t284 = t534 * t396 - t479 * t672;
t604 = t538 * pkin(1) + t536 * qJ(2);
t768 = -t538 * t671 + t354 + t604 + t607;
t651 = t283 * t265 + t284 * t768;
t769 = -t530 * t672 + t604 + t834;
t699 = (-t163 * t277 + t216 * t253 + t217 * t769 + t278 * t841) * t752 + (t267 * t620 + t309 * t835 + t651) * t753;
t561 = -t478 * t533 + t531 * t704;
t245 = t534 * t621 + t536 * t561;
t246 = t454 * t534 + t538 * t561 + t369;
t654 = t245 * t253 + t246 * t769;
t700 = (t132 * t621 + t133 * t305 + t654) * t752 + (-t202 * t394 - t203 * t396 + t651) * t753;
t778 = t699 - t700;
t776 = -t667 / 0.2e1 - t668 / 0.2e1 - t830 * t675 / 0.2e1;
t775 = t533 / 0.2e1;
t730 = -t534 / 0.2e1;
t774 = m(6) * t533;
t767 = (t703 + t705) * t538 + t604 + t400;
t764 = t163 * t851 + t267 * t849;
t411 = -rSges(4,1) * t552 - rSges(4,2) * t495;
t412 = rSges(4,1) * t496 + rSges(4,2) * t497;
t589 = (-t277 * t536 - t278 * t538) * t752 + (-t309 * t538 + t536 * t620) * t753 + (-t411 * t536 + t412 * t538) * t754;
t759 = 0.4e1 * qJD(1);
t758 = 2 * qJD(3);
t757 = 4 * qJD(3);
t756 = 2 * qJD(4);
t749 = m(5) * t202 * t652;
t243 = t352 * t672 - t354 * t673;
t627 = t397 * t672 + t398 * t673;
t154 = t627 + t243;
t258 = t394 * t672 + t396 * t673;
t63 = t154 * t258 - t202 * t283 + t203 * t284;
t748 = m(5) * t63;
t94 = t243 * t258 - t267 * t283 + t284 * t835;
t93 = m(5) * t94;
t742 = m(6) * t132 * t655;
t75 = -t132 * t673 + t133 * t672;
t741 = m(6) * ((t132 * t536 - t134 * t538) * t533 + t75);
t143 = t641 * t672 + t673 * t825;
t107 = t143 + t627;
t629 = t393 * t672 + t395 * t673;
t218 = (t453 * t538 + t454 * t536) * t533 + t629;
t40 = t107 * t218 - t132 * t245 + t133 * t246;
t739 = m(6) * t40;
t161 = t363 * t672 + t364 * t673 + t629;
t738 = m(6) * (t143 * t161 - t163 * t216 + t217 * t841);
t100 = -t163 * t673 + t672 * t841;
t735 = m(6) * ((t163 * t536 - t538 * t841) * t533 + t100);
t732 = m(6) * t75;
t729 = m(3) * (-(-rSges(3,2) * t673 - t538 * rSges(3,3) + pkin(1) * t536 - t526) * t538 + (-rSges(3,2) * t672 + rSges(3,3) * t536 + t604) * t536);
t301 = -t526 + (t705 + pkin(1) + (rSges(4,3) + pkin(6)) * t533) * t536 + t570;
t727 = m(4) * (t301 * t411 - t412 * t767);
t725 = m(4) * (-t301 * t538 + t767 * t536);
t721 = m(5) * (-t265 * t620 + t309 * t768);
t720 = m(5) * (t265 * t394 - t396 * t768);
t719 = m(5) * (-t265 * t538 + t768 * t536);
t715 = m(6) * t100;
t713 = m(6) * (t253 * t277 + t278 * t769);
t712 = m(6) * (-t253 * t621 + t305 * t769);
t711 = m(6) * (-t253 * t538 + t769 * t536);
t710 = (-t277 * t538 + t278 * t536) * t774;
t708 = (t305 * t536 + t538 * t621) * t774;
t698 = m(6) * qJD(1);
t695 = Icges(4,4) * t495;
t677 = t384 * t538;
t676 = t385 * t536;
t674 = t524 * t533;
t152 = t253 * t673 + t672 * t769;
t387 = -Icges(4,2) * t552 + Icges(4,6) * t673 + t695;
t625 = Icges(4,1) * t552 + t387 + t695;
t624 = -Icges(4,1) * t496 + t389 - t485;
t483 = Icges(4,4) * t552;
t390 = Icges(4,1) * t495 + Icges(4,5) * t673 - t483;
t623 = -Icges(4,2) * t495 + t390 - t483;
t622 = Icges(4,2) * t497 - t391 + t484;
t613 = -t486 * t672 + t487 * t673;
t612 = t534 * t487 + t531 * t518;
t611 = t465 - t493;
t610 = t466 + t492;
t403 = (t536 ^ 2 + t761) * t774;
t602 = t403 * qJD(1);
t196 = t384 * t673 - t387 * t552 + t495 * t390;
t197 = -t385 * t673 - t816;
t600 = -((t196 + t826) * t536 + ((-t676 + t677) * t533 - t197 - t791) * t538) * t533 / 0.2e1 + (t196 * t536 - t197 * t538) * t775;
t599 = (t385 * t798 - t826 * t538 + (t197 + (t676 + t677) * t533 - t791 + t816) * t536) * t775 + (t534 / 0.2e1 + t730) * ((-Icges(4,3) * t534 + (Icges(4,5) * t537 - Icges(4,6) * t535) * t533) * t672 - t465 * t496 + t466 * t497);
t590 = t143 * t218 - t163 * t245 + t246 * t841;
t585 = -t674 / 0.2e1;
t584 = t674 / 0.2e1;
t567 = ((t797 * t536 + t796 * t538) * t533 + t779) * t730 + ((-t474 * t828 + t782 * t475) * t672 + t795 * t534 + (t783 * t474 - t781 * t475 + t814) * t673) * t582 + ((t783 * t476 + t781 * t477) * t673 + t794 * t534 + (-t476 * t828 - t782 * t477 - t814) * t672) * t581;
t542 = (-t534 * t161 + (-t216 * t538 + t217 * t536) * t533) * t752;
t543 = m(6) * (-t534 * t218 + (-t245 * t538 + t246 * t536) * t533);
t54 = t542 - t543 / 0.2e1;
t550 = (-t216 * t536 - t217 * t538) * t752;
t551 = m(6) * (-t245 * t536 - t246 * t538);
t79 = t550 - t551 / 0.2e1;
t566 = t79 * qJD(2) + t54 * qJD(5);
t560 = t93 + t567;
t547 = t811 * t580 + t812 * t581 + t777 * t582;
t545 = t786 * t584 - t585 * t843 + t776;
t544 = t547 + t764;
t540 = -t584 * t843 + t786 * t585 - t776;
t505 = t531 * t516;
t498 = (-rSges(4,1) * t535 - rSges(4,2) * t537) * t533;
t406 = Icges(4,5) * t496 + Icges(4,6) * t497;
t405 = -Icges(4,5) * t552 - Icges(4,6) * t495;
t307 = t412 * t534 - t498 * t672;
t306 = -t411 * t534 - t498 * t673;
t251 = t284 + t612;
t250 = t534 * t620 + t505 - t595;
t242 = -t666 + (-t535 * t610 - t537 * t611) * t533;
t223 = t708 / 0.2e1;
t220 = t613 + t258;
t214 = m(5) * (-t283 * t536 - t284 * t538);
t208 = -t491 * t672 + t496 * t610 + t497 * t611;
t207 = t491 * t673 - t495 * t611 - t552 * t610;
t200 = t710 / 0.2e1;
t184 = t217 + t612;
t183 = t505 + t572 + (t486 + t628) * t534;
t151 = -t406 * t534 + (-t535 * t622 - t537 * t624) * t533;
t150 = -t405 * t534 + (-t535 * t623 - t537 * t625) * t533;
t146 = t161 + t613;
t99 = t715 / 0.2e1;
t74 = t732 / 0.2e1;
t72 = t711 + t719 + t725 + t729;
t66 = t214 + t551 / 0.2e1 + t550;
t56 = t545 + t712 + t720;
t55 = t543 / 0.2e1 + t542;
t48 = t735 / 0.2e1;
t41 = t545 + t727 + t721 + t713 - t844;
t25 = t741 / 0.2e1;
t22 = t99 + t223 - t735 / 0.2e1;
t21 = t99 + t48 - t708 / 0.2e1;
t20 = t223 + t48 - t715 / 0.2e1;
t17 = t649 + t656;
t15 = t74 + t200 - t741 / 0.2e1;
t14 = t74 + t25 - t710 / 0.2e1;
t13 = t200 + t25 - t732 / 0.2e1;
t10 = t598 + t601 - t589;
t9 = t589 + t785;
t8 = t589 - t785;
t7 = t560 + t738;
t6 = t567 + t739 + t748;
t4 = t547 + (t536 * t599 + t538 * t600) * t533 + t749 + t742;
t3 = t544 - t778;
t2 = t544 + t778;
t1 = t699 + t700 - t764 + t842;
t5 = [m(6) * t152 * qJD(5) + t72 * qJD(2) + t41 * qJD(3) + t56 * qJD(4), qJD(1) * t72 + qJD(3) * t9 + qJD(4) * t17, t41 * qJD(1) + t9 * qJD(2) + t1 * qJD(4) + t15 * qJD(5) + (-t749 / 0.4e1 - t742 / 0.4e1) * t757 + ((-t280 * t411 + t301 * t306 + t307 * t767 - t412 * t821) * t754 + (t202 * t620 + t203 * t309 + t250 * t265 + t251 * t768) * t753 + (-t132 * t277 + t133 * t278 + t183 * t253 + t184 * t769) * t752) * t758 + (t541 + (-t242 + t650) * t534 + ((-t151 / 0.2e1 - t208 / 0.2e1 - t600) * t538 + (t150 / 0.2e1 + t207 / 0.2e1 - t599) * t536) * t533) * qJD(3), t56 * qJD(1) + t17 * qJD(2) + t1 * qJD(3) + t842 * qJD(4) + t22 * qJD(5) + ((-t267 * t394 - t396 * t835 + t651) * t753 + (t163 * t621 + t305 * t841 + t654) * t752) * t756, t15 * qJD(3) + t22 * qJD(4) + t152 * t698; t8 * qJD(3) - t18 * qJD(4) - t403 * qJD(5) + (-t729 / 0.4e1 - t725 / 0.4e1 - t719 / 0.4e1 - t711 / 0.4e1) * t759, 0, t8 * qJD(1) + t66 * qJD(4) + ((-t306 * t536 - t307 * t538) * t754 + (-t250 * t536 - t251 * t538) * t753 + (-t183 * t536 - t184 * t538) * t752) * t758, -t854 + t66 * qJD(3) + (t214 + t551) * qJD(4), -t602; t10 * qJD(2) + t4 * qJD(3) + t3 * qJD(4) + t14 * qJD(5) + (-t727 / 0.4e1 - t721 / 0.4e1 - t713 / 0.4e1) * t759 + 0.2e1 * (-t253 * t851 - t265 * t849) * qJD(1) + (t540 + t844) * qJD(1), qJD(1) * t10 - qJD(4) * t79, t4 * qJD(1) + (m(5) * (t154 * t220 - t202 * t250 + t203 * t251) + m(4) * (-t280 * t306 + t821 * t307 + (t399 * t538 - t400 * t536) * t531 * (t411 * t538 + t412 * t536)) + (-t242 * t534 + (t150 * t536 - t151 * t538) * t533) * t730 + (-t208 * t534 + (-t405 * t672 + t496 * t623 + t497 * t625) * t673 - (-t406 * t672 + t496 * t622 + t497 * t624) * t672) * t581 + (-t207 * t534 + (t405 * t673 - t495 * t625 - t552 * t623) * t673 - (t406 * t673 - t495 * t624 - t552 * t622) * t672) * t582 + m(6) * (t107 * t146 - t132 * t183 + t133 * t184) + t567) * qJD(3) + t6 * qJD(4), t3 * qJD(1) + t6 * qJD(3) + ((t590 + t40) * t752 + (t94 + t63) * t753) * t756 - t566 + (t567 - t93 - t738) * qJD(4), qJD(1) * t14 - qJD(4) * t54; t540 * qJD(1) + t18 * qJD(2) + t2 * qJD(3) + t547 * qJD(4) + t21 * qJD(5) + (-t720 / 0.4e1 - t712 / 0.4e1) * t759, qJD(3) * t79 + t854, t2 * qJD(1) + t567 * qJD(3) + t7 * qJD(4) + (-t739 / 0.4e1 - t748 / 0.4e1) * t757 + ((t107 * t161 - t132 * t216 + t133 * t217 + t143 * t146 - t163 * t183 + t184 * t841) * t752 + (t220 * t243 - t250 * t267 + t251 * t835 + t63) * t753) * t758 + t566, t547 * qJD(1) + t7 * qJD(3) + (m(6) * t590 + t560) * qJD(4), qJD(1) * t21 + qJD(3) * t54; (-t253 * t536 - t538 * t769) * t533 * t698 + t403 * qJD(2) + t13 * qJD(3) + t20 * qJD(4), t602, t13 * qJD(1) + m(6) * (-t146 * t534 + (-t183 * t538 + t184 * t536) * t533) * qJD(3) + t55 * qJD(4), t20 * qJD(1) + t55 * qJD(3) + qJD(4) * t543, 0;];
Cq = t5;
