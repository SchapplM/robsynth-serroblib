% Calculate kinetic energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:13
% EndTime: 2018-11-23 10:29:17
% DurationCPUTime: 4.19s
% Computational Cost: add. (34197->416), mult. (34877->628), div. (0->0), fcn. (34988->30), ass. (0->202)
t702 = pkin(6) - qJ(2);
t693 = cos(t702) / 0.2e1;
t701 = pkin(6) + qJ(2);
t698 = cos(t701);
t684 = t693 + t698 / 0.2e1;
t713 = sin(qJ(2));
t714 = sin(qJ(1));
t719 = cos(qJ(1));
t670 = t684 * t719 - t713 * t714;
t691 = sin(t701) / 0.2e1;
t696 = sin(t702);
t680 = t691 - t696 / 0.2e1;
t718 = cos(qJ(2));
t671 = t680 * t719 + t714 * t718;
t705 = sin(pkin(6));
t742 = t705 * t719;
t645 = Icges(3,5) * t671 + Icges(3,6) * t670 - Icges(3,3) * t742;
t672 = -t684 * t714 - t713 * t719;
t673 = -t680 * t714 + t718 * t719;
t743 = t705 * t714;
t646 = Icges(3,5) * t673 + Icges(3,6) * t672 + Icges(3,3) * t743;
t746 = (t645 * t719 - t646 * t714) * t705;
t744 = cos(qJ(5));
t704 = sin(pkin(7));
t707 = cos(pkin(7));
t660 = -t670 * t704 - t707 * t742;
t641 = pkin(2) * t671 + pkin(11) * t660;
t661 = -t672 * t704 + t707 * t743;
t642 = pkin(2) * t673 + pkin(11) * t661;
t739 = qJD(2) * t705;
t688 = t714 * t739;
t736 = t719 * t739;
t741 = t641 * t688 + t642 * t736;
t655 = qJD(3) * t661 + t688;
t740 = qJD(1) * (pkin(1) * t714 - pkin(10) * t742);
t708 = cos(pkin(6));
t689 = qJD(2) * t708 + qJD(1);
t738 = pkin(8) - qJ(4);
t737 = pkin(8) + qJ(4);
t699 = pkin(7) + qJ(3);
t690 = sin(t699) / 0.2e1;
t700 = pkin(7) - qJ(3);
t695 = sin(t700);
t677 = t690 + t695 / 0.2e1;
t692 = cos(t700) / 0.2e1;
t697 = cos(t699);
t682 = t692 + t697 / 0.2e1;
t712 = sin(qJ(3));
t635 = t672 * t682 - t673 * t712 + t677 * t743;
t703 = sin(pkin(8));
t706 = cos(pkin(8));
t624 = -t635 * t703 + t661 * t706;
t608 = qJD(4) * t624 + t655;
t679 = t691 + t696 / 0.2e1;
t669 = -t679 * t704 + t707 * t708;
t667 = qJD(3) * t669 + t689;
t735 = cos(t737);
t734 = sin(t738);
t678 = t690 - t695 / 0.2e1;
t683 = t692 - t697 / 0.2e1;
t717 = cos(qJ(3));
t636 = t672 * t678 + t673 * t717 + t683 * t743;
t711 = sin(qJ(4));
t732 = sin(t737) / 0.2e1;
t727 = t732 + t734 / 0.2e1;
t733 = cos(t738) / 0.2e1;
t728 = t733 + t735 / 0.2e1;
t595 = -t635 * t728 + t636 * t711 - t661 * t727;
t579 = qJD(5) * t595 + t608;
t685 = t693 - t698 / 0.2e1;
t643 = t677 * t708 + t679 * t682 - t685 * t712;
t632 = -t643 * t703 + t669 * t706;
t625 = qJD(4) * t632 + t667;
t656 = qJD(3) * t660 - t736;
t644 = t678 * t679 + t683 * t708 + t685 * t717;
t611 = -t643 * t728 + t644 * t711 - t669 * t727;
t586 = qJD(5) * t611 + t625;
t633 = t670 * t682 - t671 * t712 - t677 * t742;
t623 = -t633 * t703 + t660 * t706;
t634 = t670 * t678 + t671 * t717 - t683 * t742;
t597 = pkin(3) * t634 + pkin(12) * t623;
t598 = pkin(3) * t636 + pkin(12) * t624;
t730 = t655 * t597 - t598 * t656 + t741;
t609 = qJD(4) * t623 + t656;
t657 = pkin(2) * t685 + pkin(11) * t669;
t674 = qJD(1) * (pkin(1) * t719 + pkin(10) * t743);
t729 = t689 * t642 - t657 * t688 + t674;
t593 = -t633 * t728 + t634 * t711 - t660 * t727;
t580 = qJD(5) * t593 + t609;
t676 = t732 - t734 / 0.2e1;
t681 = t733 - t735 / 0.2e1;
t716 = cos(qJ(4));
t594 = t633 * t676 + t634 * t716 + t660 * t681;
t571 = pkin(4) * t594 + pkin(13) * t593;
t596 = t635 * t676 + t636 * t716 + t661 * t681;
t572 = pkin(4) * t596 + pkin(13) * t595;
t726 = t608 * t571 - t572 * t609 + t730;
t725 = -t641 * t689 - t657 * t736 - t740;
t613 = pkin(3) * t644 + pkin(12) * t632;
t724 = t667 * t598 - t613 * t655 + t729;
t723 = -t597 * t667 + t656 * t613 + t725;
t612 = t643 * t676 + t644 * t716 + t669 * t681;
t585 = pkin(4) * t612 + pkin(13) * t611;
t722 = t625 * t572 - t585 * t608 + t724;
t721 = -t571 * t625 + t609 * t585 + t723;
t715 = cos(qJ(6));
t710 = sin(qJ(5));
t709 = sin(qJ(6));
t687 = rSges(2,1) * t719 - rSges(2,2) * t714;
t686 = rSges(2,1) * t714 + rSges(2,2) * t719;
t666 = rSges(3,1) * t685 + rSges(3,2) * t679 + rSges(3,3) * t708;
t665 = Icges(3,1) * t685 + Icges(3,4) * t679 + Icges(3,5) * t708;
t664 = Icges(3,4) * t685 + Icges(3,2) * t679 + Icges(3,6) * t708;
t663 = Icges(3,5) * t685 + Icges(3,6) * t679 + Icges(3,3) * t708;
t652 = rSges(3,1) * t673 + rSges(3,2) * t672 + rSges(3,3) * t743;
t651 = rSges(3,1) * t671 + rSges(3,2) * t670 - rSges(3,3) * t742;
t650 = Icges(3,1) * t673 + Icges(3,4) * t672 + Icges(3,5) * t743;
t649 = Icges(3,1) * t671 + Icges(3,4) * t670 - Icges(3,5) * t742;
t648 = Icges(3,4) * t673 + Icges(3,2) * t672 + Icges(3,6) * t743;
t647 = Icges(3,4) * t671 + Icges(3,2) * t670 - Icges(3,6) * t742;
t627 = t652 * t689 - t666 * t688 + t674;
t626 = -t651 * t689 - t666 * t736 - t740;
t618 = (t651 * t714 + t652 * t719) * t739;
t617 = rSges(4,1) * t644 + rSges(4,2) * t643 + rSges(4,3) * t669;
t616 = Icges(4,1) * t644 + Icges(4,4) * t643 + Icges(4,5) * t669;
t615 = Icges(4,4) * t644 + Icges(4,2) * t643 + Icges(4,6) * t669;
t614 = Icges(4,5) * t644 + Icges(4,6) * t643 + Icges(4,3) * t669;
t606 = rSges(4,1) * t636 + rSges(4,2) * t635 + rSges(4,3) * t661;
t605 = rSges(4,1) * t634 + rSges(4,2) * t633 + rSges(4,3) * t660;
t604 = Icges(4,1) * t636 + Icges(4,4) * t635 + Icges(4,5) * t661;
t603 = Icges(4,1) * t634 + Icges(4,4) * t633 + Icges(4,5) * t660;
t602 = Icges(4,4) * t636 + Icges(4,2) * t635 + Icges(4,6) * t661;
t601 = Icges(4,4) * t634 + Icges(4,2) * t633 + Icges(4,6) * t660;
t600 = Icges(4,5) * t636 + Icges(4,6) * t635 + Icges(4,3) * t661;
t599 = Icges(4,5) * t634 + Icges(4,6) * t633 + Icges(4,3) * t660;
t588 = t612 * t744 + t632 * t710;
t587 = t612 * t710 - t632 * t744;
t584 = t596 * t744 + t624 * t710;
t583 = t596 * t710 - t624 * t744;
t582 = t594 * t744 + t623 * t710;
t581 = t594 * t710 - t623 * t744;
t578 = rSges(5,1) * t612 - rSges(5,2) * t611 + rSges(5,3) * t632;
t577 = Icges(5,1) * t612 - Icges(5,4) * t611 + Icges(5,5) * t632;
t576 = Icges(5,4) * t612 - Icges(5,2) * t611 + Icges(5,6) * t632;
t575 = Icges(5,5) * t612 - Icges(5,6) * t611 + Icges(5,3) * t632;
t574 = t588 * t715 + t611 * t709;
t573 = -t588 * t709 + t611 * t715;
t569 = pkin(5) * t588 + pkin(14) * t587;
t568 = qJD(6) * t587 + t586;
t567 = t606 * t667 - t617 * t655 + t729;
t566 = -t605 * t667 + t617 * t656 + t725;
t564 = rSges(5,1) * t596 - rSges(5,2) * t595 + rSges(5,3) * t624;
t563 = rSges(5,1) * t594 - rSges(5,2) * t593 + rSges(5,3) * t623;
t562 = Icges(5,1) * t596 - Icges(5,4) * t595 + Icges(5,5) * t624;
t561 = Icges(5,1) * t594 - Icges(5,4) * t593 + Icges(5,5) * t623;
t560 = Icges(5,4) * t596 - Icges(5,2) * t595 + Icges(5,6) * t624;
t559 = Icges(5,4) * t594 - Icges(5,2) * t593 + Icges(5,6) * t623;
t558 = Icges(5,5) * t596 - Icges(5,6) * t595 + Icges(5,3) * t624;
t557 = Icges(5,5) * t594 - Icges(5,6) * t593 + Icges(5,3) * t623;
t556 = t584 * t715 + t595 * t709;
t555 = -t584 * t709 + t595 * t715;
t554 = t582 * t715 + t593 * t709;
t553 = -t582 * t709 + t593 * t715;
t551 = t605 * t655 - t606 * t656 + t741;
t550 = rSges(6,1) * t588 - rSges(6,2) * t587 + rSges(6,3) * t611;
t549 = Icges(6,1) * t588 - Icges(6,4) * t587 + Icges(6,5) * t611;
t548 = Icges(6,4) * t588 - Icges(6,2) * t587 + Icges(6,6) * t611;
t547 = Icges(6,5) * t588 - Icges(6,6) * t587 + Icges(6,3) * t611;
t546 = pkin(5) * t584 + pkin(14) * t583;
t545 = pkin(5) * t582 + pkin(14) * t581;
t544 = qJD(6) * t581 + t580;
t543 = qJD(6) * t583 + t579;
t542 = rSges(6,1) * t584 - rSges(6,2) * t583 + rSges(6,3) * t595;
t541 = rSges(6,1) * t582 - rSges(6,2) * t581 + rSges(6,3) * t593;
t540 = Icges(6,1) * t584 - Icges(6,4) * t583 + Icges(6,5) * t595;
t539 = Icges(6,1) * t582 - Icges(6,4) * t581 + Icges(6,5) * t593;
t538 = Icges(6,4) * t584 - Icges(6,2) * t583 + Icges(6,6) * t595;
t537 = Icges(6,4) * t582 - Icges(6,2) * t581 + Icges(6,6) * t593;
t536 = Icges(6,5) * t584 - Icges(6,6) * t583 + Icges(6,3) * t595;
t535 = Icges(6,5) * t582 - Icges(6,6) * t581 + Icges(6,3) * t593;
t534 = rSges(7,1) * t574 + rSges(7,2) * t573 + rSges(7,3) * t587;
t533 = Icges(7,1) * t574 + Icges(7,4) * t573 + Icges(7,5) * t587;
t532 = Icges(7,4) * t574 + Icges(7,2) * t573 + Icges(7,6) * t587;
t531 = Icges(7,5) * t574 + Icges(7,6) * t573 + Icges(7,3) * t587;
t530 = rSges(7,1) * t556 + rSges(7,2) * t555 + rSges(7,3) * t583;
t529 = rSges(7,1) * t554 + rSges(7,2) * t553 + rSges(7,3) * t581;
t528 = Icges(7,1) * t556 + Icges(7,4) * t555 + Icges(7,5) * t583;
t527 = Icges(7,1) * t554 + Icges(7,4) * t553 + Icges(7,5) * t581;
t526 = Icges(7,4) * t556 + Icges(7,2) * t555 + Icges(7,6) * t583;
t525 = Icges(7,4) * t554 + Icges(7,2) * t553 + Icges(7,6) * t581;
t524 = Icges(7,5) * t556 + Icges(7,6) * t555 + Icges(7,3) * t583;
t523 = Icges(7,5) * t554 + Icges(7,6) * t553 + Icges(7,3) * t581;
t522 = t564 * t625 - t578 * t608 + t724;
t521 = -t563 * t625 + t578 * t609 + t723;
t520 = t563 * t608 - t564 * t609 + t730;
t519 = t542 * t586 - t550 * t579 + t722;
t518 = -t541 * t586 + t550 * t580 + t721;
t517 = t541 * t579 - t542 * t580 + t726;
t516 = t530 * t568 - t534 * t543 + t546 * t586 - t569 * t579 + t722;
t515 = -t529 * t568 + t534 * t544 - t545 * t586 + t569 * t580 + t721;
t514 = t529 * t543 - t530 * t544 + t545 * t579 - t546 * t580 + t726;
t1 = m(5) * (t520 ^ 2 + t521 ^ 2 + t522 ^ 2) / 0.2e1 + m(4) * (t551 ^ 2 + t566 ^ 2 + t567 ^ 2) / 0.2e1 + m(3) * (t618 ^ 2 + t626 ^ 2 + t627 ^ 2) / 0.2e1 + t689 * ((t663 * t708 + t664 * t679 + t665 * t685) * t689 + ((t646 * t708 + t648 * t679 + t650 * t685) * t714 - (t645 * t708 + t647 * t679 + t649 * t685) * t719) * t739) / 0.2e1 + t544 * ((t524 * t581 + t526 * t553 + t528 * t554) * t543 + (t581 * t523 + t553 * t525 + t554 * t527) * t544 + (t531 * t581 + t532 * t553 + t533 * t554) * t568) / 0.2e1 + t568 * ((t524 * t587 + t526 * t573 + t528 * t574) * t543 + (t523 * t587 + t525 * t573 + t527 * t574) * t544 + (t587 * t531 + t573 * t532 + t574 * t533) * t568) / 0.2e1 + t579 * ((t595 * t536 - t583 * t538 + t584 * t540) * t579 + (t535 * t595 - t537 * t583 + t539 * t584) * t580 + (t547 * t595 - t548 * t583 + t549 * t584) * t586) / 0.2e1 + t586 * ((t536 * t611 - t538 * t587 + t540 * t588) * t579 + (t535 * t611 - t537 * t587 + t539 * t588) * t580 + (t611 * t547 - t587 * t548 + t588 * t549) * t586) / 0.2e1 + t543 * ((t583 * t524 + t555 * t526 + t556 * t528) * t543 + (t523 * t583 + t525 * t555 + t527 * t556) * t544 + (t531 * t583 + t532 * t555 + t533 * t556) * t568) / 0.2e1 + t580 * ((t536 * t593 - t538 * t581 + t540 * t582) * t579 + (t593 * t535 - t581 * t537 + t582 * t539) * t580 + (t547 * t593 - t548 * t581 + t549 * t582) * t586) / 0.2e1 + t608 * ((t624 * t558 - t595 * t560 + t596 * t562) * t608 + (t557 * t624 - t559 * t595 + t561 * t596) * t609 + (t575 * t624 - t576 * t595 + t577 * t596) * t625) / 0.2e1 + t609 * ((t558 * t623 - t560 * t593 + t562 * t594) * t608 + (t623 * t557 - t593 * t559 + t594 * t561) * t609 + (t575 * t623 - t576 * t593 + t577 * t594) * t625) / 0.2e1 + t625 * ((t558 * t632 - t560 * t611 + t562 * t612) * t608 + (t557 * t632 - t559 * t611 + t561 * t612) * t609 + (t575 * t632 - t576 * t611 + t577 * t612) * t625) / 0.2e1 + t656 * ((t600 * t660 + t602 * t633 + t604 * t634) * t655 + (t599 * t660 + t601 * t633 + t603 * t634) * t656 + (t614 * t660 + t615 * t633 + t616 * t634) * t667) / 0.2e1 + t655 * ((t600 * t661 + t602 * t635 + t604 * t636) * t655 + (t599 * t661 + t601 * t635 + t603 * t636) * t656 + (t614 * t661 + t615 * t635 + t616 * t636) * t667) / 0.2e1 + t667 * ((t600 * t669 + t602 * t643 + t604 * t644) * t655 + (t599 * t669 + t601 * t643 + t603 * t644) * t656 + (t614 * t669 + t615 * t643 + t616 * t644) * t667) / 0.2e1 + m(7) * (t514 ^ 2 + t515 ^ 2 + t516 ^ 2) / 0.2e1 + m(6) * (t517 ^ 2 + t518 ^ 2 + t519 ^ 2) / 0.2e1 - ((-t663 * t742 + t664 * t670 + t665 * t671) * t689 + ((t648 * t670 + t650 * t671) * t714 + (-t670 * t647 - t671 * t649 + t746) * t719) * t739) * t736 / 0.2e1 + ((t663 * t743 + t664 * t672 + t665 * t673) * t689 + (-(t647 * t672 + t649 * t673) * t719 + (t672 * t648 + t673 * t650 - t746) * t714) * t739) * t688 / 0.2e1 + (m(2) * (t686 ^ 2 + t687 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
