% Calculate kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:52
% EndTime: 2018-12-10 18:09:57
% DurationCPUTime: 4.53s
% Computational Cost: add. (33653->417), mult. (34205->613), div. (0->0), fcn. (34252->30), ass. (0->199)
t698 = pkin(6) - qJ(2);
t689 = cos(t698) / 0.2e1;
t697 = pkin(6) + qJ(2);
t694 = cos(t697);
t680 = t689 + t694 / 0.2e1;
t710 = sin(qJ(2));
t711 = sin(qJ(1));
t715 = cos(qJ(1));
t666 = t680 * t715 - t710 * t711;
t688 = sin(t697) / 0.2e1;
t693 = sin(t698);
t678 = t688 - t693 / 0.2e1;
t714 = cos(qJ(2));
t667 = t678 * t715 + t711 * t714;
t702 = sin(pkin(6));
t742 = t702 * t715;
t644 = Icges(3,5) * t667 + Icges(3,6) * t666 - Icges(3,3) * t742;
t668 = -t680 * t711 - t710 * t715;
t669 = -t678 * t711 + t714 * t715;
t743 = t702 * t711;
t645 = Icges(3,5) * t669 + Icges(3,6) * t668 + Icges(3,3) * t743;
t747 = t702 * (t644 * t715 - t645 * t711);
t744 = cos(qJ(5));
t695 = pkin(7) + pkin(14);
t686 = sin(t695) / 0.2e1;
t696 = pkin(7) - pkin(14);
t690 = sin(t696);
t671 = t686 + t690 / 0.2e1;
t687 = cos(t696) / 0.2e1;
t691 = cos(t695);
t673 = t687 + t691 / 0.2e1;
t699 = sin(pkin(14));
t634 = t668 * t673 - t669 * t699 + t671 * t743;
t701 = sin(pkin(7));
t705 = cos(pkin(7));
t658 = -t668 * t701 + t705 * t743;
t700 = sin(pkin(8));
t704 = cos(pkin(8));
t623 = -t634 * t700 + t658 * t704;
t740 = qJD(2) * t702;
t684 = t711 * t740;
t615 = qJD(4) * t623 + t684;
t741 = qJD(1) * (pkin(1) * t711 - pkin(10) * t742);
t706 = cos(pkin(6));
t685 = qJD(2) * t706 + qJD(1);
t739 = pkin(8) - qJ(4);
t738 = pkin(8) + qJ(4);
t641 = pkin(2) * t669 + qJ(3) * t658;
t657 = -t666 * t701 - t705 * t742;
t670 = qJD(1) * (pkin(1) * t715 + pkin(10) * t743);
t737 = qJD(3) * t657 + t685 * t641 + t670;
t640 = pkin(2) * t667 + qJ(3) * t657;
t677 = t688 + t693 / 0.2e1;
t665 = -t677 * t701 + t705 * t706;
t735 = t715 * t740;
t736 = qJD(3) * t665 + t640 * t684 + t641 * t735;
t672 = t686 - t690 / 0.2e1;
t674 = t687 - t691 / 0.2e1;
t703 = cos(pkin(14));
t635 = t668 * t672 + t669 * t703 + t674 * t743;
t709 = sin(qJ(4));
t728 = sin(t738) / 0.2e1;
t733 = sin(t739);
t721 = t728 + t733 / 0.2e1;
t729 = cos(t739) / 0.2e1;
t734 = cos(t738);
t722 = t729 + t734 / 0.2e1;
t592 = -t634 * t722 + t635 * t709 - t658 * t721;
t578 = qJD(5) * t592 + t615;
t681 = t689 - t694 / 0.2e1;
t642 = t671 * t706 + t673 * t677 - t681 * t699;
t631 = -t642 * t700 + t665 * t704;
t626 = qJD(4) * t631 + t685;
t732 = qJD(3) * t658 - t741;
t643 = t672 * t677 + t674 * t706 + t681 * t703;
t608 = -t642 * t722 + t643 * t709 - t665 * t721;
t585 = qJD(5) * t608 + t626;
t654 = pkin(2) * t681 + qJ(3) * t665;
t727 = (-pkin(3) * t643 - pkin(11) * t631 - t654) * t740;
t726 = (-rSges(4,1) * t643 - rSges(4,2) * t642 - rSges(4,3) * t665 - t654) * t740;
t632 = t666 * t673 - t667 * t699 - t671 * t742;
t622 = -t632 * t700 + t657 * t704;
t616 = qJD(4) * t622 - t735;
t633 = t666 * t672 + t667 * t703 - t674 * t742;
t597 = pkin(3) * t633 + pkin(11) * t622;
t598 = pkin(3) * t635 + pkin(11) * t623;
t725 = t597 * t684 + t598 * t735 + t736;
t590 = -t632 * t722 + t633 * t709 - t657 * t721;
t579 = qJD(5) * t590 + t616;
t676 = t728 - t733 / 0.2e1;
t679 = t729 - t734 / 0.2e1;
t713 = cos(qJ(4));
t591 = t632 * t676 + t633 * t713 + t657 * t679;
t570 = pkin(4) * t591 + pkin(12) * t590;
t593 = t634 * t676 + t635 * t713 + t658 * t679;
t571 = pkin(4) * t593 + pkin(12) * t592;
t723 = t615 * t570 - t571 * t616 + t725;
t720 = t685 * t598 + t711 * t727 + t737;
t609 = t642 * t676 + t643 * t713 + t665 * t679;
t584 = pkin(4) * t609 + pkin(12) * t608;
t719 = t626 * t571 - t584 * t615 + t720;
t718 = (-t597 - t640) * t685 + t715 * t727 + t732;
t717 = -t570 * t626 + t616 * t584 + t718;
t712 = cos(qJ(6));
t708 = sin(qJ(5));
t707 = sin(qJ(6));
t683 = rSges(2,1) * t715 - rSges(2,2) * t711;
t682 = rSges(2,1) * t711 + rSges(2,2) * t715;
t663 = rSges(3,1) * t681 + rSges(3,2) * t677 + rSges(3,3) * t706;
t662 = Icges(3,1) * t681 + Icges(3,4) * t677 + Icges(3,5) * t706;
t661 = Icges(3,4) * t681 + Icges(3,2) * t677 + Icges(3,6) * t706;
t660 = Icges(3,5) * t681 + Icges(3,6) * t677 + Icges(3,3) * t706;
t651 = rSges(3,1) * t669 + rSges(3,2) * t668 + rSges(3,3) * t743;
t650 = rSges(3,1) * t667 + rSges(3,2) * t666 - rSges(3,3) * t742;
t649 = Icges(3,1) * t669 + Icges(3,4) * t668 + Icges(3,5) * t743;
t648 = Icges(3,1) * t667 + Icges(3,4) * t666 - Icges(3,5) * t742;
t647 = Icges(3,4) * t669 + Icges(3,2) * t668 + Icges(3,6) * t743;
t646 = Icges(3,4) * t667 + Icges(3,2) * t666 - Icges(3,6) * t742;
t625 = t651 * t685 - t663 * t684 + t670;
t624 = -t650 * t685 - t663 * t735 - t741;
t619 = (t650 * t711 + t651 * t715) * t740;
t613 = Icges(4,1) * t643 + Icges(4,4) * t642 + Icges(4,5) * t665;
t612 = Icges(4,4) * t643 + Icges(4,2) * t642 + Icges(4,6) * t665;
t611 = Icges(4,5) * t643 + Icges(4,6) * t642 + Icges(4,3) * t665;
t606 = rSges(4,1) * t635 + rSges(4,2) * t634 + rSges(4,3) * t658;
t605 = rSges(4,1) * t633 + rSges(4,2) * t632 + rSges(4,3) * t657;
t604 = Icges(4,1) * t635 + Icges(4,4) * t634 + Icges(4,5) * t658;
t603 = Icges(4,1) * t633 + Icges(4,4) * t632 + Icges(4,5) * t657;
t602 = Icges(4,4) * t635 + Icges(4,2) * t634 + Icges(4,6) * t658;
t601 = Icges(4,4) * t633 + Icges(4,2) * t632 + Icges(4,6) * t657;
t600 = Icges(4,5) * t635 + Icges(4,6) * t634 + Icges(4,3) * t658;
t599 = Icges(4,5) * t633 + Icges(4,6) * t632 + Icges(4,3) * t657;
t587 = t609 * t744 + t631 * t708;
t586 = t609 * t708 - t631 * t744;
t583 = t593 * t744 + t623 * t708;
t582 = t593 * t708 - t623 * t744;
t581 = t591 * t744 + t622 * t708;
t580 = t591 * t708 - t622 * t744;
t577 = rSges(5,1) * t609 - rSges(5,2) * t608 + rSges(5,3) * t631;
t576 = Icges(5,1) * t609 - Icges(5,4) * t608 + Icges(5,5) * t631;
t575 = Icges(5,4) * t609 - Icges(5,2) * t608 + Icges(5,6) * t631;
t574 = Icges(5,5) * t609 - Icges(5,6) * t608 + Icges(5,3) * t631;
t573 = t587 * t712 + t608 * t707;
t572 = -t587 * t707 + t608 * t712;
t568 = t606 * t685 + t711 * t726 + t737;
t567 = (-t605 - t640) * t685 + t715 * t726 + t732;
t566 = pkin(5) * t587 + pkin(13) * t586;
t565 = qJD(6) * t586 + t585;
t563 = rSges(5,1) * t593 - rSges(5,2) * t592 + rSges(5,3) * t623;
t562 = rSges(5,1) * t591 - rSges(5,2) * t590 + rSges(5,3) * t622;
t561 = Icges(5,1) * t593 - Icges(5,4) * t592 + Icges(5,5) * t623;
t560 = Icges(5,1) * t591 - Icges(5,4) * t590 + Icges(5,5) * t622;
t559 = Icges(5,4) * t593 - Icges(5,2) * t592 + Icges(5,6) * t623;
t558 = Icges(5,4) * t591 - Icges(5,2) * t590 + Icges(5,6) * t622;
t557 = Icges(5,5) * t593 - Icges(5,6) * t592 + Icges(5,3) * t623;
t556 = Icges(5,5) * t591 - Icges(5,6) * t590 + Icges(5,3) * t622;
t555 = (t605 * t711 + t606 * t715) * t740 + t736;
t554 = t583 * t712 + t592 * t707;
t553 = -t583 * t707 + t592 * t712;
t552 = t581 * t712 + t590 * t707;
t551 = -t581 * t707 + t590 * t712;
t549 = rSges(6,1) * t587 - rSges(6,2) * t586 + rSges(6,3) * t608;
t548 = Icges(6,1) * t587 - Icges(6,4) * t586 + Icges(6,5) * t608;
t547 = Icges(6,4) * t587 - Icges(6,2) * t586 + Icges(6,6) * t608;
t546 = Icges(6,5) * t587 - Icges(6,6) * t586 + Icges(6,3) * t608;
t545 = pkin(5) * t583 + pkin(13) * t582;
t544 = pkin(5) * t581 + pkin(13) * t580;
t543 = qJD(6) * t580 + t579;
t542 = qJD(6) * t582 + t578;
t541 = rSges(6,1) * t583 - rSges(6,2) * t582 + rSges(6,3) * t592;
t540 = rSges(6,1) * t581 - rSges(6,2) * t580 + rSges(6,3) * t590;
t539 = Icges(6,1) * t583 - Icges(6,4) * t582 + Icges(6,5) * t592;
t538 = Icges(6,1) * t581 - Icges(6,4) * t580 + Icges(6,5) * t590;
t537 = Icges(6,4) * t583 - Icges(6,2) * t582 + Icges(6,6) * t592;
t536 = Icges(6,4) * t581 - Icges(6,2) * t580 + Icges(6,6) * t590;
t535 = Icges(6,5) * t583 - Icges(6,6) * t582 + Icges(6,3) * t592;
t534 = Icges(6,5) * t581 - Icges(6,6) * t580 + Icges(6,3) * t590;
t533 = rSges(7,1) * t573 + rSges(7,2) * t572 + rSges(7,3) * t586;
t532 = Icges(7,1) * t573 + Icges(7,4) * t572 + Icges(7,5) * t586;
t531 = Icges(7,4) * t573 + Icges(7,2) * t572 + Icges(7,6) * t586;
t530 = Icges(7,5) * t573 + Icges(7,6) * t572 + Icges(7,3) * t586;
t529 = rSges(7,1) * t554 + rSges(7,2) * t553 + rSges(7,3) * t582;
t528 = rSges(7,1) * t552 + rSges(7,2) * t551 + rSges(7,3) * t580;
t527 = Icges(7,1) * t554 + Icges(7,4) * t553 + Icges(7,5) * t582;
t526 = Icges(7,1) * t552 + Icges(7,4) * t551 + Icges(7,5) * t580;
t525 = Icges(7,4) * t554 + Icges(7,2) * t553 + Icges(7,6) * t582;
t524 = Icges(7,4) * t552 + Icges(7,2) * t551 + Icges(7,6) * t580;
t523 = Icges(7,5) * t554 + Icges(7,6) * t553 + Icges(7,3) * t582;
t522 = Icges(7,5) * t552 + Icges(7,6) * t551 + Icges(7,3) * t580;
t521 = t563 * t626 - t577 * t615 + t720;
t520 = -t562 * t626 + t577 * t616 + t718;
t519 = t562 * t615 - t563 * t616 + t725;
t518 = t541 * t585 - t549 * t578 + t719;
t517 = -t540 * t585 + t549 * t579 + t717;
t516 = t540 * t578 - t541 * t579 + t723;
t515 = t529 * t565 - t533 * t542 + t545 * t585 - t566 * t578 + t719;
t514 = -t528 * t565 + t533 * t543 - t544 * t585 + t566 * t579 + t717;
t513 = t528 * t542 - t529 * t543 + t544 * t578 - t545 * t579 + t723;
t1 = t543 * ((t523 * t580 + t525 * t551 + t527 * t552) * t542 + (t580 * t522 + t551 * t524 + t552 * t526) * t543 + (t530 * t580 + t531 * t551 + t532 * t552) * t565) / 0.2e1 + t542 * ((t582 * t523 + t553 * t525 + t554 * t527) * t542 + (t522 * t582 + t524 * t553 + t526 * t554) * t543 + (t530 * t582 + t531 * t553 + t532 * t554) * t565) / 0.2e1 + t565 * ((t523 * t586 + t525 * t572 + t527 * t573) * t542 + (t522 * t586 + t524 * t572 + t526 * t573) * t543 + (t586 * t530 + t572 * t531 + t573 * t532) * t565) / 0.2e1 + t579 * ((t535 * t590 - t537 * t580 + t539 * t581) * t578 + (t590 * t534 - t580 * t536 + t581 * t538) * t579 + (t546 * t590 - t547 * t580 + t548 * t581) * t585) / 0.2e1 + t578 * ((t592 * t535 - t582 * t537 + t583 * t539) * t578 + (t534 * t592 - t536 * t582 + t538 * t583) * t579 + (t546 * t592 - t547 * t582 + t548 * t583) * t585) / 0.2e1 + t585 * ((t535 * t608 - t537 * t586 + t539 * t587) * t578 + (t534 * t608 - t536 * t586 + t538 * t587) * t579 + (t608 * t546 - t586 * t547 + t587 * t548) * t585) / 0.2e1 + t615 * ((t623 * t557 - t592 * t559 + t593 * t561) * t615 + (t556 * t623 - t558 * t592 + t560 * t593) * t616 + (t574 * t623 - t575 * t592 + t576 * t593) * t626) / 0.2e1 + t616 * ((t557 * t622 - t559 * t590 + t561 * t591) * t615 + (t622 * t556 - t590 * t558 + t591 * t560) * t616 + (t574 * t622 - t575 * t590 + t576 * t591) * t626) / 0.2e1 + t626 * ((t557 * t631 - t559 * t608 + t561 * t609) * t615 + (t556 * t631 - t558 * t608 + t560 * t609) * t616 + (t574 * t631 - t575 * t608 + t576 * t609) * t626) / 0.2e1 + m(7) * (t513 ^ 2 + t514 ^ 2 + t515 ^ 2) / 0.2e1 + m(6) * (t516 ^ 2 + t517 ^ 2 + t518 ^ 2) / 0.2e1 + m(5) * (t519 ^ 2 + t520 ^ 2 + t521 ^ 2) / 0.2e1 + m(4) * (t555 ^ 2 + t567 ^ 2 + t568 ^ 2) / 0.2e1 + m(3) * (t619 ^ 2 + t624 ^ 2 + t625 ^ 2) / 0.2e1 + (((-t599 * t665 - t601 * t642 - t603 * t643 - t644 * t706 - t646 * t677 - t648 * t681) * t715 + (t600 * t665 + t602 * t642 + t604 * t643 + t645 * t706 + t647 * t677 + t649 * t681) * t711) * t740 + (t611 * t665 + t612 * t642 + t613 * t643 + t660 * t706 + t661 * t677 + t662 * t681) * t685) * t685 / 0.2e1 + (Icges(2,3) + m(2) * (t682 ^ 2 + t683 ^ 2)) * qJD(1) ^ 2 / 0.2e1 - (((-t599 * t657 - t601 * t632 - t603 * t633 - t666 * t646 - t667 * t648 + t747) * t715 + (t600 * t657 + t602 * t632 + t604 * t633 + t647 * t666 + t649 * t667) * t711) * t740 + (t611 * t657 + t612 * t632 + t613 * t633 - t660 * t742 + t661 * t666 + t662 * t667) * t685) * t735 / 0.2e1 + (((-t599 * t658 - t601 * t634 - t603 * t635 - t646 * t668 - t648 * t669) * t715 + (t600 * t658 + t602 * t634 + t604 * t635 + t668 * t647 + t669 * t649 - t747) * t711) * t740 + (t611 * t658 + t612 * t634 + t613 * t635 + t660 * t743 + t661 * t668 + t662 * t669) * t685) * t684 / 0.2e1;
T  = t1;
