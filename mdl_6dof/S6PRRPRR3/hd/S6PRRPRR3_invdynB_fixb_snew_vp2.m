% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR3
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 04:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:46:28
% EndTime: 2019-05-05 04:47:08
% DurationCPUTime: 39.92s
% Computational Cost: add. (639886->351), mult. (1554289->467), div. (0->0), fcn. (1238975->16), ass. (0->156)
t750 = -2 * qJD(4);
t705 = sin(pkin(13));
t709 = cos(pkin(13));
t715 = sin(qJ(3));
t719 = cos(qJ(3));
t707 = sin(pkin(7));
t738 = qJD(2) * t707;
t683 = (t705 * t715 - t709 * t719) * t738;
t706 = sin(pkin(12));
t710 = cos(pkin(12));
t695 = g(1) * t706 - g(2) * t710;
t696 = -g(1) * t710 - g(2) * t706;
t704 = -g(3) + qJDD(1);
t716 = sin(qJ(2));
t712 = cos(pkin(6));
t720 = cos(qJ(2));
t740 = t712 * t720;
t708 = sin(pkin(6));
t744 = t708 * t720;
t665 = t695 * t740 - t696 * t716 + t704 * t744;
t721 = qJD(2) ^ 2;
t749 = pkin(9) * t707;
t657 = qJDD(2) * pkin(2) + t721 * t749 + t665;
t741 = t712 * t716;
t745 = t708 * t716;
t666 = t695 * t741 + t720 * t696 + t704 * t745;
t658 = -pkin(2) * t721 + qJDD(2) * t749 + t666;
t685 = -t695 * t708 + t704 * t712;
t711 = cos(pkin(7));
t742 = t711 * t719;
t746 = t707 * t719;
t629 = t657 * t742 - t715 * t658 + t685 * t746;
t736 = qJD(2) * qJD(3);
t690 = (qJDD(2) * t715 + t719 * t736) * t707;
t700 = qJDD(2) * t711 + qJDD(3);
t701 = qJD(2) * t711 + qJD(3);
t733 = t719 * t738;
t748 = t707 ^ 2 * t721;
t621 = (t701 * t733 - t690) * qJ(4) + (t715 * t719 * t748 + t700) * pkin(3) + t629;
t743 = t711 * t715;
t747 = t707 * t715;
t630 = t657 * t743 + t719 * t658 + t685 * t747;
t734 = t715 * t738;
t686 = pkin(3) * t701 - qJ(4) * t734;
t691 = (qJDD(2) * t719 - t715 * t736) * t707;
t735 = t719 ^ 2 * t748;
t622 = -pkin(3) * t735 + qJ(4) * t691 - t686 * t701 + t630;
t684 = (t705 * t719 + t709 * t715) * t738;
t611 = t709 * t621 - t705 * t622 + t684 * t750;
t612 = t705 * t621 + t709 * t622 + t683 * t750;
t659 = mrSges(5,1) * t683 + mrSges(5,2) * t684;
t663 = -t690 * t705 + t691 * t709;
t671 = mrSges(5,1) * t701 - mrSges(5,3) * t684;
t660 = pkin(4) * t683 - pkin(10) * t684;
t699 = t701 ^ 2;
t610 = -pkin(4) * t699 + pkin(10) * t700 - t660 * t683 + t612;
t642 = -t707 * t657 + t711 * t685;
t631 = -t691 * pkin(3) - qJ(4) * t735 + t686 * t734 + qJDD(4) + t642;
t664 = t690 * t709 + t691 * t705;
t614 = (t683 * t701 - t664) * pkin(10) + (t684 * t701 - t663) * pkin(4) + t631;
t714 = sin(qJ(5));
t718 = cos(qJ(5));
t607 = t718 * t610 + t714 * t614;
t669 = t684 * t718 + t701 * t714;
t640 = -qJD(5) * t669 - t664 * t714 + t700 * t718;
t668 = -t684 * t714 + t701 * t718;
t643 = -mrSges(6,1) * t668 + mrSges(6,2) * t669;
t682 = qJD(5) + t683;
t649 = mrSges(6,1) * t682 - mrSges(6,3) * t669;
t662 = qJDD(5) - t663;
t644 = -pkin(5) * t668 - pkin(11) * t669;
t681 = t682 ^ 2;
t604 = -pkin(5) * t681 + pkin(11) * t662 + t644 * t668 + t607;
t609 = -t700 * pkin(4) - t699 * pkin(10) + t684 * t660 - t611;
t641 = qJD(5) * t668 + t664 * t718 + t700 * t714;
t605 = (-t668 * t682 - t641) * pkin(11) + (t669 * t682 - t640) * pkin(5) + t609;
t713 = sin(qJ(6));
t717 = cos(qJ(6));
t601 = -t604 * t713 + t605 * t717;
t646 = -t669 * t713 + t682 * t717;
t617 = qJD(6) * t646 + t641 * t717 + t662 * t713;
t647 = t669 * t717 + t682 * t713;
t627 = -mrSges(7,1) * t646 + mrSges(7,2) * t647;
t667 = qJD(6) - t668;
t632 = -mrSges(7,2) * t667 + mrSges(7,3) * t646;
t639 = qJDD(6) - t640;
t599 = m(7) * t601 + mrSges(7,1) * t639 - mrSges(7,3) * t617 - t627 * t647 + t632 * t667;
t602 = t604 * t717 + t605 * t713;
t616 = -qJD(6) * t647 - t641 * t713 + t662 * t717;
t633 = mrSges(7,1) * t667 - mrSges(7,3) * t647;
t600 = m(7) * t602 - mrSges(7,2) * t639 + mrSges(7,3) * t616 + t627 * t646 - t633 * t667;
t729 = -t599 * t713 + t717 * t600;
t592 = m(6) * t607 - mrSges(6,2) * t662 + mrSges(6,3) * t640 + t643 * t668 - t649 * t682 + t729;
t606 = -t610 * t714 + t614 * t718;
t648 = -mrSges(6,2) * t682 + mrSges(6,3) * t668;
t603 = -pkin(5) * t662 - pkin(11) * t681 + t644 * t669 - t606;
t724 = -m(7) * t603 + t616 * mrSges(7,1) - mrSges(7,2) * t617 + t646 * t632 - t633 * t647;
t597 = m(6) * t606 + mrSges(6,1) * t662 - mrSges(6,3) * t641 - t643 * t669 + t648 * t682 + t724;
t730 = t718 * t592 - t597 * t714;
t584 = m(5) * t612 - mrSges(5,2) * t700 + mrSges(5,3) * t663 - t659 * t683 - t671 * t701 + t730;
t670 = -mrSges(5,2) * t701 - mrSges(5,3) * t683;
t593 = t599 * t717 + t600 * t713;
t722 = -m(6) * t609 + t640 * mrSges(6,1) - mrSges(6,2) * t641 + t668 * t648 - t649 * t669 - t593;
t589 = m(5) * t611 + mrSges(5,1) * t700 - mrSges(5,3) * t664 - t659 * t684 + t670 * t701 + t722;
t579 = t705 * t584 + t709 * t589;
t688 = -mrSges(4,2) * t701 + mrSges(4,3) * t733;
t689 = (-mrSges(4,1) * t719 + mrSges(4,2) * t715) * t738;
t577 = m(4) * t629 + mrSges(4,1) * t700 - mrSges(4,3) * t690 + t688 * t701 - t689 * t734 + t579;
t687 = mrSges(4,1) * t701 - mrSges(4,3) * t734;
t731 = t709 * t584 - t589 * t705;
t578 = m(4) * t630 - mrSges(4,2) * t700 + mrSges(4,3) * t691 - t687 * t701 + t689 * t733 + t731;
t587 = t714 * t592 + t718 * t597;
t723 = m(5) * t631 - t663 * mrSges(5,1) + t664 * mrSges(5,2) + t683 * t670 + t684 * t671 + t587;
t586 = m(4) * t642 - t691 * mrSges(4,1) + t690 * mrSges(4,2) + (t687 * t715 - t688 * t719) * t738 + t723;
t564 = t577 * t742 + t578 * t743 - t586 * t707;
t560 = m(3) * t665 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t721 + t564;
t563 = t577 * t746 + t578 * t747 + t711 * t586;
t562 = m(3) * t685 + t563;
t570 = -t577 * t715 + t719 * t578;
t569 = m(3) * t666 - mrSges(3,1) * t721 - qJDD(2) * mrSges(3,2) + t570;
t550 = t560 * t740 - t562 * t708 + t569 * t741;
t548 = m(2) * t695 + t550;
t556 = -t560 * t716 + t720 * t569;
t555 = m(2) * t696 + t556;
t739 = t710 * t548 + t706 * t555;
t549 = t560 * t744 + t712 * t562 + t569 * t745;
t732 = -t548 * t706 + t710 * t555;
t623 = Ifges(7,5) * t647 + Ifges(7,6) * t646 + Ifges(7,3) * t667;
t625 = Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t667;
t594 = -mrSges(7,1) * t603 + mrSges(7,3) * t602 + Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t639 - t623 * t647 + t625 * t667;
t624 = Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t667;
t595 = mrSges(7,2) * t603 - mrSges(7,3) * t601 + Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t639 + t623 * t646 - t624 * t667;
t634 = Ifges(6,5) * t669 + Ifges(6,6) * t668 + Ifges(6,3) * t682;
t635 = Ifges(6,4) * t669 + Ifges(6,2) * t668 + Ifges(6,6) * t682;
t580 = mrSges(6,2) * t609 - mrSges(6,3) * t606 + Ifges(6,1) * t641 + Ifges(6,4) * t640 + Ifges(6,5) * t662 - pkin(11) * t593 - t594 * t713 + t595 * t717 + t634 * t668 - t635 * t682;
t636 = Ifges(6,1) * t669 + Ifges(6,4) * t668 + Ifges(6,5) * t682;
t581 = -mrSges(6,1) * t609 - mrSges(7,1) * t601 + mrSges(7,2) * t602 + mrSges(6,3) * t607 + Ifges(6,4) * t641 - Ifges(7,5) * t617 + Ifges(6,2) * t640 + Ifges(6,6) * t662 - Ifges(7,6) * t616 - Ifges(7,3) * t639 - pkin(5) * t593 - t624 * t647 + t625 * t646 - t634 * t669 + t636 * t682;
t651 = Ifges(5,4) * t684 - Ifges(5,2) * t683 + Ifges(5,6) * t701;
t652 = Ifges(5,1) * t684 - Ifges(5,4) * t683 + Ifges(5,5) * t701;
t675 = Ifges(4,6) * t701 + (Ifges(4,4) * t715 + Ifges(4,2) * t719) * t738;
t676 = Ifges(4,5) * t701 + (Ifges(4,1) * t715 + Ifges(4,4) * t719) * t738;
t557 = Ifges(4,5) * t690 + Ifges(4,6) * t691 + mrSges(4,1) * t629 - mrSges(4,2) * t630 + Ifges(5,5) * t664 + Ifges(5,6) * t663 + t684 * t651 + t683 * t652 + mrSges(5,1) * t611 - mrSges(5,2) * t612 + t714 * t580 + t718 * t581 + pkin(4) * t722 + pkin(10) * t730 + pkin(3) * t579 + (Ifges(4,3) + Ifges(5,3)) * t700 + (t675 * t715 - t676 * t719) * t738;
t650 = Ifges(5,5) * t684 - Ifges(5,6) * t683 + Ifges(5,3) * t701;
t565 = mrSges(5,2) * t631 - mrSges(5,3) * t611 + Ifges(5,1) * t664 + Ifges(5,4) * t663 + Ifges(5,5) * t700 - pkin(10) * t587 + t580 * t718 - t581 * t714 - t650 * t683 - t651 * t701;
t571 = Ifges(5,4) * t664 + Ifges(5,2) * t663 + Ifges(5,6) * t700 - t684 * t650 + t701 * t652 - mrSges(5,1) * t631 + mrSges(5,3) * t612 - Ifges(6,5) * t641 - Ifges(6,6) * t640 - Ifges(6,3) * t662 - t669 * t635 + t668 * t636 - mrSges(6,1) * t606 + mrSges(6,2) * t607 - t713 * t595 - t717 * t594 - pkin(5) * t724 - pkin(11) * t729 - pkin(4) * t587;
t674 = Ifges(4,3) * t701 + (Ifges(4,5) * t715 + Ifges(4,6) * t719) * t738;
t551 = -mrSges(4,1) * t642 + mrSges(4,3) * t630 + Ifges(4,4) * t690 + Ifges(4,2) * t691 + Ifges(4,6) * t700 - pkin(3) * t723 + qJ(4) * t731 + t705 * t565 + t709 * t571 - t674 * t734 + t701 * t676;
t552 = mrSges(4,2) * t642 - mrSges(4,3) * t629 + Ifges(4,1) * t690 + Ifges(4,4) * t691 + Ifges(4,5) * t700 - qJ(4) * t579 + t565 * t709 - t571 * t705 + t674 * t733 - t675 * t701;
t725 = pkin(9) * t570 + t551 * t719 + t552 * t715;
t545 = -mrSges(3,1) * t685 + mrSges(3,3) * t666 + t721 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t563 - t707 * t557 + t725 * t711;
t546 = mrSges(3,2) * t685 - mrSges(3,3) * t665 + Ifges(3,5) * qJDD(2) - t721 * Ifges(3,6) - t715 * t551 + t719 * t552 + (-t563 * t707 - t564 * t711) * pkin(9);
t726 = pkin(8) * t556 + t545 * t720 + t546 * t716;
t544 = mrSges(3,1) * t665 - mrSges(3,2) * t666 + Ifges(3,3) * qJDD(2) + pkin(2) * t564 + t711 * t557 + t725 * t707;
t543 = mrSges(2,2) * t704 - mrSges(2,3) * t695 - t716 * t545 + t720 * t546 + (-t549 * t708 - t550 * t712) * pkin(8);
t542 = -mrSges(2,1) * t704 + mrSges(2,3) * t696 - pkin(1) * t549 - t708 * t544 + t726 * t712;
t1 = [-m(1) * g(1) + t732; -m(1) * g(2) + t739; -m(1) * g(3) + m(2) * t704 + t549; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t739 - t706 * t542 + t710 * t543; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t732 + t710 * t542 + t706 * t543; -mrSges(1,1) * g(2) + mrSges(2,1) * t695 + mrSges(1,2) * g(1) - mrSges(2,2) * t696 + pkin(1) * t550 + t712 * t544 + t726 * t708;];
tauB  = t1;
