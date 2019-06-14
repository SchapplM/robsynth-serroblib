% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:47:55
% EndTime: 2019-05-07 20:48:45
% DurationCPUTime: 37.76s
% Computational Cost: add. (615502->385), mult. (1286818->484), div. (0->0), fcn. (955229->12), ass. (0->148)
t731 = sin(qJ(1));
t736 = cos(qJ(1));
t719 = -g(1) * t736 - g(2) * t731;
t738 = qJD(1) ^ 2;
t703 = -pkin(1) * t738 + qJDD(1) * pkin(7) + t719;
t730 = sin(qJ(2));
t735 = cos(qJ(2));
t693 = -g(3) * t730 + t735 * t703;
t711 = (-mrSges(3,1) * t735 + mrSges(3,2) * t730) * qJD(1);
t751 = qJD(1) * qJD(2);
t722 = t730 * t751;
t714 = qJDD(1) * t735 - t722;
t753 = qJD(1) * t730;
t716 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t753;
t718 = g(1) * t731 - t736 * g(2);
t702 = -qJDD(1) * pkin(1) - pkin(7) * t738 - t718;
t750 = t735 * t751;
t713 = qJDD(1) * t730 + t750;
t671 = (-t713 - t750) * pkin(8) + (-t714 + t722) * pkin(2) + t702;
t712 = (-pkin(2) * t735 - pkin(8) * t730) * qJD(1);
t737 = qJD(2) ^ 2;
t752 = qJD(1) * t735;
t674 = -pkin(2) * t737 + qJDD(2) * pkin(8) + t712 * t752 + t693;
t729 = sin(qJ(3));
t734 = cos(qJ(3));
t653 = t734 * t671 - t674 * t729;
t709 = qJD(2) * t734 - t729 * t753;
t685 = qJD(3) * t709 + qJDD(2) * t729 + t713 * t734;
t708 = qJDD(3) - t714;
t710 = qJD(2) * t729 + t734 * t753;
t721 = qJD(3) - t752;
t636 = (t709 * t721 - t685) * pkin(9) + (t709 * t710 + t708) * pkin(3) + t653;
t654 = t729 * t671 + t734 * t674;
t684 = -qJD(3) * t710 + qJDD(2) * t734 - t713 * t729;
t694 = pkin(3) * t721 - pkin(9) * t710;
t707 = t709 ^ 2;
t638 = -pkin(3) * t707 + pkin(9) * t684 - t694 * t721 + t654;
t728 = sin(qJ(4));
t733 = cos(qJ(4));
t620 = t733 * t636 - t638 * t728;
t687 = t709 * t733 - t710 * t728;
t652 = qJD(4) * t687 + t684 * t728 + t685 * t733;
t688 = t709 * t728 + t710 * t733;
t704 = qJDD(4) + t708;
t720 = qJD(4) + t721;
t613 = (t687 * t720 - t652) * qJ(5) + (t687 * t688 + t704) * pkin(4) + t620;
t621 = t728 * t636 + t733 * t638;
t651 = -qJD(4) * t688 + t684 * t733 - t685 * t728;
t676 = pkin(4) * t720 - qJ(5) * t688;
t686 = t687 ^ 2;
t619 = -pkin(4) * t686 + qJ(5) * t651 - t676 * t720 + t621;
t725 = sin(pkin(11));
t726 = cos(pkin(11));
t667 = t687 * t725 + t688 * t726;
t607 = -0.2e1 * qJD(5) * t667 + t726 * t613 - t619 * t725;
t633 = t651 * t725 + t652 * t726;
t666 = t687 * t726 - t688 * t725;
t605 = (t666 * t720 - t633) * pkin(10) + (t666 * t667 + t704) * pkin(5) + t607;
t608 = 0.2e1 * qJD(5) * t666 + t725 * t613 + t726 * t619;
t632 = t651 * t726 - t652 * t725;
t657 = pkin(5) * t720 - pkin(10) * t667;
t665 = t666 ^ 2;
t606 = -pkin(5) * t665 + pkin(10) * t632 - t657 * t720 + t608;
t727 = sin(qJ(6));
t732 = cos(qJ(6));
t603 = t605 * t732 - t606 * t727;
t646 = t666 * t732 - t667 * t727;
t617 = qJD(6) * t646 + t632 * t727 + t633 * t732;
t647 = t666 * t727 + t667 * t732;
t629 = -mrSges(7,1) * t646 + mrSges(7,2) * t647;
t715 = qJD(6) + t720;
t639 = -mrSges(7,2) * t715 + mrSges(7,3) * t646;
t698 = qJDD(6) + t704;
t599 = m(7) * t603 + mrSges(7,1) * t698 - mrSges(7,3) * t617 - t629 * t647 + t639 * t715;
t604 = t605 * t727 + t606 * t732;
t616 = -qJD(6) * t647 + t632 * t732 - t633 * t727;
t640 = mrSges(7,1) * t715 - mrSges(7,3) * t647;
t600 = m(7) * t604 - mrSges(7,2) * t698 + mrSges(7,3) * t616 + t629 * t646 - t640 * t715;
t593 = t732 * t599 + t727 * t600;
t648 = -mrSges(6,1) * t666 + mrSges(6,2) * t667;
t655 = -mrSges(6,2) * t720 + mrSges(6,3) * t666;
t591 = m(6) * t607 + mrSges(6,1) * t704 - mrSges(6,3) * t633 - t648 * t667 + t655 * t720 + t593;
t656 = mrSges(6,1) * t720 - mrSges(6,3) * t667;
t744 = -t599 * t727 + t732 * t600;
t592 = m(6) * t608 - mrSges(6,2) * t704 + mrSges(6,3) * t632 + t648 * t666 - t656 * t720 + t744;
t587 = t726 * t591 + t725 * t592;
t668 = -mrSges(5,1) * t687 + mrSges(5,2) * t688;
t675 = -mrSges(5,2) * t720 + mrSges(5,3) * t687;
t585 = m(5) * t620 + mrSges(5,1) * t704 - mrSges(5,3) * t652 - t668 * t688 + t675 * t720 + t587;
t677 = mrSges(5,1) * t720 - mrSges(5,3) * t688;
t745 = -t591 * t725 + t726 * t592;
t586 = m(5) * t621 - mrSges(5,2) * t704 + mrSges(5,3) * t651 + t668 * t687 - t677 * t720 + t745;
t579 = t733 * t585 + t728 * t586;
t689 = -mrSges(4,1) * t709 + mrSges(4,2) * t710;
t690 = -mrSges(4,2) * t721 + mrSges(4,3) * t709;
t577 = m(4) * t653 + mrSges(4,1) * t708 - mrSges(4,3) * t685 - t689 * t710 + t690 * t721 + t579;
t691 = mrSges(4,1) * t721 - mrSges(4,3) * t710;
t746 = -t585 * t728 + t733 * t586;
t578 = m(4) * t654 - mrSges(4,2) * t708 + mrSges(4,3) * t684 + t689 * t709 - t691 * t721 + t746;
t747 = -t577 * t729 + t734 * t578;
t572 = m(3) * t693 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t714 - qJD(2) * t716 + t711 * t752 + t747;
t692 = -t735 * g(3) - t730 * t703;
t717 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t752;
t673 = -qJDD(2) * pkin(2) - pkin(8) * t737 + t712 * t753 - t692;
t649 = -pkin(3) * t684 - pkin(9) * t707 + t710 * t694 + t673;
t623 = -pkin(4) * t651 - qJ(5) * t686 + t688 * t676 + qJDD(5) + t649;
t610 = -pkin(5) * t632 - pkin(10) * t665 + t657 * t667 + t623;
t743 = m(7) * t610 - t616 * mrSges(7,1) + t617 * mrSges(7,2) - t646 * t639 + t647 * t640;
t742 = m(6) * t623 - t632 * mrSges(6,1) + t633 * mrSges(6,2) - t666 * t655 + t667 * t656 + t743;
t740 = m(5) * t649 - t651 * mrSges(5,1) + t652 * mrSges(5,2) - t687 * t675 + t688 * t677 + t742;
t739 = -m(4) * t673 + t684 * mrSges(4,1) - t685 * mrSges(4,2) + t709 * t690 - t710 * t691 - t740;
t602 = m(3) * t692 + qJDD(2) * mrSges(3,1) - t713 * mrSges(3,3) + qJD(2) * t717 - t711 * t753 + t739;
t748 = t735 * t572 - t602 * t730;
t566 = m(2) * t719 - mrSges(2,1) * t738 - qJDD(1) * mrSges(2,2) + t748;
t573 = t577 * t734 + t578 * t729;
t741 = -m(3) * t702 + t714 * mrSges(3,1) - mrSges(3,2) * t713 - t716 * t753 + t717 * t752 - t573;
t569 = m(2) * t718 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t738 + t741;
t754 = t731 * t566 + t736 * t569;
t567 = t730 * t572 + t735 * t602;
t749 = t736 * t566 - t569 * t731;
t701 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t730 + Ifges(3,4) * t735) * qJD(1);
t700 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t730 + Ifges(3,2) * t735) * qJD(1);
t699 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t730 + Ifges(3,6) * t735) * qJD(1);
t680 = Ifges(4,1) * t710 + Ifges(4,4) * t709 + Ifges(4,5) * t721;
t679 = Ifges(4,4) * t710 + Ifges(4,2) * t709 + Ifges(4,6) * t721;
t678 = Ifges(4,5) * t710 + Ifges(4,6) * t709 + Ifges(4,3) * t721;
t660 = Ifges(5,1) * t688 + Ifges(5,4) * t687 + Ifges(5,5) * t720;
t659 = Ifges(5,4) * t688 + Ifges(5,2) * t687 + Ifges(5,6) * t720;
t658 = Ifges(5,5) * t688 + Ifges(5,6) * t687 + Ifges(5,3) * t720;
t643 = Ifges(6,1) * t667 + Ifges(6,4) * t666 + Ifges(6,5) * t720;
t642 = Ifges(6,4) * t667 + Ifges(6,2) * t666 + Ifges(6,6) * t720;
t641 = Ifges(6,5) * t667 + Ifges(6,6) * t666 + Ifges(6,3) * t720;
t626 = Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t715;
t625 = Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t715;
t624 = Ifges(7,5) * t647 + Ifges(7,6) * t646 + Ifges(7,3) * t715;
t595 = mrSges(7,2) * t610 - mrSges(7,3) * t603 + Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t698 + t624 * t646 - t625 * t715;
t594 = -mrSges(7,1) * t610 + mrSges(7,3) * t604 + Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t698 - t624 * t647 + t626 * t715;
t581 = mrSges(6,2) * t623 - mrSges(6,3) * t607 + Ifges(6,1) * t633 + Ifges(6,4) * t632 + Ifges(6,5) * t704 - pkin(10) * t593 - t594 * t727 + t595 * t732 + t641 * t666 - t642 * t720;
t580 = -mrSges(6,1) * t623 + mrSges(6,3) * t608 + Ifges(6,4) * t633 + Ifges(6,2) * t632 + Ifges(6,6) * t704 - pkin(5) * t743 + pkin(10) * t744 + t732 * t594 + t727 * t595 - t667 * t641 + t720 * t643;
t575 = mrSges(5,2) * t649 - mrSges(5,3) * t620 + Ifges(5,1) * t652 + Ifges(5,4) * t651 + Ifges(5,5) * t704 - qJ(5) * t587 - t580 * t725 + t581 * t726 + t658 * t687 - t659 * t720;
t574 = -mrSges(5,1) * t649 + mrSges(5,3) * t621 + Ifges(5,4) * t652 + Ifges(5,2) * t651 + Ifges(5,6) * t704 - pkin(4) * t742 + qJ(5) * t745 + t726 * t580 + t725 * t581 - t688 * t658 + t720 * t660;
t563 = t709 * t680 - t710 * t679 + Ifges(3,4) * t713 + Ifges(3,2) * t714 - Ifges(7,3) * t698 + qJD(2) * t701 - mrSges(3,1) * t702 - Ifges(4,3) * t708 - t688 * t659 + mrSges(3,3) * t693 - Ifges(4,6) * t684 - Ifges(4,5) * t685 + t687 * t660 - t667 * t642 + t666 * t643 - Ifges(5,6) * t651 - Ifges(5,5) * t652 - mrSges(4,1) * t653 + mrSges(4,2) * t654 + t646 * t626 - t647 * t625 - Ifges(6,6) * t632 - Ifges(6,5) * t633 - mrSges(5,1) * t620 + mrSges(5,2) * t621 - Ifges(7,6) * t616 - Ifges(7,5) * t617 - mrSges(6,1) * t607 + mrSges(6,2) * t608 + mrSges(7,2) * t604 - mrSges(7,1) * t603 + (-Ifges(6,3) - Ifges(5,3)) * t704 - pkin(5) * t593 - pkin(4) * t587 - pkin(3) * t579 - t699 * t753 + Ifges(3,6) * qJDD(2) - pkin(2) * t573;
t562 = mrSges(4,2) * t673 - mrSges(4,3) * t653 + Ifges(4,1) * t685 + Ifges(4,4) * t684 + Ifges(4,5) * t708 - pkin(9) * t579 - t574 * t728 + t575 * t733 + t678 * t709 - t679 * t721;
t561 = -mrSges(4,1) * t673 + mrSges(4,3) * t654 + Ifges(4,4) * t685 + Ifges(4,2) * t684 + Ifges(4,6) * t708 - pkin(3) * t740 + pkin(9) * t746 + t733 * t574 + t728 * t575 - t710 * t678 + t721 * t680;
t560 = mrSges(3,2) * t702 - mrSges(3,3) * t692 + Ifges(3,1) * t713 + Ifges(3,4) * t714 + Ifges(3,5) * qJDD(2) - pkin(8) * t573 - qJD(2) * t700 - t561 * t729 + t562 * t734 + t699 * t752;
t559 = Ifges(2,6) * qJDD(1) + t738 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t719 - Ifges(3,5) * t713 - Ifges(3,6) * t714 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t692 + mrSges(3,2) * t693 - t729 * t562 - t734 * t561 - pkin(2) * t739 - pkin(8) * t747 - pkin(1) * t567 + (-t700 * t730 + t701 * t735) * qJD(1);
t558 = -mrSges(2,2) * g(3) - mrSges(2,3) * t718 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t738 - pkin(7) * t567 + t560 * t735 - t563 * t730;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t754; (-m(1) - m(2)) * g(3) + t567; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t754 + t736 * t558 - t731 * t559; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t749 + t731 * t558 + t736 * t559; -mrSges(1,1) * g(2) + mrSges(2,1) * t718 + mrSges(1,2) * g(1) - mrSges(2,2) * t719 + Ifges(2,3) * qJDD(1) + pkin(1) * t741 + pkin(7) * t748 + t730 * t560 + t735 * t563;];
tauB  = t1;
