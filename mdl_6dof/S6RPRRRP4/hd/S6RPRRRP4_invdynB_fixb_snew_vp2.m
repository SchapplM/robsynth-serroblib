% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:24:23
% EndTime: 2019-05-06 01:24:35
% DurationCPUTime: 12.19s
% Computational Cost: add. (173474->343), mult. (416102->414), div. (0->0), fcn. (324466->10), ass. (0->141)
t758 = Ifges(6,1) + Ifges(7,1);
t752 = Ifges(6,4) + Ifges(7,4);
t751 = Ifges(6,5) + Ifges(7,5);
t757 = Ifges(6,2) + Ifges(7,2);
t756 = Ifges(6,6) + Ifges(7,6);
t755 = Ifges(6,3) + Ifges(7,3);
t718 = qJD(1) ^ 2;
t709 = cos(pkin(10));
t754 = pkin(2) * t709;
t753 = -mrSges(6,2) - mrSges(7,2);
t708 = sin(pkin(10));
t749 = mrSges(3,2) * t708;
t706 = t709 ^ 2;
t748 = t706 * t718;
t713 = sin(qJ(1));
t717 = cos(qJ(1));
t696 = -g(1) * t717 - g(2) * t713;
t692 = -pkin(1) * t718 + qJDD(1) * qJ(2) + t696;
t740 = qJD(1) * qJD(2);
t736 = -t709 * g(3) - 0.2e1 * t708 * t740;
t667 = (-pkin(7) * qJDD(1) + t718 * t754 - t692) * t708 + t736;
t683 = -g(3) * t708 + (t692 + 0.2e1 * t740) * t709;
t739 = qJDD(1) * t709;
t668 = -pkin(2) * t748 + pkin(7) * t739 + t683;
t712 = sin(qJ(3));
t716 = cos(qJ(3));
t645 = t716 * t667 - t712 * t668;
t725 = t708 * t716 + t709 * t712;
t724 = -t708 * t712 + t709 * t716;
t690 = t724 * qJD(1);
t741 = t690 * qJD(3);
t681 = qJDD(1) * t725 + t741;
t691 = t725 * qJD(1);
t615 = (-t681 + t741) * pkin(8) + (t690 * t691 + qJDD(3)) * pkin(3) + t645;
t646 = t712 * t667 + t716 * t668;
t680 = -t691 * qJD(3) + qJDD(1) * t724;
t686 = qJD(3) * pkin(3) - pkin(8) * t691;
t689 = t690 ^ 2;
t623 = -pkin(3) * t689 + pkin(8) * t680 - qJD(3) * t686 + t646;
t711 = sin(qJ(4));
t715 = cos(qJ(4));
t613 = t711 * t615 + t715 * t623;
t674 = t690 * t711 + t691 * t715;
t641 = -qJD(4) * t674 + t680 * t715 - t681 * t711;
t673 = t690 * t715 - t691 * t711;
t657 = -mrSges(5,1) * t673 + mrSges(5,2) * t674;
t707 = qJD(3) + qJD(4);
t665 = mrSges(5,1) * t707 - mrSges(5,3) * t674;
t704 = qJDD(3) + qJDD(4);
t658 = -pkin(4) * t673 - pkin(9) * t674;
t703 = t707 ^ 2;
t608 = -pkin(4) * t703 + pkin(9) * t704 + t658 * t673 + t613;
t705 = t708 ^ 2;
t695 = t713 * g(1) - t717 * g(2);
t729 = qJDD(2) - t695;
t679 = (-pkin(1) - t754) * qJDD(1) + (-qJ(2) + (-t705 - t706) * pkin(7)) * t718 + t729;
t635 = -t680 * pkin(3) - t689 * pkin(8) + t691 * t686 + t679;
t642 = qJD(4) * t673 + t680 * t711 + t681 * t715;
t611 = (-t673 * t707 - t642) * pkin(9) + (t674 * t707 - t641) * pkin(4) + t635;
t710 = sin(qJ(5));
t714 = cos(qJ(5));
t603 = -t710 * t608 + t714 * t611;
t660 = -t674 * t710 + t707 * t714;
t620 = qJD(5) * t660 + t642 * t714 + t704 * t710;
t640 = qJDD(5) - t641;
t661 = t674 * t714 + t707 * t710;
t643 = -mrSges(7,1) * t660 + mrSges(7,2) * t661;
t644 = -mrSges(6,1) * t660 + mrSges(6,2) * t661;
t669 = qJD(5) - t673;
t648 = -mrSges(6,2) * t669 + mrSges(6,3) * t660;
t600 = -0.2e1 * qJD(6) * t661 + (t660 * t669 - t620) * qJ(6) + (t660 * t661 + t640) * pkin(5) + t603;
t647 = -mrSges(7,2) * t669 + mrSges(7,3) * t660;
t738 = m(7) * t600 + t640 * mrSges(7,1) + t669 * t647;
t592 = m(6) * t603 + t640 * mrSges(6,1) + t669 * t648 + (-t643 - t644) * t661 + (-mrSges(6,3) - mrSges(7,3)) * t620 + t738;
t604 = t714 * t608 + t710 * t611;
t619 = -qJD(5) * t661 - t642 * t710 + t704 * t714;
t649 = pkin(5) * t669 - qJ(6) * t661;
t659 = t660 ^ 2;
t602 = -pkin(5) * t659 + qJ(6) * t619 + 0.2e1 * qJD(6) * t660 - t649 * t669 + t604;
t737 = m(7) * t602 + t619 * mrSges(7,3) + t660 * t643;
t650 = mrSges(7,1) * t669 - mrSges(7,3) * t661;
t743 = -mrSges(6,1) * t669 + mrSges(6,3) * t661 - t650;
t595 = m(6) * t604 + t619 * mrSges(6,3) + t640 * t753 + t660 * t644 + t669 * t743 + t737;
t731 = -t592 * t710 + t714 * t595;
t588 = m(5) * t613 - mrSges(5,2) * t704 + mrSges(5,3) * t641 + t657 * t673 - t665 * t707 + t731;
t612 = t615 * t715 - t711 * t623;
t664 = -mrSges(5,2) * t707 + mrSges(5,3) * t673;
t607 = -pkin(4) * t704 - pkin(9) * t703 + t674 * t658 - t612;
t605 = -pkin(5) * t619 - qJ(6) * t659 + t649 * t661 + qJDD(6) + t607;
t730 = m(7) * t605 - t619 * mrSges(7,1) - t660 * t647;
t720 = -m(6) * t607 + t619 * mrSges(6,1) + t620 * t753 + t660 * t648 + t661 * t743 - t730;
t597 = m(5) * t612 + t704 * mrSges(5,1) - t642 * mrSges(5,3) - t674 * t657 + t707 * t664 + t720;
t582 = t711 * t588 + t715 * t597;
t677 = -mrSges(4,1) * t690 + mrSges(4,2) * t691;
t684 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t690;
t580 = m(4) * t645 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t681 + qJD(3) * t684 - t677 * t691 + t582;
t685 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t691;
t732 = t715 * t588 - t597 * t711;
t581 = m(4) * t646 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t680 - qJD(3) * t685 + t677 * t690 + t732;
t575 = t716 * t580 + t712 * t581;
t682 = -t708 * t692 + t736;
t723 = mrSges(3,3) * qJDD(1) + t718 * (-mrSges(3,1) * t709 + t749);
t573 = m(3) * t682 - t708 * t723 + t575;
t733 = -t712 * t580 + t716 * t581;
t574 = m(3) * t683 + t709 * t723 + t733;
t734 = -t573 * t708 + t709 * t574;
t567 = m(2) * t696 - mrSges(2,1) * t718 - qJDD(1) * mrSges(2,2) + t734;
t688 = -qJDD(1) * pkin(1) - t718 * qJ(2) + t729;
t590 = t714 * t592 + t710 * t595;
t722 = m(5) * t635 - t641 * mrSges(5,1) + t642 * mrSges(5,2) - t673 * t664 + t674 * t665 + t590;
t721 = m(4) * t679 - t680 * mrSges(4,1) + t681 * mrSges(4,2) - t690 * t684 + t691 * t685 + t722;
t719 = -m(3) * t688 + mrSges(3,1) * t739 - t721 + (t705 * t718 + t748) * mrSges(3,3);
t585 = t719 + (mrSges(2,1) - t749) * qJDD(1) + m(2) * t695 - t718 * mrSges(2,2);
t747 = t713 * t567 + t717 * t585;
t568 = t709 * t573 + t708 * t574;
t746 = t660 * t756 + t661 * t751 + t669 * t755;
t745 = -t660 * t757 - t661 * t752 - t669 * t756;
t744 = t752 * t660 + t661 * t758 + t751 * t669;
t726 = Ifges(3,5) * t708 + Ifges(3,6) * t709;
t742 = t718 * t726;
t735 = t717 * t567 - t585 * t713;
t728 = Ifges(3,1) * t708 + Ifges(3,4) * t709;
t727 = Ifges(3,4) * t708 + Ifges(3,2) * t709;
t672 = Ifges(4,1) * t691 + Ifges(4,4) * t690 + Ifges(4,5) * qJD(3);
t671 = Ifges(4,4) * t691 + Ifges(4,2) * t690 + Ifges(4,6) * qJD(3);
t670 = Ifges(4,5) * t691 + Ifges(4,6) * t690 + Ifges(4,3) * qJD(3);
t654 = Ifges(5,1) * t674 + Ifges(5,4) * t673 + Ifges(5,5) * t707;
t653 = Ifges(5,4) * t674 + Ifges(5,2) * t673 + Ifges(5,6) * t707;
t652 = Ifges(5,5) * t674 + Ifges(5,6) * t673 + Ifges(5,3) * t707;
t598 = -t620 * mrSges(7,3) - t661 * t643 + t738;
t589 = mrSges(6,2) * t607 + mrSges(7,2) * t605 - mrSges(6,3) * t603 - mrSges(7,3) * t600 - qJ(6) * t598 + t752 * t619 + t620 * t758 + t751 * t640 + t746 * t660 + t745 * t669;
t583 = -mrSges(6,1) * t607 + mrSges(6,3) * t604 - mrSges(7,1) * t605 + mrSges(7,3) * t602 - pkin(5) * t730 + qJ(6) * t737 + (-qJ(6) * t650 + t744) * t669 + (-pkin(5) * t650 - t746) * t661 + (-mrSges(7,2) * qJ(6) + t756) * t640 + (-mrSges(7,2) * pkin(5) + t752) * t620 + t757 * t619;
t576 = -mrSges(5,1) * t635 - mrSges(6,1) * t603 - mrSges(7,1) * t600 + mrSges(6,2) * t604 + mrSges(7,2) * t602 + mrSges(5,3) * t613 + Ifges(5,4) * t642 + Ifges(5,2) * t641 + Ifges(5,6) * t704 - pkin(4) * t590 - pkin(5) * t598 - t674 * t652 + t707 * t654 + t745 * t661 + t744 * t660 - t755 * t640 - t751 * t620 - t756 * t619;
t569 = mrSges(5,2) * t635 - mrSges(5,3) * t612 + Ifges(5,1) * t642 + Ifges(5,4) * t641 + Ifges(5,5) * t704 - pkin(9) * t590 - t583 * t710 + t589 * t714 + t652 * t673 - t653 * t707;
t564 = mrSges(4,2) * t679 - mrSges(4,3) * t645 + Ifges(4,1) * t681 + Ifges(4,4) * t680 + Ifges(4,5) * qJDD(3) - pkin(8) * t582 - qJD(3) * t671 + t569 * t715 - t576 * t711 + t670 * t690;
t563 = -mrSges(4,1) * t679 + mrSges(4,3) * t646 + Ifges(4,4) * t681 + Ifges(4,2) * t680 + Ifges(4,6) * qJDD(3) - pkin(3) * t722 + pkin(8) * t732 + qJD(3) * t672 + t711 * t569 + t715 * t576 - t691 * t670;
t562 = -pkin(9) * t731 + (-t708 * t727 + t709 * t728 + Ifges(2,5)) * t718 + (Ifges(2,6) - t726) * qJDD(1) - pkin(3) * t582 - pkin(2) * t575 - Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) - Ifges(5,3) * t704 - pkin(1) * t568 - mrSges(5,1) * t612 + mrSges(5,2) * t613 - Ifges(4,6) * t680 - Ifges(4,5) * t681 - mrSges(3,1) * t682 + mrSges(3,2) * t683 + t690 * t672 - t691 * t671 + mrSges(2,3) * t696 - pkin(4) * t720 - t710 * t589 - t714 * t583 + t673 * t654 - t674 * t653 - Ifges(5,6) * t641 - Ifges(5,5) * t642 - mrSges(4,1) * t645 + mrSges(4,2) * t646;
t561 = mrSges(3,2) * t688 - mrSges(3,3) * t682 - pkin(7) * t575 + qJDD(1) * t728 - t712 * t563 + t716 * t564 + t709 * t742;
t560 = -mrSges(3,1) * t688 + mrSges(3,3) * t683 - pkin(2) * t721 + pkin(7) * t733 + qJDD(1) * t727 + t716 * t563 + t712 * t564 - t708 * t742;
t559 = -mrSges(2,2) * g(3) - mrSges(2,3) * t695 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t718 - qJ(2) * t568 - t560 * t708 + t561 * t709;
t1 = [-m(1) * g(1) + t735; -m(1) * g(2) + t747; (-m(1) - m(2)) * g(3) + t568; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t747 + t717 * t559 - t713 * t562; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t735 + t713 * t559 + t717 * t562; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t695 - mrSges(2,2) * t696 + t708 * t561 + t709 * t560 + pkin(1) * (-qJDD(1) * t749 + t719) + qJ(2) * t734;];
tauB  = t1;
