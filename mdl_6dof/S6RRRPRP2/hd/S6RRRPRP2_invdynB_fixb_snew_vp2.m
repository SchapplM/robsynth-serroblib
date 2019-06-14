% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:31:21
% EndTime: 2019-05-07 07:31:37
% DurationCPUTime: 13.40s
% Computational Cost: add. (193669->363), mult. (434940->448), div. (0->0), fcn. (318068->10), ass. (0->139)
t763 = Ifges(6,1) + Ifges(7,1);
t756 = Ifges(6,4) - Ifges(7,5);
t762 = -Ifges(6,5) - Ifges(7,4);
t761 = Ifges(6,2) + Ifges(7,3);
t754 = Ifges(6,6) - Ifges(7,6);
t760 = -Ifges(6,3) - Ifges(7,2);
t725 = sin(qJ(2));
t728 = cos(qJ(2));
t744 = qJD(1) * qJD(2);
t707 = qJDD(1) * t725 + t728 * t744;
t726 = sin(qJ(1));
t729 = cos(qJ(1));
t713 = -g(1) * t729 - g(2) * t726;
t730 = qJD(1) ^ 2;
t702 = -pkin(1) * t730 + qJDD(1) * pkin(7) + t713;
t753 = t702 * t725;
t758 = pkin(2) * t730;
t667 = qJDD(2) * pkin(2) - pkin(8) * t707 - t753 + (pkin(8) * t744 + t725 * t758 - g(3)) * t728;
t690 = -g(3) * t725 + t728 * t702;
t708 = qJDD(1) * t728 - t725 * t744;
t747 = qJD(1) * t725;
t711 = qJD(2) * pkin(2) - pkin(8) * t747;
t720 = t728 ^ 2;
t668 = pkin(8) * t708 - qJD(2) * t711 - t720 * t758 + t690;
t724 = sin(qJ(3));
t727 = cos(qJ(3));
t641 = t727 * t667 - t668 * t724;
t699 = (-t724 * t725 + t727 * t728) * qJD(1);
t673 = qJD(3) * t699 + t707 * t727 + t708 * t724;
t700 = (t724 * t728 + t725 * t727) * qJD(1);
t718 = qJDD(2) + qJDD(3);
t719 = qJD(2) + qJD(3);
t620 = (t699 * t719 - t673) * qJ(4) + (t699 * t700 + t718) * pkin(3) + t641;
t642 = t724 * t667 + t727 * t668;
t672 = -qJD(3) * t700 - t707 * t724 + t708 * t727;
t692 = pkin(3) * t719 - qJ(4) * t700;
t695 = t699 ^ 2;
t623 = -pkin(3) * t695 + qJ(4) * t672 - t692 * t719 + t642;
t721 = sin(pkin(10));
t722 = cos(pkin(10));
t687 = t699 * t721 + t700 * t722;
t615 = -0.2e1 * qJD(4) * t687 + t620 * t722 - t721 * t623;
t759 = cos(qJ(5));
t757 = -mrSges(6,3) - mrSges(7,2);
t686 = t699 * t722 - t700 * t721;
t616 = 0.2e1 * qJD(4) * t686 + t721 * t620 + t722 * t623;
t648 = t672 * t722 - t673 * t721;
t662 = -mrSges(5,1) * t686 + mrSges(5,2) * t687;
t676 = mrSges(5,1) * t719 - mrSges(5,3) * t687;
t663 = -pkin(4) * t686 - pkin(9) * t687;
t717 = t719 ^ 2;
t614 = -pkin(4) * t717 + pkin(9) * t718 + t663 * t686 + t616;
t712 = g(1) * t726 - t729 * g(2);
t735 = -qJDD(1) * pkin(1) - t712;
t674 = -pkin(2) * t708 + t711 * t747 + (-pkin(8) * t720 - pkin(7)) * t730 + t735;
t629 = -pkin(3) * t672 - qJ(4) * t695 + t700 * t692 + qJDD(4) + t674;
t649 = t672 * t721 + t673 * t722;
t618 = (-t686 * t719 - t649) * pkin(9) + (t687 * t719 - t648) * pkin(4) + t629;
t723 = sin(qJ(5));
t611 = t759 * t614 + t723 * t618;
t671 = t759 * t687 + t723 * t719;
t626 = qJD(5) * t671 + t649 * t723 - t759 * t718;
t647 = qJDD(5) - t648;
t680 = qJD(5) - t686;
t655 = mrSges(6,1) * t680 - mrSges(6,3) * t671;
t670 = t687 * t723 - t759 * t719;
t650 = pkin(5) * t670 - qJ(6) * t671;
t679 = t680 ^ 2;
t607 = -pkin(5) * t679 + qJ(6) * t647 + 0.2e1 * qJD(6) * t680 - t650 * t670 + t611;
t656 = -mrSges(7,1) * t680 + mrSges(7,2) * t671;
t743 = m(7) * t607 + t647 * mrSges(7,3) + t680 * t656;
t651 = mrSges(7,1) * t670 - mrSges(7,3) * t671;
t748 = -mrSges(6,1) * t670 - mrSges(6,2) * t671 - t651;
t602 = m(6) * t611 - mrSges(6,2) * t647 + t757 * t626 - t655 * t680 + t748 * t670 + t743;
t610 = -t723 * t614 + t759 * t618;
t627 = -t670 * qJD(5) + t759 * t649 + t723 * t718;
t654 = -mrSges(6,2) * t680 - mrSges(6,3) * t670;
t608 = -t647 * pkin(5) - t679 * qJ(6) + t671 * t650 + qJDD(6) - t610;
t653 = -mrSges(7,2) * t670 + mrSges(7,3) * t680;
t736 = -m(7) * t608 + t647 * mrSges(7,1) + t680 * t653;
t604 = m(6) * t610 + mrSges(6,1) * t647 + t757 * t627 + t654 * t680 + t748 * t671 + t736;
t738 = t759 * t602 - t604 * t723;
t594 = m(5) * t616 - mrSges(5,2) * t718 + mrSges(5,3) * t648 + t662 * t686 - t676 * t719 + t738;
t675 = -mrSges(5,2) * t719 + mrSges(5,3) * t686;
t613 = -pkin(4) * t718 - pkin(9) * t717 + t687 * t663 - t615;
t609 = -0.2e1 * qJD(6) * t671 + (t670 * t680 - t627) * qJ(6) + (t671 * t680 + t626) * pkin(5) + t613;
t605 = m(7) * t609 + mrSges(7,1) * t626 - t627 * mrSges(7,3) + t653 * t670 - t671 * t656;
t732 = -m(6) * t613 - t626 * mrSges(6,1) - mrSges(6,2) * t627 - t670 * t654 - t655 * t671 - t605;
t599 = m(5) * t615 + mrSges(5,1) * t718 - mrSges(5,3) * t649 - t662 * t687 + t675 * t719 + t732;
t589 = t721 * t594 + t722 * t599;
t688 = -mrSges(4,1) * t699 + mrSges(4,2) * t700;
t691 = -mrSges(4,2) * t719 + mrSges(4,3) * t699;
t587 = m(4) * t641 + mrSges(4,1) * t718 - mrSges(4,3) * t673 - t688 * t700 + t691 * t719 + t589;
t693 = mrSges(4,1) * t719 - mrSges(4,3) * t700;
t739 = t722 * t594 - t599 * t721;
t588 = m(4) * t642 - mrSges(4,2) * t718 + mrSges(4,3) * t672 + t688 * t699 - t693 * t719 + t739;
t581 = t727 * t587 + t724 * t588;
t689 = -g(3) * t728 - t753;
t706 = (-mrSges(3,1) * t728 + mrSges(3,2) * t725) * qJD(1);
t746 = qJD(1) * t728;
t710 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t746;
t579 = m(3) * t689 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t707 + qJD(2) * t710 - t706 * t747 + t581;
t709 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t747;
t740 = -t587 * t724 + t727 * t588;
t580 = m(3) * t690 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t708 - qJD(2) * t709 + t706 * t746 + t740;
t741 = -t579 * t725 + t728 * t580;
t574 = m(2) * t713 - mrSges(2,1) * t730 - qJDD(1) * mrSges(2,2) + t741;
t701 = -pkin(7) * t730 + t735;
t597 = t723 * t602 + t759 * t604;
t734 = m(5) * t629 - t648 * mrSges(5,1) + t649 * mrSges(5,2) - t686 * t675 + t687 * t676 + t597;
t733 = m(4) * t674 - t672 * mrSges(4,1) + mrSges(4,2) * t673 - t699 * t691 + t693 * t700 + t734;
t731 = -m(3) * t701 + t708 * mrSges(3,1) - mrSges(3,2) * t707 - t709 * t747 + t710 * t746 - t733;
t591 = m(2) * t712 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t730 + t731;
t752 = t726 * t574 + t729 * t591;
t575 = t728 * t579 + t725 * t580;
t751 = t761 * t670 - t756 * t671 - t754 * t680;
t750 = t754 * t670 + t762 * t671 + t760 * t680;
t749 = -t756 * t670 + t763 * t671 - t762 * t680;
t742 = t729 * t574 - t591 * t726;
t698 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t725 + Ifges(3,4) * t728) * qJD(1);
t697 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t725 + Ifges(3,2) * t728) * qJD(1);
t696 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t725 + Ifges(3,6) * t728) * qJD(1);
t683 = Ifges(4,1) * t700 + Ifges(4,4) * t699 + Ifges(4,5) * t719;
t682 = Ifges(4,4) * t700 + Ifges(4,2) * t699 + Ifges(4,6) * t719;
t681 = Ifges(4,5) * t700 + Ifges(4,6) * t699 + Ifges(4,3) * t719;
t659 = Ifges(5,1) * t687 + Ifges(5,4) * t686 + Ifges(5,5) * t719;
t658 = Ifges(5,4) * t687 + Ifges(5,2) * t686 + Ifges(5,6) * t719;
t657 = Ifges(5,5) * t687 + Ifges(5,6) * t686 + Ifges(5,3) * t719;
t596 = mrSges(6,2) * t613 + mrSges(7,2) * t608 - mrSges(6,3) * t610 - mrSges(7,3) * t609 - qJ(6) * t605 - t756 * t626 + t763 * t627 - t647 * t762 + t750 * t670 + t751 * t680;
t595 = -mrSges(6,1) * t613 - mrSges(7,1) * t609 + mrSges(7,2) * t607 + mrSges(6,3) * t611 - pkin(5) * t605 - t761 * t626 + t756 * t627 + t754 * t647 + t750 * t671 + t749 * t680;
t583 = Ifges(5,4) * t649 + Ifges(5,2) * t648 + Ifges(5,6) * t718 - t687 * t657 + t719 * t659 - mrSges(5,1) * t629 + mrSges(5,3) * t616 - mrSges(6,1) * t610 + mrSges(6,2) * t611 + mrSges(7,1) * t608 - mrSges(7,3) * t607 - pkin(5) * t736 - qJ(6) * t743 - pkin(4) * t597 + (pkin(5) * t651 + t751) * t671 + (qJ(6) * t651 - t749) * t670 + t760 * t647 + (mrSges(7,2) * pkin(5) + t762) * t627 + (mrSges(7,2) * qJ(6) + t754) * t626;
t582 = mrSges(5,2) * t629 - mrSges(5,3) * t615 + Ifges(5,1) * t649 + Ifges(5,4) * t648 + Ifges(5,5) * t718 - pkin(9) * t597 - t723 * t595 + t759 * t596 + t686 * t657 - t719 * t658;
t571 = mrSges(4,2) * t674 - mrSges(4,3) * t641 + Ifges(4,1) * t673 + Ifges(4,4) * t672 + Ifges(4,5) * t718 - qJ(4) * t589 + t582 * t722 - t583 * t721 + t681 * t699 - t682 * t719;
t570 = -mrSges(4,1) * t674 + mrSges(4,3) * t642 + Ifges(4,4) * t673 + Ifges(4,2) * t672 + Ifges(4,6) * t718 - pkin(3) * t734 + qJ(4) * t739 + t721 * t582 + t722 * t583 - t700 * t681 + t719 * t683;
t569 = -Ifges(5,6) * t648 - Ifges(5,5) * t649 - mrSges(4,1) * t641 + mrSges(4,2) * t642 + mrSges(5,2) * t616 - mrSges(5,1) * t615 - Ifges(3,3) * qJDD(2) - pkin(3) * t589 + mrSges(2,1) * g(3) - pkin(2) * t581 - t723 * t596 + t730 * Ifges(2,5) - pkin(9) * t738 - Ifges(4,6) * t672 - Ifges(4,5) * t673 + t686 * t659 - t687 * t658 - mrSges(3,1) * t689 + mrSges(3,2) * t690 + mrSges(2,3) * t713 - Ifges(3,5) * t707 - Ifges(3,6) * t708 + t699 * t683 - t700 * t682 - pkin(1) * t575 + Ifges(2,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * t718 + (-t697 * t725 + t698 * t728) * qJD(1) - pkin(4) * t732 - t759 * t595;
t568 = mrSges(3,2) * t701 - mrSges(3,3) * t689 + Ifges(3,1) * t707 + Ifges(3,4) * t708 + Ifges(3,5) * qJDD(2) - pkin(8) * t581 - qJD(2) * t697 - t570 * t724 + t571 * t727 + t696 * t746;
t567 = -mrSges(3,1) * t701 + mrSges(3,3) * t690 + Ifges(3,4) * t707 + Ifges(3,2) * t708 + Ifges(3,6) * qJDD(2) - pkin(2) * t733 + pkin(8) * t740 + qJD(2) * t698 + t727 * t570 + t724 * t571 - t696 * t747;
t566 = -mrSges(2,2) * g(3) - mrSges(2,3) * t712 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t730 - pkin(7) * t575 - t567 * t725 + t568 * t728;
t1 = [-m(1) * g(1) + t742; -m(1) * g(2) + t752; (-m(1) - m(2)) * g(3) + t575; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t752 + t729 * t566 - t726 * t569; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t742 + t726 * t566 + t729 * t569; -mrSges(1,1) * g(2) + mrSges(2,1) * t712 + mrSges(1,2) * g(1) - mrSges(2,2) * t713 + Ifges(2,3) * qJDD(1) + pkin(1) * t731 + pkin(7) * t741 + t728 * t567 + t725 * t568;];
tauB  = t1;
