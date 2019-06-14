% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:52:04
% EndTime: 2019-05-06 08:52:12
% DurationCPUTime: 4.43s
% Computational Cost: add. (40135->350), mult. (88438->405), div. (0->0), fcn. (54306->8), ass. (0->137)
t751 = sin(pkin(9));
t753 = sin(qJ(2));
t790 = qJD(1) * t753;
t798 = cos(pkin(9));
t730 = t751 * qJD(2) + t798 * t790;
t808 = -0.2e1 * t730;
t807 = 2 * qJD(4);
t806 = Ifges(4,1) + Ifges(6,3) + Ifges(5,2);
t805 = Ifges(5,1) + Ifges(6,2) + Ifges(4,3);
t785 = Ifges(6,4) - Ifges(5,5) + Ifges(4,6);
t784 = -Ifges(4,5) + Ifges(6,6) + Ifges(5,4);
t783 = -Ifges(5,6) + Ifges(6,5) - Ifges(4,4);
t804 = -Ifges(5,3) - Ifges(6,1) - Ifges(4,2);
t754 = sin(qJ(1));
t757 = cos(qJ(1));
t741 = -g(1) * t757 - g(2) * t754;
t759 = qJD(1) ^ 2;
t724 = -pkin(1) * t759 + qJDD(1) * pkin(7) + t741;
t756 = cos(qJ(2));
t694 = -t756 * g(3) - t753 * t724;
t734 = (-pkin(2) * t756 - qJ(3) * t753) * qJD(1);
t758 = qJD(2) ^ 2;
t766 = qJDD(2) * pkin(2) + t758 * qJ(3) - t734 * t790 - qJDD(3) + t694;
t729 = -t798 * qJD(2) + t751 * t790;
t789 = qJD(1) * t756;
t779 = t729 * t789;
t803 = -qJ(4) * t779 - t766;
t690 = pkin(3) * t729 - qJ(4) * t730;
t786 = qJD(1) * qJD(2);
t745 = t753 * t786;
t737 = qJDD(1) * t756 - t745;
t740 = t754 * g(1) - t757 * g(2);
t723 = -qJDD(1) * pkin(1) - t759 * pkin(7) - t740;
t777 = t756 * t786;
t736 = qJDD(1) * t753 + t777;
t662 = (-t736 - t777) * qJ(3) + (-t737 + t745) * pkin(2) + t723;
t695 = -g(3) * t753 + t756 * t724;
t666 = -pkin(2) * t758 + qJDD(2) * qJ(3) + t734 * t789 + t695;
t792 = t751 * t662 + t798 * t666;
t795 = t756 ^ 2 * t759;
t802 = pkin(3) * t795 + t737 * qJ(4) + t729 * t690 + t789 * t807 - t792;
t649 = qJD(3) * t808 + t798 * t662 - t751 * t666;
t801 = 2 * qJD(5);
t726 = t729 ^ 2;
t800 = pkin(4) * t726;
t799 = mrSges(5,2) - mrSges(6,3);
t711 = t751 * qJDD(2) + t798 * t736;
t797 = t711 * qJ(4);
t796 = t729 * t730;
t735 = (-mrSges(3,1) * t756 + mrSges(3,2) * t753) * qJD(1);
t738 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t790;
t691 = mrSges(4,1) * t729 + mrSges(4,2) * t730;
t707 = -mrSges(6,1) * t789 + mrSges(6,2) * t729;
t708 = mrSges(4,2) * t789 - mrSges(4,3) * t729;
t647 = t737 * pkin(3) - qJ(4) * t795 + t730 * t690 + qJDD(4) - t649;
t692 = -mrSges(5,2) * t729 - mrSges(5,3) * t730;
t705 = mrSges(5,1) * t729 + mrSges(5,3) * t789;
t641 = t789 * t801 + (t737 + t796) * qJ(5) + (t711 - t779) * pkin(4) + t647;
t689 = mrSges(6,1) * t730 - mrSges(6,3) * t729;
t712 = -pkin(5) * t789 - pkin(8) * t729;
t727 = t730 ^ 2;
t637 = -t727 * pkin(5) + t711 * pkin(8) + t712 * t789 + t641;
t703 = pkin(4) * t730 + qJ(5) * t789;
t710 = -t798 * qJDD(2) + t736 * t751;
t788 = qJD(3) * t729;
t718 = -0.2e1 * t788;
t643 = -t710 * pkin(4) - t726 * qJ(5) - t703 * t789 + qJDD(5) + t718 - t802;
t778 = t730 * t789;
t638 = t643 + (-t710 - t778) * pkin(8) + (-t737 + t796) * pkin(5);
t752 = sin(qJ(6));
t755 = cos(qJ(6));
t635 = -t637 * t752 + t638 * t755;
t687 = -t729 * t752 + t730 * t755;
t654 = qJD(6) * t687 + t710 * t755 + t711 * t752;
t688 = t729 * t755 + t730 * t752;
t660 = -mrSges(7,1) * t687 + mrSges(7,2) * t688;
t742 = qJD(6) - t789;
t667 = -mrSges(7,2) * t742 + mrSges(7,3) * t687;
t733 = qJDD(6) - t737;
t633 = m(7) * t635 + mrSges(7,1) * t733 - mrSges(7,3) * t654 - t660 * t688 + t667 * t742;
t636 = t637 * t755 + t638 * t752;
t653 = -qJD(6) * t688 - t710 * t752 + t711 * t755;
t668 = mrSges(7,1) * t742 - mrSges(7,3) * t688;
t634 = m(7) * t636 - mrSges(7,2) * t733 + mrSges(7,3) * t653 + t660 * t687 - t668 * t742;
t793 = -t752 * t633 + t755 * t634;
t769 = m(6) * t641 - t711 * mrSges(6,2) - t730 * t689 + t793;
t764 = -m(5) * t647 - t711 * mrSges(5,1) - t730 * t692 + t705 * t789 - t769;
t619 = m(4) * t649 - mrSges(4,3) * t711 - t691 * t730 + (-t707 - t708) * t789 + (-mrSges(4,1) + t799) * t737 + t764;
t650 = t718 + t792;
t704 = -mrSges(6,2) * t730 + mrSges(6,3) * t789;
t646 = 0.2e1 * t788 + t802;
t623 = t755 * t633 + t752 * t634;
t771 = -m(6) * t643 - t710 * mrSges(6,2) - t729 * t689 - t623;
t767 = -m(5) * t646 - t737 * mrSges(5,3) - t771;
t706 = mrSges(5,1) * t730 - mrSges(5,2) * t789;
t791 = mrSges(4,1) * t789 + mrSges(4,3) * t730 + t706;
t621 = m(4) * t650 + (-mrSges(6,1) + mrSges(4,2)) * t737 + (-t691 - t692) * t729 + (-mrSges(4,3) - mrSges(5,1)) * t710 + (-t704 - t791) * t789 + t767;
t773 = -t619 * t751 + t798 * t621;
t617 = m(3) * t695 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t737 - qJD(2) * t738 + t735 * t789 + t773;
t739 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t789;
t761 = qJD(4) * t808 + (t710 - t778) * pkin(3) + t803;
t648 = t761 - t797;
t645 = t800 + t797 - 0.2e1 * qJD(5) * t729 + (-pkin(3) - qJ(5)) * t710 + (pkin(3) * t789 + t703 + t807) * t730 - t803;
t640 = t761 + (t801 + t712) * t729 + (-pkin(5) - qJ(4)) * t711 - pkin(8) * t727 - t800 + qJ(5) * t710 - t703 * t730;
t770 = -m(7) * t640 + t653 * mrSges(7,1) - t654 * mrSges(7,2) + t687 * t667 - t688 * t668;
t765 = -m(6) * t645 - t711 * mrSges(6,1) + t710 * mrSges(6,3) - t730 * t704 + t729 * t707 - t770;
t763 = -m(5) * t648 + t710 * mrSges(5,2) + t729 * t705 - t765;
t760 = (-mrSges(4,2) + mrSges(5,3)) * t711 + t791 * t730 + m(4) * t766 - t710 * mrSges(4,1) - t729 * t708 + t763;
t625 = m(3) * t694 + qJDD(2) * mrSges(3,1) - t736 * mrSges(3,3) + qJD(2) * t739 - t735 * t790 + t760;
t774 = t756 * t617 - t625 * t753;
t611 = m(2) * t741 - mrSges(2,1) * t759 - qJDD(1) * mrSges(2,2) + t774;
t618 = t798 * t619 + t751 * t621;
t762 = -m(3) * t723 + t737 * mrSges(3,1) - t736 * mrSges(3,2) - t738 * t790 + t739 * t789 - t618;
t614 = m(2) * t740 + qJDD(1) * mrSges(2,1) - t759 * mrSges(2,2) + t762;
t794 = t754 * t611 + t757 * t614;
t612 = t753 * t617 + t756 * t625;
t782 = -t729 * t804 + t730 * t783 + t785 * t789;
t781 = -t729 * t783 - t730 * t806 - t784 * t789;
t780 = t729 * t785 + t730 * t784 + t789 * t805;
t775 = t757 * t611 - t614 * t754;
t722 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t753 + Ifges(3,4) * t756) * qJD(1);
t721 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t753 + Ifges(3,2) * t756) * qJD(1);
t720 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t753 + Ifges(3,6) * t756) * qJD(1);
t657 = Ifges(7,1) * t688 + Ifges(7,4) * t687 + Ifges(7,5) * t742;
t656 = Ifges(7,4) * t688 + Ifges(7,2) * t687 + Ifges(7,6) * t742;
t655 = Ifges(7,5) * t688 + Ifges(7,6) * t687 + Ifges(7,3) * t742;
t628 = -t711 * mrSges(5,3) - t730 * t706 - t763;
t627 = mrSges(7,2) * t640 - mrSges(7,3) * t635 + Ifges(7,1) * t654 + Ifges(7,4) * t653 + Ifges(7,5) * t733 + t655 * t687 - t656 * t742;
t626 = -mrSges(7,1) * t640 + mrSges(7,3) * t636 + Ifges(7,4) * t654 + Ifges(7,2) * t653 + Ifges(7,6) * t733 - t655 * t688 + t657 * t742;
t622 = mrSges(6,3) * t737 + t707 * t789 + t769;
t608 = t755 * t626 + t752 * t627 - mrSges(4,2) * t766 - mrSges(4,3) * t649 + pkin(5) * t770 + mrSges(6,1) * t645 + mrSges(5,1) * t647 - mrSges(5,3) * t648 - mrSges(6,2) * t641 + pkin(8) * t793 - qJ(4) * t628 + pkin(4) * t622 + t784 * t737 + t780 * t729 + t806 * t711 + t783 * t710 - t782 * t789;
t607 = -t755 * t627 - pkin(4) * t771 + t752 * t626 - qJ(5) * t765 + mrSges(4,1) * t766 + mrSges(4,3) * t650 + mrSges(6,3) * t645 - mrSges(5,1) * t646 + mrSges(5,2) * t648 - mrSges(6,2) * t643 - pkin(3) * t628 + pkin(8) * t623 + (-mrSges(6,1) * pkin(4) - t785) * t737 + t780 * t730 - t783 * t711 + t804 * t710 + (-pkin(4) * t704 + t781) * t789;
t606 = -pkin(3) * t764 - qJ(4) * t767 + (mrSges(5,1) * qJ(4) + t785) * t710 + (qJ(4) * t692 + t781) * t729 + t782 * t730 + t784 * t711 + (qJ(4) * mrSges(6,1) - pkin(3) * t799 + Ifges(3,2) + t805) * t737 - Ifges(7,3) * t733 + Ifges(3,4) * t736 + qJD(2) * t722 - mrSges(3,1) * t723 + t687 * t657 - t688 * t656 + mrSges(3,3) * t695 - mrSges(4,1) * t649 + mrSges(4,2) * t650 - Ifges(7,6) * t653 - Ifges(7,5) * t654 + mrSges(5,3) * t646 - mrSges(5,2) * t647 + mrSges(6,3) * t641 - mrSges(6,1) * t643 - mrSges(7,1) * t635 + mrSges(7,2) * t636 - pkin(5) * t623 + qJ(5) * t622 + Ifges(3,6) * qJDD(2) - pkin(2) * t618 + (-t753 * t720 + (pkin(3) * t707 - qJ(4) * (-t704 - t706)) * t756) * qJD(1);
t605 = mrSges(3,2) * t723 - mrSges(3,3) * t694 + Ifges(3,1) * t736 + Ifges(3,4) * t737 + Ifges(3,5) * qJDD(2) - qJ(3) * t618 - qJD(2) * t721 - t751 * t607 + t798 * t608 + t720 * t789;
t604 = Ifges(2,6) * qJDD(1) + t759 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t741 - Ifges(3,5) * t736 - Ifges(3,6) * t737 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t694 + mrSges(3,2) * t695 - t751 * t608 - t798 * t607 - pkin(2) * t760 - qJ(3) * t773 - pkin(1) * t612 + (-t721 * t753 + t722 * t756) * qJD(1);
t603 = -mrSges(2,2) * g(3) - mrSges(2,3) * t740 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t759 - pkin(7) * t612 + t605 * t756 - t606 * t753;
t1 = [-m(1) * g(1) + t775; -m(1) * g(2) + t794; (-m(1) - m(2)) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t794 + t757 * t603 - t754 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t775 + t754 * t603 + t757 * t604; -mrSges(1,1) * g(2) + mrSges(2,1) * t740 + mrSges(1,2) * g(1) - mrSges(2,2) * t741 + Ifges(2,3) * qJDD(1) + pkin(1) * t762 + pkin(7) * t774 + t753 * t605 + t756 * t606;];
tauB  = t1;
