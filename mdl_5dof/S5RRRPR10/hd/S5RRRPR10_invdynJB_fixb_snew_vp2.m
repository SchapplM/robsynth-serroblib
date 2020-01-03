% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:27:51
% EndTime: 2019-12-31 21:28:07
% DurationCPUTime: 15.82s
% Computational Cost: add. (249639->327), mult. (552907->429), div. (0->0), fcn. (432335->12), ass. (0->139)
t784 = sin(pkin(5));
t789 = sin(qJ(2));
t793 = cos(qJ(2));
t807 = qJD(1) * qJD(2);
t771 = (-qJDD(1) * t793 + t789 * t807) * t784;
t809 = qJD(1) * t784;
t769 = (-pkin(2) * t793 - pkin(8) * t789) * t809;
t786 = cos(pkin(5));
t779 = qJD(1) * t786 + qJD(2);
t777 = t779 ^ 2;
t778 = qJDD(1) * t786 + qJDD(2);
t808 = qJD(1) * t793;
t790 = sin(qJ(1));
t794 = cos(qJ(1));
t775 = t790 * g(1) - g(2) * t794;
t795 = qJD(1) ^ 2;
t818 = pkin(7) * t784;
t766 = qJDD(1) * pkin(1) + t795 * t818 + t775;
t776 = -g(1) * t794 - g(2) * t790;
t767 = -pkin(1) * t795 + qJDD(1) * t818 + t776;
t813 = t786 * t789;
t810 = t766 * t813 + t793 * t767;
t728 = -pkin(2) * t777 + pkin(8) * t778 + (-g(3) * t789 + t769 * t808) * t784 + t810;
t770 = (qJDD(1) * t789 + t793 * t807) * t784;
t817 = g(3) * t786;
t729 = pkin(2) * t771 - pkin(8) * t770 - t817 + (-t766 + (pkin(2) * t789 - pkin(8) * t793) * t779 * qJD(1)) * t784;
t788 = sin(qJ(3));
t792 = cos(qJ(3));
t705 = -t728 * t788 + t792 * t729;
t806 = t789 * t809;
t759 = t779 * t792 - t788 * t806;
t741 = qJD(3) * t759 + t770 * t792 + t778 * t788;
t760 = t779 * t788 + t792 * t806;
t763 = qJDD(3) + t771;
t805 = t784 * t808;
t774 = qJD(3) - t805;
t698 = (t759 * t774 - t741) * qJ(4) + (t759 * t760 + t763) * pkin(3) + t705;
t706 = t792 * t728 + t788 * t729;
t740 = -qJD(3) * t760 - t770 * t788 + t778 * t792;
t750 = pkin(3) * t774 - qJ(4) * t760;
t758 = t759 ^ 2;
t700 = -pkin(3) * t758 + qJ(4) * t740 - t750 * t774 + t706;
t783 = sin(pkin(10));
t785 = cos(pkin(10));
t746 = t759 * t785 - t760 * t783;
t819 = 2 * qJD(4);
t695 = t783 * t698 + t785 * t700 + t746 * t819;
t747 = t759 * t783 + t760 * t785;
t723 = -pkin(4) * t746 - pkin(9) * t747;
t773 = t774 ^ 2;
t693 = -pkin(4) * t773 + pkin(9) * t763 + t723 * t746 + t695;
t812 = t786 * t793;
t814 = t784 * t793;
t742 = -g(3) * t814 + t766 * t812 - t789 * t767;
t727 = -pkin(2) * t778 - pkin(8) * t777 + t769 * t806 - t742;
t701 = -pkin(3) * t740 - qJ(4) * t758 + t760 * t750 + qJDD(4) + t727;
t716 = t740 * t785 - t741 * t783;
t717 = t740 * t783 + t741 * t785;
t696 = (-t746 * t774 - t717) * pkin(9) + (t747 * t774 - t716) * pkin(4) + t701;
t787 = sin(qJ(5));
t791 = cos(qJ(5));
t690 = -t693 * t787 + t696 * t791;
t730 = -t747 * t787 + t774 * t791;
t704 = qJD(5) * t730 + t717 * t791 + t763 * t787;
t731 = t747 * t791 + t774 * t787;
t711 = -mrSges(6,1) * t730 + mrSges(6,2) * t731;
t745 = qJD(5) - t746;
t712 = -mrSges(6,2) * t745 + mrSges(6,3) * t730;
t715 = qJDD(5) - t716;
t687 = m(6) * t690 + mrSges(6,1) * t715 - t704 * mrSges(6,3) - t711 * t731 + t712 * t745;
t691 = t693 * t791 + t696 * t787;
t703 = -qJD(5) * t731 - t717 * t787 + t763 * t791;
t713 = mrSges(6,1) * t745 - mrSges(6,3) * t731;
t688 = m(6) * t691 - mrSges(6,2) * t715 + t703 * mrSges(6,3) + t711 * t730 - t713 * t745;
t679 = -t687 * t787 + t791 * t688;
t722 = -mrSges(5,1) * t746 + mrSges(5,2) * t747;
t733 = mrSges(5,1) * t774 - mrSges(5,3) * t747;
t676 = m(5) * t695 - mrSges(5,2) * t763 + mrSges(5,3) * t716 + t722 * t746 - t733 * t774 + t679;
t801 = -t698 * t785 + t700 * t783;
t692 = -pkin(4) * t763 - pkin(9) * t773 + (t819 + t723) * t747 + t801;
t689 = -m(6) * t692 + t703 * mrSges(6,1) - t704 * mrSges(6,2) + t730 * t712 - t713 * t731;
t694 = -0.2e1 * qJD(4) * t747 - t801;
t732 = -mrSges(5,2) * t774 + mrSges(5,3) * t746;
t683 = m(5) * t694 + mrSges(5,1) * t763 - mrSges(5,3) * t717 - t722 * t747 + t732 * t774 + t689;
t670 = t783 * t676 + t785 * t683;
t707 = Ifges(6,5) * t731 + Ifges(6,6) * t730 + Ifges(6,3) * t745;
t709 = Ifges(6,1) * t731 + Ifges(6,4) * t730 + Ifges(6,5) * t745;
t680 = -mrSges(6,1) * t692 + mrSges(6,3) * t691 + Ifges(6,4) * t704 + Ifges(6,2) * t703 + Ifges(6,6) * t715 - t707 * t731 + t709 * t745;
t708 = Ifges(6,4) * t731 + Ifges(6,2) * t730 + Ifges(6,6) * t745;
t681 = mrSges(6,2) * t692 - mrSges(6,3) * t690 + Ifges(6,1) * t704 + Ifges(6,4) * t703 + Ifges(6,5) * t715 + t707 * t730 - t708 * t745;
t719 = Ifges(5,4) * t747 + Ifges(5,2) * t746 + Ifges(5,6) * t774;
t720 = Ifges(5,1) * t747 + Ifges(5,4) * t746 + Ifges(5,5) * t774;
t735 = Ifges(4,4) * t760 + Ifges(4,2) * t759 + Ifges(4,6) * t774;
t736 = Ifges(4,1) * t760 + Ifges(4,4) * t759 + Ifges(4,5) * t774;
t820 = Ifges(4,5) * t741 + Ifges(4,6) * t740 + t760 * t735 - t759 * t736 + mrSges(4,1) * t705 - mrSges(4,2) * t706 + Ifges(5,5) * t717 + Ifges(5,6) * t716 + t747 * t719 - t746 * t720 + mrSges(5,1) * t694 - mrSges(5,2) * t695 + t787 * t681 + t791 * t680 + pkin(4) * t689 + pkin(9) * t679 + pkin(3) * t670 + (Ifges(4,3) + Ifges(5,3)) * t763;
t815 = t784 * t789;
t743 = -g(3) * t815 + t810;
t764 = mrSges(3,1) * t779 - mrSges(3,3) * t806;
t768 = (-mrSges(3,1) * t793 + mrSges(3,2) * t789) * t809;
t748 = -mrSges(4,1) * t759 + mrSges(4,2) * t760;
t749 = -mrSges(4,2) * t774 + mrSges(4,3) * t759;
t668 = m(4) * t705 + mrSges(4,1) * t763 - mrSges(4,3) * t741 - t748 * t760 + t749 * t774 + t670;
t751 = mrSges(4,1) * t774 - mrSges(4,3) * t760;
t802 = t785 * t676 - t683 * t783;
t669 = m(4) * t706 - mrSges(4,2) * t763 + mrSges(4,3) * t740 + t748 * t759 - t751 * t774 + t802;
t803 = -t668 * t788 + t792 * t669;
t659 = m(3) * t743 - mrSges(3,2) * t778 - mrSges(3,3) * t771 - t764 * t779 + t768 * t805 + t803;
t662 = t792 * t668 + t788 * t669;
t755 = -t766 * t784 - t817;
t765 = -mrSges(3,2) * t779 + mrSges(3,3) * t805;
t661 = m(3) * t755 + mrSges(3,1) * t771 + mrSges(3,2) * t770 + (t764 * t789 - t765 * t793) * t809 + t662;
t678 = t791 * t687 + t787 * t688;
t677 = m(5) * t701 - t716 * mrSges(5,1) + mrSges(5,2) * t717 - t746 * t732 + t733 * t747 + t678;
t797 = -m(4) * t727 + t740 * mrSges(4,1) - mrSges(4,2) * t741 + t759 * t749 - t751 * t760 - t677;
t673 = m(3) * t742 + mrSges(3,1) * t778 - mrSges(3,3) * t770 + t765 * t779 - t768 * t806 + t797;
t648 = t659 * t813 - t661 * t784 + t673 * t812;
t645 = m(2) * t775 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t795 + t648;
t655 = t793 * t659 - t673 * t789;
t653 = m(2) * t776 - mrSges(2,1) * t795 - qJDD(1) * mrSges(2,2) + t655;
t811 = t794 * t645 + t790 * t653;
t647 = t659 * t815 + t786 * t661 + t673 * t814;
t804 = -t645 * t790 + t794 * t653;
t718 = Ifges(5,5) * t747 + Ifges(5,6) * t746 + Ifges(5,3) * t774;
t663 = mrSges(5,2) * t701 - mrSges(5,3) * t694 + Ifges(5,1) * t717 + Ifges(5,4) * t716 + Ifges(5,5) * t763 - pkin(9) * t678 - t680 * t787 + t681 * t791 + t718 * t746 - t719 * t774;
t798 = mrSges(6,1) * t690 - mrSges(6,2) * t691 + Ifges(6,5) * t704 + Ifges(6,6) * t703 + Ifges(6,3) * t715 + t708 * t731 - t709 * t730;
t664 = -mrSges(5,1) * t701 + mrSges(5,3) * t695 + Ifges(5,4) * t717 + Ifges(5,2) * t716 + Ifges(5,6) * t763 - pkin(4) * t678 - t718 * t747 + t720 * t774 - t798;
t734 = Ifges(4,5) * t760 + Ifges(4,6) * t759 + Ifges(4,3) * t774;
t649 = -mrSges(4,1) * t727 + mrSges(4,3) * t706 + Ifges(4,4) * t741 + Ifges(4,2) * t740 + Ifges(4,6) * t763 - pkin(3) * t677 + qJ(4) * t802 + t783 * t663 + t785 * t664 - t760 * t734 + t774 * t736;
t650 = mrSges(4,2) * t727 - mrSges(4,3) * t705 + Ifges(4,1) * t741 + Ifges(4,4) * t740 + Ifges(4,5) * t763 - qJ(4) * t670 + t663 * t785 - t664 * t783 + t734 * t759 - t735 * t774;
t753 = Ifges(3,6) * t779 + (Ifges(3,4) * t789 + Ifges(3,2) * t793) * t809;
t754 = Ifges(3,5) * t779 + (Ifges(3,1) * t789 + Ifges(3,4) * t793) * t809;
t639 = Ifges(3,5) * t770 - Ifges(3,6) * t771 + Ifges(3,3) * t778 + mrSges(3,1) * t742 - mrSges(3,2) * t743 + t788 * t650 + t792 * t649 + pkin(2) * t797 + pkin(8) * t803 + (t753 * t789 - t754 * t793) * t809;
t752 = Ifges(3,3) * t779 + (Ifges(3,5) * t789 + Ifges(3,6) * t793) * t809;
t641 = mrSges(3,2) * t755 - mrSges(3,3) * t742 + Ifges(3,1) * t770 - Ifges(3,4) * t771 + Ifges(3,5) * t778 - pkin(8) * t662 - t649 * t788 + t650 * t792 + t752 * t805 - t753 * t779;
t643 = -mrSges(3,1) * t755 + mrSges(3,3) * t743 + Ifges(3,4) * t770 - Ifges(3,2) * t771 + Ifges(3,6) * t778 - pkin(2) * t662 - t752 * t806 + t779 * t754 - t820;
t799 = mrSges(2,1) * t775 - mrSges(2,2) * t776 + Ifges(2,3) * qJDD(1) + pkin(1) * t648 + t786 * t639 + t641 * t815 + t643 * t814 + t655 * t818;
t637 = -mrSges(2,2) * g(3) - mrSges(2,3) * t775 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t795 + t641 * t793 - t643 * t789 + (-t647 * t784 - t648 * t786) * pkin(7);
t636 = mrSges(2,1) * g(3) + mrSges(2,3) * t776 + Ifges(2,5) * t795 + Ifges(2,6) * qJDD(1) - pkin(1) * t647 - t639 * t784 + (pkin(7) * t655 + t641 * t789 + t643 * t793) * t786;
t1 = [-m(1) * g(1) + t804; -m(1) * g(2) + t811; (-m(1) - m(2)) * g(3) + t647; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t811 - t790 * t636 + t794 * t637; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t804 + t794 * t636 + t790 * t637; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t799; t799; t639; t820; t677; t798;];
tauJB = t1;
