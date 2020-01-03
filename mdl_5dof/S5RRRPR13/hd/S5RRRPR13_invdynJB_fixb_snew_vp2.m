% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:44:00
% EndTime: 2019-12-31 21:44:09
% DurationCPUTime: 6.40s
% Computational Cost: add. (90957->309), mult. (196631->385), div. (0->0), fcn. (146239->10), ass. (0->134)
t848 = Ifges(4,1) + Ifges(5,2);
t840 = Ifges(4,4) + Ifges(5,6);
t839 = Ifges(4,5) - Ifges(5,4);
t847 = -Ifges(4,2) - Ifges(5,3);
t838 = Ifges(4,6) - Ifges(5,5);
t846 = Ifges(4,3) + Ifges(5,1);
t798 = sin(pkin(5));
t802 = sin(qJ(2));
t805 = cos(qJ(2));
t823 = qJD(1) * qJD(2);
t783 = (-qJDD(1) * t805 + t802 * t823) * t798;
t799 = cos(pkin(5));
t794 = t799 * qJD(1) + qJD(2);
t801 = sin(qJ(3));
t825 = qJD(1) * t798;
t822 = t802 * t825;
t843 = cos(qJ(3));
t770 = -t843 * t794 + t801 * t822;
t824 = qJD(1) * t805;
t821 = t798 * t824;
t788 = -qJD(3) + t821;
t757 = t770 * mrSges(5,1) + t788 * mrSges(5,3);
t775 = qJDD(3) + t783;
t781 = (-pkin(2) * t805 - pkin(8) * t802) * t825;
t792 = t794 ^ 2;
t793 = t799 * qJDD(1) + qJDD(2);
t803 = sin(qJ(1));
t806 = cos(qJ(1));
t790 = t803 * g(1) - t806 * g(2);
t807 = qJD(1) ^ 2;
t842 = pkin(7) * t798;
t778 = qJDD(1) * pkin(1) + t807 * t842 + t790;
t791 = -t806 * g(1) - t803 * g(2);
t779 = -t807 * pkin(1) + qJDD(1) * t842 + t791;
t833 = t799 * t802;
t826 = t778 * t833 + t805 * t779;
t725 = -t792 * pkin(2) + t793 * pkin(8) + (-g(3) * t802 + t781 * t824) * t798 + t826;
t782 = (qJDD(1) * t802 + t805 * t823) * t798;
t841 = t799 * g(3);
t726 = t783 * pkin(2) - t782 * pkin(8) - t841 + (-t778 + (pkin(2) * t802 - pkin(8) * t805) * t794 * qJD(1)) * t798;
t709 = -t801 * t725 + t843 * t726;
t771 = t801 * t794 + t843 * t822;
t749 = t770 * pkin(3) - t771 * qJ(4);
t787 = t788 ^ 2;
t707 = -t775 * pkin(3) - t787 * qJ(4) + t771 * t749 + qJDD(4) - t709;
t746 = -t770 * qJD(3) + t843 * t782 + t801 * t793;
t836 = t770 * t788;
t702 = (t770 * t771 - t775) * pkin(9) + (t746 - t836) * pkin(4) + t707;
t745 = t771 * qJD(3) + t801 * t782 - t843 * t793;
t759 = t771 * pkin(4) + t788 * pkin(9);
t769 = t770 ^ 2;
t832 = t799 * t805;
t834 = t798 * t805;
t747 = -g(3) * t834 + t778 * t832 - t802 * t779;
t724 = -t793 * pkin(2) - t792 * pkin(8) + t781 * t822 - t747;
t844 = -2 * qJD(4);
t809 = (-t746 - t836) * qJ(4) + t724 + (-t788 * pkin(3) + t844) * t771;
t705 = -t769 * pkin(4) - t771 * t759 + (pkin(3) + pkin(9)) * t745 + t809;
t800 = sin(qJ(5));
t804 = cos(qJ(5));
t700 = t804 * t702 - t800 * t705;
t753 = t804 * t770 + t800 * t788;
t716 = t753 * qJD(5) + t800 * t745 + t804 * t775;
t754 = t800 * t770 - t804 * t788;
t727 = -t753 * mrSges(6,1) + t754 * mrSges(6,2);
t768 = qJD(5) + t771;
t730 = -t768 * mrSges(6,2) + t753 * mrSges(6,3);
t742 = qJDD(5) + t746;
t697 = m(6) * t700 + t742 * mrSges(6,1) - t716 * mrSges(6,3) - t754 * t727 + t768 * t730;
t701 = t800 * t702 + t804 * t705;
t715 = -t754 * qJD(5) + t804 * t745 - t800 * t775;
t731 = t768 * mrSges(6,1) - t754 * mrSges(6,3);
t698 = m(6) * t701 - t742 * mrSges(6,2) + t715 * mrSges(6,3) + t753 * t727 - t768 * t731;
t688 = t804 * t697 + t800 * t698;
t751 = -t770 * mrSges(5,2) - t771 * mrSges(5,3);
t815 = -m(5) * t707 - t746 * mrSges(5,1) - t771 * t751 - t688;
t687 = t775 * mrSges(5,2) - t788 * t757 - t815;
t710 = t843 * t725 + t801 * t726;
t814 = -t787 * pkin(3) + t775 * qJ(4) - t770 * t749 + t710;
t704 = -t745 * pkin(4) - t769 * pkin(9) + (t844 - t759) * t788 + t814;
t717 = Ifges(6,5) * t754 + Ifges(6,6) * t753 + Ifges(6,3) * t768;
t719 = Ifges(6,1) * t754 + Ifges(6,4) * t753 + Ifges(6,5) * t768;
t689 = -mrSges(6,1) * t704 + mrSges(6,3) * t701 + Ifges(6,4) * t716 + Ifges(6,2) * t715 + Ifges(6,6) * t742 - t754 * t717 + t768 * t719;
t718 = Ifges(6,4) * t754 + Ifges(6,2) * t753 + Ifges(6,6) * t768;
t690 = mrSges(6,2) * t704 - mrSges(6,3) * t700 + Ifges(6,1) * t716 + Ifges(6,4) * t715 + Ifges(6,5) * t742 + t753 * t717 - t768 * t718;
t706 = 0.2e1 * qJD(4) * t788 - t814;
t758 = t771 * mrSges(5,1) - t788 * mrSges(5,2);
t816 = -m(6) * t704 + t715 * mrSges(6,1) - t716 * mrSges(6,2) + t753 * t730 - t754 * t731;
t811 = -m(5) * t706 + t775 * mrSges(5,3) - t788 * t758 - t816;
t827 = -t840 * t770 + t848 * t771 - t839 * t788;
t828 = t847 * t770 + t840 * t771 - t838 * t788;
t845 = -t838 * t745 + t839 * t746 + t827 * t770 + t828 * t771 + t846 * t775 + mrSges(4,1) * t709 - mrSges(4,2) * t710 + mrSges(5,2) * t707 - mrSges(5,3) * t706 - pkin(3) * t687 - pkin(9) * t688 + qJ(4) * (-t745 * mrSges(5,1) - t770 * t751 + t811) - t800 * t689 + t804 * t690;
t835 = t798 * t802;
t748 = -g(3) * t835 + t826;
t776 = t794 * mrSges(3,1) - mrSges(3,3) * t822;
t780 = (-mrSges(3,1) * t805 + mrSges(3,2) * t802) * t825;
t750 = t770 * mrSges(4,1) + t771 * mrSges(4,2);
t755 = t788 * mrSges(4,2) - t770 * mrSges(4,3);
t685 = m(4) * t709 - t746 * mrSges(4,3) - t771 * t750 + (-t755 + t757) * t788 + (mrSges(4,1) - mrSges(5,2)) * t775 + t815;
t756 = -t788 * mrSges(4,1) - t771 * mrSges(4,3);
t693 = m(4) * t710 - t775 * mrSges(4,2) + t788 * t756 + (-t750 - t751) * t770 + (-mrSges(4,3) - mrSges(5,1)) * t745 + t811;
t819 = -t801 * t685 + t843 * t693;
t677 = m(3) * t748 - t793 * mrSges(3,2) - t783 * mrSges(3,3) - t794 * t776 + t780 * t821 + t819;
t680 = t843 * t685 + t801 * t693;
t764 = -t798 * t778 - t841;
t777 = -t794 * mrSges(3,2) + mrSges(3,3) * t821;
t679 = m(3) * t764 + t783 * mrSges(3,1) + t782 * mrSges(3,2) + (t776 * t802 - t777 * t805) * t825 + t680;
t708 = t745 * pkin(3) + t809;
t830 = -t800 * t697 + t804 * t698;
t818 = -m(5) * t708 + t745 * mrSges(5,2) + t770 * t757 - t830;
t810 = -m(4) * t724 - t745 * mrSges(4,1) - t770 * t755 + (-t756 + t758) * t771 + (-mrSges(4,2) + mrSges(5,3)) * t746 + t818;
t683 = m(3) * t747 + t793 * mrSges(3,1) - t782 * mrSges(3,3) + t794 * t777 - t780 * t822 + t810;
t666 = t677 * t833 - t798 * t679 + t683 * t832;
t663 = m(2) * t790 + qJDD(1) * mrSges(2,1) - t807 * mrSges(2,2) + t666;
t673 = t805 * t677 - t802 * t683;
t671 = m(2) * t791 - t807 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t673;
t831 = t806 * t663 + t803 * t671;
t829 = t838 * t770 - t839 * t771 + t846 * t788;
t665 = t677 * t835 + t799 * t679 + t683 * t834;
t820 = -t803 * t663 + t806 * t671;
t686 = -t746 * mrSges(5,3) - t771 * t758 - t818;
t667 = -mrSges(4,1) * t724 - mrSges(5,1) * t706 + mrSges(5,2) * t708 + mrSges(4,3) * t710 - pkin(3) * t686 - pkin(4) * t816 - pkin(9) * t830 - t804 * t689 - t800 * t690 + t847 * t745 + t840 * t746 + t829 * t771 + t838 * t775 - t827 * t788;
t812 = mrSges(6,1) * t700 - mrSges(6,2) * t701 + Ifges(6,5) * t716 + Ifges(6,6) * t715 + Ifges(6,3) * t742 + t754 * t718 - t753 * t719;
t668 = mrSges(5,1) * t707 + mrSges(4,2) * t724 - mrSges(4,3) * t709 - mrSges(5,3) * t708 + pkin(4) * t688 - qJ(4) * t686 - t840 * t745 + t848 * t746 + t829 * t770 + t839 * t775 + t828 * t788 + t812;
t762 = Ifges(3,6) * t794 + (Ifges(3,4) * t802 + Ifges(3,2) * t805) * t825;
t763 = Ifges(3,5) * t794 + (Ifges(3,1) * t802 + Ifges(3,4) * t805) * t825;
t657 = Ifges(3,5) * t782 - Ifges(3,6) * t783 + Ifges(3,3) * t793 + mrSges(3,1) * t747 - mrSges(3,2) * t748 + t801 * t668 + t843 * t667 + pkin(2) * t810 + pkin(8) * t819 + (t762 * t802 - t763 * t805) * t825;
t761 = Ifges(3,3) * t794 + (Ifges(3,5) * t802 + Ifges(3,6) * t805) * t825;
t659 = mrSges(3,2) * t764 - mrSges(3,3) * t747 + Ifges(3,1) * t782 - Ifges(3,4) * t783 + Ifges(3,5) * t793 - pkin(8) * t680 - t801 * t667 + t843 * t668 + t761 * t821 - t794 * t762;
t661 = -mrSges(3,1) * t764 + mrSges(3,3) * t748 + Ifges(3,4) * t782 - Ifges(3,2) * t783 + Ifges(3,6) * t793 - pkin(2) * t680 - t761 * t822 + t794 * t763 - t845;
t813 = mrSges(2,1) * t790 - mrSges(2,2) * t791 + Ifges(2,3) * qJDD(1) + pkin(1) * t666 + t799 * t657 + t659 * t835 + t661 * t834 + t673 * t842;
t655 = -mrSges(2,2) * g(3) - mrSges(2,3) * t790 + Ifges(2,5) * qJDD(1) - t807 * Ifges(2,6) + t805 * t659 - t802 * t661 + (-t665 * t798 - t666 * t799) * pkin(7);
t654 = mrSges(2,1) * g(3) + mrSges(2,3) * t791 + t807 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t665 - t798 * t657 + (pkin(7) * t673 + t659 * t802 + t661 * t805) * t799;
t1 = [-m(1) * g(1) + t820; -m(1) * g(2) + t831; (-m(1) - m(2)) * g(3) + t665; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t831 - t803 * t654 + t806 * t655; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t820 + t806 * t654 + t803 * t655; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t813; t813; t657; t845; t687; t812;];
tauJB = t1;
