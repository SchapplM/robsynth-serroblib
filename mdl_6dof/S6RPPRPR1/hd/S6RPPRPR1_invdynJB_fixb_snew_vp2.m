% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:55:10
% EndTime: 2019-05-05 13:55:22
% DurationCPUTime: 12.07s
% Computational Cost: add. (193290->322), mult. (446856->405), div. (0->0), fcn. (316996->12), ass. (0->141)
t854 = qJD(1) ^ 2;
t885 = cos(qJ(4));
t846 = cos(pkin(10));
t884 = pkin(3) * t846;
t843 = sin(pkin(10));
t883 = mrSges(4,2) * t843;
t838 = t846 ^ 2;
t882 = t838 * t854;
t850 = sin(qJ(1));
t852 = cos(qJ(1));
t824 = t850 * g(1) - g(2) * t852;
t821 = qJDD(1) * pkin(1) + t824;
t825 = -g(1) * t852 - g(2) * t850;
t822 = -pkin(1) * t854 + t825;
t844 = sin(pkin(9));
t847 = cos(pkin(9));
t805 = t844 * t821 + t847 * t822;
t795 = -pkin(2) * t854 + qJDD(1) * qJ(3) + t805;
t841 = -g(3) + qJDD(2);
t876 = qJD(1) * qJD(3);
t880 = t846 * t841 - 0.2e1 * t843 * t876;
t776 = (-pkin(7) * qJDD(1) + t854 * t884 - t795) * t843 + t880;
t782 = t843 * t841 + (t795 + 0.2e1 * t876) * t846;
t874 = qJDD(1) * t846;
t777 = -pkin(3) * t882 + pkin(7) * t874 + t782;
t849 = sin(qJ(4));
t761 = t849 * t776 + t885 * t777;
t873 = t846 * t885;
t879 = qJD(1) * t843;
t814 = -qJD(1) * t873 + t849 * t879;
t861 = t885 * t843 + t846 * t849;
t815 = t861 * qJD(1);
t797 = pkin(4) * t814 - qJ(5) * t815;
t853 = qJD(4) ^ 2;
t753 = -pkin(4) * t853 + qJDD(4) * qJ(5) - t797 * t814 + t761;
t837 = t843 ^ 2;
t804 = t847 * t821 - t844 * t822;
t863 = qJDD(3) - t804;
t778 = (-pkin(2) - t884) * qJDD(1) + (-qJ(3) + (-t837 - t838) * pkin(7)) * t854 + t863;
t875 = qJDD(1) * t843;
t878 = qJD(4) * t815;
t801 = -qJDD(1) * t873 + t849 * t875 + t878;
t877 = t814 * qJD(4);
t802 = t861 * qJDD(1) - t877;
t756 = (-t802 + t877) * qJ(5) + (t801 + t878) * pkin(4) + t778;
t842 = sin(pkin(11));
t845 = cos(pkin(11));
t810 = qJD(4) * t842 + t815 * t845;
t748 = -0.2e1 * qJD(5) * t810 - t842 * t753 + t845 * t756;
t789 = qJDD(4) * t842 + t802 * t845;
t809 = qJD(4) * t845 - t815 * t842;
t746 = (t809 * t814 - t789) * pkin(8) + (t809 * t810 + t801) * pkin(5) + t748;
t749 = 0.2e1 * qJD(5) * t809 + t845 * t753 + t842 * t756;
t787 = pkin(5) * t814 - pkin(8) * t810;
t788 = qJDD(4) * t845 - t802 * t842;
t808 = t809 ^ 2;
t747 = -pkin(5) * t808 + pkin(8) * t788 - t787 * t814 + t749;
t848 = sin(qJ(6));
t851 = cos(qJ(6));
t744 = t746 * t851 - t747 * t848;
t779 = t809 * t851 - t810 * t848;
t759 = qJD(6) * t779 + t788 * t848 + t789 * t851;
t780 = t809 * t848 + t810 * t851;
t766 = -mrSges(7,1) * t779 + mrSges(7,2) * t780;
t813 = qJD(6) + t814;
t767 = -mrSges(7,2) * t813 + mrSges(7,3) * t779;
t800 = qJDD(6) + t801;
t741 = m(7) * t744 + mrSges(7,1) * t800 - mrSges(7,3) * t759 - t766 * t780 + t767 * t813;
t745 = t746 * t848 + t747 * t851;
t758 = -qJD(6) * t780 + t788 * t851 - t789 * t848;
t768 = mrSges(7,1) * t813 - mrSges(7,3) * t780;
t742 = m(7) * t745 - mrSges(7,2) * t800 + mrSges(7,3) * t758 + t766 * t779 - t768 * t813;
t733 = t851 * t741 + t848 * t742;
t783 = -mrSges(6,1) * t809 + mrSges(6,2) * t810;
t867 = -mrSges(6,2) * t814 + mrSges(6,3) * t809;
t731 = m(6) * t748 + t801 * mrSges(6,1) - t789 * mrSges(6,3) - t810 * t783 + t814 * t867 + t733;
t786 = mrSges(6,1) * t814 - mrSges(6,3) * t810;
t868 = -t741 * t848 + t851 * t742;
t732 = m(6) * t749 - mrSges(6,2) * t801 + mrSges(6,3) * t788 + t783 * t809 - t786 * t814 + t868;
t727 = -t731 * t842 + t845 * t732;
t798 = mrSges(5,1) * t814 + mrSges(5,2) * t815;
t812 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t815;
t725 = m(5) * t761 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t801 - qJD(4) * t812 - t798 * t814 + t727;
t760 = t885 * t776 - t849 * t777;
t752 = -qJDD(4) * pkin(4) - t853 * qJ(5) + t815 * t797 + qJDD(5) - t760;
t750 = -t788 * pkin(5) - t808 * pkin(8) + t810 * t787 + t752;
t860 = m(7) * t750 - t758 * mrSges(7,1) + mrSges(7,2) * t759 - t779 * t767 + t768 * t780;
t743 = m(6) * t752 - t788 * mrSges(6,1) + mrSges(6,2) * t789 + t786 * t810 - t809 * t867 + t860;
t811 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t814;
t737 = m(5) * t760 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t802 + qJD(4) * t811 - t798 * t815 - t743;
t716 = t849 * t725 + t885 * t737;
t781 = -t795 * t843 + t880;
t862 = mrSges(4,3) * qJDD(1) + t854 * (-mrSges(4,1) * t846 + t883);
t714 = m(4) * t781 - t862 * t843 + t716;
t869 = t885 * t725 - t849 * t737;
t715 = m(4) * t782 + t862 * t846 + t869;
t870 = -t714 * t843 + t846 * t715;
t706 = m(3) * t805 - mrSges(3,1) * t854 - qJDD(1) * mrSges(3,2) + t870;
t791 = -qJDD(1) * pkin(2) - t854 * qJ(3) + t863;
t726 = t845 * t731 + t842 * t732;
t859 = m(5) * t778 + t801 * mrSges(5,1) + mrSges(5,2) * t802 + t814 * t811 + t812 * t815 + t726;
t857 = -m(4) * t791 + mrSges(4,1) * t874 - t859 + (t837 * t854 + t882) * mrSges(4,3);
t720 = t857 + (mrSges(3,1) - t883) * qJDD(1) + m(3) * t804 - mrSges(3,2) * t854;
t702 = t844 * t706 + t847 * t720;
t699 = m(2) * t824 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t854 + t702;
t871 = t847 * t706 - t720 * t844;
t700 = m(2) * t825 - mrSges(2,1) * t854 - qJDD(1) * mrSges(2,2) + t871;
t881 = t852 * t699 + t850 * t700;
t709 = t846 * t714 + t843 * t715;
t707 = m(3) * t841 + t709;
t872 = -t699 * t850 + t852 * t700;
t866 = Ifges(4,1) * t843 + Ifges(4,4) * t846;
t865 = Ifges(4,4) * t843 + Ifges(4,2) * t846;
t864 = Ifges(4,5) * t843 + Ifges(4,6) * t846;
t762 = Ifges(7,5) * t780 + Ifges(7,6) * t779 + Ifges(7,3) * t813;
t764 = Ifges(7,1) * t780 + Ifges(7,4) * t779 + Ifges(7,5) * t813;
t734 = -mrSges(7,1) * t750 + mrSges(7,3) * t745 + Ifges(7,4) * t759 + Ifges(7,2) * t758 + Ifges(7,6) * t800 - t762 * t780 + t764 * t813;
t763 = Ifges(7,4) * t780 + Ifges(7,2) * t779 + Ifges(7,6) * t813;
t735 = mrSges(7,2) * t750 - mrSges(7,3) * t744 + Ifges(7,1) * t759 + Ifges(7,4) * t758 + Ifges(7,5) * t800 + t762 * t779 - t763 * t813;
t770 = Ifges(6,5) * t810 + Ifges(6,6) * t809 + Ifges(6,3) * t814;
t772 = Ifges(6,1) * t810 + Ifges(6,4) * t809 + Ifges(6,5) * t814;
t717 = -mrSges(6,1) * t752 + mrSges(6,3) * t749 + Ifges(6,4) * t789 + Ifges(6,2) * t788 + Ifges(6,6) * t801 - pkin(5) * t860 + pkin(8) * t868 + t851 * t734 + t848 * t735 - t810 * t770 + t814 * t772;
t771 = Ifges(6,4) * t810 + Ifges(6,2) * t809 + Ifges(6,6) * t814;
t718 = mrSges(6,2) * t752 - mrSges(6,3) * t748 + Ifges(6,1) * t789 + Ifges(6,4) * t788 + Ifges(6,5) * t801 - pkin(8) * t733 - t734 * t848 + t735 * t851 + t770 * t809 - t771 * t814;
t792 = Ifges(5,5) * t815 - Ifges(5,6) * t814 + Ifges(5,3) * qJD(4);
t793 = Ifges(5,4) * t815 - Ifges(5,2) * t814 + Ifges(5,6) * qJD(4);
t703 = mrSges(5,2) * t778 - mrSges(5,3) * t760 + Ifges(5,1) * t802 - Ifges(5,4) * t801 + Ifges(5,5) * qJDD(4) - qJ(5) * t726 - qJD(4) * t793 - t717 * t842 + t718 * t845 - t792 * t814;
t794 = Ifges(5,1) * t815 - Ifges(5,4) * t814 + Ifges(5,5) * qJD(4);
t856 = mrSges(7,1) * t744 - mrSges(7,2) * t745 + Ifges(7,5) * t759 + Ifges(7,6) * t758 + Ifges(7,3) * t800 + t780 * t763 - t779 * t764;
t710 = -pkin(5) * t733 - pkin(4) * t726 + Ifges(5,4) * t802 + Ifges(5,6) * qJDD(4) - t815 * t792 - Ifges(6,6) * t788 - Ifges(6,5) * t789 + qJD(4) * t794 - mrSges(5,1) * t778 + mrSges(5,3) * t761 - t856 + t809 * t772 - t810 * t771 - mrSges(6,1) * t748 + mrSges(6,2) * t749 + (-Ifges(5,2) - Ifges(6,3)) * t801;
t820 = t864 * qJD(1);
t692 = -mrSges(4,1) * t791 + mrSges(4,3) * t782 - pkin(3) * t859 + pkin(7) * t869 + t865 * qJDD(1) + t849 * t703 + t885 * t710 - t820 * t879;
t695 = t846 * qJD(1) * t820 + mrSges(4,2) * t791 - mrSges(4,3) * t781 - pkin(7) * t716 + t866 * qJDD(1) + t885 * t703 - t849 * t710;
t722 = mrSges(4,2) * t875 - t857;
t858 = mrSges(2,1) * t824 + mrSges(3,1) * t804 - mrSges(2,2) * t825 - mrSges(3,2) * t805 + pkin(1) * t702 - pkin(2) * t722 + qJ(3) * t870 + t846 * t692 + t843 * t695 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t855 = mrSges(5,1) * t760 - mrSges(5,2) * t761 + Ifges(5,5) * t802 - Ifges(5,6) * t801 + Ifges(5,3) * qJDD(4) - pkin(4) * t743 + qJ(5) * t727 + t845 * t717 + t842 * t718 + t815 * t793 + t814 * t794;
t693 = mrSges(3,3) * t805 - pkin(3) * t716 - t855 - pkin(2) * t709 - mrSges(3,1) * t841 - mrSges(4,1) * t781 + mrSges(4,2) * t782 + (Ifges(3,6) - t864) * qJDD(1) + (-t843 * t865 + t846 * t866 + Ifges(3,5)) * t854;
t690 = mrSges(3,2) * t841 - mrSges(3,3) * t804 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t854 - qJ(3) * t709 - t692 * t843 + t695 * t846;
t689 = -mrSges(2,2) * g(3) - mrSges(2,3) * t824 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t854 - qJ(2) * t702 + t690 * t847 - t693 * t844;
t688 = mrSges(2,1) * g(3) + mrSges(2,3) * t825 + t854 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t707 + qJ(2) * t871 + t844 * t690 + t847 * t693;
t1 = [-m(1) * g(1) + t872; -m(1) * g(2) + t881; (-m(1) - m(2)) * g(3) + t707; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t881 - t850 * t688 + t852 * t689; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t872 + t852 * t688 + t850 * t689; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t858; t858; t707; t722; t855; t743; t856;];
tauJB  = t1;
