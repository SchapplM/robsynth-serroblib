% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:02:53
% EndTime: 2019-05-05 11:03:22
% DurationCPUTime: 27.78s
% Computational Cost: add. (493825->343), mult. (973845->437), div. (0->0), fcn. (717683->14), ass. (0->148)
t865 = sin(pkin(12));
t867 = cos(pkin(12));
t854 = g(1) * t865 - g(2) * t867;
t855 = -g(1) * t867 - g(2) * t865;
t864 = -g(3) + qJDD(1);
t866 = sin(pkin(6));
t868 = cos(pkin(6));
t873 = sin(qJ(2));
t878 = cos(qJ(2));
t814 = -t873 * t855 + (t854 * t868 + t864 * t866) * t878;
t901 = t868 * t873;
t902 = t866 * t873;
t815 = t854 * t901 + t878 * t855 + t864 * t902;
t880 = qJD(2) ^ 2;
t810 = -pkin(2) * t880 + qJDD(2) * pkin(8) + t815;
t832 = -t854 * t866 + t864 * t868;
t872 = sin(qJ(3));
t877 = cos(qJ(3));
t803 = t877 * t810 + t872 * t832;
t851 = (-pkin(3) * t877 - pkin(9) * t872) * qJD(2);
t879 = qJD(3) ^ 2;
t898 = qJD(2) * t877;
t792 = -pkin(3) * t879 + qJDD(3) * pkin(9) + t851 * t898 + t803;
t809 = -qJDD(2) * pkin(2) - t880 * pkin(8) - t814;
t897 = qJD(2) * qJD(3);
t896 = t877 * t897;
t852 = qJDD(2) * t872 + t896;
t862 = t872 * t897;
t853 = qJDD(2) * t877 - t862;
t797 = (-t852 - t896) * pkin(9) + (-t853 + t862) * pkin(3) + t809;
t871 = sin(qJ(4));
t876 = cos(qJ(4));
t775 = -t871 * t792 + t876 * t797;
t899 = qJD(2) * t872;
t848 = qJD(3) * t876 - t871 * t899;
t823 = qJD(4) * t848 + qJDD(3) * t871 + t852 * t876;
t845 = qJDD(4) - t853;
t849 = qJD(3) * t871 + t876 * t899;
t861 = qJD(4) - t898;
t771 = (t848 * t861 - t823) * pkin(10) + (t848 * t849 + t845) * pkin(4) + t775;
t776 = t876 * t792 + t871 * t797;
t822 = -qJD(4) * t849 + qJDD(3) * t876 - t852 * t871;
t831 = pkin(4) * t861 - pkin(10) * t849;
t844 = t848 ^ 2;
t773 = -pkin(4) * t844 + pkin(10) * t822 - t831 * t861 + t776;
t870 = sin(qJ(5));
t875 = cos(qJ(5));
t759 = t875 * t771 - t870 * t773;
t825 = t848 * t875 - t849 * t870;
t789 = qJD(5) * t825 + t822 * t870 + t823 * t875;
t826 = t848 * t870 + t849 * t875;
t841 = qJDD(5) + t845;
t860 = qJD(5) + t861;
t756 = (t825 * t860 - t789) * pkin(11) + (t825 * t826 + t841) * pkin(5) + t759;
t760 = t870 * t771 + t875 * t773;
t788 = -qJD(5) * t826 + t822 * t875 - t823 * t870;
t813 = pkin(5) * t860 - pkin(11) * t826;
t824 = t825 ^ 2;
t757 = -pkin(5) * t824 + pkin(11) * t788 - t813 * t860 + t760;
t869 = sin(qJ(6));
t874 = cos(qJ(6));
t755 = t756 * t869 + t757 * t874;
t802 = -t872 * t810 + t832 * t877;
t791 = -qJDD(3) * pkin(3) - pkin(9) * t879 + t851 * t899 - t802;
t777 = -pkin(4) * t822 - pkin(10) * t844 + t849 * t831 + t791;
t762 = -pkin(5) * t788 - pkin(11) * t824 + t813 * t826 + t777;
t805 = t825 * t869 + t826 * t874;
t767 = -qJD(6) * t805 + t788 * t874 - t789 * t869;
t804 = t825 * t874 - t826 * t869;
t768 = qJD(6) * t804 + t788 * t869 + t789 * t874;
t856 = qJD(6) + t860;
t778 = Ifges(7,5) * t805 + Ifges(7,6) * t804 + Ifges(7,3) * t856;
t780 = Ifges(7,1) * t805 + Ifges(7,4) * t804 + Ifges(7,5) * t856;
t836 = qJDD(6) + t841;
t743 = -mrSges(7,1) * t762 + mrSges(7,3) * t755 + Ifges(7,4) * t768 + Ifges(7,2) * t767 + Ifges(7,6) * t836 - t778 * t805 + t780 * t856;
t754 = t756 * t874 - t757 * t869;
t779 = Ifges(7,4) * t805 + Ifges(7,2) * t804 + Ifges(7,6) * t856;
t744 = mrSges(7,2) * t762 - mrSges(7,3) * t754 + Ifges(7,1) * t768 + Ifges(7,4) * t767 + Ifges(7,5) * t836 + t778 * t804 - t779 * t856;
t798 = Ifges(6,5) * t826 + Ifges(6,6) * t825 + Ifges(6,3) * t860;
t800 = Ifges(6,1) * t826 + Ifges(6,4) * t825 + Ifges(6,5) * t860;
t795 = -mrSges(7,2) * t856 + mrSges(7,3) * t804;
t796 = mrSges(7,1) * t856 - mrSges(7,3) * t805;
t890 = m(7) * t762 - t767 * mrSges(7,1) + t768 * mrSges(7,2) - t804 * t795 + t805 * t796;
t783 = -mrSges(7,1) * t804 + mrSges(7,2) * t805;
t750 = m(7) * t754 + mrSges(7,1) * t836 - mrSges(7,3) * t768 - t783 * t805 + t795 * t856;
t751 = m(7) * t755 - mrSges(7,2) * t836 + mrSges(7,3) * t767 + t783 * t804 - t796 * t856;
t892 = -t750 * t869 + t874 * t751;
t730 = -mrSges(6,1) * t777 + mrSges(6,3) * t760 + Ifges(6,4) * t789 + Ifges(6,2) * t788 + Ifges(6,6) * t841 - pkin(5) * t890 + pkin(11) * t892 + t874 * t743 + t869 * t744 - t826 * t798 + t860 * t800;
t742 = t874 * t750 + t869 * t751;
t799 = Ifges(6,4) * t826 + Ifges(6,2) * t825 + Ifges(6,6) * t860;
t731 = mrSges(6,2) * t777 - mrSges(6,3) * t759 + Ifges(6,1) * t789 + Ifges(6,4) * t788 + Ifges(6,5) * t841 - pkin(11) * t742 - t743 * t869 + t744 * t874 + t798 * t825 - t799 * t860;
t816 = Ifges(5,5) * t849 + Ifges(5,6) * t848 + Ifges(5,3) * t861;
t818 = Ifges(5,1) * t849 + Ifges(5,4) * t848 + Ifges(5,5) * t861;
t811 = -mrSges(6,2) * t860 + mrSges(6,3) * t825;
t812 = mrSges(6,1) * t860 - mrSges(6,3) * t826;
t885 = m(6) * t777 - t788 * mrSges(6,1) + t789 * mrSges(6,2) - t825 * t811 + t826 * t812 + t890;
t806 = -mrSges(6,1) * t825 + mrSges(6,2) * t826;
t739 = m(6) * t759 + mrSges(6,1) * t841 - mrSges(6,3) * t789 - t806 * t826 + t811 * t860 + t742;
t740 = m(6) * t760 - mrSges(6,2) * t841 + mrSges(6,3) * t788 + t806 * t825 - t812 * t860 + t892;
t893 = -t739 * t870 + t875 * t740;
t713 = -mrSges(5,1) * t791 + mrSges(5,3) * t776 + Ifges(5,4) * t823 + Ifges(5,2) * t822 + Ifges(5,6) * t845 - pkin(4) * t885 + pkin(10) * t893 + t875 * t730 + t870 * t731 - t849 * t816 + t861 * t818;
t735 = t875 * t739 + t870 * t740;
t817 = Ifges(5,4) * t849 + Ifges(5,2) * t848 + Ifges(5,6) * t861;
t714 = mrSges(5,2) * t791 - mrSges(5,3) * t775 + Ifges(5,1) * t823 + Ifges(5,4) * t822 + Ifges(5,5) * t845 - pkin(10) * t735 - t730 * t870 + t731 * t875 + t816 * t848 - t817 * t861;
t827 = -mrSges(5,1) * t848 + mrSges(5,2) * t849;
t829 = -mrSges(5,2) * t861 + mrSges(5,3) * t848;
t733 = m(5) * t775 + mrSges(5,1) * t845 - mrSges(5,3) * t823 - t827 * t849 + t829 * t861 + t735;
t830 = mrSges(5,1) * t861 - mrSges(5,3) * t849;
t734 = m(5) * t776 - mrSges(5,2) * t845 + mrSges(5,3) * t822 + t827 * t848 - t830 * t861 + t893;
t729 = -t733 * t871 + t876 * t734;
t752 = -m(5) * t791 + t822 * mrSges(5,1) - t823 * mrSges(5,2) + t848 * t829 - t849 * t830 - t885;
t839 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t872 + Ifges(4,2) * t877) * qJD(2);
t840 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t872 + Ifges(4,4) * t877) * qJD(2);
t904 = mrSges(4,1) * t802 - mrSges(4,2) * t803 + Ifges(4,5) * t852 + Ifges(4,6) * t853 + Ifges(4,3) * qJDD(3) + pkin(3) * t752 + pkin(9) * t729 + t876 * t713 + t871 * t714 + (t839 * t872 - t840 * t877) * qJD(2);
t728 = t733 * t876 + t734 * t871;
t857 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t899;
t858 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t898;
t884 = -m(4) * t809 + t853 * mrSges(4,1) - mrSges(4,2) * t852 - t857 * t899 + t858 * t898 - t728;
t724 = m(3) * t814 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t880 + t884;
t903 = t724 * t878;
t850 = (-mrSges(4,1) * t877 + mrSges(4,2) * t872) * qJD(2);
t727 = m(4) * t803 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t853 - qJD(3) * t857 + t850 * t898 + t729;
t746 = m(4) * t802 + qJDD(3) * mrSges(4,1) - t852 * mrSges(4,3) + qJD(3) * t858 - t850 * t899 + t752;
t894 = t877 * t727 - t746 * t872;
t718 = m(3) * t815 - mrSges(3,1) * t880 - qJDD(2) * mrSges(3,2) + t894;
t721 = t872 * t727 + t877 * t746;
t720 = m(3) * t832 + t721;
t707 = t718 * t901 - t720 * t866 + t868 * t903;
t705 = m(2) * t854 + t707;
t711 = t878 * t718 - t724 * t873;
t710 = m(2) * t855 + t711;
t900 = t867 * t705 + t865 * t710;
t706 = t718 * t902 + t868 * t720 + t866 * t903;
t895 = -t705 * t865 + t867 * t710;
t891 = m(2) * t864 + t706;
t838 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t872 + Ifges(4,6) * t877) * qJD(2);
t703 = mrSges(4,2) * t809 - mrSges(4,3) * t802 + Ifges(4,1) * t852 + Ifges(4,4) * t853 + Ifges(4,5) * qJDD(3) - pkin(9) * t728 - qJD(3) * t839 - t713 * t871 + t714 * t876 + t838 * t898;
t886 = -mrSges(7,1) * t754 + mrSges(7,2) * t755 - Ifges(7,5) * t768 - Ifges(7,6) * t767 - Ifges(7,3) * t836 - t805 * t779 + t804 * t780;
t883 = -mrSges(6,1) * t759 + mrSges(6,2) * t760 - Ifges(6,5) * t789 - Ifges(6,6) * t788 - Ifges(6,3) * t841 - pkin(5) * t742 - t826 * t799 + t825 * t800 + t886;
t881 = mrSges(5,1) * t775 - mrSges(5,2) * t776 + Ifges(5,5) * t823 + Ifges(5,6) * t822 + Ifges(5,3) * t845 + pkin(4) * t735 + t849 * t817 - t848 * t818 - t883;
t712 = -mrSges(4,1) * t809 + mrSges(4,3) * t803 + Ifges(4,4) * t852 + Ifges(4,2) * t853 + Ifges(4,6) * qJDD(3) - pkin(3) * t728 + qJD(3) * t840 - t838 * t899 - t881;
t701 = mrSges(3,2) * t832 - mrSges(3,3) * t814 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t880 - pkin(8) * t721 + t703 * t877 - t712 * t872;
t702 = -mrSges(3,1) * t832 + mrSges(3,3) * t815 + t880 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t721 - t904;
t887 = pkin(7) * t711 + t701 * t873 + t702 * t878;
t700 = mrSges(3,1) * t814 - mrSges(3,2) * t815 + Ifges(3,3) * qJDD(2) + pkin(2) * t884 + pkin(8) * t894 + t872 * t703 + t877 * t712;
t699 = mrSges(2,2) * t864 - mrSges(2,3) * t854 + t878 * t701 - t873 * t702 + (-t706 * t866 - t707 * t868) * pkin(7);
t698 = -mrSges(2,1) * t864 + mrSges(2,3) * t855 - pkin(1) * t706 - t866 * t700 + t868 * t887;
t1 = [-m(1) * g(1) + t895; -m(1) * g(2) + t900; -m(1) * g(3) + t891; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t900 - t865 * t698 + t867 * t699; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t895 + t867 * t698 + t865 * t699; -mrSges(1,1) * g(2) + mrSges(2,1) * t854 + mrSges(1,2) * g(1) - mrSges(2,2) * t855 + pkin(1) * t707 + t868 * t700 + t866 * t887; t891; t700; t904; t881; -t883; -t886;];
tauJB  = t1;
