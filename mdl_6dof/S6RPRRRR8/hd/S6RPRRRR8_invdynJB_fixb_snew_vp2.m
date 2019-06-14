% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:15:04
% EndTime: 2019-05-06 04:15:14
% DurationCPUTime: 9.41s
% Computational Cost: add. (163290->340), mult. (319575->415), div. (0->0), fcn. (219156->10), ass. (0->138)
t866 = sin(qJ(1));
t871 = cos(qJ(1));
t844 = -t871 * g(1) - t866 * g(2);
t885 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t844;
t901 = -pkin(1) - pkin(7);
t900 = mrSges(2,1) - mrSges(3,2);
t899 = Ifges(2,5) - Ifges(3,4);
t898 = (-Ifges(2,6) + Ifges(3,5));
t843 = t866 * g(1) - t871 * g(2);
t872 = qJD(1) ^ 2;
t884 = -t872 * qJ(2) + qJDD(2) - t843;
t815 = t901 * qJDD(1) + t884;
t865 = sin(qJ(3));
t870 = cos(qJ(3));
t804 = t865 * g(3) + t870 * t815;
t894 = qJD(1) * qJD(3);
t892 = t865 * t894;
t838 = t870 * qJDD(1) - t892;
t777 = (-t838 - t892) * pkin(8) + (-t865 * t870 * t872 + qJDD(3)) * pkin(3) + t804;
t805 = -t870 * g(3) + t865 * t815;
t837 = -t865 * qJDD(1) - t870 * t894;
t895 = qJD(1) * t870;
t842 = qJD(3) * pkin(3) - pkin(8) * t895;
t859 = t865 ^ 2;
t780 = -t859 * t872 * pkin(3) + t837 * pkin(8) - qJD(3) * t842 + t805;
t864 = sin(qJ(4));
t869 = cos(qJ(4));
t763 = t864 * t777 + t869 * t780;
t828 = (-t864 * t865 + t869 * t870) * qJD(1);
t792 = -t828 * qJD(4) + t869 * t837 - t864 * t838;
t896 = qJD(1) * t865;
t827 = -t864 * t895 - t869 * t896;
t801 = -t827 * mrSges(5,1) + t828 * mrSges(5,2);
t853 = qJD(3) + qJD(4);
t813 = t853 * mrSges(5,1) - t828 * mrSges(5,3);
t852 = qJDD(3) + qJDD(4);
t784 = -t837 * pkin(3) + t842 * t895 + (-pkin(8) * t859 + t901) * t872 + t885;
t793 = t827 * qJD(4) + t864 * t837 + t869 * t838;
t753 = (-t827 * t853 - t793) * pkin(9) + (t828 * t853 - t792) * pkin(4) + t784;
t802 = -t827 * pkin(4) - t828 * pkin(9);
t851 = t853 ^ 2;
t756 = -t851 * pkin(4) + t852 * pkin(9) + t827 * t802 + t763;
t863 = sin(qJ(5));
t868 = cos(qJ(5));
t742 = t868 * t753 - t863 * t756;
t807 = -t863 * t828 + t868 * t853;
t767 = t807 * qJD(5) + t868 * t793 + t863 * t852;
t791 = qJDD(5) - t792;
t808 = t868 * t828 + t863 * t853;
t823 = qJD(5) - t827;
t740 = (t807 * t823 - t767) * pkin(10) + (t807 * t808 + t791) * pkin(5) + t742;
t743 = t863 * t753 + t868 * t756;
t766 = -t808 * qJD(5) - t863 * t793 + t868 * t852;
t796 = t823 * pkin(5) - t808 * pkin(10);
t806 = t807 ^ 2;
t741 = -t806 * pkin(5) + t766 * pkin(10) - t823 * t796 + t743;
t862 = sin(qJ(6));
t867 = cos(qJ(6));
t738 = t867 * t740 - t862 * t741;
t778 = t867 * t807 - t862 * t808;
t749 = t778 * qJD(6) + t862 * t766 + t867 * t767;
t779 = t862 * t807 + t867 * t808;
t764 = -t778 * mrSges(7,1) + t779 * mrSges(7,2);
t820 = qJD(6) + t823;
t768 = -t820 * mrSges(7,2) + t778 * mrSges(7,3);
t786 = qJDD(6) + t791;
t734 = m(7) * t738 + t786 * mrSges(7,1) - t749 * mrSges(7,3) - t779 * t764 + t820 * t768;
t739 = t862 * t740 + t867 * t741;
t748 = -t779 * qJD(6) + t867 * t766 - t862 * t767;
t769 = t820 * mrSges(7,1) - t779 * mrSges(7,3);
t735 = m(7) * t739 - t786 * mrSges(7,2) + t748 * mrSges(7,3) + t778 * t764 - t820 * t769;
t726 = t867 * t734 + t862 * t735;
t781 = -t807 * mrSges(6,1) + t808 * mrSges(6,2);
t794 = -t823 * mrSges(6,2) + t807 * mrSges(6,3);
t724 = m(6) * t742 + t791 * mrSges(6,1) - t767 * mrSges(6,3) - t808 * t781 + t823 * t794 + t726;
t795 = t823 * mrSges(6,1) - t808 * mrSges(6,3);
t887 = -t862 * t734 + t867 * t735;
t725 = m(6) * t743 - t791 * mrSges(6,2) + t766 * mrSges(6,3) + t807 * t781 - t823 * t795 + t887;
t888 = -t863 * t724 + t868 * t725;
t718 = m(5) * t763 - t852 * mrSges(5,2) + t792 * mrSges(5,3) + t827 * t801 - t853 * t813 + t888;
t762 = t869 * t777 - t864 * t780;
t812 = -t853 * mrSges(5,2) + t827 * mrSges(5,3);
t755 = -t852 * pkin(4) - t851 * pkin(9) + t828 * t802 - t762;
t744 = -t766 * pkin(5) - t806 * pkin(10) + t808 * t796 + t755;
t882 = m(7) * t744 - t748 * mrSges(7,1) + t749 * mrSges(7,2) - t778 * t768 + t779 * t769;
t876 = -m(6) * t755 + t766 * mrSges(6,1) - t767 * mrSges(6,2) + t807 * t794 - t808 * t795 - t882;
t730 = m(5) * t762 + t852 * mrSges(5,1) - t793 * mrSges(5,3) - t828 * t801 + t853 * t812 + t876;
t709 = t864 * t718 + t869 * t730;
t836 = (mrSges(4,1) * t865 + mrSges(4,2) * t870) * qJD(1);
t840 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t896;
t704 = m(4) * t804 + qJDD(3) * mrSges(4,1) - t838 * mrSges(4,3) + qJD(3) * t840 - t836 * t895 + t709;
t841 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t895;
t889 = t869 * t718 - t864 * t730;
t705 = m(4) * t805 - qJDD(3) * mrSges(4,2) + t837 * mrSges(4,3) - qJD(3) * t841 - t836 * t896 + t889;
t701 = t870 * t704 + t865 * t705;
t822 = -qJDD(1) * pkin(1) + t884;
t883 = -m(3) * t822 + (t872 * mrSges(3,3)) - t701;
t697 = m(2) * t843 - (t872 * mrSges(2,2)) + t900 * qJDD(1) + t883;
t818 = t872 * pkin(1) - t885;
t814 = t901 * t872 + t885;
t720 = t868 * t724 + t863 * t725;
t881 = m(5) * t784 - t792 * mrSges(5,1) + t793 * mrSges(5,2) - t827 * t812 + t828 * t813 + t720;
t878 = -m(4) * t814 + t837 * mrSges(4,1) - t838 * mrSges(4,2) - t840 * t896 - t841 * t895 - t881;
t875 = -m(3) * t818 + (t872 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t878;
t714 = m(2) * t844 - (t872 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t875;
t897 = t871 * t697 + t866 * t714;
t891 = -t866 * t697 + t871 * t714;
t890 = -t865 * t704 + t870 * t705;
t758 = Ifges(7,4) * t779 + Ifges(7,2) * t778 + Ifges(7,6) * t820;
t759 = Ifges(7,1) * t779 + Ifges(7,4) * t778 + Ifges(7,5) * t820;
t880 = -mrSges(7,1) * t738 + mrSges(7,2) * t739 - Ifges(7,5) * t749 - Ifges(7,6) * t748 - Ifges(7,3) * t786 - t779 * t758 + t778 * t759;
t757 = Ifges(7,5) * t779 + Ifges(7,6) * t778 + Ifges(7,3) * t820;
t727 = -mrSges(7,1) * t744 + mrSges(7,3) * t739 + Ifges(7,4) * t749 + Ifges(7,2) * t748 + Ifges(7,6) * t786 - t779 * t757 + t820 * t759;
t728 = mrSges(7,2) * t744 - mrSges(7,3) * t738 + Ifges(7,1) * t749 + Ifges(7,4) * t748 + Ifges(7,5) * t786 + t778 * t757 - t820 * t758;
t770 = Ifges(6,5) * t808 + Ifges(6,6) * t807 + Ifges(6,3) * t823;
t772 = Ifges(6,1) * t808 + Ifges(6,4) * t807 + Ifges(6,5) * t823;
t707 = -mrSges(6,1) * t755 + mrSges(6,3) * t743 + Ifges(6,4) * t767 + Ifges(6,2) * t766 + Ifges(6,6) * t791 - pkin(5) * t882 + pkin(10) * t887 + t867 * t727 + t862 * t728 - t808 * t770 + t823 * t772;
t771 = Ifges(6,4) * t808 + Ifges(6,2) * t807 + Ifges(6,6) * t823;
t711 = mrSges(6,2) * t755 - mrSges(6,3) * t742 + Ifges(6,1) * t767 + Ifges(6,4) * t766 + Ifges(6,5) * t791 - pkin(10) * t726 - t862 * t727 + t867 * t728 + t807 * t770 - t823 * t771;
t798 = Ifges(5,4) * t828 + Ifges(5,2) * t827 + Ifges(5,6) * t853;
t799 = Ifges(5,1) * t828 + Ifges(5,4) * t827 + Ifges(5,5) * t853;
t879 = mrSges(5,1) * t762 - mrSges(5,2) * t763 + Ifges(5,5) * t793 + Ifges(5,6) * t792 + Ifges(5,3) * t852 + pkin(4) * t876 + pkin(9) * t888 + t868 * t707 + t863 * t711 + t828 * t798 - t827 * t799;
t797 = Ifges(5,5) * t828 + Ifges(5,6) * t827 + Ifges(5,3) * t853;
t695 = mrSges(5,2) * t784 - mrSges(5,3) * t762 + Ifges(5,1) * t793 + Ifges(5,4) * t792 + Ifges(5,5) * t852 - pkin(9) * t720 - t863 * t707 + t868 * t711 + t827 * t797 - t853 * t798;
t873 = mrSges(6,1) * t742 - mrSges(6,2) * t743 + Ifges(6,5) * t767 + Ifges(6,6) * t766 + Ifges(6,3) * t791 + pkin(5) * t726 + t808 * t771 - t807 * t772 - t880;
t702 = -mrSges(5,1) * t784 + mrSges(5,3) * t763 + Ifges(5,4) * t793 + Ifges(5,2) * t792 + Ifges(5,6) * t852 - pkin(4) * t720 - t828 * t797 + t853 * t799 - t873;
t824 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t870 - Ifges(4,6) * t865) * qJD(1);
t826 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t870 - Ifges(4,4) * t865) * qJD(1);
t692 = -mrSges(4,1) * t814 + mrSges(4,3) * t805 + Ifges(4,4) * t838 + Ifges(4,2) * t837 + Ifges(4,6) * qJDD(3) - pkin(3) * t881 + pkin(8) * t889 + qJD(3) * t826 + t864 * t695 + t869 * t702 - t824 * t895;
t825 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t870 - Ifges(4,2) * t865) * qJD(1);
t694 = mrSges(4,2) * t814 - mrSges(4,3) * t804 + Ifges(4,1) * t838 + Ifges(4,4) * t837 + Ifges(4,5) * qJDD(3) - pkin(8) * t709 - qJD(3) * t825 + t869 * t695 - t864 * t702 - t824 * t896;
t699 = qJDD(1) * mrSges(3,2) - t883;
t877 = mrSges(2,1) * t843 - mrSges(2,2) * t844 + mrSges(3,2) * t822 - mrSges(3,3) * t818 - pkin(1) * t699 - pkin(7) * t701 + qJ(2) * t875 - t865 * t692 + t870 * t694 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t874 = mrSges(4,1) * t804 - mrSges(4,2) * t805 + Ifges(4,5) * t838 + Ifges(4,6) * t837 + Ifges(4,3) * qJDD(3) + pkin(3) * t709 + t825 * t895 + t826 * t896 + t879;
t700 = -m(3) * g(3) + t890;
t691 = (t898 * t872) - mrSges(2,3) * t843 + mrSges(3,1) * t822 + t899 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + pkin(2) * t701 - qJ(2) * t700 + t874;
t690 = -mrSges(3,1) * t818 + mrSges(2,3) * t844 - pkin(1) * t700 - pkin(2) * t878 - pkin(7) * t890 + t900 * g(3) - t898 * qJDD(1) - t870 * t692 - t865 * t694 + t899 * t872;
t1 = [-m(1) * g(1) + t891; -m(1) * g(2) + t897; (-m(1) - m(2) - m(3)) * g(3) + t890; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t897 - t866 * t690 + t871 * t691; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t891 + t871 * t690 + t866 * t691; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t877; t877; t699; t874; t879; t873; -t880;];
tauJB  = t1;
