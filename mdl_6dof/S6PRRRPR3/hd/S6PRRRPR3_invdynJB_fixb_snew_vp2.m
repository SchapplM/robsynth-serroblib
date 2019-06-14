% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 07:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:28:17
% EndTime: 2019-05-05 07:28:27
% DurationCPUTime: 9.27s
% Computational Cost: add. (143144->322), mult. (285909->398), div. (0->0), fcn. (198067->12), ass. (0->141)
t938 = Ifges(5,4) + Ifges(6,6);
t949 = -Ifges(5,2) - Ifges(6,3);
t945 = Ifges(5,6) - Ifges(6,5);
t900 = sin(qJ(4));
t904 = cos(qJ(3));
t928 = qJD(2) * t904;
t901 = sin(qJ(3));
t929 = qJD(2) * t901;
t939 = cos(qJ(4));
t868 = t900 * t929 - t939 * t928;
t892 = qJD(3) + qJD(4);
t856 = t868 * mrSges(6,1) - t892 * mrSges(6,3);
t891 = qJDD(3) + qJDD(4);
t895 = sin(pkin(11));
t897 = cos(pkin(11));
t879 = t895 * g(1) - t897 * g(2);
t880 = -t897 * g(1) - t895 * g(2);
t894 = -g(3) + qJDD(1);
t905 = cos(qJ(2));
t898 = cos(pkin(6));
t902 = sin(qJ(2));
t934 = t898 * t902;
t896 = sin(pkin(6));
t935 = t896 * t902;
t842 = t879 * t934 + t905 * t880 + t894 * t935;
t906 = qJD(2) ^ 2;
t834 = -t906 * pkin(2) + qJDD(2) * pkin(8) + t842;
t859 = -t896 * t879 + t898 * t894;
t811 = -t901 * t834 + t904 * t859;
t927 = qJD(2) * qJD(3);
t926 = t904 * t927;
t877 = t901 * qJDD(2) + t926;
t798 = (-t877 + t926) * pkin(9) + (t901 * t904 * t906 + qJDD(3)) * pkin(3) + t811;
t812 = t904 * t834 + t901 * t859;
t878 = t904 * qJDD(2) - t901 * t927;
t884 = qJD(3) * pkin(3) - pkin(9) * t929;
t893 = t904 ^ 2;
t799 = -t893 * t906 * pkin(3) + t878 * pkin(9) - qJD(3) * t884 + t812;
t793 = t939 * t798 - t900 * t799;
t869 = (t900 * t904 + t939 * t901) * qJD(2);
t846 = t868 * pkin(4) - t869 * qJ(5);
t890 = t892 ^ 2;
t789 = -t891 * pkin(4) - t890 * qJ(5) + t869 * t846 + qJDD(5) - t793;
t830 = -t868 * qJD(4) + t939 * t877 + t900 * t878;
t936 = t868 * t892;
t783 = (t868 * t869 - t891) * pkin(10) + (t830 + t936) * pkin(5) + t789;
t829 = t869 * qJD(4) + t900 * t877 - t939 * t878;
t858 = t869 * pkin(5) - t892 * pkin(10);
t864 = t868 ^ 2;
t841 = -t902 * t880 + (t879 * t898 + t894 * t896) * t905;
t913 = -qJDD(2) * pkin(2) - t841;
t806 = -t878 * pkin(3) + t884 * t929 + (-pkin(9) * t893 - pkin(8)) * t906 + t913;
t940 = -2 * qJD(5);
t908 = (-t830 + t936) * qJ(5) + t806 + (t892 * pkin(4) + t940) * t869;
t786 = -t864 * pkin(5) - t869 * t858 + (pkin(4) + pkin(10)) * t829 + t908;
t899 = sin(qJ(6));
t903 = cos(qJ(6));
t781 = t903 * t783 - t899 * t786;
t850 = t903 * t868 - t899 * t892;
t805 = t850 * qJD(6) + t899 * t829 + t903 * t891;
t851 = t899 * t868 + t903 * t892;
t813 = -t850 * mrSges(7,1) + t851 * mrSges(7,2);
t827 = qJDD(6) + t830;
t862 = qJD(6) + t869;
t831 = -t862 * mrSges(7,2) + t850 * mrSges(7,3);
t778 = m(7) * t781 + t827 * mrSges(7,1) - t805 * mrSges(7,3) - t851 * t813 + t862 * t831;
t782 = t899 * t783 + t903 * t786;
t804 = -t851 * qJD(6) + t903 * t829 - t899 * t891;
t832 = t862 * mrSges(7,1) - t851 * mrSges(7,3);
t779 = m(7) * t782 - t827 * mrSges(7,2) + t804 * mrSges(7,3) + t850 * t813 - t862 * t832;
t767 = t903 * t778 + t899 * t779;
t848 = -t868 * mrSges(6,2) - t869 * mrSges(6,3);
t916 = -m(6) * t789 - t830 * mrSges(6,1) - t869 * t848 - t767;
t766 = t891 * mrSges(6,2) + t892 * t856 - t916;
t794 = t900 * t798 + t939 * t799;
t915 = -t890 * pkin(4) + t891 * qJ(5) - t868 * t846 + t794;
t785 = -t829 * pkin(5) - t864 * pkin(10) + ((2 * qJD(5)) + t858) * t892 + t915;
t807 = Ifges(7,5) * t851 + Ifges(7,6) * t850 + Ifges(7,3) * t862;
t809 = Ifges(7,1) * t851 + Ifges(7,4) * t850 + Ifges(7,5) * t862;
t769 = -mrSges(7,1) * t785 + mrSges(7,3) * t782 + Ifges(7,4) * t805 + Ifges(7,2) * t804 + Ifges(7,6) * t827 - t851 * t807 + t862 * t809;
t808 = Ifges(7,4) * t851 + Ifges(7,2) * t850 + Ifges(7,6) * t862;
t770 = mrSges(7,2) * t785 - mrSges(7,3) * t781 + Ifges(7,1) * t805 + Ifges(7,4) * t804 + Ifges(7,5) * t827 + t850 * t807 - t862 * t808;
t787 = t892 * t940 - t915;
t857 = t869 * mrSges(6,1) + t892 * mrSges(6,2);
t917 = -m(7) * t785 + t804 * mrSges(7,1) - t805 * mrSges(7,2) + t850 * t831 - t851 * t832;
t912 = -m(6) * t787 + t891 * mrSges(6,3) + t892 * t857 - t917;
t946 = Ifges(5,5) - Ifges(6,4);
t947 = Ifges(5,1) + Ifges(6,2);
t930 = -t938 * t868 + t869 * t947 + t946 * t892;
t943 = t949 * t868 + t938 * t869 + t945 * t892;
t944 = Ifges(5,3) + Ifges(6,1);
t948 = t930 * t868 - mrSges(5,2) * t794 - mrSges(6,3) * t787 - pkin(4) * t766 - pkin(10) * t767 - t899 * t769 + t903 * t770 + qJ(5) * (-t868 * t848 + t912) + mrSges(6,2) * t789 + mrSges(5,1) * t793 + t944 * t891 + t943 * t869 + t946 * t830 + (-qJ(5) * mrSges(6,1) - t945) * t829;
t847 = t868 * mrSges(5,1) + t869 * mrSges(5,2);
t854 = -t892 * mrSges(5,2) - t868 * mrSges(5,3);
t763 = m(5) * t793 - t830 * mrSges(5,3) - t869 * t847 + (t854 - t856) * t892 + (mrSges(5,1) - mrSges(6,2)) * t891 + t916;
t855 = t892 * mrSges(5,1) - t869 * mrSges(5,3);
t773 = m(5) * t794 - t891 * mrSges(5,2) - t892 * t855 + (-t847 - t848) * t868 + (-mrSges(5,3) - mrSges(6,1)) * t829 + t912;
t758 = t939 * t763 + t900 * t773;
t866 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t901 + Ifges(4,2) * t904) * qJD(2);
t867 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t901 + Ifges(4,4) * t904) * qJD(2);
t942 = mrSges(4,1) * t811 - mrSges(4,2) * t812 + Ifges(4,5) * t877 + Ifges(4,6) * t878 + Ifges(4,3) * qJDD(3) + pkin(3) * t758 + (t901 * t866 - t904 * t867) * qJD(2) + t948;
t833 = -t906 * pkin(8) + t913;
t881 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t929;
t882 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t928;
t791 = t829 * pkin(4) + t908;
t932 = -t899 * t778 + t903 * t779;
t765 = m(6) * t791 - t829 * mrSges(6,2) - t830 * mrSges(6,3) - t868 * t856 - t869 * t857 + t932;
t911 = m(5) * t806 + t829 * mrSges(5,1) + t830 * mrSges(5,2) + t868 * t854 + t869 * t855 + t765;
t909 = -m(4) * t833 + t878 * mrSges(4,1) - t877 * mrSges(4,2) - t881 * t929 + t882 * t928 - t911;
t761 = m(3) * t841 + qJDD(2) * mrSges(3,1) - t906 * mrSges(3,2) + t909;
t937 = t761 * t905;
t876 = (-mrSges(4,1) * t904 + mrSges(4,2) * t901) * qJD(2);
t756 = m(4) * t811 + qJDD(3) * mrSges(4,1) - t877 * mrSges(4,3) + qJD(3) * t882 - t876 * t929 + t758;
t922 = -t900 * t763 + t939 * t773;
t757 = m(4) * t812 - qJDD(3) * mrSges(4,2) + t878 * mrSges(4,3) - qJD(3) * t881 + t876 * t928 + t922;
t923 = -t901 * t756 + t904 * t757;
t748 = m(3) * t842 - t906 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t923;
t751 = t904 * t756 + t901 * t757;
t750 = m(3) * t859 + t751;
t739 = t748 * t934 - t896 * t750 + t898 * t937;
t737 = m(2) * t879 + t739;
t743 = t905 * t748 - t902 * t761;
t742 = m(2) * t880 + t743;
t933 = t897 * t737 + t895 * t742;
t931 = t868 * t945 - t869 * t946 - t892 * t944;
t738 = t748 * t935 + t898 * t750 + t896 * t937;
t924 = -t895 * t737 + t897 * t742;
t921 = m(2) * t894 + t738;
t744 = -mrSges(5,1) * t806 - mrSges(6,1) * t787 + mrSges(6,2) * t791 + mrSges(5,3) * t794 - pkin(4) * t765 - pkin(5) * t917 - pkin(10) * t932 - t903 * t769 - t899 * t770 + t949 * t829 + t938 * t830 + t931 * t869 + t945 * t891 + t930 * t892;
t914 = mrSges(7,1) * t781 - mrSges(7,2) * t782 + Ifges(7,5) * t805 + Ifges(7,6) * t804 + Ifges(7,3) * t827 + t851 * t808 - t850 * t809;
t752 = mrSges(6,1) * t789 + mrSges(5,2) * t806 - mrSges(5,3) * t793 - mrSges(6,3) * t791 + pkin(5) * t767 - qJ(5) * t765 - t938 * t829 + t830 * t947 + t931 * t868 + t946 * t891 - t943 * t892 + t914;
t865 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t901 + Ifges(4,6) * t904) * qJD(2);
t733 = -mrSges(4,1) * t833 + mrSges(4,3) * t812 + Ifges(4,4) * t877 + Ifges(4,2) * t878 + Ifges(4,6) * qJDD(3) - pkin(3) * t911 + pkin(9) * t922 + qJD(3) * t867 + t939 * t744 + t900 * t752 - t865 * t929;
t735 = mrSges(4,2) * t833 - mrSges(4,3) * t811 + Ifges(4,1) * t877 + Ifges(4,4) * t878 + Ifges(4,5) * qJDD(3) - pkin(9) * t758 - qJD(3) * t866 - t900 * t744 + t939 * t752 + t865 * t928;
t732 = mrSges(3,2) * t859 - mrSges(3,3) * t841 + Ifges(3,5) * qJDD(2) - t906 * Ifges(3,6) - pkin(8) * t751 - t901 * t733 + t904 * t735;
t734 = -mrSges(3,1) * t859 + mrSges(3,3) * t842 + t906 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t751 - t942;
t918 = pkin(7) * t743 + t732 * t902 + t734 * t905;
t731 = mrSges(3,1) * t841 - mrSges(3,2) * t842 + Ifges(3,3) * qJDD(2) + pkin(2) * t909 + pkin(8) * t923 + t904 * t733 + t901 * t735;
t730 = mrSges(2,2) * t894 - mrSges(2,3) * t879 + t905 * t732 - t902 * t734 + (-t738 * t896 - t739 * t898) * pkin(7);
t729 = -mrSges(2,1) * t894 + mrSges(2,3) * t880 - pkin(1) * t738 - t896 * t731 + t918 * t898;
t1 = [-m(1) * g(1) + t924; -m(1) * g(2) + t933; -m(1) * g(3) + t921; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t933 - t895 * t729 + t897 * t730; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t924 + t897 * t729 + t895 * t730; -mrSges(1,1) * g(2) + mrSges(2,1) * t879 + mrSges(1,2) * g(1) - mrSges(2,2) * t880 + pkin(1) * t739 + t898 * t731 + t918 * t896; t921; t731; t942; t948; t766; t914;];
tauJB  = t1;
