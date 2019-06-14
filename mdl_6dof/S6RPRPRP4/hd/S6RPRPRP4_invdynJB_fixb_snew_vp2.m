% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:42:55
% EndTime: 2019-05-05 17:43:01
% DurationCPUTime: 4.14s
% Computational Cost: add. (37345->304), mult. (72622->350), div. (0->0), fcn. (38033->8), ass. (0->130)
t936 = Ifges(6,1) + Ifges(7,1);
t917 = Ifges(6,4) - Ifges(7,5);
t932 = Ifges(7,4) + Ifges(6,5);
t935 = Ifges(6,2) + Ifges(7,3);
t930 = Ifges(6,6) - Ifges(7,6);
t934 = -2 * qJD(4);
t933 = Ifges(4,1) + Ifges(5,2);
t918 = Ifges(4,4) + Ifges(5,6);
t916 = Ifges(4,5) - Ifges(5,4);
t931 = Ifges(4,2) + Ifges(5,3);
t915 = Ifges(4,6) - Ifges(5,5);
t929 = Ifges(4,3) + Ifges(5,1);
t928 = Ifges(6,3) + Ifges(7,2);
t875 = sin(qJ(5));
t878 = cos(qJ(5));
t879 = cos(qJ(3));
t904 = qJD(1) * t879;
t839 = qJD(3) * t875 + t878 * t904;
t840 = qJD(3) * t878 - t875 * t904;
t876 = sin(qJ(3));
t905 = qJD(1) * t876;
t860 = qJD(5) + t905;
t927 = t935 * t839 - t917 * t840 - t930 * t860;
t926 = -t917 * t839 + t936 * t840 + t932 * t860;
t843 = (mrSges(5,2) * t879 - mrSges(5,3) * t876) * qJD(1);
t853 = -mrSges(5,1) * t904 - qJD(3) * mrSges(5,3);
t903 = qJD(1) * qJD(3);
t899 = t876 * t903;
t847 = qJDD(1) * t879 - t899;
t855 = pkin(4) * t905 - qJD(3) * pkin(8);
t871 = t879 ^ 2;
t882 = qJD(1) ^ 2;
t900 = t879 * t903;
t846 = qJDD(1) * t876 + t900;
t877 = sin(qJ(1));
t880 = cos(qJ(1));
t856 = t877 * g(1) - g(2) * t880;
t841 = qJDD(1) * pkin(1) + t856;
t857 = -g(1) * t880 - g(2) * t877;
t845 = -pkin(1) * t882 + t857;
t873 = sin(pkin(9));
t874 = cos(pkin(9));
t805 = t874 * t841 - t873 * t845;
t893 = -qJDD(1) * pkin(2) - t805;
t888 = pkin(3) * t899 + t905 * t934 + (-t846 - t900) * qJ(4) + t893;
t922 = -pkin(3) - pkin(8);
t772 = -t855 * t905 + (-pkin(4) * t871 - pkin(7)) * t882 + t922 * t847 + t888;
t872 = -g(3) + qJDD(2);
t806 = t873 * t841 + t874 * t845;
t787 = -pkin(2) * t882 + qJDD(1) * pkin(7) + t806;
t784 = t876 * t787;
t842 = (-pkin(3) * t879 - qJ(4) * t876) * qJD(1);
t881 = qJD(3) ^ 2;
t894 = -t881 * qJ(4) + t842 * t905 + qJDD(4) + t784;
t921 = pkin(8) * t882;
t776 = t846 * pkin(4) + t922 * qJDD(3) + (-pkin(4) * t903 - t876 * t921 - t872) * t879 + t894;
t770 = t878 * t772 + t875 * t776;
t802 = qJD(5) * t840 + qJDD(3) * t875 + t878 * t847;
t814 = mrSges(6,1) * t860 - mrSges(6,3) * t840;
t838 = qJDD(5) + t846;
t809 = pkin(5) * t839 - qJ(6) * t840;
t858 = t860 ^ 2;
t764 = -pkin(5) * t858 + qJ(6) * t838 + 0.2e1 * qJD(6) * t860 - t809 * t839 + t770;
t815 = -mrSges(7,1) * t860 + mrSges(7,2) * t840;
t901 = m(7) * t764 + t838 * mrSges(7,3) + t860 * t815;
t810 = mrSges(7,1) * t839 - mrSges(7,3) * t840;
t909 = -mrSges(6,1) * t839 - mrSges(6,2) * t840 - t810;
t919 = -mrSges(6,3) - mrSges(7,2);
t754 = m(6) * t770 - t838 * mrSges(6,2) + t919 * t802 - t860 * t814 + t909 * t839 + t901;
t769 = -t772 * t875 + t776 * t878;
t803 = -qJD(5) * t839 + qJDD(3) * t878 - t847 * t875;
t812 = -mrSges(6,2) * t860 - mrSges(6,3) * t839;
t765 = -pkin(5) * t838 - qJ(6) * t858 + t809 * t840 + qJDD(6) - t769;
t813 = -mrSges(7,2) * t839 + mrSges(7,3) * t860;
t895 = -m(7) * t765 + t838 * mrSges(7,1) + t860 * t813;
t756 = m(6) * t769 + t838 * mrSges(6,1) + t919 * t803 + t860 * t812 + t909 * t840 + t895;
t749 = t875 * t754 + t878 * t756;
t913 = t879 * t872;
t779 = -qJDD(3) * pkin(3) + t894 - t913;
t890 = -m(5) * t779 - t846 * mrSges(5,1) - t749;
t746 = qJDD(3) * mrSges(5,2) + qJD(3) * t853 + t843 * t905 - t890;
t783 = t879 * t787 + t876 * t872;
t778 = t881 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t934 - t842 * t904 - t783;
t775 = t847 * pkin(4) + qJD(3) * t855 - t871 * t921 - t778;
t767 = -0.2e1 * qJD(6) * t840 + (t839 * t860 - t803) * qJ(6) + (t840 * t860 + t802) * pkin(5) + t775;
t761 = m(7) * t767 + t802 * mrSges(7,1) - t803 * mrSges(7,3) + t839 * t813 - t840 * t815;
t910 = t930 * t839 - t932 * t840 - t928 * t860;
t747 = -mrSges(6,1) * t775 - mrSges(7,1) * t767 + mrSges(7,2) * t764 + mrSges(6,3) * t770 - pkin(5) * t761 - t935 * t802 + t917 * t803 + t930 * t838 + t910 * t840 + t926 * t860;
t748 = mrSges(6,2) * t775 + mrSges(7,2) * t765 - mrSges(6,3) * t769 - mrSges(7,3) * t767 - qJ(6) * t761 - t917 * t802 + t936 * t803 + t932 * t838 + t910 * t839 + t927 * t860;
t782 = -t784 + t913;
t854 = mrSges(5,1) * t905 + qJD(3) * mrSges(5,2);
t887 = m(6) * t775 + t802 * mrSges(6,1) + t803 * mrSges(6,2) + t839 * t812 + t840 * t814 + t761;
t886 = -m(5) * t778 + qJDD(3) * mrSges(5,3) + qJD(3) * t854 + t843 * t904 + t887;
t906 = t916 * qJD(3) + (t933 * t876 + t918 * t879) * qJD(1);
t907 = t915 * qJD(3) + (t918 * t876 + t931 * t879) * qJD(1);
t925 = (t907 * t876 - t906 * t879) * qJD(1) + t929 * qJDD(3) + t916 * t846 + t915 * t847 + mrSges(4,1) * t782 - mrSges(4,2) * t783 + mrSges(5,2) * t779 - mrSges(5,3) * t778 - pkin(3) * t746 - pkin(8) * t749 + qJ(4) * (mrSges(5,1) * t847 + t886) - t875 * t747 + t878 * t748;
t920 = t882 * pkin(7);
t844 = (-mrSges(4,1) * t879 + mrSges(4,2) * t876) * qJD(1);
t852 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t904;
t744 = m(4) * t782 - t846 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t852 - t853) * qJD(3) + (-t843 - t844) * t905 + t890;
t851 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t905;
t752 = t886 - qJD(3) * t851 + (mrSges(4,3) + mrSges(5,1)) * t847 - qJDD(3) * mrSges(4,2) + m(4) * t783 + t844 * t904;
t896 = -t744 * t876 + t879 * t752;
t736 = m(3) * t806 - mrSges(3,1) * t882 - qJDD(1) * mrSges(3,2) + t896;
t786 = t893 - t920;
t777 = -t847 * pkin(3) + t888 - t920;
t911 = t878 * t754 - t875 * t756;
t892 = -m(5) * t777 - t847 * mrSges(5,2) + t854 * t905 - t911;
t885 = -m(4) * t786 + t852 * t904 + t847 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t846 + (-t851 * t876 - t853 * t879) * qJD(1) + t892;
t741 = m(3) * t805 + qJDD(1) * mrSges(3,1) - t882 * mrSges(3,2) + t885;
t733 = t873 * t736 + t874 * t741;
t730 = m(2) * t856 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t882 + t733;
t897 = t874 * t736 - t741 * t873;
t731 = m(2) * t857 - mrSges(2,1) * t882 - qJDD(1) * mrSges(2,2) + t897;
t912 = t880 * t730 + t877 * t731;
t739 = t879 * t744 + t876 * t752;
t908 = t929 * qJD(3) + (t916 * t876 + t915 * t879) * qJD(1);
t737 = m(3) * t872 + t739;
t898 = -t730 * t877 + t880 * t731;
t745 = -t846 * mrSges(5,3) + t853 * t904 - t892;
t724 = -mrSges(4,1) * t786 - mrSges(5,1) * t778 + mrSges(5,2) * t777 + mrSges(4,3) * t783 - pkin(3) * t745 + pkin(4) * t887 - pkin(8) * t911 + t906 * qJD(3) + t915 * qJDD(3) - t878 * t747 - t875 * t748 + t918 * t846 + t931 * t847 - t908 * t905;
t760 = t803 * mrSges(7,2) + t840 * t810 - t895;
t883 = mrSges(6,1) * t769 - mrSges(7,1) * t765 - mrSges(6,2) * t770 + mrSges(7,3) * t764 - pkin(5) * t760 + qJ(6) * t901 - t927 * t840 + (-qJ(6) * t810 + t926) * t839 + t928 * t838 + t932 * t803 + (-mrSges(7,2) * qJ(6) - t930) * t802;
t726 = mrSges(5,1) * t779 + mrSges(4,2) * t786 - mrSges(4,3) * t782 - mrSges(5,3) * t777 + pkin(4) * t749 - qJ(4) * t745 - t907 * qJD(3) + t916 * qJDD(3) + t933 * t846 + t918 * t847 + t908 * t904 + t883;
t889 = mrSges(2,1) * t856 + mrSges(3,1) * t805 - mrSges(2,2) * t857 - mrSges(3,2) * t806 + pkin(1) * t733 + pkin(2) * t885 + pkin(7) * t896 + t879 * t724 + t876 * t726 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t722 = -mrSges(3,1) * t872 + mrSges(3,3) * t806 + t882 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t739 - t925;
t721 = mrSges(3,2) * t872 - mrSges(3,3) * t805 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t882 - pkin(7) * t739 - t724 * t876 + t726 * t879;
t720 = -mrSges(2,2) * g(3) - mrSges(2,3) * t856 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t882 - qJ(2) * t733 + t721 * t874 - t722 * t873;
t719 = mrSges(2,1) * g(3) + mrSges(2,3) * t857 + t882 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t737 + qJ(2) * t897 + t873 * t721 + t874 * t722;
t1 = [-m(1) * g(1) + t898; -m(1) * g(2) + t912; (-m(1) - m(2)) * g(3) + t737; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t912 - t877 * t719 + t880 * t720; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t898 + t880 * t719 + t877 * t720; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t889; t889; t737; t925; t746; t883; t760;];
tauJB  = t1;
