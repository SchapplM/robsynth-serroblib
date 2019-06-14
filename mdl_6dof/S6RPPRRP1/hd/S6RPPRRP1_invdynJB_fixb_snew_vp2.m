% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:44:02
% EndTime: 2019-05-05 14:44:10
% DurationCPUTime: 7.76s
% Computational Cost: add. (91950->300), mult. (203637->358), div. (0->0), fcn. (139653->10), ass. (0->133)
t884 = Ifges(6,1) + Ifges(7,1);
t878 = Ifges(6,4) + Ifges(7,4);
t877 = Ifges(6,5) + Ifges(7,5);
t883 = Ifges(6,2) + Ifges(7,2);
t876 = Ifges(6,6) + Ifges(7,6);
t882 = Ifges(6,3) + Ifges(7,3);
t842 = qJD(1) ^ 2;
t831 = sin(pkin(10));
t833 = cos(pkin(10));
t836 = sin(qJ(4));
t839 = cos(qJ(4));
t849 = t831 * t836 - t833 * t839;
t805 = t849 * qJD(1);
t850 = t831 * t839 + t833 * t836;
t806 = t850 * qJD(1);
t864 = t806 * qJD(4);
t794 = -t849 * qJDD(1) - t864;
t865 = t805 * qJD(4);
t795 = t850 * qJDD(1) - t865;
t835 = sin(qJ(5));
t838 = cos(qJ(5));
t800 = qJD(4) * t838 - t806 * t835;
t768 = qJD(5) * t800 + qJDD(4) * t835 + t795 * t838;
t801 = qJD(4) * t835 + t806 * t838;
t774 = -mrSges(7,1) * t800 + mrSges(7,2) * t801;
t837 = sin(qJ(1));
t840 = cos(qJ(1));
t815 = t837 * g(1) - g(2) * t840;
t812 = qJDD(1) * pkin(1) + t815;
t816 = -g(1) * t840 - g(2) * t837;
t813 = -pkin(1) * t842 + t816;
t832 = sin(pkin(9));
t834 = cos(pkin(9));
t798 = t832 * t812 + t834 * t813;
t787 = -pkin(2) * t842 + qJDD(1) * qJ(3) + t798;
t830 = -g(3) + qJDD(2);
t863 = qJD(1) * qJD(3);
t867 = t833 * t830 - 0.2e1 * t831 * t863;
t880 = pkin(3) * t833;
t766 = (-pkin(7) * qJDD(1) + t842 * t880 - t787) * t831 + t867;
t773 = t831 * t830 + (t787 + 0.2e1 * t863) * t833;
t862 = qJDD(1) * t833;
t827 = t833 ^ 2;
t873 = t827 * t842;
t769 = -pkin(3) * t873 + pkin(7) * t862 + t773;
t750 = t836 * t766 + t839 * t769;
t793 = pkin(4) * t805 - pkin(8) * t806;
t841 = qJD(4) ^ 2;
t745 = -pkin(4) * t841 + qJDD(4) * pkin(8) - t793 * t805 + t750;
t826 = t831 ^ 2;
t797 = t834 * t812 - t832 * t813;
t851 = qJDD(3) - t797;
t771 = (-pkin(2) - t880) * qJDD(1) + (-qJ(3) + (-t826 - t827) * pkin(7)) * t842 + t851;
t748 = (-t795 + t865) * pkin(8) + (-t794 + t864) * pkin(4) + t771;
t740 = -t835 * t745 + t838 * t748;
t792 = qJDD(5) - t794;
t804 = qJD(5) + t805;
t737 = -0.2e1 * qJD(6) * t801 + (t800 * t804 - t768) * qJ(6) + (t800 * t801 + t792) * pkin(5) + t740;
t777 = -mrSges(7,2) * t804 + mrSges(7,3) * t800;
t861 = m(7) * t737 + t792 * mrSges(7,1) + t804 * t777;
t734 = -t768 * mrSges(7,3) - t801 * t774 + t861;
t741 = t838 * t745 + t835 * t748;
t767 = -qJD(5) * t801 + qJDD(4) * t838 - t795 * t835;
t779 = pkin(5) * t804 - qJ(6) * t801;
t799 = t800 ^ 2;
t739 = -pkin(5) * t799 + qJ(6) * t767 + 0.2e1 * qJD(6) * t800 - t779 * t804 + t741;
t869 = t878 * t800 + t801 * t884 + t877 * t804;
t870 = -t800 * t883 - t801 * t878 - t804 * t876;
t881 = mrSges(6,1) * t740 + mrSges(7,1) * t737 - mrSges(6,2) * t741 - mrSges(7,2) * t739 + pkin(5) * t734 + t876 * t767 + t877 * t768 + t792 * t882 - t869 * t800 - t870 * t801;
t879 = -mrSges(6,2) - mrSges(7,2);
t874 = mrSges(4,2) * t831;
t775 = -mrSges(6,1) * t800 + mrSges(6,2) * t801;
t778 = -mrSges(6,2) * t804 + mrSges(6,3) * t800;
t727 = m(6) * t740 + t792 * mrSges(6,1) + t804 * t778 + (-t774 - t775) * t801 + (-mrSges(6,3) - mrSges(7,3)) * t768 + t861;
t860 = m(7) * t739 + t767 * mrSges(7,3) + t800 * t774;
t780 = mrSges(7,1) * t804 - mrSges(7,3) * t801;
t868 = -mrSges(6,1) * t804 + mrSges(6,3) * t801 - t780;
t730 = m(6) * t741 + t767 * mrSges(6,3) + t800 * t775 + t879 * t792 + t868 * t804 + t860;
t725 = -t727 * t835 + t838 * t730;
t790 = mrSges(5,1) * t805 + mrSges(5,2) * t806;
t803 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t806;
t722 = m(5) * t750 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t794 - qJD(4) * t803 - t790 * t805 + t725;
t749 = t766 * t839 - t836 * t769;
t744 = -qJDD(4) * pkin(4) - pkin(8) * t841 + t806 * t793 - t749;
t742 = -pkin(5) * t767 - qJ(6) * t799 + t779 * t801 + qJDD(6) + t744;
t855 = -m(7) * t742 + t767 * mrSges(7,1) + t800 * t777;
t733 = -m(6) * t744 + t767 * mrSges(6,1) + t879 * t768 + t800 * t778 + t868 * t801 + t855;
t802 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t805;
t732 = m(5) * t749 + qJDD(4) * mrSges(5,1) - t795 * mrSges(5,3) + qJD(4) * t802 - t806 * t790 + t733;
t714 = t836 * t722 + t839 * t732;
t772 = -t787 * t831 + t867;
t848 = mrSges(4,3) * qJDD(1) + t842 * (-mrSges(4,1) * t833 + t874);
t712 = m(4) * t772 - t848 * t831 + t714;
t856 = t839 * t722 - t836 * t732;
t713 = m(4) * t773 + t848 * t833 + t856;
t857 = -t712 * t831 + t833 * t713;
t704 = m(3) * t798 - mrSges(3,1) * t842 - qJDD(1) * mrSges(3,2) + t857;
t783 = -qJDD(1) * pkin(2) - t842 * qJ(3) + t851;
t724 = t838 * t727 + t835 * t730;
t847 = m(5) * t771 - t794 * mrSges(5,1) + t795 * mrSges(5,2) + t805 * t802 + t806 * t803 + t724;
t844 = -m(4) * t783 + mrSges(4,1) * t862 - t847 + (t826 * t842 + t873) * mrSges(4,3);
t717 = t844 - t842 * mrSges(3,2) + m(3) * t797 + (mrSges(3,1) - t874) * qJDD(1);
t700 = t832 * t704 + t834 * t717;
t697 = m(2) * t815 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t842 + t700;
t858 = t834 * t704 - t717 * t832;
t698 = m(2) * t816 - mrSges(2,1) * t842 - qJDD(1) * mrSges(2,2) + t858;
t872 = t840 * t697 + t837 * t698;
t707 = t833 * t712 + t831 * t713;
t871 = -t800 * t876 - t801 * t877 - t804 * t882;
t852 = Ifges(4,5) * t831 + Ifges(4,6) * t833;
t866 = t842 * t852;
t705 = m(3) * t830 + t707;
t859 = -t697 * t837 + t840 * t698;
t854 = Ifges(4,1) * t831 + Ifges(4,4) * t833;
t853 = Ifges(4,4) * t831 + Ifges(4,2) * t833;
t735 = t768 * mrSges(7,2) + t801 * t780 - t855;
t715 = -mrSges(6,1) * t744 + mrSges(6,3) * t741 - mrSges(7,1) * t742 + mrSges(7,3) * t739 - pkin(5) * t735 + qJ(6) * t860 + (-qJ(6) * t780 + t869) * t804 + t871 * t801 + (-mrSges(7,2) * qJ(6) + t876) * t792 + t878 * t768 + t883 * t767;
t723 = mrSges(6,2) * t744 + mrSges(7,2) * t742 - mrSges(6,3) * t740 - mrSges(7,3) * t737 - qJ(6) * t734 + t878 * t767 + t768 * t884 + t877 * t792 - t871 * t800 + t870 * t804;
t784 = Ifges(5,5) * t806 - Ifges(5,6) * t805 + Ifges(5,3) * qJD(4);
t785 = Ifges(5,4) * t806 - Ifges(5,2) * t805 + Ifges(5,6) * qJD(4);
t701 = mrSges(5,2) * t771 - mrSges(5,3) * t749 + Ifges(5,1) * t795 + Ifges(5,4) * t794 + Ifges(5,5) * qJDD(4) - pkin(8) * t724 - qJD(4) * t785 - t715 * t835 + t723 * t838 - t784 * t805;
t786 = Ifges(5,1) * t806 - Ifges(5,4) * t805 + Ifges(5,5) * qJD(4);
t708 = -mrSges(5,1) * t771 + mrSges(5,3) * t750 + Ifges(5,4) * t795 + Ifges(5,2) * t794 + Ifges(5,6) * qJDD(4) - pkin(4) * t724 + qJD(4) * t786 - t806 * t784 - t881;
t691 = -mrSges(4,1) * t783 + mrSges(4,3) * t773 - pkin(3) * t847 + pkin(7) * t856 + t853 * qJDD(1) + t836 * t701 + t839 * t708 - t831 * t866;
t693 = mrSges(4,2) * t783 - mrSges(4,3) * t772 - pkin(7) * t714 + t854 * qJDD(1) + t839 * t701 - t836 * t708 + t833 * t866;
t719 = qJDD(1) * t874 - t844;
t845 = mrSges(2,1) * t815 + mrSges(3,1) * t797 - mrSges(2,2) * t816 - mrSges(3,2) * t798 + pkin(1) * t700 - pkin(2) * t719 + qJ(3) * t857 + t833 * t691 + t831 * t693 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t843 = mrSges(5,1) * t749 - mrSges(5,2) * t750 + Ifges(5,5) * t795 + Ifges(5,6) * t794 + Ifges(5,3) * qJDD(4) + pkin(4) * t733 + pkin(8) * t725 + t838 * t715 + t835 * t723 + t806 * t785 + t805 * t786;
t689 = -t843 + (Ifges(3,6) - t852) * qJDD(1) - mrSges(3,1) * t830 + mrSges(3,3) * t798 + mrSges(4,2) * t773 - mrSges(4,1) * t772 - pkin(3) * t714 - pkin(2) * t707 + (-t831 * t853 + t833 * t854 + Ifges(3,5)) * t842;
t688 = mrSges(3,2) * t830 - mrSges(3,3) * t797 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t842 - qJ(3) * t707 - t691 * t831 + t693 * t833;
t687 = -mrSges(2,2) * g(3) - mrSges(2,3) * t815 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t842 - qJ(2) * t700 + t688 * t834 - t689 * t832;
t686 = mrSges(2,1) * g(3) + mrSges(2,3) * t816 + t842 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t705 + qJ(2) * t858 + t832 * t688 + t834 * t689;
t1 = [-m(1) * g(1) + t859; -m(1) * g(2) + t872; (-m(1) - m(2)) * g(3) + t705; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t872 - t837 * t686 + t840 * t687; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t859 + t840 * t686 + t837 * t687; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t845; t845; t705; t719; t843; t881; t735;];
tauJB  = t1;
