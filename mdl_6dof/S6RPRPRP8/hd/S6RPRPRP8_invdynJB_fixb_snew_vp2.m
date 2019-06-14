% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:00:54
% EndTime: 2019-05-05 18:01:01
% DurationCPUTime: 5.24s
% Computational Cost: add. (58319->314), mult. (125990->373), div. (0->0), fcn. (80399->8), ass. (0->126)
t895 = -2 * qJD(4);
t894 = Ifges(6,1) + Ifges(7,1);
t886 = Ifges(6,4) - Ifges(7,5);
t884 = -Ifges(6,5) - Ifges(7,4);
t893 = Ifges(6,2) + Ifges(7,3);
t882 = Ifges(6,6) - Ifges(7,6);
t892 = -Ifges(6,3) - Ifges(7,2);
t847 = sin(qJ(1));
t849 = cos(qJ(1));
t826 = t847 * g(1) - t849 * g(2);
t851 = qJD(1) ^ 2;
t859 = -t851 * qJ(2) + qJDD(2) - t826;
t890 = -pkin(1) - pkin(7);
t798 = qJDD(1) * t890 + t859;
t846 = sin(qJ(3));
t848 = cos(qJ(3));
t788 = t846 * g(3) + t848 * t798;
t872 = qJD(1) * qJD(3);
t869 = t846 * t872;
t821 = qJDD(1) * t848 - t869;
t760 = (-t821 - t869) * qJ(4) + (-t846 * t848 * t851 + qJDD(3)) * pkin(3) + t788;
t789 = -g(3) * t848 + t846 * t798;
t820 = -qJDD(1) * t846 - t848 * t872;
t874 = qJD(1) * t848;
t824 = qJD(3) * pkin(3) - qJ(4) * t874;
t840 = t846 ^ 2;
t761 = -pkin(3) * t840 * t851 + qJ(4) * t820 - qJD(3) * t824 + t789;
t843 = sin(pkin(9));
t844 = cos(pkin(9));
t875 = qJD(1) * t846;
t809 = -t843 * t875 + t844 * t874;
t741 = t844 * t760 - t843 * t761 + t809 * t895;
t827 = -t849 * g(1) - t847 * g(2);
t860 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t827;
t808 = (t843 * t848 + t844 * t846) * qJD(1);
t786 = t820 * t843 + t821 * t844;
t845 = sin(qJ(5));
t889 = cos(qJ(5));
t790 = -qJD(3) * t889 + t809 * t845;
t754 = -t790 * qJD(5) + t845 * qJDD(3) + t786 * t889;
t791 = t845 * qJD(3) + t809 * t889;
t766 = mrSges(7,1) * t790 - mrSges(7,3) * t791;
t742 = t843 * t760 + t844 * t761 + t808 * t895;
t778 = pkin(4) * t808 - pkin(8) * t809;
t850 = qJD(3) ^ 2;
t737 = -pkin(4) * t850 + qJDD(3) * pkin(8) - t778 * t808 + t742;
t763 = -t820 * pkin(3) + qJDD(4) + t824 * t874 + (-qJ(4) * t840 + t890) * t851 + t860;
t785 = t820 * t844 - t821 * t843;
t739 = (qJD(3) * t808 - t786) * pkin(8) + (qJD(3) * t809 - t785) * pkin(4) + t763;
t733 = -t845 * t737 + t739 * t889;
t765 = pkin(5) * t790 - qJ(6) * t791;
t784 = qJDD(5) - t785;
t806 = qJD(5) + t808;
t805 = t806 ^ 2;
t731 = -t784 * pkin(5) - t805 * qJ(6) + t791 * t765 + qJDD(6) - t733;
t769 = -mrSges(7,2) * t790 + mrSges(7,3) * t806;
t863 = -m(7) * t731 + t784 * mrSges(7,1) + t806 * t769;
t727 = t754 * mrSges(7,2) + t791 * t766 - t863;
t734 = t889 * t737 + t845 * t739;
t730 = -pkin(5) * t805 + t784 * qJ(6) + 0.2e1 * qJD(6) * t806 - t765 * t790 + t734;
t753 = qJD(5) * t791 - qJDD(3) * t889 + t786 * t845;
t772 = -mrSges(7,1) * t806 + mrSges(7,2) * t791;
t870 = m(7) * t730 + t784 * mrSges(7,3) + t806 * t772;
t877 = t790 * t886 - t791 * t894 + t806 * t884;
t878 = t790 * t893 - t791 * t886 - t806 * t882;
t891 = -t753 * t882 - t754 * t884 - t892 * t784 - t790 * t877 - t791 * t878 + mrSges(6,1) * t733 - mrSges(7,1) * t731 - mrSges(6,2) * t734 + mrSges(7,3) * t730 - pkin(5) * t727 + qJ(6) * (-t753 * mrSges(7,2) - t790 * t766 + t870);
t888 = mrSges(2,1) - mrSges(3,2);
t887 = -mrSges(6,3) - mrSges(7,2);
t885 = Ifges(2,5) - Ifges(3,4);
t883 = -Ifges(2,6) + Ifges(3,5);
t777 = mrSges(5,1) * t808 + mrSges(5,2) * t809;
t797 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t809;
t771 = mrSges(6,1) * t806 - mrSges(6,3) * t791;
t876 = -mrSges(6,1) * t790 - mrSges(6,2) * t791 - t766;
t722 = m(6) * t734 - t784 * mrSges(6,2) + t753 * t887 - t806 * t771 + t790 * t876 + t870;
t770 = -mrSges(6,2) * t806 - mrSges(6,3) * t790;
t724 = m(6) * t733 + t784 * mrSges(6,1) + t754 * t887 + t806 * t770 + t791 * t876 + t863;
t865 = t889 * t722 - t724 * t845;
t710 = m(5) * t742 - qJDD(3) * mrSges(5,2) + t785 * mrSges(5,3) - qJD(3) * t797 - t777 * t808 + t865;
t796 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t808;
t736 = -qJDD(3) * pkin(4) - t850 * pkin(8) + t809 * t778 - t741;
t732 = -0.2e1 * qJD(6) * t791 + (t790 * t806 - t754) * qJ(6) + (t791 * t806 + t753) * pkin(5) + t736;
t728 = m(7) * t732 + t753 * mrSges(7,1) - t754 * mrSges(7,3) + t769 * t790 - t791 * t772;
t854 = -m(6) * t736 - t753 * mrSges(6,1) - t754 * mrSges(6,2) - t790 * t770 - t771 * t791 - t728;
t719 = m(5) * t741 + qJDD(3) * mrSges(5,1) - t786 * mrSges(5,3) + qJD(3) * t796 - t777 * t809 + t854;
t703 = t843 * t710 + t844 * t719;
t819 = (mrSges(4,1) * t846 + mrSges(4,2) * t848) * qJD(1);
t823 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t875;
t700 = m(4) * t788 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t821 + qJD(3) * t823 - t819 * t874 + t703;
t825 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t874;
t866 = t844 * t710 - t719 * t843;
t701 = m(4) * t789 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t820 - qJD(3) * t825 - t819 * t875 + t866;
t696 = t848 * t700 + t846 * t701;
t807 = -qJDD(1) * pkin(1) + t859;
t858 = -m(3) * t807 + t851 * mrSges(3,3) - t696;
t692 = m(2) * t826 - t851 * mrSges(2,2) + qJDD(1) * t888 + t858;
t801 = t851 * pkin(1) - t860;
t717 = t845 * t722 + t889 * t724;
t715 = m(5) * t763 - t785 * mrSges(5,1) + t786 * mrSges(5,2) + t796 * t808 + t809 * t797 + t717;
t795 = t851 * t890 + t860;
t856 = -m(4) * t795 + mrSges(4,1) * t820 - t821 * mrSges(4,2) - t823 * t875 - t825 * t874 - t715;
t853 = -m(3) * t801 + t851 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t856;
t706 = m(2) * t827 - mrSges(2,1) * t851 - qJDD(1) * mrSges(2,2) + t853;
t880 = t849 * t692 + t847 * t706;
t879 = t790 * t882 + t791 * t884 + t806 * t892;
t868 = -t692 * t847 + t849 * t706;
t867 = -t846 * t700 + t848 * t701;
t712 = -mrSges(6,1) * t736 - mrSges(7,1) * t732 + mrSges(7,2) * t730 + mrSges(6,3) * t734 - pkin(5) * t728 - t753 * t893 + t886 * t754 + t882 * t784 + t879 * t791 - t877 * t806;
t714 = mrSges(6,2) * t736 + mrSges(7,2) * t731 - mrSges(6,3) * t733 - mrSges(7,3) * t732 - qJ(6) * t728 - t886 * t753 + t754 * t894 - t884 * t784 + t879 * t790 + t878 * t806;
t773 = Ifges(5,5) * t809 - Ifges(5,6) * t808 + Ifges(5,3) * qJD(3);
t774 = Ifges(5,4) * t809 - Ifges(5,2) * t808 + Ifges(5,6) * qJD(3);
t697 = mrSges(5,2) * t763 - mrSges(5,3) * t741 + Ifges(5,1) * t786 + Ifges(5,4) * t785 + Ifges(5,5) * qJDD(3) - pkin(8) * t717 - qJD(3) * t774 - t845 * t712 + t714 * t889 - t808 * t773;
t775 = Ifges(5,1) * t809 - Ifges(5,4) * t808 + Ifges(5,5) * qJD(3);
t698 = -mrSges(5,1) * t763 + mrSges(5,3) * t742 + Ifges(5,4) * t786 + Ifges(5,2) * t785 + Ifges(5,6) * qJDD(3) - pkin(4) * t717 + qJD(3) * t775 - t809 * t773 - t891;
t810 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t848 - Ifges(4,6) * t846) * qJD(1);
t812 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t848 - Ifges(4,4) * t846) * qJD(1);
t688 = -mrSges(4,1) * t795 + mrSges(4,3) * t789 + Ifges(4,4) * t821 + Ifges(4,2) * t820 + Ifges(4,6) * qJDD(3) - pkin(3) * t715 + qJ(4) * t866 + qJD(3) * t812 + t843 * t697 + t844 * t698 - t810 * t874;
t811 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t848 - Ifges(4,2) * t846) * qJD(1);
t690 = mrSges(4,2) * t795 - mrSges(4,3) * t788 + Ifges(4,1) * t821 + Ifges(4,4) * t820 + Ifges(4,5) * qJDD(3) - qJ(4) * t703 - qJD(3) * t811 + t697 * t844 - t698 * t843 - t810 * t875;
t694 = qJDD(1) * mrSges(3,2) - t858;
t855 = mrSges(2,1) * t826 - mrSges(2,2) * t827 + mrSges(3,2) * t807 - mrSges(3,3) * t801 - pkin(1) * t694 - pkin(7) * t696 + qJ(2) * t853 - t688 * t846 + t848 * t690 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t852 = mrSges(5,1) * t741 - mrSges(4,2) * t789 - mrSges(5,2) * t742 + Ifges(5,5) * t786 + Ifges(5,6) * t785 + pkin(3) * t703 + pkin(4) * t854 + pkin(8) * t865 + t889 * t712 + t845 * t714 + t809 * t774 + t808 * t775 + mrSges(4,1) * t788 + t812 * t875 + t811 * t874 + Ifges(4,6) * t820 + Ifges(4,5) * t821 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t695 = -m(3) * g(3) + t867;
t687 = t852 + t885 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t883 * t851 - mrSges(2,3) * t826 + mrSges(3,1) * t807 - qJ(2) * t695 + pkin(2) * t696;
t686 = -mrSges(3,1) * t801 + mrSges(2,3) * t827 - pkin(1) * t695 - pkin(2) * t856 - pkin(7) * t867 + g(3) * t888 - qJDD(1) * t883 - t848 * t688 - t846 * t690 + t851 * t885;
t1 = [-m(1) * g(1) + t868; -m(1) * g(2) + t880; (-m(1) - m(2) - m(3)) * g(3) + t867; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t880 - t847 * t686 + t849 * t687; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t868 + t849 * t686 + t847 * t687; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t855; t855; t694; t852; t715; t891; t727;];
tauJB  = t1;
