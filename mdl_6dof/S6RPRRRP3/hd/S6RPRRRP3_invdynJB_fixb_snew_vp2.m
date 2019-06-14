% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:18:21
% EndTime: 2019-05-06 01:18:31
% DurationCPUTime: 8.26s
% Computational Cost: add. (126897->320), mult. (244080->386), div. (0->0), fcn. (159521->10), ass. (0->131)
t876 = Ifges(6,1) + Ifges(7,1);
t866 = Ifges(6,4) - Ifges(7,5);
t874 = Ifges(7,4) + Ifges(6,5);
t875 = Ifges(6,2) + Ifges(7,3);
t873 = Ifges(6,6) - Ifges(7,6);
t872 = -Ifges(6,3) - Ifges(7,2);
t838 = sin(qJ(4));
t841 = cos(qJ(4));
t839 = sin(qJ(3));
t862 = qJD(1) * t839;
t812 = qJD(3) * t841 - t838 * t862;
t813 = qJD(3) * t838 + t841 * t862;
t837 = sin(qJ(5));
t868 = cos(qJ(5));
t787 = -t812 * t868 + t837 * t813;
t788 = t837 * t812 + t813 * t868;
t842 = cos(qJ(3));
t861 = qJD(1) * t842;
t827 = qJD(4) - t861;
t826 = qJD(5) + t827;
t871 = t875 * t787 - t866 * t788 - t873 * t826;
t870 = -t866 * t787 + t876 * t788 + t874 * t826;
t840 = sin(qJ(1));
t843 = cos(qJ(1));
t823 = t840 * g(1) - g(2) * t843;
t814 = qJDD(1) * pkin(1) + t823;
t824 = -g(1) * t843 - g(2) * t840;
t845 = qJD(1) ^ 2;
t816 = -pkin(1) * t845 + t824;
t835 = sin(pkin(10));
t836 = cos(pkin(10));
t790 = t835 * t814 + t836 * t816;
t777 = -pkin(2) * t845 + qJDD(1) * pkin(7) + t790;
t834 = -g(3) + qJDD(2);
t767 = -t839 * t777 + t842 * t834;
t817 = (-pkin(3) * t842 - pkin(8) * t839) * qJD(1);
t844 = qJD(3) ^ 2;
t760 = -qJDD(3) * pkin(3) - t844 * pkin(8) + t817 * t862 - t767;
t860 = qJD(1) * qJD(3);
t858 = t842 * t860;
t818 = qJDD(1) * t839 + t858;
t784 = -qJD(4) * t813 + qJDD(3) * t841 - t818 * t838;
t794 = pkin(4) * t827 - pkin(9) * t813;
t810 = t812 ^ 2;
t733 = -t784 * pkin(4) - t810 * pkin(9) + t813 * t794 + t760;
t785 = qJD(4) * t812 + qJDD(3) * t838 + t818 * t841;
t745 = t788 * qJD(5) - t784 * t868 + t837 * t785;
t746 = -t787 * qJD(5) + t837 * t784 + t785 * t868;
t725 = -0.2e1 * qJD(6) * t788 + (t787 * t826 - t746) * qJ(6) + (t788 * t826 + t745) * pkin(5) + t733;
t769 = -mrSges(7,2) * t787 + mrSges(7,3) * t826;
t772 = -mrSges(7,1) * t826 + mrSges(7,2) * t788;
t719 = m(7) * t725 + t745 * mrSges(7,1) - t746 * mrSges(7,3) + t787 * t769 - t788 * t772;
t789 = t836 * t814 - t835 * t816;
t776 = -qJDD(1) * pkin(2) - t845 * pkin(7) - t789;
t829 = t839 * t860;
t819 = qJDD(1) * t842 - t829;
t755 = (-t818 - t858) * pkin(8) + (-t819 + t829) * pkin(3) + t776;
t768 = t842 * t777 + t839 * t834;
t761 = -pkin(3) * t844 + qJDD(3) * pkin(8) + t817 * t861 + t768;
t734 = t841 * t755 - t838 * t761;
t811 = qJDD(4) - t819;
t730 = (t812 * t827 - t785) * pkin(9) + (t812 * t813 + t811) * pkin(4) + t734;
t735 = t838 * t755 + t841 * t761;
t732 = -pkin(4) * t810 + pkin(9) * t784 - t794 * t827 + t735;
t728 = t837 * t730 + t732 * t868;
t762 = pkin(5) * t787 - qJ(6) * t788;
t807 = qJDD(5) + t811;
t825 = t826 ^ 2;
t722 = -pkin(5) * t825 + qJ(6) * t807 + 0.2e1 * qJD(6) * t826 - t762 * t787 + t728;
t864 = t873 * t787 - t874 * t788 + t872 * t826;
t704 = -mrSges(6,1) * t733 - mrSges(7,1) * t725 + mrSges(7,2) * t722 + mrSges(6,3) * t728 - pkin(5) * t719 - t875 * t745 + t866 * t746 + t864 * t788 + t873 * t807 + t870 * t826;
t727 = t730 * t868 - t837 * t732;
t723 = -t807 * pkin(5) - t825 * qJ(6) + t788 * t762 + qJDD(6) - t727;
t705 = mrSges(6,2) * t733 + mrSges(7,2) * t723 - mrSges(6,3) * t727 - mrSges(7,3) * t725 - qJ(6) * t719 - t866 * t745 + t876 * t746 + t864 * t787 + t874 * t807 + t871 * t826;
t778 = Ifges(5,5) * t813 + Ifges(5,6) * t812 + Ifges(5,3) * t827;
t780 = Ifges(5,1) * t813 + Ifges(5,4) * t812 + Ifges(5,5) * t827;
t770 = -mrSges(6,2) * t826 - mrSges(6,3) * t787;
t771 = mrSges(6,1) * t826 - mrSges(6,3) * t788;
t850 = m(6) * t733 + t745 * mrSges(6,1) + t746 * mrSges(6,2) + t787 * t770 + t788 * t771 + t719;
t859 = m(7) * t722 + t807 * mrSges(7,3) + t826 * t772;
t763 = mrSges(7,1) * t787 - mrSges(7,3) * t788;
t863 = -mrSges(6,1) * t787 - mrSges(6,2) * t788 - t763;
t867 = -mrSges(6,3) - mrSges(7,2);
t711 = m(6) * t728 - t807 * mrSges(6,2) + t745 * t867 - t826 * t771 + t787 * t863 + t859;
t853 = -m(7) * t723 + t807 * mrSges(7,1) + t826 * t769;
t713 = m(6) * t727 + t807 * mrSges(6,1) + t746 * t867 + t826 * t770 + t788 * t863 + t853;
t854 = t711 * t868 - t713 * t837;
t685 = -mrSges(5,1) * t760 + mrSges(5,3) * t735 + Ifges(5,4) * t785 + Ifges(5,2) * t784 + Ifges(5,6) * t811 - pkin(4) * t850 + pkin(9) * t854 + t704 * t868 + t837 * t705 - t813 * t778 + t827 * t780;
t706 = t837 * t711 + t713 * t868;
t779 = Ifges(5,4) * t813 + Ifges(5,2) * t812 + Ifges(5,6) * t827;
t686 = mrSges(5,2) * t760 - mrSges(5,3) * t734 + Ifges(5,1) * t785 + Ifges(5,4) * t784 + Ifges(5,5) * t811 - pkin(9) * t706 - t837 * t704 + t705 * t868 + t812 * t778 - t827 * t779;
t791 = -mrSges(5,1) * t812 + mrSges(5,2) * t813;
t792 = -mrSges(5,2) * t827 + mrSges(5,3) * t812;
t702 = m(5) * t734 + mrSges(5,1) * t811 - mrSges(5,3) * t785 - t791 * t813 + t792 * t827 + t706;
t793 = mrSges(5,1) * t827 - mrSges(5,3) * t813;
t703 = m(5) * t735 - mrSges(5,2) * t811 + mrSges(5,3) * t784 + t791 * t812 - t793 * t827 + t854;
t700 = -t702 * t838 + t841 * t703;
t714 = -m(5) * t760 + t784 * mrSges(5,1) - t785 * mrSges(5,2) + t812 * t792 - t813 * t793 - t850;
t805 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t839 + Ifges(4,2) * t842) * qJD(1);
t806 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t839 + Ifges(4,4) * t842) * qJD(1);
t869 = mrSges(4,1) * t767 - mrSges(4,2) * t768 + Ifges(4,5) * t818 + Ifges(4,6) * t819 + Ifges(4,3) * qJDD(3) + pkin(3) * t714 + pkin(8) * t700 + t841 * t685 + t838 * t686 + (t805 * t839 - t806 * t842) * qJD(1);
t815 = (-mrSges(4,1) * t842 + mrSges(4,2) * t839) * qJD(1);
t821 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t862;
t698 = m(4) * t768 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t819 - qJD(3) * t821 + t815 * t861 + t700;
t822 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t861;
t708 = m(4) * t767 + qJDD(3) * mrSges(4,1) - t818 * mrSges(4,3) + qJD(3) * t822 - t815 * t862 + t714;
t855 = t842 * t698 - t708 * t839;
t689 = m(3) * t790 - mrSges(3,1) * t845 - qJDD(1) * mrSges(3,2) + t855;
t699 = t702 * t841 + t703 * t838;
t849 = -m(4) * t776 + t819 * mrSges(4,1) - mrSges(4,2) * t818 - t821 * t862 + t822 * t861 - t699;
t694 = m(3) * t789 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t845 + t849;
t682 = t835 * t689 + t836 * t694;
t679 = m(2) * t823 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t845 + t682;
t856 = t836 * t689 - t694 * t835;
t680 = m(2) * t824 - mrSges(2,1) * t845 - qJDD(1) * mrSges(2,2) + t856;
t865 = t843 * t679 + t840 * t680;
t692 = t839 * t698 + t842 * t708;
t690 = m(3) * t834 + t692;
t857 = -t679 * t840 + t843 * t680;
t804 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t839 + Ifges(4,6) * t842) * qJD(1);
t675 = mrSges(4,2) * t776 - mrSges(4,3) * t767 + Ifges(4,1) * t818 + Ifges(4,4) * t819 + Ifges(4,5) * qJDD(3) - pkin(8) * t699 - qJD(3) * t805 - t685 * t838 + t686 * t841 + t804 * t861;
t718 = t746 * mrSges(7,2) + t788 * t763 - t853;
t848 = -mrSges(6,1) * t727 + mrSges(7,1) * t723 + mrSges(6,2) * t728 - mrSges(7,3) * t722 + pkin(5) * t718 - qJ(6) * t859 + t872 * t807 + t871 * t788 + (qJ(6) * t763 - t870) * t787 - t874 * t746 + (qJ(6) * mrSges(7,2) + t873) * t745;
t846 = mrSges(5,1) * t734 - mrSges(5,2) * t735 + Ifges(5,5) * t785 + Ifges(5,6) * t784 + Ifges(5,3) * t811 + pkin(4) * t706 + t813 * t779 - t812 * t780 - t848;
t684 = -mrSges(4,1) * t776 + mrSges(4,3) * t768 + Ifges(4,4) * t818 + Ifges(4,2) * t819 + Ifges(4,6) * qJDD(3) - pkin(3) * t699 + qJD(3) * t806 - t804 * t862 - t846;
t851 = mrSges(2,1) * t823 + mrSges(3,1) * t789 - mrSges(2,2) * t824 - mrSges(3,2) * t790 + pkin(1) * t682 + pkin(2) * t849 + pkin(7) * t855 + t839 * t675 + t842 * t684 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t673 = -mrSges(3,1) * t834 + mrSges(3,3) * t790 + t845 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t692 - t869;
t672 = mrSges(3,2) * t834 - mrSges(3,3) * t789 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t845 - pkin(7) * t692 + t675 * t842 - t684 * t839;
t671 = -mrSges(2,2) * g(3) - mrSges(2,3) * t823 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t845 - qJ(2) * t682 + t672 * t836 - t673 * t835;
t670 = mrSges(2,1) * g(3) + mrSges(2,3) * t824 + t845 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t690 + qJ(2) * t856 + t835 * t672 + t836 * t673;
t1 = [-m(1) * g(1) + t857; -m(1) * g(2) + t865; (-m(1) - m(2)) * g(3) + t690; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t865 - t840 * t670 + t843 * t671; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t857 + t843 * t670 + t840 * t671; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t851; t851; t690; t869; t846; -t848; t718;];
tauJB  = t1;
