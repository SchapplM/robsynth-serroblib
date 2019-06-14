% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:09:01
% EndTime: 2019-05-05 23:09:13
% DurationCPUTime: 10.58s
% Computational Cost: add. (178063->338), mult. (361467->412), div. (0->0), fcn. (240902->10), ass. (0->136)
t843 = sin(qJ(1));
t847 = cos(qJ(1));
t823 = -t847 * g(1) - t843 * g(2);
t877 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t823;
t849 = qJD(1) ^ 2;
t875 = (-pkin(1) - pkin(7));
t793 = (t849 * t875) - t877;
t846 = cos(qJ(3));
t868 = qJD(1) * qJD(3);
t826 = t846 * t868;
t842 = sin(qJ(3));
t817 = -t842 * qJDD(1) - t826;
t866 = t842 * t868;
t818 = qJDD(1) * t846 - t866;
t764 = (-t818 + t866) * pkin(8) + (-t817 + t826) * pkin(3) + t793;
t822 = t843 * g(1) - t847 * g(2);
t858 = -t849 * qJ(2) + qJDD(2) - t822;
t794 = qJDD(1) * t875 + t858;
t787 = -g(3) * t846 + t842 * t794;
t816 = (pkin(3) * t842 - pkin(8) * t846) * qJD(1);
t828 = t842 * qJD(1);
t848 = qJD(3) ^ 2;
t767 = -pkin(3) * t848 + qJDD(3) * pkin(8) - t816 * t828 + t787;
t841 = sin(qJ(4));
t845 = cos(qJ(4));
t748 = t845 * t764 - t841 * t767;
t869 = qJD(1) * t846;
t813 = qJD(3) * t845 - t841 * t869;
t780 = qJD(4) * t813 + qJDD(3) * t841 + t818 * t845;
t812 = qJDD(4) - t817;
t814 = qJD(3) * t841 + t845 * t869;
t825 = t828 + qJD(4);
t738 = (t813 * t825 - t780) * qJ(5) + (t813 * t814 + t812) * pkin(4) + t748;
t749 = t841 * t764 + t845 * t767;
t779 = -qJD(4) * t814 + qJDD(3) * t845 - t818 * t841;
t789 = pkin(4) * t825 - qJ(5) * t814;
t811 = t813 ^ 2;
t740 = -pkin(4) * t811 + qJ(5) * t779 - t789 * t825 + t749;
t838 = sin(pkin(10));
t839 = cos(pkin(10));
t783 = t813 * t838 + t814 * t839;
t725 = -0.2e1 * qJD(5) * t783 + t839 * t738 - t838 * t740;
t757 = t779 * t838 + t780 * t839;
t782 = t813 * t839 - t814 * t838;
t723 = (t782 * t825 - t757) * pkin(9) + (t782 * t783 + t812) * pkin(5) + t725;
t726 = 0.2e1 * qJD(5) * t782 + t838 * t738 + t839 * t740;
t756 = t779 * t839 - t780 * t838;
t770 = pkin(5) * t825 - pkin(9) * t783;
t781 = t782 ^ 2;
t724 = -pkin(5) * t781 + pkin(9) * t756 - t770 * t825 + t726;
t840 = sin(qJ(6));
t844 = cos(qJ(6));
t721 = t723 * t844 - t724 * t840;
t759 = t782 * t844 - t783 * t840;
t734 = qJD(6) * t759 + t756 * t840 + t757 * t844;
t760 = t782 * t840 + t783 * t844;
t746 = -mrSges(7,1) * t759 + mrSges(7,2) * t760;
t824 = qJD(6) + t825;
t750 = -mrSges(7,2) * t824 + mrSges(7,3) * t759;
t805 = qJDD(6) + t812;
t715 = m(7) * t721 + mrSges(7,1) * t805 - mrSges(7,3) * t734 - t746 * t760 + t750 * t824;
t722 = t723 * t840 + t724 * t844;
t733 = -qJD(6) * t760 + t756 * t844 - t757 * t840;
t751 = mrSges(7,1) * t824 - mrSges(7,3) * t760;
t716 = m(7) * t722 - mrSges(7,2) * t805 + mrSges(7,3) * t733 + t746 * t759 - t751 * t824;
t709 = t844 * t715 + t840 * t716;
t761 = -mrSges(6,1) * t782 + mrSges(6,2) * t783;
t768 = -mrSges(6,2) * t825 + mrSges(6,3) * t782;
t707 = m(6) * t725 + mrSges(6,1) * t812 - mrSges(6,3) * t757 - t761 * t783 + t768 * t825 + t709;
t769 = mrSges(6,1) * t825 - mrSges(6,3) * t783;
t861 = -t715 * t840 + t844 * t716;
t708 = m(6) * t726 - mrSges(6,2) * t812 + mrSges(6,3) * t756 + t761 * t782 - t769 * t825 + t861;
t703 = t839 * t707 + t838 * t708;
t754 = Ifges(6,4) * t783 + Ifges(6,2) * t782 + Ifges(6,6) * t825;
t755 = Ifges(6,1) * t783 + Ifges(6,4) * t782 + Ifges(6,5) * t825;
t772 = Ifges(5,4) * t814 + Ifges(5,2) * t813 + Ifges(5,6) * t825;
t773 = Ifges(5,1) * t814 + Ifges(5,4) * t813 + Ifges(5,5) * t825;
t742 = Ifges(7,4) * t760 + Ifges(7,2) * t759 + Ifges(7,6) * t824;
t743 = Ifges(7,1) * t760 + Ifges(7,4) * t759 + Ifges(7,5) * t824;
t855 = -mrSges(7,1) * t721 + mrSges(7,2) * t722 - Ifges(7,5) * t734 - Ifges(7,6) * t733 - Ifges(7,3) * t805 - t760 * t742 + t759 * t743;
t876 = mrSges(5,1) * t748 + mrSges(6,1) * t725 - mrSges(5,2) * t749 - mrSges(6,2) * t726 + Ifges(5,5) * t780 + Ifges(6,5) * t757 + Ifges(5,6) * t779 + Ifges(6,6) * t756 + pkin(4) * t703 + pkin(5) * t709 + t783 * t754 - t782 * t755 + t814 * t772 - t813 * t773 + (Ifges(5,3) + Ifges(6,3)) * t812 - t855;
t874 = mrSges(2,1) - mrSges(3,2);
t873 = -Ifges(3,4) + Ifges(2,5);
t872 = (Ifges(3,5) - Ifges(2,6));
t815 = (mrSges(4,1) * t842 + mrSges(4,2) * t846) * qJD(1);
t821 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t869;
t785 = -mrSges(5,1) * t813 + mrSges(5,2) * t814;
t788 = -mrSges(5,2) * t825 + mrSges(5,3) * t813;
t701 = m(5) * t748 + mrSges(5,1) * t812 - mrSges(5,3) * t780 - t785 * t814 + t788 * t825 + t703;
t790 = mrSges(5,1) * t825 - mrSges(5,3) * t814;
t862 = -t707 * t838 + t839 * t708;
t702 = m(5) * t749 - mrSges(5,2) * t812 + mrSges(5,3) * t779 + t785 * t813 - t790 * t825 + t862;
t863 = -t701 * t841 + t845 * t702;
t693 = m(4) * t787 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t817 - qJD(3) * t821 - t815 * t828 + t863;
t786 = g(3) * t842 + t794 * t846;
t820 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t828;
t766 = -qJDD(3) * pkin(3) - pkin(8) * t848 + t816 * t869 - t786;
t747 = -pkin(4) * t779 - qJ(5) * t811 + t814 * t789 + qJDD(5) + t766;
t728 = -pkin(5) * t756 - pkin(9) * t781 + t770 * t783 + t747;
t860 = m(7) * t728 - t733 * mrSges(7,1) + t734 * mrSges(7,2) - t759 * t750 + t760 * t751;
t719 = m(6) * t747 - t756 * mrSges(6,1) + t757 * mrSges(6,2) - t782 * t768 + t783 * t769 + t860;
t851 = -m(5) * t766 + t779 * mrSges(5,1) - t780 * mrSges(5,2) + t813 * t788 - t814 * t790 - t719;
t717 = m(4) * t786 + qJDD(3) * mrSges(4,1) - t818 * mrSges(4,3) + qJD(3) * t820 - t815 * t869 + t851;
t687 = t842 * t693 + t846 * t717;
t799 = -qJDD(1) * pkin(1) + t858;
t857 = -m(3) * t799 + (t849 * mrSges(3,3)) - t687;
t683 = m(2) * t822 - (t849 * mrSges(2,2)) + qJDD(1) * t874 + t857;
t797 = t849 * pkin(1) + t877;
t695 = t845 * t701 + t841 * t702;
t856 = -m(4) * t793 + mrSges(4,1) * t817 - t818 * mrSges(4,2) - t820 * t828 - t821 * t869 - t695;
t853 = -m(3) * t797 + (t849 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t856;
t690 = m(2) * t823 - (mrSges(2,1) * t849) - qJDD(1) * mrSges(2,2) + t853;
t870 = t847 * t683 + t843 * t690;
t865 = -t683 * t843 + t847 * t690;
t864 = t846 * t693 - t842 * t717;
t741 = Ifges(7,5) * t760 + Ifges(7,6) * t759 + Ifges(7,3) * t824;
t710 = -mrSges(7,1) * t728 + mrSges(7,3) * t722 + Ifges(7,4) * t734 + Ifges(7,2) * t733 + Ifges(7,6) * t805 - t741 * t760 + t743 * t824;
t711 = mrSges(7,2) * t728 - mrSges(7,3) * t721 + Ifges(7,1) * t734 + Ifges(7,4) * t733 + Ifges(7,5) * t805 + t741 * t759 - t742 * t824;
t753 = Ifges(6,5) * t783 + Ifges(6,6) * t782 + Ifges(6,3) * t825;
t696 = -mrSges(6,1) * t747 + mrSges(6,3) * t726 + Ifges(6,4) * t757 + Ifges(6,2) * t756 + Ifges(6,6) * t812 - pkin(5) * t860 + pkin(9) * t861 + t844 * t710 + t840 * t711 - t783 * t753 + t825 * t755;
t697 = mrSges(6,2) * t747 - mrSges(6,3) * t725 + Ifges(6,1) * t757 + Ifges(6,4) * t756 + Ifges(6,5) * t812 - pkin(9) * t709 - t710 * t840 + t711 * t844 + t753 * t782 - t754 * t825;
t771 = Ifges(5,5) * t814 + Ifges(5,6) * t813 + Ifges(5,3) * t825;
t679 = -mrSges(5,1) * t766 + mrSges(5,3) * t749 + Ifges(5,4) * t780 + Ifges(5,2) * t779 + Ifges(5,6) * t812 - pkin(4) * t719 + qJ(5) * t862 + t839 * t696 + t838 * t697 - t814 * t771 + t825 * t773;
t681 = mrSges(5,2) * t766 - mrSges(5,3) * t748 + Ifges(5,1) * t780 + Ifges(5,4) * t779 + Ifges(5,5) * t812 - qJ(5) * t703 - t696 * t838 + t697 * t839 + t771 * t813 - t772 * t825;
t803 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t846 - Ifges(4,2) * t842) * qJD(1);
t804 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t846 - Ifges(4,4) * t842) * qJD(1);
t854 = mrSges(4,1) * t786 - mrSges(4,2) * t787 + Ifges(4,5) * t818 + Ifges(4,6) * t817 + Ifges(4,3) * qJDD(3) + pkin(3) * t851 + pkin(8) * t863 + t845 * t679 + t841 * t681 + t803 * t869 + t804 * t828;
t802 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t846 - Ifges(4,6) * t842) * qJD(1);
t676 = mrSges(4,2) * t793 - mrSges(4,3) * t786 + Ifges(4,1) * t818 + Ifges(4,4) * t817 + Ifges(4,5) * qJDD(3) - pkin(8) * t695 - qJD(3) * t803 - t679 * t841 + t681 * t845 - t802 * t828;
t677 = -mrSges(4,1) * t793 + mrSges(4,3) * t787 + Ifges(4,4) * t818 + Ifges(4,2) * t817 + Ifges(4,6) * qJDD(3) - pkin(3) * t695 + qJD(3) * t804 - t802 * t869 - t876;
t685 = qJDD(1) * mrSges(3,2) - t857;
t852 = mrSges(2,1) * t822 - mrSges(2,2) * t823 + mrSges(3,2) * t799 - mrSges(3,3) * t797 - pkin(1) * t685 - pkin(7) * t687 + qJ(2) * t853 + t846 * t676 - t677 * t842 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t686 = -m(3) * g(3) + t864;
t674 = t854 + (t872 * t849) + pkin(2) * t687 - qJ(2) * t686 - mrSges(2,3) * t822 + mrSges(3,1) * t799 + t873 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t673 = -mrSges(3,1) * t797 + mrSges(2,3) * t823 - pkin(1) * t686 - pkin(2) * t856 - pkin(7) * t864 + g(3) * t874 - qJDD(1) * t872 - t842 * t676 - t846 * t677 + t849 * t873;
t1 = [-m(1) * g(1) + t865; -m(1) * g(2) + t870; (-m(1) - m(2) - m(3)) * g(3) + t864; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t870 - t843 * t673 + t847 * t674; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t865 + t847 * t673 + t843 * t674; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t852; t852; t685; t854; t876; t719; -t855;];
tauJB  = t1;
