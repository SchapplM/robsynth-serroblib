% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:33:52
% EndTime: 2019-05-06 08:34:00
% DurationCPUTime: 4.29s
% Computational Cost: add. (40418->347), mult. (88170->413), div. (0->0), fcn. (46942->8), ass. (0->136)
t852 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t830 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t829 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t828 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t851 = -Ifges(3,3) - Ifges(4,2) - Ifges(5,3);
t850 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t803 = cos(qJ(2));
t800 = sin(qJ(2));
t833 = qJD(1) * qJD(2);
t823 = t800 * t833;
t764 = qJDD(1) * t803 - t823;
t834 = t800 * qJD(1);
t769 = -qJD(2) * pkin(3) - qJ(4) * t834;
t832 = qJD(1) * qJD(4);
t806 = qJD(1) ^ 2;
t840 = t803 ^ 2 * t806;
t801 = sin(qJ(1));
t804 = cos(qJ(1));
t777 = -g(1) * t804 - g(2) * t801;
t739 = -pkin(1) * t806 + qJDD(1) * pkin(7) + t777;
t714 = -t800 * g(3) + t739 * t803;
t758 = (-pkin(2) * t803 - qJ(3) * t800) * qJD(1);
t836 = qJD(1) * t803;
t848 = qJDD(2) * qJ(3) + t758 * t836 + t714;
t849 = pkin(3) * t840 + qJ(4) * t764 - qJD(2) * t769 + 0.2e1 * t803 * t832 - t848;
t822 = t803 * t833;
t763 = qJDD(1) * t800 + t822;
t847 = -0.2e1 * t800 * t832 + (-t763 + t822) * qJ(4);
t846 = 2 * qJD(5);
t805 = qJD(2) ^ 2;
t845 = t805 * pkin(2);
t844 = mrSges(3,3) + mrSges(4,2);
t843 = -pkin(2) - qJ(5);
t842 = pkin(3) + qJ(5);
t839 = t803 * t806;
t713 = -g(3) * t803 - t739 * t800;
t759 = (-mrSges(4,1) * t803 - mrSges(4,3) * t800) * qJD(1);
t760 = (-mrSges(3,1) * t803 + mrSges(3,2) * t800) * qJD(1);
t773 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t836;
t774 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t836;
t818 = t758 * t834 + qJDD(3) - t713;
t703 = -qJDD(2) * pkin(2) - t805 * qJ(3) + t818;
t775 = mrSges(4,2) * t836 + qJD(2) * mrSges(4,3);
t692 = (-t800 * t839 - qJDD(2)) * pkin(3) + t703 + t847;
t762 = (mrSges(5,1) * t800 - mrSges(5,2) * t803) * qJD(1);
t776 = g(1) * t801 - g(2) * t804;
t738 = -qJDD(1) * pkin(1) - pkin(7) * t806 - t776;
t815 = -t764 * pkin(2) + t738 + (-t763 - t822) * qJ(3);
t811 = -qJ(4) * t840 + qJDD(4) - t815 + ((2 * qJD(3)) + t769) * t834;
t683 = t811 + t842 * t764 + t763 * pkin(4) + (pkin(4) * t803 + t800 * t843) * t833;
t761 = (pkin(4) * t800 + qJ(5) * t803) * qJD(1);
t687 = (-pkin(4) - qJ(3)) * t805 + (-pkin(3) * t839 - qJD(1) * t761) * t800 + (-pkin(2) - t842) * qJDD(2) + t818 + t847;
t796 = sin(pkin(9));
t797 = cos(pkin(9));
t750 = qJD(2) * t796 + t797 * t836;
t677 = t683 * t797 - t796 * t687 + t750 * t846;
t719 = -qJDD(2) * t796 - t764 * t797;
t749 = -qJD(2) * t797 + t796 * t836;
t675 = (t749 * t834 - t719) * pkin(8) + (-t749 * t750 + t763) * pkin(5) + t677;
t678 = t683 * t796 + t687 * t797 + t749 * t846;
t718 = -qJDD(2) * t797 + t764 * t796;
t720 = pkin(5) * t834 + pkin(8) * t750;
t746 = t749 ^ 2;
t676 = -pkin(5) * t746 + pkin(8) * t718 - t720 * t834 + t678;
t799 = sin(qJ(6));
t802 = cos(qJ(6));
t673 = t675 * t802 - t676 * t799;
t710 = t749 * t802 + t750 * t799;
t694 = qJD(6) * t710 + t718 * t799 + t719 * t802;
t711 = t749 * t799 - t750 * t802;
t699 = -mrSges(7,1) * t710 + mrSges(7,2) * t711;
t780 = qJD(6) + t834;
t704 = -mrSges(7,2) * t780 + mrSges(7,3) * t710;
t756 = qJDD(6) + t763;
t671 = m(7) * t673 + mrSges(7,1) * t756 - mrSges(7,3) * t694 - t699 * t711 + t704 * t780;
t674 = t675 * t799 + t676 * t802;
t693 = -qJD(6) * t711 + t718 * t802 - t719 * t799;
t705 = mrSges(7,1) * t780 - mrSges(7,3) * t711;
t672 = m(7) * t674 - mrSges(7,2) * t756 + mrSges(7,3) * t693 + t699 * t710 - t705 * t780;
t662 = t671 * t802 + t672 * t799;
t712 = -mrSges(6,1) * t749 - mrSges(6,2) * t750;
t716 = -mrSges(6,2) * t834 + mrSges(6,3) * t749;
t660 = m(6) * t677 + mrSges(6,1) * t763 - mrSges(6,3) * t719 + t712 * t750 + t716 * t834 + t662;
t717 = mrSges(6,1) * t834 + mrSges(6,3) * t750;
t819 = -t671 * t799 + t672 * t802;
t661 = m(6) * t678 - mrSges(6,2) * t763 + mrSges(6,3) * t718 + t712 * t749 - t717 * t834 + t819;
t837 = -t660 * t796 + t661 * t797;
t817 = -m(5) * t692 + t762 * t834 - t837;
t813 = -m(4) * t703 + qJDD(2) * mrSges(4,1) + qJD(2) * t775 + t817;
t652 = m(3) * t713 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t773 + t774) * qJD(2) + (-t759 - t760) * t834 + (mrSges(5,3) - t844) * t763 + t813;
t771 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t834;
t831 = qJD(3) * qJD(2);
t785 = 0.2e1 * t831;
t702 = t785 - t845 + t848;
t772 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t834;
t691 = -0.2e1 * t831 + t845 + t849;
t770 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t834;
t686 = qJDD(2) * pkin(4) - t761 * t836 + t805 * t843 + qJDD(5) + t785 - t849;
t680 = -t718 * pkin(5) - t746 * pkin(8) - t750 * t720 + t686;
t814 = m(7) * t680 - t693 * mrSges(7,1) + mrSges(7,2) * t694 - t710 * t704 + t705 * t711;
t810 = -m(6) * t686 + t718 * mrSges(6,1) - mrSges(6,2) * t719 + t749 * t716 + t717 * t750 - t814;
t809 = -m(5) * t691 + qJDD(2) * mrSges(5,1) - mrSges(5,3) * t764 + qJD(2) * t770 - t810;
t807 = m(4) * t702 + qJDD(2) * mrSges(4,3) + qJD(2) * t772 + t759 * t836 + t809;
t667 = t807 + t844 * t764 - qJD(2) * t771 + m(3) * t714 - qJDD(2) * mrSges(3,2) + (t760 - t762) * t836;
t820 = -t652 * t800 + t667 * t803;
t645 = m(2) * t777 - mrSges(2,1) * t806 - qJDD(1) * mrSges(2,2) + t820;
t700 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t834 + t815;
t655 = t660 * t797 + t661 * t796;
t689 = -pkin(2) * t823 + t764 * pkin(3) + t811;
t816 = -m(5) * t689 - mrSges(5,1) * t763 + mrSges(5,2) * t764 - t770 * t834 + t773 * t836 - t655;
t653 = m(4) * t700 - mrSges(4,1) * t764 - mrSges(4,3) * t763 - t772 * t834 - t775 * t836 + t816;
t808 = -m(3) * t738 + mrSges(3,1) * t764 - mrSges(3,2) * t763 - t771 * t834 + t774 * t836 - t653;
t650 = m(2) * t776 + qJDD(1) * mrSges(2,1) - t806 * mrSges(2,2) + t808;
t838 = t645 * t801 + t650 * t804;
t646 = t652 * t803 + t800 * t667;
t835 = qJD(2) * t773;
t827 = t851 * qJD(2) + (-t829 * t800 - t828 * t803) * qJD(1);
t826 = -t828 * qJD(2) + (-t800 * t830 - t803 * t850) * qJD(1);
t825 = t829 * qJD(2) + (t800 * t852 + t803 * t830) * qJD(1);
t821 = t645 * t804 - t650 * t801;
t708 = -Ifges(6,1) * t750 + Ifges(6,4) * t749 + Ifges(6,5) * t834;
t707 = -Ifges(6,4) * t750 + Ifges(6,2) * t749 + Ifges(6,6) * t834;
t706 = -Ifges(6,5) * t750 + Ifges(6,6) * t749 + Ifges(6,3) * t834;
t697 = Ifges(7,1) * t711 + Ifges(7,4) * t710 + Ifges(7,5) * t780;
t696 = Ifges(7,4) * t711 + Ifges(7,2) * t710 + Ifges(7,6) * t780;
t695 = Ifges(7,5) * t711 + Ifges(7,6) * t710 + Ifges(7,3) * t780;
t664 = mrSges(7,2) * t680 - mrSges(7,3) * t673 + Ifges(7,1) * t694 + Ifges(7,4) * t693 + Ifges(7,5) * t756 + t695 * t710 - t696 * t780;
t663 = -mrSges(7,1) * t680 + mrSges(7,3) * t674 + Ifges(7,4) * t694 + Ifges(7,2) * t693 + Ifges(7,6) * t756 - t695 * t711 + t697 * t780;
t654 = qJDD(2) * mrSges(5,2) - t763 * mrSges(5,3) - t817 + t835;
t648 = mrSges(6,2) * t686 - mrSges(6,3) * t677 + Ifges(6,1) * t719 + Ifges(6,4) * t718 + Ifges(6,5) * t763 - pkin(8) * t662 - t663 * t799 + t664 * t802 + t706 * t749 - t707 * t834;
t647 = -mrSges(6,1) * t686 + mrSges(6,3) * t678 + Ifges(6,4) * t719 + Ifges(6,2) * t718 + Ifges(6,6) * t763 - pkin(5) * t814 + pkin(8) * t819 + t802 * t663 + t799 * t664 + t750 * t706 + t708 * t834;
t642 = -t827 * t836 + t829 * qJDD(2) + t830 * t764 + t826 * qJD(2) + (Ifges(6,3) + t852) * t763 + Ifges(7,3) * t756 - t749 * t708 - t750 * t707 + Ifges(6,6) * t718 + Ifges(6,5) * t719 + mrSges(3,2) * t738 - t710 * t697 + t711 * t696 - mrSges(3,3) * t713 - mrSges(4,3) * t700 + mrSges(4,2) * t703 + mrSges(5,1) * t689 - mrSges(5,3) * t692 + Ifges(7,6) * t693 + Ifges(7,5) * t694 + mrSges(6,1) * t677 - mrSges(6,2) * t678 + mrSges(7,1) * t673 - mrSges(7,2) * t674 + pkin(5) * t662 + pkin(4) * t655 - qJ(4) * t654 - qJ(3) * t653;
t641 = -qJ(4) * t809 + t796 * t647 - t797 * t648 - pkin(3) * t816 - mrSges(3,1) * t738 + mrSges(3,3) * t714 - mrSges(4,1) * t700 + mrSges(4,2) * t702 - mrSges(5,2) * t689 + mrSges(5,3) * t691 + qJ(5) * t655 - pkin(2) * t653 + t850 * t764 + t830 * t763 + t828 * qJDD(2) + t825 * qJD(2) + (qJ(4) * t803 * t762 + t800 * t827) * qJD(1);
t640 = -pkin(2) * (t813 - t835) - qJ(3) * t807 + Ifges(2,6) * qJDD(1) + t806 * Ifges(2,5) + t796 * t648 + t797 * t647 + mrSges(2,3) * t777 + pkin(4) * t810 - mrSges(3,1) * t713 + mrSges(3,2) * t714 - mrSges(4,3) * t702 + mrSges(4,1) * t703 + mrSges(5,1) * t691 - mrSges(5,2) * t692 + qJ(5) * t837 + pkin(3) * t654 - pkin(1) * t646 + mrSges(2,1) * g(3) + (-mrSges(4,2) * qJ(3) - t828) * t764 + (mrSges(5,2) * pkin(2) + t851) * qJDD(2) + (-pkin(2) * (-mrSges(4,2) + mrSges(5,3)) - t829) * t763 + ((qJ(3) * t762 + t825) * t803 + (pkin(2) * t759 + t826) * t800) * qJD(1);
t639 = -mrSges(2,2) * g(3) - mrSges(2,3) * t776 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t806 - pkin(7) * t646 - t641 * t800 + t642 * t803;
t1 = [-m(1) * g(1) + t821; -m(1) * g(2) + t838; (-m(1) - m(2)) * g(3) + t646; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t838 + t639 * t804 - t640 * t801; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t821 + t801 * t639 + t804 * t640; -mrSges(1,1) * g(2) + mrSges(2,1) * t776 + mrSges(1,2) * g(1) - mrSges(2,2) * t777 + Ifges(2,3) * qJDD(1) + pkin(1) * t808 + pkin(7) * t820 + t803 * t641 + t800 * t642;];
tauB  = t1;
