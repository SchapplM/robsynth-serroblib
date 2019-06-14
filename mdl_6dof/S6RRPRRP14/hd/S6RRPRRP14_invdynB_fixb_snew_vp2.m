% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:18:11
% EndTime: 2019-05-06 19:18:31
% DurationCPUTime: 9.28s
% Computational Cost: add. (127492->359), mult. (283873->436), div. (0->0), fcn. (204130->10), ass. (0->151)
t838 = -2 * qJD(3);
t837 = Ifges(3,1) + Ifges(4,2);
t836 = Ifges(6,1) + Ifges(7,1);
t824 = Ifges(3,4) + Ifges(4,6);
t823 = Ifges(6,4) - Ifges(7,5);
t821 = Ifges(3,5) - Ifges(4,4);
t835 = -Ifges(6,5) - Ifges(7,4);
t834 = Ifges(3,2) + Ifges(4,3);
t833 = Ifges(6,2) + Ifges(7,3);
t820 = Ifges(3,6) - Ifges(4,5);
t819 = Ifges(6,6) - Ifges(7,6);
t832 = Ifges(3,3) + Ifges(4,1);
t831 = -Ifges(6,3) - Ifges(7,2);
t775 = cos(pkin(6));
t769 = qJD(1) * t775 + qJD(2);
t778 = sin(qJ(2));
t774 = sin(pkin(6));
t804 = qJD(1) * t774;
t797 = t778 * t804;
t830 = (pkin(2) * t769 + t838) * t797;
t779 = sin(qJ(1));
t782 = cos(qJ(1));
t764 = t779 * g(1) - g(2) * t782;
t783 = qJD(1) ^ 2;
t748 = pkin(8) * t774 * t783 + qJDD(1) * pkin(1) + t764;
t765 = -g(1) * t782 - g(2) * t779;
t801 = qJDD(1) * t774;
t749 = -pkin(1) * t783 + pkin(8) * t801 + t765;
t781 = cos(qJ(2));
t815 = t775 * t778;
t817 = t774 * t778;
t713 = -g(3) * t817 + t748 * t815 + t781 * t749;
t750 = (-pkin(2) * t781 - qJ(3) * t778) * t804;
t767 = t769 ^ 2;
t768 = qJDD(1) * t775 + qJDD(2);
t803 = qJD(1) * t781;
t796 = t774 * t803;
t681 = t767 * pkin(2) - t768 * qJ(3) - t750 * t796 + t769 * t838 - t713;
t829 = -pkin(2) - pkin(9);
t828 = cos(qJ(5));
t827 = t775 * g(3);
t826 = mrSges(3,1) - mrSges(4,2);
t825 = -mrSges(6,3) - mrSges(7,2);
t818 = t774 ^ 2 * t783;
t816 = t774 * t781;
t814 = t775 * t781;
t726 = -t774 * t748 - t827;
t744 = mrSges(3,1) * t769 - mrSges(3,3) * t797;
t745 = -mrSges(3,2) * t769 + mrSges(3,3) * t796;
t747 = mrSges(4,1) * t797 + mrSges(4,2) * t769;
t754 = (qJD(2) * t803 + qJDD(1) * t778) * t774;
t755 = -qJD(2) * t797 + t781 * t801;
t682 = -t755 * pkin(2) + (-t769 * t796 - t754) * qJ(3) + t726 + t830;
t746 = -mrSges(4,1) * t796 - mrSges(4,3) * t769;
t753 = pkin(3) * t797 - pkin(9) * t769;
t800 = t781 ^ 2 * t818;
t673 = -pkin(3) * t800 - t827 - t754 * qJ(3) + t829 * t755 + (-t748 + (-qJ(3) * t769 * t781 - t753 * t778) * qJD(1)) * t774 + t830;
t805 = g(3) * t816 + t778 * t749;
t790 = -t767 * qJ(3) + t750 * t797 + qJDD(3) + t805;
t675 = t754 * pkin(3) + t829 * t768 + (-pkin(3) * t769 * t804 - pkin(9) * t778 * t818 - t748 * t775) * t781 + t790;
t777 = sin(qJ(4));
t780 = cos(qJ(4));
t668 = t780 * t673 + t777 * t675;
t738 = t769 * t780 - t777 * t796;
t710 = -qJD(4) * t738 - t755 * t780 - t768 * t777;
t737 = -t769 * t777 - t780 * t796;
t714 = -mrSges(5,1) * t737 + mrSges(5,2) * t738;
t760 = qJD(4) + t797;
t719 = mrSges(5,1) * t760 - mrSges(5,3) * t738;
t743 = qJDD(4) + t754;
t715 = -pkin(4) * t737 - pkin(10) * t738;
t758 = t760 ^ 2;
t664 = -pkin(4) * t758 + pkin(10) * t743 + t715 * t737 + t668;
t672 = t755 * pkin(3) - pkin(9) * t800 + t769 * t753 - t681;
t711 = qJD(4) * t737 - t755 * t777 + t768 * t780;
t666 = (-t737 * t760 - t711) * pkin(10) + (t738 * t760 - t710) * pkin(4) + t672;
t776 = sin(qJ(5));
t661 = t828 * t664 + t776 * t666;
t717 = t828 * t738 + t776 * t760;
t678 = qJD(5) * t717 + t711 * t776 - t828 * t743;
t735 = qJD(5) - t737;
t700 = mrSges(6,1) * t735 - mrSges(6,3) * t717;
t708 = qJDD(5) - t710;
t716 = t738 * t776 - t828 * t760;
t694 = pkin(5) * t716 - qJ(6) * t717;
t734 = t735 ^ 2;
t657 = -pkin(5) * t734 + qJ(6) * t708 + 0.2e1 * qJD(6) * t735 - t694 * t716 + t661;
t701 = -mrSges(7,1) * t735 + mrSges(7,2) * t717;
t798 = m(7) * t657 + t708 * mrSges(7,3) + t735 * t701;
t695 = mrSges(7,1) * t716 - mrSges(7,3) * t717;
t809 = -mrSges(6,1) * t716 - mrSges(6,2) * t717 - t695;
t652 = m(6) * t661 - t708 * mrSges(6,2) + t825 * t678 - t735 * t700 + t809 * t716 + t798;
t660 = -t776 * t664 + t828 * t666;
t679 = -t716 * qJD(5) + t828 * t711 + t776 * t743;
t699 = -mrSges(6,2) * t735 - mrSges(6,3) * t716;
t658 = -t708 * pkin(5) - t734 * qJ(6) + t717 * t694 + qJDD(6) - t660;
t698 = -mrSges(7,2) * t716 + mrSges(7,3) * t735;
t792 = -m(7) * t658 + t708 * mrSges(7,1) + t735 * t698;
t654 = m(6) * t660 + t708 * mrSges(6,1) + t825 * t679 + t735 * t699 + t809 * t717 + t792;
t793 = t828 * t652 - t654 * t776;
t645 = m(5) * t668 - mrSges(5,2) * t743 + mrSges(5,3) * t710 + t714 * t737 - t719 * t760 + t793;
t667 = -t777 * t673 + t780 * t675;
t718 = -mrSges(5,2) * t760 + mrSges(5,3) * t737;
t663 = -t743 * pkin(4) - t758 * pkin(10) + t738 * t715 - t667;
t659 = -0.2e1 * qJD(6) * t717 + (t716 * t735 - t679) * qJ(6) + (t717 * t735 + t678) * pkin(5) + t663;
t655 = m(7) * t659 + mrSges(7,1) * t678 - t679 * mrSges(7,3) + t698 * t716 - t717 * t701;
t784 = -m(6) * t663 - t678 * mrSges(6,1) - mrSges(6,2) * t679 - t716 * t699 - t700 * t717 - t655;
t649 = m(5) * t667 + mrSges(5,1) * t743 - mrSges(5,3) * t711 - t714 * t738 + t718 * t760 + t784;
t794 = t780 * t645 - t777 * t649;
t791 = m(4) * t682 - t754 * mrSges(4,3) + t746 * t796 + t794;
t636 = m(3) * t726 + t754 * mrSges(3,2) - t826 * t755 + (-t745 * t781 + (t744 - t747) * t778) * t804 + t791;
t799 = t748 * t814;
t712 = t799 - t805;
t751 = (mrSges(4,2) * t781 - mrSges(4,3) * t778) * t804;
t752 = (-mrSges(3,1) * t781 + mrSges(3,2) * t778) * t804;
t639 = t777 * t645 + t780 * t649;
t691 = -t768 * pkin(2) + t790 - t799;
t788 = -m(4) * t691 - t754 * mrSges(4,1) - t639;
t637 = m(3) * t712 - t754 * mrSges(3,3) + (t745 - t746) * t769 + t826 * t768 + (-t751 - t752) * t797 + t788;
t648 = t776 * t652 + t828 * t654;
t786 = -m(5) * t672 + mrSges(5,1) * t710 - t711 * mrSges(5,2) + t718 * t737 - t738 * t719 - t648;
t785 = -m(4) * t681 + t768 * mrSges(4,3) + t769 * t747 + t751 * t796 - t786;
t643 = m(3) * t713 - mrSges(3,2) * t768 - t744 * t769 + (mrSges(3,3) + mrSges(4,1)) * t755 + t785 + t752 * t796;
t626 = -t636 * t774 + t637 * t814 + t643 * t815;
t624 = m(2) * t764 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t783 + t626;
t631 = -t637 * t778 + t781 * t643;
t630 = m(2) * t765 - mrSges(2,1) * t783 - qJDD(1) * mrSges(2,2) + t631;
t813 = t782 * t624 + t779 * t630;
t812 = t716 * t833 - t717 * t823 - t735 * t819;
t811 = t716 * t819 + t717 * t835 + t735 * t831;
t810 = -t716 * t823 + t717 * t836 - t735 * t835;
t808 = (t778 * t821 + t781 * t820) * t804 + t832 * t769;
t807 = (-t778 * t824 - t781 * t834) * t804 - t820 * t769;
t806 = (t778 * t837 + t781 * t824) * t804 + t821 * t769;
t625 = t775 * t636 + t637 * t816 + t643 * t817;
t795 = -t624 * t779 + t782 * t630;
t646 = -mrSges(6,1) * t663 - mrSges(7,1) * t659 + mrSges(7,2) * t657 + mrSges(6,3) * t661 - pkin(5) * t655 - t678 * t833 + t823 * t679 + t819 * t708 + t811 * t717 + t810 * t735;
t647 = mrSges(6,2) * t663 + mrSges(7,2) * t658 - mrSges(6,3) * t660 - mrSges(7,3) * t659 - qJ(6) * t655 - t823 * t678 + t679 * t836 - t708 * t835 + t811 * t716 + t812 * t735;
t704 = Ifges(5,5) * t738 + Ifges(5,6) * t737 + Ifges(5,3) * t760;
t705 = Ifges(5,4) * t738 + Ifges(5,2) * t737 + Ifges(5,6) * t760;
t627 = mrSges(5,2) * t672 - mrSges(5,3) * t667 + Ifges(5,1) * t711 + Ifges(5,4) * t710 + Ifges(5,5) * t743 - pkin(10) * t648 - t776 * t646 + t828 * t647 + t737 * t704 - t760 * t705;
t706 = Ifges(5,1) * t738 + Ifges(5,4) * t737 + Ifges(5,5) * t760;
t632 = Ifges(5,4) * t711 + Ifges(5,2) * t710 + Ifges(5,6) * t743 - t738 * t704 + t760 * t706 - mrSges(5,1) * t672 + mrSges(5,3) * t668 - mrSges(6,1) * t660 + mrSges(6,2) * t661 + mrSges(7,1) * t658 - mrSges(7,3) * t657 - pkin(5) * t792 - qJ(6) * t798 - pkin(4) * t648 + (pkin(5) * t695 + t812) * t717 + (qJ(6) * t695 - t810) * t716 + t831 * t708 + (mrSges(7,2) * pkin(5) + t835) * t679 + (mrSges(7,2) * qJ(6) + t819) * t678;
t638 = t755 * mrSges(4,2) - t747 * t797 + t791;
t621 = -mrSges(3,1) * t726 - mrSges(4,1) * t681 + mrSges(4,2) * t682 + mrSges(3,3) * t713 - pkin(2) * t638 - pkin(3) * t786 - pkin(9) * t794 - t777 * t627 - t780 * t632 + t824 * t754 + t755 * t834 + t820 * t768 + t806 * t769 - t808 * t797;
t622 = t828 * t646 - qJ(3) * t638 + pkin(3) * t639 + t776 * t647 + pkin(10) * t793 + mrSges(5,1) * t667 - mrSges(5,2) * t668 - mrSges(4,3) * t682 + mrSges(4,1) * t691 + Ifges(5,6) * t710 + Ifges(5,5) * t711 - mrSges(3,3) * t712 + pkin(4) * t784 + mrSges(3,2) * t726 - t737 * t706 + t738 * t705 + Ifges(5,3) * t743 + t807 * t769 + t821 * t768 + t824 * t755 + t837 * t754 + t808 * t796;
t789 = pkin(8) * t631 + t621 * t781 + t622 * t778;
t620 = mrSges(3,1) * t712 - mrSges(3,2) * t713 + mrSges(4,2) * t691 - mrSges(4,3) * t681 + t780 * t627 - t777 * t632 - pkin(9) * t639 + pkin(2) * (-t769 * t746 + t788) + qJ(3) * t785 + (-mrSges(4,2) * pkin(2) + t832) * t768 + (mrSges(4,1) * qJ(3) + t820) * t755 + t821 * t754 + (-t806 * t781 + (-pkin(2) * t751 - t807) * t778) * t804;
t619 = -mrSges(2,2) * g(3) - mrSges(2,3) * t764 + Ifges(2,5) * qJDD(1) - t783 * Ifges(2,6) - t778 * t621 + t781 * t622 + (-t625 * t774 - t626 * t775) * pkin(8);
t618 = mrSges(2,1) * g(3) + mrSges(2,3) * t765 + t783 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t625 - t774 * t620 + t789 * t775;
t1 = [-m(1) * g(1) + t795; -m(1) * g(2) + t813; (-m(1) - m(2)) * g(3) + t625; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t813 - t779 * t618 + t782 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t795 + t782 * t618 + t779 * t619; -mrSges(1,1) * g(2) + mrSges(2,1) * t764 + mrSges(1,2) * g(1) - mrSges(2,2) * t765 + Ifges(2,3) * qJDD(1) + pkin(1) * t626 + t775 * t620 + t789 * t774;];
tauB  = t1;
