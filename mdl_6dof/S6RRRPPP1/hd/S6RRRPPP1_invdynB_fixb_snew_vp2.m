% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-05-07 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:58:17
% EndTime: 2019-05-07 03:58:34
% DurationCPUTime: 12.71s
% Computational Cost: add. (189981->365), mult. (401113->438), div. (0->0), fcn. (292645->10), ass. (0->148)
t811 = -2 * qJD(4);
t810 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t809 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t791 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t808 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t790 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t807 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t754 = sin(qJ(1));
t756 = cos(qJ(1));
t746 = t754 * g(1) - t756 * g(2);
t758 = qJD(1) ^ 2;
t731 = -qJDD(1) * pkin(1) - t758 * pkin(8) - t746;
t753 = sin(qJ(2));
t804 = cos(qJ(2));
t783 = t804 * qJD(1);
t773 = qJD(2) * t783;
t741 = t753 * qJDD(1) + t773;
t794 = qJD(1) * t753;
t782 = qJD(2) * t794;
t742 = t804 * qJDD(1) - t782;
t694 = (-t773 - t741) * pkin(9) + (-t742 + t782) * pkin(2) + t731;
t747 = -t756 * g(1) - t754 * g(2);
t732 = -t758 * pkin(1) + qJDD(1) * pkin(8) + t747;
t722 = -t753 * g(3) + t804 * t732;
t740 = (-t804 * pkin(2) - pkin(9) * t753) * qJD(1);
t757 = qJD(2) ^ 2;
t703 = -t757 * pkin(2) + qJDD(2) * pkin(9) + t740 * t783 + t722;
t752 = sin(qJ(3));
t755 = cos(qJ(3));
t669 = t755 * t694 - t752 * t703;
t737 = t755 * qJD(2) - t752 * t794;
t711 = t737 * qJD(3) + t752 * qJDD(2) + t755 * t741;
t738 = t752 * qJD(2) + t755 * t794;
t751 = sin(pkin(6));
t799 = qJ(4) * t751;
t712 = -t737 * pkin(3) - t738 * t799;
t771 = t783 - qJD(3);
t769 = t771 * t751;
t801 = cos(pkin(6));
t763 = t737 * t801 - t769;
t713 = t763 * qJ(4);
t736 = qJDD(3) - t742;
t781 = qJ(4) * t801;
t648 = t736 * pkin(3) - t711 * t781 - t738 * t712 - t771 * t713 + t669;
t670 = t752 * t694 + t755 * t703;
t718 = -t771 * pkin(3) - t738 * t781;
t710 = -t738 * qJD(3) + t755 * qJDD(2) - t752 * t741;
t768 = t801 * t710 + t736 * t751;
t649 = t768 * qJ(4) + t737 * t712 + t771 * t718 + t670;
t721 = -t804 * g(3) - t753 * t732;
t702 = -qJDD(2) * pkin(2) - t757 * pkin(9) + t740 * t794 - t721;
t652 = -t710 * pkin(3) - t711 * t799 - t737 * t713 + t738 * t718 + t702;
t750 = sin(pkin(10));
t800 = cos(pkin(10));
t700 = t800 * t738 + t763 * t750;
t772 = t801 * t800;
t779 = t751 * t800;
t642 = t648 * t772 - t750 * t649 + t652 * t779 + t700 * t811;
t699 = -t737 * t772 + t750 * t738 + t800 * t769;
t671 = -t700 * mrSges(7,2) + t699 * mrSges(7,3);
t673 = t699 * mrSges(5,1) + t700 * mrSges(5,2);
t679 = t800 * t711 + t768 * t750;
t693 = -t751 * t710 + t801 * t736;
t716 = t751 * t737 + t801 * t771;
t672 = t699 * pkin(4) - t700 * qJ(5);
t715 = t716 ^ 2;
t639 = -t693 * pkin(4) - t715 * qJ(5) + t700 * t672 + qJDD(5) - t642;
t674 = -t699 * mrSges(6,2) - t700 * mrSges(6,3);
t798 = t699 * t716;
t805 = 2 * qJD(6);
t633 = t716 * t805 + (t699 * t700 - t693) * qJ(6) + (t679 - t798) * pkin(5) + t639;
t685 = -t699 * mrSges(7,1) - t716 * mrSges(7,2);
t774 = m(7) * t633 - t693 * mrSges(7,3) + t716 * t685;
t764 = -m(6) * t639 - t679 * mrSges(6,1) - t700 * t674 - t774;
t684 = t699 * mrSges(6,1) + t716 * mrSges(6,3);
t795 = t716 * mrSges(5,2) - t699 * mrSges(5,3) - t684;
t802 = -mrSges(7,1) - mrSges(5,3);
t803 = mrSges(5,1) - mrSges(6,2);
t626 = m(5) * t642 - t795 * t716 + (-t671 - t673) * t700 + t803 * t693 + t802 * t679 + t764;
t696 = t699 * t811;
t780 = t750 * t801;
t788 = t751 * t750 * t652 + t648 * t780 + t800 * t649;
t643 = t696 + t788;
t678 = -t710 * t772 + t750 * t711 - t736 * t779;
t681 = -t716 * mrSges(5,1) - t700 * mrSges(5,3);
t765 = t715 * pkin(4) - t693 * qJ(5) - t788;
t638 = 0.2e1 * qJD(5) * t716 + ((2 * qJD(4)) + t672) * t699 + t765;
t686 = t700 * mrSges(6,1) - t716 * mrSges(6,2);
t682 = t700 * pkin(5) + t716 * qJ(6);
t698 = t699 ^ 2;
t806 = -0.2e1 * qJD(5);
t635 = -t678 * pkin(5) - t698 * qJ(6) - t699 * t672 + qJDD(6) + t696 + (t806 - t682) * t716 - t765;
t683 = t700 * mrSges(7,1) + t716 * mrSges(7,3);
t784 = -m(7) * t635 - t693 * mrSges(7,2) + t716 * t683;
t766 = -m(6) * t638 + t693 * mrSges(6,3) - t716 * t686 - t784;
t796 = -t671 - t674;
t629 = m(5) * t643 - t693 * mrSges(5,2) + t716 * t681 + (-t673 + t796) * t699 + (-mrSges(6,1) + t802) * t678 + t766;
t644 = -t751 * t648 + t801 * t652 + qJDD(4);
t761 = (-t679 - t798) * qJ(5) + t644 + (-t716 * pkin(4) + t806) * t700;
t641 = t678 * pkin(4) + t761;
t637 = -t698 * pkin(5) + t699 * t805 - t700 * t682 + (pkin(4) + qJ(6)) * t678 + t761;
t789 = m(7) * t637 + t678 * mrSges(7,3) + t699 * t685;
t770 = m(6) * t641 - t679 * mrSges(6,3) - t700 * t686 + t789;
t630 = m(5) * t644 + (t681 - t683) * t700 + t795 * t699 + (mrSges(5,2) - mrSges(7,2)) * t679 + t803 * t678 + t770;
t619 = (t800 * t626 + t629 * t750) * t751 + t801 * t630;
t739 = (-mrSges(3,1) * t804 + mrSges(3,2) * t753) * qJD(1);
t744 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t794;
t620 = t626 * t772 + t629 * t780 - t751 * t630;
t714 = -t737 * mrSges(4,1) + t738 * mrSges(4,2);
t719 = t771 * mrSges(4,2) + t737 * mrSges(4,3);
t618 = m(4) * t669 + t736 * mrSges(4,1) - t711 * mrSges(4,3) - t738 * t714 - t771 * t719 + t620;
t624 = -t750 * t626 + t800 * t629;
t720 = -t771 * mrSges(4,1) - t738 * mrSges(4,3);
t623 = m(4) * t670 - t736 * mrSges(4,2) + t710 * mrSges(4,3) + t737 * t714 + t771 * t720 + t624;
t775 = -t752 * t618 + t755 * t623;
t612 = m(3) * t722 - qJDD(2) * mrSges(3,2) + t742 * mrSges(3,3) - qJD(2) * t744 + t739 * t783 + t775;
t745 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t783;
t759 = -m(4) * t702 + t710 * mrSges(4,1) - t711 * mrSges(4,2) + t737 * t719 - t738 * t720 - t619;
t617 = m(3) * t721 + qJDD(2) * mrSges(3,1) - t741 * mrSges(3,3) + qJD(2) * t745 - t739 * t794 + t759;
t776 = t804 * t612 - t753 * t617;
t606 = m(2) * t747 - t758 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t776;
t613 = t755 * t618 + t752 * t623;
t760 = -m(3) * t731 + t742 * mrSges(3,1) - t741 * mrSges(3,2) - t744 * t794 + t745 * t783 - t613;
t609 = m(2) * t746 + qJDD(1) * mrSges(2,1) - t758 * mrSges(2,2) + t760;
t797 = t754 * t606 + t756 * t609;
t607 = t753 * t612 + t804 * t617;
t787 = t699 * t790 - t700 * t791 + t716 * t807;
t786 = t699 * t808 + t700 * t809 - t716 * t790;
t785 = t699 * t809 - t700 * t810 + t716 * t791;
t777 = t756 * t606 - t754 * t609;
t730 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t753 + t804 * Ifges(3,4)) * qJD(1);
t729 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t753 + t804 * Ifges(3,2)) * qJD(1);
t728 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t753 + t804 * Ifges(3,6)) * qJD(1);
t707 = Ifges(4,1) * t738 + Ifges(4,4) * t737 - Ifges(4,5) * t771;
t706 = Ifges(4,4) * t738 + Ifges(4,2) * t737 - Ifges(4,6) * t771;
t705 = Ifges(4,5) * t738 + Ifges(4,6) * t737 - Ifges(4,3) * t771;
t632 = t679 * mrSges(7,1) + t700 * t671 + t774;
t631 = -t678 * mrSges(6,2) - t679 * mrSges(7,2) - t700 * t683 - t699 * t684 + t770;
t621 = mrSges(6,1) * t639 + mrSges(7,1) * t633 + mrSges(5,2) * t644 - mrSges(7,2) * t637 - mrSges(5,3) * t642 - mrSges(6,3) * t641 + pkin(5) * t632 - qJ(5) * t631 - t678 * t809 + t679 * t810 + t791 * t693 + t787 * t699 + t786 * t716;
t615 = -mrSges(5,1) * t644 + mrSges(5,3) * t643 - mrSges(6,1) * t638 + mrSges(6,2) * t641 + mrSges(7,1) * t635 - mrSges(7,3) * t637 - pkin(5) * (t699 * t671 + t784) - qJ(6) * t789 - pkin(4) * t631 + t785 * t716 + (qJ(6) * t683 + t787) * t700 + t790 * t693 + (qJ(6) * mrSges(7,2) + t809) * t679 + (-pkin(5) * mrSges(7,1) + t808) * t678;
t614 = mrSges(5,1) * t642 - mrSges(5,2) * t643 + mrSges(6,2) * t639 - mrSges(6,3) * t638 + mrSges(7,2) * t635 - mrSges(7,3) * t633 - qJ(6) * t632 + pkin(4) * (t716 * t684 + t764) + qJ(5) * t766 + (-pkin(4) * t671 + t786) * t700 + (-pkin(4) * mrSges(6,2) + t807) * t693 + (-pkin(4) * mrSges(7,1) + t791) * t679 + (qJ(5) * t796 - t785) * t699 + (qJ(5) * (-mrSges(6,1) - mrSges(7,1)) - t790) * t678;
t603 = Ifges(4,1) * t711 + Ifges(4,4) * t710 + Ifges(4,5) * t736 + t737 * t705 + t771 * t706 + mrSges(4,2) * t702 - mrSges(4,3) * t669 + t800 * t621 - t750 * t615 + (-t751 * t619 - t801 * t620) * qJ(4);
t602 = -mrSges(4,1) * t702 + mrSges(4,3) * t670 + Ifges(4,4) * t711 + Ifges(4,2) * t710 + Ifges(4,6) * t736 - pkin(3) * t619 - t751 * t614 + t615 * t772 + t621 * t780 + t624 * t781 - t738 * t705 - t771 * t707;
t601 = Ifges(3,4) * t741 + Ifges(3,2) * t742 + Ifges(3,6) * qJDD(2) - t728 * t794 + qJD(2) * t730 - mrSges(3,1) * t731 + mrSges(3,3) * t722 - Ifges(4,5) * t711 - Ifges(4,6) * t710 - Ifges(4,3) * t736 - t738 * t706 + t737 * t707 - mrSges(4,1) * t669 + mrSges(4,2) * t670 - t801 * t614 - pkin(3) * t620 - pkin(2) * t613 + (-qJ(4) * t624 - t800 * t615 - t750 * t621) * t751;
t600 = mrSges(3,2) * t731 - mrSges(3,3) * t721 + Ifges(3,1) * t741 + Ifges(3,4) * t742 + Ifges(3,5) * qJDD(2) - pkin(9) * t613 - qJD(2) * t729 - t752 * t602 + t755 * t603 + t728 * t783;
t599 = Ifges(2,6) * qJDD(1) + t758 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t747 - Ifges(3,5) * t741 - Ifges(3,6) * t742 - Ifges(3,3) * qJDD(2) - t729 * t794 + t730 * t783 - mrSges(3,1) * t721 + mrSges(3,2) * t722 - t752 * t603 - t755 * t602 - pkin(2) * t759 - pkin(9) * t775 - pkin(1) * t607;
t598 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - t758 * Ifges(2,6) - pkin(8) * t607 + t804 * t600 - t753 * t601;
t1 = [-m(1) * g(1) + t777; -m(1) * g(2) + t797; (-m(1) - m(2)) * g(3) + t607; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t797 + t756 * t598 - t754 * t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t777 + t754 * t598 + t756 * t599; -mrSges(1,1) * g(2) + mrSges(2,1) * t746 + mrSges(1,2) * g(1) - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) + pkin(1) * t760 + pkin(8) * t776 + t753 * t600 + t804 * t601;];
tauB  = t1;
