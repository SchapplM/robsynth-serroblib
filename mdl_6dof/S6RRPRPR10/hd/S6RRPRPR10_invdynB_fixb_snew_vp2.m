% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:38:56
% EndTime: 2019-05-06 15:39:24
% DurationCPUTime: 19.74s
% Computational Cost: add. (305350->378), mult. (696483->473), div. (0->0), fcn. (553273->12), ass. (0->155)
t811 = Ifges(5,1) + Ifges(6,2);
t805 = Ifges(5,4) + Ifges(6,6);
t804 = Ifges(5,5) - Ifges(6,4);
t810 = -Ifges(5,2) - Ifges(6,3);
t803 = Ifges(5,6) - Ifges(6,5);
t809 = Ifges(5,3) + Ifges(6,1);
t808 = -2 * qJD(5);
t807 = cos(qJ(4));
t767 = cos(pkin(6));
t806 = t767 * g(3);
t761 = t767 * qJD(1) + qJD(2);
t764 = sin(pkin(11));
t766 = cos(pkin(11));
t770 = sin(qJ(2));
t765 = sin(pkin(6));
t791 = qJD(1) * t765;
t788 = t770 * t791;
t739 = t766 * t761 - t764 * t788;
t740 = t764 * t761 + t766 * t788;
t769 = sin(qJ(4));
t719 = -t807 * t739 + t769 * t740;
t773 = cos(qJ(2));
t790 = qJD(1) * t773;
t787 = t765 * t790;
t754 = -qJD(4) + t787;
t802 = t719 * t754;
t801 = t765 * t770;
t800 = t765 * t773;
t799 = t767 * t770;
t798 = t767 * t773;
t771 = sin(qJ(1));
t774 = cos(qJ(1));
t756 = t771 * g(1) - t774 * g(2);
t775 = qJD(1) ^ 2;
t747 = t775 * t765 * pkin(8) + qJDD(1) * pkin(1) + t756;
t757 = -t774 * g(1) - t771 * g(2);
t789 = qJDD(1) * t765;
t748 = -t775 * pkin(1) + pkin(8) * t789 + t757;
t792 = t747 * t799 + t773 * t748;
t716 = -g(3) * t801 + t792;
t745 = t761 * mrSges(3,1) - mrSges(3,3) * t788;
t750 = (-mrSges(3,1) * t773 + mrSges(3,2) * t770) * t791;
t752 = -qJD(2) * t788 + t773 * t789;
t760 = t767 * qJDD(1) + qJDD(2);
t749 = (-pkin(2) * t773 - qJ(3) * t770) * t791;
t759 = t761 ^ 2;
t698 = -t759 * pkin(2) + t760 * qJ(3) + (-g(3) * t770 + t749 * t790) * t765 + t792;
t751 = (qJD(2) * t790 + qJDD(1) * t770) * t765;
t699 = -t752 * pkin(2) - t806 - t751 * qJ(3) + (-t747 + (pkin(2) * t770 - qJ(3) * t773) * t761 * qJD(1)) * t765;
t660 = -0.2e1 * qJD(3) * t740 - t764 * t698 + t766 * t699;
t727 = t766 * t751 + t764 * t760;
t656 = (-t739 * t787 - t727) * pkin(9) + (t739 * t740 - t752) * pkin(3) + t660;
t661 = 0.2e1 * qJD(3) * t739 + t766 * t698 + t764 * t699;
t726 = -t764 * t751 + t766 * t760;
t728 = -pkin(3) * t787 - t740 * pkin(9);
t738 = t739 ^ 2;
t659 = -t738 * pkin(3) + t726 * pkin(9) + t728 * t787 + t661;
t651 = t807 * t656 - t769 * t659;
t679 = -t719 * qJD(4) + t769 * t726 + t807 * t727;
t720 = t769 * t739 + t807 * t740;
t692 = t719 * mrSges(5,1) + t720 * mrSges(5,2);
t703 = t719 * mrSges(6,1) + t754 * mrSges(6,3);
t705 = t754 * mrSges(5,2) - t719 * mrSges(5,3);
t744 = qJDD(4) - t752;
t691 = t719 * pkin(4) - t720 * qJ(5);
t753 = t754 ^ 2;
t650 = -t744 * pkin(4) - t753 * qJ(5) + t720 * t691 + qJDD(5) - t651;
t645 = (t719 * t720 - t744) * pkin(10) + (t679 - t802) * pkin(5) + t650;
t678 = t720 * qJD(4) - t807 * t726 + t769 * t727;
t707 = t720 * pkin(5) + t754 * pkin(10);
t718 = t719 ^ 2;
t715 = -g(3) * t800 + t747 * t798 - t770 * t748;
t697 = -t760 * pkin(2) - t759 * qJ(3) + t749 * t788 + qJDD(3) - t715;
t665 = -t726 * pkin(3) - t738 * pkin(9) + t740 * t728 + t697;
t776 = (-t679 - t802) * qJ(5) + t665 + (-t754 * pkin(4) + t808) * t720;
t648 = -t720 * t707 - t718 * pkin(5) + t776 + (pkin(4) + pkin(10)) * t678;
t768 = sin(qJ(6));
t772 = cos(qJ(6));
t643 = t772 * t645 - t768 * t648;
t701 = t772 * t719 + t768 * t754;
t664 = t701 * qJD(6) + t768 * t678 + t772 * t744;
t702 = t768 * t719 - t772 * t754;
t672 = -t701 * mrSges(7,1) + t702 * mrSges(7,2);
t677 = qJDD(6) + t679;
t717 = qJD(6) + t720;
t680 = -t717 * mrSges(7,2) + t701 * mrSges(7,3);
t641 = m(7) * t643 + t677 * mrSges(7,1) - t664 * mrSges(7,3) - t702 * t672 + t717 * t680;
t644 = t768 * t645 + t772 * t648;
t663 = -t702 * qJD(6) + t772 * t678 - t768 * t744;
t681 = t717 * mrSges(7,1) - t702 * mrSges(7,3);
t642 = m(7) * t644 - t677 * mrSges(7,2) + t663 * mrSges(7,3) + t701 * t672 - t717 * t681;
t633 = t772 * t641 + t768 * t642;
t693 = -t719 * mrSges(6,2) - t720 * mrSges(6,3);
t781 = -m(6) * t650 - t679 * mrSges(6,1) - t720 * t693 - t633;
t628 = m(5) * t651 - t679 * mrSges(5,3) - t720 * t692 + (t703 - t705) * t754 + (mrSges(5,1) - mrSges(6,2)) * t744 + t781;
t652 = t769 * t656 + t807 * t659;
t706 = -t754 * mrSges(5,1) - t720 * mrSges(5,3);
t780 = -t753 * pkin(4) + t744 * qJ(5) - t719 * t691 + t652;
t649 = 0.2e1 * qJD(5) * t754 - t780;
t704 = t720 * mrSges(6,1) - t754 * mrSges(6,2);
t647 = -t678 * pkin(5) - t718 * pkin(10) + (t808 - t707) * t754 + t780;
t782 = -m(7) * t647 + t663 * mrSges(7,1) - t664 * mrSges(7,2) + t701 * t680 - t702 * t681;
t779 = -m(6) * t649 + t744 * mrSges(6,3) - t754 * t704 - t782;
t638 = m(5) * t652 - t744 * mrSges(5,2) + t754 * t706 + (-t692 - t693) * t719 + (-mrSges(5,3) - mrSges(6,1)) * t678 + t779;
t626 = t807 * t628 + t769 * t638;
t721 = -t739 * mrSges(4,1) + t740 * mrSges(4,2);
t724 = mrSges(4,2) * t787 + t739 * mrSges(4,3);
t624 = m(4) * t660 - t752 * mrSges(4,1) - t727 * mrSges(4,3) - t740 * t721 - t724 * t787 + t626;
t725 = -mrSges(4,1) * t787 - t740 * mrSges(4,3);
t784 = -t769 * t628 + t807 * t638;
t625 = m(4) * t661 + t752 * mrSges(4,2) + t726 * mrSges(4,3) + t739 * t721 + t725 * t787 + t784;
t785 = -t764 * t624 + t766 * t625;
t616 = m(3) * t716 - t760 * mrSges(3,2) + t752 * mrSges(3,3) - t761 * t745 + t750 * t787 + t785;
t619 = t766 * t624 + t764 * t625;
t732 = -t765 * t747 - t806;
t746 = -t761 * mrSges(3,2) + mrSges(3,3) * t787;
t618 = m(3) * t732 - t752 * mrSges(3,1) + t751 * mrSges(3,2) + (t745 * t770 - t746 * t773) * t791 + t619;
t654 = t678 * pkin(4) + t776;
t796 = -t768 * t641 + t772 * t642;
t632 = m(6) * t654 - t678 * mrSges(6,2) - t679 * mrSges(6,3) - t719 * t703 - t720 * t704 + t796;
t778 = m(5) * t665 + t678 * mrSges(5,1) + t679 * mrSges(5,2) + t719 * t705 + t720 * t706 + t632;
t777 = -m(4) * t697 + t726 * mrSges(4,1) - t727 * mrSges(4,2) + t739 * t724 - t740 * t725 - t778;
t631 = m(3) * t715 + t760 * mrSges(3,1) - t751 * mrSges(3,3) + t761 * t746 - t750 * t788 + t777;
t607 = t616 * t799 - t765 * t618 + t631 * t798;
t605 = m(2) * t756 + qJDD(1) * mrSges(2,1) - t775 * mrSges(2,2) + t607;
t611 = t773 * t616 - t770 * t631;
t610 = m(2) * t757 - t775 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t611;
t797 = t774 * t605 + t771 * t610;
t795 = t803 * t719 - t804 * t720 + t809 * t754;
t794 = t810 * t719 + t805 * t720 - t803 * t754;
t793 = t805 * t719 - t811 * t720 + t804 * t754;
t606 = t616 * t801 + t767 * t618 + t631 * t800;
t786 = -t771 * t605 + t774 * t610;
t666 = Ifges(7,5) * t702 + Ifges(7,6) * t701 + Ifges(7,3) * t717;
t668 = Ifges(7,1) * t702 + Ifges(7,4) * t701 + Ifges(7,5) * t717;
t634 = -mrSges(7,1) * t647 + mrSges(7,3) * t644 + Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t677 - t702 * t666 + t717 * t668;
t667 = Ifges(7,4) * t702 + Ifges(7,2) * t701 + Ifges(7,6) * t717;
t635 = mrSges(7,2) * t647 - mrSges(7,3) * t643 + Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t677 + t701 * t666 - t717 * t667;
t612 = -mrSges(5,1) * t665 - mrSges(6,1) * t649 + mrSges(6,2) * t654 + mrSges(5,3) * t652 - pkin(4) * t632 - pkin(5) * t782 - pkin(10) * t796 - t772 * t634 - t768 * t635 + t810 * t678 + t805 * t679 + t795 * t720 + t803 * t744 + t793 * t754;
t620 = mrSges(6,1) * t650 + mrSges(7,1) * t643 + mrSges(5,2) * t665 - mrSges(7,2) * t644 - mrSges(5,3) * t651 - mrSges(6,3) * t654 + Ifges(7,5) * t664 + Ifges(7,6) * t663 + Ifges(7,3) * t677 + pkin(5) * t633 - qJ(5) * t632 + t702 * t667 - t701 * t668 + t794 * t754 + t804 * t744 + t795 * t719 + t811 * t679 - t805 * t678;
t709 = Ifges(4,5) * t740 + Ifges(4,6) * t739 - Ifges(4,3) * t787;
t711 = Ifges(4,1) * t740 + Ifges(4,4) * t739 - Ifges(4,5) * t787;
t601 = -mrSges(4,1) * t697 + mrSges(4,3) * t661 + Ifges(4,4) * t727 + Ifges(4,2) * t726 - Ifges(4,6) * t752 - pkin(3) * t778 + pkin(9) * t784 + t807 * t612 + t769 * t620 - t740 * t709 - t711 * t787;
t710 = Ifges(4,4) * t740 + Ifges(4,2) * t739 - Ifges(4,6) * t787;
t603 = mrSges(4,2) * t697 - mrSges(4,3) * t660 + Ifges(4,1) * t727 + Ifges(4,4) * t726 - Ifges(4,5) * t752 - pkin(9) * t626 - t769 * t612 + t807 * t620 + t739 * t709 + t710 * t787;
t729 = Ifges(3,3) * t761 + (Ifges(3,5) * t770 + Ifges(3,6) * t773) * t791;
t730 = Ifges(3,6) * t761 + (Ifges(3,4) * t770 + Ifges(3,2) * t773) * t791;
t600 = mrSges(3,2) * t732 - mrSges(3,3) * t715 + Ifges(3,1) * t751 + Ifges(3,4) * t752 + Ifges(3,5) * t760 - qJ(3) * t619 - t764 * t601 + t766 * t603 + t729 * t787 - t761 * t730;
t731 = Ifges(3,5) * t761 + (Ifges(3,1) * t770 + Ifges(3,4) * t773) * t791;
t602 = (qJ(5) * mrSges(6,1) + t803) * t678 - t804 * t679 - pkin(4) * (t754 * t703 + t781) + (qJ(5) * t693 + t793) * t719 - t794 * t720 - t729 * t788 - qJ(5) * t779 - t772 * t635 + t768 * t634 + Ifges(3,6) * t760 + t761 * t731 + Ifges(3,4) * t751 + t739 * t711 - t740 * t710 - Ifges(4,6) * t726 - Ifges(4,5) * t727 - mrSges(3,1) * t732 + mrSges(3,3) * t716 + mrSges(4,2) * t661 - mrSges(4,1) * t660 - mrSges(5,1) * t651 + mrSges(5,2) * t652 + mrSges(6,3) * t649 - mrSges(6,2) * t650 + pkin(10) * t633 + (Ifges(3,2) + Ifges(4,3)) * t752 - pkin(3) * t626 + (pkin(4) * mrSges(6,2) - t809) * t744 - pkin(2) * t619;
t783 = pkin(8) * t611 + t600 * t770 + t602 * t773;
t599 = Ifges(3,5) * t751 + Ifges(3,6) * t752 + Ifges(3,3) * t760 + mrSges(3,1) * t715 - mrSges(3,2) * t716 + t764 * t603 + t766 * t601 + pkin(2) * t777 + qJ(3) * t785 + (t730 * t770 - t731 * t773) * t791;
t598 = -mrSges(2,2) * g(3) - mrSges(2,3) * t756 + Ifges(2,5) * qJDD(1) - t775 * Ifges(2,6) + t773 * t600 - t770 * t602 + (-t606 * t765 - t607 * t767) * pkin(8);
t597 = mrSges(2,1) * g(3) + mrSges(2,3) * t757 + t775 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t606 - t765 * t599 + t783 * t767;
t1 = [-m(1) * g(1) + t786; -m(1) * g(2) + t797; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t797 - t771 * t597 + t774 * t598; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t786 + t774 * t597 + t771 * t598; -mrSges(1,1) * g(2) + mrSges(2,1) * t756 + mrSges(1,2) * g(1) - mrSges(2,2) * t757 + Ifges(2,3) * qJDD(1) + pkin(1) * t607 + t767 * t599 + t783 * t765;];
tauB  = t1;
