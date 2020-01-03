% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:38:20
% EndTime: 2019-12-31 21:38:36
% DurationCPUTime: 15.99s
% Computational Cost: add. (261171->326), mult. (572723->429), div. (0->0), fcn. (445203->12), ass. (0->137)
t776 = sin(pkin(5));
t781 = sin(qJ(2));
t784 = cos(qJ(2));
t799 = qJD(1) * qJD(2);
t761 = (-qJDD(1) * t784 + t781 * t799) * t776;
t810 = cos(qJ(3));
t809 = pkin(7) * t776;
t778 = cos(pkin(5));
t808 = t778 * g(3);
t807 = t776 * t781;
t806 = t776 * t784;
t805 = t778 * t781;
t804 = t778 * t784;
t782 = sin(qJ(1));
t785 = cos(qJ(1));
t767 = t782 * g(1) - t785 * g(2);
t786 = qJD(1) ^ 2;
t756 = qJDD(1) * pkin(1) + t786 * t809 + t767;
t768 = -t785 * g(1) - t782 * g(2);
t757 = -t786 * pkin(1) + qJDD(1) * t809 + t768;
t802 = t756 * t805 + t784 * t757;
t732 = -g(3) * t807 + t802;
t771 = t778 * qJD(1) + qJD(2);
t801 = qJD(1) * t776;
t798 = t781 * t801;
t754 = t771 * mrSges(3,1) - mrSges(3,3) * t798;
t758 = (-mrSges(3,1) * t784 + mrSges(3,2) * t781) * t801;
t770 = t778 * qJDD(1) + qJDD(2);
t759 = (-pkin(2) * t784 - pkin(8) * t781) * t801;
t769 = t771 ^ 2;
t800 = qJD(1) * t784;
t712 = -t769 * pkin(2) + t770 * pkin(8) + (-g(3) * t781 + t759 * t800) * t776 + t802;
t760 = (qJDD(1) * t781 + t784 * t799) * t776;
t713 = t761 * pkin(2) - t760 * pkin(8) - t808 + (-t756 + (pkin(2) * t781 - pkin(8) * t784) * t771 * qJD(1)) * t776;
t780 = sin(qJ(3));
t696 = t810 * t712 + t780 * t713;
t749 = -t810 * t771 + t780 * t798;
t750 = t780 * t771 + t810 * t798;
t733 = t749 * pkin(3) - t750 * qJ(4);
t753 = qJDD(3) + t761;
t797 = t776 * t800;
t766 = qJD(3) - t797;
t765 = t766 ^ 2;
t688 = -t765 * pkin(3) + t753 * qJ(4) - t749 * t733 + t696;
t731 = -g(3) * t806 + t756 * t804 - t781 * t757;
t711 = -t770 * pkin(2) - t769 * pkin(8) + t759 * t798 - t731;
t729 = t750 * qJD(3) + t780 * t760 - t810 * t770;
t730 = -t749 * qJD(3) + t810 * t760 + t780 * t770;
t691 = (t749 * t766 - t730) * qJ(4) + (t750 * t766 + t729) * pkin(3) + t711;
t775 = sin(pkin(10));
t777 = cos(pkin(10));
t739 = t777 * t750 + t775 * t766;
t683 = -0.2e1 * qJD(4) * t739 - t775 * t688 + t777 * t691;
t718 = t777 * t730 + t775 * t753;
t738 = -t775 * t750 + t777 * t766;
t681 = (t749 * t738 - t718) * pkin(9) + (t738 * t739 + t729) * pkin(4) + t683;
t684 = 0.2e1 * qJD(4) * t738 + t777 * t688 + t775 * t691;
t717 = -t775 * t730 + t777 * t753;
t722 = t749 * pkin(4) - t739 * pkin(9);
t737 = t738 ^ 2;
t682 = -t737 * pkin(4) + t717 * pkin(9) - t749 * t722 + t684;
t779 = sin(qJ(5));
t783 = cos(qJ(5));
t679 = t783 * t681 - t779 * t682;
t714 = t783 * t738 - t779 * t739;
t694 = t714 * qJD(5) + t779 * t717 + t783 * t718;
t715 = t779 * t738 + t783 * t739;
t701 = -t714 * mrSges(6,1) + t715 * mrSges(6,2);
t748 = qJD(5) + t749;
t702 = -t748 * mrSges(6,2) + t714 * mrSges(6,3);
t727 = qJDD(5) + t729;
t676 = m(6) * t679 + t727 * mrSges(6,1) - t694 * mrSges(6,3) - t715 * t701 + t748 * t702;
t680 = t779 * t681 + t783 * t682;
t693 = -t715 * qJD(5) + t783 * t717 - t779 * t718;
t703 = t748 * mrSges(6,1) - t715 * mrSges(6,3);
t677 = m(6) * t680 - t727 * mrSges(6,2) + t693 * mrSges(6,3) + t714 * t701 - t748 * t703;
t668 = t783 * t676 + t779 * t677;
t719 = -t738 * mrSges(5,1) + t739 * mrSges(5,2);
t793 = -t749 * mrSges(5,2) + t738 * mrSges(5,3);
t666 = m(5) * t683 + t729 * mrSges(5,1) - t718 * mrSges(5,3) - t739 * t719 + t749 * t793 + t668;
t721 = t749 * mrSges(5,1) - t739 * mrSges(5,3);
t794 = -t779 * t676 + t783 * t677;
t667 = m(5) * t684 - t729 * mrSges(5,2) + t717 * mrSges(5,3) + t738 * t719 - t749 * t721 + t794;
t664 = -t775 * t666 + t777 * t667;
t734 = t749 * mrSges(4,1) + t750 * mrSges(4,2);
t741 = t766 * mrSges(4,1) - t750 * mrSges(4,3);
t662 = m(4) * t696 - t753 * mrSges(4,2) - t729 * mrSges(4,3) - t749 * t734 - t766 * t741 + t664;
t695 = -t780 * t712 + t810 * t713;
t687 = -t753 * pkin(3) - t765 * qJ(4) + t750 * t733 + qJDD(4) - t695;
t685 = -t717 * pkin(4) - t737 * pkin(9) + t739 * t722 + t687;
t791 = m(6) * t685 - t693 * mrSges(6,1) + t694 * mrSges(6,2) - t714 * t702 + t715 * t703;
t678 = m(5) * t687 - t717 * mrSges(5,1) + t718 * mrSges(5,2) + t739 * t721 - t738 * t793 + t791;
t740 = -t766 * mrSges(4,2) - t749 * mrSges(4,3);
t672 = m(4) * t695 + t753 * mrSges(4,1) - t730 * mrSges(4,3) - t750 * t734 + t766 * t740 - t678;
t795 = t810 * t662 - t780 * t672;
t651 = m(3) * t732 - t770 * mrSges(3,2) - t761 * mrSges(3,3) - t771 * t754 + t758 * t797 + t795;
t654 = t780 * t662 + t810 * t672;
t745 = -t776 * t756 - t808;
t755 = -t771 * mrSges(3,2) + mrSges(3,3) * t797;
t653 = m(3) * t745 + t761 * mrSges(3,1) + t760 * mrSges(3,2) + (t754 * t781 - t755 * t784) * t801 + t654;
t663 = t777 * t666 + t775 * t667;
t789 = -m(4) * t711 - t729 * mrSges(4,1) - t730 * mrSges(4,2) - t749 * t740 - t750 * t741 - t663;
t659 = m(3) * t731 + t770 * mrSges(3,1) - t760 * mrSges(3,3) + t771 * t755 - t758 * t798 + t789;
t640 = t651 * t805 - t776 * t653 + t659 * t804;
t637 = m(2) * t767 + qJDD(1) * mrSges(2,1) - t786 * mrSges(2,2) + t640;
t646 = t784 * t651 - t781 * t659;
t644 = m(2) * t768 - t786 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t646;
t803 = t785 * t637 + t782 * t644;
t639 = t651 * t807 + t778 * t653 + t659 * t806;
t796 = -t782 * t637 + t785 * t644;
t697 = Ifges(6,5) * t715 + Ifges(6,6) * t714 + Ifges(6,3) * t748;
t699 = Ifges(6,1) * t715 + Ifges(6,4) * t714 + Ifges(6,5) * t748;
t669 = -mrSges(6,1) * t685 + mrSges(6,3) * t680 + Ifges(6,4) * t694 + Ifges(6,2) * t693 + Ifges(6,6) * t727 - t715 * t697 + t748 * t699;
t698 = Ifges(6,4) * t715 + Ifges(6,2) * t714 + Ifges(6,6) * t748;
t670 = mrSges(6,2) * t685 - mrSges(6,3) * t679 + Ifges(6,1) * t694 + Ifges(6,4) * t693 + Ifges(6,5) * t727 + t714 * t697 - t748 * t698;
t704 = Ifges(5,5) * t739 + Ifges(5,6) * t738 + Ifges(5,3) * t749;
t706 = Ifges(5,1) * t739 + Ifges(5,4) * t738 + Ifges(5,5) * t749;
t655 = -mrSges(5,1) * t687 + mrSges(5,3) * t684 + Ifges(5,4) * t718 + Ifges(5,2) * t717 + Ifges(5,6) * t729 - pkin(4) * t791 + pkin(9) * t794 + t783 * t669 + t779 * t670 - t739 * t704 + t749 * t706;
t705 = Ifges(5,4) * t739 + Ifges(5,2) * t738 + Ifges(5,6) * t749;
t656 = mrSges(5,2) * t687 - mrSges(5,3) * t683 + Ifges(5,1) * t718 + Ifges(5,4) * t717 + Ifges(5,5) * t729 - pkin(9) * t668 - t779 * t669 + t783 * t670 + t738 * t704 - t749 * t705;
t723 = Ifges(4,5) * t750 - Ifges(4,6) * t749 + Ifges(4,3) * t766;
t724 = Ifges(4,4) * t750 - Ifges(4,2) * t749 + Ifges(4,6) * t766;
t641 = mrSges(4,2) * t711 - mrSges(4,3) * t695 + Ifges(4,1) * t730 - Ifges(4,4) * t729 + Ifges(4,5) * t753 - qJ(4) * t663 - t775 * t655 + t777 * t656 - t749 * t723 - t766 * t724;
t725 = Ifges(4,1) * t750 - Ifges(4,4) * t749 + Ifges(4,5) * t766;
t788 = mrSges(6,1) * t679 - mrSges(6,2) * t680 + Ifges(6,5) * t694 + Ifges(6,6) * t693 + Ifges(6,3) * t727 + t715 * t698 - t714 * t699;
t647 = -t788 + (-Ifges(4,2) - Ifges(5,3)) * t729 + t766 * t725 - t750 * t723 + Ifges(4,6) * t753 + t738 * t706 - t739 * t705 + Ifges(4,4) * t730 - Ifges(5,6) * t717 - Ifges(5,5) * t718 - mrSges(4,1) * t711 + mrSges(4,3) * t696 + mrSges(5,2) * t684 - mrSges(5,1) * t683 - pkin(4) * t668 - pkin(3) * t663;
t743 = Ifges(3,6) * t771 + (Ifges(3,4) * t781 + Ifges(3,2) * t784) * t801;
t744 = Ifges(3,5) * t771 + (Ifges(3,1) * t781 + Ifges(3,4) * t784) * t801;
t631 = Ifges(3,5) * t760 - Ifges(3,6) * t761 + Ifges(3,3) * t770 + mrSges(3,1) * t731 - mrSges(3,2) * t732 + t780 * t641 + t810 * t647 + pkin(2) * t789 + pkin(8) * t795 + (t743 * t781 - t744 * t784) * t801;
t742 = Ifges(3,3) * t771 + (Ifges(3,5) * t781 + Ifges(3,6) * t784) * t801;
t633 = mrSges(3,2) * t745 - mrSges(3,3) * t731 + Ifges(3,1) * t760 - Ifges(3,4) * t761 + Ifges(3,5) * t770 - pkin(8) * t654 + t810 * t641 - t780 * t647 + t742 * t797 - t771 * t743;
t787 = mrSges(4,1) * t695 - mrSges(4,2) * t696 + Ifges(4,5) * t730 - Ifges(4,6) * t729 + Ifges(4,3) * t753 - pkin(3) * t678 + qJ(4) * t664 + t777 * t655 + t775 * t656 + t750 * t724 + t749 * t725;
t635 = -mrSges(3,1) * t745 + mrSges(3,3) * t732 + Ifges(3,4) * t760 - Ifges(3,2) * t761 + Ifges(3,6) * t770 - pkin(2) * t654 - t742 * t798 + t771 * t744 - t787;
t790 = mrSges(2,1) * t767 - mrSges(2,2) * t768 + Ifges(2,3) * qJDD(1) + pkin(1) * t640 + t778 * t631 + t633 * t807 + t635 * t806 + t646 * t809;
t629 = -mrSges(2,2) * g(3) - mrSges(2,3) * t767 + Ifges(2,5) * qJDD(1) - t786 * Ifges(2,6) + t784 * t633 - t781 * t635 + (-t639 * t776 - t640 * t778) * pkin(7);
t628 = mrSges(2,1) * g(3) + mrSges(2,3) * t768 + t786 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t639 - t776 * t631 + (pkin(7) * t646 + t633 * t781 + t635 * t784) * t778;
t1 = [-m(1) * g(1) + t796; -m(1) * g(2) + t803; (-m(1) - m(2)) * g(3) + t639; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t803 - t782 * t628 + t785 * t629; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t796 + t785 * t628 + t782 * t629; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t790; t790; t631; t787; t678; t788;];
tauJB = t1;
