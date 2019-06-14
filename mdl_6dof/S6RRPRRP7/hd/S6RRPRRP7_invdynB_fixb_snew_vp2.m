% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:13:00
% EndTime: 2019-05-06 18:13:09
% DurationCPUTime: 5.58s
% Computational Cost: add. (55478->340), mult. (112939->404), div. (0->0), fcn. (69951->8), ass. (0->132)
t790 = Ifges(3,1) + Ifges(4,1);
t789 = Ifges(6,1) + Ifges(7,1);
t779 = Ifges(3,4) - Ifges(4,5);
t778 = Ifges(6,4) - Ifges(7,5);
t776 = Ifges(3,5) + Ifges(4,4);
t788 = -Ifges(6,5) - Ifges(7,4);
t787 = Ifges(3,2) + Ifges(4,3);
t786 = Ifges(6,2) + Ifges(7,3);
t775 = Ifges(3,6) - Ifges(4,6);
t774 = Ifges(6,6) - Ifges(7,6);
t785 = Ifges(3,3) + Ifges(4,2);
t784 = -Ifges(6,3) - Ifges(7,2);
t739 = sin(qJ(4));
t740 = sin(qJ(2));
t742 = cos(qJ(4));
t743 = cos(qJ(2));
t694 = (t739 * t740 + t742 * t743) * qJD(1);
t783 = 2 * qJD(3);
t782 = cos(qJ(5));
t781 = mrSges(3,3) + mrSges(4,2);
t780 = -mrSges(6,3) - mrSges(7,2);
t746 = qJD(1) ^ 2;
t773 = t743 ^ 2 * t746;
t741 = sin(qJ(1));
t744 = cos(qJ(1));
t719 = -g(1) * t744 - g(2) * t741;
t697 = -pkin(1) * t746 + qJDD(1) * pkin(7) + t719;
t678 = -g(3) * t740 + t743 * t697;
t708 = (-mrSges(3,1) * t743 + mrSges(3,2) * t740) * qJD(1);
t762 = qJD(1) * qJD(2);
t759 = t740 * t762;
t710 = qJDD(1) * t743 - t759;
t764 = qJD(1) * t740;
t713 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t764;
t706 = (-pkin(2) * t743 - qJ(3) * t740) * qJD(1);
t745 = qJD(2) ^ 2;
t763 = qJD(1) * t743;
t657 = -pkin(2) * t745 + qJDD(2) * qJ(3) + qJD(2) * t783 + t706 * t763 + t678;
t707 = (-mrSges(4,1) * t743 - mrSges(4,3) * t740) * qJD(1);
t714 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t764;
t717 = -qJD(2) * pkin(3) - pkin(8) * t764;
t637 = -pkin(3) * t773 - pkin(8) * t710 + qJD(2) * t717 + t657;
t677 = -t743 * g(3) - t740 * t697;
t661 = -qJDD(2) * pkin(2) - t745 * qJ(3) + t706 * t764 + qJDD(3) - t677;
t760 = t743 * t762;
t709 = qJDD(1) * t740 + t760;
t638 = (-t709 + t760) * pkin(8) + (-t740 * t743 * t746 - qJDD(2)) * pkin(3) + t661;
t627 = t742 * t637 + t739 * t638;
t695 = (-t739 * t743 + t740 * t742) * qJD(1);
t662 = -qJD(4) * t695 - t709 * t739 - t710 * t742;
t673 = mrSges(5,1) * t694 + mrSges(5,2) * t695;
t731 = -qJD(2) + qJD(4);
t680 = mrSges(5,1) * t731 - mrSges(5,3) * t695;
t730 = -qJDD(2) + qJDD(4);
t718 = t741 * g(1) - t744 * g(2);
t696 = -qJDD(1) * pkin(1) - t746 * pkin(7) - t718;
t751 = -t710 * pkin(2) + t696 + (-t709 - t760) * qJ(3);
t629 = -pkin(2) * t759 + t710 * pkin(3) - pkin(8) * t773 - t751 + (t717 + t783) * t764;
t663 = -qJD(4) * t694 + t709 * t742 - t710 * t739;
t622 = t629 + (t694 * t731 - t663) * pkin(9) + (t695 * t731 - t662) * pkin(4);
t674 = pkin(4) * t694 - pkin(9) * t695;
t729 = t731 ^ 2;
t625 = -pkin(4) * t729 + pkin(9) * t730 - t674 * t694 + t627;
t738 = sin(qJ(5));
t620 = t738 * t622 + t782 * t625;
t676 = t695 * t782 + t738 * t731;
t632 = t676 * qJD(5) + t738 * t663 - t730 * t782;
t660 = qJDD(5) - t662;
t687 = qJD(5) + t694;
t666 = mrSges(6,1) * t687 - mrSges(6,3) * t676;
t675 = t738 * t695 - t731 * t782;
t650 = pkin(5) * t675 - qJ(6) * t676;
t683 = t687 ^ 2;
t616 = -pkin(5) * t683 + qJ(6) * t660 + 0.2e1 * qJD(6) * t687 - t650 * t675 + t620;
t667 = -mrSges(7,1) * t687 + mrSges(7,2) * t676;
t761 = m(7) * t616 + t660 * mrSges(7,3) + t687 * t667;
t651 = mrSges(7,1) * t675 - mrSges(7,3) * t676;
t768 = -mrSges(6,1) * t675 - mrSges(6,2) * t676 - t651;
t611 = m(6) * t620 - t660 * mrSges(6,2) + t632 * t780 - t687 * t666 + t675 * t768 + t761;
t619 = t622 * t782 - t738 * t625;
t633 = -t675 * qJD(5) + t663 * t782 + t738 * t730;
t665 = -mrSges(6,2) * t687 - mrSges(6,3) * t675;
t617 = -t660 * pkin(5) - t683 * qJ(6) + t676 * t650 + qJDD(6) - t619;
t664 = -mrSges(7,2) * t675 + mrSges(7,3) * t687;
t754 = -m(7) * t617 + t660 * mrSges(7,1) + t687 * t664;
t613 = m(6) * t619 + t660 * mrSges(6,1) + t633 * t780 + t687 * t665 + t676 * t768 + t754;
t755 = t782 * t611 - t613 * t738;
t604 = m(5) * t627 - mrSges(5,2) * t730 + mrSges(5,3) * t662 - t673 * t694 - t680 * t731 + t755;
t626 = -t739 * t637 + t742 * t638;
t679 = -mrSges(5,2) * t731 - mrSges(5,3) * t694;
t624 = -t730 * pkin(4) - t729 * pkin(9) + t695 * t674 - t626;
t618 = -0.2e1 * qJD(6) * t676 + (t675 * t687 - t633) * qJ(6) + (t676 * t687 + t632) * pkin(5) + t624;
t614 = m(7) * t618 + mrSges(7,1) * t632 - t633 * mrSges(7,3) + t664 * t675 - t676 * t667;
t748 = -m(6) * t624 - t632 * mrSges(6,1) - mrSges(6,2) * t633 - t675 * t665 - t666 * t676 - t614;
t608 = m(5) * t626 + mrSges(5,1) * t730 - mrSges(5,3) * t663 - t673 * t695 + t679 * t731 + t748;
t756 = t742 * t604 - t739 * t608;
t750 = m(4) * t657 + qJDD(2) * mrSges(4,3) + qJD(2) * t714 + t707 * t763 + t756;
t597 = m(3) * t678 - qJDD(2) * mrSges(3,2) - qJD(2) * t713 + t708 * t763 + t710 * t781 + t750;
t715 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t763;
t599 = t739 * t604 + t742 * t608;
t716 = mrSges(4,2) * t763 + qJD(2) * mrSges(4,3);
t749 = -m(4) * t661 + qJDD(2) * mrSges(4,1) + qJD(2) * t716 - t599;
t598 = m(3) * t677 + qJDD(2) * mrSges(3,1) + qJD(2) * t715 - t781 * t709 + (-t707 - t708) * t764 + t749;
t757 = t743 * t597 - t598 * t740;
t590 = m(2) * t719 - mrSges(2,1) * t746 - qJDD(1) * mrSges(2,2) + t757;
t649 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t764 + t751;
t607 = t738 * t611 + t782 * t613;
t752 = -m(5) * t629 + t662 * mrSges(5,1) - t663 * mrSges(5,2) - t694 * t679 - t695 * t680 - t607;
t602 = m(4) * t649 - mrSges(4,1) * t710 - t709 * mrSges(4,3) - t714 * t764 - t716 * t763 + t752;
t747 = -m(3) * t696 + t710 * mrSges(3,1) - mrSges(3,2) * t709 - t713 * t764 + t715 * t763 - t602;
t601 = m(2) * t718 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t746 + t747;
t772 = t741 * t590 + t744 * t601;
t591 = t740 * t597 + t743 * t598;
t771 = t786 * t675 - t778 * t676 - t774 * t687;
t770 = t774 * t675 + t788 * t676 + t784 * t687;
t769 = -t778 * t675 + t789 * t676 - t788 * t687;
t767 = t785 * qJD(2) + (t776 * t740 + t775 * t743) * qJD(1);
t766 = -t775 * qJD(2) + (-t779 * t740 - t787 * t743) * qJD(1);
t765 = t776 * qJD(2) + (t790 * t740 + t779 * t743) * qJD(1);
t758 = t744 * t590 - t601 * t741;
t670 = Ifges(5,1) * t695 - Ifges(5,4) * t694 + Ifges(5,5) * t731;
t669 = Ifges(5,4) * t695 - Ifges(5,2) * t694 + Ifges(5,6) * t731;
t668 = Ifges(5,5) * t695 - Ifges(5,6) * t694 + Ifges(5,3) * t731;
t606 = mrSges(6,2) * t624 + mrSges(7,2) * t617 - mrSges(6,3) * t619 - mrSges(7,3) * t618 - qJ(6) * t614 - t778 * t632 + t789 * t633 - t660 * t788 + t770 * t675 + t771 * t687;
t605 = -mrSges(6,1) * t624 - mrSges(7,1) * t618 + mrSges(7,2) * t616 + mrSges(6,3) * t620 - pkin(5) * t614 - t786 * t632 + t778 * t633 + t774 * t660 + t770 * t676 + t769 * t687;
t593 = Ifges(5,4) * t663 + Ifges(5,2) * t662 + Ifges(5,6) * t730 - t695 * t668 + t731 * t670 - mrSges(5,1) * t629 + mrSges(5,3) * t627 - mrSges(6,1) * t619 + mrSges(6,2) * t620 + mrSges(7,1) * t617 - mrSges(7,3) * t616 - pkin(5) * t754 - qJ(6) * t761 - pkin(4) * t607 + (pkin(5) * t651 + t771) * t676 + (qJ(6) * t651 - t769) * t675 + t784 * t660 + (mrSges(7,2) * pkin(5) + t788) * t633 + (mrSges(7,2) * qJ(6) + t774) * t632;
t592 = mrSges(5,2) * t629 - mrSges(5,3) * t626 + Ifges(5,1) * t663 + Ifges(5,4) * t662 + Ifges(5,5) * t730 - pkin(9) * t607 - t738 * t605 + t606 * t782 - t694 * t668 - t731 * t669;
t587 = mrSges(3,2) * t696 + mrSges(4,2) * t661 - mrSges(3,3) * t677 - mrSges(4,3) * t649 - pkin(8) * t599 - qJ(3) * t602 + t766 * qJD(2) + t776 * qJDD(2) + t742 * t592 - t739 * t593 + t790 * t709 + t779 * t710 + t767 * t763;
t586 = -mrSges(3,1) * t696 - mrSges(4,1) * t649 + mrSges(4,2) * t657 + mrSges(3,3) * t678 - pkin(2) * t602 - pkin(3) * t752 - pkin(8) * t756 + t765 * qJD(2) + t775 * qJDD(2) - t739 * t592 - t742 * t593 + t779 * t709 + t787 * t710 - t767 * t764;
t585 = (t765 * t743 + (pkin(2) * t707 + t766) * t740) * qJD(1) + pkin(9) * t755 - pkin(1) * t591 + t782 * t605 + pkin(3) * t599 - qJ(3) * t750 + pkin(4) * t748 - pkin(2) * t749 + Ifges(2,6) * qJDD(1) - t785 * qJDD(2) + mrSges(2,1) * g(3) + (-mrSges(4,2) * qJ(3) - t775) * t710 + (mrSges(4,2) * pkin(2) - t776) * t709 + t746 * Ifges(2,5) + t738 * t606 + Ifges(5,3) * t730 + mrSges(2,3) * t719 + t694 * t670 + t695 * t669 - mrSges(3,1) * t677 + mrSges(3,2) * t678 + mrSges(4,1) * t661 + Ifges(5,6) * t662 + Ifges(5,5) * t663 - mrSges(4,3) * t657 + mrSges(5,1) * t626 - mrSges(5,2) * t627;
t584 = -mrSges(2,2) * g(3) - mrSges(2,3) * t718 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t746 - pkin(7) * t591 - t586 * t740 + t587 * t743;
t1 = [-m(1) * g(1) + t758; -m(1) * g(2) + t772; (-m(1) - m(2)) * g(3) + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t772 + t744 * t584 - t741 * t585; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t758 + t741 * t584 + t744 * t585; -mrSges(1,1) * g(2) + mrSges(2,1) * t718 + mrSges(1,2) * g(1) - mrSges(2,2) * t719 + Ifges(2,3) * qJDD(1) + pkin(1) * t747 + pkin(7) * t757 + t743 * t586 + t740 * t587;];
tauB  = t1;
