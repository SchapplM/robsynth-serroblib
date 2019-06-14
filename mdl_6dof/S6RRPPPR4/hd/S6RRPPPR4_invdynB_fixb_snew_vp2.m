% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-05-06 08:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:41:59
% EndTime: 2019-05-06 08:42:05
% DurationCPUTime: 4.34s
% Computational Cost: add. (39901->350), mult. (88140->406), div. (0->0), fcn. (49938->8), ass. (0->139)
t801 = -2 * qJD(4);
t800 = Ifges(5,1) + Ifges(6,1);
t790 = Ifges(3,4) + Ifges(4,6);
t789 = Ifges(5,4) - Ifges(6,5);
t788 = Ifges(3,5) - Ifges(4,4);
t787 = Ifges(5,5) + Ifges(6,4);
t799 = Ifges(3,2) + Ifges(4,3);
t798 = Ifges(4,2) + Ifges(3,1);
t797 = Ifges(5,2) + Ifges(6,3);
t796 = -Ifges(6,2) - Ifges(5,3);
t786 = Ifges(3,6) - Ifges(4,5);
t785 = Ifges(5,6) - Ifges(6,6);
t795 = Ifges(3,3) + Ifges(4,1);
t746 = cos(qJ(2));
t743 = sin(qJ(2));
t771 = qJD(1) * qJD(2);
t768 = t743 * t771;
t714 = t746 * qJDD(1) - t768;
t739 = sin(pkin(9));
t740 = cos(pkin(9));
t680 = t740 * qJDD(2) - t739 * t714;
t774 = qJD(1) * t746;
t702 = t739 * qJD(2) + t740 * t774;
t772 = t743 * qJD(1);
t769 = t702 * t772;
t794 = (-t680 + t769) * qJ(5);
t719 = pkin(3) * t772 - qJD(2) * qJ(4);
t738 = t746 ^ 2;
t749 = qJD(1) ^ 2;
t767 = t746 * t771;
t713 = t743 * qJDD(1) + t767;
t744 = sin(qJ(1));
t747 = cos(qJ(1));
t722 = t744 * g(1) - t747 * g(2);
t760 = -qJDD(1) * pkin(1) - t722;
t753 = pkin(2) * t768 - 0.2e1 * qJD(3) * t772 + (-t713 - t767) * qJ(3) + t760;
t636 = -t719 * t772 + (-pkin(3) * t738 - pkin(7)) * t749 + (-pkin(2) - qJ(4)) * t714 + t753;
t723 = -t747 * g(1) - t744 * g(2);
t697 = -t749 * pkin(1) + qJDD(1) * pkin(7) + t723;
t671 = -t746 * g(3) - t743 * t697;
t710 = (-pkin(2) * t746 - qJ(3) * t743) * qJD(1);
t748 = qJD(2) ^ 2;
t651 = -qJDD(2) * pkin(2) - t748 * qJ(3) + t710 * t772 + qJDD(3) - t671;
t646 = (-t743 * t746 * t749 - qJDD(2)) * qJ(4) + (t713 - t767) * pkin(3) + t651;
t703 = t740 * qJD(2) - t739 * t774;
t630 = -t739 * t636 + t740 * t646 + t703 * t801;
t793 = 2 * qJD(5);
t792 = t749 * pkin(7);
t791 = -mrSges(5,3) - mrSges(6,2);
t784 = t743 ^ 2 * t749;
t711 = (mrSges(4,2) * t746 - mrSges(4,3) * t743) * qJD(1);
t712 = (-mrSges(3,1) * t746 + mrSges(3,2) * t743) * qJD(1);
t718 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t774;
t720 = -mrSges(4,1) * t774 - qJD(2) * mrSges(4,3);
t631 = t740 * t636 + t739 * t646 + t702 * t801;
t677 = mrSges(5,1) * t772 - t703 * mrSges(5,3);
t679 = t739 * qJDD(2) + t740 * t714;
t666 = t702 * pkin(4) - t703 * qJ(5);
t627 = -pkin(4) * t784 + t713 * qJ(5) - t702 * t666 + t772 * t793 + t631;
t678 = -mrSges(6,1) * t772 + t703 * mrSges(6,2);
t628 = -t713 * pkin(4) - qJ(5) * t784 + t703 * t666 + qJDD(5) - t630;
t624 = (-t680 - t769) * pkin(8) + (t702 * t703 - t713) * pkin(5) + t628;
t681 = -pkin(5) * t772 - t703 * pkin(8);
t700 = t702 ^ 2;
t625 = -t700 * pkin(5) + t679 * pkin(8) + t681 * t772 + t627;
t742 = sin(qJ(6));
t745 = cos(qJ(6));
t622 = t745 * t624 - t742 * t625;
t664 = t745 * t702 - t742 * t703;
t639 = t664 * qJD(6) + t742 * t679 + t745 * t680;
t665 = t742 * t702 + t745 * t703;
t648 = -t664 * mrSges(7,1) + t665 * mrSges(7,2);
t726 = qJD(6) - t772;
t652 = -t726 * mrSges(7,2) + t664 * mrSges(7,3);
t709 = qJDD(6) - t713;
t619 = m(7) * t622 + t709 * mrSges(7,1) - t639 * mrSges(7,3) - t665 * t648 + t726 * t652;
t623 = t742 * t624 + t745 * t625;
t638 = -t665 * qJD(6) + t745 * t679 - t742 * t680;
t653 = t726 * mrSges(7,1) - t665 * mrSges(7,3);
t620 = m(7) * t623 - t709 * mrSges(7,2) + t638 * mrSges(7,3) + t664 * t648 - t726 * t653;
t764 = -t742 * t619 + t745 * t620;
t758 = m(6) * t627 + t713 * mrSges(6,3) + t678 * t772 + t764;
t667 = t702 * mrSges(6,1) - t703 * mrSges(6,3);
t778 = -t702 * mrSges(5,1) - t703 * mrSges(5,2) - t667;
t609 = m(5) * t631 - t713 * mrSges(5,2) - t677 * t772 + t791 * t679 + t778 * t702 + t758;
t676 = -mrSges(5,2) * t772 - t702 * mrSges(5,3);
t612 = t745 * t619 + t742 * t620;
t675 = -t702 * mrSges(6,2) + mrSges(6,3) * t772;
t755 = -m(6) * t628 + t713 * mrSges(6,1) + t675 * t772 - t612;
t611 = m(5) * t630 + t713 * mrSges(5,1) + t676 * t772 + t791 * t680 + t778 * t703 + t755;
t607 = t739 * t609 + t740 * t611;
t757 = -m(4) * t651 - t713 * mrSges(4,1) - t607;
t605 = m(3) * t671 - t713 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t718 - t720) * qJD(2) + (-t711 - t712) * t772 + t757;
t672 = -t743 * g(3) + t746 * t697;
t717 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t772;
t770 = qJD(3) * qJD(2);
t731 = -0.2e1 * t770;
t761 = t748 * pkin(2) - qJDD(2) * qJ(3) - t710 * t774 - t672;
t650 = t731 + t761;
t721 = mrSges(4,1) * t772 + qJD(2) * mrSges(4,2);
t754 = t738 * t749 * qJ(4) - t714 * pkin(3) - qJD(2) * t719 - qJDD(4) + t761;
t642 = 0.2e1 * t770 - t754;
t633 = -0.2e1 * qJD(5) * t703 + t794 + (t703 * t772 + t679) * pkin(4) + t642;
t629 = -t700 * pkin(8) + t731 + (-pkin(4) - pkin(5)) * t679 - t794 + (-pkin(4) * t772 + t681 + t793) * t703 + t754;
t756 = -m(7) * t629 + t638 * mrSges(7,1) - t639 * mrSges(7,2) + t664 * t652 - t665 * t653;
t621 = m(6) * t633 + t679 * mrSges(6,1) - t680 * mrSges(6,3) + t702 * t675 - t703 * t678 + t756;
t751 = -m(5) * t642 - t679 * mrSges(5,1) - t680 * mrSges(5,2) - t702 * t676 - t703 * t677 - t621;
t750 = -m(4) * t650 + qJDD(2) * mrSges(4,3) + qJD(2) * t721 + t711 * t774 - t751;
t617 = t712 * t774 + t750 + (mrSges(3,3) + mrSges(4,1)) * t714 - qJD(2) * t717 + m(3) * t672 - qJDD(2) * mrSges(3,2);
t765 = -t743 * t605 + t746 * t617;
t599 = m(2) * t723 - t749 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t765;
t696 = t760 - t792;
t649 = -t714 * pkin(2) + t753 - t792;
t782 = t740 * t609 - t739 * t611;
t759 = -m(4) * t649 - t714 * mrSges(4,2) + t721 * t772 - t782;
t752 = -m(3) * t696 + t718 * t774 + t714 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t713 + (-t717 * t743 - t720 * t746) * qJD(1) + t759;
t603 = m(2) * t722 + qJDD(1) * mrSges(2,1) - t749 * mrSges(2,2) + t752;
t783 = t744 * t599 + t747 * t603;
t600 = t746 * t605 + t743 * t617;
t781 = t797 * t702 - t789 * t703 - t785 * t772;
t780 = t785 * t702 - t787 * t703 + t796 * t772;
t779 = -t789 * t702 + t800 * t703 + t787 * t772;
t777 = t795 * qJD(2) + (t788 * t743 + t786 * t746) * qJD(1);
t776 = -t786 * qJD(2) + (-t790 * t743 - t799 * t746) * qJD(1);
t775 = t788 * qJD(2) + (t798 * t743 + t790 * t746) * qJD(1);
t766 = t747 * t599 - t744 * t603;
t645 = Ifges(7,1) * t665 + Ifges(7,4) * t664 + Ifges(7,5) * t726;
t644 = Ifges(7,4) * t665 + Ifges(7,2) * t664 + Ifges(7,6) * t726;
t643 = Ifges(7,5) * t665 + Ifges(7,6) * t664 + Ifges(7,3) * t726;
t614 = mrSges(7,2) * t629 - mrSges(7,3) * t622 + Ifges(7,1) * t639 + Ifges(7,4) * t638 + Ifges(7,5) * t709 + t664 * t643 - t726 * t644;
t613 = -mrSges(7,1) * t629 + mrSges(7,3) * t623 + Ifges(7,4) * t639 + Ifges(7,2) * t638 + Ifges(7,6) * t709 - t665 * t643 + t726 * t645;
t606 = -t713 * mrSges(4,3) + t720 * t774 - t759;
t601 = mrSges(5,2) * t642 + mrSges(6,2) * t628 - mrSges(5,3) * t630 - mrSges(6,3) * t633 - pkin(8) * t612 - qJ(5) * t621 - t742 * t613 + t745 * t614 - t789 * t679 + t800 * t680 + t780 * t702 + t787 * t713 + t781 * t772;
t596 = -mrSges(5,1) * t642 - mrSges(6,1) * t633 + mrSges(6,2) * t627 + mrSges(5,3) * t631 - pkin(4) * t621 - pkin(5) * t756 - pkin(8) * t764 - t745 * t613 - t742 * t614 - t797 * t679 + t789 * t680 + t780 * t703 + t785 * t713 + t779 * t772;
t595 = (-qJ(5) * mrSges(6,2) - t785) * t679 + (-pkin(4) * mrSges(6,2) + t787) * t680 + t788 * qJDD(2) + t790 * t714 + t776 * qJD(2) + t777 * t774 + (-qJ(5) * t667 + t779) * t702 + (-pkin(4) * t667 - t781) * t703 + qJ(5) * t758 + (-t796 + t798) * t713 - Ifges(7,3) * t709 + mrSges(3,2) * t696 - mrSges(3,3) * t671 + t664 * t645 - t665 * t644 - mrSges(4,3) * t649 + mrSges(4,1) * t651 - Ifges(7,6) * t638 - Ifges(7,5) * t639 + mrSges(6,3) * t627 - mrSges(6,1) * t628 + mrSges(5,1) * t630 - mrSges(5,2) * t631 + mrSges(7,2) * t623 - mrSges(7,1) * t622 - pkin(5) * t612 + pkin(3) * t607 - qJ(3) * t606 + pkin(4) * t755;
t594 = -mrSges(3,1) * t696 - mrSges(4,1) * t650 + mrSges(4,2) * t649 + mrSges(3,3) * t672 - pkin(2) * t606 - pkin(3) * t751 - qJ(4) * t782 + t775 * qJD(2) + t786 * qJDD(2) - t740 * t596 - t739 * t601 + t790 * t713 + t799 * t714 - t777 * t772;
t593 = -pkin(1) * t600 + mrSges(2,3) * t723 - pkin(2) * (-qJD(2) * t720 + t757) - qJ(3) * t750 + t739 * t596 + qJ(4) * t607 - mrSges(3,1) * t671 + mrSges(3,2) * t672 - mrSges(4,2) * t651 + mrSges(4,3) * t650 - t740 * t601 + mrSges(2,1) * g(3) + t749 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t786) * t714 - t788 * t713 + (pkin(2) * mrSges(4,2) - t795) * qJDD(2) + (t775 * t746 + (pkin(2) * t711 + t776) * t743) * qJD(1);
t592 = -mrSges(2,2) * g(3) - mrSges(2,3) * t722 + Ifges(2,5) * qJDD(1) - t749 * Ifges(2,6) - pkin(7) * t600 - t743 * t594 + t746 * t595;
t1 = [-m(1) * g(1) + t766; -m(1) * g(2) + t783; (-m(1) - m(2)) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t783 + t747 * t592 - t744 * t593; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t766 + t744 * t592 + t747 * t593; -mrSges(1,1) * g(2) + mrSges(2,1) * t722 + mrSges(1,2) * g(1) - mrSges(2,2) * t723 + Ifges(2,3) * qJDD(1) + pkin(1) * t752 + pkin(7) * t765 + t746 * t594 + t743 * t595;];
tauB  = t1;
