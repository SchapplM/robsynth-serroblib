% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:16
% EndTime: 2019-12-31 19:10:23
% DurationCPUTime: 6.84s
% Computational Cost: add. (95789->293), mult. (226044->368), div. (0->0), fcn. (165474->10), ass. (0->130)
t753 = qJD(1) ^ 2;
t743 = cos(pkin(9));
t781 = pkin(2) * t743;
t742 = sin(pkin(9));
t780 = mrSges(3,2) * t742;
t740 = t743 ^ 2;
t779 = t740 * t753;
t747 = sin(qJ(1));
t751 = cos(qJ(1));
t730 = -t751 * g(1) - t747 * g(2);
t725 = -t753 * pkin(1) + qJDD(1) * qJ(2) + t730;
t774 = qJD(1) * qJD(2);
t771 = -t743 * g(3) - 0.2e1 * t742 * t774;
t696 = (-pkin(6) * qJDD(1) + t753 * t781 - t725) * t742 + t771;
t712 = -t742 * g(3) + (t725 + 0.2e1 * t774) * t743;
t772 = qJDD(1) * t743;
t697 = -pkin(2) * t779 + pkin(6) * t772 + t712;
t746 = sin(qJ(3));
t750 = cos(qJ(3));
t676 = t746 * t696 + t750 * t697;
t775 = t743 * qJD(1);
t776 = t742 * qJD(1);
t723 = -t746 * t776 + t750 * t775;
t762 = t742 * t750 + t743 * t746;
t724 = t762 * qJD(1);
t707 = -t723 * pkin(3) - t724 * pkin(7);
t752 = qJD(3) ^ 2;
t666 = -t752 * pkin(3) + qJDD(3) * pkin(7) + t723 * t707 + t676;
t739 = t742 ^ 2;
t729 = t747 * g(1) - t751 * g(2);
t766 = qJDD(2) - t729;
t708 = (-pkin(1) - t781) * qJDD(1) + (-qJ(2) + (-t739 - t740) * pkin(6)) * t753 + t766;
t720 = t724 * qJD(3);
t773 = qJDD(1) * t742;
t709 = -t746 * t773 + t750 * t772 - t720;
t777 = t723 * qJD(3);
t710 = t762 * qJDD(1) + t777;
t669 = (-t710 - t777) * pkin(7) + (-t709 + t720) * pkin(3) + t708;
t745 = sin(qJ(4));
t749 = cos(qJ(4));
t656 = -t745 * t666 + t749 * t669;
t714 = t749 * qJD(3) - t745 * t724;
t685 = t714 * qJD(4) + t745 * qJDD(3) + t749 * t710;
t706 = qJDD(4) - t709;
t715 = t745 * qJD(3) + t749 * t724;
t721 = qJD(4) - t723;
t653 = (t714 * t721 - t685) * pkin(8) + (t714 * t715 + t706) * pkin(4) + t656;
t657 = t749 * t666 + t745 * t669;
t684 = -t715 * qJD(4) + t749 * qJDD(3) - t745 * t710;
t695 = t721 * pkin(4) - t715 * pkin(8);
t713 = t714 ^ 2;
t654 = -t713 * pkin(4) + t684 * pkin(8) - t721 * t695 + t657;
t744 = sin(qJ(5));
t748 = cos(qJ(5));
t651 = t748 * t653 - t744 * t654;
t686 = t748 * t714 - t744 * t715;
t662 = t686 * qJD(5) + t744 * t684 + t748 * t685;
t687 = t744 * t714 + t748 * t715;
t674 = -t686 * mrSges(6,1) + t687 * mrSges(6,2);
t719 = qJD(5) + t721;
t677 = -t719 * mrSges(6,2) + t686 * mrSges(6,3);
t703 = qJDD(5) + t706;
t647 = m(6) * t651 + t703 * mrSges(6,1) - t662 * mrSges(6,3) - t687 * t674 + t719 * t677;
t652 = t744 * t653 + t748 * t654;
t661 = -t687 * qJD(5) + t748 * t684 - t744 * t685;
t678 = t719 * mrSges(6,1) - t687 * mrSges(6,3);
t648 = m(6) * t652 - t703 * mrSges(6,2) + t661 * mrSges(6,3) + t686 * t674 - t719 * t678;
t639 = t748 * t647 + t744 * t648;
t688 = -t714 * mrSges(5,1) + t715 * mrSges(5,2);
t691 = -t721 * mrSges(5,2) + t714 * mrSges(5,3);
t637 = m(5) * t656 + t706 * mrSges(5,1) - t685 * mrSges(5,3) - t715 * t688 + t721 * t691 + t639;
t692 = t721 * mrSges(5,1) - t715 * mrSges(5,3);
t767 = -t744 * t647 + t748 * t648;
t638 = m(5) * t657 - t706 * mrSges(5,2) + t684 * mrSges(5,3) + t714 * t688 - t721 * t692 + t767;
t633 = -t745 * t637 + t749 * t638;
t704 = -t723 * mrSges(4,1) + t724 * mrSges(4,2);
t717 = qJD(3) * mrSges(4,1) - t724 * mrSges(4,3);
t631 = m(4) * t676 - qJDD(3) * mrSges(4,2) + t709 * mrSges(4,3) - qJD(3) * t717 + t723 * t704 + t633;
t675 = t750 * t696 - t746 * t697;
t665 = -qJDD(3) * pkin(3) - t752 * pkin(7) + t724 * t707 - t675;
t655 = -t684 * pkin(4) - t713 * pkin(8) + t715 * t695 + t665;
t759 = m(6) * t655 - t661 * mrSges(6,1) + t662 * mrSges(6,2) - t686 * t677 + t687 * t678;
t649 = -m(5) * t665 + t684 * mrSges(5,1) - t685 * mrSges(5,2) + t714 * t691 - t715 * t692 - t759;
t716 = -qJD(3) * mrSges(4,2) + t723 * mrSges(4,3);
t643 = m(4) * t675 + qJDD(3) * mrSges(4,1) - t710 * mrSges(4,3) + qJD(3) * t716 - t724 * t704 + t649;
t623 = t746 * t631 + t750 * t643;
t711 = -t742 * t725 + t771;
t761 = mrSges(3,3) * qJDD(1) + t753 * (-mrSges(3,1) * t743 + t780);
t621 = m(3) * t711 - t761 * t742 + t623;
t768 = t750 * t631 - t746 * t643;
t622 = m(3) * t712 + t761 * t743 + t768;
t769 = -t742 * t621 + t743 * t622;
t613 = m(2) * t730 - t753 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t769;
t722 = -qJDD(1) * pkin(1) - t753 * qJ(2) + t766;
t632 = t749 * t637 + t745 * t638;
t757 = m(4) * t708 - t709 * mrSges(4,1) + t710 * mrSges(4,2) - t723 * t716 + t724 * t717 + t632;
t756 = -m(3) * t722 + mrSges(3,1) * t772 - t757 + (t739 * t753 + t779) * mrSges(3,3);
t626 = (mrSges(2,1) - t780) * qJDD(1) - t753 * mrSges(2,2) + m(2) * t729 + t756;
t778 = t747 * t613 + t751 * t626;
t615 = t743 * t621 + t742 * t622;
t770 = t751 * t613 - t747 * t626;
t765 = Ifges(3,1) * t742 + Ifges(3,4) * t743;
t764 = Ifges(3,4) * t742 + Ifges(3,2) * t743;
t763 = Ifges(3,5) * t742 + Ifges(3,6) * t743;
t670 = Ifges(6,5) * t687 + Ifges(6,6) * t686 + Ifges(6,3) * t719;
t672 = Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t719;
t640 = -mrSges(6,1) * t655 + mrSges(6,3) * t652 + Ifges(6,4) * t662 + Ifges(6,2) * t661 + Ifges(6,6) * t703 - t687 * t670 + t719 * t672;
t671 = Ifges(6,4) * t687 + Ifges(6,2) * t686 + Ifges(6,6) * t719;
t641 = mrSges(6,2) * t655 - mrSges(6,3) * t651 + Ifges(6,1) * t662 + Ifges(6,4) * t661 + Ifges(6,5) * t703 + t686 * t670 - t719 * t671;
t679 = Ifges(5,5) * t715 + Ifges(5,6) * t714 + Ifges(5,3) * t721;
t681 = Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * t721;
t620 = -mrSges(5,1) * t665 + mrSges(5,3) * t657 + Ifges(5,4) * t685 + Ifges(5,2) * t684 + Ifges(5,6) * t706 - pkin(4) * t759 + pkin(8) * t767 + t748 * t640 + t744 * t641 - t715 * t679 + t721 * t681;
t680 = Ifges(5,4) * t715 + Ifges(5,2) * t714 + Ifges(5,6) * t721;
t624 = mrSges(5,2) * t665 - mrSges(5,3) * t656 + Ifges(5,1) * t685 + Ifges(5,4) * t684 + Ifges(5,5) * t706 - pkin(8) * t639 - t744 * t640 + t748 * t641 + t714 * t679 - t721 * t680;
t698 = Ifges(4,5) * t724 + Ifges(4,6) * t723 + Ifges(4,3) * qJD(3);
t699 = Ifges(4,4) * t724 + Ifges(4,2) * t723 + Ifges(4,6) * qJD(3);
t610 = mrSges(4,2) * t708 - mrSges(4,3) * t675 + Ifges(4,1) * t710 + Ifges(4,4) * t709 + Ifges(4,5) * qJDD(3) - pkin(7) * t632 - qJD(3) * t699 - t745 * t620 + t749 * t624 + t723 * t698;
t700 = Ifges(4,1) * t724 + Ifges(4,4) * t723 + Ifges(4,5) * qJD(3);
t758 = -mrSges(6,1) * t651 + mrSges(6,2) * t652 - Ifges(6,5) * t662 - Ifges(6,6) * t661 - Ifges(6,3) * t703 - t687 * t671 + t686 * t672;
t754 = mrSges(5,1) * t656 - mrSges(5,2) * t657 + Ifges(5,5) * t685 + Ifges(5,6) * t684 + Ifges(5,3) * t706 + pkin(4) * t639 + t715 * t680 - t714 * t681 - t758;
t616 = -mrSges(4,1) * t708 + mrSges(4,3) * t676 + Ifges(4,4) * t710 + Ifges(4,2) * t709 + Ifges(4,6) * qJDD(3) - pkin(3) * t632 + qJD(3) * t700 - t724 * t698 - t754;
t727 = t763 * qJD(1);
t606 = -mrSges(3,1) * t722 + mrSges(3,3) * t712 - pkin(2) * t757 + pkin(6) * t768 + t764 * qJDD(1) + t746 * t610 + t750 * t616 - t727 * t776;
t609 = mrSges(3,2) * t722 - mrSges(3,3) * t711 - pkin(6) * t623 + t765 * qJDD(1) + t750 * t610 - t746 * t616 + t727 * t775;
t628 = mrSges(3,2) * t773 - t756;
t760 = mrSges(2,1) * t729 - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) - pkin(1) * t628 + qJ(2) * t769 + t743 * t606 + t742 * t609;
t755 = mrSges(4,1) * t675 - mrSges(4,2) * t676 + Ifges(4,5) * t710 + Ifges(4,6) * t709 + Ifges(4,3) * qJDD(3) + pkin(3) * t649 + pkin(7) * t633 + t749 * t620 + t745 * t624 + t724 * t699 - t723 * t700;
t607 = -pkin(1) * t615 + mrSges(2,1) * g(3) + (Ifges(2,6) - t763) * qJDD(1) + mrSges(2,3) * t730 - mrSges(3,1) * t711 + mrSges(3,2) * t712 - t755 - pkin(2) * t623 + (-t742 * t764 + t743 * t765 + Ifges(2,5)) * t753;
t604 = -mrSges(2,2) * g(3) - mrSges(2,3) * t729 + Ifges(2,5) * qJDD(1) - t753 * Ifges(2,6) - qJ(2) * t615 - t742 * t606 + t743 * t609;
t1 = [-m(1) * g(1) + t770; -m(1) * g(2) + t778; (-m(1) - m(2)) * g(3) + t615; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t778 + t751 * t604 - t747 * t607; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t770 + t747 * t604 + t751 * t607; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t760; t760; t628; t755; t754; -t758;];
tauJB = t1;
