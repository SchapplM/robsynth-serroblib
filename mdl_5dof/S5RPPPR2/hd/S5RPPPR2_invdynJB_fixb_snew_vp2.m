% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:48
% EndTime: 2020-01-03 11:22:53
% DurationCPUTime: 5.67s
% Computational Cost: add. (44592->280), mult. (124905->382), div. (0->0), fcn. (84227->10), ass. (0->133)
t754 = sin(pkin(7));
t757 = cos(pkin(7));
t759 = sin(qJ(1));
t761 = cos(qJ(1));
t735 = -t759 * g(2) + t761 * g(3);
t762 = qJD(1) ^ 2;
t809 = -t762 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t735;
t706 = -t754 * g(1) + t809 * t757;
t775 = -pkin(2) * t757 - qJ(3) * t754;
t729 = t775 * qJD(1);
t795 = t757 * qJD(1);
t697 = t729 * t795 + t706;
t753 = sin(pkin(8));
t756 = cos(pkin(8));
t736 = -t761 * g(2) - t759 * g(3);
t770 = -t762 * qJ(2) + qJDD(2) - t736;
t796 = t754 * qJD(1);
t808 = (-pkin(1) + t775) * qJDD(1) + t770 - 0.2e1 * qJD(3) * t796;
t677 = -t753 * t697 + t808 * t756;
t705 = -t757 * g(1) - t809 * t754;
t752 = sin(pkin(9));
t755 = cos(pkin(9));
t800 = t754 * t756;
t769 = t752 * t800 + t755 * t757;
t718 = t769 * qJD(1);
t716 = t769 * qJDD(1);
t807 = 2 * qJD(4);
t806 = mrSges(3,2) * t754;
t805 = Ifges(4,4) * t756;
t804 = Ifges(4,6) * t757;
t803 = t754 ^ 2 * t762;
t802 = t757 ^ 2 * t762;
t801 = t753 * t754;
t799 = t757 * t762;
t730 = (-mrSges(3,1) * t757 + t806) * qJD(1);
t678 = t756 * t697 + t808 * t753;
t779 = mrSges(4,1) * t753 + mrSges(4,2) * t756;
t723 = t779 * t796;
t772 = -mrSges(4,1) * t757 - mrSges(4,3) * t800;
t727 = t772 * qJD(1);
t771 = mrSges(4,2) * t757 - mrSges(4,3) * t801;
t722 = (pkin(3) * t753 - qJ(4) * t756) * t796;
t790 = t753 * t796;
t792 = qJDD(1) * t757;
t675 = -pkin(3) * t802 - qJ(4) * t792 - t722 * t790 + t678;
t696 = t729 * t796 + qJDD(3) - t705;
t683 = ((-qJDD(1) * t756 - t753 * t799) * qJ(4) + (qJDD(1) * t753 - t756 * t799) * pkin(3)) * t754 + t696;
t671 = t755 * t675 + t752 * t683 - t718 * t807;
t789 = t756 * t796;
t719 = -t752 * t795 + t755 * t789;
t698 = t718 * mrSges(5,1) + t719 * mrSges(5,2);
t704 = mrSges(5,1) * t790 - t719 * mrSges(5,3);
t699 = t718 * pkin(4) - t719 * pkin(6);
t793 = qJDD(1) * t754;
t786 = t753 * t793;
t791 = t753 ^ 2 * t803;
t669 = -pkin(4) * t791 + pkin(6) * t786 - t718 * t699 + t671;
t674 = pkin(3) * t792 - qJ(4) * t802 + t722 * t789 + qJDD(4) - t677;
t717 = (-t752 * t757 + t755 * t800) * qJDD(1);
t672 = (t718 * t790 - t717) * pkin(6) + (t719 * t790 + t716) * pkin(4) + t674;
t758 = sin(qJ(5));
t760 = cos(qJ(5));
t666 = -t758 * t669 + t760 * t672;
t700 = -t758 * t719 + t760 * t790;
t701 = t760 * t719 + t758 * t790;
t685 = -t700 * mrSges(6,1) + t701 * mrSges(6,2);
t687 = t700 * qJD(5) + t760 * t717 + t758 * t786;
t715 = qJD(5) + t718;
t688 = -t715 * mrSges(6,2) + t700 * mrSges(6,3);
t714 = qJDD(5) + t716;
t664 = m(6) * t666 + t714 * mrSges(6,1) - t687 * mrSges(6,3) - t701 * t685 + t715 * t688;
t667 = t760 * t669 + t758 * t672;
t686 = -t701 * qJD(5) - t758 * t717 + t760 * t786;
t689 = t715 * mrSges(6,1) - t701 * mrSges(6,3);
t665 = m(6) * t667 - t714 * mrSges(6,2) + t686 * mrSges(6,3) + t700 * t685 - t715 * t689;
t781 = -t758 * t664 + t760 * t665;
t657 = m(5) * t671 - t716 * mrSges(5,3) - t718 * t698 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t704) * t801 + t781;
t774 = t752 * t675 - t755 * t683;
t670 = -0.2e1 * qJD(4) * t719 - t774;
t703 = -mrSges(5,2) * t790 - t718 * mrSges(5,3);
t668 = -pkin(4) * t786 - pkin(6) * t791 + (t807 + t699) * t719 + t774;
t766 = -m(6) * t668 + t686 * mrSges(6,1) - t687 * mrSges(6,2) + t700 * t688 - t701 * t689;
t662 = m(5) * t670 - t717 * mrSges(5,3) - t719 * t698 + (mrSges(5,1) * qJDD(1) + qJD(1) * t703) * t801 + t766;
t782 = t755 * t657 - t752 * t662;
t652 = m(4) * t678 + t771 * qJDD(1) + (-t723 * t801 + t727 * t757) * qJD(1) + t782;
t659 = t760 * t664 + t758 * t665;
t658 = m(5) * t674 + t716 * mrSges(5,1) + t717 * mrSges(5,2) + t718 * t703 + t719 * t704 + t659;
t726 = t771 * qJD(1);
t655 = m(4) * t677 + t772 * qJDD(1) + (-t723 * t800 - t726 * t757) * qJD(1) - t658;
t783 = t756 * t652 - t753 * t655;
t643 = m(3) * t706 + (qJDD(1) * mrSges(3,3) + qJD(1) * t730) * t757 + t783;
t654 = t752 * t657 + t755 * t662;
t768 = m(4) * t696 + t654;
t773 = t726 * t753 + t727 * t756;
t650 = m(3) * t705 + ((-mrSges(3,3) - t779) * qJDD(1) + (-t730 - t773) * qJD(1)) * t754 - t768;
t784 = t757 * t643 - t754 * t650;
t635 = m(2) * t735 - t762 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t784;
t646 = t753 * t652 + t756 * t655;
t725 = -qJDD(1) * pkin(1) + t770;
t764 = -m(3) * t725 + mrSges(3,1) * t792 - t646 + (t802 + t803) * mrSges(3,3);
t640 = m(2) * t736 - t762 * mrSges(2,2) + (mrSges(2,1) - t806) * qJDD(1) + t764;
t798 = t759 * t635 + t761 * t640;
t637 = t754 * t643 + t757 * t650;
t785 = -t761 * t635 + t759 * t640;
t778 = Ifges(3,1) * t754 + Ifges(3,4) * t757;
t777 = Ifges(3,5) * t754 + Ifges(3,6) * t757;
t776 = Ifges(4,5) * t756 - Ifges(4,6) * t753;
t679 = Ifges(6,5) * t701 + Ifges(6,6) * t700 + Ifges(6,3) * t715;
t681 = Ifges(6,1) * t701 + Ifges(6,4) * t700 + Ifges(6,5) * t715;
t660 = -mrSges(6,1) * t668 + mrSges(6,3) * t667 + Ifges(6,4) * t687 + Ifges(6,2) * t686 + Ifges(6,6) * t714 - t701 * t679 + t715 * t681;
t680 = Ifges(6,4) * t701 + Ifges(6,2) * t700 + Ifges(6,6) * t715;
t661 = mrSges(6,2) * t668 - mrSges(6,3) * t666 + Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t714 + t700 * t679 - t715 * t680;
t690 = Ifges(5,5) * t719 - Ifges(5,6) * t718 + Ifges(5,3) * t790;
t691 = Ifges(5,4) * t719 - Ifges(5,2) * t718 + Ifges(5,6) * t790;
t647 = mrSges(5,2) * t674 - mrSges(5,3) * t670 + Ifges(5,1) * t717 - Ifges(5,4) * t716 - pkin(6) * t659 - t758 * t660 + t760 * t661 - t718 * t690 + (Ifges(5,5) * qJDD(1) - qJD(1) * t691) * t801;
t692 = Ifges(5,1) * t719 - Ifges(5,4) * t718 + Ifges(5,5) * t790;
t763 = mrSges(6,1) * t666 - mrSges(6,2) * t667 + Ifges(6,5) * t687 + Ifges(6,6) * t686 + Ifges(6,3) * t714 + t701 * t680 - t700 * t681;
t648 = -mrSges(5,1) * t674 + mrSges(5,3) * t671 + Ifges(5,4) * t717 - Ifges(5,2) * t716 - pkin(4) * t659 - t719 * t690 + (Ifges(5,6) * qJDD(1) + qJD(1) * t692) * t801 - t763;
t709 = (-Ifges(4,3) * t757 + t776 * t754) * qJD(1);
t710 = (-t804 + (-Ifges(4,2) * t753 + t805) * t754) * qJD(1);
t765 = -Ifges(4,5) * t757 + (Ifges(4,1) * t756 - Ifges(4,4) * t753) * t754;
t632 = mrSges(4,2) * t696 - mrSges(4,3) * t677 - qJ(4) * t654 + t755 * t647 - t752 * t648 + (-t709 * t801 + t710 * t757) * qJD(1) + t765 * qJDD(1);
t711 = t765 * qJD(1);
t633 = -mrSges(4,1) * t696 + mrSges(4,3) * t678 - Ifges(5,5) * t717 + Ifges(5,6) * t716 - t719 * t691 - t718 * t692 - mrSges(5,1) * t670 + mrSges(5,2) * t671 - t758 * t661 - t760 * t660 - pkin(4) * t766 - pkin(6) * t781 - pkin(3) * t654 + (-t709 * t800 - t757 * t711) * qJD(1) + (-t804 + (t805 + (-Ifges(4,2) - Ifges(5,3)) * t753) * t754) * qJDD(1);
t731 = t777 * qJD(1);
t629 = mrSges(3,2) * t725 - mrSges(3,3) * t705 - qJ(3) * t646 + t778 * qJDD(1) + t756 * t632 - t753 * t633 + t731 * t795;
t631 = -mrSges(3,1) * t725 + mrSges(3,3) * t706 - mrSges(4,1) * t677 + mrSges(4,2) * t678 - t752 * t647 - t755 * t648 + pkin(3) * t658 - qJ(4) * t782 - pkin(2) * t646 + (Ifges(3,2) + Ifges(4,3)) * t792 + ((Ifges(3,4) - t776) * qJDD(1) + (-t710 * t756 - t711 * t753 - t731) * qJD(1)) * t754;
t645 = mrSges(3,2) * t793 - t764;
t767 = mrSges(2,1) * t736 - mrSges(2,2) * t735 + Ifges(2,3) * qJDD(1) - pkin(1) * t645 + qJ(2) * t784 + t754 * t629 + t757 * t631;
t653 = (t773 * qJD(1) + t779 * qJDD(1)) * t754 + t768;
t627 = mrSges(2,1) * g(1) + mrSges(2,3) * t735 - mrSges(3,1) * t705 + mrSges(3,2) * t706 - t753 * t632 - t756 * t633 + pkin(2) * t653 - qJ(3) * t783 - pkin(1) * t637 + (Ifges(2,6) - t777) * qJDD(1) + (Ifges(2,5) - t754 * (Ifges(3,4) * t754 + Ifges(3,2) * t757) + t757 * t778) * t762;
t626 = -mrSges(2,2) * g(1) - mrSges(2,3) * t736 + Ifges(2,5) * qJDD(1) - t762 * Ifges(2,6) - qJ(2) * t637 + t757 * t629 - t754 * t631;
t1 = [(-m(1) - m(2)) * g(1) + t637; -m(1) * g(2) + t798; -m(1) * g(3) + t785; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t767; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t785 + t759 * t626 + t761 * t627; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t798 - t761 * t626 + t759 * t627; t767; t645; t653; t658; t763;];
tauJB = t1;
