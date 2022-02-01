% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:11
% EndTime: 2022-01-23 09:32:17
% DurationCPUTime: 4.19s
% Computational Cost: add. (35957->268), mult. (86495->338), div. (0->0), fcn. (56474->8), ass. (0->118)
t806 = Ifges(5,4) + Ifges(6,4);
t816 = Ifges(5,2) + Ifges(6,2);
t811 = Ifges(5,6) + Ifges(6,6);
t765 = sin(qJ(1));
t768 = cos(qJ(1));
t747 = -t768 * g(1) - t765 * g(2);
t769 = qJD(1) ^ 2;
t815 = -t769 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t747;
t763 = sin(qJ(4));
t764 = sin(qJ(3));
t766 = cos(qJ(4));
t767 = cos(qJ(3));
t761 = sin(pkin(8));
t796 = t761 * qJD(1);
t724 = (-t763 * t767 - t764 * t766) * t796;
t793 = qJD(1) * qJD(3);
t733 = (-qJDD(1) * t764 - t767 * t793) * t761;
t734 = (qJDD(1) * t767 - t764 * t793) * t761;
t693 = t724 * qJD(4) + t763 * t733 + t766 * t734;
t725 = (-t763 * t764 + t766 * t767) * t796;
t706 = -t724 * mrSges(6,1) + t725 * mrSges(6,2);
t762 = cos(pkin(8));
t717 = -t761 * g(3) + t815 * t762;
t780 = -pkin(2) * t762 - pkin(6) * t761;
t740 = t780 * qJD(1);
t795 = t762 * qJD(1);
t705 = t740 * t795 + t717;
t746 = t765 * g(1) - t768 * g(2);
t776 = -t769 * qJ(2) + qJDD(2) - t746;
t718 = (-pkin(1) + t780) * qJDD(1) + t776;
t715 = t767 * t718;
t792 = t762 * qJDD(1);
t749 = qJDD(3) - t792;
t750 = qJD(3) - t795;
t802 = t761 ^ 2 * t769;
t676 = t749 * pkin(3) - t734 * pkin(7) + t715 + (-pkin(3) * t767 * t802 - pkin(7) * t750 * t796 - t705) * t764;
t680 = t767 * t705 + t764 * t718;
t787 = t767 * t796;
t732 = t750 * pkin(3) - pkin(7) * t787;
t791 = t764 ^ 2 * t802;
t677 = -pkin(3) * t791 + t733 * pkin(7) - t750 * t732 + t680;
t669 = t766 * t676 - t763 * t677;
t745 = qJDD(4) + t749;
t748 = qJD(4) + t750;
t665 = -0.2e1 * qJD(5) * t725 + (t724 * t748 - t693) * qJ(5) + (t724 * t725 + t745) * pkin(4) + t669;
t709 = -t748 * mrSges(6,2) + t724 * mrSges(6,3);
t790 = m(6) * t665 + t745 * mrSges(6,1) + t748 * t709;
t661 = -t693 * mrSges(6,3) - t725 * t706 + t790;
t670 = t763 * t676 + t766 * t677;
t692 = -t725 * qJD(4) + t766 * t733 - t763 * t734;
t711 = t748 * pkin(4) - t725 * qJ(5);
t723 = t724 ^ 2;
t667 = -t723 * pkin(4) + t692 * qJ(5) + 0.2e1 * qJD(5) * t724 - t748 * t711 + t670;
t812 = Ifges(5,5) + Ifges(6,5);
t813 = Ifges(5,1) + Ifges(6,1);
t799 = -t806 * t724 - t813 * t725 - t812 * t748;
t809 = t816 * t724 + t806 * t725 + t811 * t748;
t810 = Ifges(5,3) + Ifges(6,3);
t814 = mrSges(5,1) * t669 + mrSges(6,1) * t665 - mrSges(5,2) * t670 - mrSges(6,2) * t667 + pkin(4) * t661 + t811 * t692 + t812 * t693 + t799 * t724 + t809 * t725 + t810 * t745;
t707 = -t724 * mrSges(5,1) + t725 * mrSges(5,2);
t710 = -t748 * mrSges(5,2) + t724 * mrSges(5,3);
t654 = m(5) * t669 + t745 * mrSges(5,1) + t748 * t710 + (-t706 - t707) * t725 + (-mrSges(5,3) - mrSges(6,3)) * t693 + t790;
t712 = t748 * mrSges(6,1) - t725 * mrSges(6,3);
t713 = t748 * mrSges(5,1) - t725 * mrSges(5,3);
t789 = m(6) * t667 + t692 * mrSges(6,3) + t724 * t706;
t657 = m(5) * t670 + t692 * mrSges(5,3) + t724 * t707 + (-t712 - t713) * t748 + (-mrSges(5,2) - mrSges(6,2)) * t745 + t789;
t652 = t766 * t654 + t763 * t657;
t679 = -t764 * t705 + t715;
t808 = -mrSges(4,1) * t679 + mrSges(4,2) * t680 - Ifges(4,5) * t734 - Ifges(4,6) * t733 - Ifges(4,3) * t749 - pkin(3) * t652 - t814;
t716 = -t762 * g(3) - t815 * t761;
t788 = t764 * t796;
t729 = -t750 * mrSges(4,2) - mrSges(4,3) * t788;
t730 = t750 * mrSges(4,1) - mrSges(4,3) * t787;
t807 = -t729 * t764 - t730 * t767;
t805 = mrSges(3,2) * t761;
t738 = (-mrSges(3,1) * t762 + t805) * qJD(1);
t731 = (mrSges(4,1) * t764 + mrSges(4,2) * t767) * t796;
t649 = m(4) * t679 + t749 * mrSges(4,1) - t734 * mrSges(4,3) + t750 * t729 - t731 * t787 + t652;
t781 = -t763 * t654 + t766 * t657;
t650 = m(4) * t680 - t749 * mrSges(4,2) + t733 * mrSges(4,3) - t750 * t730 - t731 * t788 + t781;
t782 = -t764 * t649 + t767 * t650;
t797 = qJDD(1) * mrSges(3,3);
t643 = m(3) * t717 + (qJD(1) * t738 + t797) * t762 + t782;
t704 = t740 * t796 - t716;
t678 = -t733 * pkin(3) - pkin(7) * t791 + t732 * t787 + t704;
t672 = -t692 * pkin(4) - t723 * qJ(5) + t725 * t711 + qJDD(5) + t678;
t662 = m(6) * t672 - t692 * mrSges(6,1) + t693 * mrSges(6,2) - t724 * t709 + t725 * t712;
t773 = m(5) * t678 - t692 * mrSges(5,1) + t693 * mrSges(5,2) - t724 * t710 + t725 * t713 + t662;
t771 = -m(4) * t704 + t733 * mrSges(4,1) - t734 * mrSges(4,2) - t773;
t659 = t771 + m(3) * t716 + (-t797 + (-t738 + t807) * qJD(1)) * t761;
t783 = t762 * t643 - t761 * t659;
t636 = m(2) * t747 - t769 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t783;
t646 = t767 * t649 + t764 * t650;
t736 = -qJDD(1) * pkin(1) + t776;
t774 = -m(3) * t736 + mrSges(3,1) * t792 - t646 + (t762 ^ 2 * t769 + t802) * mrSges(3,3);
t640 = m(2) * t746 - t769 * mrSges(2,2) + (mrSges(2,1) - t805) * qJDD(1) + t774;
t801 = t765 * t636 + t768 * t640;
t638 = t761 * t643 + t762 * t659;
t800 = -t811 * t724 - t812 * t725 - t810 * t748;
t784 = t768 * t636 - t765 * t640;
t779 = Ifges(3,1) * t761 + Ifges(3,4) * t762;
t778 = Ifges(3,5) * t761 + Ifges(3,6) * t762;
t720 = Ifges(4,6) * t750 + (Ifges(4,4) * t767 - Ifges(4,2) * t764) * t796;
t721 = Ifges(4,5) * t750 + (Ifges(4,1) * t767 - Ifges(4,4) * t764) * t796;
t777 = t720 * t767 + t721 * t764;
t647 = -mrSges(5,1) * t678 + mrSges(5,3) * t670 - mrSges(6,1) * t672 + mrSges(6,3) * t667 - pkin(4) * t662 + qJ(5) * t789 + (-qJ(5) * t712 - t799) * t748 + (-qJ(5) * mrSges(6,2) + t811) * t745 + t800 * t725 + t806 * t693 + t816 * t692;
t651 = mrSges(5,2) * t678 + mrSges(6,2) * t672 - mrSges(5,3) * t669 - mrSges(6,3) * t665 - qJ(5) * t661 + t806 * t692 + t813 * t693 - t800 * t724 + t812 * t745 - t809 * t748;
t719 = Ifges(4,3) * t750 + (Ifges(4,5) * t767 - Ifges(4,6) * t764) * t796;
t632 = -mrSges(4,1) * t704 + mrSges(4,3) * t680 + Ifges(4,4) * t734 + Ifges(4,2) * t733 + Ifges(4,6) * t749 - pkin(3) * t773 + pkin(7) * t781 + t766 * t647 + t763 * t651 - t719 * t787 + t750 * t721;
t633 = mrSges(4,2) * t704 - mrSges(4,3) * t679 + Ifges(4,1) * t734 + Ifges(4,4) * t733 + Ifges(4,5) * t749 - pkin(7) * t652 - t763 * t647 + t766 * t651 - t719 * t788 - t750 * t720;
t739 = t778 * qJD(1);
t629 = mrSges(3,2) * t736 - mrSges(3,3) * t716 - pkin(6) * t646 + qJDD(1) * t779 - t764 * t632 + t767 * t633 + t739 * t795;
t631 = (Ifges(3,4) * qJDD(1) + (-t739 - t777) * qJD(1)) * t761 + mrSges(3,3) * t717 - pkin(2) * t646 + Ifges(3,2) * t792 - mrSges(3,1) * t736 + t808;
t645 = qJDD(1) * t805 - t774;
t775 = mrSges(2,1) * t746 - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) - pkin(1) * t645 + qJ(2) * t783 + t761 * t629 + t762 * t631;
t627 = t769 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t747 - mrSges(3,1) * t716 + mrSges(3,2) * t717 - t764 * t633 - t767 * t632 - pkin(2) * t771 - pkin(6) * t782 - pkin(1) * t638 + (Ifges(2,6) - t778) * qJDD(1) + (-pkin(2) * t807 * t761 + (-t761 * (Ifges(3,4) * t761 + Ifges(3,2) * t762) + t762 * t779) * qJD(1)) * qJD(1);
t626 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - t769 * Ifges(2,6) - qJ(2) * t638 + t762 * t629 - t761 * t631;
t1 = [-m(1) * g(1) + t784; -m(1) * g(2) + t801; (-m(1) - m(2)) * g(3) + t638; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t801 + t768 * t626 - t765 * t627; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t784 + t765 * t626 + t768 * t627; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t775; t775; t645; t777 * t796 - t808; t814; t662;];
tauJB = t1;
