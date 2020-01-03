% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:32
% EndTime: 2019-12-31 21:19:38
% DurationCPUTime: 4.50s
% Computational Cost: add. (42892->294), mult. (89100->355), div. (0->0), fcn. (57798->8), ass. (0->121)
t799 = Ifges(4,4) + Ifges(5,6);
t811 = -Ifges(4,2) - Ifges(5,3);
t807 = Ifges(4,6) - Ifges(5,5);
t766 = sin(qJ(3));
t770 = cos(qJ(2));
t791 = qJD(1) * t770;
t767 = sin(qJ(2));
t792 = qJD(1) * t767;
t801 = cos(qJ(3));
t738 = t766 * t792 - t801 * t791;
t762 = qJD(2) + qJD(3);
t727 = t738 * mrSges(5,1) - t762 * mrSges(5,3);
t761 = qJDD(2) + qJDD(3);
t739 = (t766 * t770 + t801 * t767) * qJD(1);
t790 = qJD(1) * qJD(2);
t747 = t767 * qJDD(1) + t770 * t790;
t748 = t770 * qJDD(1) - t767 * t790;
t703 = t739 * qJD(3) + t766 * t747 - t801 * t748;
t729 = t739 * pkin(4) - t762 * pkin(8);
t734 = t738 ^ 2;
t704 = -t738 * qJD(3) + t801 * t747 + t766 * t748;
t752 = qJD(2) * pkin(2) - pkin(7) * t792;
t764 = t770 ^ 2;
t772 = qJD(1) ^ 2;
t768 = sin(qJ(1));
t771 = cos(qJ(1));
t753 = t768 * g(1) - t771 * g(2);
t784 = -qJDD(1) * pkin(1) - t753;
t705 = -t748 * pkin(2) + t752 * t792 + (-pkin(7) * t764 - pkin(6)) * t772 + t784;
t798 = t738 * t762;
t802 = -2 * qJD(4);
t776 = (-t704 + t798) * qJ(4) + t705 + (t762 * pkin(3) + t802) * t739;
t661 = -t734 * pkin(4) - t739 * t729 + (pkin(3) + pkin(8)) * t703 + t776;
t754 = -t771 * g(1) - t768 * g(2);
t741 = -t772 * pkin(1) + qJDD(1) * pkin(6) + t754;
t797 = t767 * t741;
t800 = pkin(2) * t772;
t689 = qJDD(2) * pkin(2) - t747 * pkin(7) - t797 + (pkin(7) * t790 + t767 * t800 - g(3)) * t770;
t724 = -t767 * g(3) + t770 * t741;
t690 = t748 * pkin(7) - qJD(2) * t752 - t764 * t800 + t724;
t671 = t801 * t689 - t766 * t690;
t717 = t738 * pkin(3) - t739 * qJ(4);
t760 = t762 ^ 2;
t669 = -t761 * pkin(3) - t760 * qJ(4) + t739 * t717 + qJDD(4) - t671;
t662 = (t738 * t739 - t761) * pkin(8) + (t704 + t798) * pkin(4) + t669;
t765 = sin(qJ(5));
t769 = cos(qJ(5));
t659 = -t765 * t661 + t769 * t662;
t721 = t769 * t738 - t765 * t762;
t678 = t721 * qJD(5) + t765 * t703 + t769 * t761;
t722 = t765 * t738 + t769 * t762;
t686 = -t721 * mrSges(6,1) + t722 * mrSges(6,2);
t702 = qJDD(5) + t704;
t733 = qJD(5) + t739;
t706 = -t733 * mrSges(6,2) + t721 * mrSges(6,3);
t656 = m(6) * t659 + t702 * mrSges(6,1) - t678 * mrSges(6,3) - t722 * t686 + t733 * t706;
t660 = t769 * t661 + t765 * t662;
t677 = -t722 * qJD(5) + t769 * t703 - t765 * t761;
t707 = t733 * mrSges(6,1) - t722 * mrSges(6,3);
t657 = m(6) * t660 - t702 * mrSges(6,2) + t677 * mrSges(6,3) + t721 * t686 - t733 * t707;
t645 = t769 * t656 + t765 * t657;
t719 = -t738 * mrSges(5,2) - t739 * mrSges(5,3);
t781 = -m(5) * t669 - t704 * mrSges(5,1) - t739 * t719 - t645;
t644 = t761 * mrSges(5,2) + t762 * t727 - t781;
t672 = t766 * t689 + t801 * t690;
t780 = -t760 * pkin(3) + t761 * qJ(4) - t738 * t717 + t672;
t664 = -t703 * pkin(4) - t734 * pkin(8) + ((2 * qJD(4)) + t729) * t762 + t780;
t679 = Ifges(6,5) * t722 + Ifges(6,6) * t721 + Ifges(6,3) * t733;
t681 = Ifges(6,1) * t722 + Ifges(6,4) * t721 + Ifges(6,5) * t733;
t647 = -mrSges(6,1) * t664 + mrSges(6,3) * t660 + Ifges(6,4) * t678 + Ifges(6,2) * t677 + Ifges(6,6) * t702 - t722 * t679 + t733 * t681;
t680 = Ifges(6,4) * t722 + Ifges(6,2) * t721 + Ifges(6,6) * t733;
t648 = mrSges(6,2) * t664 - mrSges(6,3) * t659 + Ifges(6,1) * t678 + Ifges(6,4) * t677 + Ifges(6,5) * t702 + t721 * t679 - t733 * t680;
t667 = t762 * t802 - t780;
t728 = t739 * mrSges(5,1) + t762 * mrSges(5,2);
t783 = -m(6) * t664 + t677 * mrSges(6,1) - t678 * mrSges(6,2) + t721 * t706 - t722 * t707;
t778 = -m(5) * t667 + t761 * mrSges(5,3) + t762 * t728 - t783;
t808 = Ifges(4,5) - Ifges(5,4);
t809 = Ifges(4,1) + Ifges(5,2);
t793 = -t799 * t738 + t809 * t739 + t808 * t762;
t805 = t811 * t738 + t799 * t739 + t807 * t762;
t806 = Ifges(4,3) + Ifges(5,1);
t810 = -mrSges(4,2) * t672 - mrSges(5,3) * t667 - pkin(3) * t644 - pkin(8) * t645 - t765 * t647 + t769 * t648 + qJ(4) * (-t738 * t719 + t778) + mrSges(5,2) * t669 + mrSges(4,1) * t671 + t806 * t761 + t805 * t739 + t808 * t704 + (-qJ(4) * mrSges(5,1) - t807) * t703 + t793 * t738;
t718 = t738 * mrSges(4,1) + t739 * mrSges(4,2);
t725 = -t762 * mrSges(4,2) - t738 * mrSges(4,3);
t641 = m(4) * t671 - t704 * mrSges(4,3) - t739 * t718 + (t725 - t727) * t762 + (mrSges(4,1) - mrSges(5,2)) * t761 + t781;
t726 = t762 * mrSges(4,1) - t739 * mrSges(4,3);
t651 = m(4) * t672 - t761 * mrSges(4,2) - t762 * t726 + (-t718 - t719) * t738 + (-mrSges(4,3) - mrSges(5,1)) * t703 + t778;
t636 = t801 * t641 + t766 * t651;
t723 = -t770 * g(3) - t797;
t736 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t767 + Ifges(3,2) * t770) * qJD(1);
t737 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t767 + Ifges(3,4) * t770) * qJD(1);
t804 = mrSges(3,1) * t723 - mrSges(3,2) * t724 + Ifges(3,5) * t747 + Ifges(3,6) * t748 + Ifges(3,3) * qJDD(2) + pkin(2) * t636 + (t767 * t736 - t770 * t737) * qJD(1) + t810;
t746 = (-mrSges(3,1) * t770 + mrSges(3,2) * t767) * qJD(1);
t751 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t791;
t634 = m(3) * t723 + qJDD(2) * mrSges(3,1) - t747 * mrSges(3,3) + qJD(2) * t751 - t746 * t792 + t636;
t750 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t792;
t786 = -t766 * t641 + t801 * t651;
t635 = m(3) * t724 - qJDD(2) * mrSges(3,2) + t748 * mrSges(3,3) - qJD(2) * t750 + t746 * t791 + t786;
t787 = -t767 * t634 + t770 * t635;
t627 = m(2) * t754 - t772 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t787;
t740 = -t772 * pkin(6) + t784;
t666 = t703 * pkin(3) + t776;
t795 = -t765 * t656 + t769 * t657;
t642 = m(5) * t666 - t703 * mrSges(5,2) - t704 * mrSges(5,3) - t738 * t727 - t739 * t728 + t795;
t777 = m(4) * t705 + t703 * mrSges(4,1) + t704 * mrSges(4,2) + t738 * t725 + t739 * t726 + t642;
t774 = -m(3) * t740 + t748 * mrSges(3,1) - t747 * mrSges(3,2) - t750 * t792 + t751 * t791 - t777;
t638 = m(2) * t753 + qJDD(1) * mrSges(2,1) - t772 * mrSges(2,2) + t774;
t796 = t768 * t627 + t771 * t638;
t629 = t770 * t634 + t767 * t635;
t794 = t807 * t738 - t808 * t739 - t806 * t762;
t788 = t771 * t627 - t768 * t638;
t624 = -mrSges(4,1) * t705 - mrSges(5,1) * t667 + mrSges(5,2) * t666 + mrSges(4,3) * t672 - pkin(3) * t642 - pkin(4) * t783 - pkin(8) * t795 - t769 * t647 - t765 * t648 + t811 * t703 + t799 * t704 + t794 * t739 + t807 * t761 + t793 * t762;
t779 = mrSges(6,1) * t659 - mrSges(6,2) * t660 + Ifges(6,5) * t678 + Ifges(6,6) * t677 + Ifges(6,3) * t702 + t722 * t680 - t721 * t681;
t630 = mrSges(5,1) * t669 + mrSges(4,2) * t705 - mrSges(4,3) * t671 - mrSges(5,3) * t666 + pkin(4) * t645 - qJ(4) * t642 - t799 * t703 + t809 * t704 + t794 * t738 + t808 * t761 - t805 * t762 + t779;
t735 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t767 + Ifges(3,6) * t770) * qJD(1);
t620 = -mrSges(3,1) * t740 + mrSges(3,3) * t724 + Ifges(3,4) * t747 + Ifges(3,2) * t748 + Ifges(3,6) * qJDD(2) - pkin(2) * t777 + pkin(7) * t786 + qJD(2) * t737 + t801 * t624 + t766 * t630 - t735 * t792;
t623 = mrSges(3,2) * t740 - mrSges(3,3) * t723 + Ifges(3,1) * t747 + Ifges(3,4) * t748 + Ifges(3,5) * qJDD(2) - pkin(7) * t636 - qJD(2) * t736 - t766 * t624 + t801 * t630 + t735 * t791;
t782 = mrSges(2,1) * t753 - mrSges(2,2) * t754 + Ifges(2,3) * qJDD(1) + pkin(1) * t774 + pkin(6) * t787 + t770 * t620 + t767 * t623;
t621 = mrSges(2,1) * g(3) + mrSges(2,3) * t754 + t772 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t629 - t804;
t618 = -mrSges(2,2) * g(3) - mrSges(2,3) * t753 + Ifges(2,5) * qJDD(1) - t772 * Ifges(2,6) - pkin(6) * t629 - t767 * t620 + t770 * t623;
t1 = [-m(1) * g(1) + t788; -m(1) * g(2) + t796; (-m(1) - m(2)) * g(3) + t629; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t796 + t771 * t618 - t768 * t621; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t788 + t768 * t618 + t771 * t621; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t782; t782; t804; t810; t644; t779;];
tauJB = t1;
