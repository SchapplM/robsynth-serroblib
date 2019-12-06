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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:31:24
% EndTime: 2019-12-05 17:31:30
% DurationCPUTime: 5.42s
% Computational Cost: add. (44592->280), mult. (124905->382), div. (0->0), fcn. (84227->10), ass. (0->133)
t740 = sin(pkin(7));
t743 = cos(pkin(7));
t745 = sin(qJ(1));
t747 = cos(qJ(1));
t723 = t745 * g(2) - t747 * g(3);
t748 = qJD(1) ^ 2;
t795 = -t748 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t723;
t694 = -t740 * g(1) + t795 * t743;
t762 = -pkin(2) * t743 - qJ(3) * t740;
t717 = t762 * qJD(1);
t782 = t743 * qJD(1);
t685 = t717 * t782 + t694;
t739 = sin(pkin(8));
t742 = cos(pkin(8));
t724 = t747 * g(2) + t745 * g(3);
t754 = -t748 * qJ(2) + qJDD(2) - t724;
t783 = t740 * qJD(1);
t794 = (-pkin(1) + t762) * qJDD(1) + t754 - 0.2e1 * qJD(3) * t783;
t665 = -t739 * t685 + t794 * t742;
t693 = -t743 * g(1) - t795 * t740;
t738 = sin(pkin(9));
t741 = cos(pkin(9));
t786 = t740 * t742;
t756 = t738 * t786 + t741 * t743;
t706 = t756 * qJD(1);
t704 = t756 * qJDD(1);
t793 = 2 * qJD(4);
t792 = mrSges(3,2) * t740;
t791 = Ifges(4,4) * t742;
t790 = Ifges(4,6) * t743;
t789 = t740 ^ 2 * t748;
t788 = t743 ^ 2 * t748;
t787 = t739 * t740;
t785 = t743 * t748;
t718 = (-mrSges(3,1) * t743 + t792) * qJD(1);
t666 = t742 * t685 + t794 * t739;
t766 = mrSges(4,1) * t739 + mrSges(4,2) * t742;
t711 = t766 * t783;
t758 = -mrSges(4,1) * t743 - mrSges(4,3) * t786;
t715 = t758 * qJD(1);
t757 = mrSges(4,2) * t743 - mrSges(4,3) * t787;
t710 = (pkin(3) * t739 - qJ(4) * t742) * t783;
t777 = t739 * t783;
t779 = qJDD(1) * t743;
t663 = -pkin(3) * t788 - qJ(4) * t779 - t710 * t777 + t666;
t684 = t717 * t783 + qJDD(3) - t693;
t671 = ((-qJDD(1) * t742 - t739 * t785) * qJ(4) + (qJDD(1) * t739 - t742 * t785) * pkin(3)) * t740 + t684;
t659 = t741 * t663 + t738 * t671 - t706 * t793;
t776 = t742 * t783;
t707 = -t738 * t782 + t741 * t776;
t686 = t706 * mrSges(5,1) + t707 * mrSges(5,2);
t692 = mrSges(5,1) * t777 - t707 * mrSges(5,3);
t687 = t706 * pkin(4) - t707 * pkin(6);
t780 = qJDD(1) * t740;
t773 = t739 * t780;
t778 = t739 ^ 2 * t789;
t657 = -pkin(4) * t778 + pkin(6) * t773 - t706 * t687 + t659;
t662 = pkin(3) * t779 - qJ(4) * t788 + t710 * t776 + qJDD(4) - t665;
t705 = (-t738 * t743 + t741 * t786) * qJDD(1);
t660 = (t706 * t777 - t705) * pkin(6) + (t707 * t777 + t704) * pkin(4) + t662;
t744 = sin(qJ(5));
t746 = cos(qJ(5));
t654 = -t744 * t657 + t746 * t660;
t688 = -t744 * t707 + t746 * t777;
t689 = t746 * t707 + t744 * t777;
t673 = -t688 * mrSges(6,1) + t689 * mrSges(6,2);
t675 = t688 * qJD(5) + t746 * t705 + t744 * t773;
t703 = qJD(5) + t706;
t676 = -t703 * mrSges(6,2) + t688 * mrSges(6,3);
t702 = qJDD(5) + t704;
t652 = m(6) * t654 + t702 * mrSges(6,1) - t675 * mrSges(6,3) - t689 * t673 + t703 * t676;
t655 = t746 * t657 + t744 * t660;
t674 = -t689 * qJD(5) - t744 * t705 + t746 * t773;
t677 = t703 * mrSges(6,1) - t689 * mrSges(6,3);
t653 = m(6) * t655 - t702 * mrSges(6,2) + t674 * mrSges(6,3) + t688 * t673 - t703 * t677;
t768 = -t744 * t652 + t746 * t653;
t645 = m(5) * t659 - t704 * mrSges(5,3) - t706 * t686 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t692) * t787 + t768;
t760 = t738 * t663 - t741 * t671;
t658 = -0.2e1 * qJD(4) * t707 - t760;
t691 = -mrSges(5,2) * t777 - t706 * mrSges(5,3);
t656 = -pkin(4) * t773 - pkin(6) * t778 + (t793 + t687) * t707 + t760;
t752 = -m(6) * t656 + t674 * mrSges(6,1) - t675 * mrSges(6,2) + t688 * t676 - t689 * t677;
t650 = m(5) * t658 - t705 * mrSges(5,3) - t707 * t686 + (mrSges(5,1) * qJDD(1) + qJD(1) * t691) * t787 + t752;
t769 = t741 * t645 - t738 * t650;
t640 = m(4) * t666 + t757 * qJDD(1) + (-t711 * t787 + t715 * t743) * qJD(1) + t769;
t647 = t746 * t652 + t744 * t653;
t646 = m(5) * t662 + t704 * mrSges(5,1) + t705 * mrSges(5,2) + t706 * t691 + t707 * t692 + t647;
t714 = t757 * qJD(1);
t643 = m(4) * t665 + t758 * qJDD(1) + (-t711 * t786 - t714 * t743) * qJD(1) - t646;
t770 = t742 * t640 - t739 * t643;
t631 = m(3) * t694 + (qJDD(1) * mrSges(3,3) + qJD(1) * t718) * t743 + t770;
t642 = t738 * t645 + t741 * t650;
t755 = m(4) * t684 + t642;
t759 = t714 * t739 + t715 * t742;
t638 = m(3) * t693 + ((-mrSges(3,3) - t766) * qJDD(1) + (-t718 - t759) * qJD(1)) * t740 - t755;
t627 = t740 * t631 + t743 * t638;
t771 = t743 * t631 - t740 * t638;
t625 = m(2) * t723 - t748 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t771;
t634 = t739 * t640 + t742 * t643;
t713 = -qJDD(1) * pkin(1) + t754;
t750 = -m(3) * t713 + mrSges(3,1) * t779 - t634 + (t788 + t789) * mrSges(3,3);
t628 = m(2) * t724 - t748 * mrSges(2,2) + (mrSges(2,1) - t792) * qJDD(1) + t750;
t772 = t747 * t625 - t745 * t628;
t765 = Ifges(3,1) * t740 + Ifges(3,4) * t743;
t764 = Ifges(3,5) * t740 + Ifges(3,6) * t743;
t763 = Ifges(4,5) * t742 - Ifges(4,6) * t739;
t761 = -t745 * t625 - t747 * t628;
t667 = Ifges(6,5) * t689 + Ifges(6,6) * t688 + Ifges(6,3) * t703;
t669 = Ifges(6,1) * t689 + Ifges(6,4) * t688 + Ifges(6,5) * t703;
t648 = -mrSges(6,1) * t656 + mrSges(6,3) * t655 + Ifges(6,4) * t675 + Ifges(6,2) * t674 + Ifges(6,6) * t702 - t689 * t667 + t703 * t669;
t668 = Ifges(6,4) * t689 + Ifges(6,2) * t688 + Ifges(6,6) * t703;
t649 = mrSges(6,2) * t656 - mrSges(6,3) * t654 + Ifges(6,1) * t675 + Ifges(6,4) * t674 + Ifges(6,5) * t702 + t688 * t667 - t703 * t668;
t678 = Ifges(5,5) * t707 - Ifges(5,6) * t706 + Ifges(5,3) * t777;
t679 = Ifges(5,4) * t707 - Ifges(5,2) * t706 + Ifges(5,6) * t777;
t635 = mrSges(5,2) * t662 - mrSges(5,3) * t658 + Ifges(5,1) * t705 - Ifges(5,4) * t704 - pkin(6) * t647 - t744 * t648 + t746 * t649 - t706 * t678 + (Ifges(5,5) * qJDD(1) - qJD(1) * t679) * t787;
t680 = Ifges(5,1) * t707 - Ifges(5,4) * t706 + Ifges(5,5) * t777;
t749 = mrSges(6,1) * t654 - mrSges(6,2) * t655 + Ifges(6,5) * t675 + Ifges(6,6) * t674 + Ifges(6,3) * t702 + t689 * t668 - t688 * t669;
t636 = -mrSges(5,1) * t662 + mrSges(5,3) * t659 + Ifges(5,4) * t705 - Ifges(5,2) * t704 - pkin(4) * t647 - t707 * t678 + (Ifges(5,6) * qJDD(1) + qJD(1) * t680) * t787 - t749;
t697 = (-Ifges(4,3) * t743 + t763 * t740) * qJD(1);
t698 = (-t790 + (-Ifges(4,2) * t739 + t791) * t740) * qJD(1);
t751 = -Ifges(4,5) * t743 + (Ifges(4,1) * t742 - Ifges(4,4) * t739) * t740;
t622 = mrSges(4,2) * t684 - mrSges(4,3) * t665 - qJ(4) * t642 + t741 * t635 - t738 * t636 + (-t697 * t787 + t698 * t743) * qJD(1) + t751 * qJDD(1);
t699 = t751 * qJD(1);
t623 = -mrSges(4,1) * t684 + mrSges(4,3) * t666 - Ifges(5,5) * t705 + Ifges(5,6) * t704 - t707 * t679 - t706 * t680 - mrSges(5,1) * t658 + mrSges(5,2) * t659 - t744 * t649 - t746 * t648 - pkin(4) * t752 - pkin(6) * t768 - pkin(3) * t642 + (-t697 * t786 - t743 * t699) * qJD(1) + (-t790 + (t791 + (-Ifges(4,2) - Ifges(5,3)) * t739) * t740) * qJDD(1);
t719 = t764 * qJD(1);
t619 = mrSges(3,2) * t713 - mrSges(3,3) * t693 - qJ(3) * t634 + t765 * qJDD(1) + t742 * t622 - t739 * t623 + t719 * t782;
t621 = -mrSges(3,1) * t713 + mrSges(3,3) * t694 - mrSges(4,1) * t665 + mrSges(4,2) * t666 - t738 * t635 - t741 * t636 + pkin(3) * t646 - qJ(4) * t769 - pkin(2) * t634 + (Ifges(3,2) + Ifges(4,3)) * t779 + ((Ifges(3,4) - t763) * qJDD(1) + (-t698 * t742 - t699 * t739 - t719) * qJD(1)) * t740;
t633 = mrSges(3,2) * t780 - t750;
t753 = mrSges(2,1) * t724 - mrSges(2,2) * t723 + Ifges(2,3) * qJDD(1) - pkin(1) * t633 + qJ(2) * t771 + t740 * t619 + t743 * t621;
t641 = (t759 * qJD(1) + t766 * qJDD(1)) * t740 + t755;
t617 = mrSges(2,1) * g(1) + mrSges(2,3) * t723 - mrSges(3,1) * t693 + mrSges(3,2) * t694 - t739 * t622 - t742 * t623 + pkin(2) * t641 - qJ(3) * t770 - pkin(1) * t627 + (Ifges(2,6) - t764) * qJDD(1) + (Ifges(2,5) - t740 * (Ifges(3,4) * t740 + Ifges(3,2) * t743) + t743 * t765) * t748;
t616 = -mrSges(2,2) * g(1) - mrSges(2,3) * t724 + Ifges(2,5) * qJDD(1) - t748 * Ifges(2,6) - qJ(2) * t627 + t743 * t619 - t740 * t621;
t1 = [(-m(1) - m(2)) * g(1) + t627; -m(1) * g(2) + t761; -m(1) * g(3) + t772; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t753; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t772 - t745 * t616 - t747 * t617; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t761 + t747 * t616 - t745 * t617; t753; t633; t641; t646; t749;];
tauJB = t1;
