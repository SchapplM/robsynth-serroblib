% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:56:38
% EndTime: 2019-12-05 17:56:47
% DurationCPUTime: 8.05s
% Computational Cost: add. (80930->290), mult. (207955->381), div. (0->0), fcn. (141748->10), ass. (0->126)
t751 = sin(qJ(1));
t754 = cos(qJ(1));
t730 = t751 * g(2) - t754 * g(3);
t755 = qJD(1) ^ 2;
t791 = -t755 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t730;
t746 = sin(pkin(8));
t748 = cos(pkin(8));
t700 = -t748 * g(1) - t791 * t746;
t780 = t748 * qJD(1);
t734 = qJD(3) - t780;
t750 = sin(qJ(3));
t781 = t746 * qJD(1);
t775 = t750 * t781;
t715 = -t734 * mrSges(4,2) - mrSges(4,3) * t775;
t753 = cos(qJ(3));
t774 = t753 * t781;
t717 = t734 * mrSges(4,1) - mrSges(4,3) * t774;
t790 = -t715 * t750 - t717 * t753;
t701 = -t746 * g(1) + t791 * t748;
t767 = -pkin(2) * t748 - pkin(6) * t746;
t726 = t767 * qJD(1);
t690 = t726 * t780 + t701;
t731 = t754 * g(2) + t751 * g(3);
t761 = -t755 * qJ(2) + qJDD(2) - t731;
t702 = (-pkin(1) + t767) * qJDD(1) + t761;
t699 = t753 * t702;
t778 = qJD(1) * qJD(3);
t720 = (qJDD(1) * t753 - t750 * t778) * t746;
t777 = t748 * qJDD(1);
t733 = qJDD(3) - t777;
t784 = t746 ^ 2 * t755;
t669 = t733 * pkin(3) - t720 * qJ(4) + t699 + (-pkin(3) * t753 * t784 - qJ(4) * t734 * t781 - t690) * t750;
t679 = t753 * t690 + t750 * t702;
t716 = t734 * pkin(3) - qJ(4) * t774;
t719 = (-qJDD(1) * t750 - t753 * t778) * t746;
t776 = t750 ^ 2 * t784;
t670 = -pkin(3) * t776 + t719 * qJ(4) - t734 * t716 + t679;
t745 = sin(pkin(9));
t747 = cos(pkin(9));
t711 = (-t745 * t750 + t747 * t753) * t781;
t655 = -0.2e1 * qJD(4) * t711 + t747 * t669 - t745 * t670;
t694 = t745 * t719 + t747 * t720;
t710 = (-t745 * t753 - t747 * t750) * t781;
t653 = (t710 * t734 - t694) * pkin(7) + (t710 * t711 + t733) * pkin(4) + t655;
t656 = 0.2e1 * qJD(4) * t710 + t745 * t669 + t747 * t670;
t693 = t747 * t719 - t745 * t720;
t697 = t734 * pkin(4) - t711 * pkin(7);
t709 = t710 ^ 2;
t654 = -t709 * pkin(4) + t693 * pkin(7) - t734 * t697 + t656;
t749 = sin(qJ(5));
t752 = cos(qJ(5));
t651 = t752 * t653 - t749 * t654;
t687 = t752 * t710 - t749 * t711;
t665 = t687 * qJD(5) + t749 * t693 + t752 * t694;
t688 = t749 * t710 + t752 * t711;
t676 = -t687 * mrSges(6,1) + t688 * mrSges(6,2);
t732 = qJD(5) + t734;
t680 = -t732 * mrSges(6,2) + t687 * mrSges(6,3);
t729 = qJDD(5) + t733;
t647 = m(6) * t651 + t729 * mrSges(6,1) - t665 * mrSges(6,3) - t688 * t676 + t732 * t680;
t652 = t749 * t653 + t752 * t654;
t664 = -t688 * qJD(5) + t752 * t693 - t749 * t694;
t681 = t732 * mrSges(6,1) - t688 * mrSges(6,3);
t648 = m(6) * t652 - t729 * mrSges(6,2) + t664 * mrSges(6,3) + t687 * t676 - t732 * t681;
t639 = t752 * t647 + t749 * t648;
t691 = -t710 * mrSges(5,1) + t711 * mrSges(5,2);
t695 = -t734 * mrSges(5,2) + t710 * mrSges(5,3);
t637 = m(5) * t655 + t733 * mrSges(5,1) - t694 * mrSges(5,3) - t711 * t691 + t734 * t695 + t639;
t696 = t734 * mrSges(5,1) - t711 * mrSges(5,3);
t768 = -t749 * t647 + t752 * t648;
t638 = m(5) * t656 - t733 * mrSges(5,2) + t693 * mrSges(5,3) + t710 * t691 - t734 * t696 + t768;
t633 = t747 * t637 + t745 * t638;
t678 = -t750 * t690 + t699;
t683 = Ifges(5,4) * t711 + Ifges(5,2) * t710 + Ifges(5,6) * t734;
t684 = Ifges(5,1) * t711 + Ifges(5,4) * t710 + Ifges(5,5) * t734;
t672 = Ifges(6,4) * t688 + Ifges(6,2) * t687 + Ifges(6,6) * t732;
t673 = Ifges(6,1) * t688 + Ifges(6,4) * t687 + Ifges(6,5) * t732;
t759 = -mrSges(6,1) * t651 + mrSges(6,2) * t652 - Ifges(6,5) * t665 - Ifges(6,6) * t664 - Ifges(6,3) * t729 - t688 * t672 + t687 * t673;
t789 = -mrSges(4,1) * t678 - mrSges(5,1) * t655 + mrSges(4,2) * t679 + mrSges(5,2) * t656 - Ifges(4,5) * t720 - Ifges(5,5) * t694 - Ifges(4,6) * t719 - Ifges(5,6) * t693 - pkin(3) * t633 - pkin(4) * t639 - t711 * t683 + t710 * t684 - (Ifges(4,3) + Ifges(5,3)) * t733 + t759;
t787 = mrSges(3,2) * t746;
t724 = (-mrSges(3,1) * t748 + t787) * qJD(1);
t718 = (mrSges(4,1) * t750 + mrSges(4,2) * t753) * t781;
t631 = m(4) * t678 + t733 * mrSges(4,1) - t720 * mrSges(4,3) + t734 * t715 - t718 * t774 + t633;
t769 = -t745 * t637 + t747 * t638;
t632 = m(4) * t679 - t733 * mrSges(4,2) + t719 * mrSges(4,3) - t734 * t717 - t718 * t775 + t769;
t770 = -t750 * t631 + t753 * t632;
t782 = qJDD(1) * mrSges(3,3);
t624 = m(3) * t701 + (qJD(1) * t724 + t782) * t748 + t770;
t689 = t726 * t781 - t700;
t677 = -t719 * pkin(3) - qJ(4) * t776 + t716 * t774 + qJDD(4) + t689;
t658 = -t693 * pkin(4) - t709 * pkin(7) + t711 * t697 + t677;
t764 = m(6) * t658 - t664 * mrSges(6,1) + t665 * mrSges(6,2) - t687 * t680 + t688 * t681;
t649 = m(5) * t677 - t693 * mrSges(5,1) + t694 * mrSges(5,2) - t710 * t695 + t711 * t696 + t764;
t757 = -m(4) * t689 + t719 * mrSges(4,1) - t720 * mrSges(4,2) - t649;
t646 = t757 + m(3) * t700 + (-t782 + (-t724 + t790) * qJD(1)) * t746;
t620 = t746 * t624 + t748 * t646;
t771 = t748 * t624 - t746 * t646;
t618 = m(2) * t730 - t755 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t771;
t627 = t753 * t631 + t750 * t632;
t722 = -qJDD(1) * pkin(1) + t761;
t758 = -m(3) * t722 + mrSges(3,1) * t777 - t627 + (t748 ^ 2 * t755 + t784) * mrSges(3,3);
t621 = m(2) * t731 - t755 * mrSges(2,2) + (mrSges(2,1) - t787) * qJDD(1) + t758;
t772 = t754 * t618 - t751 * t621;
t766 = Ifges(3,1) * t746 + Ifges(3,4) * t748;
t765 = Ifges(3,5) * t746 + Ifges(3,6) * t748;
t763 = -t751 * t618 - t754 * t621;
t704 = Ifges(4,6) * t734 + (Ifges(4,4) * t753 - Ifges(4,2) * t750) * t781;
t705 = Ifges(4,5) * t734 + (Ifges(4,1) * t753 - Ifges(4,4) * t750) * t781;
t762 = t704 * t753 + t705 * t750;
t671 = Ifges(6,5) * t688 + Ifges(6,6) * t687 + Ifges(6,3) * t732;
t640 = -mrSges(6,1) * t658 + mrSges(6,3) * t652 + Ifges(6,4) * t665 + Ifges(6,2) * t664 + Ifges(6,6) * t729 - t688 * t671 + t732 * t673;
t641 = mrSges(6,2) * t658 - mrSges(6,3) * t651 + Ifges(6,1) * t665 + Ifges(6,4) * t664 + Ifges(6,5) * t729 + t687 * t671 - t732 * t672;
t682 = Ifges(5,5) * t711 + Ifges(5,6) * t710 + Ifges(5,3) * t734;
t628 = -mrSges(5,1) * t677 + mrSges(5,3) * t656 + Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t733 - pkin(4) * t764 + pkin(7) * t768 + t752 * t640 + t749 * t641 - t711 * t682 + t734 * t684;
t629 = mrSges(5,2) * t677 - mrSges(5,3) * t655 + Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t733 - pkin(7) * t639 - t749 * t640 + t752 * t641 + t710 * t682 - t734 * t683;
t703 = Ifges(4,3) * t734 + (Ifges(4,5) * t753 - Ifges(4,6) * t750) * t781;
t615 = -mrSges(4,1) * t689 + mrSges(4,3) * t679 + Ifges(4,4) * t720 + Ifges(4,2) * t719 + Ifges(4,6) * t733 - pkin(3) * t649 + qJ(4) * t769 + t747 * t628 + t745 * t629 - t703 * t774 + t734 * t705;
t616 = mrSges(4,2) * t689 - mrSges(4,3) * t678 + Ifges(4,1) * t720 + Ifges(4,4) * t719 + Ifges(4,5) * t733 - qJ(4) * t633 - t745 * t628 + t747 * t629 - t703 * t775 - t734 * t704;
t725 = t765 * qJD(1);
t612 = mrSges(3,2) * t722 - mrSges(3,3) * t700 - pkin(6) * t627 + t766 * qJDD(1) - t750 * t615 + t753 * t616 + t725 * t780;
t614 = (Ifges(3,4) * qJDD(1) + (-t725 - t762) * qJD(1)) * t746 - pkin(2) * t627 + mrSges(3,3) * t701 - mrSges(3,1) * t722 + Ifges(3,2) * t777 + t789;
t626 = qJDD(1) * t787 - t758;
t760 = mrSges(2,1) * t731 - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) - pkin(1) * t626 + qJ(2) * t771 + t746 * t612 + t748 * t614;
t610 = t755 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t730 - mrSges(3,1) * t700 + mrSges(3,2) * t701 - t750 * t616 - t753 * t615 - pkin(2) * t757 - pkin(6) * t770 - pkin(1) * t620 + (Ifges(2,6) - t765) * qJDD(1) + (-pkin(2) * t790 * t746 + (-t746 * (Ifges(3,4) * t746 + Ifges(3,2) * t748) + t748 * t766) * qJD(1)) * qJD(1);
t609 = -mrSges(2,2) * g(1) - mrSges(2,3) * t731 + Ifges(2,5) * qJDD(1) - t755 * Ifges(2,6) - qJ(2) * t620 + t748 * t612 - t746 * t614;
t1 = [(-m(1) - m(2)) * g(1) + t620; -m(1) * g(2) + t763; -m(1) * g(3) + t772; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t760; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t772 - t751 * t609 - t754 * t610; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t763 + t754 * t609 - t751 * t610; t760; t626; t762 * t781 - t789; t649; -t759;];
tauJB = t1;
