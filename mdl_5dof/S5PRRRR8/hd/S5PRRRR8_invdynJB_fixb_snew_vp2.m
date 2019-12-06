% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:15:42
% EndTime: 2019-12-05 17:15:50
% DurationCPUTime: 7.03s
% Computational Cost: add. (87721->269), mult. (171760->352), div. (0->0), fcn. (119898->12), ass. (0->119)
t689 = sin(qJ(4));
t690 = sin(qJ(3));
t693 = cos(qJ(4));
t694 = cos(qJ(3));
t661 = (t689 * t690 - t693 * t694) * qJD(2);
t684 = sin(pkin(10));
t686 = cos(pkin(10));
t671 = t684 * g(1) - t686 * g(2);
t672 = -t686 * g(1) - t684 * g(2);
t683 = -g(3) + qJDD(1);
t685 = sin(pkin(5));
t687 = cos(pkin(5));
t691 = sin(qJ(2));
t695 = cos(qJ(2));
t644 = -t691 * t672 + (t671 * t687 + t683 * t685) * t695;
t718 = t687 * t691;
t719 = t685 * t691;
t645 = t671 * t718 + t695 * t672 + t683 * t719;
t696 = qJD(2) ^ 2;
t640 = -t696 * pkin(2) + qJDD(2) * pkin(7) + t645;
t655 = -t685 * t671 + t687 * t683;
t626 = -t690 * t640 + t694 * t655;
t714 = qJD(2) * qJD(3);
t713 = t694 * t714;
t669 = t690 * qJDD(2) + t713;
t616 = (-t669 + t713) * pkin(8) + (t690 * t694 * t696 + qJDD(3)) * pkin(3) + t626;
t627 = t694 * t640 + t690 * t655;
t670 = t694 * qJDD(2) - t690 * t714;
t716 = qJD(2) * t690;
t676 = qJD(3) * pkin(3) - pkin(8) * t716;
t682 = t694 ^ 2;
t617 = -t682 * t696 * pkin(3) + t670 * pkin(8) - qJD(3) * t676 + t627;
t613 = t689 * t616 + t693 * t617;
t662 = (t689 * t694 + t690 * t693) * qJD(2);
t635 = -t662 * qJD(4) - t689 * t669 + t693 * t670;
t647 = t661 * mrSges(5,1) + t662 * mrSges(5,2);
t681 = qJD(3) + qJD(4);
t654 = t681 * mrSges(5,1) - t662 * mrSges(5,3);
t680 = qJDD(3) + qJDD(4);
t648 = t661 * pkin(4) - t662 * pkin(9);
t679 = t681 ^ 2;
t609 = -t679 * pkin(4) + t680 * pkin(9) - t661 * t648 + t613;
t702 = -qJDD(2) * pkin(2) - t644;
t621 = -t670 * pkin(3) + t676 * t716 + (-pkin(8) * t682 - pkin(7)) * t696 + t702;
t636 = -t661 * qJD(4) + t693 * t669 + t689 * t670;
t610 = (t661 * t681 - t636) * pkin(9) + (t662 * t681 - t635) * pkin(4) + t621;
t688 = sin(qJ(5));
t692 = cos(qJ(5));
t606 = -t688 * t609 + t692 * t610;
t649 = -t688 * t662 + t692 * t681;
t620 = t649 * qJD(5) + t692 * t636 + t688 * t680;
t650 = t692 * t662 + t688 * t681;
t628 = -t649 * mrSges(6,1) + t650 * mrSges(6,2);
t633 = qJDD(5) - t635;
t656 = qJD(5) + t661;
t637 = -t656 * mrSges(6,2) + t649 * mrSges(6,3);
t602 = m(6) * t606 + t633 * mrSges(6,1) - t620 * mrSges(6,3) - t650 * t628 + t656 * t637;
t607 = t692 * t609 + t688 * t610;
t619 = -t650 * qJD(5) - t688 * t636 + t692 * t680;
t638 = t656 * mrSges(6,1) - t650 * mrSges(6,3);
t603 = m(6) * t607 - t633 * mrSges(6,2) + t619 * mrSges(6,3) + t649 * t628 - t656 * t638;
t709 = -t688 * t602 + t692 * t603;
t589 = m(5) * t613 - t680 * mrSges(5,2) + t635 * mrSges(5,3) - t661 * t647 - t681 * t654 + t709;
t612 = t693 * t616 - t689 * t617;
t653 = -t681 * mrSges(5,2) - t661 * mrSges(5,3);
t608 = -t680 * pkin(4) - t679 * pkin(9) + t662 * t648 - t612;
t703 = -m(6) * t608 + t619 * mrSges(6,1) - t620 * mrSges(6,2) + t649 * t637 - t650 * t638;
t598 = m(5) * t612 + t680 * mrSges(5,1) - t636 * mrSges(5,3) - t662 * t647 + t681 * t653 + t703;
t583 = t689 * t589 + t693 * t598;
t659 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t690 + Ifges(4,2) * t694) * qJD(2);
t660 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t690 + Ifges(4,4) * t694) * qJD(2);
t622 = Ifges(6,5) * t650 + Ifges(6,6) * t649 + Ifges(6,3) * t656;
t624 = Ifges(6,1) * t650 + Ifges(6,4) * t649 + Ifges(6,5) * t656;
t595 = -mrSges(6,1) * t608 + mrSges(6,3) * t607 + Ifges(6,4) * t620 + Ifges(6,2) * t619 + Ifges(6,6) * t633 - t650 * t622 + t656 * t624;
t623 = Ifges(6,4) * t650 + Ifges(6,2) * t649 + Ifges(6,6) * t656;
t596 = mrSges(6,2) * t608 - mrSges(6,3) * t606 + Ifges(6,1) * t620 + Ifges(6,4) * t619 + Ifges(6,5) * t633 + t649 * t622 - t656 * t623;
t642 = Ifges(5,4) * t662 - Ifges(5,2) * t661 + Ifges(5,6) * t681;
t643 = Ifges(5,1) * t662 - Ifges(5,4) * t661 + Ifges(5,5) * t681;
t700 = -mrSges(5,1) * t612 + mrSges(5,2) * t613 - Ifges(5,5) * t636 - Ifges(5,6) * t635 - Ifges(5,3) * t680 - pkin(4) * t703 - pkin(9) * t709 - t692 * t595 - t688 * t596 - t662 * t642 - t661 * t643;
t721 = mrSges(4,1) * t626 - mrSges(4,2) * t627 + Ifges(4,5) * t669 + Ifges(4,6) * t670 + Ifges(4,3) * qJDD(3) + pkin(3) * t583 + (t690 * t659 - t694 * t660) * qJD(2) - t700;
t639 = -t696 * pkin(7) + t702;
t673 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t716;
t715 = qJD(2) * t694;
t674 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t715;
t591 = t692 * t602 + t688 * t603;
t701 = m(5) * t621 - t635 * mrSges(5,1) + t636 * mrSges(5,2) + t661 * t653 + t662 * t654 + t591;
t698 = -m(4) * t639 + t670 * mrSges(4,1) - t669 * mrSges(4,2) - t673 * t716 + t674 * t715 - t701;
t586 = m(3) * t644 + qJDD(2) * mrSges(3,1) - t696 * mrSges(3,2) + t698;
t720 = t586 * t695;
t668 = (-mrSges(4,1) * t694 + mrSges(4,2) * t690) * qJD(2);
t581 = m(4) * t626 + qJDD(3) * mrSges(4,1) - t669 * mrSges(4,3) + qJD(3) * t674 - t668 * t716 + t583;
t710 = t693 * t589 - t689 * t598;
t582 = m(4) * t627 - qJDD(3) * mrSges(4,2) + t670 * mrSges(4,3) - qJD(3) * t673 + t668 * t715 + t710;
t711 = -t690 * t581 + t694 * t582;
t572 = m(3) * t645 - t696 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t711;
t575 = t694 * t581 + t690 * t582;
t574 = m(3) * t655 + t575;
t562 = t572 * t718 - t685 * t574 + t687 * t720;
t560 = m(2) * t671 + t562;
t568 = t695 * t572 - t691 * t586;
t567 = m(2) * t672 + t568;
t717 = t686 * t560 + t684 * t567;
t561 = t572 * t719 + t687 * t574 + t685 * t720;
t712 = -t684 * t560 + t686 * t567;
t708 = m(2) * t683 + t561;
t641 = Ifges(5,5) * t662 - Ifges(5,6) * t661 + Ifges(5,3) * t681;
t576 = mrSges(5,2) * t621 - mrSges(5,3) * t612 + Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * t680 - pkin(9) * t591 - t688 * t595 + t692 * t596 - t661 * t641 - t681 * t642;
t699 = mrSges(6,1) * t606 - mrSges(6,2) * t607 + Ifges(6,5) * t620 + Ifges(6,6) * t619 + Ifges(6,3) * t633 + t650 * t623 - t649 * t624;
t577 = -mrSges(5,1) * t621 + mrSges(5,3) * t613 + Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * t680 - pkin(4) * t591 - t662 * t641 + t681 * t643 - t699;
t658 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t690 + Ifges(4,6) * t694) * qJD(2);
t563 = -mrSges(4,1) * t639 + mrSges(4,3) * t627 + Ifges(4,4) * t669 + Ifges(4,2) * t670 + Ifges(4,6) * qJDD(3) - pkin(3) * t701 + pkin(8) * t710 + qJD(3) * t660 + t689 * t576 + t693 * t577 - t658 * t716;
t564 = mrSges(4,2) * t639 - mrSges(4,3) * t626 + Ifges(4,1) * t669 + Ifges(4,4) * t670 + Ifges(4,5) * qJDD(3) - pkin(8) * t583 - qJD(3) * t659 + t693 * t576 - t689 * t577 + t658 * t715;
t557 = mrSges(3,2) * t655 - mrSges(3,3) * t644 + Ifges(3,5) * qJDD(2) - t696 * Ifges(3,6) - pkin(7) * t575 - t690 * t563 + t694 * t564;
t558 = -mrSges(3,1) * t655 + mrSges(3,3) * t645 + t696 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t575 - t721;
t704 = pkin(6) * t568 + t557 * t691 + t558 * t695;
t556 = mrSges(3,1) * t644 - mrSges(3,2) * t645 + Ifges(3,3) * qJDD(2) + pkin(2) * t698 + pkin(7) * t711 + t694 * t563 + t690 * t564;
t555 = mrSges(2,2) * t683 - mrSges(2,3) * t671 + t695 * t557 - t691 * t558 + (-t561 * t685 - t562 * t687) * pkin(6);
t554 = -mrSges(2,1) * t683 + mrSges(2,3) * t672 - pkin(1) * t561 - t685 * t556 + t704 * t687;
t1 = [-m(1) * g(1) + t712; -m(1) * g(2) + t717; -m(1) * g(3) + t708; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t717 - t684 * t554 + t686 * t555; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t712 + t686 * t554 + t684 * t555; -mrSges(1,1) * g(2) + mrSges(2,1) * t671 + mrSges(1,2) * g(1) - mrSges(2,2) * t672 + pkin(1) * t562 + t687 * t556 + t704 * t685; t708; t556; t721; -t700; t699;];
tauJB = t1;
