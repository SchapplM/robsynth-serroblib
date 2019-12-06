% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR5
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:26:12
% EndTime: 2019-12-05 16:26:21
% DurationCPUTime: 6.85s
% Computational Cost: add. (77385->268), mult. (163577->353), div. (0->0), fcn. (112699->12), ass. (0->118)
t667 = sin(pkin(10));
t670 = cos(pkin(10));
t674 = sin(qJ(3));
t677 = cos(qJ(3));
t644 = (t667 * t674 - t670 * t677) * qJD(2);
t668 = sin(pkin(9));
t671 = cos(pkin(9));
t658 = t668 * g(1) - t671 * g(2);
t659 = -t671 * g(1) - t668 * g(2);
t666 = -g(3) + qJDD(1);
t669 = sin(pkin(5));
t672 = cos(pkin(5));
t675 = sin(qJ(2));
t678 = cos(qJ(2));
t626 = -t675 * t659 + (t658 * t672 + t666 * t669) * t678;
t699 = t672 * t675;
t700 = t669 * t675;
t627 = t658 * t699 + t678 * t659 + t666 * t700;
t680 = qJD(2) ^ 2;
t622 = -t680 * pkin(2) + qJDD(2) * pkin(7) + t627;
t641 = -t669 * t658 + t672 * t666;
t608 = -t674 * t622 + t677 * t641;
t695 = qJD(2) * qJD(3);
t694 = t677 * t695;
t656 = t674 * qJDD(2) + t694;
t605 = (-t656 + t694) * qJ(4) + (t674 * t677 * t680 + qJDD(3)) * pkin(3) + t608;
t609 = t677 * t622 + t674 * t641;
t657 = t677 * qJDD(2) - t674 * t695;
t697 = qJD(2) * t674;
t660 = qJD(3) * pkin(3) - qJ(4) * t697;
t665 = t677 ^ 2;
t606 = -t665 * t680 * pkin(3) + t657 * qJ(4) - qJD(3) * t660 + t609;
t703 = 2 * qJD(4);
t601 = t667 * t605 + t670 * t606 - t644 * t703;
t645 = (t667 * t677 + t670 * t674) * qJD(2);
t630 = t644 * pkin(4) - t645 * pkin(8);
t679 = qJD(3) ^ 2;
t599 = -t679 * pkin(4) + qJDD(3) * pkin(8) - t644 * t630 + t601;
t684 = -qJDD(2) * pkin(2) - t626;
t607 = -t657 * pkin(3) + qJDD(4) + t660 * t697 + (-qJ(4) * t665 - pkin(7)) * t680 + t684;
t633 = -t667 * t656 + t670 * t657;
t634 = t670 * t656 + t667 * t657;
t602 = (qJD(3) * t644 - t634) * pkin(8) + (qJD(3) * t645 - t633) * pkin(4) + t607;
t673 = sin(qJ(5));
t676 = cos(qJ(5));
t596 = -t673 * t599 + t676 * t602;
t635 = t676 * qJD(3) - t673 * t645;
t616 = t635 * qJD(5) + t673 * qJDD(3) + t676 * t634;
t636 = t673 * qJD(3) + t676 * t645;
t617 = -t635 * mrSges(6,1) + t636 * mrSges(6,2);
t643 = qJD(5) + t644;
t619 = -t643 * mrSges(6,2) + t635 * mrSges(6,3);
t632 = qJDD(5) - t633;
t593 = m(6) * t596 + t632 * mrSges(6,1) - t616 * mrSges(6,3) - t636 * t617 + t643 * t619;
t597 = t676 * t599 + t673 * t602;
t615 = -t636 * qJD(5) + t676 * qJDD(3) - t673 * t634;
t620 = t643 * mrSges(6,1) - t636 * mrSges(6,3);
t594 = m(6) * t597 - t632 * mrSges(6,2) + t615 * mrSges(6,3) + t635 * t617 - t643 * t620;
t585 = -t673 * t593 + t676 * t594;
t629 = t644 * mrSges(5,1) + t645 * mrSges(5,2);
t640 = qJD(3) * mrSges(5,1) - t645 * mrSges(5,3);
t582 = m(5) * t601 - qJDD(3) * mrSges(5,2) + t633 * mrSges(5,3) - qJD(3) * t640 - t644 * t629 + t585;
t689 = -t670 * t605 + t667 * t606;
t598 = -qJDD(3) * pkin(4) - t679 * pkin(8) + (t703 + t630) * t645 + t689;
t595 = -m(6) * t598 + t615 * mrSges(6,1) - t616 * mrSges(6,2) + t635 * t619 - t636 * t620;
t600 = -0.2e1 * qJD(4) * t645 - t689;
t639 = -qJD(3) * mrSges(5,2) - t644 * mrSges(5,3);
t589 = m(5) * t600 + qJDD(3) * mrSges(5,1) - t634 * mrSges(5,3) + qJD(3) * t639 - t645 * t629 + t595;
t576 = t667 * t582 + t670 * t589;
t610 = Ifges(6,5) * t636 + Ifges(6,6) * t635 + Ifges(6,3) * t643;
t612 = Ifges(6,1) * t636 + Ifges(6,4) * t635 + Ifges(6,5) * t643;
t586 = -mrSges(6,1) * t598 + mrSges(6,3) * t597 + Ifges(6,4) * t616 + Ifges(6,2) * t615 + Ifges(6,6) * t632 - t636 * t610 + t643 * t612;
t611 = Ifges(6,4) * t636 + Ifges(6,2) * t635 + Ifges(6,6) * t643;
t587 = mrSges(6,2) * t598 - mrSges(6,3) * t596 + Ifges(6,1) * t616 + Ifges(6,4) * t615 + Ifges(6,5) * t632 + t635 * t610 - t643 * t611;
t624 = Ifges(5,4) * t645 - Ifges(5,2) * t644 + Ifges(5,6) * qJD(3);
t625 = Ifges(5,1) * t645 - Ifges(5,4) * t644 + Ifges(5,5) * qJD(3);
t648 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t674 + Ifges(4,2) * t677) * qJD(2);
t649 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t674 + Ifges(4,4) * t677) * qJD(2);
t704 = (t674 * t648 - t677 * t649) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t608 + mrSges(5,1) * t600 - mrSges(4,2) * t609 - mrSges(5,2) * t601 + Ifges(4,5) * t656 + Ifges(5,5) * t634 + Ifges(4,6) * t657 + Ifges(5,6) * t633 + pkin(3) * t576 + pkin(4) * t595 + pkin(8) * t585 + t676 * t586 + t673 * t587 + t645 * t624 + t644 * t625;
t584 = t676 * t593 + t673 * t594;
t583 = m(5) * t607 - t633 * mrSges(5,1) + t634 * mrSges(5,2) + t644 * t639 + t645 * t640 + t584;
t621 = -t680 * pkin(7) + t684;
t661 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t697;
t696 = qJD(2) * t677;
t662 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t696;
t682 = -m(4) * t621 + t657 * mrSges(4,1) - t656 * mrSges(4,2) - t661 * t697 + t662 * t696 - t583;
t579 = m(3) * t626 + qJDD(2) * mrSges(3,1) - t680 * mrSges(3,2) + t682;
t701 = t579 * t678;
t655 = (-mrSges(4,1) * t677 + mrSges(4,2) * t674) * qJD(2);
t574 = m(4) * t608 + qJDD(3) * mrSges(4,1) - t656 * mrSges(4,3) + qJD(3) * t662 - t655 * t697 + t576;
t691 = t670 * t582 - t667 * t589;
t575 = m(4) * t609 - qJDD(3) * mrSges(4,2) + t657 * mrSges(4,3) - qJD(3) * t661 + t655 * t696 + t691;
t692 = -t674 * t574 + t677 * t575;
t565 = m(3) * t627 - t680 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t692;
t568 = t677 * t574 + t674 * t575;
t567 = m(3) * t641 + t568;
t555 = t565 * t699 - t669 * t567 + t672 * t701;
t553 = m(2) * t658 + t555;
t561 = t678 * t565 - t675 * t579;
t560 = m(2) * t659 + t561;
t698 = t671 * t553 + t668 * t560;
t554 = t565 * t700 + t672 * t567 + t669 * t701;
t693 = -t668 * t553 + t671 * t560;
t690 = m(2) * t666 + t554;
t623 = Ifges(5,5) * t645 - Ifges(5,6) * t644 + Ifges(5,3) * qJD(3);
t569 = mrSges(5,2) * t607 - mrSges(5,3) * t600 + Ifges(5,1) * t634 + Ifges(5,4) * t633 + Ifges(5,5) * qJDD(3) - pkin(8) * t584 - qJD(3) * t624 - t673 * t586 + t676 * t587 - t644 * t623;
t683 = mrSges(6,1) * t596 - mrSges(6,2) * t597 + Ifges(6,5) * t616 + Ifges(6,6) * t615 + Ifges(6,3) * t632 + t636 * t611 - t635 * t612;
t570 = -mrSges(5,1) * t607 + mrSges(5,3) * t601 + Ifges(5,4) * t634 + Ifges(5,2) * t633 + Ifges(5,6) * qJDD(3) - pkin(4) * t584 + qJD(3) * t625 - t645 * t623 - t683;
t647 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t674 + Ifges(4,6) * t677) * qJD(2);
t556 = -mrSges(4,1) * t621 + mrSges(4,3) * t609 + Ifges(4,4) * t656 + Ifges(4,2) * t657 + Ifges(4,6) * qJDD(3) - pkin(3) * t583 + qJ(4) * t691 + qJD(3) * t649 + t667 * t569 + t670 * t570 - t647 * t697;
t557 = mrSges(4,2) * t621 - mrSges(4,3) * t608 + Ifges(4,1) * t656 + Ifges(4,4) * t657 + Ifges(4,5) * qJDD(3) - qJ(4) * t576 - qJD(3) * t648 + t670 * t569 - t667 * t570 + t647 * t696;
t550 = mrSges(3,2) * t641 - mrSges(3,3) * t626 + Ifges(3,5) * qJDD(2) - t680 * Ifges(3,6) - pkin(7) * t568 - t674 * t556 + t677 * t557;
t551 = -mrSges(3,1) * t641 + mrSges(3,3) * t627 + t680 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t568 - t704;
t685 = pkin(6) * t561 + t550 * t675 + t551 * t678;
t549 = mrSges(3,1) * t626 - mrSges(3,2) * t627 + Ifges(3,3) * qJDD(2) + pkin(2) * t682 + pkin(7) * t692 + t677 * t556 + t674 * t557;
t548 = mrSges(2,2) * t666 - mrSges(2,3) * t658 + t678 * t550 - t675 * t551 + (-t554 * t669 - t555 * t672) * pkin(6);
t547 = -mrSges(2,1) * t666 + mrSges(2,3) * t659 - pkin(1) * t554 - t669 * t549 + t672 * t685;
t1 = [-m(1) * g(1) + t693; -m(1) * g(2) + t698; -m(1) * g(3) + t690; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t698 - t668 * t547 + t671 * t548; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t693 + t671 * t547 + t668 * t548; -mrSges(1,1) * g(2) + mrSges(2,1) * t658 + mrSges(1,2) * g(1) - mrSges(2,2) * t659 + pkin(1) * t555 + t672 * t549 + t669 * t685; t690; t549; t704; t583; t683;];
tauJB = t1;
