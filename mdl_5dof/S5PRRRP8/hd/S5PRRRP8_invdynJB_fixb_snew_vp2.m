% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:42
% EndTime: 2019-12-05 16:58:49
% DurationCPUTime: 3.54s
% Computational Cost: add. (38819->243), mult. (72931->305), div. (0->0), fcn. (47248->10), ass. (0->109)
t697 = Ifges(5,1) + Ifges(6,1);
t690 = Ifges(5,4) - Ifges(6,5);
t689 = -Ifges(5,5) - Ifges(6,4);
t696 = Ifges(5,2) + Ifges(6,3);
t688 = Ifges(5,6) - Ifges(6,6);
t695 = -Ifges(5,3) - Ifges(6,2);
t652 = sin(pkin(9));
t654 = cos(pkin(9));
t641 = t652 * g(1) - t654 * g(2);
t642 = -t654 * g(1) - t652 * g(2);
t651 = -g(3) + qJDD(1);
t653 = sin(pkin(5));
t655 = cos(pkin(5));
t658 = sin(qJ(2));
t660 = cos(qJ(2));
t598 = -t658 * t642 + (t641 * t655 + t651 * t653) * t660;
t684 = t655 * t658;
t685 = t653 * t658;
t599 = t641 * t684 + t660 * t642 + t651 * t685;
t662 = qJD(2) ^ 2;
t597 = -t662 * pkin(2) + qJDD(2) * pkin(7) + t599;
t622 = -t653 * t641 + t655 * t651;
t657 = sin(qJ(3));
t659 = cos(qJ(3));
t592 = -t657 * t597 + t659 * t622;
t638 = (-pkin(3) * t659 - pkin(8) * t657) * qJD(2);
t661 = qJD(3) ^ 2;
t678 = qJD(2) * t657;
t588 = -qJDD(3) * pkin(3) - t661 * pkin(8) + t638 * t678 - t592;
t656 = sin(qJ(4));
t692 = cos(qJ(4));
t636 = t656 * qJD(3) + t692 * t678;
t676 = qJD(2) * qJD(3);
t673 = t659 * t676;
t639 = t657 * qJDD(2) + t673;
t610 = t636 * qJD(4) - t692 * qJDD(3) + t656 * t639;
t635 = -t692 * qJD(3) + t656 * t678;
t611 = -t635 * qJD(4) + t656 * qJDD(3) + t692 * t639;
t677 = t659 * qJD(2);
t647 = qJD(4) - t677;
t584 = -0.2e1 * qJD(5) * t636 + (t635 * t647 - t611) * qJ(5) + (t636 * t647 + t610) * pkin(4) + t588;
t620 = -t647 * mrSges(6,1) + t636 * mrSges(6,2);
t621 = -t635 * mrSges(6,2) + t647 * mrSges(6,3);
t580 = m(6) * t584 + t610 * mrSges(6,1) - t611 * mrSges(6,3) - t636 * t620 + t635 * t621;
t593 = t659 * t597 + t657 * t622;
t589 = -t661 * pkin(3) + qJDD(3) * pkin(8) + t638 * t677 + t593;
t596 = -qJDD(2) * pkin(2) - t662 * pkin(7) - t598;
t674 = t657 * t676;
t640 = t659 * qJDD(2) - t674;
t591 = (-t639 - t673) * pkin(8) + (-t640 + t674) * pkin(3) + t596;
t586 = t692 * t589 + t656 * t591;
t614 = t635 * pkin(4) - t636 * qJ(5);
t632 = qJDD(4) - t640;
t646 = t647 ^ 2;
t582 = -t646 * pkin(4) + t632 * qJ(5) + 0.2e1 * qJD(5) * t647 - t635 * t614 + t586;
t680 = t690 * t635 - t697 * t636 + t689 * t647;
t682 = t688 * t635 + t689 * t636 + t695 * t647;
t568 = -mrSges(5,1) * t588 - mrSges(6,1) * t584 + mrSges(6,2) * t582 + mrSges(5,3) * t586 - pkin(4) * t580 - t696 * t610 + t690 * t611 + t688 * t632 + t682 * t636 - t680 * t647;
t585 = -t656 * t589 + t692 * t591;
t583 = -t632 * pkin(4) - t646 * qJ(5) + t636 * t614 + qJDD(5) - t585;
t681 = t696 * t635 - t690 * t636 - t688 * t647;
t569 = mrSges(5,2) * t588 + mrSges(6,2) * t583 - mrSges(5,3) * t585 - mrSges(6,3) * t584 - qJ(5) * t580 - t690 * t610 + t697 * t611 - t689 * t632 + t682 * t635 + t681 * t647;
t619 = t647 * mrSges(5,1) - t636 * mrSges(5,3);
t675 = m(6) * t582 + t632 * mrSges(6,3) + t647 * t620;
t615 = t635 * mrSges(6,1) - t636 * mrSges(6,3);
t679 = -t635 * mrSges(5,1) - t636 * mrSges(5,2) - t615;
t691 = -mrSges(5,3) - mrSges(6,2);
t575 = m(5) * t586 - t632 * mrSges(5,2) + t691 * t610 - t647 * t619 + t679 * t635 + t675;
t618 = -t647 * mrSges(5,2) - t635 * mrSges(5,3);
t669 = -m(6) * t583 + t632 * mrSges(6,1) + t647 * t621;
t576 = m(5) * t585 + t632 * mrSges(5,1) + t691 * t611 + t647 * t618 + t679 * t636 + t669;
t571 = t692 * t575 - t656 * t576;
t577 = -m(5) * t588 - t610 * mrSges(5,1) - t611 * mrSges(5,2) - t635 * t618 - t636 * t619 - t580;
t626 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t657 + Ifges(4,2) * t659) * qJD(2);
t627 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t657 + Ifges(4,4) * t659) * qJD(2);
t694 = mrSges(4,1) * t592 - mrSges(4,2) * t593 + Ifges(4,5) * t639 + Ifges(4,6) * t640 + Ifges(4,3) * qJDD(3) + pkin(3) * t577 + pkin(8) * t571 + (t626 * t657 - t627 * t659) * qJD(2) + t692 * t568 + t656 * t569;
t579 = t611 * mrSges(6,2) + t636 * t615 - t669;
t693 = -t688 * t610 - t689 * t611 - t695 * t632 - t680 * t635 - t681 * t636 + mrSges(5,1) * t585 - mrSges(6,1) * t583 - mrSges(5,2) * t586 + mrSges(6,3) * t582 - pkin(4) * t579 + qJ(5) * (-t610 * mrSges(6,2) - t635 * t615 + t675);
t570 = t656 * t575 + t692 * t576;
t643 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t678;
t644 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t677;
t664 = -m(4) * t596 + t640 * mrSges(4,1) - t639 * mrSges(4,2) - t643 * t678 + t644 * t677 - t570;
t564 = m(3) * t598 + qJDD(2) * mrSges(3,1) - t662 * mrSges(3,2) + t664;
t686 = t564 * t660;
t637 = (-mrSges(4,1) * t659 + mrSges(4,2) * t657) * qJD(2);
t567 = m(4) * t593 - qJDD(3) * mrSges(4,2) + t640 * mrSges(4,3) - qJD(3) * t643 + t637 * t677 + t571;
t573 = m(4) * t592 + qJDD(3) * mrSges(4,1) - t639 * mrSges(4,3) + qJD(3) * t644 - t637 * t678 + t577;
t671 = t659 * t567 - t657 * t573;
t558 = m(3) * t599 - t662 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t671;
t561 = t657 * t567 + t659 * t573;
t560 = m(3) * t622 + t561;
t548 = t558 * t684 - t653 * t560 + t655 * t686;
t546 = m(2) * t641 + t548;
t553 = t660 * t558 - t658 * t564;
t552 = m(2) * t642 + t553;
t683 = t654 * t546 + t652 * t552;
t547 = t558 * t685 + t655 * t560 + t653 * t686;
t672 = -t652 * t546 + t654 * t552;
t670 = m(2) * t651 + t547;
t625 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t657 + Ifges(4,6) * t659) * qJD(2);
t549 = mrSges(4,2) * t596 - mrSges(4,3) * t592 + Ifges(4,1) * t639 + Ifges(4,4) * t640 + Ifges(4,5) * qJDD(3) - pkin(8) * t570 - qJD(3) * t626 - t656 * t568 + t692 * t569 + t625 * t677;
t554 = -mrSges(4,1) * t596 + mrSges(4,3) * t593 + Ifges(4,4) * t639 + Ifges(4,2) * t640 + Ifges(4,6) * qJDD(3) - pkin(3) * t570 + qJD(3) * t627 - t625 * t678 - t693;
t543 = mrSges(3,2) * t622 - mrSges(3,3) * t598 + Ifges(3,5) * qJDD(2) - t662 * Ifges(3,6) - pkin(7) * t561 + t659 * t549 - t657 * t554;
t544 = -mrSges(3,1) * t622 + mrSges(3,3) * t599 + t662 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t561 - t694;
t666 = pkin(6) * t553 + t543 * t658 + t544 * t660;
t542 = mrSges(3,1) * t598 - mrSges(3,2) * t599 + Ifges(3,3) * qJDD(2) + pkin(2) * t664 + pkin(7) * t671 + t657 * t549 + t659 * t554;
t541 = mrSges(2,2) * t651 - mrSges(2,3) * t641 + t660 * t543 - t658 * t544 + (-t547 * t653 - t548 * t655) * pkin(6);
t540 = -mrSges(2,1) * t651 + mrSges(2,3) * t642 - pkin(1) * t547 - t653 * t542 + t655 * t666;
t1 = [-m(1) * g(1) + t672; -m(1) * g(2) + t683; -m(1) * g(3) + t670; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t683 - t652 * t540 + t654 * t541; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t672 + t654 * t540 + t652 * t541; -mrSges(1,1) * g(2) + mrSges(2,1) * t641 + mrSges(1,2) * g(1) - mrSges(2,2) * t642 + pkin(1) * t548 + t655 * t542 + t653 * t666; t670; t542; t694; t693; t579;];
tauJB = t1;
