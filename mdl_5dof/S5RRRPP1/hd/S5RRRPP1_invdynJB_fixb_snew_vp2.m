% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:39
% EndTime: 2019-12-31 20:49:42
% DurationCPUTime: 2.92s
% Computational Cost: add. (40754->249), mult. (54199->304), div. (0->0), fcn. (31035->8), ass. (0->104)
t692 = Ifges(5,1) + Ifges(6,1);
t685 = Ifges(5,4) - Ifges(6,5);
t684 = Ifges(5,5) + Ifges(6,4);
t691 = -Ifges(5,2) - Ifges(6,3);
t690 = -Ifges(6,2) - Ifges(5,3);
t683 = Ifges(5,6) - Ifges(6,6);
t644 = qJDD(1) + qJDD(2);
t651 = sin(qJ(3));
t654 = cos(qJ(3));
t645 = qJD(1) + qJD(2);
t673 = qJD(3) * t645;
t627 = t651 * t644 + t654 * t673;
t653 = sin(qJ(1));
t656 = cos(qJ(1));
t639 = t653 * g(1) - t656 * g(2);
t632 = qJDD(1) * pkin(1) + t639;
t640 = -t656 * g(1) - t653 * g(2);
t658 = qJD(1) ^ 2;
t633 = -t658 * pkin(1) + t640;
t652 = sin(qJ(2));
t655 = cos(qJ(2));
t609 = t652 * t632 + t655 * t633;
t643 = t645 ^ 2;
t601 = -t643 * pkin(2) + t644 * pkin(7) + t609;
t679 = t651 * t601;
t687 = pkin(3) * t643;
t582 = qJDD(3) * pkin(3) - t627 * qJ(4) - t679 + (qJ(4) * t673 + t651 * t687 - g(3)) * t654;
t586 = -t651 * g(3) + t654 * t601;
t628 = t654 * t644 - t651 * t673;
t681 = t645 * t651;
t634 = qJD(3) * pkin(3) - qJ(4) * t681;
t649 = t654 ^ 2;
t583 = t628 * qJ(4) - qJD(3) * t634 - t649 * t687 + t586;
t650 = sin(pkin(8));
t680 = t645 * t654;
t682 = cos(pkin(8));
t617 = t650 * t681 - t682 * t680;
t688 = -2 * qJD(4);
t579 = t650 * t582 + t682 * t583 + t617 * t688;
t605 = t650 * t627 - t682 * t628;
t618 = (t650 * t654 + t682 * t651) * t645;
t613 = qJD(3) * mrSges(5,1) - t618 * mrSges(5,3);
t597 = t617 * pkin(4) - t618 * qJ(5);
t657 = qJD(3) ^ 2;
t574 = -t657 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t617 * t597 + t579;
t614 = -qJD(3) * mrSges(6,1) + t618 * mrSges(6,2);
t671 = m(6) * t574 + qJDD(3) * mrSges(6,3) + qJD(3) * t614;
t598 = t617 * mrSges(6,1) - t618 * mrSges(6,3);
t674 = -t617 * mrSges(5,1) - t618 * mrSges(5,2) - t598;
t686 = -mrSges(5,3) - mrSges(6,2);
t567 = m(5) * t579 - qJDD(3) * mrSges(5,2) - qJD(3) * t613 + t686 * t605 + t674 * t617 + t671;
t663 = t682 * t582 - t650 * t583;
t578 = t618 * t688 + t663;
t606 = t682 * t627 + t650 * t628;
t612 = -qJD(3) * mrSges(5,2) - t617 * mrSges(5,3);
t575 = -qJDD(3) * pkin(4) - t657 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t597) * t618 - t663;
t615 = -t617 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t666 = -m(6) * t575 + qJDD(3) * mrSges(6,1) + qJD(3) * t615;
t568 = m(5) * t578 + qJDD(3) * mrSges(5,1) + qJD(3) * t612 + t686 * t606 + t674 * t618 + t666;
t560 = t650 * t567 + t682 * t568;
t571 = t606 * mrSges(6,2) + t618 * t598 - t666;
t585 = -t654 * g(3) - t679;
t621 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t651 + Ifges(4,2) * t654) * t645;
t622 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t651 + Ifges(4,4) * t654) * t645;
t675 = t684 * qJD(3) - t685 * t617 + t692 * t618;
t676 = t683 * qJD(3) + t691 * t617 + t685 * t618;
t689 = (Ifges(4,3) - t690) * qJDD(3) - t683 * t605 + t684 * t606 + t675 * t617 + t676 * t618 + (t651 * t621 - t654 * t622) * t645 + mrSges(4,1) * t585 + mrSges(5,1) * t578 - mrSges(6,1) * t575 - mrSges(4,2) * t586 - mrSges(5,2) * t579 + mrSges(6,3) * t574 + Ifges(4,5) * t627 + Ifges(4,6) * t628 + pkin(3) * t560 - pkin(4) * t571 + qJ(5) * (-t605 * mrSges(6,2) - t617 * t598 + t671);
t626 = (-mrSges(4,1) * t654 + mrSges(4,2) * t651) * t645;
t636 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t680;
t556 = m(4) * t585 + qJDD(3) * mrSges(4,1) - t627 * mrSges(4,3) + qJD(3) * t636 - t626 * t681 + t560;
t635 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t681;
t667 = t682 * t567 - t650 * t568;
t557 = m(4) * t586 - qJDD(3) * mrSges(4,2) + t628 * mrSges(4,3) - qJD(3) * t635 + t626 * t680 + t667;
t668 = -t651 * t556 + t654 * t557;
t550 = m(3) * t609 - t643 * mrSges(3,1) - t644 * mrSges(3,2) + t668;
t608 = t655 * t632 - t652 * t633;
t664 = -t644 * pkin(2) - t608;
t584 = -t628 * pkin(3) + qJDD(4) + t634 * t681 + (-qJ(4) * t649 - pkin(7)) * t643 + t664;
t577 = -0.2e1 * qJD(5) * t618 + (qJD(3) * t617 - t606) * qJ(5) + (qJD(3) * t618 + t605) * pkin(4) + t584;
t572 = m(6) * t577 + t605 * mrSges(6,1) - t606 * mrSges(6,3) - t618 * t614 + t617 * t615;
t569 = m(5) * t584 + t605 * mrSges(5,1) + t606 * mrSges(5,2) + t617 * t612 + t618 * t613 + t572;
t600 = -t643 * pkin(7) + t664;
t660 = -m(4) * t600 + t628 * mrSges(4,1) - t627 * mrSges(4,2) - t635 * t681 + t636 * t680 - t569;
t562 = m(3) * t608 + t644 * mrSges(3,1) - t643 * mrSges(3,2) + t660;
t547 = t652 * t550 + t655 * t562;
t544 = m(2) * t639 + qJDD(1) * mrSges(2,1) - t658 * mrSges(2,2) + t547;
t669 = t655 * t550 - t652 * t562;
t545 = m(2) * t640 - t658 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t669;
t678 = t656 * t544 + t653 * t545;
t552 = t654 * t556 + t651 * t557;
t677 = t690 * qJD(3) + t683 * t617 - t684 * t618;
t670 = -t653 * t544 + t656 * t545;
t558 = -mrSges(5,1) * t584 - mrSges(6,1) * t577 + mrSges(6,2) * t574 + mrSges(5,3) * t579 - pkin(4) * t572 + t675 * qJD(3) + t683 * qJDD(3) + t691 * t605 + t685 * t606 + t677 * t618;
t559 = mrSges(5,2) * t584 + mrSges(6,2) * t575 - mrSges(5,3) * t578 - mrSges(6,3) * t577 - qJ(5) * t572 - t676 * qJD(3) + t684 * qJDD(3) - t685 * t605 + t692 * t606 + t677 * t617;
t620 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t651 + Ifges(4,6) * t654) * t645;
t538 = -mrSges(4,1) * t600 + mrSges(4,3) * t586 + Ifges(4,4) * t627 + Ifges(4,2) * t628 + Ifges(4,6) * qJDD(3) - pkin(3) * t569 + qJ(4) * t667 + qJD(3) * t622 + t682 * t558 + t650 * t559 - t620 * t681;
t540 = mrSges(4,2) * t600 - mrSges(4,3) * t585 + Ifges(4,1) * t627 + Ifges(4,4) * t628 + Ifges(4,5) * qJDD(3) - qJ(4) * t560 - qJD(3) * t621 - t650 * t558 + t682 * t559 + t620 * t680;
t662 = mrSges(3,1) * t608 - mrSges(3,2) * t609 + Ifges(3,3) * t644 + pkin(2) * t660 + pkin(7) * t668 + t654 * t538 + t651 * t540;
t661 = mrSges(2,1) * t639 - mrSges(2,2) * t640 + Ifges(2,3) * qJDD(1) + pkin(1) * t547 + t662;
t536 = mrSges(3,1) * g(3) + mrSges(3,3) * t609 + t643 * Ifges(3,5) + Ifges(3,6) * t644 - pkin(2) * t552 - t689;
t535 = -mrSges(3,2) * g(3) - mrSges(3,3) * t608 + Ifges(3,5) * t644 - t643 * Ifges(3,6) - pkin(7) * t552 - t651 * t538 + t654 * t540;
t534 = -mrSges(2,2) * g(3) - mrSges(2,3) * t639 + Ifges(2,5) * qJDD(1) - t658 * Ifges(2,6) - pkin(6) * t547 + t655 * t535 - t652 * t536;
t533 = Ifges(2,6) * qJDD(1) + t658 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t640 + t652 * t535 + t655 * t536 - pkin(1) * (-m(3) * g(3) + t552) + pkin(6) * t669;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t678; (-m(1) - m(2) - m(3)) * g(3) + t552; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t678 - t653 * t533 + t656 * t534; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t670 + t656 * t533 + t653 * t534; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t661; t661; t662; t689; t569; t571;];
tauJB = t1;
