% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:17
% EndTime: 2019-12-31 20:53:19
% DurationCPUTime: 1.48s
% Computational Cost: add. (14825->229), mult. (18638->265), div. (0->0), fcn. (8001->6), ass. (0->96)
t695 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t679 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t678 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t694 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t677 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t693 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t642 = qJD(1) + qJD(2);
t651 = cos(qJ(3));
t683 = t642 * t651;
t628 = -mrSges(5,1) * t683 - qJD(3) * mrSges(5,3);
t650 = sin(qJ(1));
t653 = cos(qJ(1));
t632 = t650 * g(1) - t653 * g(2);
t622 = qJDD(1) * pkin(1) + t632;
t633 = -t653 * g(1) - t650 * g(2);
t655 = qJD(1) ^ 2;
t623 = -t655 * pkin(1) + t633;
t649 = sin(qJ(2));
t652 = cos(qJ(2));
t588 = t649 * t622 + t652 * t623;
t640 = t642 ^ 2;
t641 = qJDD(1) + qJDD(2);
t585 = -t640 * pkin(2) + t641 * pkin(7) + t588;
t648 = sin(qJ(3));
t580 = -t651 * g(3) - t648 * t585;
t610 = (-pkin(3) * t651 - qJ(4) * t648) * t642;
t654 = qJD(3) ^ 2;
t684 = t642 * t648;
t579 = -qJDD(3) * pkin(3) - t654 * qJ(4) + t610 * t684 + qJDD(4) - t580;
t680 = qJD(3) * t642;
t671 = t651 * t680;
t614 = t648 * t641 + t671;
t687 = -2 * qJD(5);
t575 = qJD(3) * t687 + (-t640 * t648 * t651 - qJDD(3)) * qJ(5) + (t614 - t671) * pkin(4) + t579;
t629 = mrSges(6,1) * t683 + qJD(3) * mrSges(6,2);
t666 = -m(6) * t575 + qJDD(3) * mrSges(6,3) + qJD(3) * t629;
t663 = m(5) * t579 + t614 * mrSges(5,1) - t666;
t611 = (mrSges(5,2) * t651 - mrSges(5,3) * t648) * t642;
t613 = (-mrSges(6,2) * t648 - mrSges(6,3) * t651) * t642;
t681 = t611 + t613;
t685 = t614 * mrSges(6,1);
t569 = qJDD(3) * mrSges(5,2) + qJD(3) * t628 + t681 * t684 + t663 + t685;
t570 = t613 * t684 - t666 + t685;
t672 = t648 * t680;
t615 = t651 * t641 - t672;
t626 = pkin(4) * t684 - qJD(3) * qJ(5);
t647 = t651 ^ 2;
t581 = -t648 * g(3) + t651 * t585;
t659 = -t654 * pkin(3) + qJDD(3) * qJ(4) + t610 * t683 + t581;
t576 = -t647 * t640 * qJ(5) + t615 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t626) * qJD(3) + t659;
t688 = -2 * qJD(4);
t578 = qJD(3) * t688 - t659;
t630 = mrSges(5,1) * t684 + qJD(3) * mrSges(5,2);
t627 = mrSges(6,1) * t684 - qJD(3) * mrSges(6,3);
t667 = m(6) * t576 + qJDD(3) * mrSges(6,2) + qJD(3) * t627 + t613 * t683;
t660 = -m(5) * t578 + qJDD(3) * mrSges(5,3) + qJD(3) * t630 + t611 * t683 + t667;
t673 = (t695 * t648 + t679 * t651) * t642 + t678 * qJD(3);
t674 = (t679 * t648 + t694 * t651) * t642 + t677 * qJD(3);
t692 = (t674 * t648 - t673 * t651) * t642 + t693 * qJDD(3) + t678 * t614 + t677 * t615 + mrSges(4,1) * t580 - mrSges(4,2) * t581 + mrSges(5,2) * t579 + mrSges(6,2) * t576 - mrSges(5,3) * t578 - mrSges(6,3) * t575 - pkin(3) * t569 + qJ(4) * ((mrSges(5,1) + mrSges(6,1)) * t615 + t660) - qJ(5) * t570;
t691 = pkin(3) * t672 + t684 * t688;
t686 = -mrSges(6,1) - mrSges(4,3);
t612 = (-mrSges(4,1) * t651 + mrSges(4,2) * t648) * t642;
t624 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t684;
t565 = t612 * t683 + m(4) * t581 - qJDD(3) * mrSges(4,2) - qJD(3) * t624 + (mrSges(5,1) - t686) * t615 + t660;
t625 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t683;
t566 = m(4) * t580 + t686 * t614 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t625 - t628) * qJD(3) + (-t612 - t681) * t684 - t663;
t668 = t651 * t565 - t648 * t566;
t556 = m(3) * t588 - t640 * mrSges(3,1) - t641 * mrSges(3,2) + t668;
t587 = t652 * t622 - t649 * t623;
t664 = -t641 * pkin(2) - t587;
t584 = -t640 * pkin(7) + t664;
t577 = -t615 * pkin(3) + (-t614 - t671) * qJ(4) + t584 + t691;
t573 = -t614 * qJ(4) + (-pkin(4) * t647 - pkin(7)) * t640 + (-pkin(3) - qJ(5)) * t615 + (-t626 * t648 + (-qJ(4) * qJD(3) + t687) * t651) * t642 + t664 + t691;
t665 = m(6) * t573 - t614 * mrSges(6,2) - t615 * mrSges(6,3) - t627 * t684 - t629 * t683;
t661 = -m(5) * t577 - t615 * mrSges(5,2) + t630 * t684 - t665;
t656 = -m(4) * t584 + t625 * t683 + t615 * mrSges(4,1) + (-t624 * t648 - t628 * t651) * t642 + (-mrSges(4,2) + mrSges(5,3)) * t614 + t661;
t560 = m(3) * t587 + t641 * mrSges(3,1) - t640 * mrSges(3,2) + t656;
t551 = t649 * t556 + t652 * t560;
t548 = m(2) * t632 + qJDD(1) * mrSges(2,1) - t655 * mrSges(2,2) + t551;
t669 = t652 * t556 - t649 * t560;
t549 = m(2) * t633 - t655 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t669;
t682 = t653 * t548 + t650 * t549;
t558 = t648 * t565 + t651 * t566;
t675 = (t678 * t648 + t677 * t651) * t642 + t693 * qJD(3);
t670 = -t650 * t548 + t653 * t549;
t567 = -t614 * mrSges(5,3) + t628 * t683 - t661;
t571 = t615 * mrSges(6,1) + t667;
t544 = -mrSges(4,1) * t584 - mrSges(5,1) * t578 + mrSges(6,1) * t576 + mrSges(5,2) * t577 + mrSges(4,3) * t581 - mrSges(6,3) * t573 - pkin(3) * t567 + pkin(4) * t571 - qJ(5) * t665 + t673 * qJD(3) + t677 * qJDD(3) + t679 * t614 + t694 * t615 - t675 * t684;
t553 = mrSges(5,1) * t579 + mrSges(6,1) * t575 + mrSges(4,2) * t584 - mrSges(6,2) * t573 - mrSges(4,3) * t580 - mrSges(5,3) * t577 + pkin(4) * t570 - qJ(4) * t567 - t674 * qJD(3) + t678 * qJDD(3) + t695 * t614 + t679 * t615 + t675 * t683;
t662 = mrSges(3,1) * t587 - mrSges(3,2) * t588 + Ifges(3,3) * t641 + pkin(2) * t656 + pkin(7) * t668 + t651 * t544 + t648 * t553;
t658 = mrSges(2,1) * t632 - mrSges(2,2) * t633 + Ifges(2,3) * qJDD(1) + pkin(1) * t551 + t662;
t542 = mrSges(3,1) * g(3) + mrSges(3,3) * t588 + t640 * Ifges(3,5) + Ifges(3,6) * t641 - pkin(2) * t558 - t692;
t541 = -mrSges(3,2) * g(3) - mrSges(3,3) * t587 + Ifges(3,5) * t641 - t640 * Ifges(3,6) - pkin(7) * t558 - t648 * t544 + t651 * t553;
t540 = -mrSges(2,2) * g(3) - mrSges(2,3) * t632 + Ifges(2,5) * qJDD(1) - t655 * Ifges(2,6) - pkin(6) * t551 + t652 * t541 - t649 * t542;
t539 = Ifges(2,6) * qJDD(1) + t655 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t633 + t649 * t541 + t652 * t542 - pkin(1) * (-m(3) * g(3) + t558) + pkin(6) * t669;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t682; (-m(1) - m(2) - m(3)) * g(3) + t558; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t682 - t650 * t539 + t653 * t540; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t670 + t653 * t539 + t650 * t540; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t658; t658; t662; t692; t569; t571;];
tauJB = t1;
