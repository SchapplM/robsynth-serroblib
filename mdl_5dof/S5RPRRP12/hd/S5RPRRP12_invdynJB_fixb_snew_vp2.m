% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:28
% EndTime: 2019-12-31 18:56:30
% DurationCPUTime: 1.63s
% Computational Cost: add. (14248->242), mult. (26952->279), div. (0->0), fcn. (14667->6), ass. (0->101)
t677 = Ifges(5,1) + Ifges(6,1);
t668 = Ifges(5,4) + Ifges(6,4);
t666 = Ifges(5,5) + Ifges(6,5);
t676 = Ifges(5,2) + Ifges(6,2);
t665 = Ifges(5,6) + Ifges(6,6);
t675 = Ifges(5,3) + Ifges(6,3);
t632 = sin(qJ(1));
t635 = cos(qJ(1));
t619 = -t635 * g(1) - t632 * g(2);
t674 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t619;
t630 = sin(qJ(4));
t633 = cos(qJ(4));
t634 = cos(qJ(3));
t658 = qJD(1) * t634;
t609 = t633 * qJD(3) - t630 * t658;
t631 = sin(qJ(3));
t656 = qJD(1) * qJD(3);
t652 = t631 * t656;
t614 = t634 * qJDD(1) - t652;
t576 = t609 * qJD(4) + t630 * qJDD(3) + t633 * t614;
t610 = t630 * qJD(3) + t633 * t658;
t579 = -t609 * mrSges(6,1) + t610 * mrSges(6,2);
t637 = qJD(1) ^ 2;
t672 = (-pkin(1) - pkin(6));
t590 = (t672 * t637) - t674;
t651 = t634 * t656;
t613 = -t631 * qJDD(1) - t651;
t559 = (-t614 + t652) * pkin(7) + (-t613 + t651) * pkin(3) + t590;
t618 = t632 * g(1) - t635 * g(2);
t645 = -t637 * qJ(2) + qJDD(2) - t618;
t591 = t672 * qJDD(1) + t645;
t582 = -t634 * g(3) + t631 * t591;
t612 = (pkin(3) * t631 - pkin(7) * t634) * qJD(1);
t636 = qJD(3) ^ 2;
t657 = t631 * qJD(1);
t562 = -t636 * pkin(3) + qJDD(3) * pkin(7) - t612 * t657 + t582;
t555 = t633 * t559 - t630 * t562;
t608 = qJDD(4) - t613;
t620 = qJD(4) + t657;
t551 = -0.2e1 * qJD(5) * t610 + (t609 * t620 - t576) * qJ(5) + (t609 * t610 + t608) * pkin(4) + t555;
t583 = -t620 * mrSges(6,2) + t609 * mrSges(6,3);
t654 = m(6) * t551 + t608 * mrSges(6,1) + t620 * t583;
t548 = -t576 * mrSges(6,3) - t610 * t579 + t654;
t556 = t630 * t559 + t633 * t562;
t575 = -t610 * qJD(4) + t633 * qJDD(3) - t630 * t614;
t585 = t620 * pkin(4) - t610 * qJ(5);
t607 = t609 ^ 2;
t553 = -t607 * pkin(4) + t575 * qJ(5) + 0.2e1 * qJD(5) * t609 - t620 * t585 + t556;
t660 = t668 * t609 + t677 * t610 + t666 * t620;
t661 = -t676 * t609 - t668 * t610 - t665 * t620;
t673 = mrSges(5,1) * t555 + mrSges(6,1) * t551 - mrSges(5,2) * t556 - mrSges(6,2) * t553 + pkin(4) * t548 + t665 * t575 + t666 * t576 + t675 * t608 - t660 * t609 - t661 * t610;
t671 = mrSges(2,1) - mrSges(3,2);
t670 = -mrSges(5,2) - mrSges(6,2);
t669 = -Ifges(3,4) + Ifges(2,5);
t667 = (Ifges(3,5) - Ifges(2,6));
t611 = (mrSges(4,1) * t631 + mrSges(4,2) * t634) * qJD(1);
t617 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t658;
t580 = -t609 * mrSges(5,1) + t610 * mrSges(5,2);
t584 = -t620 * mrSges(5,2) + t609 * mrSges(5,3);
t542 = m(5) * t555 + t608 * mrSges(5,1) + t620 * t584 + (-t579 - t580) * t610 + (-mrSges(5,3) - mrSges(6,3)) * t576 + t654;
t653 = m(6) * t553 + t575 * mrSges(6,3) + t609 * t579;
t586 = t620 * mrSges(6,1) - t610 * mrSges(6,3);
t659 = -t620 * mrSges(5,1) + t610 * mrSges(5,3) - t586;
t545 = m(5) * t556 + t575 * mrSges(5,3) + t609 * t580 + t670 * t608 + t659 * t620 + t653;
t648 = -t630 * t542 + t633 * t545;
t538 = m(4) * t582 - qJDD(3) * mrSges(4,2) + t613 * mrSges(4,3) - qJD(3) * t617 - t611 * t657 + t648;
t581 = t631 * g(3) + t634 * t591;
t616 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t657;
t561 = -qJDD(3) * pkin(3) - t636 * pkin(7) + t612 * t658 - t581;
t554 = -t575 * pkin(4) - t607 * qJ(5) + t610 * t585 + qJDD(5) + t561;
t647 = -m(6) * t554 + t575 * mrSges(6,1) + t609 * t583;
t639 = -m(5) * t561 + t575 * mrSges(5,1) + t670 * t576 + t609 * t584 + t659 * t610 + t647;
t546 = m(4) * t581 + qJDD(3) * mrSges(4,1) - t614 * mrSges(4,3) + qJD(3) * t616 - t611 * t658 + t639;
t528 = t631 * t538 + t634 * t546;
t596 = -qJDD(1) * pkin(1) + t645;
t644 = -m(3) * t596 + (t637 * mrSges(3,3)) - t528;
t524 = m(2) * t618 - (t637 * mrSges(2,2)) + t671 * qJDD(1) + t644;
t594 = t637 * pkin(1) + t674;
t540 = t633 * t542 + t630 * t545;
t643 = -m(4) * t590 + t613 * mrSges(4,1) - t614 * mrSges(4,2) - t616 * t657 - t617 * t658 - t540;
t640 = -m(3) * t594 + (t637 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t643;
t533 = m(2) * t619 - (t637 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t640;
t663 = t635 * t524 + t632 * t533;
t662 = -t665 * t609 - t666 * t610 - t675 * t620;
t650 = -t632 * t524 + t635 * t533;
t649 = t634 * t538 - t631 * t546;
t549 = t576 * mrSges(6,2) + t610 * t586 - t647;
t530 = -mrSges(5,1) * t561 + mrSges(5,3) * t556 - mrSges(6,1) * t554 + mrSges(6,3) * t553 - pkin(4) * t549 + qJ(5) * t653 + (-qJ(5) * t586 + t660) * t620 + t662 * t610 + (-qJ(5) * mrSges(6,2) + t665) * t608 + t668 * t576 + t676 * t575;
t536 = mrSges(5,2) * t561 + mrSges(6,2) * t554 - mrSges(5,3) * t555 - mrSges(6,3) * t551 - qJ(5) * t548 + t668 * t575 + t677 * t576 + t666 * t608 - t662 * t609 + t661 * t620;
t599 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t634 - Ifges(4,2) * t631) * qJD(1);
t600 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t634 - Ifges(4,4) * t631) * qJD(1);
t641 = mrSges(4,1) * t581 - mrSges(4,2) * t582 + Ifges(4,5) * t614 + Ifges(4,6) * t613 + Ifges(4,3) * qJDD(3) + pkin(3) * t639 + pkin(7) * t648 + t633 * t530 + t630 * t536 + t599 * t658 + t600 * t657;
t598 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t634 - Ifges(4,6) * t631) * qJD(1);
t521 = mrSges(4,2) * t590 - mrSges(4,3) * t581 + Ifges(4,1) * t614 + Ifges(4,4) * t613 + Ifges(4,5) * qJDD(3) - pkin(7) * t540 - qJD(3) * t599 - t630 * t530 + t633 * t536 - t598 * t657;
t522 = -mrSges(4,1) * t590 + mrSges(4,3) * t582 + Ifges(4,4) * t614 + Ifges(4,2) * t613 + Ifges(4,6) * qJDD(3) - pkin(3) * t540 + qJD(3) * t600 - t598 * t658 - t673;
t526 = qJDD(1) * mrSges(3,2) - t644;
t638 = mrSges(2,1) * t618 - mrSges(2,2) * t619 + mrSges(3,2) * t596 - mrSges(3,3) * t594 - pkin(1) * t526 - pkin(6) * t528 + qJ(2) * t640 + t634 * t521 - t631 * t522 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t527 = -m(3) * g(3) + t649;
t519 = (t667 * t637) + t641 - mrSges(2,3) * t618 + mrSges(3,1) * t596 + pkin(2) * t528 - qJ(2) * t527 + t669 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t518 = -mrSges(3,1) * t594 + mrSges(2,3) * t619 - pkin(1) * t527 - pkin(2) * t643 - pkin(6) * t649 + t671 * g(3) - t667 * qJDD(1) - t631 * t521 - t634 * t522 + t669 * t637;
t1 = [-m(1) * g(1) + t650; -m(1) * g(2) + t663; (-m(1) - m(2) - m(3)) * g(3) + t649; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t663 - t632 * t518 + t635 * t519; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t650 + t635 * t518 + t632 * t519; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t638; t638; t526; t641; t673; t549;];
tauJB = t1;
