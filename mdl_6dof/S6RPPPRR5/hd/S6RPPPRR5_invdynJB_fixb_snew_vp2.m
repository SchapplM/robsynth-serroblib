% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:19
% EndTime: 2019-05-05 13:50:21
% DurationCPUTime: 2.08s
% Computational Cost: add. (22999->251), mult. (38662->288), div. (0->0), fcn. (15847->8), ass. (0->106)
t658 = -2 * qJD(1);
t621 = sin(qJ(1));
t624 = cos(qJ(1));
t598 = -t624 * g(1) - t621 * g(2);
t657 = -qJDD(1) * qJ(2) + (qJD(2) * t658) - t598;
t626 = qJD(1) ^ 2;
t597 = t621 * g(1) - t624 * g(2);
t636 = qJDD(2) - t597;
t653 = -pkin(1) - qJ(3);
t631 = (qJD(3) * t658) + t653 * qJDD(1) + t636;
t571 = (-pkin(3) - qJ(2)) * t626 + t631;
t576 = t653 * t626 + qJDD(3) - t657;
t572 = qJDD(1) * pkin(3) + t576;
t617 = sin(pkin(9));
t618 = cos(pkin(9));
t558 = -t617 * t571 + t618 * t572;
t555 = -qJDD(1) * pkin(4) - (t626 * pkin(7)) - t558;
t620 = sin(qJ(5));
t623 = cos(qJ(5));
t645 = qJD(1) * qJD(5);
t640 = t623 * t645;
t592 = t620 * qJDD(1) + t640;
t641 = t620 * t645;
t593 = t623 * qJDD(1) - t641;
t549 = (-t592 - t640) * pkin(8) + (-t593 + t641) * pkin(5) + t555;
t559 = t618 * t571 + t617 * t572;
t556 = -t626 * pkin(4) + qJDD(1) * pkin(7) + t559;
t613 = -g(3) + qJDD(4);
t553 = t623 * t556 + t620 * t613;
t591 = (-pkin(5) * t623 - pkin(8) * t620) * qJD(1);
t625 = qJD(5) ^ 2;
t646 = t623 * qJD(1);
t551 = -t625 * pkin(5) + qJDD(5) * pkin(8) + t591 * t646 + t553;
t619 = sin(qJ(6));
t622 = cos(qJ(6));
t547 = t622 * t549 - t619 * t551;
t647 = qJD(1) * t620;
t588 = t622 * qJD(5) - t619 * t647;
t566 = t588 * qJD(6) + t619 * qJDD(5) + t622 * t592;
t589 = t619 * qJD(5) + t622 * t647;
t570 = -t588 * mrSges(7,1) + t589 * mrSges(7,2);
t599 = qJD(6) - t646;
t577 = -t599 * mrSges(7,2) + t588 * mrSges(7,3);
t587 = qJDD(6) - t593;
t544 = m(7) * t547 + t587 * mrSges(7,1) - t566 * mrSges(7,3) - t589 * t570 + t599 * t577;
t548 = t619 * t549 + t622 * t551;
t565 = -t589 * qJD(6) + t622 * qJDD(5) - t619 * t592;
t578 = t599 * mrSges(7,1) - t589 * mrSges(7,3);
t545 = m(7) * t548 - t587 * mrSges(7,2) + t565 * mrSges(7,3) + t588 * t570 - t599 * t578;
t538 = -t619 * t544 + t622 * t545;
t651 = t623 * t613;
t550 = -qJDD(5) * pkin(5) - t625 * pkin(8) - t651 + (qJD(1) * t591 + t556) * t620;
t560 = Ifges(7,5) * t589 + Ifges(7,6) * t588 + Ifges(7,3) * t599;
t562 = Ifges(7,1) * t589 + Ifges(7,4) * t588 + Ifges(7,5) * t599;
t539 = -mrSges(7,1) * t550 + mrSges(7,3) * t548 + Ifges(7,4) * t566 + Ifges(7,2) * t565 + Ifges(7,6) * t587 - t589 * t560 + t599 * t562;
t561 = Ifges(7,4) * t589 + Ifges(7,2) * t588 + Ifges(7,6) * t599;
t540 = mrSges(7,2) * t550 - mrSges(7,3) * t547 + Ifges(7,1) * t566 + Ifges(7,4) * t565 + Ifges(7,5) * t587 + t588 * t560 - t599 * t561;
t546 = -m(7) * t550 + t565 * mrSges(7,1) - t566 * mrSges(7,2) + t588 * t577 - t589 * t578;
t552 = -t620 * t556 + t651;
t583 = (Ifges(6,6) * qJD(5)) + (Ifges(6,4) * t620 + Ifges(6,2) * t623) * qJD(1);
t584 = (Ifges(6,5) * qJD(5)) + (Ifges(6,1) * t620 + Ifges(6,4) * t623) * qJD(1);
t656 = mrSges(6,1) * t552 - mrSges(6,2) * t553 + Ifges(6,5) * t592 + Ifges(6,6) * t593 + Ifges(6,3) * qJDD(5) + pkin(5) * t546 + pkin(8) * t538 + t622 * t539 + t619 * t540 + (t620 * t583 - t623 * t584) * qJD(1);
t655 = -m(3) - m(4);
t654 = mrSges(2,1) - mrSges(3,2);
t652 = t626 * mrSges(4,3);
t650 = t626 * qJ(2);
t579 = (t626 * pkin(1)) + t657;
t590 = (-mrSges(6,1) * t623 + mrSges(6,2) * t620) * qJD(1);
t595 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t647;
t536 = m(6) * t553 - qJDD(5) * mrSges(6,2) + t593 * mrSges(6,3) - qJD(5) * t595 + t590 * t646 + t538;
t596 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t646;
t542 = m(6) * t552 + qJDD(5) * mrSges(6,1) - t592 * mrSges(6,3) + qJD(5) * t596 - t590 * t647 + t546;
t638 = t623 * t536 - t620 * t542;
t525 = m(5) * t559 - (t626 * mrSges(5,1)) - qJDD(1) * mrSges(5,2) + t638;
t537 = t622 * t544 + t619 * t545;
t630 = -m(6) * t555 + t593 * mrSges(6,1) - t592 * mrSges(6,2) - t595 * t647 + t596 * t646 - t537;
t532 = m(5) * t558 + qJDD(1) * mrSges(5,1) - t626 * mrSges(5,2) + t630;
t518 = t617 * t525 + t618 * t532;
t637 = m(4) * t576 + qJDD(1) * mrSges(4,1) + t518;
t633 = -m(3) * t579 + (t626 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t637;
t510 = m(2) * t598 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t626) + t633;
t575 = t631 - t650;
t648 = t618 * t525 - t617 * t532;
t515 = m(4) * t575 - t626 * mrSges(4,1) - qJDD(1) * mrSges(4,3) + t648;
t581 = -qJDD(1) * pkin(1) + t636 - t650;
t632 = -m(3) * t581 + t626 * mrSges(3,3) - t515;
t511 = m(2) * t597 - t626 * mrSges(2,2) + t654 * qJDD(1) + t632;
t649 = t621 * t510 + t624 * t511;
t529 = t620 * t536 + t623 * t542;
t643 = -Ifges(3,4) + Ifges(2,5) - Ifges(4,6);
t642 = (Ifges(4,4) - Ifges(3,5) + Ifges(2,6));
t527 = m(5) * t613 + t529;
t639 = t624 * t510 - t621 * t511;
t629 = mrSges(7,1) * t547 - mrSges(7,2) * t548 + Ifges(7,5) * t566 + Ifges(7,6) * t565 + Ifges(7,3) * t587 + t589 * t561 - t588 * t562;
t514 = qJDD(1) * mrSges(3,2) - t632;
t582 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t620 + Ifges(6,6) * t623) * qJD(1);
t520 = mrSges(6,2) * t555 - mrSges(6,3) * t552 + Ifges(6,1) * t592 + Ifges(6,4) * t593 + Ifges(6,5) * qJDD(5) - pkin(8) * t537 - qJD(5) * t583 - t619 * t539 + t622 * t540 + t582 * t646;
t522 = -mrSges(6,1) * t555 + mrSges(6,3) * t553 + Ifges(6,4) * t592 + Ifges(6,2) * t593 + Ifges(6,6) * qJDD(5) - pkin(5) * t537 + qJD(5) * t584 - t582 * t647 - t629;
t627 = -mrSges(2,2) * t598 - mrSges(5,2) * t559 - mrSges(3,3) * t579 - mrSges(4,3) * t575 - qJ(3) * t515 + qJ(2) * (t633 - t652) - pkin(1) * t514 + pkin(3) * t518 + t620 * t520 + t623 * t522 + pkin(7) * t638 + pkin(4) * t630 + mrSges(5,1) * t558 + mrSges(4,1) * t576 + mrSges(3,2) * t581 + mrSges(2,1) * t597 + (Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + Ifges(5,3)) * qJDD(1);
t526 = t655 * g(3) + t527;
t516 = t637 - t652;
t506 = -mrSges(5,1) * t613 + mrSges(5,3) * t559 + (t626 * Ifges(5,5)) + Ifges(5,6) * qJDD(1) - pkin(4) * t529 - t656;
t505 = mrSges(5,2) * t613 - mrSges(5,3) * t558 + Ifges(5,5) * qJDD(1) - t626 * Ifges(5,6) - pkin(7) * t529 + t623 * t520 - t620 * t522;
t504 = -qJ(2) * t526 - mrSges(2,3) * t597 + pkin(2) * t515 + mrSges(3,1) * t581 + t617 * t505 + t618 * t506 - pkin(3) * t527 + qJ(4) * t648 - mrSges(4,2) * t575 - (t642 * t626) + t643 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,1)) * g(3);
t503 = -pkin(1) * t526 + mrSges(2,3) * t598 + pkin(2) * t516 - qJ(3) * t527 - t618 * t505 + t617 * t506 + qJ(4) * t518 - mrSges(3,1) * t579 - mrSges(4,2) * t576 + t643 * t626 + t642 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t654) * g(3);
t1 = [-m(1) * g(1) + t639; -m(1) * g(2) + t649; (-m(1) - m(2) + t655) * g(3) + t527; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t649 - t621 * t503 + t624 * t504; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t639 + t624 * t503 + t621 * t504; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t627; t627; t514; t516; t527; t656; t629;];
tauJB  = t1;
