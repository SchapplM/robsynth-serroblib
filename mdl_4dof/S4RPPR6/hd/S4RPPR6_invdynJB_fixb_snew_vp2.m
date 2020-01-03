% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:42
% EndTime: 2019-12-31 16:40:43
% DurationCPUTime: 0.97s
% Computational Cost: add. (5230->182), mult. (12074->225), div. (0->0), fcn. (6755->6), ass. (0->84)
t523 = sin(qJ(1));
t525 = cos(qJ(1));
t502 = -t525 * g(1) - t523 * g(2);
t526 = qJD(1) ^ 2;
t563 = -t526 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t502;
t520 = sin(pkin(6));
t521 = cos(pkin(6));
t562 = (Ifges(3,6) - Ifges(4,6)) * t521 + (Ifges(4,4) + Ifges(3,5)) * t520;
t516 = t520 ^ 2;
t517 = t521 ^ 2;
t551 = t517 * t526;
t559 = t516 * t526 + t551;
t556 = Ifges(3,4) - Ifges(4,5);
t558 = t556 * t520;
t482 = -t521 * g(3) - t563 * t520;
t557 = Ifges(3,1) + Ifges(4,1);
t555 = Ifges(3,2) + Ifges(4,3);
t554 = mrSges(3,2) * t520;
t553 = qJ(3) * t520;
t550 = t526 * qJ(2);
t496 = (-mrSges(4,1) * t521 - mrSges(4,3) * t520) * qJD(1);
t497 = (-mrSges(3,1) * t521 + t554) * qJD(1);
t536 = -pkin(2) * t521 - t553;
t495 = t536 * qJD(1);
t546 = t520 * qJD(1);
t466 = t495 * t546 + qJDD(3) - t482;
t461 = (-pkin(3) * t521 * t526 - pkin(5) * qJDD(1)) * t520 + t466;
t483 = -t520 * g(3) + t563 * t521;
t545 = t521 * qJD(1);
t468 = t495 * t545 + t483;
t542 = qJDD(1) * t521;
t462 = -pkin(3) * t551 - pkin(5) * t542 + t468;
t522 = sin(qJ(4));
t524 = cos(qJ(4));
t459 = t524 * t461 - t522 * t462;
t533 = -t520 * t522 - t521 * t524;
t492 = t533 * qJD(1);
t534 = t520 * t524 - t521 * t522;
t493 = t534 * qJD(1);
t474 = -t492 * mrSges(5,1) + t493 * mrSges(5,2);
t481 = t492 * qJD(4) + t534 * qJDD(1);
t484 = -qJD(4) * mrSges(5,2) + t492 * mrSges(5,3);
t456 = m(5) * t459 + qJDD(4) * mrSges(5,1) - t481 * mrSges(5,3) + qJD(4) * t484 - t493 * t474;
t460 = t522 * t461 + t524 * t462;
t480 = -t493 * qJD(4) + t533 * qJDD(1);
t485 = qJD(4) * mrSges(5,1) - t493 * mrSges(5,3);
t457 = m(5) * t460 - qJDD(4) * mrSges(5,2) + t480 * mrSges(5,3) - qJD(4) * t485 + t492 * t474;
t447 = t524 * t456 + t522 * t457;
t530 = m(4) * t466 + t447;
t444 = m(3) * t482 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t496 - t497) * qJD(1)) * t520 - t530;
t537 = -t522 * t456 + t524 * t457;
t532 = m(4) * t468 + mrSges(4,2) * t542 + t496 * t545 + t537;
t445 = m(3) * t483 + (qJDD(1) * mrSges(3,3) + qJD(1) * t497) * t521 + t532;
t538 = -t520 * t444 + t521 * t445;
t438 = m(2) * t502 - t526 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t538;
t501 = t523 * g(1) - t525 * g(2);
t540 = -qJDD(2) + t501;
t531 = -0.2e1 * qJD(3) * t546 - t540;
t479 = -t550 + (-pkin(1) + t536) * qJDD(1) + t531;
t465 = (qJ(2) + (-t516 - t517) * pkin(5)) * t526 + (t553 + pkin(1) + (pkin(2) + pkin(3)) * t521) * qJDD(1) - t531;
t535 = -m(5) * t465 + t480 * mrSges(5,1) - t481 * mrSges(5,2) + t492 * t484 - t493 * t485;
t543 = qJDD(1) * t520;
t454 = m(4) * t479 - mrSges(4,1) * t542 - t559 * mrSges(4,2) - mrSges(4,3) * t543 + t535;
t491 = -qJDD(1) * pkin(1) - t540 - t550;
t527 = -m(3) * t491 + mrSges(3,1) * t542 + t559 * mrSges(3,3) - t454;
t451 = (mrSges(2,1) - t554) * qJDD(1) - t526 * mrSges(2,2) + t527 + m(2) * t501;
t549 = t523 * t438 + t525 * t451;
t440 = t521 * t444 + t520 * t445;
t547 = t562 * qJD(1);
t539 = t525 * t438 - t523 * t451;
t469 = Ifges(5,5) * t493 + Ifges(5,6) * t492 + Ifges(5,3) * qJD(4);
t471 = Ifges(5,1) * t493 + Ifges(5,4) * t492 + Ifges(5,5) * qJD(4);
t448 = -mrSges(5,1) * t465 + mrSges(5,3) * t460 + Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * qJDD(4) + qJD(4) * t471 - t493 * t469;
t470 = Ifges(5,4) * t493 + Ifges(5,2) * t492 + Ifges(5,6) * qJD(4);
t449 = mrSges(5,2) * t465 - mrSges(5,3) * t459 + Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * qJDD(4) - qJD(4) * t470 + t492 * t469;
t433 = -mrSges(3,1) * t491 + mrSges(3,3) * t483 - mrSges(4,1) * t479 + mrSges(4,2) * t468 - t522 * t449 - t524 * t448 - pkin(3) * t535 - pkin(5) * t537 - pkin(2) * t454 - t547 * t546 + (t555 * t521 + t558) * qJDD(1);
t435 = mrSges(3,2) * t491 + mrSges(4,2) * t466 - mrSges(3,3) * t482 - mrSges(4,3) * t479 - pkin(5) * t447 - qJ(3) * t454 - t522 * t448 + t524 * t449 + t547 * t545 + (t557 * t520 + t556 * t521) * qJDD(1);
t453 = mrSges(3,2) * t543 - t527;
t529 = mrSges(2,1) * t501 - mrSges(2,2) * t502 + Ifges(2,3) * qJDD(1) - pkin(1) * t453 + qJ(2) * t538 + t521 * t433 + t520 * t435;
t528 = mrSges(5,1) * t459 - mrSges(5,2) * t460 + Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * qJDD(4) + t493 * t470 - t492 * t471;
t446 = (qJDD(1) * mrSges(4,2) + qJD(1) * t496) * t520 + t530;
t431 = pkin(3) * t447 + pkin(2) * t446 - pkin(1) * t440 - mrSges(3,1) * t482 + mrSges(3,2) * t483 + t528 + mrSges(4,1) * t466 - mrSges(4,3) * t468 + mrSges(2,3) * t502 + mrSges(2,1) * g(3) + (Ifges(2,6) - t562) * qJDD(1) - qJ(3) * t532 + (Ifges(2,5) + t556 * t517 + (-t558 + (-t555 + t557) * t521) * t520) * t526;
t430 = -mrSges(2,2) * g(3) - mrSges(2,3) * t501 + Ifges(2,5) * qJDD(1) - t526 * Ifges(2,6) - qJ(2) * t440 - t520 * t433 + t521 * t435;
t1 = [-m(1) * g(1) + t539; -m(1) * g(2) + t549; (-m(1) - m(2)) * g(3) + t440; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t549 + t525 * t430 - t523 * t431; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t539 + t523 * t430 + t525 * t431; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t529; t529; t453; t446; t528;];
tauJB = t1;
