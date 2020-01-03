% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:44
% EndTime: 2019-12-31 19:52:45
% DurationCPUTime: 1.31s
% Computational Cost: add. (10971->200), mult. (13016->234), div. (0->0), fcn. (5490->6), ass. (0->87)
t590 = Ifges(5,1) + Ifges(6,1);
t577 = Ifges(5,4) - Ifges(6,5);
t588 = Ifges(5,5) + Ifges(6,4);
t589 = Ifges(5,2) + Ifges(6,3);
t587 = Ifges(5,6) - Ifges(6,6);
t586 = Ifges(5,3) + Ifges(6,2);
t541 = qJD(1) + qJD(2);
t547 = sin(qJ(4));
t550 = cos(qJ(4));
t585 = (t589 * t547 - t577 * t550) * t541 - t587 * qJD(4);
t584 = (-t577 * t547 + t590 * t550) * t541 + t588 * qJD(4);
t549 = sin(qJ(1));
t552 = cos(qJ(1));
t529 = t549 * g(1) - t552 * g(2);
t522 = qJDD(1) * pkin(1) + t529;
t530 = -t552 * g(1) - t549 * g(2);
t554 = qJD(1) ^ 2;
t523 = -t554 * pkin(1) + t530;
t548 = sin(qJ(2));
t551 = cos(qJ(2));
t493 = t548 * t522 + t551 * t523;
t540 = qJDD(1) + qJDD(2);
t583 = -t540 * qJ(3) - 0.2e1 * qJD(3) * t541 - t493;
t582 = -m(3) - m(4);
t581 = -pkin(2) - pkin(7);
t580 = t547 * g(3);
t579 = mrSges(3,1) - mrSges(4,2);
t578 = -mrSges(5,3) - mrSges(6,2);
t576 = Ifges(3,5) - Ifges(4,4);
t575 = -Ifges(3,6) + Ifges(4,5);
t574 = t541 * t547;
t573 = t541 * t550;
t492 = t551 * t522 - t548 * t523;
t539 = t541 ^ 2;
t561 = -t539 * qJ(3) + qJDD(3) - t492;
t487 = t581 * t540 + t561;
t483 = -t550 * g(3) + t547 * t487;
t570 = qJD(4) * t541;
t515 = t547 * t540 + t550 * t570;
t525 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t573;
t513 = (mrSges(6,1) * t547 - mrSges(6,3) * t550) * t541;
t566 = t541 * (-t513 - (mrSges(5,1) * t547 + mrSges(5,2) * t550) * t541);
t512 = (pkin(4) * t547 - qJ(5) * t550) * t541;
t553 = qJD(4) ^ 2;
t479 = -t553 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t512 * t574 + t483;
t526 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t573;
t569 = m(6) * t479 + qJDD(4) * mrSges(6,3) + qJD(4) * t526;
t469 = m(5) * t483 - qJDD(4) * mrSges(5,2) - qJD(4) * t525 + t578 * t515 + t547 * t566 + t569;
t482 = t550 * t487 + t580;
t516 = t550 * t540 - t547 * t570;
t524 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t574;
t480 = -qJDD(4) * pkin(4) - t580 - t553 * qJ(5) + qJDD(5) + (t512 * t541 - t487) * t550;
t527 = -mrSges(6,2) * t574 + qJD(4) * mrSges(6,3);
t562 = -m(6) * t480 + qJDD(4) * mrSges(6,1) + qJD(4) * t527;
t470 = m(5) * t482 + qJDD(4) * mrSges(5,1) + qJD(4) * t524 + t578 * t516 + t550 * t566 + t562;
t463 = t547 * t469 + t550 * t470;
t490 = -t540 * pkin(2) + t561;
t560 = -m(4) * t490 + t539 * mrSges(4,3) - t463;
t459 = m(3) * t492 - t539 * mrSges(3,2) + t579 * t540 + t560;
t488 = t539 * pkin(2) + t583;
t485 = t581 * t539 - t583;
t476 = t515 * pkin(4) - t516 * qJ(5) + (-0.2e1 * qJD(5) * t550 + (pkin(4) * t550 + qJ(5) * t547) * qJD(4)) * t541 + t485;
t471 = m(6) * t476 + t515 * mrSges(6,1) - t516 * mrSges(6,3) - t526 * t573 + t527 * t574;
t559 = -m(5) * t485 - t515 * mrSges(5,1) - t516 * mrSges(5,2) - t524 * t574 - t525 * t573 - t471;
t557 = -m(4) * t488 + t539 * mrSges(4,2) + t540 * mrSges(4,3) - t559;
t466 = m(3) * t493 - t539 * mrSges(3,1) - t540 * mrSges(3,2) + t557;
t454 = t551 * t459 + t548 * t466;
t451 = m(2) * t529 + qJDD(1) * mrSges(2,1) - t554 * mrSges(2,2) + t454;
t564 = -t548 * t459 + t551 * t466;
t452 = m(2) * t530 - t554 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t564;
t572 = t552 * t451 + t549 * t452;
t567 = t541 * ((t587 * t547 - t588 * t550) * t541 - t586 * qJD(4));
t565 = -t549 * t451 + t552 * t452;
t563 = -t550 * t469 + t547 * t470;
t456 = -mrSges(5,1) * t485 - mrSges(6,1) * t476 + mrSges(6,2) * t479 + mrSges(5,3) * t483 - pkin(4) * t471 + t584 * qJD(4) + t587 * qJDD(4) - t589 * t515 + t577 * t516 + t550 * t567;
t457 = mrSges(5,2) * t485 + mrSges(6,2) * t480 - mrSges(5,3) * t482 - mrSges(6,3) * t476 - qJ(5) * t471 + t585 * qJD(4) + t588 * qJDD(4) - t577 * t515 + t590 * t516 + t547 * t567;
t461 = t540 * mrSges(4,2) - t560;
t558 = mrSges(3,1) * t492 - mrSges(3,2) * t493 + mrSges(4,2) * t490 - mrSges(4,3) * t488 - pkin(2) * t461 - pkin(7) * t463 + qJ(3) * t557 - t547 * t456 + t550 * t457 + (Ifges(3,3) + Ifges(4,1)) * t540;
t474 = t516 * mrSges(6,2) + t513 * t573 - t562;
t556 = mrSges(5,1) * t482 - mrSges(6,1) * t480 - mrSges(5,2) * t483 + mrSges(6,3) * t479 - pkin(4) * t474 + qJ(5) * t569 + (-qJ(5) * t513 + t584) * t574 - t585 * t573 + t588 * t516 + (-qJ(5) * mrSges(6,2) - t587) * t515 + t586 * qJDD(4);
t555 = mrSges(2,1) * t529 - mrSges(2,2) * t530 + Ifges(2,3) * qJDD(1) + pkin(1) * t454 + t558;
t462 = -m(4) * g(3) - t563;
t447 = (-mrSges(3,2) + mrSges(4,3)) * g(3) + t556 + t576 * t540 + t575 * t539 - qJ(3) * t462 + pkin(3) * t463 + mrSges(4,1) * t490 - mrSges(3,3) * t492;
t446 = -mrSges(4,1) * t488 + mrSges(3,3) * t493 - pkin(2) * t462 - pkin(3) * t559 + pkin(7) * t563 + t579 * g(3) - t550 * t456 - t547 * t457 + t576 * t539 - t575 * t540;
t445 = -mrSges(2,2) * g(3) - mrSges(2,3) * t529 + Ifges(2,5) * qJDD(1) - t554 * Ifges(2,6) - pkin(6) * t454 - t548 * t446 + t551 * t447;
t444 = Ifges(2,6) * qJDD(1) + t554 * Ifges(2,5) + mrSges(2,3) * t530 + t548 * t447 + t551 * t446 + pkin(1) * t563 + pkin(6) * t564 + (-pkin(1) * t582 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t565; -m(1) * g(2) + t572; (-m(1) - m(2) + t582) * g(3) - t563; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t572 - t549 * t444 + t552 * t445; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t565 + t552 * t444 + t549 * t445; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t555; t555; t558; t461; t556; t474;];
tauJB = t1;
