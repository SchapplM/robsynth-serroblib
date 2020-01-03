% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:18
% EndTime: 2019-12-31 17:11:19
% DurationCPUTime: 1.19s
% Computational Cost: add. (7746->222), mult. (15797->267), div. (0->0), fcn. (7773->6), ass. (0->93)
t591 = Ifges(3,1) + Ifges(4,2);
t583 = Ifges(3,4) + Ifges(4,6);
t582 = Ifges(3,5) - Ifges(4,4);
t590 = Ifges(3,2) + Ifges(4,3);
t581 = Ifges(3,6) - Ifges(4,5);
t589 = Ifges(3,3) + Ifges(4,1);
t550 = sin(qJ(2));
t553 = cos(qJ(2));
t527 = (mrSges(4,2) * t553 - mrSges(4,3) * t550) * qJD(1);
t574 = qJD(1) * t553;
t535 = -mrSges(4,1) * t574 - qJD(2) * mrSges(4,3);
t572 = qJD(1) * qJD(2);
t571 = t550 * t572;
t530 = t553 * qJDD(1) - t571;
t573 = t550 * qJD(1);
t537 = pkin(3) * t573 - qJD(2) * pkin(6);
t548 = t553 ^ 2;
t556 = qJD(1) ^ 2;
t570 = t553 * t572;
t529 = t550 * qJDD(1) + t570;
t551 = sin(qJ(1));
t554 = cos(qJ(1));
t538 = t551 * g(1) - t554 * g(2);
t567 = -qJDD(1) * pkin(1) - t538;
t585 = -2 * qJD(3);
t559 = pkin(2) * t571 + t573 * t585 + (-t529 - t570) * qJ(3) + t567;
t483 = -t537 * t573 + (-pkin(3) * t548 - pkin(5)) * t556 + (-pkin(2) - pkin(6)) * t530 + t559;
t539 = -t554 * g(1) - t551 * g(2);
t516 = -t556 * pkin(1) + qJDD(1) * pkin(5) + t539;
t503 = -t553 * g(3) - t550 * t516;
t526 = (-pkin(2) * t553 - qJ(3) * t550) * qJD(1);
t555 = qJD(2) ^ 2;
t490 = -qJDD(2) * pkin(2) - t555 * qJ(3) + t526 * t573 + qJDD(3) - t503;
t486 = (-t550 * t553 * t556 - qJDD(2)) * pkin(6) + (t529 - t570) * pkin(3) + t490;
t549 = sin(qJ(4));
t552 = cos(qJ(4));
t481 = -t549 * t483 + t552 * t486;
t524 = -t549 * qJD(2) - t552 * t574;
t499 = t524 * qJD(4) + t552 * qJDD(2) - t549 * t530;
t525 = t552 * qJD(2) - t549 * t574;
t500 = -t524 * mrSges(5,1) + t525 * mrSges(5,2);
t541 = qJD(4) + t573;
t501 = -t541 * mrSges(5,2) + t524 * mrSges(5,3);
t523 = qJDD(4) + t529;
t478 = m(5) * t481 + t523 * mrSges(5,1) - t499 * mrSges(5,3) - t525 * t500 + t541 * t501;
t482 = t552 * t483 + t549 * t486;
t498 = -t525 * qJD(4) - t549 * qJDD(2) - t552 * t530;
t502 = t541 * mrSges(5,1) - t525 * mrSges(5,3);
t479 = m(5) * t482 - t523 * mrSges(5,2) + t498 * mrSges(5,3) + t524 * t500 - t541 * t502;
t469 = t552 * t478 + t549 * t479;
t564 = -m(4) * t490 - t529 * mrSges(4,1) - t469;
t468 = qJDD(2) * mrSges(4,2) + qJD(2) * t535 + t527 * t573 - t564;
t504 = -t550 * g(3) + t553 * t516;
t561 = -t555 * pkin(2) + qJDD(2) * qJ(3) + t526 * t574 + t504;
t485 = -t548 * t556 * pkin(6) + t530 * pkin(3) + ((2 * qJD(3)) + t537) * qJD(2) + t561;
t491 = Ifges(5,5) * t525 + Ifges(5,6) * t524 + Ifges(5,3) * t541;
t493 = Ifges(5,1) * t525 + Ifges(5,4) * t524 + Ifges(5,5) * t541;
t470 = -mrSges(5,1) * t485 + mrSges(5,3) * t482 + Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t523 - t525 * t491 + t541 * t493;
t492 = Ifges(5,4) * t525 + Ifges(5,2) * t524 + Ifges(5,6) * t541;
t471 = mrSges(5,2) * t485 - mrSges(5,3) * t481 + Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t523 + t524 * t491 - t541 * t492;
t489 = qJD(2) * t585 - t561;
t536 = mrSges(4,1) * t573 + qJD(2) * mrSges(4,2);
t565 = -m(5) * t485 + t498 * mrSges(5,1) - t499 * mrSges(5,2) + t524 * t501 - t525 * t502;
t560 = -m(4) * t489 + qJDD(2) * mrSges(4,3) + qJD(2) * t536 + t527 * t574 - t565;
t575 = t582 * qJD(2) + (t591 * t550 + t583 * t553) * qJD(1);
t576 = t581 * qJD(2) + (t583 * t550 + t590 * t553) * qJD(1);
t588 = (t550 * t576 - t575 * t553) * qJD(1) + t589 * qJDD(2) + t582 * t529 + t581 * t530 + mrSges(3,1) * t503 - mrSges(3,2) * t504 + mrSges(4,2) * t490 - mrSges(4,3) * t489 - pkin(2) * t468 - pkin(6) * t469 + qJ(3) * (t530 * mrSges(4,1) + t560) - t549 * t470 + t552 * t471;
t584 = t556 * pkin(5);
t528 = (-mrSges(3,1) * t553 + mrSges(3,2) * t550) * qJD(1);
t534 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t574;
t466 = m(3) * t503 - t529 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t534 - t535) * qJD(2) + (-t527 - t528) * t573 + t564;
t533 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t573;
t474 = t528 * t574 + m(3) * t504 - qJDD(2) * mrSges(3,2) - qJD(2) * t533 + (mrSges(3,3) + mrSges(4,1)) * t530 + t560;
t568 = -t550 * t466 + t553 * t474;
t459 = m(2) * t539 - t556 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t568;
t515 = t567 - t584;
t487 = -t530 * pkin(2) + t559 - t584;
t578 = -t549 * t478 + t552 * t479;
t566 = -m(4) * t487 - t530 * mrSges(4,2) + t536 * t573 - t578;
t558 = -m(3) * t515 + t534 * t574 + t530 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t529 + (-t533 * t550 - t535 * t553) * qJD(1) + t566;
t463 = m(2) * t538 + qJDD(1) * mrSges(2,1) - t556 * mrSges(2,2) + t558;
t579 = t551 * t459 + t554 * t463;
t461 = t553 * t466 + t550 * t474;
t577 = t589 * qJD(2) + (t582 * t550 + t581 * t553) * qJD(1);
t569 = t554 * t459 - t551 * t463;
t467 = -t529 * mrSges(4,3) + t535 * t574 - t566;
t454 = -mrSges(3,1) * t515 - mrSges(4,1) * t489 + mrSges(4,2) * t487 + mrSges(3,3) * t504 - pkin(2) * t467 - pkin(3) * t565 - pkin(6) * t578 + t575 * qJD(2) + t581 * qJDD(2) - t552 * t470 - t549 * t471 + t583 * t529 + t590 * t530 - t577 * t573;
t562 = mrSges(5,1) * t481 - mrSges(5,2) * t482 + Ifges(5,5) * t499 + Ifges(5,6) * t498 + Ifges(5,3) * t523 + t525 * t492 - t524 * t493;
t456 = mrSges(4,1) * t490 + mrSges(3,2) * t515 - mrSges(3,3) * t503 - mrSges(4,3) * t487 + pkin(3) * t469 - qJ(3) * t467 - t576 * qJD(2) + t582 * qJDD(2) + t591 * t529 + t583 * t530 + t577 * t574 + t562;
t563 = mrSges(2,1) * t538 - mrSges(2,2) * t539 + Ifges(2,3) * qJDD(1) + pkin(1) * t558 + pkin(5) * t568 + t553 * t454 + t550 * t456;
t452 = mrSges(2,1) * g(3) + mrSges(2,3) * t539 + t556 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t461 - t588;
t451 = -mrSges(2,2) * g(3) - mrSges(2,3) * t538 + Ifges(2,5) * qJDD(1) - t556 * Ifges(2,6) - pkin(5) * t461 - t550 * t454 + t553 * t456;
t1 = [-m(1) * g(1) + t569; -m(1) * g(2) + t579; (-m(1) - m(2)) * g(3) + t461; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t579 + t554 * t451 - t551 * t452; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t569 + t551 * t451 + t554 * t452; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t563; t563; t588; t468; t562;];
tauJB = t1;
