% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:38
% EndTime: 2019-12-31 18:05:40
% DurationCPUTime: 0.96s
% Computational Cost: add. (7999->216), mult. (14506->251), div. (0->0), fcn. (6761->6), ass. (0->87)
t555 = sin(qJ(1));
t558 = cos(qJ(1));
t532 = t555 * g(1) - t558 * g(2);
t560 = qJD(1) ^ 2;
t512 = -qJDD(1) * pkin(1) - t560 * qJ(2) + qJDD(2) - t532;
t504 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t512;
t533 = -t558 * g(1) - t555 * g(2);
t588 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t533;
t587 = -m(3) - m(4);
t554 = sin(qJ(4));
t586 = t554 * g(3);
t585 = mrSges(2,1) - mrSges(3,2);
t584 = t560 * mrSges(4,3);
t510 = t560 * pkin(1) - t588;
t505 = qJDD(3) + (-pkin(1) - qJ(3)) * t560 + t588;
t501 = -qJDD(1) * pkin(6) + t505;
t557 = cos(qJ(4));
t496 = -t557 * g(3) + t554 * t501;
t525 = (mrSges(5,1) * t554 + mrSges(5,2) * t557) * qJD(1);
t580 = qJD(1) * qJD(4);
t574 = t557 * t580;
t527 = -t554 * qJDD(1) - t574;
t582 = qJD(1) * t557;
t531 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t582;
t500 = -t560 * pkin(6) - t504;
t575 = t554 * t580;
t528 = t557 * qJDD(1) - t575;
t484 = (-t528 + t575) * pkin(7) + (-t527 + t574) * pkin(4) + t500;
t526 = (pkin(4) * t554 - pkin(7) * t557) * qJD(1);
t559 = qJD(4) ^ 2;
t581 = t554 * qJD(1);
t486 = -t559 * pkin(4) + qJDD(4) * pkin(7) - t526 * t581 + t496;
t553 = sin(qJ(5));
t556 = cos(qJ(5));
t482 = t556 * t484 - t553 * t486;
t523 = t556 * qJD(4) - t553 * t582;
t494 = t523 * qJD(5) + t553 * qJDD(4) + t556 * t528;
t524 = t553 * qJD(4) + t556 * t582;
t498 = -t523 * mrSges(6,1) + t524 * mrSges(6,2);
t534 = qJD(5) + t581;
t506 = -t534 * mrSges(6,2) + t523 * mrSges(6,3);
t522 = qJDD(5) - t527;
t479 = m(6) * t482 + t522 * mrSges(6,1) - t494 * mrSges(6,3) - t524 * t498 + t534 * t506;
t483 = t553 * t484 + t556 * t486;
t493 = -t524 * qJD(5) + t556 * qJDD(4) - t553 * t528;
t507 = t534 * mrSges(6,1) - t524 * mrSges(6,3);
t480 = m(6) * t483 - t522 * mrSges(6,2) + t493 * mrSges(6,3) + t523 * t498 - t534 * t507;
t571 = -t553 * t479 + t556 * t480;
t467 = m(5) * t496 - qJDD(4) * mrSges(5,2) + t527 * mrSges(5,3) - qJD(4) * t531 - t525 * t581 + t571;
t495 = t557 * t501 + t586;
t530 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t581;
t485 = -qJDD(4) * pkin(4) - t559 * pkin(7) - t586 + (qJD(1) * t526 - t501) * t557;
t565 = -m(6) * t485 + t493 * mrSges(6,1) - t494 * mrSges(6,2) + t523 * t506 - t524 * t507;
t475 = m(5) * t495 + qJDD(4) * mrSges(5,1) - t528 * mrSges(5,3) + qJD(4) * t530 - t525 * t582 + t565;
t459 = t554 * t467 + t557 * t475;
t570 = m(4) * t505 + qJDD(1) * mrSges(4,2) + t459;
t566 = -m(3) * t510 + t560 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t570;
t455 = m(2) * t533 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t560 + t566;
t469 = t556 * t479 + t553 * t480;
t567 = -m(5) * t500 + t527 * mrSges(5,1) - t528 * mrSges(5,2) - t530 * t581 - t531 * t582 - t469;
t464 = m(4) * t504 - t560 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t567;
t563 = -m(3) * t512 + t560 * mrSges(3,3) - t464;
t461 = m(2) * t532 - t560 * mrSges(2,2) + t585 * qJDD(1) + t563;
t583 = t555 * t455 + t558 * t461;
t577 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t576 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t573 = t558 * t455 - t555 * t461;
t572 = t557 * t467 - t554 * t475;
t487 = Ifges(6,5) * t524 + Ifges(6,6) * t523 + Ifges(6,3) * t534;
t489 = Ifges(6,1) * t524 + Ifges(6,4) * t523 + Ifges(6,5) * t534;
t472 = -mrSges(6,1) * t485 + mrSges(6,3) * t483 + Ifges(6,4) * t494 + Ifges(6,2) * t493 + Ifges(6,6) * t522 - t524 * t487 + t534 * t489;
t488 = Ifges(6,4) * t524 + Ifges(6,2) * t523 + Ifges(6,6) * t534;
t473 = mrSges(6,2) * t485 - mrSges(6,3) * t482 + Ifges(6,1) * t494 + Ifges(6,4) * t493 + Ifges(6,5) * t522 + t523 * t487 - t534 * t488;
t514 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t557 - Ifges(5,2) * t554) * qJD(1);
t515 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t557 - Ifges(5,4) * t554) * qJD(1);
t564 = mrSges(5,1) * t495 - mrSges(5,2) * t496 + Ifges(5,5) * t528 + Ifges(5,6) * t527 + Ifges(5,3) * qJDD(4) + pkin(4) * t565 + pkin(7) * t571 + t556 * t472 + t553 * t473 + t514 * t582 + t515 * t581;
t562 = mrSges(6,1) * t482 - mrSges(6,2) * t483 + Ifges(6,5) * t494 + Ifges(6,6) * t493 + Ifges(6,3) * t522 + t524 * t488 - t523 * t489;
t513 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t557 - Ifges(5,6) * t554) * qJD(1);
t451 = mrSges(5,2) * t500 - mrSges(5,3) * t495 + Ifges(5,1) * t528 + Ifges(5,4) * t527 + Ifges(5,5) * qJDD(4) - pkin(7) * t469 - qJD(4) * t514 - t553 * t472 + t556 * t473 - t513 * t581;
t452 = -mrSges(5,1) * t500 + mrSges(5,3) * t496 + Ifges(5,4) * t528 + Ifges(5,2) * t527 + Ifges(5,6) * qJDD(4) - pkin(4) * t469 + qJD(4) * t515 - t513 * t582 - t562;
t463 = qJDD(1) * mrSges(3,2) - t563;
t561 = -mrSges(2,2) * t533 - mrSges(3,3) * t510 - mrSges(4,3) * t504 - pkin(6) * t459 - qJ(3) * t464 + t557 * t451 - t554 * t452 + qJ(2) * (t566 - t584) - pkin(1) * t463 + mrSges(4,2) * t505 + mrSges(3,2) * t512 + mrSges(2,1) * t532 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t458 = t587 * g(3) + t572;
t457 = t570 - t584;
t449 = t564 + t577 * t560 - qJ(3) * t572 + t576 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t585) * g(3) + mrSges(2,3) * t533 + mrSges(4,1) * t505 - mrSges(3,1) * t510 + pkin(3) * t459 + pkin(2) * t457 - pkin(1) * t458;
t448 = -qJ(2) * t458 - mrSges(2,3) * t532 + pkin(2) * t464 + mrSges(3,1) * t512 + t554 * t451 + t557 * t452 + pkin(3) * t567 + pkin(6) * t572 + mrSges(4,1) * t504 - t576 * t560 + t577 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t573; -m(1) * g(2) + t583; (-m(1) - m(2) + t587) * g(3) + t572; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t583 + t558 * t448 - t555 * t449; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t573 + t555 * t448 + t558 * t449; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t561; t561; t463; t457; t564; t562;];
tauJB = t1;
