% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:02
% EndTime: 2020-01-03 12:02:05
% DurationCPUTime: 3.08s
% Computational Cost: add. (53053->230), mult. (68329->290), div. (0->0), fcn. (38611->10), ass. (0->99)
t585 = sin(qJ(1));
t589 = cos(qJ(1));
t564 = -t589 * g(2) - t585 * g(3);
t557 = qJDD(1) * pkin(1) + t564;
t563 = -t585 * g(2) + t589 * g(3);
t590 = qJD(1) ^ 2;
t558 = -t590 * pkin(1) + t563;
t584 = sin(qJ(2));
t588 = cos(qJ(2));
t540 = t588 * t557 - t584 * t558;
t574 = qJDD(1) + qJDD(2);
t537 = t574 * pkin(2) + t540;
t541 = t584 * t557 + t588 * t558;
t576 = qJD(1) + qJD(2);
t572 = t576 ^ 2;
t538 = -t572 * pkin(2) + t541;
t580 = sin(pkin(9));
t581 = cos(pkin(9));
t522 = t580 * t537 + t581 * t538;
t519 = -t572 * pkin(3) + t574 * pkin(7) + t522;
t579 = -g(1) + qJDD(3);
t583 = sin(qJ(4));
t587 = cos(qJ(4));
t515 = -t583 * t519 + t587 * t579;
t605 = qJD(4) * t576;
t604 = t587 * t605;
t552 = t583 * t574 + t604;
t512 = (-t552 + t604) * pkin(8) + (t572 * t583 * t587 + qJDD(4)) * pkin(4) + t515;
t516 = t587 * t519 + t583 * t579;
t553 = t587 * t574 - t583 * t605;
t608 = t576 * t583;
t561 = qJD(4) * pkin(4) - pkin(8) * t608;
t578 = t587 ^ 2;
t513 = -t578 * t572 * pkin(4) + t553 * pkin(8) - qJD(4) * t561 + t516;
t582 = sin(qJ(5));
t586 = cos(qJ(5));
t510 = t586 * t512 - t582 * t513;
t547 = (-t582 * t583 + t586 * t587) * t576;
t528 = t547 * qJD(5) + t586 * t552 + t582 * t553;
t548 = (t582 * t587 + t583 * t586) * t576;
t533 = -t547 * mrSges(6,1) + t548 * mrSges(6,2);
t575 = qJD(4) + qJD(5);
t542 = -t575 * mrSges(6,2) + t547 * mrSges(6,3);
t573 = qJDD(4) + qJDD(5);
t507 = m(6) * t510 + t573 * mrSges(6,1) - t528 * mrSges(6,3) - t548 * t533 + t575 * t542;
t511 = t582 * t512 + t586 * t513;
t527 = -t548 * qJD(5) - t582 * t552 + t586 * t553;
t543 = t575 * mrSges(6,1) - t548 * mrSges(6,3);
t508 = m(6) * t511 - t573 * mrSges(6,2) + t527 * mrSges(6,3) + t547 * t533 - t575 * t543;
t498 = t586 * t507 + t582 * t508;
t545 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t583 + Ifges(5,2) * t587) * t576;
t546 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t583 + Ifges(5,4) * t587) * t576;
t530 = Ifges(6,4) * t548 + Ifges(6,2) * t547 + Ifges(6,6) * t575;
t531 = Ifges(6,1) * t548 + Ifges(6,4) * t547 + Ifges(6,5) * t575;
t595 = -mrSges(6,1) * t510 + mrSges(6,2) * t511 - Ifges(6,5) * t528 - Ifges(6,6) * t527 - Ifges(6,3) * t573 - t548 * t530 + t547 * t531;
t609 = mrSges(5,1) * t515 - mrSges(5,2) * t516 + Ifges(5,5) * t552 + Ifges(5,6) * t553 + Ifges(5,3) * qJDD(4) + pkin(4) * t498 + (t583 * t545 - t587 * t546) * t576 - t595;
t607 = t576 * t587;
t551 = (-mrSges(5,1) * t587 + mrSges(5,2) * t583) * t576;
t560 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t607;
t496 = m(5) * t515 + qJDD(4) * mrSges(5,1) - t552 * mrSges(5,3) + qJD(4) * t560 - t551 * t608 + t498;
t559 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t608;
t599 = -t582 * t507 + t586 * t508;
t497 = m(5) * t516 - qJDD(4) * mrSges(5,2) + t553 * mrSges(5,3) - qJD(4) * t559 + t551 * t607 + t599;
t600 = -t583 * t496 + t587 * t497;
t489 = m(4) * t522 - t572 * mrSges(4,1) - t574 * mrSges(4,2) + t600;
t521 = t581 * t537 - t580 * t538;
t597 = -t574 * pkin(3) - t521;
t518 = -t572 * pkin(7) + t597;
t514 = t561 * t608 - t553 * pkin(4) + (-pkin(8) * t578 - pkin(7)) * t572 + t597;
t596 = m(6) * t514 - t527 * mrSges(6,1) + t528 * mrSges(6,2) - t547 * t542 + t548 * t543;
t592 = -m(5) * t518 + t553 * mrSges(5,1) - t552 * mrSges(5,2) - t559 * t608 + t560 * t607 - t596;
t502 = m(4) * t521 + t574 * mrSges(4,1) - t572 * mrSges(4,2) + t592;
t484 = t580 * t489 + t581 * t502;
t481 = m(3) * t540 + t574 * mrSges(3,1) - t572 * mrSges(3,2) + t484;
t601 = t581 * t489 - t580 * t502;
t482 = m(3) * t541 - t572 * mrSges(3,1) - t574 * mrSges(3,2) + t601;
t602 = -t584 * t481 + t588 * t482;
t471 = m(2) * t563 - t590 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t602;
t474 = t588 * t481 + t584 * t482;
t472 = m(2) * t564 + qJDD(1) * mrSges(2,1) - t590 * mrSges(2,2) + t474;
t606 = t585 * t471 + t589 * t472;
t492 = t587 * t496 + t583 * t497;
t490 = m(4) * t579 + t492;
t603 = -t589 * t471 + t585 * t472;
t529 = Ifges(6,5) * t548 + Ifges(6,6) * t547 + Ifges(6,3) * t575;
t499 = -mrSges(6,1) * t514 + mrSges(6,3) * t511 + Ifges(6,4) * t528 + Ifges(6,2) * t527 + Ifges(6,6) * t573 - t548 * t529 + t575 * t531;
t500 = mrSges(6,2) * t514 - mrSges(6,3) * t510 + Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * t573 + t547 * t529 - t575 * t530;
t544 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t583 + Ifges(5,6) * t587) * t576;
t477 = -mrSges(5,1) * t518 + mrSges(5,3) * t516 + Ifges(5,4) * t552 + Ifges(5,2) * t553 + Ifges(5,6) * qJDD(4) - pkin(4) * t596 + pkin(8) * t599 + qJD(4) * t546 + t586 * t499 + t582 * t500 - t544 * t608;
t486 = mrSges(5,2) * t518 - mrSges(5,3) * t515 + Ifges(5,1) * t552 + Ifges(5,4) * t553 + Ifges(5,5) * qJDD(4) - pkin(8) * t498 - qJD(4) * t545 - t582 * t499 + t586 * t500 + t544 * t607;
t594 = mrSges(3,1) * t540 + mrSges(4,1) * t521 - mrSges(3,2) * t541 - mrSges(4,2) * t522 + pkin(2) * t484 + pkin(3) * t592 + pkin(7) * t600 + t587 * t477 + t583 * t486 + (Ifges(4,3) + Ifges(3,3)) * t574;
t593 = mrSges(2,1) * t564 - mrSges(2,2) * t563 + Ifges(2,3) * qJDD(1) + pkin(1) * t474 + t594;
t475 = -mrSges(4,1) * t579 + mrSges(4,3) * t522 + t572 * Ifges(4,5) + Ifges(4,6) * t574 - pkin(3) * t492 - t609;
t467 = mrSges(4,2) * t579 - mrSges(4,3) * t521 + Ifges(4,5) * t574 - t572 * Ifges(4,6) - pkin(7) * t492 - t583 * t477 + t587 * t486;
t466 = -mrSges(3,2) * g(1) - mrSges(3,3) * t540 + Ifges(3,5) * t574 - t572 * Ifges(3,6) - qJ(3) * t484 + t581 * t467 - t580 * t475;
t465 = mrSges(3,1) * g(1) + mrSges(3,3) * t541 + t572 * Ifges(3,5) + Ifges(3,6) * t574 - pkin(2) * t490 + qJ(3) * t601 + t580 * t467 + t581 * t475;
t464 = -mrSges(2,2) * g(1) - mrSges(2,3) * t564 + Ifges(2,5) * qJDD(1) - t590 * Ifges(2,6) - pkin(6) * t474 - t584 * t465 + t588 * t466;
t463 = Ifges(2,6) * qJDD(1) + t590 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t563 + t584 * t466 + t588 * t465 - pkin(1) * (-m(3) * g(1) + t490) + pkin(6) * t602;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t490; -m(1) * g(2) + t606; -m(1) * g(3) + t603; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t593; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t603 + t589 * t463 + t585 * t464; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t606 + t585 * t463 - t589 * t464; t593; t594; t490; t609; -t595;];
tauJB = t1;
