% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:52
% EndTime: 2019-12-31 17:25:55
% DurationCPUTime: 2.34s
% Computational Cost: add. (24265->241), mult. (49557->309), div. (0->0), fcn. (31857->8), ass. (0->99)
t554 = sin(qJ(3));
t555 = sin(qJ(2));
t558 = cos(qJ(3));
t559 = cos(qJ(2));
t530 = (t554 * t555 - t558 * t559) * qJD(1);
t576 = qJD(1) * qJD(2);
t538 = t555 * qJDD(1) + t559 * t576;
t556 = sin(qJ(1));
t560 = cos(qJ(1));
t545 = -t560 * g(1) - t556 * g(2);
t561 = qJD(1) ^ 2;
t533 = -t561 * pkin(1) + qJDD(1) * pkin(5) + t545;
t580 = t555 * t533;
t581 = pkin(2) * t561;
t502 = qJDD(2) * pkin(2) - t538 * pkin(6) - t580 + (pkin(6) * t576 + t555 * t581 - g(3)) * t559;
t522 = -t555 * g(3) + t559 * t533;
t539 = t559 * qJDD(1) - t555 * t576;
t578 = qJD(1) * t555;
t543 = qJD(2) * pkin(2) - pkin(6) * t578;
t552 = t559 ^ 2;
t503 = t539 * pkin(6) - qJD(2) * t543 - t552 * t581 + t522;
t490 = t554 * t502 + t558 * t503;
t531 = (t554 * t559 + t555 * t558) * qJD(1);
t508 = -t531 * qJD(3) - t554 * t538 + t558 * t539;
t517 = t530 * mrSges(4,1) + t531 * mrSges(4,2);
t550 = qJD(2) + qJD(3);
t524 = t550 * mrSges(4,1) - t531 * mrSges(4,3);
t549 = qJDD(2) + qJDD(3);
t509 = -t530 * qJD(3) + t558 * t538 + t554 * t539;
t544 = t556 * g(1) - t560 * g(2);
t569 = -qJDD(1) * pkin(1) - t544;
t510 = -t539 * pkin(2) + t543 * t578 + (-pkin(6) * t552 - pkin(5)) * t561 + t569;
t485 = (t530 * t550 - t509) * pkin(7) + (t531 * t550 - t508) * pkin(3) + t510;
t518 = t530 * pkin(3) - t531 * pkin(7);
t548 = t550 ^ 2;
t487 = -t548 * pkin(3) + t549 * pkin(7) - t530 * t518 + t490;
t553 = sin(qJ(4));
t557 = cos(qJ(4));
t483 = t557 * t485 - t553 * t487;
t519 = -t553 * t531 + t557 * t550;
t493 = t519 * qJD(4) + t557 * t509 + t553 * t549;
t520 = t557 * t531 + t553 * t550;
t500 = -t519 * mrSges(5,1) + t520 * mrSges(5,2);
t507 = qJDD(4) - t508;
t526 = qJD(4) + t530;
t511 = -t526 * mrSges(5,2) + t519 * mrSges(5,3);
t479 = m(5) * t483 + t507 * mrSges(5,1) - t493 * mrSges(5,3) - t520 * t500 + t526 * t511;
t484 = t553 * t485 + t557 * t487;
t492 = -t520 * qJD(4) - t553 * t509 + t557 * t549;
t512 = t526 * mrSges(5,1) - t520 * mrSges(5,3);
t480 = m(5) * t484 - t507 * mrSges(5,2) + t492 * mrSges(5,3) + t519 * t500 - t526 * t512;
t572 = -t553 * t479 + t557 * t480;
t466 = m(4) * t490 - t549 * mrSges(4,2) + t508 * mrSges(4,3) - t530 * t517 - t550 * t524 + t572;
t489 = t558 * t502 - t554 * t503;
t523 = -t550 * mrSges(4,2) - t530 * mrSges(4,3);
t486 = -t549 * pkin(3) - t548 * pkin(7) + t531 * t518 - t489;
t567 = -m(5) * t486 + t492 * mrSges(5,1) - t493 * mrSges(5,2) + t519 * t511 - t520 * t512;
t475 = m(4) * t489 + t549 * mrSges(4,1) - t509 * mrSges(4,3) - t531 * t517 + t550 * t523 + t567;
t460 = t554 * t466 + t558 * t475;
t521 = -t559 * g(3) - t580;
t528 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t555 + Ifges(3,2) * t559) * qJD(1);
t529 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t555 + Ifges(3,4) * t559) * qJD(1);
t494 = Ifges(5,5) * t520 + Ifges(5,6) * t519 + Ifges(5,3) * t526;
t496 = Ifges(5,1) * t520 + Ifges(5,4) * t519 + Ifges(5,5) * t526;
t472 = -mrSges(5,1) * t486 + mrSges(5,3) * t484 + Ifges(5,4) * t493 + Ifges(5,2) * t492 + Ifges(5,6) * t507 - t520 * t494 + t526 * t496;
t495 = Ifges(5,4) * t520 + Ifges(5,2) * t519 + Ifges(5,6) * t526;
t473 = mrSges(5,2) * t486 - mrSges(5,3) * t483 + Ifges(5,1) * t493 + Ifges(5,4) * t492 + Ifges(5,5) * t507 + t519 * t494 - t526 * t495;
t514 = Ifges(4,4) * t531 - Ifges(4,2) * t530 + Ifges(4,6) * t550;
t515 = Ifges(4,1) * t531 - Ifges(4,4) * t530 + Ifges(4,5) * t550;
t565 = -mrSges(4,1) * t489 + mrSges(4,2) * t490 - Ifges(4,5) * t509 - Ifges(4,6) * t508 - Ifges(4,3) * t549 - pkin(3) * t567 - pkin(7) * t572 - t557 * t472 - t553 * t473 - t531 * t514 - t530 * t515;
t582 = mrSges(3,1) * t521 - mrSges(3,2) * t522 + Ifges(3,5) * t538 + Ifges(3,6) * t539 + Ifges(3,3) * qJDD(2) + pkin(2) * t460 + (t555 * t528 - t559 * t529) * qJD(1) - t565;
t537 = (-mrSges(3,1) * t559 + mrSges(3,2) * t555) * qJD(1);
t577 = qJD(1) * t559;
t542 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t577;
t458 = m(3) * t521 + qJDD(2) * mrSges(3,1) - t538 * mrSges(3,3) + qJD(2) * t542 - t537 * t578 + t460;
t541 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t578;
t573 = t558 * t466 - t554 * t475;
t459 = m(3) * t522 - qJDD(2) * mrSges(3,2) + t539 * mrSges(3,3) - qJD(2) * t541 + t537 * t577 + t573;
t574 = -t555 * t458 + t559 * t459;
t450 = m(2) * t545 - t561 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t574;
t532 = -t561 * pkin(5) + t569;
t468 = t557 * t479 + t553 * t480;
t566 = m(4) * t510 - t508 * mrSges(4,1) + t509 * mrSges(4,2) + t530 * t523 + t531 * t524 + t468;
t563 = -m(3) * t532 + t539 * mrSges(3,1) - t538 * mrSges(3,2) - t541 * t578 + t542 * t577 - t566;
t462 = m(2) * t544 + qJDD(1) * mrSges(2,1) - t561 * mrSges(2,2) + t563;
t579 = t556 * t450 + t560 * t462;
t452 = t559 * t458 + t555 * t459;
t575 = t560 * t450 - t556 * t462;
t513 = Ifges(4,5) * t531 - Ifges(4,6) * t530 + Ifges(4,3) * t550;
t453 = mrSges(4,2) * t510 - mrSges(4,3) * t489 + Ifges(4,1) * t509 + Ifges(4,4) * t508 + Ifges(4,5) * t549 - pkin(7) * t468 - t553 * t472 + t557 * t473 - t530 * t513 - t550 * t514;
t564 = mrSges(5,1) * t483 - mrSges(5,2) * t484 + Ifges(5,5) * t493 + Ifges(5,6) * t492 + Ifges(5,3) * t507 + t520 * t495 - t519 * t496;
t454 = -mrSges(4,1) * t510 + mrSges(4,3) * t490 + Ifges(4,4) * t509 + Ifges(4,2) * t508 + Ifges(4,6) * t549 - pkin(3) * t468 - t531 * t513 + t550 * t515 - t564;
t527 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t555 + Ifges(3,6) * t559) * qJD(1);
t445 = -mrSges(3,1) * t532 + mrSges(3,3) * t522 + Ifges(3,4) * t538 + Ifges(3,2) * t539 + Ifges(3,6) * qJDD(2) - pkin(2) * t566 + pkin(6) * t573 + qJD(2) * t529 + t554 * t453 + t558 * t454 - t527 * t578;
t447 = mrSges(3,2) * t532 - mrSges(3,3) * t521 + Ifges(3,1) * t538 + Ifges(3,4) * t539 + Ifges(3,5) * qJDD(2) - pkin(6) * t460 - qJD(2) * t528 + t558 * t453 - t554 * t454 + t527 * t577;
t568 = mrSges(2,1) * t544 - mrSges(2,2) * t545 + Ifges(2,3) * qJDD(1) + pkin(1) * t563 + pkin(5) * t574 + t559 * t445 + t555 * t447;
t443 = mrSges(2,1) * g(3) + mrSges(2,3) * t545 + t561 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t452 - t582;
t442 = -mrSges(2,2) * g(3) - mrSges(2,3) * t544 + Ifges(2,5) * qJDD(1) - t561 * Ifges(2,6) - pkin(5) * t452 - t555 * t445 + t559 * t447;
t1 = [-m(1) * g(1) + t575; -m(1) * g(2) + t579; (-m(1) - m(2)) * g(3) + t452; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t579 + t560 * t442 - t556 * t443; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t575 + t556 * t442 + t560 * t443; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t568; t568; t582; -t565; t564;];
tauJB = t1;
