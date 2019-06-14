% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:52:30
% EndTime: 2019-05-05 18:52:34
% DurationCPUTime: 3.52s
% Computational Cost: add. (30561->302), mult. (74084->372), div. (0->0), fcn. (54978->10), ass. (0->131)
t585 = Ifges(4,1) + Ifges(5,1);
t579 = Ifges(4,4) - Ifges(5,5);
t578 = Ifges(4,5) + Ifges(5,4);
t584 = Ifges(4,2) + Ifges(5,3);
t577 = Ifges(4,6) - Ifges(5,6);
t583 = Ifges(4,3) + Ifges(5,2);
t548 = qJD(1) ^ 2;
t582 = 2 * qJD(4);
t581 = cos(qJ(3));
t580 = -mrSges(4,3) - mrSges(5,2);
t539 = cos(pkin(10));
t531 = t539 ^ 2;
t576 = t531 * t548;
t543 = sin(qJ(1));
t546 = cos(qJ(1));
t560 = -g(1) * t546 - g(2) * t543;
t515 = -pkin(1) * t548 + qJDD(1) * qJ(2) + t560;
t538 = sin(pkin(10));
t566 = qJD(1) * qJD(2);
t563 = -t539 * g(3) - 0.2e1 * t538 * t566;
t477 = (pkin(2) * t539 * t548 - pkin(7) * qJDD(1) - t515) * t538 + t563;
t500 = -g(3) * t538 + (t515 + 0.2e1 * t566) * t539;
t565 = qJDD(1) * t539;
t478 = -pkin(2) * t576 + pkin(7) * t565 + t500;
t542 = sin(qJ(3));
t459 = t542 * t477 + t581 * t478;
t564 = t539 * t581;
t557 = t581 * t538 + t539 * t542;
t514 = t557 * qJD(1);
t567 = t514 * qJD(3);
t497 = t567 + (t538 * t542 - t564) * qJDD(1);
t504 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t514;
t569 = qJD(1) * t538;
t513 = -qJD(1) * t564 + t542 * t569;
t491 = pkin(3) * t513 - qJ(4) * t514;
t547 = qJD(3) ^ 2;
t447 = -pkin(3) * t547 + qJDD(3) * qJ(4) + qJD(3) * t582 - t513 * t491 + t459;
t505 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t514;
t458 = t581 * t477 - t542 * t478;
t448 = -qJDD(3) * pkin(3) - t547 * qJ(4) + t514 * t491 + qJDD(4) - t458;
t568 = t513 * qJD(3);
t498 = t557 * qJDD(1) - t568;
t438 = (-t498 - t568) * pkin(8) + (t513 * t514 - qJDD(3)) * pkin(4) + t448;
t507 = -qJD(3) * pkin(4) - pkin(8) * t514;
t512 = t513 ^ 2;
t440 = -pkin(4) * t512 + pkin(8) * t497 + qJD(3) * t507 + t447;
t541 = sin(qJ(5));
t545 = cos(qJ(5));
t435 = t541 * t438 + t545 * t440;
t486 = t513 * t545 - t514 * t541;
t487 = t513 * t541 + t514 * t545;
t467 = -pkin(5) * t486 - pkin(9) * t487;
t532 = -qJD(3) + qJD(5);
t528 = t532 ^ 2;
t529 = -qJDD(3) + qJDD(5);
t432 = -pkin(5) * t528 + pkin(9) * t529 + t467 * t486 + t435;
t570 = t543 * g(1) - t546 * g(2);
t511 = -qJDD(1) * pkin(1) - t548 * qJ(2) + qJDD(2) - t570;
t530 = t538 ^ 2;
t496 = -pkin(2) * t565 + t511 + (-t530 * t548 - t576) * pkin(7);
t554 = t497 * pkin(3) + t496 + (-t498 + t568) * qJ(4);
t436 = -pkin(3) * t567 - t497 * pkin(4) - t512 * pkin(8) - t554 + (t507 + t582) * t514;
t455 = -qJD(5) * t487 + t497 * t545 - t498 * t541;
t456 = qJD(5) * t486 + t497 * t541 + t498 * t545;
t433 = t436 + (t487 * t532 - t455) * pkin(5) + (-t486 * t532 - t456) * pkin(9);
t540 = sin(qJ(6));
t544 = cos(qJ(6));
t429 = -t432 * t540 + t433 * t544;
t470 = -t487 * t540 + t532 * t544;
t443 = qJD(6) * t470 + t456 * t544 + t529 * t540;
t454 = qJDD(6) - t455;
t471 = t487 * t544 + t532 * t540;
t457 = -mrSges(7,1) * t470 + mrSges(7,2) * t471;
t479 = qJD(6) - t486;
t460 = -mrSges(7,2) * t479 + mrSges(7,3) * t470;
t426 = m(7) * t429 + mrSges(7,1) * t454 - mrSges(7,3) * t443 - t457 * t471 + t460 * t479;
t430 = t432 * t544 + t433 * t540;
t442 = -qJD(6) * t471 - t456 * t540 + t529 * t544;
t461 = mrSges(7,1) * t479 - mrSges(7,3) * t471;
t427 = m(7) * t430 - mrSges(7,2) * t454 + mrSges(7,3) * t442 + t457 * t470 - t461 * t479;
t419 = -t426 * t540 + t544 * t427;
t466 = -mrSges(6,1) * t486 + mrSges(6,2) * t487;
t474 = mrSges(6,1) * t532 - mrSges(6,3) * t487;
t417 = m(6) * t435 - mrSges(6,2) * t529 + mrSges(6,3) * t455 + t466 * t486 - t474 * t532 + t419;
t434 = t438 * t545 - t440 * t541;
t431 = -pkin(5) * t529 - pkin(9) * t528 + t467 * t487 - t434;
t428 = -m(7) * t431 + t442 * mrSges(7,1) - mrSges(7,2) * t443 + t470 * t460 - t461 * t471;
t473 = -mrSges(6,2) * t532 + mrSges(6,3) * t486;
t422 = m(6) * t434 + mrSges(6,1) * t529 - mrSges(6,3) * t456 - t466 * t487 + t473 * t532 + t428;
t561 = t545 * t417 - t541 * t422;
t556 = m(5) * t447 + qJDD(3) * mrSges(5,3) + qJD(3) * t505 + t561;
t492 = mrSges(5,1) * t513 - mrSges(5,3) * t514;
t571 = -mrSges(4,1) * t513 - mrSges(4,2) * t514 - t492;
t410 = m(4) * t459 - qJDD(3) * mrSges(4,2) - qJD(3) * t504 + t580 * t497 + t571 * t513 + t556;
t503 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t513;
t413 = t541 * t417 + t545 * t422;
t506 = -mrSges(5,2) * t513 + qJD(3) * mrSges(5,3);
t555 = -m(5) * t448 + qJDD(3) * mrSges(5,1) + qJD(3) * t506 - t413;
t411 = m(4) * t458 + qJDD(3) * mrSges(4,1) + qJD(3) * t503 + t580 * t498 + t571 * t514 + t555;
t575 = t542 * t410 + t581 * t411;
t418 = t544 * t426 + t540 * t427;
t574 = -t577 * qJD(3) + t584 * t513 - t579 * t514;
t573 = -t583 * qJD(3) + t577 * t513 - t578 * t514;
t572 = t578 * qJD(3) - t579 * t513 + t585 * t514;
t562 = t581 * t410 - t542 * t411;
t559 = -mrSges(3,1) * t539 + mrSges(3,2) * t538;
t558 = mrSges(3,3) * qJDD(1) + t548 * t559;
t553 = -m(6) * t436 + t455 * mrSges(6,1) - t456 * mrSges(6,2) + t486 * t473 - t487 * t474 - t418;
t445 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t514 + t554;
t552 = m(5) * t445 + t497 * mrSges(5,1) + t513 * t506 + t553;
t450 = Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t479;
t451 = Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t479;
t551 = mrSges(7,1) * t429 - mrSges(7,2) * t430 + Ifges(7,5) * t443 + Ifges(7,6) * t442 + Ifges(7,3) * t454 + t450 * t471 - t451 * t470;
t449 = Ifges(7,5) * t471 + Ifges(7,6) * t470 + Ifges(7,3) * t479;
t420 = -mrSges(7,1) * t431 + mrSges(7,3) * t430 + Ifges(7,4) * t443 + Ifges(7,2) * t442 + Ifges(7,6) * t454 - t449 * t471 + t451 * t479;
t421 = mrSges(7,2) * t431 - mrSges(7,3) * t429 + Ifges(7,1) * t443 + Ifges(7,4) * t442 + Ifges(7,5) * t454 + t449 * t470 - t450 * t479;
t463 = Ifges(6,4) * t487 + Ifges(6,2) * t486 + Ifges(6,6) * t532;
t464 = Ifges(6,1) * t487 + Ifges(6,4) * t486 + Ifges(6,5) * t532;
t550 = mrSges(6,1) * t434 - mrSges(6,2) * t435 + Ifges(6,5) * t456 + Ifges(6,6) * t455 + Ifges(6,3) * t529 + pkin(5) * t428 + pkin(9) * t419 + t544 * t420 + t540 * t421 + t487 * t463 - t486 * t464;
t549 = t552 + t513 * t503 + m(4) * t496 + t497 * mrSges(4,1) + (t504 - t505) * t514 + (mrSges(4,2) - mrSges(5,3)) * t498;
t517 = (Ifges(3,5) * t538 + Ifges(3,6) * t539) * qJD(1);
t499 = -t538 * t515 + t563;
t462 = Ifges(6,5) * t487 + Ifges(6,6) * t486 + Ifges(6,3) * t532;
t415 = -t498 * mrSges(5,3) - t514 * t505 + t552;
t414 = t549 + m(3) * t511 + t559 * qJDD(1) + (-t530 - t531) * t548 * mrSges(3,3);
t412 = t498 * mrSges(5,2) + t514 * t492 - t555;
t406 = -mrSges(6,1) * t436 + mrSges(6,3) * t435 + Ifges(6,4) * t456 + Ifges(6,2) * t455 + Ifges(6,6) * t529 - pkin(5) * t418 - t462 * t487 + t464 * t532 - t551;
t405 = mrSges(6,2) * t436 - mrSges(6,3) * t434 + Ifges(6,1) * t456 + Ifges(6,4) * t455 + Ifges(6,5) * t529 - pkin(9) * t418 - t420 * t540 + t421 * t544 + t462 * t486 - t463 * t532;
t404 = mrSges(4,2) * t496 + mrSges(5,2) * t448 - mrSges(4,3) * t458 - mrSges(5,3) * t445 - pkin(8) * t413 - qJ(4) * t415 + t574 * qJD(3) + t578 * qJDD(3) + t545 * t405 - t541 * t406 - t579 * t497 + t585 * t498 + t573 * t513;
t403 = -mrSges(4,1) * t496 - mrSges(5,1) * t445 + mrSges(5,2) * t447 + mrSges(4,3) * t459 - pkin(3) * t415 - pkin(4) * t553 - pkin(8) * t561 + t572 * qJD(3) + t577 * qJDD(3) - t541 * t405 - t545 * t406 - t584 * t497 + t579 * t498 + t573 * t514;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t570 - mrSges(2,2) * t560 + t538 * (t539 * qJD(1) * t517 + mrSges(3,2) * t511 - mrSges(3,3) * t499 + t581 * t404 - t542 * t403 - pkin(7) * t575 + (Ifges(3,1) * t538 + Ifges(3,4) * t539) * qJDD(1)) + t539 * (-t517 * t569 - mrSges(3,1) * t511 + mrSges(3,3) * t500 + t542 * t404 + t581 * t403 - pkin(2) * t549 + pkin(7) * t562 + (Ifges(3,4) * t538 + Ifges(3,2) * t539) * qJDD(1)) - pkin(1) * t414 + qJ(2) * ((m(3) * t500 + t558 * t539 + t562) * t539 + (-m(3) * t499 + t558 * t538 - t575) * t538); t414; -t550 + mrSges(4,1) * t458 - mrSges(4,2) * t459 + mrSges(5,3) * t447 - mrSges(5,1) * t448 + t583 * qJDD(3) - pkin(4) * t413 - pkin(3) * t412 - t574 * t514 + t578 * t498 + (-qJ(4) * t492 + t572) * t513 + qJ(4) * t556 + (-mrSges(5,2) * qJ(4) - t577) * t497; t412; t550; t551;];
tauJ  = t1;
