% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:40:58
% EndTime: 2019-05-05 13:41:03
% DurationCPUTime: 4.04s
% Computational Cost: add. (51026->275), mult. (102914->331), div. (0->0), fcn. (58887->10), ass. (0->120)
t552 = qJD(1) ^ 2;
t541 = sin(pkin(10));
t543 = cos(pkin(10));
t546 = sin(qJ(5));
t549 = cos(qJ(5));
t563 = t541 * t546 - t543 * t549;
t517 = t563 * qJD(1);
t564 = -t541 * t549 - t543 * t546;
t518 = t564 * qJD(1);
t575 = t518 * qJD(5);
t502 = t563 * qJDD(1) - t575;
t586 = -pkin(1) - pkin(2);
t585 = mrSges(2,1) + mrSges(3,1);
t584 = Ifges(3,4) + Ifges(2,5);
t583 = Ifges(2,6) - Ifges(3,6);
t582 = mrSges(5,1) * t543;
t581 = mrSges(5,2) * t541;
t534 = t543 ^ 2;
t580 = t534 * t552;
t547 = sin(qJ(1));
t550 = cos(qJ(1));
t522 = -t550 * g(1) - t547 * g(2);
t559 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t522;
t514 = -pkin(1) * t552 + t559;
t509 = t586 * t552 + t559;
t521 = t547 * g(1) - t550 * g(2);
t558 = -t552 * qJ(2) + qJDD(2) - t521;
t513 = t586 * qJDD(1) + t558;
t542 = sin(pkin(9));
t544 = cos(pkin(9));
t491 = t544 * t509 + t542 * t513;
t488 = -pkin(3) * t552 - qJDD(1) * qJ(4) + t491;
t538 = g(3) + qJDD(3);
t574 = qJD(1) * qJD(4);
t578 = t543 * t538 + 0.2e1 * t541 * t574;
t474 = (pkin(4) * t543 * t552 + pkin(7) * qJDD(1) - t488) * t541 + t578;
t478 = t541 * t538 + (t488 - 0.2e1 * t574) * t543;
t573 = qJDD(1) * t543;
t475 = -pkin(4) * t580 - pkin(7) * t573 + t478;
t471 = t546 * t474 + t549 * t475;
t498 = -mrSges(6,1) * t517 + mrSges(6,2) * t518;
t511 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t518;
t501 = -pkin(5) * t517 - pkin(8) * t518;
t551 = qJD(5) ^ 2;
t468 = -pkin(5) * t551 + qJDD(5) * pkin(8) + t501 * t517 + t471;
t533 = t541 ^ 2;
t490 = -t542 * t509 + t544 * t513;
t560 = qJDD(1) * pkin(3) + qJDD(4) - t490;
t476 = pkin(4) * t573 + (-qJ(4) + (-t533 - t534) * pkin(7)) * t552 + t560;
t576 = t517 * qJD(5);
t503 = t564 * qJDD(1) + t576;
t469 = (-t503 - t576) * pkin(8) + (-t502 + t575) * pkin(5) + t476;
t545 = sin(qJ(6));
t548 = cos(qJ(6));
t465 = -t468 * t545 + t469 * t548;
t506 = qJD(5) * t548 - t518 * t545;
t485 = qJD(6) * t506 + qJDD(5) * t545 + t503 * t548;
t507 = qJD(5) * t545 + t518 * t548;
t489 = -mrSges(7,1) * t506 + mrSges(7,2) * t507;
t515 = qJD(6) - t517;
t492 = -mrSges(7,2) * t515 + mrSges(7,3) * t506;
t500 = qJDD(6) - t502;
t463 = m(7) * t465 + mrSges(7,1) * t500 - mrSges(7,3) * t485 - t489 * t507 + t492 * t515;
t466 = t468 * t548 + t469 * t545;
t484 = -qJD(6) * t507 + qJDD(5) * t548 - t503 * t545;
t493 = mrSges(7,1) * t515 - mrSges(7,3) * t507;
t464 = m(7) * t466 - mrSges(7,2) * t500 + mrSges(7,3) * t484 + t489 * t506 - t493 * t515;
t568 = -t463 * t545 + t548 * t464;
t454 = m(6) * t471 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t502 - qJD(5) * t511 + t498 * t517 + t568;
t470 = t474 * t549 - t475 * t546;
t510 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t517;
t467 = -qJDD(5) * pkin(5) - pkin(8) * t551 + t501 * t518 - t470;
t556 = -m(7) * t467 + t484 * mrSges(7,1) - mrSges(7,2) * t485 + t506 * t492 - t493 * t507;
t459 = m(6) * t470 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t503 + qJD(5) * t510 - t498 * t518 + t556;
t450 = t546 * t454 + t549 * t459;
t477 = -t488 * t541 + t578;
t562 = mrSges(5,3) * qJDD(1) + t552 * (-t581 + t582);
t448 = m(5) * t477 + t562 * t541 + t450;
t569 = t549 * t454 - t546 * t459;
t449 = m(5) * t478 - t562 * t543 + t569;
t570 = -t448 * t541 + t543 * t449;
t442 = m(4) * t491 - mrSges(4,1) * t552 + qJDD(1) * mrSges(4,2) + t570;
t487 = -t552 * qJ(4) + t560;
t455 = t548 * t463 + t545 * t464;
t554 = m(6) * t476 - t502 * mrSges(6,1) + mrSges(6,2) * t503 - t517 * t510 + t511 * t518 + t455;
t553 = -m(5) * t487 + qJDD(1) * t581 - t554 + (t533 * t552 + t580) * mrSges(5,3);
t451 = t553 + (-mrSges(4,1) - t582) * qJDD(1) + m(4) * t490 - mrSges(4,2) * t552;
t571 = t544 * t442 - t542 * t451;
t561 = m(3) * t514 + qJDD(1) * mrSges(3,3) + t571;
t438 = m(2) * t522 - qJDD(1) * mrSges(2,2) - t585 * t552 + t561;
t440 = t442 * t542 + t451 * t544;
t516 = -qJDD(1) * pkin(1) + t558;
t555 = -m(3) * t516 + qJDD(1) * mrSges(3,1) + t552 * mrSges(3,3) - t440;
t439 = m(2) * t521 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t552 + t555;
t579 = t547 * t438 + t550 * t439;
t565 = -Ifges(5,5) * t541 - Ifges(5,6) * t543;
t577 = t552 * t565;
t572 = t550 * t438 - t439 * t547;
t567 = -Ifges(5,1) * t541 - Ifges(5,4) * t543;
t566 = -Ifges(5,4) * t541 - Ifges(5,2) * t543;
t444 = t543 * t448 + t541 * t449;
t557 = -m(4) * t538 - t444;
t496 = Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * qJD(5);
t495 = Ifges(6,4) * t518 + Ifges(6,2) * t517 + Ifges(6,6) * qJD(5);
t494 = Ifges(6,5) * t518 + Ifges(6,6) * t517 + Ifges(6,3) * qJD(5);
t481 = Ifges(7,1) * t507 + Ifges(7,4) * t506 + Ifges(7,5) * t515;
t480 = Ifges(7,4) * t507 + Ifges(7,2) * t506 + Ifges(7,6) * t515;
t479 = Ifges(7,5) * t507 + Ifges(7,6) * t506 + Ifges(7,3) * t515;
t457 = mrSges(7,2) * t467 - mrSges(7,3) * t465 + Ifges(7,1) * t485 + Ifges(7,4) * t484 + Ifges(7,5) * t500 + t479 * t506 - t480 * t515;
t456 = -mrSges(7,1) * t467 + mrSges(7,3) * t466 + Ifges(7,4) * t485 + Ifges(7,2) * t484 + Ifges(7,6) * t500 - t479 * t507 + t481 * t515;
t446 = -mrSges(6,1) * t476 - mrSges(7,1) * t465 + mrSges(7,2) * t466 + mrSges(6,3) * t471 + Ifges(6,4) * t503 - Ifges(7,5) * t485 + Ifges(6,2) * t502 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t484 - Ifges(7,3) * t500 - pkin(5) * t455 + qJD(5) * t496 - t480 * t507 + t481 * t506 - t494 * t518;
t445 = mrSges(6,2) * t476 - mrSges(6,3) * t470 + Ifges(6,1) * t503 + Ifges(6,4) * t502 + Ifges(6,5) * qJDD(5) - pkin(8) * t455 - qJD(5) * t495 - t456 * t545 + t457 * t548 + t494 * t517;
t443 = -m(3) * g(3) + t557;
t434 = mrSges(5,2) * t487 - mrSges(5,3) * t477 - pkin(7) * t450 + t567 * qJDD(1) + t549 * t445 - t546 * t446 - t543 * t577;
t433 = -mrSges(5,1) * t487 + mrSges(5,3) * t478 - pkin(4) * t554 + pkin(7) * t569 + t566 * qJDD(1) + t546 * t445 + t549 * t446 + t541 * t577;
t432 = -pkin(3) * t444 + mrSges(4,3) * t491 - mrSges(4,1) * t538 - pkin(4) * t450 - mrSges(5,1) * t477 + mrSges(5,2) * t478 - t545 * t457 - t548 * t456 - pkin(5) * t556 - pkin(8) * t568 - Ifges(6,5) * t503 - Ifges(6,6) * t502 - Ifges(6,3) * qJDD(5) - t518 * t495 + t517 * t496 - mrSges(6,1) * t470 + mrSges(6,2) * t471 + (-Ifges(4,6) - t565) * qJDD(1) + (t541 * t566 - t543 * t567 + Ifges(4,5)) * t552;
t431 = mrSges(4,2) * t538 - mrSges(4,3) * t490 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t552 - qJ(4) * t444 - t433 * t541 + t434 * t543;
t430 = mrSges(3,2) * t516 - mrSges(2,3) * t521 - qJ(2) * t443 - qJ(3) * t440 + t544 * t431 - t542 * t432 - t583 * t552 + t584 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t429 = mrSges(3,2) * t514 + mrSges(2,3) * t522 - pkin(1) * t443 - pkin(2) * t557 + t585 * g(3) - qJ(3) * t571 + t583 * qJDD(1) - t542 * t431 - t544 * t432 + t584 * t552;
t1 = [-m(1) * g(1) + t572; -m(1) * g(2) + t579; (-m(1) - m(2) - m(3)) * g(3) + t557; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t579 - t547 * t429 + t550 * t430; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t572 + t550 * t429 + t547 * t430; qJ(2) * (-mrSges(3,1) * t552 + t561) + mrSges(2,1) * t521 - mrSges(2,2) * t522 + pkin(1) * t555 - pkin(2) * t440 - mrSges(3,1) * t516 + mrSges(3,3) * t514 - qJ(4) * t570 - t541 * t434 - t543 * t433 - pkin(3) * t553 - mrSges(4,1) * t490 + mrSges(4,2) * t491 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (pkin(3) * t582 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB  = t1;
