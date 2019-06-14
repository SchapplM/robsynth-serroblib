% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:47:58
% EndTime: 2019-05-08 04:48:12
% DurationCPUTime: 6.84s
% Computational Cost: add. (73122->327), mult. (145341->402), div. (0->0), fcn. (105595->10), ass. (0->131)
t570 = Ifges(6,1) + Ifges(7,1);
t560 = Ifges(6,4) - Ifges(7,5);
t568 = Ifges(7,4) + Ifges(6,5);
t569 = Ifges(6,2) + Ifges(7,3);
t567 = Ifges(6,6) - Ifges(7,6);
t566 = -Ifges(6,3) - Ifges(7,2);
t532 = sin(qJ(3));
t533 = sin(qJ(2));
t536 = cos(qJ(3));
t537 = cos(qJ(2));
t512 = (t532 * t537 + t533 * t536) * qJD(1);
t528 = qJD(2) + qJD(3);
t531 = sin(qJ(4));
t535 = cos(qJ(4));
t497 = -t512 * t531 + t528 * t535;
t498 = t512 * t535 + t528 * t531;
t530 = sin(qJ(5));
t563 = cos(qJ(5));
t466 = -t563 * t497 + t530 * t498;
t467 = t530 * t497 + t563 * t498;
t555 = qJD(1) * t537;
t556 = qJD(1) * t533;
t511 = -t532 * t556 + t536 * t555;
t507 = qJD(4) - t511;
t505 = qJD(5) + t507;
t565 = t466 * t569 - t467 * t560 - t505 * t567;
t564 = -t560 * t466 + t467 * t570 + t568 * t505;
t539 = qJD(1) ^ 2;
t562 = pkin(2) * t539;
t561 = -mrSges(6,3) - mrSges(7,2);
t534 = sin(qJ(1));
t538 = cos(qJ(1));
t547 = -g(1) * t538 - g(2) * t534;
t514 = -pkin(1) * t539 + qJDD(1) * pkin(7) + t547;
t559 = t533 * t514;
t554 = qJD(1) * qJD(2);
t518 = qJDD(1) * t533 + t537 * t554;
t473 = qJDD(2) * pkin(2) - t518 * pkin(8) - t559 + (pkin(8) * t554 + t533 * t562 - g(3)) * t537;
t500 = -g(3) * t533 + t537 * t514;
t519 = qJDD(1) * t537 - t533 * t554;
t522 = qJD(2) * pkin(2) - pkin(8) * t556;
t529 = t537 ^ 2;
t474 = pkin(8) * t519 - qJD(2) * t522 - t529 * t562 + t500;
t450 = t532 * t473 + t536 * t474;
t484 = -t512 * qJD(3) - t532 * t518 + t519 * t536;
t493 = -mrSges(4,1) * t511 + mrSges(4,2) * t512;
t502 = mrSges(4,1) * t528 - mrSges(4,3) * t512;
t527 = qJDD(2) + qJDD(3);
t485 = qJD(3) * t511 + t518 * t536 + t519 * t532;
t552 = t534 * g(1) - t538 * g(2);
t546 = -qJDD(1) * pkin(1) - t552;
t486 = -t519 * pkin(2) + t522 * t556 + (-pkin(8) * t529 - pkin(7)) * t539 + t546;
t432 = (-t511 * t528 - t485) * pkin(9) + (t512 * t528 - t484) * pkin(3) + t486;
t494 = -pkin(3) * t511 - pkin(9) * t512;
t526 = t528 ^ 2;
t435 = -pkin(3) * t526 + pkin(9) * t527 + t494 * t511 + t450;
t416 = t535 * t432 - t531 * t435;
t456 = qJD(4) * t497 + t485 * t535 + t527 * t531;
t483 = qJDD(4) - t484;
t413 = (t497 * t507 - t456) * pkin(10) + (t497 * t498 + t483) * pkin(4) + t416;
t417 = t531 * t432 + t535 * t435;
t455 = -qJD(4) * t498 - t485 * t531 + t527 * t535;
t489 = pkin(4) * t507 - pkin(10) * t498;
t496 = t497 ^ 2;
t415 = -pkin(4) * t496 + pkin(10) * t455 - t489 * t507 + t417;
t409 = t530 * t413 + t563 * t415;
t426 = t467 * qJD(5) - t563 * t455 + t530 * t456;
t459 = mrSges(6,1) * t505 - mrSges(6,3) * t467;
t479 = qJDD(5) + t483;
t445 = pkin(5) * t466 - qJ(6) * t467;
t503 = t505 ^ 2;
t405 = -pkin(5) * t503 + qJ(6) * t479 + 0.2e1 * qJD(6) * t505 - t445 * t466 + t409;
t460 = -mrSges(7,1) * t505 + mrSges(7,2) * t467;
t553 = m(7) * t405 + t479 * mrSges(7,3) + t505 * t460;
t446 = mrSges(7,1) * t466 - mrSges(7,3) * t467;
t557 = -mrSges(6,1) * t466 - mrSges(6,2) * t467 - t446;
t392 = m(6) * t409 - t479 * mrSges(6,2) + t561 * t426 - t505 * t459 + t557 * t466 + t553;
t408 = t563 * t413 - t530 * t415;
t427 = -t466 * qJD(5) + t530 * t455 + t563 * t456;
t458 = -mrSges(6,2) * t505 - mrSges(6,3) * t466;
t406 = -t479 * pkin(5) - t503 * qJ(6) + t467 * t445 + qJDD(6) - t408;
t457 = -mrSges(7,2) * t466 + mrSges(7,3) * t505;
t548 = -m(7) * t406 + t479 * mrSges(7,1) + t505 * t457;
t394 = m(6) * t408 + t479 * mrSges(6,1) + t561 * t427 + t505 * t458 + t557 * t467 + t548;
t389 = t530 * t392 + t563 * t394;
t471 = -mrSges(5,1) * t497 + mrSges(5,2) * t498;
t487 = -mrSges(5,2) * t507 + mrSges(5,3) * t497;
t385 = m(5) * t416 + mrSges(5,1) * t483 - mrSges(5,3) * t456 - t471 * t498 + t487 * t507 + t389;
t488 = mrSges(5,1) * t507 - mrSges(5,3) * t498;
t549 = t563 * t392 - t394 * t530;
t386 = m(5) * t417 - mrSges(5,2) * t483 + mrSges(5,3) * t455 + t471 * t497 - t488 * t507 + t549;
t550 = -t385 * t531 + t535 * t386;
t379 = m(4) * t450 - mrSges(4,2) * t527 + mrSges(4,3) * t484 + t493 * t511 - t502 * t528 + t550;
t449 = t536 * t473 - t532 * t474;
t501 = -mrSges(4,2) * t528 + mrSges(4,3) * t511;
t434 = -t527 * pkin(3) - t526 * pkin(9) + t512 * t494 - t449;
t418 = -t455 * pkin(4) - t496 * pkin(10) + t498 * t489 + t434;
t411 = -0.2e1 * qJD(6) * t467 + (t466 * t505 - t427) * qJ(6) + (t467 * t505 + t426) * pkin(5) + t418;
t402 = m(7) * t411 + t426 * mrSges(7,1) - t427 * mrSges(7,3) + t466 * t457 - t467 * t460;
t543 = m(6) * t418 + t426 * mrSges(6,1) + mrSges(6,2) * t427 + t466 * t458 + t459 * t467 + t402;
t541 = -m(5) * t434 + t455 * mrSges(5,1) - mrSges(5,2) * t456 + t497 * t487 - t488 * t498 - t543;
t396 = m(4) * t449 + mrSges(4,1) * t527 - mrSges(4,3) * t485 - t493 * t512 + t501 * t528 + t541;
t376 = t532 * t379 + t536 * t396;
t381 = t535 * t385 + t531 * t386;
t558 = t567 * t466 - t467 * t568 + t566 * t505;
t551 = t536 * t379 - t396 * t532;
t387 = -mrSges(6,1) * t418 - mrSges(7,1) * t411 + mrSges(7,2) * t405 + mrSges(6,3) * t409 - pkin(5) * t402 - t426 * t569 + t560 * t427 + t558 * t467 + t567 * t479 + t564 * t505;
t388 = mrSges(6,2) * t418 + mrSges(7,2) * t406 - mrSges(6,3) * t408 - mrSges(7,3) * t411 - qJ(6) * t402 - t560 * t426 + t427 * t570 + t558 * t466 + t568 * t479 + t565 * t505;
t461 = Ifges(5,5) * t498 + Ifges(5,6) * t497 + Ifges(5,3) * t507;
t463 = Ifges(5,1) * t498 + Ifges(5,4) * t497 + Ifges(5,5) * t507;
t373 = -mrSges(5,1) * t434 + mrSges(5,3) * t417 + Ifges(5,4) * t456 + Ifges(5,2) * t455 + Ifges(5,6) * t483 - pkin(4) * t543 + pkin(10) * t549 + t563 * t387 + t530 * t388 - t498 * t461 + t507 * t463;
t462 = Ifges(5,4) * t498 + Ifges(5,2) * t497 + Ifges(5,6) * t507;
t375 = mrSges(5,2) * t434 - mrSges(5,3) * t416 + Ifges(5,1) * t456 + Ifges(5,4) * t455 + Ifges(5,5) * t483 - pkin(10) * t389 - t530 * t387 + t563 * t388 + t497 * t461 - t507 * t462;
t491 = Ifges(4,4) * t512 + Ifges(4,2) * t511 + Ifges(4,6) * t528;
t492 = Ifges(4,1) * t512 + Ifges(4,4) * t511 + Ifges(4,5) * t528;
t545 = mrSges(4,1) * t449 - mrSges(4,2) * t450 + Ifges(4,5) * t485 + Ifges(4,6) * t484 + Ifges(4,3) * t527 + pkin(3) * t541 + pkin(9) * t550 + t535 * t373 + t531 * t375 + t512 * t491 - t511 * t492;
t544 = m(4) * t486 - t484 * mrSges(4,1) + t485 * mrSges(4,2) - t511 * t501 + t512 * t502 + t381;
t401 = t427 * mrSges(7,2) + t467 * t446 - t548;
t542 = -mrSges(6,1) * t408 + mrSges(7,1) * t406 + mrSges(6,2) * t409 - mrSges(7,3) * t405 + pkin(5) * t401 - qJ(6) * t553 + t566 * t479 + t565 * t467 + (qJ(6) * t446 - t564) * t466 - t568 * t427 + (mrSges(7,2) * qJ(6) + t567) * t426;
t540 = mrSges(5,1) * t416 - mrSges(5,2) * t417 + Ifges(5,5) * t456 + Ifges(5,6) * t455 + Ifges(5,3) * t483 + pkin(4) * t389 + t498 * t462 - t497 * t463 - t542;
t521 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t555;
t520 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t556;
t517 = (-mrSges(3,1) * t537 + mrSges(3,2) * t533) * qJD(1);
t513 = -t539 * pkin(7) + t546;
t510 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t533 + Ifges(3,4) * t537) * qJD(1);
t509 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t533 + Ifges(3,2) * t537) * qJD(1);
t499 = -t537 * g(3) - t559;
t490 = Ifges(4,5) * t512 + Ifges(4,6) * t511 + Ifges(4,3) * t528;
t371 = -mrSges(4,1) * t486 + mrSges(4,3) * t450 + Ifges(4,4) * t485 + Ifges(4,2) * t484 + Ifges(4,6) * t527 - pkin(3) * t381 - t512 * t490 + t528 * t492 - t540;
t370 = mrSges(4,2) * t486 - mrSges(4,3) * t449 + Ifges(4,1) * t485 + Ifges(4,4) * t484 + Ifges(4,5) * t527 - pkin(9) * t381 - t373 * t531 + t375 * t535 + t490 * t511 - t491 * t528;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t552 - mrSges(2,2) * t547 + t533 * (mrSges(3,2) * t513 - mrSges(3,3) * t499 + Ifges(3,1) * t518 + Ifges(3,4) * t519 + Ifges(3,5) * qJDD(2) - pkin(8) * t376 - qJD(2) * t509 + t536 * t370 - t532 * t371) + t537 * (-mrSges(3,1) * t513 + mrSges(3,3) * t500 + Ifges(3,4) * t518 + Ifges(3,2) * t519 + Ifges(3,6) * qJDD(2) - pkin(2) * t544 + pkin(8) * t551 + qJD(2) * t510 + t532 * t370 + t536 * t371) + pkin(1) * (-m(3) * t513 + t519 * mrSges(3,1) - t518 * mrSges(3,2) + (-t520 * t533 + t521 * t537) * qJD(1) - t544) + pkin(7) * (t537 * (m(3) * t500 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t519 - qJD(2) * t520 + t517 * t555 + t551) - t533 * (m(3) * t499 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t518 + qJD(2) * t521 - t517 * t556 + t376)); (t533 * t509 - t537 * t510) * qJD(1) + Ifges(3,5) * t518 + Ifges(3,6) * t519 + mrSges(3,1) * t499 - mrSges(3,2) * t500 + pkin(2) * t376 + t545 + Ifges(3,3) * qJDD(2); t545; t540; -t542; t401;];
tauJ  = t1;
