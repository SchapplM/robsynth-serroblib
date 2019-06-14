% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:24:13
% EndTime: 2019-05-07 18:24:22
% DurationCPUTime: 4.28s
% Computational Cost: add. (27801->302), mult. (55327->357), div. (0->0), fcn. (38134->8), ass. (0->116)
t553 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t578 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t568 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t528 = sin(qJ(3));
t531 = cos(qJ(3));
t529 = sin(qJ(2));
t556 = qJD(1) * t529;
t511 = qJD(2) * t531 - t528 * t556;
t512 = qJD(2) * t528 + t531 * t556;
t527 = sin(qJ(4));
t562 = cos(qJ(4));
t486 = -t562 * t511 + t527 * t512;
t487 = t527 * t511 + t562 * t512;
t532 = cos(qJ(2));
t555 = qJD(1) * t532;
t523 = qJD(3) - t555;
t521 = -qJD(4) - t523;
t577 = t578 * t486 + t553 * t487 - t568 * t521;
t470 = -mrSges(7,2) * t521 + mrSges(7,3) * t486;
t554 = qJD(1) * qJD(2);
t524 = t529 * t554;
t516 = qJDD(1) * t532 - t524;
t510 = qJDD(3) - t516;
t508 = -qJDD(4) - t510;
t535 = qJD(1) ^ 2;
t530 = sin(qJ(1));
t533 = cos(qJ(1));
t549 = t530 * g(1) - t533 * g(2);
t506 = -qJDD(1) * pkin(1) - t535 * pkin(7) - t549;
t550 = t532 * t554;
t515 = qJDD(1) * t529 + t550;
t463 = (-t515 - t550) * pkin(8) + (-t516 + t524) * pkin(2) + t506;
t543 = -g(1) * t533 - g(2) * t530;
t507 = -pkin(1) * t535 + qJDD(1) * pkin(7) + t543;
t492 = -g(3) * t529 + t532 * t507;
t514 = (-pkin(2) * t532 - pkin(8) * t529) * qJD(1);
t534 = qJD(2) ^ 2;
t468 = -pkin(2) * t534 + qJDD(2) * pkin(8) + t514 * t555 + t492;
t437 = t531 * t463 - t528 * t468;
t484 = qJD(3) * t511 + qJDD(2) * t528 + t515 * t531;
t416 = (t511 * t523 - t484) * pkin(9) + (t511 * t512 + t510) * pkin(3) + t437;
t438 = t528 * t463 + t531 * t468;
t483 = -qJD(3) * t512 + qJDD(2) * t531 - t515 * t528;
t493 = pkin(3) * t523 - pkin(9) * t512;
t509 = t511 ^ 2;
t419 = -pkin(3) * t509 + pkin(9) * t483 - t493 * t523 + t438;
t413 = t562 * t416 - t527 * t419;
t457 = pkin(4) * t486 - qJ(5) * t487;
t520 = t521 ^ 2;
t409 = t508 * pkin(4) - t520 * qJ(5) + t487 * t457 + qJDD(5) - t413;
t436 = -t486 * qJD(4) + t527 * t483 + t562 * t484;
t559 = t486 * t521;
t400 = -0.2e1 * qJD(6) * t487 + (-t436 + t559) * qJ(6) + (t486 * t487 + t508) * pkin(5) + t409;
t459 = -mrSges(7,1) * t486 + mrSges(7,2) * t487;
t544 = -m(7) * t400 + t436 * mrSges(7,3) + t487 * t459;
t398 = t508 * mrSges(7,1) + t521 * t470 - t544;
t458 = mrSges(6,1) * t486 - mrSges(6,3) * t487;
t469 = -mrSges(6,2) * t486 - mrSges(6,3) * t521;
t565 = -m(6) * t409 - t508 * mrSges(6,1) - t521 * t469;
t395 = t436 * mrSges(6,2) + t487 * t458 + t398 - t565;
t414 = t527 * t416 + t562 * t419;
t563 = -2 * qJD(5);
t408 = -pkin(4) * t520 - t508 * qJ(5) - t486 * t457 + t521 * t563 + t414;
t435 = t487 * qJD(4) - t562 * t483 + t527 * t484;
t472 = pkin(5) * t521 - qJ(6) * t487;
t485 = t486 ^ 2;
t403 = -pkin(5) * t485 + qJ(6) * t435 + 0.2e1 * qJD(6) * t486 - t472 * t521 + t408;
t473 = mrSges(7,1) * t521 - mrSges(7,3) * t487;
t475 = mrSges(6,1) * t521 + mrSges(6,2) * t487;
t552 = m(7) * t403 + t435 * mrSges(7,3) + t486 * t459;
t542 = m(6) * t408 - t508 * mrSges(6,3) - t521 * t475 + t552;
t570 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t574 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t567 = -t486 * t553 + t574 * t487 - t570 * t521;
t569 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t576 = -mrSges(6,1) * t409 - mrSges(7,1) * t400 - mrSges(5,2) * t414 - pkin(4) * t395 - pkin(5) * t398 + qJ(5) * (-t521 * t473 + t542) + mrSges(7,2) * t403 + mrSges(6,3) * t408 + mrSges(5,1) * t413 + (-mrSges(7,2) * qJ(5) - t569) * t508 + (-qJ(5) * t458 + t567) * t486 + t570 * t436 + (-mrSges(6,2) * qJ(5) - t568) * t435 + t577 * t487;
t471 = mrSges(5,2) * t521 - mrSges(5,3) * t486;
t557 = -mrSges(5,1) * t486 - mrSges(5,2) * t487 - t458;
t560 = -mrSges(5,3) - mrSges(6,2);
t389 = m(5) * t413 + (-t470 - t471) * t521 + (-mrSges(5,1) - mrSges(7,1)) * t508 + t557 * t487 + t560 * t436 + t544 + t565;
t474 = -mrSges(5,1) * t521 - mrSges(5,3) * t487;
t392 = m(5) * t414 + (-t473 + t474) * t521 + (mrSges(5,2) - mrSges(7,2)) * t508 + t557 * t486 + t560 * t435 + t542;
t387 = t562 * t389 + t527 * t392;
t478 = Ifges(4,4) * t512 + Ifges(4,2) * t511 + Ifges(4,6) * t523;
t479 = Ifges(4,1) * t512 + Ifges(4,4) * t511 + Ifges(4,5) * t523;
t575 = mrSges(4,1) * t437 - mrSges(4,2) * t438 + Ifges(4,5) * t484 + Ifges(4,6) * t483 + Ifges(4,3) * t510 + pkin(3) * t387 + t512 * t478 - t511 * t479 + t576;
t491 = -t532 * g(3) - t529 * t507;
t541 = qJDD(2) * pkin(2) + t534 * pkin(8) - t514 * t556 + t491;
t540 = t483 * pkin(3) + t509 * pkin(9) - t512 * t493 + t541;
t564 = -(t436 + t559) * qJ(5) - t540;
t551 = t568 * t486 - t570 * t487 + t569 * t521;
t488 = -mrSges(4,1) * t511 + mrSges(4,2) * t512;
t489 = -mrSges(4,2) * t523 + mrSges(4,3) * t511;
t385 = m(4) * t437 + mrSges(4,1) * t510 - mrSges(4,3) * t484 - t488 * t512 + t489 * t523 + t387;
t490 = mrSges(4,1) * t523 - mrSges(4,3) * t512;
t546 = -t389 * t527 + t562 * t392;
t386 = m(4) * t438 - mrSges(4,2) * t510 + mrSges(4,3) * t483 + t488 * t511 - t490 * t523 + t546;
t547 = -t385 * t528 + t531 * t386;
t405 = -t485 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t435 + (pkin(4) * t521 + (2 * qJD(5)) + t472) * t487 - t564;
t399 = m(7) * t405 - t435 * mrSges(7,1) + t436 * mrSges(7,2) - t486 * t470 + t487 * t473;
t381 = t531 * t385 + t528 * t386;
t411 = t487 * t563 + (-t487 * t521 + t435) * pkin(4) + t564;
t393 = m(6) * t411 + t435 * mrSges(6,1) - t436 * mrSges(6,3) + t486 * t469 - t487 * t475 - t399;
t539 = -m(5) * t540 + t435 * mrSges(5,1) + t436 * mrSges(5,2) + t486 * t471 + t487 * t474 + t393;
t538 = m(4) * t541 + t483 * mrSges(4,1) - t484 * mrSges(4,2) + t511 * t489 - t512 * t490 - t539;
t519 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t555;
t518 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t556;
t513 = (-mrSges(3,1) * t532 + mrSges(3,2) * t529) * qJD(1);
t505 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t529 + Ifges(3,4) * t532) * qJD(1);
t504 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t529 + Ifges(3,2) * t532) * qJD(1);
t477 = Ifges(4,5) * t512 + Ifges(4,6) * t511 + Ifges(4,3) * t523;
t383 = -mrSges(5,2) * t540 + mrSges(6,2) * t409 + mrSges(7,2) * t405 - mrSges(5,3) * t413 - mrSges(6,3) * t411 - mrSges(7,3) * t400 - qJ(5) * t393 - qJ(6) * t398 - t553 * t435 + t574 * t436 + t551 * t486 - t570 * t508 + t577 * t521;
t382 = mrSges(5,1) * t540 + mrSges(5,3) * t414 - mrSges(6,1) * t411 + mrSges(6,2) * t408 + mrSges(7,1) * t405 - mrSges(7,3) * t403 + pkin(5) * t399 - qJ(6) * t552 - pkin(4) * t393 + (qJ(6) * t473 - t567) * t521 + (mrSges(7,2) * qJ(6) - t568) * t508 + t551 * t487 + t553 * t436 + t578 * t435;
t380 = -mrSges(4,2) * t541 - mrSges(4,3) * t437 + Ifges(4,1) * t484 + Ifges(4,4) * t483 + Ifges(4,5) * t510 - pkin(9) * t387 - t527 * t382 + t562 * t383 + t511 * t477 - t523 * t478;
t379 = mrSges(4,1) * t541 + mrSges(4,3) * t438 + Ifges(4,4) * t484 + Ifges(4,2) * t483 + Ifges(4,6) * t510 - pkin(3) * t539 + pkin(9) * t546 + t562 * t382 + t527 * t383 - t512 * t477 + t523 * t479;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t549 - mrSges(2,2) * t543 + t529 * (mrSges(3,2) * t506 - mrSges(3,3) * t491 + Ifges(3,1) * t515 + Ifges(3,4) * t516 + Ifges(3,5) * qJDD(2) - pkin(8) * t381 - qJD(2) * t504 - t528 * t379 + t531 * t380) + t532 * (-mrSges(3,1) * t506 + mrSges(3,3) * t492 + Ifges(3,4) * t515 + Ifges(3,2) * t516 + Ifges(3,6) * qJDD(2) - pkin(2) * t381 + qJD(2) * t505 - t575) + pkin(1) * (-m(3) * t506 + t516 * mrSges(3,1) - t515 * mrSges(3,2) + (-t518 * t529 + t519 * t532) * qJD(1) - t381) + pkin(7) * (t532 * (m(3) * t492 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t516 - qJD(2) * t518 + t513 * t555 + t547) - t529 * (m(3) * t491 + qJDD(2) * mrSges(3,1) - t515 * mrSges(3,3) + qJD(2) * t519 - t513 * t556 + t538)); Ifges(3,5) * t515 + Ifges(3,6) * t516 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t491 - mrSges(3,2) * t492 + t528 * t380 + t531 * t379 + pkin(2) * t538 + pkin(8) * t547 + (t529 * t504 - t532 * t505) * qJD(1); t575; t576; t395; t399;];
tauJ  = t1;
