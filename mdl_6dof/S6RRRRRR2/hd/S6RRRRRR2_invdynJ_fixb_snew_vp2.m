% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 08:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:26:05
% EndTime: 2019-05-08 08:26:23
% DurationCPUTime: 10.29s
% Computational Cost: add. (159706->351), mult. (344174->445), div. (0->0), fcn. (262950->12), ass. (0->141)
t547 = qJD(1) ^ 2;
t567 = pkin(2) * t547;
t540 = sin(qJ(1));
t546 = cos(qJ(1));
t557 = -g(1) * t546 - g(2) * t540;
t518 = -pkin(1) * t547 + qJDD(1) * pkin(7) + t557;
t539 = sin(qJ(2));
t566 = t518 * t539;
t545 = cos(qJ(2));
t563 = qJD(1) * qJD(2);
t521 = qJDD(1) * t539 + t545 * t563;
t481 = qJDD(2) * pkin(2) - pkin(8) * t521 - t566 + (pkin(8) * t563 + t539 * t567 - g(3)) * t545;
t505 = -g(3) * t539 + t518 * t545;
t522 = qJDD(1) * t545 - t539 * t563;
t565 = qJD(1) * t539;
t525 = qJD(2) * pkin(2) - pkin(8) * t565;
t534 = t545 ^ 2;
t482 = pkin(8) * t522 - qJD(2) * t525 - t534 * t567 + t505;
t538 = sin(qJ(3));
t544 = cos(qJ(3));
t463 = t481 * t544 - t482 * t538;
t515 = (-t538 * t539 + t544 * t545) * qJD(1);
t490 = qJD(3) * t515 + t521 * t544 + t522 * t538;
t516 = (t538 * t545 + t539 * t544) * qJD(1);
t532 = qJDD(2) + qJDD(3);
t533 = qJD(2) + qJD(3);
t436 = (t515 * t533 - t490) * pkin(9) + (t515 * t516 + t532) * pkin(3) + t463;
t464 = t481 * t538 + t482 * t544;
t489 = -qJD(3) * t516 - t521 * t538 + t522 * t544;
t508 = pkin(3) * t533 - pkin(9) * t516;
t511 = t515 ^ 2;
t441 = -pkin(3) * t511 + pkin(9) * t489 - t508 * t533 + t464;
t537 = sin(qJ(4));
t543 = cos(qJ(4));
t429 = t436 * t537 + t441 * t543;
t502 = t515 * t537 + t516 * t543;
t458 = -qJD(4) * t502 + t489 * t543 - t490 * t537;
t501 = t515 * t543 - t516 * t537;
t475 = -mrSges(5,1) * t501 + mrSges(5,2) * t502;
t530 = qJD(4) + t533;
t493 = mrSges(5,1) * t530 - mrSges(5,3) * t502;
t529 = qJDD(4) + t532;
t476 = -pkin(4) * t501 - pkin(10) * t502;
t528 = t530 ^ 2;
t418 = -pkin(4) * t528 + pkin(10) * t529 + t476 * t501 + t429;
t562 = g(1) * t540 - g(2) * t546;
t556 = -qJDD(1) * pkin(1) - t562;
t491 = -pkin(2) * t522 + t525 * t565 + (-pkin(8) * t534 - pkin(7)) * t547 + t556;
t445 = -pkin(3) * t489 - pkin(9) * t511 + t508 * t516 + t491;
t459 = qJD(4) * t501 + t489 * t537 + t490 * t543;
t426 = (-t501 * t530 - t459) * pkin(10) + (t502 * t530 - t458) * pkin(4) + t445;
t536 = sin(qJ(5));
t542 = cos(qJ(5));
t413 = -t418 * t536 + t426 * t542;
t485 = -t502 * t536 + t530 * t542;
t443 = qJD(5) * t485 + t459 * t542 + t529 * t536;
t456 = qJDD(5) - t458;
t486 = t502 * t542 + t530 * t536;
t499 = qJD(5) - t501;
t411 = (t485 * t499 - t443) * pkin(11) + (t485 * t486 + t456) * pkin(5) + t413;
t414 = t418 * t542 + t426 * t536;
t442 = -qJD(5) * t486 - t459 * t536 + t529 * t542;
t470 = pkin(5) * t499 - pkin(11) * t486;
t484 = t485 ^ 2;
t412 = -pkin(5) * t484 + pkin(11) * t442 - t470 * t499 + t414;
t535 = sin(qJ(6));
t541 = cos(qJ(6));
t409 = t411 * t541 - t412 * t535;
t465 = t485 * t541 - t486 * t535;
t425 = qJD(6) * t465 + t442 * t535 + t443 * t541;
t466 = t485 * t535 + t486 * t541;
t437 = -mrSges(7,1) * t465 + mrSges(7,2) * t466;
t494 = qJD(6) + t499;
t446 = -mrSges(7,2) * t494 + mrSges(7,3) * t465;
t449 = qJDD(6) + t456;
t405 = m(7) * t409 + mrSges(7,1) * t449 - mrSges(7,3) * t425 - t437 * t466 + t446 * t494;
t410 = t411 * t535 + t412 * t541;
t424 = -qJD(6) * t466 + t442 * t541 - t443 * t535;
t447 = mrSges(7,1) * t494 - mrSges(7,3) * t466;
t406 = m(7) * t410 - mrSges(7,2) * t449 + mrSges(7,3) * t424 + t437 * t465 - t447 * t494;
t397 = t405 * t541 + t406 * t535;
t467 = -mrSges(6,1) * t485 + mrSges(6,2) * t486;
t468 = -mrSges(6,2) * t499 + mrSges(6,3) * t485;
t395 = m(6) * t413 + mrSges(6,1) * t456 - mrSges(6,3) * t443 - t467 * t486 + t468 * t499 + t397;
t469 = mrSges(6,1) * t499 - mrSges(6,3) * t486;
t558 = -t405 * t535 + t406 * t541;
t396 = m(6) * t414 - mrSges(6,2) * t456 + mrSges(6,3) * t442 + t467 * t485 - t469 * t499 + t558;
t559 = -t395 * t536 + t396 * t542;
t389 = m(5) * t429 - mrSges(5,2) * t529 + mrSges(5,3) * t458 + t475 * t501 - t493 * t530 + t559;
t428 = t436 * t543 - t441 * t537;
t492 = -mrSges(5,2) * t530 + mrSges(5,3) * t501;
t417 = -pkin(4) * t529 - pkin(10) * t528 + t476 * t502 - t428;
t415 = -pkin(5) * t442 - pkin(11) * t484 + t470 * t486 + t417;
t555 = m(7) * t415 - mrSges(7,1) * t424 + mrSges(7,2) * t425 - t446 * t465 + t447 * t466;
t551 = -m(6) * t417 + mrSges(6,1) * t442 - mrSges(6,2) * t443 + t468 * t485 - t469 * t486 - t555;
t401 = m(5) * t428 + mrSges(5,1) * t529 - mrSges(5,3) * t459 - t475 * t502 + t492 * t530 + t551;
t382 = t389 * t537 + t401 * t543;
t503 = -mrSges(4,1) * t515 + mrSges(4,2) * t516;
t506 = -mrSges(4,2) * t533 + mrSges(4,3) * t515;
t379 = m(4) * t463 + mrSges(4,1) * t532 - mrSges(4,3) * t490 - t503 * t516 + t506 * t533 + t382;
t507 = mrSges(4,1) * t533 - mrSges(4,3) * t516;
t560 = t389 * t543 - t401 * t537;
t380 = m(4) * t464 - mrSges(4,2) * t532 + mrSges(4,3) * t489 + t503 * t515 - t507 * t533 + t560;
t374 = t379 * t544 + t380 * t538;
t391 = t395 * t542 + t396 * t536;
t564 = qJD(1) * t545;
t561 = -t379 * t538 + t380 * t544;
t554 = -m(5) * t445 + mrSges(5,1) * t458 - mrSges(5,2) * t459 + t492 * t501 - t493 * t502 - t391;
t432 = Ifges(7,4) * t466 + Ifges(7,2) * t465 + Ifges(7,6) * t494;
t433 = Ifges(7,1) * t466 + Ifges(7,4) * t465 + Ifges(7,5) * t494;
t553 = -mrSges(7,1) * t409 + mrSges(7,2) * t410 - Ifges(7,5) * t425 - Ifges(7,6) * t424 - Ifges(7,3) * t449 - t432 * t466 + t433 * t465;
t431 = Ifges(7,5) * t466 + Ifges(7,6) * t465 + Ifges(7,3) * t494;
t398 = -mrSges(7,1) * t415 + mrSges(7,3) * t410 + Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t449 - t431 * t466 + t433 * t494;
t399 = mrSges(7,2) * t415 - mrSges(7,3) * t409 + Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t449 + t431 * t465 - t432 * t494;
t450 = Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * t499;
t452 = Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * t499;
t384 = -mrSges(6,1) * t417 + mrSges(6,3) * t414 + Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t456 - pkin(5) * t555 + pkin(11) * t558 + t541 * t398 + t535 * t399 - t486 * t450 + t499 * t452;
t451 = Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * t499;
t386 = mrSges(6,2) * t417 - mrSges(6,3) * t413 + Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t456 - pkin(11) * t397 - t398 * t535 + t399 * t541 + t450 * t485 - t451 * t499;
t472 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t530;
t473 = Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t530;
t552 = mrSges(5,1) * t428 - mrSges(5,2) * t429 + Ifges(5,5) * t459 + Ifges(5,6) * t458 + Ifges(5,3) * t529 + pkin(4) * t551 + pkin(10) * t559 + t384 * t542 + t386 * t536 + t472 * t502 - t501 * t473;
t550 = m(4) * t491 - mrSges(4,1) * t489 + mrSges(4,2) * t490 - t506 * t515 + t507 * t516 - t554;
t496 = Ifges(4,4) * t516 + Ifges(4,2) * t515 + Ifges(4,6) * t533;
t497 = Ifges(4,1) * t516 + Ifges(4,4) * t515 + Ifges(4,5) * t533;
t549 = mrSges(4,1) * t463 - mrSges(4,2) * t464 + Ifges(4,5) * t490 + Ifges(4,6) * t489 + Ifges(4,3) * t532 + pkin(3) * t382 + t496 * t516 - t515 * t497 + t552;
t548 = mrSges(6,1) * t413 - mrSges(6,2) * t414 + Ifges(6,5) * t443 + Ifges(6,6) * t442 + Ifges(6,3) * t456 + pkin(5) * t397 + t451 * t486 - t452 * t485 - t553;
t524 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t564;
t523 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t565;
t520 = (-mrSges(3,1) * t545 + mrSges(3,2) * t539) * qJD(1);
t517 = -pkin(7) * t547 + t556;
t514 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t539 + Ifges(3,4) * t545) * qJD(1);
t513 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t539 + Ifges(3,2) * t545) * qJD(1);
t504 = -g(3) * t545 - t566;
t495 = Ifges(4,5) * t516 + Ifges(4,6) * t515 + Ifges(4,3) * t533;
t471 = Ifges(5,5) * t502 + Ifges(5,6) * t501 + Ifges(5,3) * t530;
t375 = -mrSges(5,1) * t445 + mrSges(5,3) * t429 + Ifges(5,4) * t459 + Ifges(5,2) * t458 + Ifges(5,6) * t529 - pkin(4) * t391 - t471 * t502 + t473 * t530 - t548;
t373 = mrSges(5,2) * t445 - mrSges(5,3) * t428 + Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * t529 - pkin(10) * t391 - t384 * t536 + t386 * t542 + t471 * t501 - t472 * t530;
t372 = mrSges(4,2) * t491 - mrSges(4,3) * t463 + Ifges(4,1) * t490 + Ifges(4,4) * t489 + Ifges(4,5) * t532 - pkin(9) * t382 + t373 * t543 - t375 * t537 + t495 * t515 - t496 * t533;
t371 = -mrSges(4,1) * t491 + mrSges(4,3) * t464 + Ifges(4,4) * t490 + Ifges(4,2) * t489 + Ifges(4,6) * t532 + pkin(3) * t554 + pkin(9) * t560 + t537 * t373 + t543 * t375 - t516 * t495 + t533 * t497;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t562 - mrSges(2,2) * t557 + t539 * (mrSges(3,2) * t517 - mrSges(3,3) * t504 + Ifges(3,1) * t521 + Ifges(3,4) * t522 + Ifges(3,5) * qJDD(2) - pkin(8) * t374 - qJD(2) * t513 - t538 * t371 + t544 * t372) + t545 * (-mrSges(3,1) * t517 + mrSges(3,3) * t505 + Ifges(3,4) * t521 + Ifges(3,2) * t522 + Ifges(3,6) * qJDD(2) - pkin(2) * t550 + pkin(8) * t561 + qJD(2) * t514 + t544 * t371 + t538 * t372) + pkin(1) * ((-t523 * t539 + t524 * t545) * qJD(1) - t550 - m(3) * t517 + mrSges(3,1) * t522 - mrSges(3,2) * t521) + pkin(7) * (t545 * (m(3) * t505 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t522 - qJD(2) * t523 + t520 * t564 + t561) - t539 * (m(3) * t504 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t521 + qJD(2) * t524 - t520 * t565 + t374)); t549 + pkin(2) * t374 + Ifges(3,3) * qJDD(2) + (t539 * t513 - t545 * t514) * qJD(1) + Ifges(3,6) * t522 + Ifges(3,5) * t521 + mrSges(3,1) * t504 - mrSges(3,2) * t505; t549; t552; t548; -t553;];
tauJ  = t1;
