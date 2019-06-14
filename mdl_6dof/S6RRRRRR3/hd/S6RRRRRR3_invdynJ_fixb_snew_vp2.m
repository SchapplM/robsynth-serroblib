% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR3
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
% Datum: 2019-05-08 09:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:54:01
% EndTime: 2019-05-08 08:54:17
% DurationCPUTime: 11.37s
% Computational Cost: add. (188746->352), mult. (375789->444), div. (0->0), fcn. (280589->12), ass. (0->141)
t549 = qJD(1) ^ 2;
t569 = pkin(2) * t549;
t542 = sin(qJ(1));
t548 = cos(qJ(1));
t559 = -g(1) * t548 - g(2) * t542;
t521 = -pkin(1) * t549 + qJDD(1) * pkin(7) + t559;
t541 = sin(qJ(2));
t568 = t521 * t541;
t547 = cos(qJ(2));
t565 = qJD(1) * qJD(2);
t525 = qJDD(1) * t541 + t547 * t565;
t482 = qJDD(2) * pkin(2) - pkin(8) * t525 - t568 + (pkin(8) * t565 + t541 * t569 - g(3)) * t547;
t507 = -g(3) * t541 + t547 * t521;
t526 = qJDD(1) * t547 - t541 * t565;
t567 = qJD(1) * t541;
t529 = qJD(2) * pkin(2) - pkin(8) * t567;
t536 = t547 ^ 2;
t484 = pkin(8) * t526 - qJD(2) * t529 - t536 * t569 + t507;
t540 = sin(qJ(3));
t546 = cos(qJ(3));
t462 = t540 * t482 + t546 * t484;
t519 = (t540 * t547 + t541 * t546) * qJD(1);
t492 = -t519 * qJD(3) - t540 * t525 + t526 * t546;
t566 = qJD(1) * t547;
t518 = -t540 * t567 + t546 * t566;
t501 = -mrSges(4,1) * t518 + mrSges(4,2) * t519;
t535 = qJD(2) + qJD(3);
t509 = mrSges(4,1) * t535 - mrSges(4,3) * t519;
t534 = qJDD(2) + qJDD(3);
t493 = qJD(3) * t518 + t525 * t546 + t526 * t540;
t564 = g(1) * t542 - t548 * g(2);
t557 = -qJDD(1) * pkin(1) - t564;
t494 = -pkin(2) * t526 + t529 * t567 + (-pkin(8) * t536 - pkin(7)) * t549 + t557;
t447 = (-t518 * t535 - t493) * pkin(9) + (t519 * t535 - t492) * pkin(3) + t494;
t502 = -pkin(3) * t518 - pkin(9) * t519;
t533 = t535 ^ 2;
t450 = -pkin(3) * t533 + pkin(9) * t534 + t502 * t518 + t462;
t539 = sin(qJ(4));
t545 = cos(qJ(4));
t430 = t545 * t447 - t450 * t539;
t504 = -t519 * t539 + t535 * t545;
t465 = qJD(4) * t504 + t493 * t545 + t534 * t539;
t491 = qJDD(4) - t492;
t505 = t519 * t545 + t535 * t539;
t514 = qJD(4) - t518;
t426 = (t504 * t514 - t465) * pkin(10) + (t504 * t505 + t491) * pkin(4) + t430;
t431 = t539 * t447 + t545 * t450;
t464 = -qJD(4) * t505 - t493 * t539 + t534 * t545;
t497 = pkin(4) * t514 - pkin(10) * t505;
t503 = t504 ^ 2;
t428 = -pkin(4) * t503 + pkin(10) * t464 - t497 * t514 + t431;
t538 = sin(qJ(5));
t544 = cos(qJ(5));
t414 = t544 * t426 - t428 * t538;
t475 = t504 * t544 - t505 * t538;
t443 = qJD(5) * t475 + t464 * t538 + t465 * t544;
t476 = t504 * t538 + t505 * t544;
t487 = qJDD(5) + t491;
t512 = qJD(5) + t514;
t411 = (t475 * t512 - t443) * pkin(11) + (t475 * t476 + t487) * pkin(5) + t414;
t415 = t538 * t426 + t544 * t428;
t442 = -qJD(5) * t476 + t464 * t544 - t465 * t538;
t468 = pkin(5) * t512 - pkin(11) * t476;
t474 = t475 ^ 2;
t412 = -pkin(5) * t474 + pkin(11) * t442 - t468 * t512 + t415;
t537 = sin(qJ(6));
t543 = cos(qJ(6));
t409 = t411 * t543 - t412 * t537;
t457 = t475 * t543 - t476 * t537;
t423 = qJD(6) * t457 + t442 * t537 + t443 * t543;
t458 = t475 * t537 + t476 * t543;
t438 = -mrSges(7,1) * t457 + mrSges(7,2) * t458;
t510 = qJD(6) + t512;
t451 = -mrSges(7,2) * t510 + mrSges(7,3) * t457;
t486 = qJDD(6) + t487;
t405 = m(7) * t409 + mrSges(7,1) * t486 - mrSges(7,3) * t423 - t438 * t458 + t451 * t510;
t410 = t411 * t537 + t412 * t543;
t422 = -qJD(6) * t458 + t442 * t543 - t443 * t537;
t452 = mrSges(7,1) * t510 - mrSges(7,3) * t458;
t406 = m(7) * t410 - mrSges(7,2) * t486 + mrSges(7,3) * t422 + t438 * t457 - t452 * t510;
t397 = t543 * t405 + t537 * t406;
t459 = -mrSges(6,1) * t475 + mrSges(6,2) * t476;
t466 = -mrSges(6,2) * t512 + mrSges(6,3) * t475;
t394 = m(6) * t414 + mrSges(6,1) * t487 - mrSges(6,3) * t443 - t459 * t476 + t466 * t512 + t397;
t467 = mrSges(6,1) * t512 - mrSges(6,3) * t476;
t560 = -t405 * t537 + t543 * t406;
t395 = m(6) * t415 - mrSges(6,2) * t487 + mrSges(6,3) * t442 + t459 * t475 - t467 * t512 + t560;
t390 = t544 * t394 + t538 * t395;
t480 = -mrSges(5,1) * t504 + mrSges(5,2) * t505;
t495 = -mrSges(5,2) * t514 + mrSges(5,3) * t504;
t388 = m(5) * t430 + mrSges(5,1) * t491 - mrSges(5,3) * t465 - t480 * t505 + t495 * t514 + t390;
t496 = mrSges(5,1) * t514 - mrSges(5,3) * t505;
t561 = -t394 * t538 + t544 * t395;
t389 = m(5) * t431 - mrSges(5,2) * t491 + mrSges(5,3) * t464 + t480 * t504 - t496 * t514 + t561;
t562 = -t388 * t539 + t545 * t389;
t380 = m(4) * t462 - mrSges(4,2) * t534 + mrSges(4,3) * t492 + t501 * t518 - t509 * t535 + t562;
t461 = t482 * t546 - t540 * t484;
t508 = -mrSges(4,2) * t535 + mrSges(4,3) * t518;
t449 = -pkin(3) * t534 - pkin(9) * t533 + t519 * t502 - t461;
t432 = -pkin(4) * t464 - pkin(10) * t503 + t505 * t497 + t449;
t417 = -pkin(5) * t442 - pkin(11) * t474 + t468 * t476 + t432;
t558 = m(7) * t417 - t422 * mrSges(7,1) + t423 * mrSges(7,2) - t457 * t451 + t458 * t452;
t553 = m(6) * t432 - t442 * mrSges(6,1) + mrSges(6,2) * t443 - t475 * t466 + t467 * t476 + t558;
t551 = -m(5) * t449 + t464 * mrSges(5,1) - mrSges(5,2) * t465 + t504 * t495 - t496 * t505 - t553;
t401 = m(4) * t461 + mrSges(4,1) * t534 - mrSges(4,3) * t493 - t501 * t519 + t508 * t535 + t551;
t377 = t540 * t380 + t546 * t401;
t382 = t545 * t388 + t539 * t389;
t563 = t546 * t380 - t401 * t540;
t434 = Ifges(7,4) * t458 + Ifges(7,2) * t457 + Ifges(7,6) * t510;
t435 = Ifges(7,1) * t458 + Ifges(7,4) * t457 + Ifges(7,5) * t510;
t556 = -mrSges(7,1) * t409 + mrSges(7,2) * t410 - Ifges(7,5) * t423 - Ifges(7,6) * t422 - Ifges(7,3) * t486 - t458 * t434 + t457 * t435;
t433 = Ifges(7,5) * t458 + Ifges(7,6) * t457 + Ifges(7,3) * t510;
t398 = -mrSges(7,1) * t417 + mrSges(7,3) * t410 + Ifges(7,4) * t423 + Ifges(7,2) * t422 + Ifges(7,6) * t486 - t433 * t458 + t435 * t510;
t399 = mrSges(7,2) * t417 - mrSges(7,3) * t409 + Ifges(7,1) * t423 + Ifges(7,4) * t422 + Ifges(7,5) * t486 + t433 * t457 - t434 * t510;
t453 = Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t512;
t455 = Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t512;
t383 = -mrSges(6,1) * t432 + mrSges(6,3) * t415 + Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t487 - pkin(5) * t558 + pkin(11) * t560 + t543 * t398 + t537 * t399 - t476 * t453 + t512 * t455;
t454 = Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t512;
t384 = mrSges(6,2) * t432 - mrSges(6,3) * t414 + Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t487 - pkin(11) * t397 - t398 * t537 + t399 * t543 + t453 * t475 - t454 * t512;
t469 = Ifges(5,5) * t505 + Ifges(5,6) * t504 + Ifges(5,3) * t514;
t471 = Ifges(5,1) * t505 + Ifges(5,4) * t504 + Ifges(5,5) * t514;
t374 = -mrSges(5,1) * t449 + mrSges(5,3) * t431 + Ifges(5,4) * t465 + Ifges(5,2) * t464 + Ifges(5,6) * t491 - pkin(4) * t553 + pkin(10) * t561 + t544 * t383 + t538 * t384 - t505 * t469 + t514 * t471;
t470 = Ifges(5,4) * t505 + Ifges(5,2) * t504 + Ifges(5,6) * t514;
t376 = mrSges(5,2) * t449 - mrSges(5,3) * t430 + Ifges(5,1) * t465 + Ifges(5,4) * t464 + Ifges(5,5) * t491 - pkin(10) * t390 - t383 * t538 + t384 * t544 + t469 * t504 - t470 * t514;
t499 = Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t535;
t500 = Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t535;
t555 = mrSges(4,1) * t461 - mrSges(4,2) * t462 + Ifges(4,5) * t493 + Ifges(4,6) * t492 + Ifges(4,3) * t534 + pkin(3) * t551 + pkin(9) * t562 + t545 * t374 + t539 * t376 + t519 * t499 - t500 * t518;
t554 = m(4) * t494 - mrSges(4,1) * t492 + mrSges(4,2) * t493 - t508 * t518 + t509 * t519 + t382;
t552 = -mrSges(6,1) * t414 + mrSges(6,2) * t415 - Ifges(6,5) * t443 - Ifges(6,6) * t442 - Ifges(6,3) * t487 - pkin(5) * t397 - t476 * t454 + t475 * t455 + t556;
t550 = mrSges(5,1) * t430 - mrSges(5,2) * t431 + Ifges(5,5) * t465 + Ifges(5,6) * t464 + Ifges(5,3) * t491 + pkin(4) * t390 + t505 * t470 - t504 * t471 - t552;
t528 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t566;
t527 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t567;
t524 = (-t547 * mrSges(3,1) + t541 * mrSges(3,2)) * qJD(1);
t520 = -pkin(7) * t549 + t557;
t517 = Ifges(3,5) * qJD(2) + (t541 * Ifges(3,1) + t547 * Ifges(3,4)) * qJD(1);
t516 = Ifges(3,6) * qJD(2) + (t541 * Ifges(3,4) + t547 * Ifges(3,2)) * qJD(1);
t506 = -g(3) * t547 - t568;
t498 = Ifges(4,5) * t519 + Ifges(4,6) * t518 + Ifges(4,3) * t535;
t372 = -mrSges(4,1) * t494 + mrSges(4,3) * t462 + Ifges(4,4) * t493 + Ifges(4,2) * t492 + Ifges(4,6) * t534 - pkin(3) * t382 - t519 * t498 + t535 * t500 - t550;
t371 = mrSges(4,2) * t494 - mrSges(4,3) * t461 + Ifges(4,1) * t493 + Ifges(4,4) * t492 + Ifges(4,5) * t534 - pkin(9) * t382 - t374 * t539 + t376 * t545 + t498 * t518 - t499 * t535;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t564 - mrSges(2,2) * t559 + t541 * (mrSges(3,2) * t520 - mrSges(3,3) * t506 + Ifges(3,1) * t525 + Ifges(3,4) * t526 + Ifges(3,5) * qJDD(2) - pkin(8) * t377 - qJD(2) * t516 + t546 * t371 - t540 * t372) + t547 * (-mrSges(3,1) * t520 + mrSges(3,3) * t507 + Ifges(3,4) * t525 + Ifges(3,2) * t526 + Ifges(3,6) * qJDD(2) - pkin(2) * t554 + pkin(8) * t563 + qJD(2) * t517 + t540 * t371 + t546 * t372) + pkin(1) * (-m(3) * t520 + mrSges(3,1) * t526 - mrSges(3,2) * t525 + (-t527 * t541 + t528 * t547) * qJD(1) - t554) + pkin(7) * (t547 * (m(3) * t507 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t526 - qJD(2) * t527 + t524 * t566 + t563) - t541 * (m(3) * t506 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t525 + qJD(2) * t528 - t524 * t567 + t377)); Ifges(3,3) * qJDD(2) + t555 + mrSges(3,1) * t506 - mrSges(3,2) * t507 + Ifges(3,5) * t525 + Ifges(3,6) * t526 + pkin(2) * t377 + (t541 * t516 - t547 * t517) * qJD(1); t555; t550; -t552; -t556;];
tauJ  = t1;
