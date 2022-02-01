% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:07
% EndTime: 2022-01-20 12:08:14
% DurationCPUTime: 5.56s
% Computational Cost: add. (93727->273), mult. (120226->346), div. (0->0), fcn. (77730->10), ass. (0->110)
t515 = qJD(1) + qJD(2);
t511 = t515 ^ 2;
t542 = pkin(3) * t511;
t519 = sin(qJ(3));
t541 = t515 * t519;
t524 = cos(qJ(3));
t540 = t515 * t524;
t521 = sin(qJ(1));
t526 = cos(qJ(1));
t506 = t521 * g(1) - t526 * g(2);
t501 = qJDD(1) * pkin(1) + t506;
t507 = -t526 * g(1) - t521 * g(2);
t527 = qJD(1) ^ 2;
t502 = -t527 * pkin(1) + t507;
t520 = sin(qJ(2));
t525 = cos(qJ(2));
t482 = t520 * t501 + t525 * t502;
t513 = qJDD(1) + qJDD(2);
t480 = -t511 * pkin(2) + t513 * pkin(7) + t482;
t539 = t519 * t480;
t537 = qJD(3) * t515;
t496 = t519 * t513 + t524 * t537;
t461 = qJDD(3) * pkin(3) - t496 * pkin(8) - t539 + (pkin(8) * t537 + t519 * t542 - g(3)) * t524;
t470 = -t519 * g(3) + t524 * t480;
t497 = t524 * t513 - t519 * t537;
t505 = qJD(3) * pkin(3) - pkin(8) * t541;
t516 = t524 ^ 2;
t462 = t497 * pkin(8) - qJD(3) * t505 - t516 * t542 + t470;
t518 = sin(qJ(4));
t523 = cos(qJ(4));
t446 = t523 * t461 - t518 * t462;
t490 = (-t518 * t519 + t523 * t524) * t515;
t466 = t490 * qJD(4) + t523 * t496 + t518 * t497;
t491 = (t518 * t524 + t519 * t523) * t515;
t512 = qJDD(3) + qJDD(4);
t514 = qJD(3) + qJD(4);
t442 = (t490 * t514 - t466) * pkin(9) + (t490 * t491 + t512) * pkin(4) + t446;
t447 = t518 * t461 + t523 * t462;
t465 = -t491 * qJD(4) - t518 * t496 + t523 * t497;
t485 = t514 * pkin(4) - t491 * pkin(9);
t486 = t490 ^ 2;
t443 = -t486 * pkin(4) + t465 * pkin(9) - t514 * t485 + t447;
t517 = sin(qJ(5));
t522 = cos(qJ(5));
t440 = t522 * t442 - t517 * t443;
t475 = t522 * t490 - t517 * t491;
t451 = t475 * qJD(5) + t517 * t465 + t522 * t466;
t476 = t517 * t490 + t522 * t491;
t457 = -t475 * mrSges(6,1) + t476 * mrSges(6,2);
t509 = qJD(5) + t514;
t467 = -t509 * mrSges(6,2) + t475 * mrSges(6,3);
t508 = qJDD(5) + t512;
t438 = m(6) * t440 + t508 * mrSges(6,1) - t451 * mrSges(6,3) - t476 * t457 + t509 * t467;
t441 = t517 * t442 + t522 * t443;
t450 = -t476 * qJD(5) + t522 * t465 - t517 * t466;
t468 = t509 * mrSges(6,1) - t476 * mrSges(6,3);
t439 = m(6) * t441 - t508 * mrSges(6,2) + t450 * mrSges(6,3) + t475 * t457 - t509 * t468;
t430 = t522 * t438 + t517 * t439;
t478 = -t490 * mrSges(5,1) + t491 * mrSges(5,2);
t483 = -t514 * mrSges(5,2) + t490 * mrSges(5,3);
t428 = m(5) * t446 + t512 * mrSges(5,1) - t466 * mrSges(5,3) - t491 * t478 + t514 * t483 + t430;
t484 = t514 * mrSges(5,1) - t491 * mrSges(5,3);
t532 = -t517 * t438 + t522 * t439;
t429 = m(5) * t447 - t512 * mrSges(5,2) + t465 * mrSges(5,3) + t490 * t478 - t514 * t484 + t532;
t424 = t523 * t428 + t518 * t429;
t469 = -t524 * g(3) - t539;
t495 = (-mrSges(4,1) * t524 + mrSges(4,2) * t519) * t515;
t504 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t540;
t422 = m(4) * t469 + qJDD(3) * mrSges(4,1) - t496 * mrSges(4,3) + qJD(3) * t504 - t495 * t541 + t424;
t503 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t541;
t533 = -t518 * t428 + t523 * t429;
t423 = m(4) * t470 - qJDD(3) * mrSges(4,2) + t497 * mrSges(4,3) - qJD(3) * t503 + t495 * t540 + t533;
t534 = -t519 * t422 + t524 * t423;
t415 = m(3) * t482 - t511 * mrSges(3,1) - t513 * mrSges(3,2) + t534;
t481 = t525 * t501 - t520 * t502;
t530 = -t513 * pkin(2) - t481;
t479 = -t511 * pkin(7) + t530;
t463 = -t497 * pkin(3) + t505 * t541 + (-pkin(8) * t516 - pkin(7)) * t511 + t530;
t445 = -t465 * pkin(4) - t486 * pkin(9) + t491 * t485 + t463;
t531 = m(6) * t445 - t450 * mrSges(6,1) + t451 * mrSges(6,2) - t475 * t467 + t476 * t468;
t529 = m(5) * t463 - t465 * mrSges(5,1) + t466 * mrSges(5,2) - t490 * t483 + t491 * t484 + t531;
t528 = -m(4) * t479 + t497 * mrSges(4,1) - t496 * mrSges(4,2) - t503 * t541 + t504 * t540 - t529;
t434 = m(3) * t481 + t513 * mrSges(3,1) - t511 * mrSges(3,2) + t528;
t412 = t520 * t415 + t525 * t434;
t410 = m(2) * t506 + qJDD(1) * mrSges(2,1) - t527 * mrSges(2,2) + t412;
t535 = t525 * t415 - t520 * t434;
t411 = m(2) * t507 - t527 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t535;
t538 = t526 * t410 + t521 * t411;
t416 = t524 * t422 + t519 * t423;
t536 = -t521 * t410 + t526 * t411;
t489 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t519 + Ifges(4,4) * t524) * t515;
t488 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t519 + Ifges(4,2) * t524) * t515;
t487 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t519 + Ifges(4,6) * t524) * t515;
t473 = Ifges(5,1) * t491 + Ifges(5,4) * t490 + Ifges(5,5) * t514;
t472 = Ifges(5,4) * t491 + Ifges(5,2) * t490 + Ifges(5,6) * t514;
t471 = Ifges(5,5) * t491 + Ifges(5,6) * t490 + Ifges(5,3) * t514;
t454 = Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t509;
t453 = Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t509;
t452 = Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t509;
t432 = mrSges(6,2) * t445 - mrSges(6,3) * t440 + Ifges(6,1) * t451 + Ifges(6,4) * t450 + Ifges(6,5) * t508 + t475 * t452 - t509 * t453;
t431 = -mrSges(6,1) * t445 + mrSges(6,3) * t441 + Ifges(6,4) * t451 + Ifges(6,2) * t450 + Ifges(6,6) * t508 - t476 * t452 + t509 * t454;
t418 = mrSges(5,2) * t463 - mrSges(5,3) * t446 + Ifges(5,1) * t466 + Ifges(5,4) * t465 + Ifges(5,5) * t512 - pkin(9) * t430 - t517 * t431 + t522 * t432 + t490 * t471 - t514 * t472;
t417 = -mrSges(5,1) * t463 + mrSges(5,3) * t447 + Ifges(5,4) * t466 + Ifges(5,2) * t465 + Ifges(5,6) * t512 - pkin(4) * t531 + pkin(9) * t532 + t522 * t431 + t517 * t432 - t491 * t471 + t514 * t473;
t406 = mrSges(4,2) * t479 - mrSges(4,3) * t469 + Ifges(4,1) * t496 + Ifges(4,4) * t497 + Ifges(4,5) * qJDD(3) - pkin(8) * t424 - qJD(3) * t488 - t518 * t417 + t523 * t418 + t487 * t540;
t405 = -mrSges(4,1) * t479 + mrSges(4,3) * t470 + Ifges(4,4) * t496 + Ifges(4,2) * t497 + Ifges(4,6) * qJDD(3) - pkin(3) * t529 + pkin(8) * t533 + qJD(3) * t489 + t523 * t417 + t518 * t418 - t487 * t541;
t404 = mrSges(3,1) * g(3) - Ifges(4,3) * qJDD(3) + (-t519 * t488 + t524 * t489) * t515 - Ifges(5,3) * t512 + Ifges(3,6) * t513 - Ifges(6,3) * t508 + t511 * Ifges(3,5) + t490 * t473 - t491 * t472 - Ifges(4,5) * t496 - Ifges(4,6) * t497 - t476 * t453 + mrSges(3,3) * t482 - Ifges(5,6) * t465 - Ifges(5,5) * t466 - mrSges(4,1) * t469 + mrSges(4,2) * t470 + t475 * t454 - mrSges(5,1) * t446 + mrSges(5,2) * t447 - Ifges(6,6) * t450 - Ifges(6,5) * t451 + mrSges(6,2) * t441 - mrSges(6,1) * t440 - pkin(4) * t430 - pkin(3) * t424 - pkin(2) * t416;
t403 = -mrSges(3,2) * g(3) - mrSges(3,3) * t481 + Ifges(3,5) * t513 - t511 * Ifges(3,6) - pkin(7) * t416 - t519 * t405 + t524 * t406;
t402 = -mrSges(2,2) * g(3) - mrSges(2,3) * t506 + Ifges(2,5) * qJDD(1) - t527 * Ifges(2,6) - pkin(6) * t412 + t525 * t403 - t520 * t404;
t401 = Ifges(2,6) * qJDD(1) + t527 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t507 + t520 * t403 + t525 * t404 - pkin(1) * (-m(3) * g(3) + t416) + pkin(6) * t535;
t1 = [-m(1) * g(1) + t536; -m(1) * g(2) + t538; (-m(1) - m(2) - m(3)) * g(3) + t416; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t538 - t521 * t401 + t526 * t402; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t536 + t526 * t401 + t521 * t402; -mrSges(1,1) * g(2) + mrSges(2,1) * t506 + mrSges(3,1) * t481 + mrSges(1,2) * g(1) - mrSges(2,2) * t507 - mrSges(3,2) * t482 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t513 + pkin(1) * t412 + pkin(2) * t528 + pkin(7) * t534 + t524 * t405 + t519 * t406;];
tauB = t1;
