% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:21:02
% EndTime: 2019-07-18 17:21:09
% DurationCPUTime: 2.34s
% Computational Cost: add. (16235->280), mult. (34981->341), div. (0->0), fcn. (21098->8), ass. (0->111)
t562 = Ifges(3,1) + Ifges(4,1);
t554 = Ifges(3,4) + Ifges(4,4);
t553 = Ifges(3,5) + Ifges(4,5);
t561 = -Ifges(3,2) - Ifges(4,2);
t560 = Ifges(3,6) + Ifges(4,6);
t559 = Ifges(3,3) + Ifges(4,3);
t520 = sin(qJ(4));
t521 = sin(qJ(2));
t524 = cos(qJ(4));
t525 = cos(qJ(2));
t489 = (t520 * t521 - t524 * t525) * qJD(1);
t558 = -pkin(1) - pkin(2);
t557 = pkin(1) * t525;
t556 = t525 * g(3);
t555 = -mrSges(3,2) - mrSges(4,2);
t551 = -pkin(3) - qJ(3);
t527 = qJD(1) ^ 2;
t550 = t525 ^ 2 * t527;
t522 = sin(qJ(1));
t526 = cos(qJ(1));
t508 = -t526 * g(1) - t522 * g(2);
t549 = t521 * t508;
t548 = t521 * t527;
t540 = qJD(1) * qJD(2);
t538 = t525 * t540;
t498 = t521 * qJDD(1) + t538;
t539 = qJD(1) * qJD(3);
t530 = qJDD(2) * pkin(1) + qJ(3) * t538 - 0.2e1 * t521 * t539 + t548 * t557 - t549;
t451 = qJDD(2) * pkin(2) + t551 * t498 + (pkin(2) * t548 + pkin(3) * t540 - g(3)) * t525 + t530;
t499 = t525 * qJDD(1) - t521 * t540;
t542 = qJD(1) * t521;
t501 = qJD(2) * pkin(1) - qJ(3) * t542;
t506 = qJD(2) * pkin(2) - pkin(3) * t542;
t480 = -t521 * g(3) + t525 * t508;
t533 = t499 * qJ(3) + 0.2e1 * t525 * t539 + t480;
t452 = t499 * pkin(3) + t558 * t550 + (-t501 - t506) * qJD(2) + t533;
t445 = t520 * t451 + t524 * t452;
t490 = (t520 * t525 + t521 * t524) * qJD(1);
t465 = -t490 * qJD(4) - t520 * t498 + t524 * t499;
t474 = t489 * mrSges(5,1) + t490 * mrSges(5,2);
t515 = qJD(2) + qJD(4);
t478 = t515 * mrSges(5,1) - t490 * mrSges(5,3);
t514 = qJDD(2) + qJDD(4);
t442 = (t489 * t490 + t514) * pkin(4) + t445;
t507 = t522 * g(1) - t526 * g(2);
t532 = t501 * t542 + qJDD(3) - t507;
t454 = t499 * t558 + t506 * t542 + t551 * t550 + t532;
t466 = -t489 * qJD(4) + t524 * t498 + t520 * t499;
t446 = (t489 * t515 - t466) * pkin(4) + t454;
t519 = sin(qJ(5));
t523 = cos(qJ(5));
t440 = -t519 * t442 + t523 * t446;
t475 = -t519 * t490 + t523 * t515;
t450 = t475 * qJD(5) + t523 * t466 + t519 * t514;
t476 = t523 * t490 + t519 * t515;
t459 = -t475 * mrSges(6,1) + t476 * mrSges(6,2);
t464 = qJDD(5) - t465;
t482 = qJD(5) + t489;
t467 = -t482 * mrSges(6,2) + t475 * mrSges(6,3);
t438 = m(6) * t440 + t464 * mrSges(6,1) - t450 * mrSges(6,3) - t476 * t459 + t482 * t467;
t441 = t523 * t442 + t519 * t446;
t449 = -t476 * qJD(5) - t519 * t466 + t523 * t514;
t468 = t482 * mrSges(6,1) - t476 * mrSges(6,3);
t439 = m(6) * t441 - t464 * mrSges(6,2) + t449 * mrSges(6,3) + t475 * t459 - t482 * t468;
t536 = -t519 * t438 + t523 * t439;
t429 = m(5) * t445 - t514 * mrSges(5,2) + t465 * mrSges(5,3) - t489 * t474 - t515 * t478 + t536;
t444 = t524 * t451 - t520 * t452;
t443 = (-t490 ^ 2 - t515 ^ 2) * pkin(4) - t444;
t477 = -t515 * mrSges(5,2) - t489 * mrSges(5,3);
t434 = m(5) * t444 - m(6) * t443 + t514 * mrSges(5,1) + t449 * mrSges(6,1) - t450 * mrSges(6,2) - t466 * mrSges(5,3) + t475 * t467 - t476 * t468 - t490 * t474 + t515 * t477;
t425 = t520 * t429 + t524 * t434;
t547 = t523 * t438 + t519 * t439;
t546 = t559 * qJD(2) + (t521 * t553 + t525 * t560) * qJD(1);
t545 = -t560 * qJD(2) + (-t521 * t554 + t525 * t561) * qJD(1);
t544 = t553 * qJD(2) + (t521 * t562 + t554 * t525) * qJD(1);
t502 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t542;
t543 = -qJD(2) * mrSges(3,1) + mrSges(3,3) * t542 - t502;
t541 = qJD(1) * t525;
t537 = t524 * t429 - t520 * t434;
t461 = -t498 * qJ(3) + t530 - t556;
t504 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t541;
t535 = m(4) * t461 + qJDD(2) * mrSges(4,1) + qJD(2) * t504 + t425;
t462 = -pkin(1) * t550 - qJD(2) * t501 + t533;
t496 = (-mrSges(4,1) * t525 + mrSges(4,2) * t521) * qJD(1);
t531 = m(4) * t462 + t499 * mrSges(4,3) + t496 * t541 + t537;
t529 = -m(5) * t454 + t465 * mrSges(5,1) - t466 * mrSges(5,2) - t489 * t477 - t490 * t478 - t547;
t469 = -t499 * pkin(1) - qJ(3) * t550 + t532;
t528 = m(4) * t469 - t529;
t505 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t541;
t497 = (-mrSges(3,1) * t525 + mrSges(3,2) * t521) * qJD(1);
t479 = -t549 - t556;
t472 = Ifges(5,1) * t490 - Ifges(5,4) * t489 + Ifges(5,5) * t515;
t471 = Ifges(5,4) * t490 - Ifges(5,2) * t489 + Ifges(5,6) * t515;
t470 = Ifges(5,5) * t490 - Ifges(5,6) * t489 + Ifges(5,3) * t515;
t457 = Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t482;
t456 = Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t482;
t455 = Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t482;
t432 = mrSges(6,2) * t443 - mrSges(6,3) * t440 + Ifges(6,1) * t450 + Ifges(6,4) * t449 + Ifges(6,5) * t464 + t475 * t455 - t482 * t456;
t431 = -mrSges(6,1) * t443 + mrSges(6,3) * t441 + Ifges(6,4) * t450 + Ifges(6,2) * t449 + Ifges(6,6) * t464 - t476 * t455 + t482 * t457;
t430 = -mrSges(5,1) * t454 - mrSges(6,1) * t440 + mrSges(6,2) * t441 + mrSges(5,3) * t445 + Ifges(5,4) * t466 - Ifges(6,5) * t450 + Ifges(5,2) * t465 + Ifges(5,6) * t514 - Ifges(6,6) * t449 - Ifges(6,3) * t464 - t476 * t456 + t475 * t457 - t490 * t470 + t515 * t472;
t426 = qJDD(1) * mrSges(2,1) - t527 * mrSges(2,2) + (m(2) + m(3)) * t507 + (mrSges(3,1) + mrSges(4,1)) * t499 + t555 * t498 + ((t504 + t505) * t525 + t543 * t521) * qJD(1) - t528;
t424 = -t498 * mrSges(4,3) - t496 * t542 + t535;
t423 = m(3) * t480 + t499 * mrSges(3,3) + qJD(2) * t543 + qJDD(2) * t555 + t497 * t541 + t531;
t422 = m(3) * t479 + qJDD(2) * mrSges(3,1) + qJD(2) * t505 + (-mrSges(3,3) - mrSges(4,3)) * t498 + (-t496 - t497) * t542 + t535;
t421 = mrSges(5,2) * t454 - mrSges(5,3) * t444 + Ifges(5,1) * t466 + Ifges(5,4) * t465 + Ifges(5,5) * t514 - pkin(4) * t547 - t519 * t431 + t523 * t432 - t489 * t470 - t515 * t471;
t420 = m(2) * t508 - t527 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t521 * t422 + t525 * t423;
t419 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + t527 * Ifges(2,5) - pkin(4) * t536 - t519 * t432 - t523 * t431 + mrSges(2,3) * t508 - Ifges(5,3) * t514 - mrSges(3,1) * t479 + mrSges(3,2) * t480 - t489 * t472 - t490 * t471 - mrSges(4,1) * t461 + mrSges(4,2) * t462 - Ifges(5,6) * t465 - Ifges(5,5) * t466 - mrSges(5,1) * t444 + mrSges(5,2) * t445 - pkin(2) * t425 - pkin(1) * t424 - t560 * t499 - t553 * t498 - t559 * qJDD(2) + (t521 * t545 + t525 * t544) * qJD(1);
t418 = -mrSges(3,2) * t507 + mrSges(4,2) * t469 - mrSges(3,3) * t479 - mrSges(4,3) * t461 - pkin(3) * t425 - qJ(3) * t424 + t545 * qJD(2) + t553 * qJDD(2) + t524 * t421 - t520 * t430 + t498 * t562 + t554 * t499 + t546 * t541;
t417 = mrSges(3,1) * t507 + mrSges(3,3) * t480 - mrSges(4,1) * t469 + mrSges(4,3) * t462 + t520 * t421 + t524 * t430 + pkin(2) * t529 + pkin(3) * t537 - pkin(1) * t528 + qJ(3) * t531 + (pkin(1) * mrSges(4,1) - t561) * t499 + (-pkin(1) * mrSges(4,2) + t554) * t498 + (-qJ(3) * mrSges(4,2) + t560) * qJDD(2) + (-qJ(3) * t502 + t544) * qJD(2) + (t504 * t557 + (-pkin(1) * t502 - t546) * t521) * qJD(1);
t416 = -mrSges(2,2) * g(3) - mrSges(2,3) * t507 + Ifges(2,5) * qJDD(1) - t527 * Ifges(2,6) - t521 * t417 + t525 * t418;
t1 = [-m(1) * g(1) + t526 * t420 - t522 * t426; -m(1) * g(2) + t522 * t420 + t526 * t426; t525 * t422 + t521 * t423 + (-m(1) - m(2)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t526 * t416 - t522 * t419; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t522 * t416 + t526 * t419; -mrSges(1,1) * g(2) + mrSges(2,1) * t507 + mrSges(1,2) * g(1) - mrSges(2,2) * t508 + Ifges(2,3) * qJDD(1) + t525 * t417 + t521 * t418;];
tauB  = t1;
