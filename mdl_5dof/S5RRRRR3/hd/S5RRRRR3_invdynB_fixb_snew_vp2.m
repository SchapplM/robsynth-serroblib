% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:37
% EndTime: 2019-07-18 17:17:45
% DurationCPUTime: 4.44s
% Computational Cost: add. (50668->290), mult. (102774->372), div. (0->0), fcn. (77087->10), ass. (0->110)
t510 = sin(qJ(1));
t515 = cos(qJ(1));
t499 = -t515 * g(1) - t510 * g(2);
t509 = sin(qJ(2));
t514 = cos(qJ(2));
t483 = -t514 * g(3) - t509 * t499;
t493 = (-mrSges(3,1) * t514 + mrSges(3,2) * t509) * qJD(1);
t523 = qJD(1) * qJD(2);
t494 = t509 * qJDD(1) + t514 * t523;
t524 = qJD(1) * t514;
t497 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t524;
t525 = qJD(1) * t509;
t516 = qJD(1) ^ 2;
t475 = (t509 * t514 * t516 + qJDD(2)) * pkin(1) + t483;
t484 = -t509 * g(3) + t514 * t499;
t478 = (-t514 ^ 2 * t516 - qJD(2) ^ 2) * pkin(1) + t484;
t508 = sin(qJ(3));
t513 = cos(qJ(3));
t455 = t508 * t475 + t513 * t478;
t490 = (t508 * t514 + t509 * t513) * qJD(1);
t522 = t509 * t523;
t495 = t514 * qJDD(1) - t522;
t462 = -t490 * qJD(3) - t508 * t494 + t513 * t495;
t489 = -t508 * t525 + t513 * t524;
t469 = -t489 * mrSges(4,1) + t490 * mrSges(4,2);
t505 = qJD(2) + qJD(3);
t480 = t505 * mrSges(4,1) - t490 * mrSges(4,3);
t504 = qJDD(2) + qJDD(3);
t463 = t489 * qJD(3) + t513 * t494 + t508 * t495;
t498 = t510 * g(1) - t515 * g(2);
t472 = -t498 + (-t495 + t522) * pkin(1);
t436 = (-t489 * t505 - t463) * pkin(5) + (t490 * t505 - t462) * pkin(2) + t472;
t470 = -t489 * pkin(2) - t490 * pkin(5);
t503 = t505 ^ 2;
t444 = -t503 * pkin(2) + t504 * pkin(5) + t489 * t470 + t455;
t507 = sin(qJ(4));
t512 = cos(qJ(4));
t428 = t512 * t436 - t507 * t444;
t476 = -t507 * t490 + t512 * t505;
t447 = t476 * qJD(4) + t512 * t463 + t507 * t504;
t477 = t512 * t490 + t507 * t505;
t458 = -t476 * mrSges(5,1) + t477 * mrSges(5,2);
t461 = qJDD(4) - t462;
t485 = qJD(4) - t489;
t464 = -t485 * mrSges(5,2) + t476 * mrSges(5,3);
t426 = (t476 * t477 + t461) * pkin(3) + t428;
t429 = t507 * t436 + t512 * t444;
t427 = (-t476 ^ 2 - t485 ^ 2) * pkin(3) + t429;
t506 = sin(qJ(5));
t511 = cos(qJ(5));
t424 = t511 * t426 - t506 * t427;
t446 = -t477 * qJD(4) - t507 * t463 + t512 * t504;
t456 = t511 * t476 - t506 * t477;
t433 = t456 * qJD(5) + t506 * t446 + t511 * t447;
t457 = t506 * t476 + t511 * t477;
t442 = -t456 * mrSges(6,1) + t457 * mrSges(6,2);
t481 = qJD(5) + t485;
t448 = -t481 * mrSges(6,2) + t456 * mrSges(6,3);
t459 = qJDD(5) + t461;
t422 = m(6) * t424 + t459 * mrSges(6,1) - t433 * mrSges(6,3) - t457 * t442 + t481 * t448;
t425 = t506 * t426 + t511 * t427;
t432 = -t457 * qJD(5) + t511 * t446 - t506 * t447;
t449 = t481 * mrSges(6,1) - t457 * mrSges(6,3);
t423 = m(6) * t425 - t459 * mrSges(6,2) + t432 * mrSges(6,3) + t456 * t442 - t481 * t449;
t526 = t511 * t422 + t506 * t423;
t414 = m(5) * t428 + t461 * mrSges(5,1) - t447 * mrSges(5,3) - t477 * t458 + t485 * t464 + t526;
t465 = t485 * mrSges(5,1) - t477 * mrSges(5,3);
t415 = m(5) * t429 - t461 * mrSges(5,2) + t446 * mrSges(5,3) - t506 * t422 + t511 * t423 + t476 * t458 - t485 * t465;
t520 = -t507 * t414 + t512 * t415;
t407 = m(4) * t455 - t504 * mrSges(4,2) + t462 * mrSges(4,3) + t489 * t469 - t505 * t480 + t520;
t454 = t513 * t475 - t508 * t478;
t479 = -t505 * mrSges(4,2) + t489 * mrSges(4,3);
t443 = -t504 * pkin(2) - t503 * pkin(5) + t490 * t470 - t454;
t430 = (t477 * t485 - t446) * pkin(3) + t443;
t519 = m(6) * t430 - t432 * mrSges(6,1) + t433 * mrSges(6,2) - t456 * t448 + t457 * t449;
t517 = -m(5) * t443 + t446 * mrSges(5,1) - t447 * mrSges(5,2) + t476 * t464 - t477 * t465 - t519;
t419 = m(4) * t454 + t504 * mrSges(4,1) - t463 * mrSges(4,3) - t490 * t469 + t505 * t479 + t517;
t527 = t508 * t407 + t513 * t419;
t402 = m(3) * t483 + qJDD(2) * mrSges(3,1) - t494 * mrSges(3,3) + qJD(2) * t497 - t493 * t525 + t527;
t496 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t525;
t403 = m(3) * t484 - qJDD(2) * mrSges(3,2) + t495 * mrSges(3,3) - qJD(2) * t496 + t513 * t407 - t508 * t419 + t493 * t524;
t400 = m(2) * t499 - t516 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t509 * t402 + t514 * t403;
t408 = t512 * t414 + t507 * t415;
t518 = m(4) * t472 - t462 * mrSges(4,1) + t463 * mrSges(4,2) - t489 * t479 + t490 * t480 + t408;
t405 = qJDD(1) * mrSges(2,1) + t495 * mrSges(3,1) - t516 * mrSges(2,2) - t494 * mrSges(3,2) + (m(2) + m(3)) * t498 + (-t496 * t509 + t497 * t514) * qJD(1) - t518;
t528 = t510 * t400 + t515 * t405;
t521 = t515 * t400 - t510 * t405;
t488 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t509 + Ifges(3,4) * t514) * qJD(1);
t487 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t509 + Ifges(3,2) * t514) * qJD(1);
t486 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t509 + Ifges(3,6) * t514) * qJD(1);
t468 = Ifges(4,1) * t490 + Ifges(4,4) * t489 + Ifges(4,5) * t505;
t467 = Ifges(4,4) * t490 + Ifges(4,2) * t489 + Ifges(4,6) * t505;
t466 = Ifges(4,5) * t490 + Ifges(4,6) * t489 + Ifges(4,3) * t505;
t452 = Ifges(5,1) * t477 + Ifges(5,4) * t476 + Ifges(5,5) * t485;
t451 = Ifges(5,4) * t477 + Ifges(5,2) * t476 + Ifges(5,6) * t485;
t450 = Ifges(5,5) * t477 + Ifges(5,6) * t476 + Ifges(5,3) * t485;
t439 = Ifges(6,1) * t457 + Ifges(6,4) * t456 + Ifges(6,5) * t481;
t438 = Ifges(6,4) * t457 + Ifges(6,2) * t456 + Ifges(6,6) * t481;
t437 = Ifges(6,5) * t457 + Ifges(6,6) * t456 + Ifges(6,3) * t481;
t417 = mrSges(6,2) * t430 - mrSges(6,3) * t424 + Ifges(6,1) * t433 + Ifges(6,4) * t432 + Ifges(6,5) * t459 + t456 * t437 - t481 * t438;
t416 = -mrSges(6,1) * t430 + mrSges(6,3) * t425 + Ifges(6,4) * t433 + Ifges(6,2) * t432 + Ifges(6,6) * t459 - t457 * t437 + t481 * t439;
t410 = mrSges(5,2) * t443 - mrSges(5,3) * t428 + Ifges(5,1) * t447 + Ifges(5,4) * t446 + Ifges(5,5) * t461 - t506 * t416 + t511 * t417 + t476 * t450 - t485 * t451;
t409 = -mrSges(5,1) * t443 + mrSges(5,3) * t429 + Ifges(5,4) * t447 + Ifges(5,2) * t446 + Ifges(5,6) * t461 - pkin(3) * t519 + t511 * t416 + t506 * t417 - t477 * t450 + t485 * t452;
t401 = Ifges(4,4) * t463 + Ifges(4,2) * t462 + Ifges(4,6) * t504 - t490 * t466 + t505 * t468 - mrSges(4,1) * t472 + mrSges(4,3) * t455 - Ifges(5,5) * t447 - Ifges(5,6) * t446 - Ifges(5,3) * t461 - t477 * t451 + t476 * t452 - mrSges(5,1) * t428 + mrSges(5,2) * t429 - Ifges(6,5) * t433 - Ifges(6,6) * t432 - Ifges(6,3) * t459 - t457 * t438 + t456 * t439 - mrSges(6,1) * t424 + mrSges(6,2) * t425 - pkin(3) * t526 - pkin(2) * t408;
t397 = mrSges(4,2) * t472 - mrSges(4,3) * t454 + Ifges(4,1) * t463 + Ifges(4,4) * t462 + Ifges(4,5) * t504 - pkin(5) * t408 - t507 * t409 + t512 * t410 + t489 * t466 - t505 * t467;
t396 = t516 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t499 - pkin(1) * t527 - Ifges(3,5) * t494 - Ifges(3,6) * t495 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t483 + mrSges(3,2) * t484 - t512 * t409 - pkin(2) * t517 - pkin(5) * t520 - t507 * t410 - Ifges(4,5) * t463 - Ifges(4,6) * t462 - Ifges(4,3) * t504 - t490 * t467 + t489 * t468 - mrSges(4,1) * t454 + mrSges(4,2) * t455 + Ifges(2,6) * qJDD(1) + (-t509 * t487 + t514 * t488) * qJD(1);
t395 = -mrSges(3,2) * t498 - mrSges(3,3) * t483 + Ifges(3,1) * t494 + Ifges(3,4) * t495 + Ifges(3,5) * qJDD(2) - qJD(2) * t487 + t513 * t397 - t508 * t401 + t486 * t524;
t394 = mrSges(3,1) * t498 + mrSges(3,3) * t484 + Ifges(3,4) * t494 + Ifges(3,2) * t495 + Ifges(3,6) * qJDD(2) - pkin(1) * t518 + qJD(2) * t488 + t508 * t397 + t513 * t401 - t486 * t525;
t393 = -mrSges(2,2) * g(3) - mrSges(2,3) * t498 + Ifges(2,5) * qJDD(1) - t516 * Ifges(2,6) - t509 * t394 + t514 * t395;
t1 = [-m(1) * g(1) + t521; -m(1) * g(2) + t528; t514 * t402 + t509 * t403 + (-m(1) - m(2)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t528 + t515 * t393 - t510 * t396; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t521 + t510 * t393 + t515 * t396; -mrSges(1,1) * g(2) + mrSges(2,1) * t498 + mrSges(1,2) * g(1) - mrSges(2,2) * t499 + Ifges(2,3) * qJDD(1) + t514 * t394 + t509 * t395;];
tauB  = t1;
