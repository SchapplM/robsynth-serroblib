% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 08:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:22:39
% EndTime: 2019-05-05 08:22:43
% DurationCPUTime: 2.22s
% Computational Cost: add. (15385->270), mult. (29358->328), div. (0->0), fcn. (19973->12), ass. (0->117)
t533 = Ifges(5,1) + Ifges(6,1);
t522 = Ifges(5,4) - Ifges(6,5);
t521 = Ifges(5,5) + Ifges(6,4);
t532 = -Ifges(5,2) - Ifges(6,3);
t520 = Ifges(5,6) - Ifges(6,6);
t529 = Ifges(5,3) + Ifges(6,2);
t482 = sin(pkin(11));
t484 = cos(pkin(11));
t469 = g(1) * t482 - g(2) * t484;
t481 = -g(3) + qJDD(1);
t483 = sin(pkin(6));
t485 = cos(pkin(6));
t531 = t469 * t485 + t481 * t483;
t487 = sin(qJ(4));
t488 = sin(qJ(3));
t512 = qJD(2) * t488;
t524 = cos(qJ(4));
t463 = -qJD(3) * t524 + t487 * t512;
t491 = cos(qJ(3));
t510 = qJD(2) * qJD(3);
t508 = t491 * t510;
t467 = qJDD(2) * t488 + t508;
t431 = -t463 * qJD(4) + t487 * qJDD(3) + t467 * t524;
t464 = t487 * qJD(3) + t512 * t524;
t437 = mrSges(6,1) * t463 - mrSges(6,3) * t464;
t470 = -g(1) * t484 - g(2) * t482;
t489 = sin(qJ(2));
t492 = cos(qJ(2));
t419 = t492 * t470 + t489 * t531;
t494 = qJD(2) ^ 2;
t413 = -pkin(2) * t494 + qJDD(2) * pkin(8) + t419;
t446 = -t469 * t483 + t481 * t485;
t408 = t491 * t413 + t488 * t446;
t466 = (-pkin(3) * t491 - pkin(9) * t488) * qJD(2);
t493 = qJD(3) ^ 2;
t511 = qJD(2) * t491;
t399 = -pkin(3) * t493 + qJDD(3) * pkin(9) + t466 * t511 + t408;
t418 = -t489 * t470 + t492 * t531;
t412 = -qJDD(2) * pkin(2) - t494 * pkin(8) - t418;
t509 = t488 * t510;
t468 = t491 * qJDD(2) - t509;
t401 = (-t467 - t508) * pkin(9) + (-t468 + t509) * pkin(3) + t412;
t387 = -t487 * t399 + t524 * t401;
t436 = pkin(4) * t463 - qJ(5) * t464;
t460 = qJDD(4) - t468;
t477 = qJD(4) - t511;
t476 = t477 ^ 2;
t385 = -t460 * pkin(4) - t476 * qJ(5) + t464 * t436 + qJDD(5) - t387;
t519 = t463 * t477;
t379 = (-t431 - t519) * pkin(10) + (t463 * t464 - t460) * pkin(5) + t385;
t388 = t524 * t399 + t487 * t401;
t525 = 2 * qJD(5);
t384 = -pkin(4) * t476 + t460 * qJ(5) - t463 * t436 + t477 * t525 + t388;
t430 = qJD(4) * t464 - qJDD(3) * t524 + t467 * t487;
t445 = -pkin(5) * t477 - pkin(10) * t464;
t459 = t463 ^ 2;
t380 = -pkin(5) * t459 + pkin(10) * t430 + t445 * t477 + t384;
t486 = sin(qJ(6));
t490 = cos(qJ(6));
t377 = t379 * t490 - t380 * t486;
t432 = t463 * t490 - t464 * t486;
t395 = qJD(6) * t432 + t430 * t486 + t431 * t490;
t433 = t463 * t486 + t464 * t490;
t409 = -mrSges(7,1) * t432 + mrSges(7,2) * t433;
t475 = qJD(6) - t477;
t414 = -mrSges(7,2) * t475 + mrSges(7,3) * t432;
t456 = qJDD(6) - t460;
t374 = m(7) * t377 + mrSges(7,1) * t456 - mrSges(7,3) * t395 - t409 * t433 + t414 * t475;
t378 = t379 * t486 + t380 * t490;
t394 = -qJD(6) * t433 + t430 * t490 - t431 * t486;
t415 = mrSges(7,1) * t475 - mrSges(7,3) * t433;
t375 = m(7) * t378 - mrSges(7,2) * t456 + mrSges(7,3) * t394 + t409 * t432 - t415 * t475;
t368 = t490 * t374 + t486 * t375;
t444 = -mrSges(6,2) * t463 + mrSges(6,3) * t477;
t499 = -m(6) * t385 + t460 * mrSges(6,1) + t477 * t444 - t368;
t367 = t431 * mrSges(6,2) + t464 * t437 - t499;
t403 = Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t475;
t404 = Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t475;
t498 = -mrSges(7,1) * t377 + mrSges(7,2) * t378 - Ifges(7,5) * t395 - Ifges(7,6) * t394 - Ifges(7,3) * t456 - t433 * t403 + t432 * t404;
t443 = -mrSges(6,1) * t477 + mrSges(6,2) * t464;
t505 = -t486 * t374 + t490 * t375;
t501 = m(6) * t384 + t460 * mrSges(6,3) + t477 * t443 + t505;
t514 = -t522 * t463 + t533 * t464 + t521 * t477;
t515 = t532 * t463 + t522 * t464 + t520 * t477;
t530 = -t520 * t430 + t521 * t431 + t529 * t460 + mrSges(5,1) * t387 - mrSges(6,1) * t385 - mrSges(5,2) * t388 + mrSges(6,3) * t384 - pkin(4) * t367 - pkin(5) * t368 + qJ(5) * (-t430 * mrSges(6,2) - t463 * t437 + t501) + t498 + t515 * t464 + t514 * t463;
t407 = -t488 * t413 + t491 * t446;
t500 = qJDD(3) * pkin(3) + t493 * pkin(9) - t466 * t512 + t407;
t526 = (-t431 + t519) * qJ(5) - t500;
t523 = -mrSges(5,3) - mrSges(6,2);
t516 = t520 * t463 - t464 * t521 - t529 * t477;
t513 = -mrSges(5,1) * t463 - mrSges(5,2) * t464 - t437;
t465 = (-mrSges(4,1) * t491 + mrSges(4,2) * t488) * qJD(2);
t471 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t512;
t442 = mrSges(5,1) * t477 - mrSges(5,3) * t464;
t364 = m(5) * t388 - t460 * mrSges(5,2) + t430 * t523 - t477 * t442 + t463 * t513 + t501;
t441 = -mrSges(5,2) * t477 - mrSges(5,3) * t463;
t365 = m(5) * t387 + t460 * mrSges(5,1) + t431 * t523 + t477 * t441 + t464 * t513 + t499;
t506 = t524 * t364 - t365 * t487;
t361 = m(4) * t408 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t468 - qJD(3) * t471 + t465 * t511 + t506;
t472 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t511;
t386 = -0.2e1 * qJD(5) * t464 + (t464 * t477 + t430) * pkin(4) + t526;
t382 = -t459 * pkin(10) + (-pkin(4) - pkin(5)) * t430 + (-pkin(4) * t477 + t445 + t525) * t464 - t526;
t503 = -m(7) * t382 + t394 * mrSges(7,1) - t395 * mrSges(7,2) + t432 * t414 - t433 * t415;
t372 = m(6) * t386 + t430 * mrSges(6,1) - t431 * mrSges(6,3) - t464 * t443 + t463 * t444 + t503;
t496 = m(5) * t500 - t430 * mrSges(5,1) - t431 * mrSges(5,2) - t463 * t441 - t464 * t442 - t372;
t371 = m(4) * t407 + qJDD(3) * mrSges(4,1) - t467 * mrSges(4,3) + qJD(3) * t472 - t465 * t512 + t496;
t507 = t491 * t361 - t371 * t488;
t362 = t487 * t364 + t524 * t365;
t497 = -m(4) * t412 + t468 * mrSges(4,1) - t467 * mrSges(4,2) - t471 * t512 + t472 * t511 - t362;
t452 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t488 + Ifges(4,4) * t491) * qJD(2);
t451 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t488 + Ifges(4,2) * t491) * qJD(2);
t402 = Ifges(7,5) * t433 + Ifges(7,6) * t432 + Ifges(7,3) * t475;
t370 = mrSges(7,2) * t382 - mrSges(7,3) * t377 + Ifges(7,1) * t395 + Ifges(7,4) * t394 + Ifges(7,5) * t456 + t402 * t432 - t403 * t475;
t369 = -mrSges(7,1) * t382 + mrSges(7,3) * t378 + Ifges(7,4) * t395 + Ifges(7,2) * t394 + Ifges(7,6) * t456 - t402 * t433 + t404 * t475;
t359 = -mrSges(5,2) * t500 + mrSges(6,2) * t385 - mrSges(5,3) * t387 - mrSges(6,3) * t386 - pkin(10) * t368 - qJ(5) * t372 - t486 * t369 + t490 * t370 - t522 * t430 + t533 * t431 + t521 * t460 + t516 * t463 - t515 * t477;
t358 = mrSges(5,1) * t500 - mrSges(6,1) * t386 + mrSges(6,2) * t384 + mrSges(5,3) * t388 - pkin(4) * t372 - pkin(5) * t503 - pkin(10) * t505 - t490 * t369 - t486 * t370 + t532 * t430 + t522 * t431 + t520 * t460 + t516 * t464 + t514 * t477;
t1 = [m(2) * t481 + t485 * (m(3) * t446 + t361 * t488 + t371 * t491) + (t489 * (m(3) * t419 - mrSges(3,1) * t494 - qJDD(2) * mrSges(3,2) + t507) + t492 * (m(3) * t418 + qJDD(2) * mrSges(3,1) - t494 * mrSges(3,2) + t497)) * t483; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t418 - mrSges(3,2) * t419 + t488 * (mrSges(4,2) * t412 - mrSges(4,3) * t407 + Ifges(4,1) * t467 + Ifges(4,4) * t468 + Ifges(4,5) * qJDD(3) - pkin(9) * t362 - qJD(3) * t451 - t487 * t358 + t359 * t524) + t491 * (-mrSges(4,1) * t412 + mrSges(4,3) * t408 + Ifges(4,4) * t467 + Ifges(4,2) * t468 + Ifges(4,6) * qJDD(3) - pkin(3) * t362 + qJD(3) * t452 - t530) + pkin(2) * t497 + pkin(8) * t507; Ifges(4,5) * t467 + Ifges(4,6) * t468 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t407 - mrSges(4,2) * t408 + t487 * t359 + t524 * t358 + pkin(3) * t496 + pkin(9) * t506 + (t451 * t488 - t452 * t491) * qJD(2); t530; t367; -t498;];
tauJ  = t1;
