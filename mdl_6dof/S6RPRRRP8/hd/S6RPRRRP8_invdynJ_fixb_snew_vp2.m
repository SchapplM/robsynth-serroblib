% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:46:35
% EndTime: 2019-05-06 01:46:39
% DurationCPUTime: 2.09s
% Computational Cost: add. (17203->266), mult. (33600->325), div. (0->0), fcn. (21988->8), ass. (0->107)
t503 = Ifges(6,1) + Ifges(7,1);
t496 = Ifges(6,4) - Ifges(7,5);
t495 = -Ifges(6,5) - Ifges(7,4);
t502 = Ifges(6,2) + Ifges(7,3);
t494 = Ifges(6,6) - Ifges(7,6);
t501 = -Ifges(6,3) - Ifges(7,2);
t464 = sin(qJ(1));
t467 = cos(qJ(1));
t478 = -g(1) * t467 - g(2) * t464;
t474 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t478;
t462 = sin(qJ(4));
t463 = sin(qJ(3));
t465 = cos(qJ(4));
t466 = cos(qJ(3));
t441 = (t466 * t462 + t463 * t465) * qJD(1);
t486 = qJD(1) * qJD(3);
t445 = -qJDD(1) * t463 - t466 * t486;
t483 = t463 * t486;
t446 = qJDD(1) * t466 - t483;
t413 = -qJD(4) * t441 + t445 * t462 + t446 * t465;
t442 = (-t463 * t462 + t466 * t465) * qJD(1);
t458 = qJD(3) + qJD(4);
t461 = sin(qJ(5));
t498 = cos(qJ(5));
t426 = t442 * t461 - t498 * t458;
t457 = qJDD(3) + qJDD(4);
t383 = -t426 * qJD(5) + t498 * t413 + t461 * t457;
t427 = t498 * t442 + t461 * t458;
t400 = mrSges(7,1) * t426 - mrSges(7,3) * t427;
t487 = t466 * qJD(1);
t449 = (qJD(3) * pkin(3)) - pkin(8) * t487;
t460 = t463 ^ 2;
t468 = qJD(1) ^ 2;
t499 = -pkin(1) - pkin(7);
t404 = -pkin(3) * t445 + t449 * t487 + (-pkin(8) * t460 + t499) * t468 + t474;
t412 = -qJD(4) * t442 + t445 * t465 - t446 * t462;
t373 = (t441 * t458 - t413) * pkin(9) + (t442 * t458 - t412) * pkin(4) + t404;
t482 = g(1) * t464 - t467 * g(2);
t473 = -qJ(2) * t468 + qJDD(2) - t482;
t433 = t499 * qJDD(1) + t473;
t424 = t463 * g(3) + t466 * t433;
t397 = (-t446 - t483) * pkin(8) + (-t463 * t466 * t468 + qJDD(3)) * pkin(3) + t424;
t425 = -g(3) * t466 + t463 * t433;
t398 = -pkin(3) * t460 * t468 + pkin(8) * t445 - qJD(3) * t449 + t425;
t379 = t462 * t397 + t465 * t398;
t423 = pkin(4) * t441 - pkin(9) * t442;
t456 = t458 ^ 2;
t376 = -pkin(4) * t456 + pkin(9) * t457 - t423 * t441 + t379;
t370 = t498 * t373 - t461 * t376;
t399 = pkin(5) * t426 - qJ(6) * t427;
t411 = qJDD(5) - t412;
t437 = qJD(5) + t441;
t435 = t437 ^ 2;
t368 = -t411 * pkin(5) - t435 * qJ(6) + t427 * t399 + qJDD(6) - t370;
t414 = -mrSges(7,2) * t426 + mrSges(7,3) * t437;
t479 = -m(7) * t368 + t411 * mrSges(7,1) + t437 * t414;
t364 = t383 * mrSges(7,2) + t400 * t427 - t479;
t371 = t461 * t373 + t498 * t376;
t367 = -pkin(5) * t435 + qJ(6) * t411 + 0.2e1 * qJD(6) * t437 - t399 * t426 + t371;
t382 = qJD(5) * t427 + t413 * t461 - t498 * t457;
t417 = -mrSges(7,1) * t437 + mrSges(7,2) * t427;
t484 = m(7) * t367 + t411 * mrSges(7,3) + t437 * t417;
t490 = t496 * t426 - t503 * t427 + t495 * t437;
t491 = t502 * t426 - t496 * t427 - t494 * t437;
t500 = -t494 * t382 - t495 * t383 - t501 * t411 - t490 * t426 - t491 * t427 + mrSges(6,1) * t370 - mrSges(7,1) * t368 - mrSges(6,2) * t371 + mrSges(7,3) * t367 - pkin(5) * t364 + qJ(6) * (-t382 * mrSges(7,2) - t400 * t426 + t484);
t497 = -mrSges(6,3) - mrSges(7,2);
t422 = mrSges(5,1) * t441 + mrSges(5,2) * t442;
t431 = mrSges(5,1) * t458 - mrSges(5,3) * t442;
t416 = mrSges(6,1) * t437 - mrSges(6,3) * t427;
t489 = -mrSges(6,1) * t426 - mrSges(6,2) * t427 - t400;
t359 = m(6) * t371 - mrSges(6,2) * t411 + t497 * t382 - t416 * t437 + t489 * t426 + t484;
t415 = -mrSges(6,2) * t437 - mrSges(6,3) * t426;
t361 = m(6) * t370 + mrSges(6,1) * t411 + t497 * t383 + t415 * t437 + t489 * t427 + t479;
t480 = t498 * t359 - t361 * t461;
t348 = m(5) * t379 - mrSges(5,2) * t457 + mrSges(5,3) * t412 - t422 * t441 - t431 * t458 + t480;
t378 = t397 * t465 - t462 * t398;
t430 = -mrSges(5,2) * t458 - mrSges(5,3) * t441;
t375 = -pkin(4) * t457 - pkin(9) * t456 + t442 * t423 - t378;
t369 = -0.2e1 * qJD(6) * t427 + (t426 * t437 - t383) * qJ(6) + (t427 * t437 + t382) * pkin(5) + t375;
t365 = m(7) * t369 + t382 * mrSges(7,1) - t383 * mrSges(7,3) + t414 * t426 - t427 * t417;
t469 = -m(6) * t375 - t382 * mrSges(6,1) - t383 * mrSges(6,2) - t426 * t415 - t416 * t427 - t365;
t356 = m(5) * t378 + mrSges(5,1) * t457 - mrSges(5,3) * t413 - t422 * t442 + t430 * t458 + t469;
t345 = t462 * t348 + t465 * t356;
t354 = t461 * t359 + t498 * t361;
t492 = t494 * t426 + t495 * t427 + t501 * t437;
t488 = qJD(1) * t463;
t481 = t465 * t348 - t462 * t356;
t444 = (t463 * mrSges(4,1) + t466 * mrSges(4,2)) * qJD(1);
t447 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t488;
t448 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t487;
t477 = (m(4) * t424 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t446 + qJD(3) * t447 - t444 * t487 + t345) * t466 + (m(4) * t425 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t445 - qJD(3) * t448 - t444 * t488 + t481) * t463;
t472 = m(5) * t404 - mrSges(5,1) * t412 + t413 * mrSges(5,2) + t441 * t430 + t442 * t431 + t354;
t350 = -mrSges(6,1) * t375 - mrSges(7,1) * t369 + mrSges(7,2) * t367 + mrSges(6,3) * t371 - pkin(5) * t365 - t502 * t382 + t496 * t383 + t494 * t411 + t492 * t427 - t490 * t437;
t352 = mrSges(6,2) * t375 + mrSges(7,2) * t368 - mrSges(6,3) * t370 - mrSges(7,3) * t369 - qJ(6) * t365 - t496 * t382 + t503 * t383 - t495 * t411 + t492 * t426 + t491 * t437;
t419 = Ifges(5,4) * t442 - Ifges(5,2) * t441 + Ifges(5,6) * t458;
t420 = Ifges(5,1) * t442 - Ifges(5,4) * t441 + Ifges(5,5) * t458;
t471 = mrSges(5,1) * t378 - mrSges(5,2) * t379 + Ifges(5,5) * t413 + Ifges(5,6) * t412 + Ifges(5,3) * t457 + pkin(4) * t469 + pkin(9) * t480 + t498 * t350 + t461 * t352 + t442 * t419 + t420 * t441;
t440 = (Ifges(4,5) * qJD(3)) + (t466 * Ifges(4,1) - t463 * Ifges(4,4)) * qJD(1);
t439 = (Ifges(4,6) * qJD(3)) + (t466 * Ifges(4,4) - t463 * Ifges(4,2)) * qJD(1);
t436 = -qJDD(1) * pkin(1) + t473;
t434 = pkin(1) * t468 - t474;
t432 = t499 * t468 + t474;
t418 = Ifges(5,5) * t442 - Ifges(5,6) * t441 + Ifges(5,3) * t458;
t342 = -mrSges(5,1) * t404 + mrSges(5,3) * t379 + Ifges(5,4) * t413 + Ifges(5,2) * t412 + Ifges(5,6) * t457 - pkin(4) * t354 - t442 * t418 + t458 * t420 - t500;
t341 = m(3) * t436 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t468) + t477;
t340 = mrSges(5,2) * t404 - mrSges(5,3) * t378 + Ifges(5,1) * t413 + Ifges(5,4) * t412 + Ifges(5,5) * t457 - pkin(9) * t354 - t461 * t350 + t498 * t352 - t441 * t418 - t458 * t419;
t1 = [mrSges(2,1) * t482 - mrSges(2,2) * t478 + mrSges(3,2) * t436 - mrSges(3,3) * t434 + t466 * (mrSges(4,2) * t432 - mrSges(4,3) * t424 + Ifges(4,1) * t446 + Ifges(4,4) * t445 + Ifges(4,5) * qJDD(3) - pkin(8) * t345 - qJD(3) * t439 + t340 * t465 - t342 * t462) - t463 * (-mrSges(4,1) * t432 + mrSges(4,3) * t425 + Ifges(4,4) * t446 + Ifges(4,2) * t445 + Ifges(4,6) * qJDD(3) - pkin(3) * t472 + pkin(8) * t481 + qJD(3) * t440 + t462 * t340 + t465 * t342) - pkin(7) * t477 - pkin(1) * t341 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t434 + m(4) * t432 - mrSges(4,1) * t445 + mrSges(3,2) * t468 + mrSges(4,2) * t446 + t472 + qJDD(1) * mrSges(3,3) + (t447 * t463 + t448 * t466) * qJD(1)) * qJ(2); t341; t471 + (t466 * t439 + t463 * t440) * qJD(1) + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t424 - mrSges(4,2) * t425 + Ifges(4,5) * t446 + Ifges(4,6) * t445 + pkin(3) * t345; t471; t500; t364;];
tauJ  = t1;
