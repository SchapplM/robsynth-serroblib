% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:20:50
% EndTime: 2019-05-05 15:20:53
% DurationCPUTime: 3.04s
% Computational Cost: add. (31128->271), mult. (69696->345), div. (0->0), fcn. (49977->12), ass. (0->123)
t480 = qJD(1) ^ 2;
t469 = cos(pkin(11));
t506 = pkin(3) * t469;
t467 = sin(pkin(11));
t505 = mrSges(4,2) * t467;
t465 = t469 ^ 2;
t504 = t465 * t480;
t474 = sin(qJ(1));
t478 = cos(qJ(1));
t495 = t474 * g(1) - g(2) * t478;
t452 = qJDD(1) * pkin(1) + t495;
t490 = -g(1) * t478 - g(2) * t474;
t453 = -pkin(1) * t480 + t490;
t468 = sin(pkin(10));
t470 = cos(pkin(10));
t436 = t468 * t452 + t470 * t453;
t426 = -pkin(2) * t480 + qJDD(1) * qJ(3) + t436;
t466 = -g(3) + qJDD(2);
t498 = qJD(1) * qJD(3);
t502 = t469 * t466 - 0.2e1 * t467 * t498;
t406 = (-pkin(7) * qJDD(1) + t480 * t506 - t426) * t467 + t502;
t414 = t467 * t466 + (t426 + 0.2e1 * t498) * t469;
t496 = qJDD(1) * t469;
t409 = -pkin(3) * t504 + pkin(7) * t496 + t414;
t473 = sin(qJ(4));
t477 = cos(qJ(4));
t390 = t473 * t406 + t477 * t409;
t500 = qJD(1) * t469;
t501 = qJD(1) * t467;
t445 = -t473 * t501 + t477 * t500;
t488 = t467 * t477 + t469 * t473;
t446 = t488 * qJD(1);
t429 = -mrSges(5,1) * t445 + mrSges(5,2) * t446;
t443 = t446 * qJD(4);
t497 = qJDD(1) * t467;
t433 = -t473 * t497 + t477 * t496 - t443;
t441 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t446;
t432 = -pkin(4) * t445 - pkin(8) * t446;
t479 = qJD(4) ^ 2;
t384 = -pkin(4) * t479 + qJDD(4) * pkin(8) + t432 * t445 + t390;
t464 = t467 ^ 2;
t435 = t470 * t452 - t468 * t453;
t489 = qJDD(3) - t435;
t410 = (-pkin(2) - t506) * qJDD(1) + (-qJ(3) + (-t464 - t465) * pkin(7)) * t480 + t489;
t499 = t445 * qJD(4);
t434 = qJDD(1) * t488 + t499;
t388 = (-t434 - t499) * pkin(8) + (-t433 + t443) * pkin(4) + t410;
t472 = sin(qJ(5));
t476 = cos(qJ(5));
t374 = -t472 * t384 + t476 * t388;
t438 = qJD(4) * t476 - t446 * t472;
t408 = qJD(5) * t438 + qJDD(4) * t472 + t434 * t476;
t431 = qJDD(5) - t433;
t439 = qJD(4) * t472 + t446 * t476;
t444 = qJD(5) - t445;
t372 = (t438 * t444 - t408) * pkin(9) + (t438 * t439 + t431) * pkin(5) + t374;
t375 = t476 * t384 + t472 * t388;
t407 = -qJD(5) * t439 + qJDD(4) * t476 - t434 * t472;
t419 = pkin(5) * t444 - pkin(9) * t439;
t437 = t438 ^ 2;
t373 = -pkin(5) * t437 + pkin(9) * t407 - t419 * t444 + t375;
t471 = sin(qJ(6));
t475 = cos(qJ(6));
t370 = t372 * t475 - t373 * t471;
t411 = t438 * t475 - t439 * t471;
t382 = qJD(6) * t411 + t407 * t471 + t408 * t475;
t412 = t438 * t471 + t439 * t475;
t395 = -mrSges(7,1) * t411 + mrSges(7,2) * t412;
t442 = qJD(6) + t444;
t396 = -mrSges(7,2) * t442 + mrSges(7,3) * t411;
t428 = qJDD(6) + t431;
t367 = m(7) * t370 + mrSges(7,1) * t428 - mrSges(7,3) * t382 - t395 * t412 + t396 * t442;
t371 = t372 * t471 + t373 * t475;
t381 = -qJD(6) * t412 + t407 * t475 - t408 * t471;
t397 = mrSges(7,1) * t442 - mrSges(7,3) * t412;
t368 = m(7) * t371 - mrSges(7,2) * t428 + mrSges(7,3) * t381 + t395 * t411 - t397 * t442;
t359 = t475 * t367 + t471 * t368;
t415 = -mrSges(6,1) * t438 + mrSges(6,2) * t439;
t417 = -mrSges(6,2) * t444 + mrSges(6,3) * t438;
t357 = m(6) * t374 + mrSges(6,1) * t431 - mrSges(6,3) * t408 - t415 * t439 + t417 * t444 + t359;
t418 = mrSges(6,1) * t444 - mrSges(6,3) * t439;
t491 = -t367 * t471 + t475 * t368;
t358 = m(6) * t375 - mrSges(6,2) * t431 + mrSges(6,3) * t407 + t415 * t438 - t418 * t444 + t491;
t492 = -t357 * t472 + t476 * t358;
t352 = m(5) * t390 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t433 - qJD(4) * t441 + t429 * t445 + t492;
t389 = t406 * t477 - t473 * t409;
t440 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t445;
t383 = -qJDD(4) * pkin(4) - pkin(8) * t479 + t446 * t432 - t389;
t376 = -pkin(5) * t407 - pkin(9) * t437 + t419 * t439 + t383;
t486 = m(7) * t376 - t381 * mrSges(7,1) + mrSges(7,2) * t382 - t411 * t396 + t397 * t412;
t482 = -m(6) * t383 + t407 * mrSges(6,1) - mrSges(6,2) * t408 + t438 * t417 - t418 * t439 - t486;
t363 = m(5) * t389 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t434 + qJD(4) * t440 - t429 * t446 + t482;
t503 = t473 * t352 + t477 * t363;
t353 = t476 * t357 + t472 * t358;
t413 = -t426 * t467 + t502;
t487 = mrSges(4,3) * qJDD(1) + t480 * (-mrSges(4,1) * t469 + t505);
t345 = m(4) * t413 - t467 * t487 + t503;
t493 = t477 * t352 - t473 * t363;
t346 = m(4) * t414 + t469 * t487 + t493;
t494 = -t345 * t467 + t469 * t346;
t392 = Ifges(7,4) * t412 + Ifges(7,2) * t411 + Ifges(7,6) * t442;
t393 = Ifges(7,1) * t412 + Ifges(7,4) * t411 + Ifges(7,5) * t442;
t485 = -mrSges(7,1) * t370 + mrSges(7,2) * t371 - Ifges(7,5) * t382 - Ifges(7,6) * t381 - Ifges(7,3) * t428 - t412 * t392 + t411 * t393;
t484 = m(5) * t410 - t433 * mrSges(5,1) + t434 * mrSges(5,2) - t445 * t440 + t446 * t441 + t353;
t421 = -qJDD(1) * pkin(2) - t480 * qJ(3) + t489;
t483 = -m(4) * t421 + mrSges(4,1) * t496 - t484 + (t464 * t480 + t504) * mrSges(4,3);
t400 = Ifges(6,4) * t439 + Ifges(6,2) * t438 + Ifges(6,6) * t444;
t401 = Ifges(6,1) * t439 + Ifges(6,4) * t438 + Ifges(6,5) * t444;
t481 = mrSges(6,1) * t374 - mrSges(6,2) * t375 + Ifges(6,5) * t408 + Ifges(6,6) * t407 + Ifges(6,3) * t431 + pkin(5) * t359 + t439 * t400 - t438 * t401 - t485;
t451 = (Ifges(4,5) * t467 + Ifges(4,6) * t469) * qJD(1);
t424 = Ifges(5,1) * t446 + Ifges(5,4) * t445 + Ifges(5,5) * qJD(4);
t423 = Ifges(5,4) * t446 + Ifges(5,2) * t445 + Ifges(5,6) * qJD(4);
t422 = Ifges(5,5) * t446 + Ifges(5,6) * t445 + Ifges(5,3) * qJD(4);
t399 = Ifges(6,5) * t439 + Ifges(6,6) * t438 + Ifges(6,3) * t444;
t391 = Ifges(7,5) * t412 + Ifges(7,6) * t411 + Ifges(7,3) * t442;
t361 = mrSges(7,2) * t376 - mrSges(7,3) * t370 + Ifges(7,1) * t382 + Ifges(7,4) * t381 + Ifges(7,5) * t428 + t391 * t411 - t392 * t442;
t360 = -mrSges(7,1) * t376 + mrSges(7,3) * t371 + Ifges(7,4) * t382 + Ifges(7,2) * t381 + Ifges(7,6) * t428 - t391 * t412 + t393 * t442;
t349 = mrSges(4,2) * t497 - t483;
t348 = mrSges(6,2) * t383 - mrSges(6,3) * t374 + Ifges(6,1) * t408 + Ifges(6,4) * t407 + Ifges(6,5) * t431 - pkin(9) * t359 - t360 * t471 + t361 * t475 + t399 * t438 - t400 * t444;
t347 = -mrSges(6,1) * t383 + mrSges(6,3) * t375 + Ifges(6,4) * t408 + Ifges(6,2) * t407 + Ifges(6,6) * t431 - pkin(5) * t486 + pkin(9) * t491 + t475 * t360 + t471 * t361 - t439 * t399 + t444 * t401;
t343 = -mrSges(5,1) * t410 + mrSges(5,3) * t390 + Ifges(5,4) * t434 + Ifges(5,2) * t433 + Ifges(5,6) * qJDD(4) - pkin(4) * t353 + qJD(4) * t424 - t446 * t422 - t481;
t342 = mrSges(5,2) * t410 - mrSges(5,3) * t389 + Ifges(5,1) * t434 + Ifges(5,4) * t433 + Ifges(5,5) * qJDD(4) - pkin(8) * t353 - qJD(4) * t423 - t347 * t472 + t348 * t476 + t422 * t445;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t495 - mrSges(2,2) * t490 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t435 - mrSges(3,2) * t436 + t467 * (t451 * t500 + mrSges(4,2) * t421 - mrSges(4,3) * t413 + t477 * t342 - t473 * t343 - pkin(7) * t503 + (Ifges(4,1) * t467 + Ifges(4,4) * t469) * qJDD(1)) + t469 * (-t451 * t501 - mrSges(4,1) * t421 + mrSges(4,3) * t414 + t473 * t342 + t477 * t343 - pkin(3) * t484 + pkin(7) * t493 + (Ifges(4,4) * t467 + Ifges(4,2) * t469) * qJDD(1)) - pkin(2) * t349 + qJ(3) * t494 + pkin(1) * (t468 * (m(3) * t436 - mrSges(3,1) * t480 - qJDD(1) * mrSges(3,2) + t494) + t470 * (-t480 * mrSges(3,2) + m(3) * t435 + t483 + (mrSges(3,1) - t505) * qJDD(1))); m(3) * t466 + t345 * t469 + t346 * t467; t349; mrSges(5,1) * t389 - mrSges(5,2) * t390 + Ifges(5,5) * t434 + Ifges(5,6) * t433 + Ifges(5,3) * qJDD(4) + pkin(4) * t482 + pkin(8) * t492 + t476 * t347 + t472 * t348 + t446 * t423 - t445 * t424; t481; -t485;];
tauJ  = t1;
