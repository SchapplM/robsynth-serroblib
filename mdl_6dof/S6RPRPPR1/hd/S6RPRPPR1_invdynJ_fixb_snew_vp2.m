% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:24:18
% EndTime: 2019-05-05 16:24:22
% DurationCPUTime: 3.30s
% Computational Cost: add. (32207->291), mult. (72025->374), div. (0->0), fcn. (48449->12), ass. (0->117)
t497 = -2 * qJD(4);
t474 = sin(qJ(1));
t477 = cos(qJ(1));
t490 = t474 * g(1) - g(2) * t477;
t453 = qJDD(1) * pkin(1) + t490;
t479 = qJD(1) ^ 2;
t485 = -g(1) * t477 - g(2) * t474;
t455 = -pkin(1) * t479 + t485;
t469 = sin(pkin(9));
t471 = cos(pkin(9));
t428 = t469 * t453 + t471 * t455;
t421 = -pkin(2) * t479 + qJDD(1) * pkin(7) + t428;
t466 = -g(3) + qJDD(2);
t473 = sin(qJ(3));
t476 = cos(qJ(3));
t408 = -t473 * t421 + t476 * t466;
t492 = qJD(1) * qJD(3);
t491 = t476 * t492;
t456 = qJDD(1) * t473 + t491;
t395 = (-t456 + t491) * qJ(4) + (t473 * t476 * t479 + qJDD(3)) * pkin(3) + t408;
t409 = t476 * t421 + t473 * t466;
t457 = qJDD(1) * t476 - t473 * t492;
t495 = qJD(1) * t473;
t458 = qJD(3) * pkin(3) - qJ(4) * t495;
t465 = t476 ^ 2;
t398 = -pkin(3) * t465 * t479 + qJ(4) * t457 - qJD(3) * t458 + t409;
t468 = sin(pkin(10));
t496 = cos(pkin(10));
t442 = (t468 * t476 + t496 * t473) * qJD(1);
t379 = t496 * t395 - t468 * t398 + t442 * t497;
t494 = qJD(1) * t476;
t441 = t468 * t495 - t496 * t494;
t380 = t468 * t395 + t496 * t398 + t441 * t497;
t424 = mrSges(5,1) * t441 + mrSges(5,2) * t442;
t429 = t456 * t468 - t496 * t457;
t437 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t442;
t423 = pkin(4) * t441 - qJ(5) * t442;
t478 = qJD(3) ^ 2;
t378 = -pkin(4) * t478 + qJDD(3) * qJ(5) - t423 * t441 + t380;
t427 = t471 * t453 - t469 * t455;
t483 = -qJDD(1) * pkin(2) - t427;
t399 = -t457 * pkin(3) + qJDD(4) + t458 * t495 + (-qJ(4) * t465 - pkin(7)) * t479 + t483;
t430 = t496 * t456 + t468 * t457;
t383 = (qJD(3) * t441 - t430) * qJ(5) + (qJD(3) * t442 + t429) * pkin(4) + t399;
t467 = sin(pkin(11));
t470 = cos(pkin(11));
t435 = qJD(3) * t467 + t442 * t470;
t373 = -0.2e1 * qJD(5) * t435 - t467 * t378 + t470 * t383;
t416 = qJDD(3) * t467 + t430 * t470;
t434 = qJD(3) * t470 - t442 * t467;
t371 = (t434 * t441 - t416) * pkin(8) + (t434 * t435 + t429) * pkin(5) + t373;
t374 = 0.2e1 * qJD(5) * t434 + t470 * t378 + t467 * t383;
t412 = pkin(5) * t441 - pkin(8) * t435;
t415 = qJDD(3) * t470 - t430 * t467;
t433 = t434 ^ 2;
t372 = -pkin(5) * t433 + pkin(8) * t415 - t412 * t441 + t374;
t472 = sin(qJ(6));
t475 = cos(qJ(6));
t369 = t371 * t475 - t372 * t472;
t404 = t434 * t475 - t435 * t472;
t386 = qJD(6) * t404 + t415 * t472 + t416 * t475;
t405 = t434 * t472 + t435 * t475;
t391 = -mrSges(7,1) * t404 + mrSges(7,2) * t405;
t440 = qJD(6) + t441;
t396 = -mrSges(7,2) * t440 + mrSges(7,3) * t404;
t426 = qJDD(6) + t429;
t366 = m(7) * t369 + mrSges(7,1) * t426 - mrSges(7,3) * t386 - t391 * t405 + t396 * t440;
t370 = t371 * t472 + t372 * t475;
t385 = -qJD(6) * t405 + t415 * t475 - t416 * t472;
t397 = mrSges(7,1) * t440 - mrSges(7,3) * t405;
t367 = m(7) * t370 - mrSges(7,2) * t426 + mrSges(7,3) * t385 + t391 * t404 - t397 * t440;
t358 = t475 * t366 + t472 * t367;
t406 = -mrSges(6,1) * t434 + mrSges(6,2) * t435;
t410 = -mrSges(6,2) * t441 + mrSges(6,3) * t434;
t356 = m(6) * t373 + mrSges(6,1) * t429 - mrSges(6,3) * t416 - t406 * t435 + t410 * t441 + t358;
t411 = mrSges(6,1) * t441 - mrSges(6,3) * t435;
t486 = -t366 * t472 + t475 * t367;
t357 = m(6) * t374 - mrSges(6,2) * t429 + mrSges(6,3) * t415 + t406 * t434 - t411 * t441 + t486;
t487 = -t356 * t467 + t470 * t357;
t350 = m(5) * t380 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t429 - qJD(3) * t437 - t424 * t441 + t487;
t377 = -qJDD(3) * pkin(4) - t478 * qJ(5) + t442 * t423 + qJDD(5) - t379;
t375 = -t415 * pkin(5) - t433 * pkin(8) + t435 * t412 + t377;
t482 = m(7) * t375 - t385 * mrSges(7,1) + mrSges(7,2) * t386 - t404 * t396 + t397 * t405;
t368 = m(6) * t377 - t415 * mrSges(6,1) + mrSges(6,2) * t416 - t434 * t410 + t411 * t435 + t482;
t436 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t441;
t362 = m(5) * t379 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t430 + qJD(3) * t436 - t424 * t442 - t368;
t345 = t468 * t350 + t496 * t362;
t352 = t470 * t356 + t467 * t357;
t454 = (-mrSges(4,1) * t476 + mrSges(4,2) * t473) * qJD(1);
t460 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t494;
t343 = m(4) * t408 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t456 + qJD(3) * t460 - t454 * t495 + t345;
t459 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t495;
t488 = t496 * t350 - t362 * t468;
t344 = m(4) * t409 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t457 - qJD(3) * t459 + t454 * t494 + t488;
t489 = -t343 * t473 + t476 * t344;
t351 = m(5) * t399 + t429 * mrSges(5,1) + mrSges(5,2) * t430 + t441 * t436 + t437 * t442 + t352;
t388 = Ifges(7,4) * t405 + Ifges(7,2) * t404 + Ifges(7,6) * t440;
t389 = Ifges(7,1) * t405 + Ifges(7,4) * t404 + Ifges(7,5) * t440;
t481 = mrSges(7,1) * t369 - mrSges(7,2) * t370 + Ifges(7,5) * t386 + Ifges(7,6) * t385 + Ifges(7,3) * t426 + t405 * t388 - t404 * t389;
t420 = -t479 * pkin(7) + t483;
t480 = -m(4) * t420 + t457 * mrSges(4,1) - mrSges(4,2) * t456 - t459 * t495 + t460 * t494 - t351;
t448 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t473 + Ifges(4,4) * t476) * qJD(1);
t447 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t473 + Ifges(4,2) * t476) * qJD(1);
t419 = Ifges(5,1) * t442 - Ifges(5,4) * t441 + Ifges(5,5) * qJD(3);
t418 = Ifges(5,4) * t442 - Ifges(5,2) * t441 + Ifges(5,6) * qJD(3);
t417 = Ifges(5,5) * t442 - Ifges(5,6) * t441 + Ifges(5,3) * qJD(3);
t402 = Ifges(6,1) * t435 + Ifges(6,4) * t434 + Ifges(6,5) * t441;
t401 = Ifges(6,4) * t435 + Ifges(6,2) * t434 + Ifges(6,6) * t441;
t400 = Ifges(6,5) * t435 + Ifges(6,6) * t434 + Ifges(6,3) * t441;
t387 = Ifges(7,5) * t405 + Ifges(7,6) * t404 + Ifges(7,3) * t440;
t360 = mrSges(7,2) * t375 - mrSges(7,3) * t369 + Ifges(7,1) * t386 + Ifges(7,4) * t385 + Ifges(7,5) * t426 + t387 * t404 - t388 * t440;
t359 = -mrSges(7,1) * t375 + mrSges(7,3) * t370 + Ifges(7,4) * t386 + Ifges(7,2) * t385 + Ifges(7,6) * t426 - t387 * t405 + t389 * t440;
t347 = mrSges(6,2) * t377 - mrSges(6,3) * t373 + Ifges(6,1) * t416 + Ifges(6,4) * t415 + Ifges(6,5) * t429 - pkin(8) * t358 - t359 * t472 + t360 * t475 + t400 * t434 - t401 * t441;
t346 = -mrSges(6,1) * t377 + mrSges(6,3) * t374 + Ifges(6,4) * t416 + Ifges(6,2) * t415 + Ifges(6,6) * t429 - pkin(5) * t482 + pkin(8) * t486 + t475 * t359 + t472 * t360 - t435 * t400 + t441 * t402;
t341 = Ifges(5,6) * qJDD(3) + (-Ifges(5,2) - Ifges(6,3)) * t429 - t442 * t417 - t435 * t401 + Ifges(5,4) * t430 + t434 * t402 + qJD(3) * t419 - Ifges(6,6) * t415 - Ifges(6,5) * t416 - mrSges(5,1) * t399 + mrSges(5,3) * t380 - mrSges(6,1) * t373 + mrSges(6,2) * t374 - pkin(5) * t358 - pkin(4) * t352 - t481;
t340 = mrSges(5,2) * t399 - mrSges(5,3) * t379 + Ifges(5,1) * t430 - Ifges(5,4) * t429 + Ifges(5,5) * qJDD(3) - qJ(5) * t352 - qJD(3) * t418 - t346 * t467 + t347 * t470 - t417 * t441;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t490 - mrSges(2,2) * t485 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t427 - mrSges(3,2) * t428 + t473 * (mrSges(4,2) * t420 - mrSges(4,3) * t408 + Ifges(4,1) * t456 + Ifges(4,4) * t457 + Ifges(4,5) * qJDD(3) - qJ(4) * t345 - qJD(3) * t447 + t496 * t340 - t468 * t341) + t476 * (-mrSges(4,1) * t420 + mrSges(4,3) * t409 + Ifges(4,4) * t456 + Ifges(4,2) * t457 + Ifges(4,6) * qJDD(3) - pkin(3) * t351 + qJ(4) * t488 + qJD(3) * t448 + t468 * t340 + t496 * t341) + pkin(2) * t480 + pkin(7) * t489 + pkin(1) * (t469 * (m(3) * t428 - mrSges(3,1) * t479 - qJDD(1) * mrSges(3,2) + t489) + t471 * (m(3) * t427 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t479 + t480)); m(3) * t466 + t343 * t476 + t344 * t473; Ifges(4,5) * t456 + Ifges(4,6) * t457 + mrSges(4,1) * t408 - mrSges(4,2) * t409 + Ifges(5,5) * t430 - Ifges(5,6) * t429 + t442 * t418 + t441 * t419 + mrSges(5,1) * t379 - mrSges(5,2) * t380 + t467 * t347 + t470 * t346 - pkin(4) * t368 + qJ(5) * t487 + pkin(3) * t345 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t447 * t473 - t448 * t476) * qJD(1); t351; t368; t481;];
tauJ  = t1;
