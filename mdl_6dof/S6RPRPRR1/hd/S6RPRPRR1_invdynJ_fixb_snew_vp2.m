% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:13:58
% EndTime: 2019-05-05 18:14:01
% DurationCPUTime: 3.28s
% Computational Cost: add. (36643->292), mult. (80830->375), div. (0->0), fcn. (56480->12), ass. (0->119)
t476 = sin(qJ(1));
t480 = cos(qJ(1));
t493 = t476 * g(1) - g(2) * t480;
t452 = qJDD(1) * pkin(1) + t493;
t481 = qJD(1) ^ 2;
t488 = -g(1) * t480 - g(2) * t476;
t454 = -pkin(1) * t481 + t488;
t470 = sin(pkin(10));
t472 = cos(pkin(10));
t432 = t470 * t452 + t472 * t454;
t426 = -pkin(2) * t481 + qJDD(1) * pkin(7) + t432;
t468 = -g(3) + qJDD(2);
t475 = sin(qJ(3));
t479 = cos(qJ(3));
t414 = -t426 * t475 + t479 * t468;
t495 = qJD(1) * qJD(3);
t494 = t479 * t495;
t455 = qJDD(1) * t475 + t494;
t409 = (-t455 + t494) * qJ(4) + (t475 * t479 * t481 + qJDD(3)) * pkin(3) + t414;
t415 = t479 * t426 + t475 * t468;
t456 = qJDD(1) * t479 - t475 * t495;
t497 = qJD(1) * t475;
t457 = qJD(3) * pkin(3) - qJ(4) * t497;
t467 = t479 ^ 2;
t410 = -pkin(3) * t467 * t481 + qJ(4) * t456 - qJD(3) * t457 + t415;
t469 = sin(pkin(11));
t471 = cos(pkin(11));
t442 = (t469 * t479 + t471 * t475) * qJD(1);
t377 = -0.2e1 * qJD(4) * t442 + t471 * t409 - t410 * t469;
t434 = t455 * t471 + t456 * t469;
t441 = (-t469 * t475 + t471 * t479) * qJD(1);
t374 = (qJD(3) * t441 - t434) * pkin(8) + (t441 * t442 + qJDD(3)) * pkin(4) + t377;
t378 = 0.2e1 * qJD(4) * t441 + t469 * t409 + t471 * t410;
t433 = -t455 * t469 + t456 * t471;
t437 = qJD(3) * pkin(4) - pkin(8) * t442;
t440 = t441 ^ 2;
t376 = -pkin(4) * t440 + pkin(8) * t433 - qJD(3) * t437 + t378;
t474 = sin(qJ(5));
t478 = cos(qJ(5));
t371 = t474 * t374 + t478 * t376;
t424 = t441 * t474 + t442 * t478;
t394 = -qJD(5) * t424 + t433 * t478 - t434 * t474;
t423 = t441 * t478 - t442 * t474;
t407 = -mrSges(6,1) * t423 + mrSges(6,2) * t424;
t466 = qJD(3) + qJD(5);
t417 = mrSges(6,1) * t466 - mrSges(6,3) * t424;
t465 = qJDD(3) + qJDD(5);
t408 = -pkin(5) * t423 - pkin(9) * t424;
t464 = t466 ^ 2;
t368 = -pkin(5) * t464 + pkin(9) * t465 + t408 * t423 + t371;
t431 = t452 * t472 - t470 * t454;
t487 = -qJDD(1) * pkin(2) - t431;
t411 = -pkin(3) * t456 + qJDD(4) + t457 * t497 + (-qJ(4) * t467 - pkin(7)) * t481 + t487;
t383 = -pkin(4) * t433 - pkin(8) * t440 + t442 * t437 + t411;
t395 = qJD(5) * t423 + t433 * t474 + t434 * t478;
t372 = (-t423 * t466 - t395) * pkin(9) + (t424 * t466 - t394) * pkin(5) + t383;
t473 = sin(qJ(6));
t477 = cos(qJ(6));
t365 = -t368 * t473 + t372 * t477;
t412 = -t424 * t473 + t466 * t477;
t381 = qJD(6) * t412 + t395 * t477 + t465 * t473;
t393 = qJDD(6) - t394;
t413 = t424 * t477 + t466 * t473;
t396 = -mrSges(7,1) * t412 + mrSges(7,2) * t413;
t419 = qJD(6) - t423;
t397 = -mrSges(7,2) * t419 + mrSges(7,3) * t412;
t362 = m(7) * t365 + mrSges(7,1) * t393 - mrSges(7,3) * t381 - t396 * t413 + t397 * t419;
t366 = t368 * t477 + t372 * t473;
t380 = -qJD(6) * t413 - t395 * t473 + t465 * t477;
t398 = mrSges(7,1) * t419 - mrSges(7,3) * t413;
t363 = m(7) * t366 - mrSges(7,2) * t393 + mrSges(7,3) * t380 + t396 * t412 - t398 * t419;
t489 = -t362 * t473 + t477 * t363;
t349 = m(6) * t371 - mrSges(6,2) * t465 + mrSges(6,3) * t394 + t407 * t423 - t417 * t466 + t489;
t370 = t374 * t478 - t376 * t474;
t416 = -mrSges(6,2) * t466 + mrSges(6,3) * t423;
t367 = -pkin(5) * t465 - pkin(9) * t464 + t408 * t424 - t370;
t485 = -m(7) * t367 + t380 * mrSges(7,1) - mrSges(7,2) * t381 + t412 * t397 - t398 * t413;
t358 = m(6) * t370 + mrSges(6,1) * t465 - mrSges(6,3) * t395 - t407 * t424 + t416 * t466 + t485;
t346 = t474 * t349 + t478 * t358;
t429 = -mrSges(5,1) * t441 + mrSges(5,2) * t442;
t435 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t441;
t344 = m(5) * t377 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t434 + qJD(3) * t435 - t429 * t442 + t346;
t436 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t442;
t490 = t478 * t349 - t358 * t474;
t345 = m(5) * t378 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t433 - qJD(3) * t436 + t429 * t441 + t490;
t338 = t471 * t344 + t469 * t345;
t352 = t477 * t362 + t473 * t363;
t496 = qJD(1) * t479;
t453 = (-mrSges(4,1) * t479 + mrSges(4,2) * t475) * qJD(1);
t459 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t496;
t336 = m(4) * t414 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t455 + qJD(3) * t459 - t453 * t497 + t338;
t458 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t497;
t491 = -t344 * t469 + t471 * t345;
t337 = m(4) * t415 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t456 - qJD(3) * t458 + t453 * t496 + t491;
t492 = -t336 * t475 + t479 * t337;
t486 = m(6) * t383 - t394 * mrSges(6,1) + t395 * mrSges(6,2) - t423 * t416 + t424 * t417 + t352;
t384 = Ifges(7,5) * t413 + Ifges(7,6) * t412 + Ifges(7,3) * t419;
t386 = Ifges(7,1) * t413 + Ifges(7,4) * t412 + Ifges(7,5) * t419;
t355 = -mrSges(7,1) * t367 + mrSges(7,3) * t366 + Ifges(7,4) * t381 + Ifges(7,2) * t380 + Ifges(7,6) * t393 - t384 * t413 + t386 * t419;
t385 = Ifges(7,4) * t413 + Ifges(7,2) * t412 + Ifges(7,6) * t419;
t356 = mrSges(7,2) * t367 - mrSges(7,3) * t365 + Ifges(7,1) * t381 + Ifges(7,4) * t380 + Ifges(7,5) * t393 + t384 * t412 - t385 * t419;
t400 = Ifges(6,4) * t424 + Ifges(6,2) * t423 + Ifges(6,6) * t466;
t401 = Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t466;
t484 = mrSges(6,1) * t370 - mrSges(6,2) * t371 + Ifges(6,5) * t395 + Ifges(6,6) * t394 + Ifges(6,3) * t465 + pkin(5) * t485 + pkin(9) * t489 + t477 * t355 + t473 * t356 + t424 * t400 - t423 * t401;
t350 = m(5) * t411 - t433 * mrSges(5,1) + mrSges(5,2) * t434 - t441 * t435 + t436 * t442 + t486;
t483 = mrSges(7,1) * t365 - mrSges(7,2) * t366 + Ifges(7,5) * t381 + Ifges(7,6) * t380 + Ifges(7,3) * t393 + t385 * t413 - t386 * t412;
t425 = -pkin(7) * t481 + t487;
t482 = -m(4) * t425 + t456 * mrSges(4,1) - mrSges(4,2) * t455 - t458 * t497 + t459 * t496 - t350;
t448 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t475 + Ifges(4,4) * t479) * qJD(1);
t447 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t475 + Ifges(4,2) * t479) * qJD(1);
t422 = Ifges(5,1) * t442 + Ifges(5,4) * t441 + Ifges(5,5) * qJD(3);
t421 = Ifges(5,4) * t442 + Ifges(5,2) * t441 + Ifges(5,6) * qJD(3);
t420 = Ifges(5,5) * t442 + Ifges(5,6) * t441 + Ifges(5,3) * qJD(3);
t399 = Ifges(6,5) * t424 + Ifges(6,6) * t423 + Ifges(6,3) * t466;
t340 = -mrSges(6,1) * t383 + mrSges(6,3) * t371 + Ifges(6,4) * t395 + Ifges(6,2) * t394 + Ifges(6,6) * t465 - pkin(5) * t352 - t399 * t424 + t401 * t466 - t483;
t339 = mrSges(6,2) * t383 - mrSges(6,3) * t370 + Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t465 - pkin(9) * t352 - t355 * t473 + t356 * t477 + t399 * t423 - t400 * t466;
t334 = mrSges(5,2) * t411 - mrSges(5,3) * t377 + Ifges(5,1) * t434 + Ifges(5,4) * t433 + Ifges(5,5) * qJDD(3) - pkin(8) * t346 - qJD(3) * t421 + t339 * t478 - t340 * t474 + t420 * t441;
t333 = -mrSges(5,1) * t411 + mrSges(5,3) * t378 + Ifges(5,4) * t434 + Ifges(5,2) * t433 + Ifges(5,6) * qJDD(3) - pkin(4) * t486 + pkin(8) * t490 + qJD(3) * t422 + t474 * t339 + t478 * t340 - t442 * t420;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t493 - mrSges(2,2) * t488 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t431 - mrSges(3,2) * t432 + t475 * (mrSges(4,2) * t425 - mrSges(4,3) * t414 + Ifges(4,1) * t455 + Ifges(4,4) * t456 + Ifges(4,5) * qJDD(3) - qJ(4) * t338 - qJD(3) * t447 - t469 * t333 + t471 * t334) + t479 * (-mrSges(4,1) * t425 + mrSges(4,3) * t415 + Ifges(4,4) * t455 + Ifges(4,2) * t456 + Ifges(4,6) * qJDD(3) - pkin(3) * t350 + qJ(4) * t491 + qJD(3) * t448 + t471 * t333 + t469 * t334) + pkin(2) * t482 + pkin(7) * t492 + pkin(1) * (t470 * (m(3) * t432 - mrSges(3,1) * t481 - qJDD(1) * mrSges(3,2) + t492) + t472 * (m(3) * t431 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t481 + t482)); m(3) * t468 + t336 * t479 + t337 * t475; Ifges(4,5) * t455 + Ifges(4,6) * t456 - t441 * t422 + t442 * t421 + Ifges(5,6) * t433 + Ifges(5,5) * t434 + mrSges(4,1) * t414 - mrSges(4,2) * t415 + mrSges(5,1) * t377 - mrSges(5,2) * t378 + pkin(4) * t346 + pkin(3) * t338 + t484 + (t475 * t447 - t479 * t448) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3); t350; t484; t483;];
tauJ  = t1;
