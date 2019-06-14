% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 20:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:01:59
% EndTime: 2019-05-05 20:02:04
% DurationCPUTime: 3.76s
% Computational Cost: add. (39043->289), mult. (82546->363), div. (0->0), fcn. (55236->10), ass. (0->115)
t466 = sin(qJ(1));
t470 = cos(qJ(1));
t480 = -g(1) * t470 - g(2) * t466;
t491 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t480;
t490 = (-pkin(1) - pkin(7));
t472 = qJD(1) ^ 2;
t428 = (t472 * t490) - t491;
t465 = sin(qJ(3));
t469 = cos(qJ(3));
t488 = qJD(1) * qJD(3);
t485 = t469 * t488;
t451 = qJDD(1) * t465 + t485;
t486 = t465 * t488;
t452 = qJDD(1) * t469 - t486;
t408 = (-t452 + t486) * qJ(4) + (t451 + t485) * pkin(3) + t428;
t484 = g(1) * t466 - t470 * g(2);
t476 = -qJ(2) * t472 + qJDD(2) - t484;
t434 = qJDD(1) * t490 + t476;
t425 = -g(3) * t469 + t465 * t434;
t449 = (t465 * pkin(3) - t469 * qJ(4)) * qJD(1);
t457 = t465 * qJD(1);
t471 = qJD(3) ^ 2;
t411 = -pkin(3) * t471 + qJDD(3) * qJ(4) - t449 * t457 + t425;
t461 = sin(pkin(10));
t462 = cos(pkin(10));
t489 = qJD(1) * t469;
t447 = qJD(3) * t461 + t462 * t489;
t391 = -0.2e1 * qJD(4) * t447 + t462 * t408 - t411 * t461;
t432 = qJDD(3) * t461 + t452 * t462;
t446 = qJD(3) * t462 - t461 * t489;
t382 = (t446 * t457 - t432) * pkin(8) + (t446 * t447 + t451) * pkin(4) + t391;
t392 = 0.2e1 * qJD(4) * t446 + t461 * t408 + t462 * t411;
t431 = qJDD(3) * t462 - t452 * t461;
t433 = pkin(4) * t457 - pkin(8) * t447;
t445 = t446 ^ 2;
t384 = -pkin(4) * t445 + pkin(8) * t431 - t433 * t457 + t392;
t464 = sin(qJ(5));
t468 = cos(qJ(5));
t369 = t468 * t382 - t384 * t464;
t421 = t446 * t468 - t447 * t464;
t398 = qJD(5) * t421 + t431 * t464 + t432 * t468;
t422 = t446 * t464 + t447 * t468;
t448 = qJDD(5) + t451;
t456 = t457 + qJD(5);
t367 = (t421 * t456 - t398) * pkin(9) + (t421 * t422 + t448) * pkin(5) + t369;
t370 = t464 * t382 + t468 * t384;
t397 = -qJD(5) * t422 + t431 * t468 - t432 * t464;
t414 = pkin(5) * t456 - pkin(9) * t422;
t420 = t421 ^ 2;
t368 = -pkin(5) * t420 + pkin(9) * t397 - t414 * t456 + t370;
t463 = sin(qJ(6));
t467 = cos(qJ(6));
t365 = t367 * t467 - t368 * t463;
t403 = t421 * t467 - t422 * t463;
t378 = qJD(6) * t403 + t397 * t463 + t398 * t467;
t404 = t421 * t463 + t422 * t467;
t390 = -mrSges(7,1) * t403 + mrSges(7,2) * t404;
t455 = qJD(6) + t456;
t394 = -mrSges(7,2) * t455 + mrSges(7,3) * t403;
t444 = qJDD(6) + t448;
t361 = m(7) * t365 + mrSges(7,1) * t444 - mrSges(7,3) * t378 - t390 * t404 + t394 * t455;
t366 = t367 * t463 + t368 * t467;
t377 = -qJD(6) * t404 + t397 * t467 - t398 * t463;
t395 = mrSges(7,1) * t455 - mrSges(7,3) * t404;
t362 = m(7) * t366 - mrSges(7,2) * t444 + mrSges(7,3) * t377 + t390 * t403 - t395 * t455;
t354 = t467 * t361 + t463 * t362;
t405 = -mrSges(6,1) * t421 + mrSges(6,2) * t422;
t412 = -mrSges(6,2) * t456 + mrSges(6,3) * t421;
t352 = m(6) * t369 + mrSges(6,1) * t448 - mrSges(6,3) * t398 - t405 * t422 + t412 * t456 + t354;
t413 = mrSges(6,1) * t456 - mrSges(6,3) * t422;
t481 = -t361 * t463 + t467 * t362;
t353 = m(6) * t370 - mrSges(6,2) * t448 + mrSges(6,3) * t397 + t405 * t421 - t413 * t456 + t481;
t348 = t468 * t352 + t464 * t353;
t423 = -mrSges(5,1) * t446 + mrSges(5,2) * t447;
t429 = -mrSges(5,2) * t457 + mrSges(5,3) * t446;
t346 = m(5) * t391 + mrSges(5,1) * t451 - mrSges(5,3) * t432 - t423 * t447 + t429 * t457 + t348;
t430 = mrSges(5,1) * t457 - mrSges(5,3) * t447;
t482 = -t352 * t464 + t468 * t353;
t347 = m(5) * t392 - mrSges(5,2) * t451 + mrSges(5,3) * t431 + t423 * t446 - t430 * t457 + t482;
t340 = t462 * t346 + t461 * t347;
t483 = -t346 * t461 + t462 * t347;
t424 = g(3) * t465 + t434 * t469;
t410 = -qJDD(3) * pkin(3) - qJ(4) * t471 + t449 * t489 + qJDD(4) - t424;
t393 = -pkin(4) * t431 - pkin(8) * t445 + t447 * t433 + t410;
t372 = -pkin(5) * t397 - pkin(9) * t420 + t414 * t422 + t393;
t479 = m(7) * t372 - t377 * mrSges(7,1) + t378 * mrSges(7,2) - t403 * t394 + t404 * t395;
t474 = m(6) * t393 - t397 * mrSges(6,1) + t398 * mrSges(6,2) - t421 * t412 + t422 * t413 + t479;
t363 = m(5) * t410 - t431 * mrSges(5,1) + t432 * mrSges(5,2) - t446 * t429 + t447 * t430 + t474;
t450 = (t465 * mrSges(4,1) + t469 * mrSges(4,2)) * qJD(1);
t453 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t457;
t454 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t489;
t478 = (m(4) * t425 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t451 - qJD(3) * t454 - t450 * t457 + t483) * t465 + (m(4) * t424 + qJDD(3) * mrSges(4,1) - t452 * mrSges(4,3) + qJD(3) * t453 - t450 * t489 - t363) * t469;
t386 = Ifges(7,4) * t404 + Ifges(7,2) * t403 + Ifges(7,6) * t455;
t387 = Ifges(7,1) * t404 + Ifges(7,4) * t403 + Ifges(7,5) * t455;
t475 = -mrSges(7,1) * t365 + mrSges(7,2) * t366 - Ifges(7,5) * t378 - Ifges(7,6) * t377 - Ifges(7,3) * t444 - t404 * t386 + t403 * t387;
t400 = Ifges(6,4) * t422 + Ifges(6,2) * t421 + Ifges(6,6) * t456;
t401 = Ifges(6,1) * t422 + Ifges(6,4) * t421 + Ifges(6,5) * t456;
t473 = mrSges(6,1) * t369 - mrSges(6,2) * t370 + Ifges(6,5) * t398 + Ifges(6,6) * t397 + Ifges(6,3) * t448 + pkin(5) * t354 + t422 * t400 - t421 * t401 - t475;
t443 = (Ifges(4,5) * qJD(3)) + (t469 * Ifges(4,1) - t465 * Ifges(4,4)) * qJD(1);
t442 = (Ifges(4,6) * qJD(3)) + (t469 * Ifges(4,4) - Ifges(4,2) * t465) * qJD(1);
t436 = -qJDD(1) * pkin(1) + t476;
t435 = pkin(1) * t472 + t491;
t417 = Ifges(5,1) * t447 + Ifges(5,4) * t446 + Ifges(5,5) * t457;
t416 = Ifges(5,4) * t447 + Ifges(5,2) * t446 + Ifges(5,6) * t457;
t415 = Ifges(5,5) * t447 + Ifges(5,6) * t446 + Ifges(5,3) * t457;
t399 = Ifges(6,5) * t422 + Ifges(6,6) * t421 + Ifges(6,3) * t456;
t385 = Ifges(7,5) * t404 + Ifges(7,6) * t403 + Ifges(7,3) * t455;
t356 = mrSges(7,2) * t372 - mrSges(7,3) * t365 + Ifges(7,1) * t378 + Ifges(7,4) * t377 + Ifges(7,5) * t444 + t385 * t403 - t386 * t455;
t355 = -mrSges(7,1) * t372 + mrSges(7,3) * t366 + Ifges(7,4) * t378 + Ifges(7,2) * t377 + Ifges(7,6) * t444 - t385 * t404 + t387 * t455;
t342 = mrSges(6,2) * t393 - mrSges(6,3) * t369 + Ifges(6,1) * t398 + Ifges(6,4) * t397 + Ifges(6,5) * t448 - pkin(9) * t354 - t355 * t463 + t356 * t467 + t399 * t421 - t400 * t456;
t341 = -mrSges(6,1) * t393 + mrSges(6,3) * t370 + Ifges(6,4) * t398 + Ifges(6,2) * t397 + Ifges(6,6) * t448 - pkin(5) * t479 + pkin(9) * t481 + t467 * t355 + t463 * t356 - t422 * t399 + t456 * t401;
t338 = m(3) * t436 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t472) + t478;
t337 = mrSges(5,2) * t410 - mrSges(5,3) * t391 + Ifges(5,1) * t432 + Ifges(5,4) * t431 + Ifges(5,5) * t451 - pkin(8) * t348 - t341 * t464 + t342 * t468 + t415 * t446 - t416 * t457;
t336 = -mrSges(5,1) * t410 + mrSges(5,3) * t392 + Ifges(5,4) * t432 + Ifges(5,2) * t431 + Ifges(5,6) * t451 - pkin(4) * t474 + pkin(8) * t482 + t468 * t341 + t464 * t342 - t447 * t415 + t417 * t457;
t1 = [mrSges(2,1) * t484 - mrSges(2,2) * t480 + mrSges(3,2) * t436 - mrSges(3,3) * t435 + t469 * (mrSges(4,2) * t428 - mrSges(4,3) * t424 + Ifges(4,1) * t452 - Ifges(4,4) * t451 + Ifges(4,5) * qJDD(3) - qJ(4) * t340 - qJD(3) * t442 - t461 * t336 + t462 * t337) - t465 * (-mrSges(4,1) * t428 - mrSges(5,1) * t391 + mrSges(5,2) * t392 + mrSges(4,3) * t425 + Ifges(4,4) * t452 - Ifges(5,5) * t432 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t431 - pkin(3) * t340 - pkin(4) * t348 + qJD(3) * t443 - t447 * t416 + t446 * t417 - t473 + (-Ifges(5,3) - Ifges(4,2)) * t451) - pkin(7) * t478 - pkin(1) * t338 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t435 + m(4) * t428 + mrSges(4,1) * t451 + mrSges(3,2) * t472 + mrSges(4,2) * t452 + t340 + qJDD(1) * mrSges(3,3) + (t453 * t465 + t454 * t469) * qJD(1)) * qJ(2); t338; Ifges(4,5) * t452 - Ifges(4,6) * t451 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t424 - mrSges(4,2) * t425 + t461 * t337 + t462 * t336 - pkin(3) * t363 + qJ(4) * t483 + (t469 * t442 + t465 * t443) * qJD(1); t363; t473; -t475;];
tauJ  = t1;
