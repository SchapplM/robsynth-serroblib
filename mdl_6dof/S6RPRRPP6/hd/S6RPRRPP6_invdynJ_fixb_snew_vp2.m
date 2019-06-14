% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-05-05 21:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:41:41
% EndTime: 2019-05-05 21:41:46
% DurationCPUTime: 2.27s
% Computational Cost: add. (16684->265), mult. (33264->321), div. (0->0), fcn. (21033->8), ass. (0->107)
t493 = Ifges(6,1) + Ifges(7,1);
t482 = Ifges(6,4) - Ifges(7,5);
t481 = Ifges(6,5) + Ifges(7,4);
t492 = -Ifges(6,2) - Ifges(7,3);
t480 = Ifges(6,6) - Ifges(7,6);
t491 = -Ifges(7,2) - Ifges(6,3);
t455 = qJD(1) ^ 2;
t484 = -pkin(1) - pkin(7);
t450 = sin(qJ(1));
t453 = cos(qJ(1));
t462 = -g(1) * t453 - g(2) * t450;
t487 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t462;
t421 = t484 * t455 - t487;
t449 = sin(qJ(3));
t452 = cos(qJ(3));
t472 = qJD(1) * qJD(3);
t467 = t452 * t472;
t438 = -qJDD(1) * t449 - t467;
t468 = t449 * t472;
t439 = qJDD(1) * t452 - t468;
t392 = (-t439 + t468) * pkin(8) + (-t438 + t467) * pkin(3) + t421;
t466 = g(1) * t450 - t453 * g(2);
t458 = -qJ(2) * t455 + qJDD(2) - t466;
t422 = t484 * qJDD(1) + t458;
t416 = -g(3) * t452 + t449 * t422;
t437 = (pkin(3) * t449 - pkin(8) * t452) * qJD(1);
t454 = qJD(3) ^ 2;
t474 = qJD(1) * t449;
t397 = -pkin(3) * t454 + qJDD(3) * pkin(8) - t437 * t474 + t416;
t448 = sin(qJ(4));
t451 = cos(qJ(4));
t371 = t451 * t392 - t397 * t448;
t473 = qJD(1) * t452;
t434 = qJD(3) * t451 - t448 * t473;
t411 = qJD(4) * t434 + qJDD(3) * t448 + t439 * t451;
t433 = qJDD(4) - t438;
t435 = qJD(3) * t448 + t451 * t473;
t443 = qJD(4) + t474;
t367 = (t434 * t443 - t411) * qJ(5) + (t434 * t435 + t433) * pkin(4) + t371;
t372 = t448 * t392 + t451 * t397;
t410 = -qJD(4) * t435 + qJDD(3) * t451 - t439 * t448;
t418 = pkin(4) * t443 - qJ(5) * t435;
t432 = t434 ^ 2;
t369 = -pkin(4) * t432 + qJ(5) * t410 - t418 * t443 + t372;
t447 = sin(pkin(9));
t479 = cos(pkin(9));
t412 = -t479 * t434 + t435 * t447;
t485 = -2 * qJD(5);
t363 = t447 * t367 + t479 * t369 + t412 * t485;
t382 = -t479 * t410 + t411 * t447;
t413 = t447 * t434 + t479 * t435;
t400 = mrSges(6,1) * t443 - mrSges(6,3) * t413;
t387 = pkin(5) * t412 - qJ(6) * t413;
t442 = t443 ^ 2;
t360 = -pkin(5) * t442 + qJ(6) * t433 + 0.2e1 * qJD(6) * t443 - t387 * t412 + t363;
t401 = -mrSges(7,1) * t443 + mrSges(7,2) * t413;
t469 = m(7) * t360 + t433 * mrSges(7,3) + t443 * t401;
t388 = mrSges(7,1) * t412 - mrSges(7,3) * t413;
t475 = -mrSges(6,1) * t412 - mrSges(6,2) * t413 - t388;
t483 = -mrSges(6,3) - mrSges(7,2);
t351 = m(6) * t363 - mrSges(6,2) * t433 + t483 * t382 - t400 * t443 + t475 * t412 + t469;
t459 = t479 * t367 - t447 * t369;
t362 = t413 * t485 + t459;
t383 = t447 * t410 + t479 * t411;
t398 = -mrSges(6,2) * t443 - mrSges(6,3) * t412;
t361 = -t433 * pkin(5) - t442 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t387) * t413 - t459;
t399 = -mrSges(7,2) * t412 + mrSges(7,3) * t443;
t463 = -m(7) * t361 + t433 * mrSges(7,1) + t443 * t399;
t353 = m(6) * t362 + mrSges(6,1) * t433 + t483 * t383 + t398 * t443 + t475 * t413 + t463;
t346 = t447 * t351 + t479 * t353;
t357 = mrSges(7,2) * t383 + t388 * t413 - t463;
t404 = Ifges(5,4) * t435 + Ifges(5,2) * t434 + Ifges(5,6) * t443;
t405 = Ifges(5,1) * t435 + Ifges(5,4) * t434 + Ifges(5,5) * t443;
t476 = -t482 * t412 + t413 * t493 + t481 * t443;
t477 = t412 * t492 + t413 * t482 + t443 * t480;
t490 = -t480 * t382 + t481 * t383 + (Ifges(5,3) - t491) * t433 + mrSges(5,1) * t371 + mrSges(6,1) * t362 - mrSges(7,1) * t361 - mrSges(5,2) * t372 - mrSges(6,2) * t363 + mrSges(7,3) * t360 + Ifges(5,5) * t411 + Ifges(5,6) * t410 + pkin(4) * t346 - pkin(5) * t357 + qJ(6) * (-mrSges(7,2) * t382 - t388 * t412 + t469) + t435 * t404 - t434 * t405 + t477 * t413 + t476 * t412;
t414 = -mrSges(5,1) * t434 + mrSges(5,2) * t435;
t417 = -mrSges(5,2) * t443 + mrSges(5,3) * t434;
t344 = m(5) * t371 + mrSges(5,1) * t433 - mrSges(5,3) * t411 - t414 * t435 + t417 * t443 + t346;
t419 = mrSges(5,1) * t443 - mrSges(5,3) * t435;
t464 = t479 * t351 - t353 * t447;
t345 = m(5) * t372 - mrSges(5,2) * t433 + mrSges(5,3) * t410 + t414 * t434 - t419 * t443 + t464;
t340 = t451 * t344 + t448 * t345;
t478 = t480 * t412 - t481 * t413 + t491 * t443;
t465 = -t344 * t448 + t451 * t345;
t415 = g(3) * t449 + t422 * t452;
t396 = -qJDD(3) * pkin(3) - pkin(8) * t454 + t437 * t473 - t415;
t370 = -pkin(4) * t410 - qJ(5) * t432 + t435 * t418 + qJDD(5) + t396;
t365 = -0.2e1 * qJD(6) * t413 + (t412 * t443 - t383) * qJ(6) + (t413 * t443 + t382) * pkin(5) + t370;
t358 = m(7) * t365 + t382 * mrSges(7,1) - t383 * mrSges(7,3) + t412 * t399 - t413 * t401;
t436 = (mrSges(4,1) * t449 + mrSges(4,2) * t452) * qJD(1);
t440 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t474;
t441 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t473;
t355 = m(6) * t370 + t382 * mrSges(6,1) + t383 * mrSges(6,2) + t412 * t398 + t413 * t400 + t358;
t457 = -m(5) * t396 + t410 * mrSges(5,1) - t411 * mrSges(5,2) + t434 * t417 - t435 * t419 - t355;
t461 = (m(4) * t416 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t438 - qJD(3) * t441 - t436 * t474 + t465) * t449 + (m(4) * t415 + qJDD(3) * mrSges(4,1) - t439 * mrSges(4,3) + qJD(3) * t440 - t436 * t473 + t457) * t452;
t429 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t452 - Ifges(4,4) * t449) * qJD(1);
t428 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t452 - Ifges(4,2) * t449) * qJD(1);
t424 = -qJDD(1) * pkin(1) + t458;
t423 = pkin(1) * t455 + t487;
t403 = Ifges(5,5) * t435 + Ifges(5,6) * t434 + Ifges(5,3) * t443;
t348 = mrSges(6,2) * t370 + mrSges(7,2) * t361 - mrSges(6,3) * t362 - mrSges(7,3) * t365 - qJ(6) * t358 - t482 * t382 + t383 * t493 + t478 * t412 + t481 * t433 - t477 * t443;
t347 = -mrSges(6,1) * t370 - mrSges(7,1) * t365 + mrSges(7,2) * t360 + mrSges(6,3) * t363 - pkin(5) * t358 + t382 * t492 + t482 * t383 + t478 * t413 + t480 * t433 + t476 * t443;
t338 = m(3) * t424 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t455 + t461;
t337 = mrSges(5,2) * t396 - mrSges(5,3) * t371 + Ifges(5,1) * t411 + Ifges(5,4) * t410 + Ifges(5,5) * t433 - qJ(5) * t346 - t447 * t347 + t479 * t348 + t434 * t403 - t443 * t404;
t336 = -mrSges(5,1) * t396 + mrSges(5,3) * t372 + Ifges(5,4) * t411 + Ifges(5,2) * t410 + Ifges(5,6) * t433 - pkin(4) * t355 + qJ(5) * t464 + t479 * t347 + t447 * t348 - t435 * t403 + t443 * t405;
t1 = [mrSges(2,1) * t466 - mrSges(2,2) * t462 + mrSges(3,2) * t424 - mrSges(3,3) * t423 + t452 * (mrSges(4,2) * t421 - mrSges(4,3) * t415 + Ifges(4,1) * t439 + Ifges(4,4) * t438 + Ifges(4,5) * qJDD(3) - pkin(8) * t340 - qJD(3) * t428 - t336 * t448 + t337 * t451) - pkin(7) * t461 - pkin(1) * t338 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (mrSges(4,1) * t421 - mrSges(4,3) * t416 - Ifges(4,4) * t439 - Ifges(4,2) * t438 - Ifges(4,6) * qJDD(3) + pkin(3) * t340 - qJD(3) * t429 + t490) * t449 + (-m(3) * t423 + m(4) * t421 - mrSges(4,1) * t438 + mrSges(3,2) * t455 + mrSges(4,2) * t439 + t340 + qJDD(1) * mrSges(3,3) + (t440 * t449 + t441 * t452) * qJD(1)) * qJ(2); t338; Ifges(4,5) * t439 + Ifges(4,6) * t438 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t415 - mrSges(4,2) * t416 + t448 * t337 + t451 * t336 + pkin(3) * t457 + pkin(8) * t465 + (t428 * t452 + t429 * t449) * qJD(1); t490; t355; t357;];
tauJ  = t1;
