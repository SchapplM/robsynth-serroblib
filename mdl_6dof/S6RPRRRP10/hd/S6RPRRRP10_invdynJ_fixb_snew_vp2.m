% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP10
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
% Datum: 2019-05-06 01:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:55:52
% EndTime: 2019-05-06 01:55:57
% DurationCPUTime: 2.44s
% Computational Cost: add. (18155->266), mult. (35095->320), div. (0->0), fcn. (22445->8), ass. (0->108)
t502 = Ifges(6,1) + Ifges(7,1);
t491 = Ifges(6,4) - Ifges(7,5);
t500 = Ifges(7,4) + Ifges(6,5);
t501 = Ifges(6,2) + Ifges(7,3);
t499 = Ifges(6,6) - Ifges(7,6);
t498 = -Ifges(6,3) - Ifges(7,2);
t464 = sin(qJ(4));
t467 = cos(qJ(4));
t468 = cos(qJ(3));
t488 = t468 * qJD(1);
t446 = qJD(3) * t467 - t464 * t488;
t447 = qJD(3) * t464 + t467 * t488;
t463 = sin(qJ(5));
t493 = cos(qJ(5));
t421 = -t493 * t446 + t447 * t463;
t422 = t463 * t446 + t493 * t447;
t465 = sin(qJ(3));
t459 = t465 * qJD(1);
t456 = t459 + qJD(4);
t455 = qJD(5) + t456;
t497 = t501 * t421 - t491 * t422 - t499 * t455;
t496 = -t491 * t421 + t502 * t422 + t500 * t455;
t466 = sin(qJ(1));
t469 = cos(qJ(1));
t479 = -g(1) * t469 - g(2) * t466;
t495 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t479;
t494 = -pkin(1) - pkin(7);
t492 = -mrSges(6,3) - mrSges(7,2);
t471 = qJD(1) ^ 2;
t430 = t494 * t471 - t495;
t487 = qJD(1) * qJD(3);
t457 = t468 * t487;
t450 = -t465 * qJDD(1) - t457;
t484 = t465 * t487;
t451 = qJDD(1) * t468 - t484;
t402 = (-t451 + t484) * pkin(8) + (-t450 + t457) * pkin(3) + t430;
t483 = g(1) * t466 - t469 * g(2);
t476 = -qJ(2) * t471 + qJDD(2) - t483;
t431 = t494 * qJDD(1) + t476;
t425 = -g(3) * t468 + t465 * t431;
t449 = (t465 * pkin(3) - t468 * pkin(8)) * qJD(1);
t470 = qJD(3) ^ 2;
t407 = -pkin(3) * t470 + qJDD(3) * pkin(8) - t449 * t459 + t425;
t375 = t467 * t402 - t407 * t464;
t420 = qJD(4) * t446 + qJDD(3) * t464 + t451 * t467;
t445 = qJDD(4) - t450;
t371 = (t446 * t456 - t420) * pkin(9) + (t446 * t447 + t445) * pkin(4) + t375;
t376 = t464 * t402 + t467 * t407;
t419 = -qJD(4) * t447 + qJDD(3) * t467 - t451 * t464;
t429 = pkin(4) * t456 - pkin(9) * t447;
t444 = t446 ^ 2;
t373 = -pkin(4) * t444 + pkin(9) * t419 - t429 * t456 + t376;
t369 = t463 * t371 + t493 * t373;
t386 = qJD(5) * t422 - t493 * t419 + t420 * t463;
t410 = mrSges(6,1) * t455 - mrSges(6,3) * t422;
t443 = qJDD(5) + t445;
t397 = pkin(5) * t421 - qJ(6) * t422;
t454 = t455 ^ 2;
t363 = -pkin(5) * t454 + qJ(6) * t443 + 0.2e1 * qJD(6) * t455 - t397 * t421 + t369;
t411 = -mrSges(7,1) * t455 + mrSges(7,2) * t422;
t485 = m(7) * t363 + t443 * mrSges(7,3) + t455 * t411;
t398 = mrSges(7,1) * t421 - mrSges(7,3) * t422;
t489 = -mrSges(6,1) * t421 - mrSges(6,2) * t422 - t398;
t353 = m(6) * t369 - mrSges(6,2) * t443 + t492 * t386 - t410 * t455 + t489 * t421 + t485;
t368 = t493 * t371 - t463 * t373;
t387 = -t421 * qJD(5) + t463 * t419 + t493 * t420;
t409 = -mrSges(6,2) * t455 - mrSges(6,3) * t421;
t364 = -t443 * pkin(5) - t454 * qJ(6) + t422 * t397 + qJDD(6) - t368;
t408 = -mrSges(7,2) * t421 + mrSges(7,3) * t455;
t480 = -m(7) * t364 + t443 * mrSges(7,1) + t455 * t408;
t355 = m(6) * t368 + mrSges(6,1) * t443 + t492 * t387 + t409 * t455 + t489 * t422 + t480;
t349 = t463 * t353 + t493 * t355;
t423 = -mrSges(5,1) * t446 + mrSges(5,2) * t447;
t426 = -mrSges(5,2) * t456 + mrSges(5,3) * t446;
t345 = m(5) * t375 + mrSges(5,1) * t445 - mrSges(5,3) * t420 - t423 * t447 + t426 * t456 + t349;
t427 = mrSges(5,1) * t456 - mrSges(5,3) * t447;
t481 = t493 * t353 - t355 * t463;
t346 = m(5) * t376 - mrSges(5,2) * t445 + mrSges(5,3) * t419 + t423 * t446 - t427 * t456 + t481;
t341 = t467 * t345 + t464 * t346;
t490 = t499 * t421 - t500 * t422 + t498 * t455;
t482 = -t345 * t464 + t467 * t346;
t424 = g(3) * t465 + t431 * t468;
t406 = -qJDD(3) * pkin(3) - pkin(8) * t470 + t449 * t488 - t424;
t374 = -pkin(4) * t419 - pkin(9) * t444 + t447 * t429 + t406;
t366 = -0.2e1 * qJD(6) * t422 + (t421 * t455 - t387) * qJ(6) + (t422 * t455 + t386) * pkin(5) + t374;
t356 = m(7) * t366 + t386 * mrSges(7,1) - t387 * mrSges(7,3) + t421 * t408 - t422 * t411;
t448 = (t465 * mrSges(4,1) + t468 * mrSges(4,2)) * qJD(1);
t452 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t459;
t453 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t488;
t475 = m(6) * t374 + t386 * mrSges(6,1) + t387 * mrSges(6,2) + t421 * t409 + t422 * t410 + t356;
t473 = -m(5) * t406 + t419 * mrSges(5,1) - t420 * mrSges(5,2) + t446 * t426 - t447 * t427 - t475;
t478 = (m(4) * t425 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t450 - qJD(3) * t453 - t448 * t459 + t482) * t465 + (m(4) * t424 + qJDD(3) * mrSges(4,1) - t451 * mrSges(4,3) + qJD(3) * t452 - t448 * t488 + t473) * t468;
t360 = mrSges(7,2) * t387 + t398 * t422 - t480;
t474 = -mrSges(6,1) * t368 + mrSges(7,1) * t364 + mrSges(6,2) * t369 - mrSges(7,3) * t363 + pkin(5) * t360 - qJ(6) * t485 + t498 * t443 + t497 * t422 + (qJ(6) * t398 - t496) * t421 - t500 * t387 + (qJ(6) * mrSges(7,2) + t499) * t386;
t414 = Ifges(5,4) * t447 + Ifges(5,2) * t446 + Ifges(5,6) * t456;
t415 = Ifges(5,1) * t447 + Ifges(5,4) * t446 + Ifges(5,5) * t456;
t472 = mrSges(5,1) * t375 - mrSges(5,2) * t376 + Ifges(5,5) * t420 + Ifges(5,6) * t419 + Ifges(5,3) * t445 + pkin(4) * t349 + t447 * t414 - t446 * t415 - t474;
t442 = Ifges(4,5) * qJD(3) + (t468 * Ifges(4,1) - t465 * Ifges(4,4)) * qJD(1);
t441 = Ifges(4,6) * qJD(3) + (t468 * Ifges(4,4) - t465 * Ifges(4,2)) * qJD(1);
t433 = -qJDD(1) * pkin(1) + t476;
t432 = pkin(1) * t471 + t495;
t413 = Ifges(5,5) * t447 + Ifges(5,6) * t446 + Ifges(5,3) * t456;
t348 = mrSges(6,2) * t374 + mrSges(7,2) * t364 - mrSges(6,3) * t368 - mrSges(7,3) * t366 - qJ(6) * t356 - t491 * t386 + t502 * t387 + t490 * t421 + t500 * t443 + t497 * t455;
t347 = -mrSges(6,1) * t374 - mrSges(7,1) * t366 + mrSges(7,2) * t363 + mrSges(6,3) * t369 - pkin(5) * t356 - t501 * t386 + t491 * t387 + t490 * t422 + t499 * t443 + t496 * t455;
t339 = m(3) * t433 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t471 + t478;
t338 = mrSges(5,2) * t406 - mrSges(5,3) * t375 + Ifges(5,1) * t420 + Ifges(5,4) * t419 + Ifges(5,5) * t445 - pkin(9) * t349 - t463 * t347 + t493 * t348 + t446 * t413 - t456 * t414;
t337 = -mrSges(5,1) * t406 + mrSges(5,3) * t376 + Ifges(5,4) * t420 + Ifges(5,2) * t419 + Ifges(5,6) * t445 - pkin(4) * t475 + pkin(9) * t481 + t493 * t347 + t463 * t348 - t447 * t413 + t456 * t415;
t1 = [mrSges(2,1) * t483 - mrSges(2,2) * t479 + mrSges(3,2) * t433 - mrSges(3,3) * t432 + t468 * (mrSges(4,2) * t430 - mrSges(4,3) * t424 + Ifges(4,1) * t451 + Ifges(4,4) * t450 + Ifges(4,5) * qJDD(3) - pkin(8) * t341 - qJD(3) * t441 - t337 * t464 + t467 * t338) - t465 * (-mrSges(4,1) * t430 + mrSges(4,3) * t425 + Ifges(4,4) * t451 + Ifges(4,2) * t450 + Ifges(4,6) * qJDD(3) - pkin(3) * t341 + qJD(3) * t442 - t472) - pkin(7) * t478 - pkin(1) * t339 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t432 + m(4) * t430 - mrSges(4,1) * t450 + mrSges(3,2) * t471 + mrSges(4,2) * t451 + t341 + qJDD(1) * mrSges(3,3) + (t452 * t465 + t453 * t468) * qJD(1)) * qJ(2); t339; Ifges(4,5) * t451 + Ifges(4,6) * t450 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t424 - mrSges(4,2) * t425 + t464 * t338 + t467 * t337 + pkin(3) * t473 + pkin(8) * t482 + (t468 * t441 + t465 * t442) * qJD(1); t472; -t474; t360;];
tauJ  = t1;
