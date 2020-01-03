% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:35
% EndTime: 2019-12-31 18:21:39
% DurationCPUTime: 3.85s
% Computational Cost: add. (45674->268), mult. (93332->339), div. (0->0), fcn. (57621->10), ass. (0->106)
t480 = sin(qJ(1));
t483 = cos(qJ(1));
t466 = t480 * g(1) - t483 * g(2);
t458 = qJDD(1) * pkin(1) + t466;
t467 = -t483 * g(1) - t480 * g(2);
t485 = qJD(1) ^ 2;
t461 = -t485 * pkin(1) + t467;
t475 = sin(pkin(8));
t477 = cos(pkin(8));
t435 = t475 * t458 + t477 * t461;
t427 = -t485 * pkin(2) + qJDD(1) * pkin(6) + t435;
t473 = -g(3) + qJDD(2);
t479 = sin(qJ(3));
t482 = cos(qJ(3));
t421 = t482 * t427 + t479 * t473;
t460 = (-mrSges(4,1) * t482 + mrSges(4,2) * t479) * qJD(1);
t496 = qJD(1) * qJD(3);
t470 = t479 * t496;
t463 = t482 * qJDD(1) - t470;
t498 = qJD(1) * t479;
t464 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t498;
t434 = t477 * t458 - t475 * t461;
t426 = -qJDD(1) * pkin(2) - t485 * pkin(6) - t434;
t494 = t482 * t496;
t462 = t479 * qJDD(1) + t494;
t415 = (-t462 - t494) * qJ(4) + (-t463 + t470) * pkin(3) + t426;
t459 = (-pkin(3) * t482 - qJ(4) * t479) * qJD(1);
t484 = qJD(3) ^ 2;
t497 = t482 * qJD(1);
t419 = -t484 * pkin(3) + qJDD(3) * qJ(4) + t459 * t497 + t421;
t474 = sin(pkin(9));
t476 = cos(pkin(9));
t455 = t474 * qJD(3) + t476 * t498;
t403 = -0.2e1 * qJD(4) * t455 + t476 * t415 - t474 * t419;
t441 = t474 * qJDD(3) + t476 * t462;
t454 = t476 * qJD(3) - t474 * t498;
t401 = (-t454 * t497 - t441) * pkin(7) + (t454 * t455 - t463) * pkin(4) + t403;
t404 = 0.2e1 * qJD(4) * t454 + t474 * t415 + t476 * t419;
t440 = t476 * qJDD(3) - t474 * t462;
t442 = -pkin(4) * t497 - t455 * pkin(7);
t453 = t454 ^ 2;
t402 = -t453 * pkin(4) + t440 * pkin(7) + t442 * t497 + t404;
t478 = sin(qJ(5));
t481 = cos(qJ(5));
t399 = t481 * t401 - t478 * t402;
t432 = t481 * t454 - t478 * t455;
t408 = t432 * qJD(5) + t478 * t440 + t481 * t441;
t433 = t478 * t454 + t481 * t455;
t418 = -t432 * mrSges(6,1) + t433 * mrSges(6,2);
t468 = qJD(5) - t497;
t422 = -t468 * mrSges(6,2) + t432 * mrSges(6,3);
t457 = qJDD(5) - t463;
t397 = m(6) * t399 + t457 * mrSges(6,1) - t408 * mrSges(6,3) - t433 * t418 + t468 * t422;
t400 = t478 * t401 + t481 * t402;
t407 = -t433 * qJD(5) + t481 * t440 - t478 * t441;
t423 = t468 * mrSges(6,1) - t433 * mrSges(6,3);
t398 = m(6) * t400 - t457 * mrSges(6,2) + t407 * mrSges(6,3) + t432 * t418 - t468 * t423;
t389 = t481 * t397 + t478 * t398;
t436 = -t454 * mrSges(5,1) + t455 * mrSges(5,2);
t438 = mrSges(5,2) * t497 + t454 * mrSges(5,3);
t387 = m(5) * t403 - t463 * mrSges(5,1) - t441 * mrSges(5,3) - t455 * t436 - t438 * t497 + t389;
t439 = -mrSges(5,1) * t497 - t455 * mrSges(5,3);
t489 = -t478 * t397 + t481 * t398;
t388 = m(5) * t404 + t463 * mrSges(5,2) + t440 * mrSges(5,3) + t454 * t436 + t439 * t497 + t489;
t490 = -t474 * t387 + t476 * t388;
t384 = m(4) * t421 - qJDD(3) * mrSges(4,2) + t463 * mrSges(4,3) - qJD(3) * t464 + t460 * t497 + t490;
t420 = -t479 * t427 + t482 * t473;
t465 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t497;
t417 = -qJDD(3) * pkin(3) - t484 * qJ(4) + t459 * t498 + qJDD(4) - t420;
t405 = -t440 * pkin(4) - t453 * pkin(7) + t455 * t442 + t417;
t488 = m(6) * t405 - t407 * mrSges(6,1) + t408 * mrSges(6,2) - t432 * t422 + t433 * t423;
t486 = -m(5) * t417 + t440 * mrSges(5,1) - t441 * mrSges(5,2) + t454 * t438 - t455 * t439 - t488;
t393 = m(4) * t420 + qJDD(3) * mrSges(4,1) - t462 * mrSges(4,3) + qJD(3) * t465 - t460 * t498 + t486;
t491 = t482 * t384 - t479 * t393;
t376 = m(3) * t435 - t485 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t491;
t385 = t476 * t387 + t474 * t388;
t487 = -m(4) * t426 + t463 * mrSges(4,1) - t462 * mrSges(4,2) - t464 * t498 + t465 * t497 - t385;
t381 = m(3) * t434 + qJDD(1) * mrSges(3,1) - t485 * mrSges(3,2) + t487;
t372 = t475 * t376 + t477 * t381;
t370 = m(2) * t466 + qJDD(1) * mrSges(2,1) - t485 * mrSges(2,2) + t372;
t492 = t477 * t376 - t475 * t381;
t371 = m(2) * t467 - t485 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t492;
t499 = t483 * t370 + t480 * t371;
t378 = t479 * t384 + t482 * t393;
t495 = m(3) * t473 + t378;
t493 = -t480 * t370 + t483 * t371;
t451 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t479 + Ifges(4,4) * t482) * qJD(1);
t450 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t479 + Ifges(4,2) * t482) * qJD(1);
t449 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t479 + Ifges(4,6) * t482) * qJD(1);
t430 = Ifges(5,1) * t455 + Ifges(5,4) * t454 - Ifges(5,5) * t497;
t429 = Ifges(5,4) * t455 + Ifges(5,2) * t454 - Ifges(5,6) * t497;
t428 = Ifges(5,5) * t455 + Ifges(5,6) * t454 - Ifges(5,3) * t497;
t413 = Ifges(6,1) * t433 + Ifges(6,4) * t432 + Ifges(6,5) * t468;
t412 = Ifges(6,4) * t433 + Ifges(6,2) * t432 + Ifges(6,6) * t468;
t411 = Ifges(6,5) * t433 + Ifges(6,6) * t432 + Ifges(6,3) * t468;
t391 = mrSges(6,2) * t405 - mrSges(6,3) * t399 + Ifges(6,1) * t408 + Ifges(6,4) * t407 + Ifges(6,5) * t457 + t432 * t411 - t468 * t412;
t390 = -mrSges(6,1) * t405 + mrSges(6,3) * t400 + Ifges(6,4) * t408 + Ifges(6,2) * t407 + Ifges(6,6) * t457 - t433 * t411 + t468 * t413;
t379 = mrSges(5,2) * t417 - mrSges(5,3) * t403 + Ifges(5,1) * t441 + Ifges(5,4) * t440 - Ifges(5,5) * t463 - pkin(7) * t389 - t478 * t390 + t481 * t391 + t454 * t428 + t429 * t497;
t377 = -mrSges(5,1) * t417 + mrSges(5,3) * t404 + Ifges(5,4) * t441 + Ifges(5,2) * t440 - Ifges(5,6) * t463 - pkin(4) * t488 + pkin(7) * t489 + t481 * t390 + t478 * t391 - t455 * t428 - t430 * t497;
t373 = Ifges(4,4) * t462 + Ifges(4,6) * qJDD(3) - t449 * t498 + qJD(3) * t451 - mrSges(4,1) * t426 + mrSges(4,3) * t421 - Ifges(5,5) * t441 - Ifges(5,6) * t440 - t455 * t429 + t454 * t430 - mrSges(5,1) * t403 + mrSges(5,2) * t404 - Ifges(6,5) * t408 - Ifges(6,6) * t407 - Ifges(6,3) * t457 - t433 * t412 + t432 * t413 - mrSges(6,1) * t399 + mrSges(6,2) * t400 - pkin(4) * t389 - pkin(3) * t385 + (Ifges(4,2) + Ifges(5,3)) * t463;
t366 = mrSges(4,2) * t426 - mrSges(4,3) * t420 + Ifges(4,1) * t462 + Ifges(4,4) * t463 + Ifges(4,5) * qJDD(3) - qJ(4) * t385 - qJD(3) * t450 - t474 * t377 + t476 * t379 + t449 * t497;
t365 = Ifges(3,6) * qJDD(1) + t485 * Ifges(3,5) - mrSges(3,1) * t473 + mrSges(3,3) * t435 - Ifges(4,5) * t462 - Ifges(4,6) * t463 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t420 + mrSges(4,2) * t421 - t474 * t379 - t476 * t377 - pkin(3) * t486 - qJ(4) * t490 - pkin(2) * t378 + (-t479 * t450 + t482 * t451) * qJD(1);
t364 = mrSges(3,2) * t473 - mrSges(3,3) * t434 + Ifges(3,5) * qJDD(1) - t485 * Ifges(3,6) - pkin(6) * t378 + t482 * t366 - t479 * t373;
t363 = -mrSges(2,2) * g(3) - mrSges(2,3) * t466 + Ifges(2,5) * qJDD(1) - t485 * Ifges(2,6) - qJ(2) * t372 + t477 * t364 - t475 * t365;
t362 = mrSges(2,1) * g(3) + mrSges(2,3) * t467 + t485 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t495 + qJ(2) * t492 + t475 * t364 + t477 * t365;
t1 = [-m(1) * g(1) + t493; -m(1) * g(2) + t499; (-m(1) - m(2)) * g(3) + t495; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t499 - t480 * t362 + t483 * t363; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t493 + t483 * t362 + t480 * t363; pkin(1) * t372 + mrSges(2,1) * t466 - mrSges(2,2) * t467 + t482 * t373 + pkin(2) * t487 + pkin(6) * t491 + t479 * t366 + mrSges(3,1) * t434 - mrSges(3,2) * t435 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
