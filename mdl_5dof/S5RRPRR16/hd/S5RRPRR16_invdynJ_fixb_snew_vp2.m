% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR16_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:45:05
% EndTime: 2019-12-31 20:45:09
% DurationCPUTime: 2.34s
% Computational Cost: add. (15980->272), mult. (36158->342), div. (0->0), fcn. (25270->10), ass. (0->123)
t499 = -2 * qJD(3);
t498 = Ifges(3,1) + Ifges(4,2);
t490 = Ifges(3,4) + Ifges(4,6);
t489 = Ifges(3,5) - Ifges(4,4);
t497 = Ifges(3,2) + Ifges(4,3);
t488 = Ifges(3,6) - Ifges(4,5);
t496 = Ifges(3,3) + Ifges(4,1);
t449 = cos(pkin(5));
t443 = qJD(1) * t449 + qJD(2);
t452 = sin(qJ(2));
t448 = sin(pkin(5));
t479 = qJD(1) * t448;
t473 = t452 * t479;
t495 = (pkin(2) * t443 + t499) * t473;
t458 = qJD(1) ^ 2;
t453 = sin(qJ(1));
t457 = cos(qJ(1));
t468 = -g(1) * t457 - g(2) * t453;
t476 = qJDD(1) * t448;
t426 = -pkin(1) * t458 + pkin(7) * t476 + t468;
t456 = cos(qJ(2));
t485 = t448 * t452;
t471 = t453 * g(1) - g(2) * t457;
t493 = pkin(7) * t448;
t425 = qJDD(1) * pkin(1) + t458 * t493 + t471;
t487 = t425 * t449;
t391 = -g(3) * t485 + t456 * t426 + t452 * t487;
t427 = (-t456 * pkin(2) - t452 * qJ(3)) * t479;
t441 = t443 ^ 2;
t442 = qJDD(1) * t449 + qJDD(2);
t478 = qJD(1) * t456;
t472 = t448 * t478;
t369 = pkin(2) * t441 - t442 * qJ(3) - t427 * t472 + t443 * t499 - t391;
t494 = -pkin(2) - pkin(8);
t492 = g(3) * t449;
t491 = mrSges(3,1) - mrSges(4,2);
t486 = t448 ^ 2 * t458;
t484 = t448 * t456;
t430 = pkin(3) * t473 - pkin(8) * t443;
t431 = (qJD(2) * t478 + qJDD(1) * t452) * t448;
t432 = -qJD(2) * t473 + t456 * t476;
t474 = t456 ^ 2 * t486;
t362 = -pkin(3) * t474 - t492 - qJ(3) * t431 + t494 * t432 + (-t425 + (-qJ(3) * t443 * t456 - t430 * t452) * qJD(1)) * t448 + t495;
t480 = g(3) * t484 + t452 * t426;
t466 = -qJ(3) * t441 + t427 * t473 + qJDD(3) + t480;
t364 = pkin(3) * t431 + t494 * t442 + (-pkin(3) * t443 * t479 - pkin(8) * t452 * t486 - t487) * t456 + t466;
t451 = sin(qJ(4));
t455 = cos(qJ(4));
t358 = t455 * t362 + t451 * t364;
t413 = -t443 * t451 - t455 * t472;
t414 = t443 * t455 - t451 * t472;
t393 = -pkin(4) * t413 - pkin(9) * t414;
t420 = qJDD(4) + t431;
t436 = qJD(4) + t473;
t434 = t436 ^ 2;
t354 = -pkin(4) * t434 + pkin(9) * t420 + t393 * t413 + t358;
t361 = pkin(3) * t432 - pkin(8) * t474 + t443 * t430 - t369;
t388 = -qJD(4) * t414 - t432 * t455 - t442 * t451;
t389 = qJD(4) * t413 - t432 * t451 + t442 * t455;
t355 = (-t413 * t436 - t389) * pkin(9) + (t414 * t436 - t388) * pkin(4) + t361;
t450 = sin(qJ(5));
t454 = cos(qJ(5));
t351 = -t354 * t450 + t355 * t454;
t394 = -t414 * t450 + t436 * t454;
t367 = qJD(5) * t394 + t389 * t454 + t420 * t450;
t395 = t414 * t454 + t436 * t450;
t376 = -mrSges(6,1) * t394 + mrSges(6,2) * t395;
t412 = qJD(5) - t413;
t378 = -mrSges(6,2) * t412 + mrSges(6,3) * t394;
t386 = qJDD(5) - t388;
t348 = m(6) * t351 + mrSges(6,1) * t386 - mrSges(6,3) * t367 - t376 * t395 + t378 * t412;
t352 = t354 * t454 + t355 * t450;
t366 = -qJD(5) * t395 - t389 * t450 + t420 * t454;
t379 = mrSges(6,1) * t412 - mrSges(6,3) * t395;
t349 = m(6) * t352 - mrSges(6,2) * t386 + mrSges(6,3) * t366 + t376 * t394 - t379 * t412;
t339 = t454 * t348 + t450 * t349;
t483 = (t452 * t489 + t456 * t488) * t479 + t496 * t443;
t482 = (t452 * t490 + t456 * t497) * t479 + t488 * t443;
t481 = (t452 * t498 + t456 * t490) * t479 + t489 * t443;
t475 = t456 * t487;
t392 = -mrSges(5,1) * t413 + mrSges(5,2) * t414;
t397 = mrSges(5,1) * t436 - mrSges(5,3) * t414;
t469 = -t348 * t450 + t454 * t349;
t337 = m(5) * t358 - mrSges(5,2) * t420 + mrSges(5,3) * t388 + t392 * t413 - t397 * t436 + t469;
t357 = -t362 * t451 + t364 * t455;
t396 = -mrSges(5,2) * t436 + mrSges(5,3) * t413;
t353 = -pkin(4) * t420 - pkin(9) * t434 + t393 * t414 - t357;
t464 = -m(6) * t353 + t366 * mrSges(6,1) - mrSges(6,2) * t367 + t394 * t378 - t379 * t395;
t344 = m(5) * t357 + mrSges(5,1) * t420 - mrSges(5,3) * t389 - t392 * t414 + t396 * t436 + t464;
t470 = t455 * t337 - t451 * t344;
t404 = -t425 * t448 - t492;
t334 = t337 * t451 + t344 * t455;
t370 = -pkin(2) * t432 + (-t443 * t472 - t431) * qJ(3) + t404 + t495;
t423 = -mrSges(4,1) * t472 - mrSges(4,3) * t443;
t467 = -m(4) * t370 + t431 * mrSges(4,3) - t423 * t472 - t470;
t375 = -pkin(2) * t442 + t466 - t475;
t465 = -m(4) * t375 - t431 * mrSges(4,1) - t334;
t462 = -m(5) * t361 + mrSges(5,1) * t388 - t389 * mrSges(5,2) + t396 * t413 - t414 * t397 - t339;
t371 = Ifges(6,5) * t395 + Ifges(6,6) * t394 + Ifges(6,3) * t412;
t373 = Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t412;
t342 = -mrSges(6,1) * t353 + mrSges(6,3) * t352 + Ifges(6,4) * t367 + Ifges(6,2) * t366 + Ifges(6,6) * t386 - t371 * t395 + t373 * t412;
t372 = Ifges(6,4) * t395 + Ifges(6,2) * t394 + Ifges(6,6) * t412;
t343 = mrSges(6,2) * t353 - mrSges(6,3) * t351 + Ifges(6,1) * t367 + Ifges(6,4) * t366 + Ifges(6,5) * t386 + t371 * t394 - t372 * t412;
t381 = Ifges(5,4) * t414 + Ifges(5,2) * t413 + Ifges(5,6) * t436;
t382 = Ifges(5,1) * t414 + Ifges(5,4) * t413 + Ifges(5,5) * t436;
t461 = mrSges(5,1) * t357 - mrSges(5,2) * t358 + Ifges(5,5) * t389 + Ifges(5,6) * t388 + Ifges(5,3) * t420 + pkin(4) * t464 + pkin(9) * t469 + t454 * t342 + t450 * t343 + t414 * t381 - t413 * t382;
t424 = mrSges(4,1) * t473 + mrSges(4,2) * t443;
t428 = (t456 * mrSges(4,2) - t452 * mrSges(4,3)) * t479;
t460 = -m(4) * t369 + t442 * mrSges(4,3) + t443 * t424 + t428 * t472 - t462;
t459 = mrSges(6,1) * t351 - mrSges(6,2) * t352 + Ifges(6,5) * t367 + Ifges(6,6) * t366 + Ifges(6,3) * t386 + t372 * t395 - t373 * t394;
t429 = (-t456 * mrSges(3,1) + t452 * mrSges(3,2)) * t479;
t422 = -mrSges(3,2) * t443 + mrSges(3,3) * t472;
t421 = mrSges(3,1) * t443 - mrSges(3,3) * t473;
t390 = t475 - t480;
t380 = Ifges(5,5) * t414 + Ifges(5,6) * t413 + Ifges(5,3) * t436;
t335 = t429 * t472 + m(3) * t391 - mrSges(3,2) * t442 - t421 * t443 + (mrSges(3,3) + mrSges(4,1)) * t432 + t460;
t333 = mrSges(4,2) * t442 + t423 * t443 + t428 * t473 - t465;
t332 = t432 * mrSges(4,2) - t424 * t473 - t467;
t331 = m(3) * t390 - mrSges(3,3) * t431 + (t422 - t423) * t443 + t491 * t442 + (-t428 - t429) * t473 + t465;
t330 = -mrSges(5,1) * t361 + mrSges(5,3) * t358 + Ifges(5,4) * t389 + Ifges(5,2) * t388 + Ifges(5,6) * t420 - pkin(4) * t339 - t380 * t414 + t382 * t436 - t459;
t329 = mrSges(5,2) * t361 - mrSges(5,3) * t357 + Ifges(5,1) * t389 + Ifges(5,4) * t388 + Ifges(5,5) * t420 - pkin(9) * t339 - t342 * t450 + t343 * t454 + t380 * t413 - t381 * t436;
t328 = mrSges(3,1) * t390 - mrSges(3,2) * t391 + mrSges(4,2) * t375 - mrSges(4,3) * t369 + t455 * t329 - t451 * t330 - pkin(8) * t334 - pkin(2) * t333 + qJ(3) * t460 + t496 * t442 + (qJ(3) * mrSges(4,1) + t488) * t432 + t489 * t431 + (t482 * t452 - t481 * t456) * t479;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t471 - mrSges(2,2) * t468 + (mrSges(4,1) * t375 + mrSges(3,2) * t404 - mrSges(3,3) * t390 - mrSges(4,3) * t370 + pkin(3) * t334 - qJ(3) * t332 + t498 * t431 + t490 * t432 + t489 * t442 - t482 * t443 + t483 * t472 + t461) * t485 + (-mrSges(3,1) * t404 - mrSges(4,1) * t369 + mrSges(4,2) * t370 + mrSges(3,3) * t391 - pkin(2) * t332 - pkin(3) * t462 - pkin(8) * t470 - t451 * t329 - t455 * t330 + t490 * t431 + t497 * t432 + t488 * t442 + t481 * t443 - t483 * t473) * t484 + t449 * t328 + pkin(1) * ((t331 * t456 + t335 * t452) * t449 + (-m(3) * t404 - t431 * mrSges(3,2) + t491 * t432 + (t422 * t456 + (-t421 + t424) * t452) * t479 + t467) * t448) + (-t331 * t452 + t456 * t335) * t493; t328; t333; t461; t459;];
tauJ = t1;
