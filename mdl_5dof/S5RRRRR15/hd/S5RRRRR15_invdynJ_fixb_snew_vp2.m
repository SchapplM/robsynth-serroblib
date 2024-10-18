% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR15
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR15_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:25:17
% EndTime: 2024-09-27 22:25:23
% DurationCPUTime: 2.28s
% Computational Cost: add. (125645->303), mult. (353460->413), div. (0->0), fcn. (286297->14), ass. (0->137)
t490 = sin(pkin(6));
t530 = pkin(11) * t490;
t492 = cos(pkin(6));
t529 = pkin(11) * t492;
t491 = sin(pkin(5));
t528 = t491 * pkin(8);
t504 = qJD(1) ^ 2;
t498 = sin(qJ(1));
t503 = cos(qJ(1));
t516 = t498 * g(1) - g(2) * t503;
t471 = qJDD(1) * pkin(1) + t504 * t528 + t516;
t493 = cos(pkin(5));
t527 = t471 * t493;
t526 = t491 ^ 2 * t504;
t494 = sin(qJ(5));
t525 = t490 * t494;
t499 = cos(qJ(5));
t524 = t490 * t499;
t497 = sin(qJ(2));
t523 = t491 * t497;
t502 = cos(qJ(2));
t522 = t491 * t502;
t520 = qJD(1) * qJD(2);
t475 = (qJDD(1) * t497 + t502 * t520) * t491;
t483 = t493 * qJDD(1) + qJDD(2);
t484 = t493 * qJD(1) + qJD(2);
t512 = -g(1) * t503 - g(2) * t498;
t472 = -pkin(1) * t504 + qJDD(1) * t528 + t512;
t513 = -t472 * t497 + t502 * t527;
t428 = pkin(2) * t483 - pkin(9) * t475 + (pkin(2) * t497 * t526 + (pkin(9) * qJD(1) * t484 - g(3)) * t491) * t502 + t513;
t451 = -g(3) * t523 + t502 * t472 + t497 * t527;
t521 = qJD(1) * t491;
t517 = t497 * t521;
t474 = pkin(2) * t484 - pkin(9) * t517;
t476 = (qJDD(1) * t502 - t497 * t520) * t491;
t519 = t502 ^ 2 * t526;
t429 = -pkin(2) * t519 + pkin(9) * t476 - t474 * t484 + t451;
t496 = sin(qJ(3));
t501 = cos(qJ(3));
t408 = t501 * t428 - t429 * t496;
t467 = (-t497 * t496 + t502 * t501) * t521;
t442 = qJD(3) * t467 + t475 * t501 + t476 * t496;
t468 = (t502 * t496 + t497 * t501) * t521;
t481 = qJDD(3) + t483;
t482 = qJD(3) + t484;
t398 = (t467 * t482 - t442) * pkin(10) + (t467 * t468 + t481) * pkin(3) + t408;
t409 = t496 * t428 + t501 * t429;
t441 = -qJD(3) * t468 - t475 * t496 + t476 * t501;
t457 = pkin(3) * t482 - pkin(10) * t468;
t465 = t467 ^ 2;
t400 = -pkin(3) * t465 + pkin(10) * t441 - t457 * t482 + t409;
t495 = sin(qJ(4));
t500 = cos(qJ(4));
t393 = t495 * t398 + t500 * t400;
t452 = t467 * t500 - t468 * t495;
t453 = t467 * t495 + t468 * t500;
t430 = -pkin(4) * t452 - t453 * t530;
t480 = qJD(4) + t482;
t440 = pkin(4) * t480 - t453 * t529;
t416 = -qJD(4) * t453 + t441 * t500 - t442 * t495;
t479 = qJDD(4) + t481;
t510 = t416 * t492 + t479 * t490;
t389 = t510 * pkin(11) + t430 * t452 - t440 * t480 + t393;
t392 = t500 * t398 - t400 * t495;
t417 = qJD(4) * t452 + t441 * t495 + t442 * t500;
t509 = t452 * t492 + t480 * t490;
t435 = t509 * pkin(11);
t388 = pkin(4) * t479 - t417 * t529 - t430 * t453 + t435 * t480 + t392;
t461 = -g(3) * t493 - t471 * t491;
t433 = -pkin(2) * t476 - pkin(9) * t519 + t474 * t517 + t461;
t406 = -pkin(3) * t441 - pkin(10) * t465 + t468 * t457 + t433;
t390 = -pkin(4) * t416 - t417 * t530 - t435 * t452 + t440 * t453 + t406;
t511 = t388 * t492 + t390 * t490;
t385 = -t389 * t494 + t511 * t499;
t419 = -t453 * t494 + t509 * t499;
t395 = qJD(5) * t419 + t417 * t499 + t510 * t494;
t420 = t453 * t499 + t509 * t494;
t404 = -mrSges(6,1) * t419 + mrSges(6,2) * t420;
t410 = -t416 * t490 + t479 * t492 + qJDD(5);
t436 = -t452 * t490 + t480 * t492 + qJD(5);
t411 = -mrSges(6,2) * t436 + mrSges(6,3) * t419;
t381 = m(6) * t385 + mrSges(6,1) * t410 - mrSges(6,3) * t395 - t404 * t420 + t411 * t436;
t386 = t389 * t499 + t511 * t494;
t394 = -qJD(5) * t420 - t417 * t494 + t510 * t499;
t412 = mrSges(6,1) * t436 - mrSges(6,3) * t420;
t382 = m(6) * t386 - mrSges(6,2) * t410 + mrSges(6,3) * t394 + t404 * t419 - t412 * t436;
t387 = -t388 * t490 + t390 * t492;
t384 = m(6) * t387 - mrSges(6,1) * t394 + mrSges(6,2) * t395 - t411 * t419 + t412 * t420;
t364 = -t384 * t490 + (t381 * t499 + t382 * t494) * t492;
t431 = -mrSges(5,1) * t452 + mrSges(5,2) * t453;
t443 = -mrSges(5,2) * t480 + mrSges(5,3) * t452;
t361 = m(5) * t392 + mrSges(5,1) * t479 - mrSges(5,3) * t417 - t431 * t453 + t443 * t480 + t364;
t369 = -t381 * t494 + t499 * t382;
t444 = mrSges(5,1) * t480 - mrSges(5,3) * t453;
t367 = m(5) * t393 - mrSges(5,2) * t479 + mrSges(5,3) * t416 + t431 * t452 - t444 * t480 + t369;
t359 = t500 * t361 + t495 * t367;
t454 = -mrSges(4,1) * t467 + mrSges(4,2) * t468;
t455 = -mrSges(4,2) * t482 + mrSges(4,3) * t467;
t356 = m(4) * t408 + mrSges(4,1) * t481 - mrSges(4,3) * t442 - t454 * t468 + t455 * t482 + t359;
t456 = mrSges(4,1) * t482 - mrSges(4,3) * t468;
t514 = -t361 * t495 + t500 * t367;
t357 = m(4) * t409 - mrSges(4,2) * t481 + mrSges(4,3) * t441 + t454 * t467 - t456 * t482 + t514;
t350 = t501 * t356 + t496 * t357;
t363 = t381 * t524 + t382 * t525 + t492 * t384;
t518 = t502 * t521;
t515 = -t356 * t496 + t501 * t357;
t508 = -m(5) * t406 + t416 * mrSges(5,1) - t417 * mrSges(5,2) + t452 * t443 - t453 * t444 - t363;
t402 = Ifges(6,4) * t420 + Ifges(6,2) * t419 + Ifges(6,6) * t436;
t403 = Ifges(6,1) * t420 + Ifges(6,4) * t419 + Ifges(6,5) * t436;
t371 = mrSges(6,1) * t385 - mrSges(6,2) * t386 + Ifges(6,5) * t395 + Ifges(6,6) * t394 + Ifges(6,3) * t410 + t402 * t420 - t403 * t419;
t401 = Ifges(6,5) * t420 + Ifges(6,6) * t419 + Ifges(6,3) * t436;
t374 = -mrSges(6,1) * t387 + mrSges(6,3) * t386 + Ifges(6,4) * t395 + Ifges(6,2) * t394 + Ifges(6,6) * t410 - t401 * t420 + t403 * t436;
t375 = mrSges(6,2) * t387 - mrSges(6,3) * t385 + Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t410 + t401 * t419 - t402 * t436;
t422 = Ifges(5,4) * t453 + Ifges(5,2) * t452 + Ifges(5,6) * t480;
t423 = Ifges(5,1) * t453 + Ifges(5,4) * t452 + Ifges(5,5) * t480;
t507 = mrSges(5,1) * t392 - mrSges(5,2) * t393 + Ifges(5,5) * t417 + Ifges(5,6) * t416 + Ifges(5,3) * t479 + pkin(4) * t364 + t369 * t530 + t492 * t371 + t374 * t524 + t375 * t525 + t453 * t422 - t452 * t423;
t506 = -m(4) * t433 + t441 * mrSges(4,1) - t442 * mrSges(4,2) + t467 * t455 - t468 * t456 + t508;
t446 = Ifges(4,4) * t468 + Ifges(4,2) * t467 + Ifges(4,6) * t482;
t447 = Ifges(4,1) * t468 + Ifges(4,4) * t467 + Ifges(4,5) * t482;
t505 = mrSges(4,1) * t408 - mrSges(4,2) * t409 + Ifges(4,5) * t442 + Ifges(4,6) * t441 + Ifges(4,3) * t481 + pkin(3) * t359 + t468 * t446 - t467 * t447 + t507;
t473 = (-t502 * mrSges(3,1) + t497 * mrSges(3,2)) * t521;
t470 = -mrSges(3,2) * t484 + mrSges(3,3) * t518;
t469 = mrSges(3,1) * t484 - mrSges(3,3) * t517;
t460 = Ifges(3,5) * t484 + (t497 * Ifges(3,1) + t502 * Ifges(3,4)) * t521;
t459 = Ifges(3,6) * t484 + (t497 * Ifges(3,4) + t502 * Ifges(3,2)) * t521;
t458 = Ifges(3,3) * t484 + (t497 * Ifges(3,5) + t502 * Ifges(3,6)) * t521;
t450 = -g(3) * t522 + t513;
t445 = Ifges(4,5) * t468 + Ifges(4,6) * t467 + Ifges(4,3) * t482;
t421 = Ifges(5,5) * t453 + Ifges(5,6) * t452 + Ifges(5,3) * t480;
t352 = mrSges(5,2) * t406 - mrSges(5,3) * t392 + Ifges(5,1) * t417 + Ifges(5,4) * t416 + Ifges(5,5) * t479 - t374 * t494 + t375 * t499 + t421 * t452 - t422 * t480 + (-t363 * t490 - t364 * t492) * pkin(11);
t351 = -mrSges(5,1) * t406 + mrSges(5,3) * t393 + Ifges(5,4) * t417 + Ifges(5,2) * t416 + Ifges(5,6) * t479 - pkin(4) * t363 - t371 * t490 - t421 * t453 + t423 * t480 + (pkin(11) * t369 + t374 * t499 + t375 * t494) * t492;
t349 = m(3) * t451 - mrSges(3,2) * t483 + mrSges(3,3) * t476 - t469 * t484 + t473 * t518 + t515;
t348 = m(3) * t450 + mrSges(3,1) * t483 - mrSges(3,3) * t475 + t470 * t484 - t473 * t517 + t350;
t347 = mrSges(4,2) * t433 - mrSges(4,3) * t408 + Ifges(4,1) * t442 + Ifges(4,4) * t441 + Ifges(4,5) * t481 - pkin(10) * t359 - t351 * t495 + t352 * t500 + t445 * t467 - t446 * t482;
t346 = (t497 * t459 - t502 * t460) * t521 + Ifges(3,3) * t483 + Ifges(3,5) * t475 + Ifges(3,6) * t476 + mrSges(3,1) * t450 - mrSges(3,2) * t451 + pkin(2) * t350 + t505;
t345 = -mrSges(4,1) * t433 + mrSges(4,3) * t409 + Ifges(4,4) * t442 + Ifges(4,2) * t441 + Ifges(4,6) * t481 + pkin(3) * t508 + pkin(10) * t514 + t500 * t351 + t495 * t352 - t468 * t445 + t482 * t447;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t516 - mrSges(2,2) * t512 + (mrSges(3,2) * t461 - mrSges(3,3) * t450 + Ifges(3,1) * t475 + Ifges(3,4) * t476 + Ifges(3,5) * t483 - pkin(9) * t350 - t345 * t496 + t347 * t501 + t458 * t518 - t484 * t459) * t523 + (-mrSges(3,1) * t461 + mrSges(3,3) * t451 + Ifges(3,4) * t475 + Ifges(3,2) * t476 + Ifges(3,6) * t483 + pkin(2) * t506 + pkin(9) * t515 + t501 * t345 + t496 * t347 - t458 * t517 + t484 * t460) * t522 + t493 * t346 + pkin(1) * ((t502 * t348 + t497 * t349) * t493 + ((-t469 * t497 + t470 * t502) * t521 - t475 * mrSges(3,2) + t476 * mrSges(3,1) - m(3) * t461 + t506) * t491) + (-t348 * t497 + t349 * t502) * t528; t346; t505; t507; t371;];
tauJ = t1;
