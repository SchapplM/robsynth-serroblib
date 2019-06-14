% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRP5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:01:52
% EndTime: 2019-05-05 10:01:56
% DurationCPUTime: 3.66s
% Computational Cost: add. (38149->279), mult. (78298->354), div. (0->0), fcn. (62173->14), ass. (0->126)
t526 = Ifges(6,1) + Ifges(7,1);
t519 = Ifges(6,4) + Ifges(7,4);
t518 = Ifges(6,5) + Ifges(7,5);
t525 = Ifges(6,2) + Ifges(7,2);
t517 = Ifges(6,6) + Ifges(7,6);
t524 = Ifges(6,3) + Ifges(7,3);
t477 = sin(pkin(12));
t480 = cos(pkin(12));
t469 = g(1) * t477 - g(2) * t480;
t476 = -g(3) + qJDD(1);
t479 = sin(pkin(6));
t482 = cos(pkin(6));
t523 = t469 * t482 + t476 * t479;
t478 = sin(pkin(7));
t485 = sin(qJ(3));
t489 = cos(qJ(3));
t504 = qJD(2) * qJD(3);
t463 = (-qJDD(2) * t489 + t485 * t504) * t478;
t470 = -g(1) * t480 - g(2) * t477;
t486 = sin(qJ(2));
t490 = cos(qJ(2));
t440 = -t470 * t486 + t523 * t490;
t491 = qJD(2) ^ 2;
t521 = pkin(9) * t478;
t436 = qJDD(2) * pkin(2) + t491 * t521 + t440;
t441 = t490 * t470 + t523 * t486;
t437 = -pkin(2) * t491 + qJDD(2) * t521 + t441;
t481 = cos(pkin(7));
t456 = -t469 * t479 + t476 * t482;
t514 = t456 * t478;
t398 = -t485 * t437 + t489 * (t436 * t481 + t514);
t475 = qJD(2) * t481 + qJD(3);
t484 = sin(qJ(4));
t488 = cos(qJ(4));
t505 = qJD(2) * t478;
t501 = t485 * t505;
t454 = t475 * t488 - t484 * t501;
t462 = (qJDD(2) * t485 + t489 * t504) * t478;
t474 = qJDD(2) * t481 + qJDD(3);
t432 = qJD(4) * t454 + t462 * t488 + t474 * t484;
t455 = t475 * t484 + t488 * t501;
t500 = t489 * t505;
t468 = qJD(4) - t500;
t483 = sin(qJ(5));
t487 = cos(qJ(5));
t443 = -t455 * t483 + t468 * t487;
t457 = qJDD(4) + t463;
t404 = qJD(5) * t443 + t432 * t487 + t457 * t483;
t444 = t455 * t487 + t468 * t483;
t415 = -mrSges(7,1) * t443 + mrSges(7,2) * t444;
t511 = t481 * t485;
t399 = t436 * t511 + t489 * t437 + t485 * t514;
t461 = (-pkin(3) * t489 - pkin(10) * t485) * t505;
t473 = t475 ^ 2;
t395 = -pkin(3) * t473 + pkin(10) * t474 + t461 * t500 + t399;
t452 = t481 * t456;
t397 = t463 * pkin(3) - t462 * pkin(10) + t452 + (-t436 + (pkin(3) * t485 - pkin(10) * t489) * t475 * qJD(2)) * t478;
t391 = t488 * t395 + t484 * t397;
t439 = -pkin(4) * t454 - pkin(11) * t455;
t467 = t468 ^ 2;
t386 = -pkin(4) * t467 + pkin(11) * t457 + t439 * t454 + t391;
t394 = -t474 * pkin(3) - t473 * pkin(10) + t461 * t501 - t398;
t431 = -qJD(4) * t455 - t462 * t484 + t474 * t488;
t389 = (-t454 * t468 - t432) * pkin(11) + (t455 * t468 - t431) * pkin(4) + t394;
t381 = -t483 * t386 + t487 * t389;
t429 = qJDD(5) - t431;
t453 = qJD(5) - t454;
t378 = -0.2e1 * qJD(6) * t444 + (t443 * t453 - t404) * qJ(6) + (t443 * t444 + t429) * pkin(5) + t381;
t419 = -mrSges(7,2) * t453 + mrSges(7,3) * t443;
t503 = m(7) * t378 + t429 * mrSges(7,1) + t453 * t419;
t375 = -t404 * mrSges(7,3) - t444 * t415 + t503;
t382 = t487 * t386 + t483 * t389;
t403 = -qJD(5) * t444 - t432 * t483 + t457 * t487;
t421 = pkin(5) * t453 - qJ(6) * t444;
t442 = t443 ^ 2;
t380 = -pkin(5) * t442 + qJ(6) * t403 + 0.2e1 * qJD(6) * t443 - t421 * t453 + t382;
t507 = t519 * t443 + t526 * t444 + t518 * t453;
t508 = -t525 * t443 - t519 * t444 - t517 * t453;
t522 = mrSges(6,1) * t381 + mrSges(7,1) * t378 - mrSges(6,2) * t382 - mrSges(7,2) * t380 + pkin(5) * t375 + t403 * t517 + t404 * t518 + t524 * t429 - t443 * t507 - t444 * t508;
t520 = -mrSges(6,2) - mrSges(7,2);
t459 = -mrSges(4,2) * t475 + mrSges(4,3) * t500;
t460 = (-mrSges(4,1) * t489 + mrSges(4,2) * t485) * t505;
t416 = -mrSges(6,1) * t443 + mrSges(6,2) * t444;
t420 = -mrSges(6,2) * t453 + mrSges(6,3) * t443;
t369 = m(6) * t381 + t429 * mrSges(6,1) + t453 * t420 + (-t415 - t416) * t444 + (-mrSges(6,3) - mrSges(7,3)) * t404 + t503;
t502 = m(7) * t380 + t403 * mrSges(7,3) + t443 * t415;
t422 = mrSges(7,1) * t453 - mrSges(7,3) * t444;
t506 = -mrSges(6,1) * t453 + mrSges(6,3) * t444 - t422;
t371 = m(6) * t382 + t403 * mrSges(6,3) + t443 * t416 + t429 * t520 + t453 * t506 + t502;
t367 = t369 * t487 + t371 * t483;
t445 = -mrSges(5,2) * t468 + mrSges(5,3) * t454;
t446 = mrSges(5,1) * t468 - mrSges(5,3) * t455;
t493 = -m(5) * t394 + t431 * mrSges(5,1) - mrSges(5,2) * t432 + t454 * t445 - t446 * t455 - t367;
t362 = m(4) * t398 + mrSges(4,1) * t474 - mrSges(4,3) * t462 + t459 * t475 - t460 * t501 + t493;
t515 = t362 * t489;
t458 = mrSges(4,1) * t475 - mrSges(4,3) * t501;
t368 = -t369 * t483 + t487 * t371;
t438 = -mrSges(5,1) * t454 + mrSges(5,2) * t455;
t365 = m(5) * t391 - mrSges(5,2) * t457 + mrSges(5,3) * t431 + t438 * t454 - t446 * t468 + t368;
t390 = -t484 * t395 + t397 * t488;
t385 = -pkin(4) * t457 - pkin(11) * t467 + t455 * t439 - t390;
t383 = -pkin(5) * t403 - qJ(6) * t442 + t421 * t444 + qJDD(6) + t385;
t497 = -m(7) * t383 + t403 * mrSges(7,1) + t443 * t419;
t374 = -m(6) * t385 + t403 * mrSges(6,1) + t404 * t520 + t443 * t420 + t444 * t506 + t497;
t373 = m(5) * t390 + t457 * mrSges(5,1) - t432 * mrSges(5,3) - t455 * t438 + t468 * t445 + t374;
t498 = t488 * t365 - t373 * t484;
t357 = m(4) * t399 - mrSges(4,2) * t474 - mrSges(4,3) * t463 - t458 * t475 + t460 * t500 + t498;
t510 = t357 * t511 + t481 * t515;
t359 = t484 * t365 + t488 * t373;
t509 = -t517 * t443 - t518 * t444 - t524 * t453;
t499 = t489 * t357 - t362 * t485;
t377 = t404 * mrSges(7,2) + t444 * t422 - t497;
t360 = -mrSges(6,1) * t385 + mrSges(6,3) * t382 - mrSges(7,1) * t383 + mrSges(7,3) * t380 - pkin(5) * t377 + qJ(6) * t502 + (-qJ(6) * t422 + t507) * t453 + t509 * t444 + (-mrSges(7,2) * qJ(6) + t517) * t429 + t519 * t404 + t525 * t403;
t366 = mrSges(6,2) * t385 + mrSges(7,2) * t383 - mrSges(6,3) * t381 - mrSges(7,3) * t378 - qJ(6) * t375 + t519 * t403 + t526 * t404 + t518 * t429 - t509 * t443 + t508 * t453;
t426 = Ifges(5,4) * t455 + Ifges(5,2) * t454 + Ifges(5,6) * t468;
t427 = Ifges(5,1) * t455 + Ifges(5,4) * t454 + Ifges(5,5) * t468;
t492 = mrSges(5,1) * t390 - mrSges(5,2) * t391 + Ifges(5,5) * t432 + Ifges(5,6) * t431 + Ifges(5,3) * t457 + pkin(4) * t374 + pkin(11) * t368 + t487 * t360 + t483 * t366 + t455 * t426 - t454 * t427;
t450 = Ifges(4,5) * t475 + (Ifges(4,1) * t485 + Ifges(4,4) * t489) * t505;
t449 = Ifges(4,6) * t475 + (Ifges(4,4) * t485 + Ifges(4,2) * t489) * t505;
t425 = Ifges(5,5) * t455 + Ifges(5,6) * t454 + Ifges(5,3) * t468;
t417 = -t478 * t436 + t452;
t358 = m(4) * t417 + t463 * mrSges(4,1) + t462 * mrSges(4,2) + (t458 * t485 - t459 * t489) * t505 + t359;
t354 = -mrSges(5,1) * t394 + mrSges(5,3) * t391 + Ifges(5,4) * t432 + Ifges(5,2) * t431 + Ifges(5,6) * t457 - pkin(4) * t367 - t455 * t425 + t468 * t427 - t522;
t353 = mrSges(5,2) * t394 - mrSges(5,3) * t390 + Ifges(5,1) * t432 + Ifges(5,4) * t431 + Ifges(5,5) * t457 - pkin(11) * t367 - t360 * t483 + t366 * t487 + t425 * t454 - t426 * t468;
t352 = Ifges(4,5) * t462 - Ifges(4,6) * t463 + Ifges(4,3) * t474 + mrSges(4,1) * t398 - mrSges(4,2) * t399 + t484 * t353 + t488 * t354 + pkin(3) * t493 + pkin(10) * t498 + (t449 * t485 - t450 * t489) * t505;
t1 = [m(2) * t476 + t482 * (m(3) * t456 + t481 * t358 + (t357 * t485 + t515) * t478) + (t486 * (m(3) * t441 - mrSges(3,1) * t491 - qJDD(2) * mrSges(3,2) + t499) + t490 * (m(3) * t440 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t491 - t358 * t478 + t510)) * t479; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t440 - mrSges(3,2) * t441 + t481 * t352 + pkin(2) * t510 + (t485 * (mrSges(4,2) * t417 - mrSges(4,3) * t398 + Ifges(4,1) * t462 - Ifges(4,4) * t463 + Ifges(4,5) * t474 - pkin(10) * t359 + t353 * t488 - t354 * t484 - t449 * t475) + t489 * (-mrSges(4,1) * t417 + mrSges(4,3) * t399 + Ifges(4,4) * t462 - Ifges(4,2) * t463 + Ifges(4,6) * t474 - pkin(3) * t359 + t475 * t450 - t492) - pkin(2) * t358 + pkin(9) * t499) * t478; t352; t492; t522; t377;];
tauJ  = t1;
