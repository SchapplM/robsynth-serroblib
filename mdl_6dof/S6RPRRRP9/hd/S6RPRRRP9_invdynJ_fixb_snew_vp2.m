% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP9
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
% Datum: 2019-05-06 01:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:51:09
% EndTime: 2019-05-06 01:51:15
% DurationCPUTime: 2.72s
% Computational Cost: add. (18897->267), mult. (36902->320), div. (0->0), fcn. (23829->8), ass. (0->107)
t472 = sin(qJ(4));
t476 = cos(qJ(4));
t477 = cos(qJ(3));
t498 = qJD(1) * t477;
t456 = qJD(3) * t472 + t476 * t498;
t473 = sin(qJ(3));
t497 = qJD(1) * qJD(3);
t493 = t473 * t497;
t460 = qJDD(1) * t477 - t493;
t429 = -qJD(4) * t456 + qJDD(3) * t476 - t460 * t472;
t455 = qJD(3) * t476 - t472 * t498;
t430 = qJD(4) * t455 + qJDD(3) * t472 + t460 * t476;
t471 = sin(qJ(5));
t475 = cos(qJ(5));
t432 = t455 * t475 - t456 * t471;
t398 = qJD(5) * t432 + t429 * t471 + t430 * t475;
t433 = t455 * t471 + t456 * t475;
t409 = -mrSges(7,1) * t432 + mrSges(7,2) * t433;
t480 = qJD(1) ^ 2;
t502 = -pkin(1) - pkin(7);
t474 = sin(qJ(1));
t478 = cos(qJ(1));
t488 = -t478 * g(1) - t474 * g(2);
t503 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t488;
t441 = t480 * t502 - t503;
t465 = t477 * t497;
t459 = -t473 * qJDD(1) - t465;
t413 = (-t460 + t493) * pkin(8) + (-t459 + t465) * pkin(3) + t441;
t492 = t474 * g(1) - t478 * g(2);
t485 = -t480 * qJ(2) + qJDD(2) - t492;
t442 = qJDD(1) * t502 + t485;
t436 = -g(3) * t477 + t473 * t442;
t458 = (pkin(3) * t473 - pkin(8) * t477) * qJD(1);
t467 = t473 * qJD(1);
t479 = qJD(3) ^ 2;
t417 = -pkin(3) * t479 + qJDD(3) * pkin(8) - t458 * t467 + t436;
t385 = t476 * t413 - t472 * t417;
t454 = qJDD(4) - t459;
t464 = t467 + qJD(4);
t381 = (t455 * t464 - t430) * pkin(9) + (t455 * t456 + t454) * pkin(4) + t385;
t386 = t472 * t413 + t476 * t417;
t440 = pkin(4) * t464 - pkin(9) * t456;
t453 = t455 ^ 2;
t383 = -pkin(4) * t453 + pkin(9) * t429 - t440 * t464 + t386;
t375 = t475 * t381 - t471 * t383;
t452 = qJDD(5) + t454;
t463 = qJD(5) + t464;
t371 = -0.2e1 * qJD(6) * t433 + (t432 * t463 - t398) * qJ(6) + (t432 * t433 + t452) * pkin(5) + t375;
t418 = -mrSges(7,2) * t463 + mrSges(7,3) * t432;
t495 = m(7) * t371 + t452 * mrSges(7,1) + t463 * t418;
t367 = -t398 * mrSges(7,3) - t433 * t409 + t495;
t376 = t471 * t381 + t475 * t383;
t397 = -qJD(5) * t433 + t429 * t475 - t430 * t471;
t420 = pkin(5) * t463 - qJ(6) * t433;
t431 = t432 ^ 2;
t373 = -pkin(5) * t431 + qJ(6) * t397 + 0.2e1 * qJD(6) * t432 - t420 * t463 + t376;
t501 = Ifges(6,4) + Ifges(7,4);
t508 = Ifges(6,5) + Ifges(7,5);
t509 = Ifges(6,1) + Ifges(7,1);
t504 = t501 * t432 + t509 * t433 + t508 * t463;
t507 = Ifges(6,6) + Ifges(7,6);
t511 = Ifges(6,2) + Ifges(7,2);
t505 = t511 * t432 + t501 * t433 + t507 * t463;
t506 = Ifges(6,3) + Ifges(7,3);
t512 = mrSges(6,1) * t375 + mrSges(7,1) * t371 - mrSges(6,2) * t376 - mrSges(7,2) * t373 + pkin(5) * t367 + t507 * t397 + t508 * t398 - t432 * t504 + t505 * t433 + t506 * t452;
t410 = -mrSges(6,1) * t432 + mrSges(6,2) * t433;
t419 = -mrSges(6,2) * t463 + mrSges(6,3) * t432;
t361 = m(6) * t375 + t452 * mrSges(6,1) + t463 * t419 + (-t409 - t410) * t433 + (-mrSges(6,3) - mrSges(7,3)) * t398 + t495;
t421 = mrSges(7,1) * t463 - mrSges(7,3) * t433;
t422 = mrSges(6,1) * t463 - mrSges(6,3) * t433;
t494 = m(7) * t373 + t397 * mrSges(7,3) + t432 * t409;
t364 = m(6) * t376 + t397 * mrSges(6,3) + t432 * t410 + (-t421 - t422) * t463 + (-mrSges(6,2) - mrSges(7,2)) * t452 + t494;
t359 = t475 * t361 + t471 * t364;
t424 = Ifges(5,4) * t456 + Ifges(5,2) * t455 + Ifges(5,6) * t464;
t425 = Ifges(5,1) * t456 + Ifges(5,4) * t455 + Ifges(5,5) * t464;
t510 = mrSges(5,1) * t385 - mrSges(5,2) * t386 + Ifges(5,5) * t430 + Ifges(5,6) * t429 + Ifges(5,3) * t454 + pkin(4) * t359 + t456 * t424 - t455 * t425 + t512;
t434 = -mrSges(5,1) * t455 + mrSges(5,2) * t456;
t437 = -mrSges(5,2) * t464 + mrSges(5,3) * t455;
t356 = m(5) * t385 + mrSges(5,1) * t454 - mrSges(5,3) * t430 - t434 * t456 + t437 * t464 + t359;
t438 = mrSges(5,1) * t464 - mrSges(5,3) * t456;
t489 = -t361 * t471 + t475 * t364;
t357 = m(5) * t386 - mrSges(5,2) * t454 + mrSges(5,3) * t429 + t434 * t455 - t438 * t464 + t489;
t351 = t476 * t356 + t472 * t357;
t500 = -t507 * t432 - t508 * t433 - t506 * t463;
t490 = -t356 * t472 + t476 * t357;
t435 = g(3) * t473 + t442 * t477;
t416 = -qJDD(3) * pkin(3) - pkin(8) * t479 + t458 * t498 - t435;
t384 = -pkin(4) * t429 - pkin(9) * t453 + t456 * t440 + t416;
t378 = -pkin(5) * t397 - qJ(6) * t431 + t420 * t433 + qJDD(6) + t384;
t368 = m(7) * t378 - t397 * mrSges(7,1) + t398 * mrSges(7,2) - t432 * t418 + t433 * t421;
t457 = (mrSges(4,1) * t473 + mrSges(4,2) * t477) * qJD(1);
t461 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t467;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t498;
t484 = m(6) * t384 - t397 * mrSges(6,1) + t398 * mrSges(6,2) - t432 * t419 + t433 * t422 + t368;
t482 = -m(5) * t416 + t429 * mrSges(5,1) - t430 * mrSges(5,2) + t455 * t437 - t456 * t438 - t484;
t487 = (m(4) * t436 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t459 - qJD(3) * t462 - t457 * t467 + t490) * t473 + (m(4) * t435 + qJDD(3) * mrSges(4,1) - t460 * mrSges(4,3) + qJD(3) * t461 - t457 * t498 + t482) * t477;
t451 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t477 - Ifges(4,4) * t473) * qJD(1);
t450 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t477 - Ifges(4,2) * t473) * qJD(1);
t444 = -qJDD(1) * pkin(1) + t485;
t443 = t480 * pkin(1) + t503;
t423 = Ifges(5,5) * t456 + Ifges(5,6) * t455 + Ifges(5,3) * t464;
t358 = mrSges(6,2) * t384 + mrSges(7,2) * t378 - mrSges(6,3) * t375 - mrSges(7,3) * t371 - qJ(6) * t367 + t501 * t397 + t509 * t398 - t500 * t432 + t508 * t452 - t505 * t463;
t352 = -mrSges(6,1) * t384 + mrSges(6,3) * t376 - mrSges(7,1) * t378 + mrSges(7,3) * t373 - pkin(5) * t368 + qJ(6) * t494 + (-qJ(6) * t421 + t504) * t463 + (-mrSges(7,2) * qJ(6) + t507) * t452 + t500 * t433 + t501 * t398 + t511 * t397;
t349 = m(3) * t444 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t480 + t487;
t348 = mrSges(5,2) * t416 - mrSges(5,3) * t385 + Ifges(5,1) * t430 + Ifges(5,4) * t429 + Ifges(5,5) * t454 - pkin(9) * t359 - t352 * t471 + t358 * t475 + t423 * t455 - t424 * t464;
t347 = -mrSges(5,1) * t416 + mrSges(5,3) * t386 + Ifges(5,4) * t430 + Ifges(5,2) * t429 + Ifges(5,6) * t454 - pkin(4) * t484 + pkin(9) * t489 + t475 * t352 + t471 * t358 - t456 * t423 + t464 * t425;
t1 = [mrSges(2,1) * t492 - mrSges(2,2) * t488 + mrSges(3,2) * t444 - mrSges(3,3) * t443 + t477 * (mrSges(4,2) * t441 - mrSges(4,3) * t435 + Ifges(4,1) * t460 + Ifges(4,4) * t459 + Ifges(4,5) * qJDD(3) - pkin(8) * t351 - qJD(3) * t450 - t347 * t472 + t348 * t476) - t473 * (-mrSges(4,1) * t441 + mrSges(4,3) * t436 + Ifges(4,4) * t460 + Ifges(4,2) * t459 + Ifges(4,6) * qJDD(3) - pkin(3) * t351 + qJD(3) * t451 - t510) - pkin(7) * t487 - pkin(1) * t349 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t443 + m(4) * t441 - t459 * mrSges(4,1) + t480 * mrSges(3,2) + t460 * mrSges(4,2) + t351 + qJDD(1) * mrSges(3,3) + (t461 * t473 + t462 * t477) * qJD(1)) * qJ(2); t349; Ifges(4,5) * t460 + Ifges(4,6) * t459 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t435 - mrSges(4,2) * t436 + t472 * t348 + t476 * t347 + pkin(3) * t482 + pkin(8) * t490 + (t450 * t477 + t451 * t473) * qJD(1); t510; t512; t368;];
tauJ  = t1;
