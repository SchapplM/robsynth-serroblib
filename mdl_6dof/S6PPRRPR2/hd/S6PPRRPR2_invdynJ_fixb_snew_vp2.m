% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:16:05
% EndTime: 2019-05-04 20:16:06
% DurationCPUTime: 1.42s
% Computational Cost: add. (7467->211), mult. (13641->262), div. (0->0), fcn. (9829->14), ass. (0->106)
t428 = sin(qJ(4));
t431 = cos(qJ(4));
t463 = Ifges(5,4) + Ifges(6,6);
t475 = t428 * t463 + t431 * (Ifges(5,2) + Ifges(6,3));
t474 = t428 * (Ifges(5,1) + Ifges(6,2)) + t431 * t463;
t462 = Ifges(5,5) - Ifges(6,4);
t461 = Ifges(5,6) - Ifges(6,5);
t420 = sin(pkin(11));
t424 = cos(pkin(11));
t406 = -g(1) * t424 - g(2) * t420;
t419 = sin(pkin(12));
t423 = cos(pkin(12));
t405 = g(1) * t420 - g(2) * t424;
t418 = -g(3) + qJDD(1);
t422 = sin(pkin(6));
t426 = cos(pkin(6));
t446 = t405 * t426 + t418 * t422;
t366 = -t406 * t419 + t446 * t423;
t381 = -t405 * t422 + t418 * t426 + qJDD(2);
t421 = sin(pkin(7));
t425 = cos(pkin(7));
t473 = t366 * t425 + t381 * t421;
t367 = t406 * t423 + t446 * t419;
t429 = sin(qJ(3));
t432 = cos(qJ(3));
t358 = -t429 * t367 + t473 * t432;
t438 = -qJDD(3) * pkin(3) - t358;
t434 = qJD(3) ^ 2;
t466 = pkin(9) * t434;
t356 = t438 - t466;
t452 = qJD(3) * qJD(4);
t451 = t431 * t452;
t402 = qJDD(3) * t428 + t451;
t450 = t428 * t452;
t403 = qJDD(3) * t431 - t450;
t454 = qJD(3) * t428;
t407 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t454;
t453 = qJD(3) * t431;
t408 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t453;
t409 = -mrSges(6,1) * t453 - qJD(4) * mrSges(6,3);
t468 = -2 * qJD(5);
t435 = pkin(4) * t450 + t454 * t468 + (-t402 - t451) * qJ(5) + t438;
t353 = -pkin(4) * t403 + t435 - t466;
t410 = mrSges(6,1) * t454 + qJD(4) * mrSges(6,2);
t361 = -t366 * t421 + t381 * t425;
t359 = t432 * t367 + t473 * t429;
t357 = -pkin(3) * t434 + qJDD(3) * pkin(9) + t359;
t354 = t428 * t357;
t399 = (-pkin(4) * t431 - qJ(5) * t428) * qJD(3);
t433 = qJD(4) ^ 2;
t445 = -qJ(5) * t433 + t399 * t454 + qJDD(5) + t354;
t465 = pkin(10) * t434;
t467 = -pkin(4) - pkin(10);
t347 = pkin(5) * t402 + t467 * qJDD(4) + (-pkin(5) * t452 - t428 * t465 - t361) * t431 + t445;
t411 = pkin(5) * t454 - qJD(4) * pkin(10);
t417 = t431 ^ 2;
t350 = -t411 * t454 + (-pkin(5) * t417 - pkin(9)) * t434 + t467 * t403 + t435;
t427 = sin(qJ(6));
t430 = cos(qJ(6));
t343 = t347 * t430 - t350 * t427;
t397 = -qJD(4) * t427 - t430 * t453;
t376 = qJD(6) * t397 + qJDD(4) * t430 - t403 * t427;
t398 = qJD(4) * t430 - t427 * t453;
t377 = -mrSges(7,1) * t397 + mrSges(7,2) * t398;
t413 = qJD(6) + t454;
t379 = -mrSges(7,2) * t413 + mrSges(7,3) * t397;
t396 = qJDD(6) + t402;
t340 = m(7) * t343 + mrSges(7,1) * t396 - mrSges(7,3) * t376 - t377 * t398 + t379 * t413;
t344 = t347 * t427 + t350 * t430;
t375 = -qJD(6) * t398 - qJDD(4) * t427 - t403 * t430;
t380 = mrSges(7,1) * t413 - mrSges(7,3) * t398;
t341 = m(7) * t344 - mrSges(7,2) * t396 + mrSges(7,3) * t375 + t377 * t397 - t380 * t413;
t457 = -t427 * t340 + t430 * t341;
t444 = m(6) * t353 + t403 * mrSges(6,2) - t410 * t454 + t457;
t470 = -m(5) * t356 + t403 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t402 + t408 * t453 + (-t407 * t428 - t409 * t431) * qJD(3) - t444;
t469 = (t475 * qJD(3) + t461 * qJD(4)) * t428 - t431 * (t474 * qJD(3) + t462 * qJD(4));
t460 = t361 * t431;
t352 = t431 * t357 + t428 * t361;
t351 = -t354 + t460;
t400 = (mrSges(6,2) * t431 - mrSges(6,3) * t428) * qJD(3);
t401 = (-mrSges(5,1) * t431 + mrSges(5,2) * t428) * qJD(3);
t333 = t340 * t430 + t341 * t427;
t349 = -qJDD(4) * pkin(4) + t445 - t460;
t441 = -m(6) * t349 - t402 * mrSges(6,1) - t333;
t330 = m(5) * t351 - mrSges(5,3) * t402 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t408 - t409) * qJD(4) + (-t400 - t401) * t454 + t441;
t437 = -pkin(4) * t433 + qJDD(4) * qJ(5) + t399 * t453 + t352;
t348 = qJD(4) * t468 - t437;
t346 = -t417 * t465 + pkin(5) * t403 + ((2 * qJD(5)) + t411) * qJD(4) + t437;
t442 = -m(7) * t346 + mrSges(7,1) * t375 - t376 * mrSges(7,2) + t379 * t397 - t398 * t380;
t436 = -m(6) * t348 + qJDD(4) * mrSges(6,3) + qJD(4) * t410 + t400 * t453 - t442;
t335 = t401 * t453 + m(5) * t352 - qJDD(4) * mrSges(5,2) - qJD(4) * t407 + (mrSges(5,3) + mrSges(6,1)) * t403 + t436;
t449 = -t330 * t428 + t431 * t335;
t327 = m(4) * t359 - mrSges(4,1) * t434 - qJDD(3) * mrSges(4,2) + t449;
t329 = m(4) * t358 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t434 + t470;
t448 = t327 * t429 + t329 * t432;
t369 = Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t413;
t370 = Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t413;
t440 = mrSges(7,1) * t343 - mrSges(7,2) * t344 + Ifges(7,5) * t376 + Ifges(7,6) * t375 + Ifges(7,3) * t396 + t398 * t369 - t397 * t370;
t368 = Ifges(7,5) * t398 + Ifges(7,6) * t397 + Ifges(7,3) * t413;
t337 = mrSges(7,2) * t346 - mrSges(7,3) * t343 + Ifges(7,1) * t376 + Ifges(7,4) * t375 + Ifges(7,5) * t396 + t368 * t397 - t369 * t413;
t336 = -mrSges(7,1) * t346 + mrSges(7,3) * t344 + Ifges(7,4) * t376 + Ifges(7,2) * t375 + Ifges(7,6) * t396 - t368 * t398 + t370 * t413;
t332 = -mrSges(6,3) * t402 + t409 * t453 + t444;
t331 = qJDD(4) * mrSges(6,2) + qJD(4) * t409 + t400 * t454 - t441;
t328 = m(4) * t361 + t330 * t431 + t335 * t428;
t326 = m(3) * t381 + t328 * t425 + t448 * t421;
t1 = [m(2) * t418 + t426 * t326 + (t419 * (m(3) * t367 + t327 * t432 - t329 * t429) + t423 * (m(3) * t366 - t328 * t421 + t448 * t425)) * t422; t326; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t358 - mrSges(4,2) * t359 + t428 * (mrSges(6,1) * t349 + mrSges(5,2) * t356 - mrSges(5,3) * t351 - mrSges(6,3) * t353 + pkin(5) * t333 - qJ(5) * t332 + t440) + t431 * (-mrSges(5,1) * t356 - mrSges(6,1) * t348 + mrSges(6,2) * t353 + mrSges(5,3) * t352 - pkin(4) * t332 - pkin(5) * t442 - pkin(10) * t457 - t430 * t336 - t427 * t337) + pkin(9) * t449 + t475 * t403 + (t428 * t462 + t431 * t461) * qJDD(4) - t469 * qJD(4) + t474 * t402 + t470 * pkin(3); mrSges(5,1) * t351 - mrSges(5,2) * t352 + mrSges(6,2) * t349 - mrSges(6,3) * t348 + t430 * t337 - t427 * t336 - pkin(10) * t333 - pkin(4) * t331 + qJ(5) * t436 + (mrSges(6,1) * qJ(5) + t461) * t403 + t462 * t402 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + t469 * qJD(3); t331; t440;];
tauJ  = t1;
