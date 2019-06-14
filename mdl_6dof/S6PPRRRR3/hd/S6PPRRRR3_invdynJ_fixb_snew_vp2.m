% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 21:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:13:44
% EndTime: 2019-05-04 21:13:47
% DurationCPUTime: 3.54s
% Computational Cost: add. (43175->239), mult. (82112->325), div. (0->0), fcn. (69591->18), ass. (0->119)
t426 = sin(pkin(13));
t431 = cos(pkin(13));
t420 = -g(1) * t431 - g(2) * t426;
t425 = sin(pkin(14));
t430 = cos(pkin(14));
t419 = g(1) * t426 - g(2) * t431;
t424 = -g(3) + qJDD(1);
t429 = sin(pkin(6));
t434 = cos(pkin(6));
t448 = t419 * t434 + t424 * t429;
t394 = -t425 * t420 + t430 * t448;
t407 = -t419 * t429 + t424 * t434 + qJDD(2);
t428 = sin(pkin(7));
t433 = cos(pkin(7));
t464 = t394 * t433 + t407 * t428;
t427 = sin(pkin(8));
t437 = sin(qJ(4));
t441 = cos(qJ(4));
t455 = qJD(3) * qJD(4);
t416 = (-qJDD(3) * t441 + t437 * t455) * t427;
t395 = t430 * t420 + t425 * t448;
t438 = sin(qJ(3));
t442 = cos(qJ(3));
t369 = -t395 * t438 + t464 * t442;
t443 = qJD(3) ^ 2;
t463 = pkin(10) * t427;
t364 = qJDD(3) * pkin(3) + t443 * t463 + t369;
t370 = t442 * t395 + t464 * t438;
t365 = -pkin(3) * t443 + qJDD(3) * t463 + t370;
t432 = cos(pkin(8));
t378 = -t394 * t428 + t407 * t433;
t461 = t378 * t427;
t356 = -t437 * t365 + t441 * (t364 * t432 + t461);
t423 = qJD(3) * t432 + qJD(4);
t456 = qJD(3) * t427;
t453 = t441 * t456;
t412 = -mrSges(5,2) * t423 + mrSges(5,3) * t453;
t413 = (-mrSges(5,1) * t441 + mrSges(5,2) * t437) * t456;
t415 = (qJDD(3) * t437 + t441 * t455) * t427;
t422 = qJDD(3) * t432 + qJDD(4);
t458 = t432 * t437;
t357 = t364 * t458 + t441 * t365 + t437 * t461;
t414 = (-pkin(4) * t441 - pkin(11) * t437) * t456;
t421 = t423 ^ 2;
t355 = -pkin(4) * t421 + pkin(11) * t422 + t414 * t453 + t357;
t377 = t432 * t378;
t359 = t416 * pkin(4) - t415 * pkin(11) + t377 + (-t364 + (pkin(4) * t437 - pkin(11) * t441) * t423 * qJD(3)) * t427;
t436 = sin(qJ(5));
t440 = cos(qJ(5));
t351 = t440 * t355 + t436 * t359;
t454 = t437 * t456;
t408 = t423 * t440 - t436 * t454;
t409 = t423 * t436 + t440 * t454;
t390 = -pkin(5) * t408 - pkin(12) * t409;
t410 = qJDD(5) + t416;
t418 = qJD(5) - t453;
t417 = t418 ^ 2;
t349 = -pkin(5) * t417 + pkin(12) * t410 + t390 * t408 + t351;
t354 = -t422 * pkin(4) - t421 * pkin(11) + t414 * t454 - t356;
t387 = -qJD(5) * t409 - t415 * t436 + t422 * t440;
t388 = qJD(5) * t408 + t415 * t440 + t422 * t436;
t352 = (-t408 * t418 - t388) * pkin(12) + (t409 * t418 - t387) * pkin(5) + t354;
t435 = sin(qJ(6));
t439 = cos(qJ(6));
t345 = -t349 * t435 + t352 * t439;
t396 = -t409 * t435 + t418 * t439;
t368 = qJD(6) * t396 + t388 * t439 + t410 * t435;
t397 = t409 * t439 + t418 * t435;
t375 = -mrSges(7,1) * t396 + mrSges(7,2) * t397;
t406 = qJD(6) - t408;
t379 = -mrSges(7,2) * t406 + mrSges(7,3) * t396;
t385 = qJDD(6) - t387;
t343 = m(7) * t345 + mrSges(7,1) * t385 - mrSges(7,3) * t368 - t375 * t397 + t379 * t406;
t346 = t349 * t439 + t352 * t435;
t367 = -qJD(6) * t397 - t388 * t435 + t410 * t439;
t380 = mrSges(7,1) * t406 - mrSges(7,3) * t397;
t344 = m(7) * t346 - mrSges(7,2) * t385 + mrSges(7,3) * t367 + t375 * t396 - t380 * t406;
t336 = t343 * t439 + t344 * t435;
t398 = -mrSges(6,2) * t418 + mrSges(6,3) * t408;
t399 = mrSges(6,1) * t418 - mrSges(6,3) * t409;
t446 = -m(6) * t354 + t387 * mrSges(6,1) - mrSges(6,2) * t388 + t408 * t398 - t399 * t409 - t336;
t332 = m(5) * t356 + mrSges(5,1) * t422 - mrSges(5,3) * t415 + t412 * t423 - t413 * t454 + t446;
t462 = t332 * t441;
t411 = mrSges(5,1) * t423 - mrSges(5,3) * t454;
t337 = -t343 * t435 + t439 * t344;
t389 = -mrSges(6,1) * t408 + mrSges(6,2) * t409;
t335 = m(6) * t351 - mrSges(6,2) * t410 + mrSges(6,3) * t387 + t389 * t408 - t399 * t418 + t337;
t350 = -t355 * t436 + t359 * t440;
t348 = -pkin(5) * t410 - pkin(12) * t417 + t390 * t409 - t350;
t347 = -m(7) * t348 + t367 * mrSges(7,1) - mrSges(7,2) * t368 + t396 * t379 - t380 * t397;
t341 = m(6) * t350 + mrSges(6,1) * t410 - mrSges(6,3) * t388 - t389 * t409 + t398 * t418 + t347;
t451 = t440 * t335 - t341 * t436;
t328 = m(5) * t357 - mrSges(5,2) * t422 - mrSges(5,3) * t416 - t411 * t423 + t413 * t453 + t451;
t457 = t328 * t458 + t432 * t462;
t330 = t436 * t335 + t440 * t341;
t452 = t441 * t328 - t332 * t437;
t360 = -t427 * t364 + t377;
t329 = m(5) * t360 + t416 * mrSges(5,1) + t415 * mrSges(5,2) + (t411 * t437 - t412 * t441) * t456 + t330;
t321 = m(4) * t369 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t443 - t329 * t427 + t457;
t323 = m(4) * t370 - mrSges(4,1) * t443 - qJDD(3) * mrSges(4,2) + t452;
t450 = t321 * t442 + t323 * t438;
t372 = Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t406;
t373 = Ifges(7,1) * t397 + Ifges(7,4) * t396 + Ifges(7,5) * t406;
t445 = mrSges(7,1) * t345 - mrSges(7,2) * t346 + Ifges(7,5) * t368 + Ifges(7,6) * t367 + Ifges(7,3) * t385 + t372 * t397 - t373 * t396;
t371 = Ifges(7,5) * t397 + Ifges(7,6) * t396 + Ifges(7,3) * t406;
t338 = -mrSges(7,1) * t348 + mrSges(7,3) * t346 + Ifges(7,4) * t368 + Ifges(7,2) * t367 + Ifges(7,6) * t385 - t371 * t397 + t373 * t406;
t339 = mrSges(7,2) * t348 - mrSges(7,3) * t345 + Ifges(7,1) * t368 + Ifges(7,4) * t367 + Ifges(7,5) * t385 + t371 * t396 - t372 * t406;
t382 = Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * t418;
t383 = Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * t418;
t444 = mrSges(6,1) * t350 - mrSges(6,2) * t351 + Ifges(6,5) * t388 + Ifges(6,6) * t387 + Ifges(6,3) * t410 + pkin(5) * t347 + pkin(12) * t337 + t439 * t338 + t435 * t339 + t409 * t382 - t408 * t383;
t404 = Ifges(5,5) * t423 + (Ifges(5,1) * t437 + Ifges(5,4) * t441) * t456;
t403 = Ifges(5,6) * t423 + (Ifges(5,4) * t437 + Ifges(5,2) * t441) * t456;
t381 = Ifges(6,5) * t409 + Ifges(6,6) * t408 + Ifges(6,3) * t418;
t325 = -mrSges(6,1) * t354 + mrSges(6,3) * t351 + Ifges(6,4) * t388 + Ifges(6,2) * t387 + Ifges(6,6) * t410 - pkin(5) * t336 - t381 * t409 + t383 * t418 - t445;
t324 = mrSges(6,2) * t354 - mrSges(6,3) * t350 + Ifges(6,1) * t388 + Ifges(6,4) * t387 + Ifges(6,5) * t410 - pkin(12) * t336 - t338 * t435 + t339 * t439 + t381 * t408 - t382 * t418;
t322 = m(4) * t378 + t432 * t329 + (t328 * t437 + t462) * t427;
t320 = Ifges(5,5) * t415 - Ifges(5,6) * t416 + Ifges(5,3) * t422 + mrSges(5,1) * t356 - mrSges(5,2) * t357 + t436 * t324 + t440 * t325 + pkin(4) * t446 + pkin(11) * t451 + (t403 * t437 - t404 * t441) * t456;
t319 = m(3) * t407 + t433 * t322 + t428 * t450;
t1 = [m(2) * t424 + t434 * t319 + (t425 * (m(3) * t395 - t321 * t438 + t323 * t442) + t430 * (m(3) * t394 - t428 * t322 + t433 * t450)) * t429; t319; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t369 - mrSges(4,2) * t370 + t432 * t320 + pkin(3) * t457 + (t437 * (mrSges(5,2) * t360 - mrSges(5,3) * t356 + Ifges(5,1) * t415 - Ifges(5,4) * t416 + Ifges(5,5) * t422 - pkin(11) * t330 + t324 * t440 - t325 * t436 - t403 * t423) + t441 * (-mrSges(5,1) * t360 + mrSges(5,3) * t357 + Ifges(5,4) * t415 - Ifges(5,2) * t416 + Ifges(5,6) * t422 - pkin(4) * t330 + t423 * t404 - t444) - pkin(3) * t329 + pkin(10) * t452) * t427; t320; t444; t445;];
tauJ  = t1;
