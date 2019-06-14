% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:08:09
% EndTime: 2019-05-04 22:08:12
% DurationCPUTime: 1.58s
% Computational Cost: add. (13325->236), mult. (26603->308), div. (0->0), fcn. (18803->14), ass. (0->107)
t425 = sin(pkin(12));
t429 = cos(pkin(12));
t434 = sin(qJ(4));
t437 = cos(qJ(4));
t401 = (t425 * t434 - t429 * t437) * qJD(2);
t460 = 2 * qJD(5);
t427 = sin(pkin(10));
t431 = cos(pkin(10));
t415 = g(1) * t427 - g(2) * t431;
t432 = cos(pkin(6));
t459 = t415 * t432;
t428 = sin(pkin(6));
t435 = sin(qJ(2));
t458 = t428 * t435;
t438 = cos(qJ(2));
t457 = t428 * t438;
t416 = -g(1) * t431 - g(2) * t427;
t424 = -g(3) + qJDD(1);
t382 = -t416 * t435 + t424 * t457 + t438 * t459;
t377 = qJDD(2) * pkin(2) + t382;
t383 = t438 * t416 + t424 * t458 + t435 * t459;
t440 = qJD(2) ^ 2;
t378 = -pkin(2) * t440 + t383;
t426 = sin(pkin(11));
t430 = cos(pkin(11));
t363 = t426 * t377 + t430 * t378;
t361 = -pkin(3) * t440 + qJDD(2) * pkin(8) + t363;
t447 = -t415 * t428 + t432 * t424;
t396 = qJDD(3) + t447;
t357 = -t434 * t361 + t437 * t396;
t453 = qJD(2) * qJD(4);
t451 = t437 * t453;
t413 = qJDD(2) * t434 + t451;
t354 = (-t413 + t451) * qJ(5) + (t434 * t437 * t440 + qJDD(4)) * pkin(4) + t357;
t358 = t437 * t361 + t434 * t396;
t414 = qJDD(2) * t437 - t434 * t453;
t455 = qJD(2) * t434;
t417 = qJD(4) * pkin(4) - qJ(5) * t455;
t423 = t437 ^ 2;
t355 = -pkin(4) * t423 * t440 + qJ(5) * t414 - qJD(4) * t417 + t358;
t350 = t425 * t354 + t429 * t355 - t401 * t460;
t402 = (t425 * t437 + t429 * t434) * qJD(2);
t385 = mrSges(6,1) * t401 + mrSges(6,2) * t402;
t389 = -t413 * t425 + t414 * t429;
t398 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t402;
t386 = pkin(5) * t401 - pkin(9) * t402;
t439 = qJD(4) ^ 2;
t348 = -pkin(5) * t439 + qJDD(4) * pkin(9) - t386 * t401 + t350;
t362 = t430 * t377 - t426 * t378;
t444 = -qJDD(2) * pkin(3) - t362;
t356 = -t414 * pkin(4) + qJDD(5) + t417 * t455 + (-qJ(5) * t423 - pkin(8)) * t440 + t444;
t390 = t413 * t429 + t414 * t425;
t351 = (qJD(4) * t401 - t390) * pkin(9) + (qJD(4) * t402 - t389) * pkin(5) + t356;
t433 = sin(qJ(6));
t436 = cos(qJ(6));
t345 = -t348 * t433 + t351 * t436;
t393 = qJD(4) * t436 - t402 * t433;
t370 = qJD(6) * t393 + qJDD(4) * t433 + t390 * t436;
t394 = qJD(4) * t433 + t402 * t436;
t371 = -mrSges(7,1) * t393 + mrSges(7,2) * t394;
t400 = qJD(6) + t401;
t372 = -mrSges(7,2) * t400 + mrSges(7,3) * t393;
t388 = qJDD(6) - t389;
t343 = m(7) * t345 + mrSges(7,1) * t388 - mrSges(7,3) * t370 - t371 * t394 + t372 * t400;
t346 = t348 * t436 + t351 * t433;
t369 = -qJD(6) * t394 + qJDD(4) * t436 - t390 * t433;
t373 = mrSges(7,1) * t400 - mrSges(7,3) * t394;
t344 = m(7) * t346 - mrSges(7,2) * t388 + mrSges(7,3) * t369 + t371 * t393 - t373 * t400;
t448 = -t343 * t433 + t436 * t344;
t333 = m(6) * t350 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t389 - qJD(4) * t398 - t385 * t401 + t448;
t446 = -t429 * t354 + t425 * t355;
t349 = -0.2e1 * qJD(5) * t402 - t446;
t397 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t401;
t347 = -qJDD(4) * pkin(5) - t439 * pkin(9) + (t460 + t386) * t402 + t446;
t443 = -m(7) * t347 + t369 * mrSges(7,1) - mrSges(7,2) * t370 + t393 * t372 - t373 * t394;
t339 = m(6) * t349 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t390 + qJD(4) * t397 - t385 * t402 + t443;
t328 = t425 * t333 + t429 * t339;
t412 = (-mrSges(5,1) * t437 + mrSges(5,2) * t434) * qJD(2);
t454 = qJD(2) * t437;
t419 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t454;
t326 = m(5) * t357 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t413 + qJD(4) * t419 - t412 * t455 + t328;
t418 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t455;
t449 = t429 * t333 - t339 * t425;
t327 = m(5) * t358 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t414 - qJD(4) * t418 + t412 * t454 + t449;
t450 = -t326 * t434 + t437 * t327;
t320 = m(4) * t363 - mrSges(4,1) * t440 - qJDD(2) * mrSges(4,2) + t450;
t335 = t436 * t343 + t433 * t344;
t334 = m(6) * t356 - t389 * mrSges(6,1) + mrSges(6,2) * t390 + t401 * t397 + t398 * t402 + t335;
t360 = -t440 * pkin(8) + t444;
t441 = -m(5) * t360 + t414 * mrSges(5,1) - mrSges(5,2) * t413 - t418 * t455 + t419 * t454 - t334;
t330 = m(4) * t362 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t440 + t441;
t456 = t426 * t320 + t430 * t330;
t452 = m(4) * t396 + t437 * t326 + t434 * t327;
t365 = Ifges(7,4) * t394 + Ifges(7,2) * t393 + Ifges(7,6) * t400;
t366 = Ifges(7,1) * t394 + Ifges(7,4) * t393 + Ifges(7,5) * t400;
t442 = mrSges(7,1) * t345 - mrSges(7,2) * t346 + Ifges(7,5) * t370 + Ifges(7,6) * t369 + Ifges(7,3) * t388 + t365 * t394 - t366 * t393;
t407 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t434 + Ifges(5,4) * t437) * qJD(2);
t406 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t434 + Ifges(5,2) * t437) * qJD(2);
t381 = Ifges(6,1) * t402 - Ifges(6,4) * t401 + Ifges(6,5) * qJD(4);
t380 = Ifges(6,4) * t402 - Ifges(6,2) * t401 + Ifges(6,6) * qJD(4);
t379 = Ifges(6,5) * t402 - Ifges(6,6) * t401 + Ifges(6,3) * qJD(4);
t364 = Ifges(7,5) * t394 + Ifges(7,6) * t393 + Ifges(7,3) * t400;
t337 = mrSges(7,2) * t347 - mrSges(7,3) * t345 + Ifges(7,1) * t370 + Ifges(7,4) * t369 + Ifges(7,5) * t388 + t364 * t393 - t365 * t400;
t336 = -mrSges(7,1) * t347 + mrSges(7,3) * t346 + Ifges(7,4) * t370 + Ifges(7,2) * t369 + Ifges(7,6) * t388 - t364 * t394 + t366 * t400;
t322 = -mrSges(6,1) * t356 + mrSges(6,3) * t350 + Ifges(6,4) * t390 + Ifges(6,2) * t389 + Ifges(6,6) * qJDD(4) - pkin(5) * t335 + qJD(4) * t381 - t379 * t402 - t442;
t321 = mrSges(6,2) * t356 - mrSges(6,3) * t349 + Ifges(6,1) * t390 + Ifges(6,4) * t389 + Ifges(6,5) * qJDD(4) - pkin(9) * t335 - qJD(4) * t380 - t336 * t433 + t337 * t436 - t379 * t401;
t1 = [m(2) * t424 + (m(3) * t383 - mrSges(3,1) * t440 - qJDD(2) * mrSges(3,2) + t320 * t430 - t330 * t426) * t458 + (m(3) * t382 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t440 + t456) * t457 + t432 * (m(3) * t447 + t452); mrSges(3,1) * t382 - mrSges(3,2) * t383 + mrSges(4,1) * t362 - mrSges(4,2) * t363 + t434 * (mrSges(5,2) * t360 - mrSges(5,3) * t357 + Ifges(5,1) * t413 + Ifges(5,4) * t414 + Ifges(5,5) * qJDD(4) - qJ(5) * t328 - qJD(4) * t406 + t321 * t429 - t322 * t425) + t437 * (-mrSges(5,1) * t360 + mrSges(5,3) * t358 + Ifges(5,4) * t413 + Ifges(5,2) * t414 + Ifges(5,6) * qJDD(4) - pkin(4) * t334 + qJ(5) * t449 + qJD(4) * t407 + t425 * t321 + t429 * t322) + pkin(3) * t441 + pkin(8) * t450 + pkin(2) * t456 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t452; Ifges(5,5) * t413 + Ifges(5,6) * t414 + mrSges(5,1) * t357 - mrSges(5,2) * t358 + Ifges(6,5) * t390 + Ifges(6,6) * t389 + t402 * t380 + t401 * t381 + mrSges(6,1) * t349 - mrSges(6,2) * t350 + t433 * t337 + t436 * t336 + pkin(5) * t443 + pkin(9) * t448 + pkin(4) * t328 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t406 * t434 - t407 * t437) * qJD(2); t334; t442;];
tauJ  = t1;
