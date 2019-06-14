% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR2
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
% Datum: 2019-05-04 22:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:20:44
% EndTime: 2019-05-04 22:20:46
% DurationCPUTime: 1.55s
% Computational Cost: add. (15161->234), mult. (29520->303), div. (0->0), fcn. (20853->14), ass. (0->105)
t416 = sin(pkin(10));
t420 = cos(pkin(10));
t404 = t416 * g(1) - t420 * g(2);
t421 = cos(pkin(6));
t445 = t404 * t421;
t417 = sin(pkin(6));
t424 = sin(qJ(2));
t444 = t417 * t424;
t427 = cos(qJ(2));
t443 = t417 * t427;
t405 = -t420 * g(1) - t416 * g(2);
t413 = -g(3) + qJDD(1);
t367 = -t424 * t405 + t413 * t443 + t427 * t445;
t365 = qJDD(2) * pkin(2) + t367;
t368 = t427 * t405 + t413 * t444 + t424 * t445;
t429 = qJD(2) ^ 2;
t366 = -t429 * pkin(2) + t368;
t415 = sin(pkin(11));
t419 = cos(pkin(11));
t351 = t415 * t365 + t419 * t366;
t349 = -t429 * pkin(3) + qJDD(2) * pkin(8) + t351;
t433 = -t417 * t404 + t421 * t413;
t379 = qJDD(3) + t433;
t423 = sin(qJ(4));
t426 = cos(qJ(4));
t345 = t426 * t349 + t423 * t379;
t401 = (-t426 * mrSges(5,1) + t423 * mrSges(5,2)) * qJD(2);
t439 = qJD(2) * qJD(4);
t412 = t423 * t439;
t403 = t426 * qJDD(2) - t412;
t441 = t423 * qJD(2);
t406 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t441;
t400 = (-t426 * pkin(4) - t423 * qJ(5)) * qJD(2);
t428 = qJD(4) ^ 2;
t440 = t426 * qJD(2);
t340 = -t428 * pkin(4) + qJDD(4) * qJ(5) + t400 * t440 + t345;
t350 = t419 * t365 - t415 * t366;
t348 = -qJDD(2) * pkin(3) - t429 * pkin(8) - t350;
t437 = t426 * t439;
t402 = t423 * qJDD(2) + t437;
t343 = (-t402 - t437) * qJ(5) + (-t403 + t412) * pkin(4) + t348;
t414 = sin(pkin(12));
t418 = cos(pkin(12));
t396 = t414 * qJD(4) + t418 * t441;
t335 = -0.2e1 * qJD(5) * t396 - t414 * t340 + t418 * t343;
t383 = t414 * qJDD(4) + t418 * t402;
t395 = t418 * qJD(4) - t414 * t441;
t333 = (-t395 * t440 - t383) * pkin(9) + (t395 * t396 - t403) * pkin(5) + t335;
t336 = 0.2e1 * qJD(5) * t395 + t418 * t340 + t414 * t343;
t382 = t418 * qJDD(4) - t414 * t402;
t384 = -pkin(5) * t440 - t396 * pkin(9);
t394 = t395 ^ 2;
t334 = -t394 * pkin(5) + t382 * pkin(9) + t384 * t440 + t336;
t422 = sin(qJ(6));
t425 = cos(qJ(6));
t331 = t425 * t333 - t422 * t334;
t373 = t425 * t395 - t422 * t396;
t354 = t373 * qJD(6) + t422 * t382 + t425 * t383;
t374 = t422 * t395 + t425 * t396;
t359 = -t373 * mrSges(7,1) + t374 * mrSges(7,2);
t410 = qJD(6) - t440;
t363 = -t410 * mrSges(7,2) + t373 * mrSges(7,3);
t398 = qJDD(6) - t403;
t326 = m(7) * t331 + t398 * mrSges(7,1) - t354 * mrSges(7,3) - t374 * t359 + t410 * t363;
t332 = t422 * t333 + t425 * t334;
t353 = -t374 * qJD(6) + t425 * t382 - t422 * t383;
t364 = t410 * mrSges(7,1) - t374 * mrSges(7,3);
t327 = m(7) * t332 - t398 * mrSges(7,2) + t353 * mrSges(7,3) + t373 * t359 - t410 * t364;
t320 = t425 * t326 + t422 * t327;
t375 = -t395 * mrSges(6,1) + t396 * mrSges(6,2);
t380 = mrSges(6,2) * t440 + t395 * mrSges(6,3);
t318 = m(6) * t335 - t403 * mrSges(6,1) - t383 * mrSges(6,3) - t396 * t375 - t380 * t440 + t320;
t381 = -mrSges(6,1) * t440 - t396 * mrSges(6,3);
t434 = -t422 * t326 + t425 * t327;
t319 = m(6) * t336 + t403 * mrSges(6,2) + t382 * mrSges(6,3) + t395 * t375 + t381 * t440 + t434;
t435 = -t414 * t318 + t418 * t319;
t315 = m(5) * t345 - qJDD(4) * mrSges(5,2) + t403 * mrSges(5,3) - qJD(4) * t406 + t401 * t440 + t435;
t344 = -t423 * t349 + t426 * t379;
t339 = -qJDD(4) * pkin(4) - t428 * qJ(5) + t400 * t441 + qJDD(5) - t344;
t337 = -t382 * pkin(5) - t394 * pkin(9) + t396 * t384 + t339;
t432 = m(7) * t337 - t353 * mrSges(7,1) + t354 * mrSges(7,2) - t373 * t363 + t374 * t364;
t330 = m(6) * t339 - t382 * mrSges(6,1) + t383 * mrSges(6,2) - t395 * t380 + t396 * t381 + t432;
t407 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t440;
t329 = m(5) * t344 + qJDD(4) * mrSges(5,1) - t402 * mrSges(5,3) + qJD(4) * t407 - t401 * t441 - t330;
t436 = t426 * t315 - t423 * t329;
t308 = m(4) * t351 - t429 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t436;
t316 = t418 * t318 + t414 * t319;
t431 = -m(5) * t348 + t403 * mrSges(5,1) - t402 * mrSges(5,2) - t406 * t441 + t407 * t440 - t316;
t312 = m(4) * t350 + qJDD(2) * mrSges(4,1) - t429 * mrSges(4,2) + t431;
t442 = t415 * t308 + t419 * t312;
t438 = m(4) * t379 + t423 * t315 + t426 * t329;
t356 = Ifges(7,4) * t374 + Ifges(7,2) * t373 + Ifges(7,6) * t410;
t357 = Ifges(7,1) * t374 + Ifges(7,4) * t373 + Ifges(7,5) * t410;
t430 = mrSges(7,1) * t331 - mrSges(7,2) * t332 + Ifges(7,5) * t354 + Ifges(7,6) * t353 + Ifges(7,3) * t398 + t374 * t356 - t373 * t357;
t392 = Ifges(5,5) * qJD(4) + (t423 * Ifges(5,1) + t426 * Ifges(5,4)) * qJD(2);
t391 = Ifges(5,6) * qJD(4) + (t423 * Ifges(5,4) + Ifges(5,2) * t426) * qJD(2);
t371 = Ifges(6,1) * t396 + Ifges(6,4) * t395 - Ifges(6,5) * t440;
t370 = Ifges(6,4) * t396 + Ifges(6,2) * t395 - Ifges(6,6) * t440;
t369 = Ifges(6,5) * t396 + Ifges(6,6) * t395 - Ifges(6,3) * t440;
t355 = Ifges(7,5) * t374 + Ifges(7,6) * t373 + Ifges(7,3) * t410;
t322 = mrSges(7,2) * t337 - mrSges(7,3) * t331 + Ifges(7,1) * t354 + Ifges(7,4) * t353 + Ifges(7,5) * t398 + t373 * t355 - t410 * t356;
t321 = -mrSges(7,1) * t337 + mrSges(7,3) * t332 + Ifges(7,4) * t354 + Ifges(7,2) * t353 + Ifges(7,6) * t398 - t374 * t355 + t410 * t357;
t310 = mrSges(6,2) * t339 - mrSges(6,3) * t335 + Ifges(6,1) * t383 + Ifges(6,4) * t382 - Ifges(6,5) * t403 - pkin(9) * t320 - t422 * t321 + t425 * t322 + t395 * t369 + t370 * t440;
t309 = -mrSges(6,1) * t339 + mrSges(6,3) * t336 + Ifges(6,4) * t383 + Ifges(6,2) * t382 - Ifges(6,6) * t403 - pkin(5) * t432 + pkin(9) * t434 + t425 * t321 + t422 * t322 - t396 * t369 - t371 * t440;
t1 = [m(2) * t413 + (m(3) * t368 - t429 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t419 * t308 - t415 * t312) * t444 + (m(3) * t367 + qJDD(2) * mrSges(3,1) - t429 * mrSges(3,2) + t442) * t443 + t421 * (m(3) * t433 + t438); mrSges(3,1) * t367 - mrSges(3,2) * t368 + mrSges(4,1) * t350 - mrSges(4,2) * t351 + t423 * (mrSges(5,2) * t348 - mrSges(5,3) * t344 + Ifges(5,1) * t402 + Ifges(5,4) * t403 + Ifges(5,5) * qJDD(4) - qJ(5) * t316 - qJD(4) * t391 - t414 * t309 + t418 * t310) + t426 * (-mrSges(5,1) * t348 - mrSges(6,1) * t335 + mrSges(6,2) * t336 + mrSges(5,3) * t345 + Ifges(5,4) * t402 - Ifges(6,5) * t383 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t382 - pkin(4) * t316 - pkin(5) * t320 + qJD(4) * t392 - t396 * t370 + t395 * t371 - t430 + (Ifges(5,2) + Ifges(6,3)) * t403) + pkin(3) * t431 + pkin(8) * t436 + pkin(2) * t442 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t438; Ifges(5,5) * t402 + Ifges(5,6) * t403 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t344 - mrSges(5,2) * t345 + t414 * t310 + t418 * t309 - pkin(4) * t330 + qJ(5) * t435 + (t423 * t391 - t426 * t392) * qJD(2); t330; t430;];
tauJ  = t1;
