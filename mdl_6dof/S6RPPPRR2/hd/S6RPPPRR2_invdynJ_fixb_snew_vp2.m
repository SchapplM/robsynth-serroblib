% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPPRR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 13:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:36:01
% EndTime: 2019-05-05 13:36:02
% DurationCPUTime: 1.33s
% Computational Cost: add. (7982->205), mult. (16643->258), div. (0->0), fcn. (10589->10), ass. (0->98)
t404 = sin(pkin(10));
t397 = t404 ^ 2;
t406 = cos(pkin(10));
t398 = t406 ^ 2;
t439 = -t397 - t398;
t448 = t439 * mrSges(5,3);
t410 = sin(qJ(1));
t413 = cos(qJ(1));
t430 = t410 * g(1) - t413 * g(2);
t387 = qJDD(1) * pkin(1) + t430;
t415 = qJD(1) ^ 2;
t427 = -t413 * g(1) - t410 * g(2);
t388 = -t415 * pkin(1) + t427;
t405 = sin(pkin(9));
t407 = cos(pkin(9));
t373 = t407 * t387 - t405 * t388;
t421 = -t415 * qJ(3) + qJDD(3) - t373;
t442 = -pkin(2) - qJ(4);
t447 = -(2 * qJD(1) * qJD(4)) + t442 * qJDD(1) + t421;
t374 = t405 * t387 + t407 * t388;
t446 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t374;
t409 = sin(qJ(5));
t412 = cos(qJ(5));
t425 = t404 * t412 + t406 * t409;
t380 = t425 * qJD(1);
t361 = t415 * pkin(2) - t446;
t445 = -m(4) * t361 + t415 * mrSges(4,2) + qJDD(1) * mrSges(4,3);
t424 = -t404 * t409 + t406 * t412;
t381 = t424 * qJD(1);
t436 = t381 * qJD(5);
t371 = -t425 * qJDD(1) - t436;
t443 = pkin(4) * t415;
t401 = -g(3) + qJDD(2);
t433 = t406 * qJDD(1);
t440 = t447 * t406;
t340 = -pkin(7) * t433 + (-t406 * t443 - t401) * t404 + t440;
t352 = t406 * t401 + t447 * t404;
t434 = t404 * qJDD(1);
t341 = -pkin(7) * t434 - t397 * t443 + t352;
t337 = t409 * t340 + t412 * t341;
t367 = t380 * mrSges(6,1) + t381 * mrSges(6,2);
t378 = qJD(5) * mrSges(6,1) - t381 * mrSges(6,3);
t370 = t380 * pkin(5) - t381 * pkin(8);
t414 = qJD(5) ^ 2;
t334 = -t414 * pkin(5) + qJDD(5) * pkin(8) - t380 * t370 + t337;
t423 = qJDD(4) + t446;
t350 = pkin(4) * t434 + (t439 * pkin(7) + t442) * t415 + t423;
t437 = t380 * qJD(5);
t372 = t424 * qJDD(1) - t437;
t335 = (-t372 + t437) * pkin(8) + (-t371 + t436) * pkin(5) + t350;
t408 = sin(qJ(6));
t411 = cos(qJ(6));
t331 = -t408 * t334 + t411 * t335;
t375 = t411 * qJD(5) - t408 * t381;
t348 = t375 * qJD(6) + t408 * qJDD(5) + t411 * t372;
t376 = t408 * qJD(5) + t411 * t381;
t353 = -mrSges(7,1) * t375 + mrSges(7,2) * t376;
t379 = qJD(6) + t380;
t358 = -t379 * mrSges(7,2) + t375 * mrSges(7,3);
t369 = qJDD(6) - t371;
t329 = m(7) * t331 + t369 * mrSges(7,1) - t348 * mrSges(7,3) - t376 * t353 + t379 * t358;
t332 = t411 * t334 + t408 * t335;
t347 = -t376 * qJD(6) + t411 * qJDD(5) - t408 * t372;
t359 = t379 * mrSges(7,1) - t376 * mrSges(7,3);
t330 = m(7) * t332 - t369 * mrSges(7,2) + t347 * mrSges(7,3) + t375 * t353 - t379 * t359;
t428 = -t408 * t329 + t411 * t330;
t320 = m(6) * t337 - qJDD(5) * mrSges(6,2) + t371 * mrSges(6,3) - qJD(5) * t378 - t380 * t367 + t428;
t336 = t412 * t340 - t409 * t341;
t377 = -qJD(5) * mrSges(6,2) - t380 * mrSges(6,3);
t333 = -qJDD(5) * pkin(5) - t414 * pkin(8) + t381 * t370 - t336;
t420 = -m(7) * t333 + t347 * mrSges(7,1) - t348 * mrSges(7,2) + t375 * t358 - t376 * t359;
t325 = m(6) * t336 + qJDD(5) * mrSges(6,1) - t372 * mrSges(6,3) + qJD(5) * t377 - t381 * t367 + t420;
t441 = t409 * t320 + t412 * t325;
t321 = t411 * t329 + t408 * t330;
t429 = t412 * t320 - t409 * t325;
t351 = -t404 * t401 + t440;
t422 = -qJDD(1) * mrSges(5,3) - t415 * (t404 * mrSges(5,1) + t406 * mrSges(5,2));
t316 = m(5) * t351 + t422 * t406 + t441;
t317 = m(5) * t352 + t422 * t404 + t429;
t426 = t406 * t316 + t404 * t317;
t419 = m(6) * t350 - t371 * mrSges(6,1) + t372 * mrSges(6,2) + t380 * t377 + t381 * t378 + t321;
t362 = -qJDD(1) * pkin(2) + t421;
t313 = m(4) * t362 + qJDD(1) * mrSges(4,2) - t415 * mrSges(4,3) + t426;
t357 = t442 * t415 + t423;
t418 = m(5) * t357 + mrSges(5,1) * t434 + mrSges(5,2) * t433 + t419;
t343 = Ifges(7,4) * t376 + Ifges(7,2) * t375 + Ifges(7,6) * t379;
t344 = Ifges(7,1) * t376 + Ifges(7,4) * t375 + Ifges(7,5) * t379;
t417 = mrSges(7,1) * t331 - mrSges(7,2) * t332 + Ifges(7,5) * t348 + Ifges(7,6) * t347 + Ifges(7,3) * t369 + t376 * t343 - t375 * t344;
t416 = t415 * t448 + t418;
t365 = Ifges(6,1) * t381 - Ifges(6,4) * t380 + Ifges(6,5) * qJD(5);
t364 = Ifges(6,4) * t381 - Ifges(6,2) * t380 + Ifges(6,6) * qJD(5);
t363 = Ifges(6,5) * t381 - Ifges(6,6) * t380 + Ifges(6,3) * qJD(5);
t342 = Ifges(7,5) * t376 + Ifges(7,6) * t375 + Ifges(7,3) * t379;
t323 = mrSges(7,2) * t333 - mrSges(7,3) * t331 + Ifges(7,1) * t348 + Ifges(7,4) * t347 + Ifges(7,5) * t369 + t375 * t342 - t379 * t343;
t322 = -mrSges(7,1) * t333 + mrSges(7,3) * t332 + Ifges(7,4) * t348 + Ifges(7,2) * t347 + Ifges(7,6) * t369 - t376 * t342 + t379 * t344;
t315 = -mrSges(6,1) * t350 + mrSges(6,3) * t337 + Ifges(6,4) * t372 + Ifges(6,2) * t371 + Ifges(6,6) * qJDD(5) - pkin(5) * t321 + qJD(5) * t365 - t381 * t363 - t417;
t314 = mrSges(6,2) * t350 - mrSges(6,3) * t336 + Ifges(6,1) * t372 + Ifges(6,4) * t371 + Ifges(6,5) * qJDD(5) - pkin(8) * t321 - qJD(5) * t364 - t408 * t322 + t411 * t323 - t380 * t363;
t1 = [pkin(1) * (t405 * (m(3) * t374 + t418 + t445) + t407 * (m(3) * t373 - t313) + (t405 * (-mrSges(3,1) + t448) - t407 * mrSges(3,2)) * t415) - mrSges(2,2) * t427 + mrSges(2,1) * t430 - pkin(2) * t313 + qJ(3) * (t416 + t445) - t404 * (-mrSges(5,1) * t357 + mrSges(5,3) * t352 - pkin(4) * t419 + pkin(7) * t429 + t409 * t314 + t412 * t315) - qJ(4) * t426 + mrSges(3,1) * t373 - mrSges(3,2) * t374 + t406 * (mrSges(5,2) * t357 - mrSges(5,3) * t351 - pkin(7) * t441 + t412 * t314 - t409 * t315) + mrSges(4,2) * t362 - mrSges(4,3) * t361 + (pkin(1) * (t407 * mrSges(3,1) - t405 * mrSges(3,2)) + Ifges(5,1) * t398 + Ifges(2,3) + Ifges(3,3) + Ifges(4,1) + (-0.2e1 * Ifges(5,4) * t406 + Ifges(5,2) * t404) * t404) * qJDD(1); -t404 * t316 + t406 * t317 + (m(3) + m(4)) * t401; t313; t416; mrSges(6,1) * t336 - mrSges(6,2) * t337 + Ifges(6,5) * t372 + Ifges(6,6) * t371 + Ifges(6,3) * qJDD(5) + pkin(5) * t420 + pkin(8) * t428 + t411 * t322 + t408 * t323 + t381 * t364 + t380 * t365; t417;];
tauJ  = t1;
