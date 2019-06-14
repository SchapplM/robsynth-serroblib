% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:29:27
% EndTime: 2019-05-04 23:29:29
% DurationCPUTime: 1.37s
% Computational Cost: add. (7340->213), mult. (13454->260), div. (0->0), fcn. (9142->12), ass. (0->99)
t456 = Ifges(6,1) + Ifges(7,1);
t448 = Ifges(6,4) + Ifges(7,4);
t447 = Ifges(6,5) + Ifges(7,5);
t455 = Ifges(6,2) + Ifges(7,2);
t446 = Ifges(6,6) + Ifges(7,6);
t415 = sin(qJ(5));
t418 = cos(qJ(5));
t416 = sin(qJ(4));
t437 = t416 * qJD(2);
t394 = t418 * qJD(4) - t415 * t437;
t419 = cos(qJ(4));
t435 = qJD(2) * qJD(4);
t430 = t419 * t435;
t398 = t416 * qJDD(2) + t430;
t370 = t394 * qJD(5) + t415 * qJDD(4) + t418 * t398;
t395 = t415 * qJD(4) + t418 * t437;
t372 = -t394 * mrSges(7,1) + t395 * mrSges(7,2);
t410 = sin(pkin(10));
t413 = cos(pkin(10));
t401 = -t413 * g(1) - t410 * g(2);
t408 = -g(3) + qJDD(1);
t417 = sin(qJ(2));
t420 = cos(qJ(2));
t411 = sin(pkin(6));
t443 = t411 * t420;
t400 = t410 * g(1) - t413 * g(2);
t414 = cos(pkin(6));
t445 = t400 * t414;
t356 = -t417 * t401 + t408 * t443 + t420 * t445;
t354 = qJDD(2) * pkin(2) + t356;
t444 = t411 * t417;
t357 = t420 * t401 + t408 * t444 + t417 * t445;
t422 = qJD(2) ^ 2;
t355 = -t422 * pkin(2) + t357;
t409 = sin(pkin(11));
t412 = cos(pkin(11));
t349 = t409 * t354 + t412 * t355;
t347 = -t422 * pkin(3) + qJDD(2) * pkin(8) + t349;
t427 = -t411 * t400 + t414 * t408;
t381 = qJDD(3) + t427;
t343 = t419 * t347 + t416 * t381;
t397 = (-t419 * pkin(4) - t416 * pkin(9)) * qJD(2);
t421 = qJD(4) ^ 2;
t436 = t419 * qJD(2);
t338 = -t421 * pkin(4) + qJDD(4) * pkin(9) + t397 * t436 + t343;
t348 = t412 * t354 - t409 * t355;
t346 = -qJDD(2) * pkin(3) - t422 * pkin(8) - t348;
t431 = t416 * t435;
t399 = t419 * qJDD(2) - t431;
t341 = (-t398 - t430) * pkin(9) + (-t399 + t431) * pkin(4) + t346;
t333 = -t415 * t338 + t418 * t341;
t392 = qJDD(5) - t399;
t406 = qJD(5) - t436;
t329 = -0.2e1 * qJD(6) * t395 + (t394 * t406 - t370) * qJ(6) + (t394 * t395 + t392) * pkin(5) + t333;
t376 = -t406 * mrSges(7,2) + t394 * mrSges(7,3);
t434 = m(7) * t329 + t392 * mrSges(7,1) + t406 * t376;
t327 = -t370 * mrSges(7,3) - t395 * t372 + t434;
t334 = t418 * t338 + t415 * t341;
t369 = -t395 * qJD(5) + t418 * qJDD(4) - t415 * t398;
t378 = t406 * pkin(5) - t395 * qJ(6);
t391 = t394 ^ 2;
t331 = -t391 * pkin(5) + t369 * qJ(6) + 0.2e1 * qJD(6) * t394 - t406 * t378 + t334;
t440 = -t455 * t394 - t448 * t395 - t446 * t406;
t450 = t448 * t394 + t456 * t395 + t447 * t406;
t452 = Ifges(6,3) + Ifges(7,3);
t454 = mrSges(6,1) * t333 + mrSges(7,1) * t329 - mrSges(6,2) * t334 - mrSges(7,2) * t331 + pkin(5) * t327 + t446 * t369 + t447 * t370 + t452 * t392 - t450 * t394 - t440 * t395;
t449 = -mrSges(6,2) - mrSges(7,2);
t396 = (-t419 * mrSges(5,1) + t416 * mrSges(5,2)) * qJD(2);
t402 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t437;
t373 = -t394 * mrSges(6,1) + t395 * mrSges(6,2);
t377 = -t406 * mrSges(6,2) + t394 * mrSges(6,3);
t322 = m(6) * t333 + t392 * mrSges(6,1) + t406 * t377 + (-t372 - t373) * t395 + (-mrSges(6,3) - mrSges(7,3)) * t370 + t434;
t433 = m(7) * t331 + t369 * mrSges(7,3) + t394 * t372;
t379 = t406 * mrSges(7,1) - t395 * mrSges(7,3);
t438 = -t406 * mrSges(6,1) + t395 * mrSges(6,3) - t379;
t324 = m(6) * t334 + t369 * mrSges(6,3) + t394 * t373 + t392 * t449 + t438 * t406 + t433;
t428 = -t415 * t322 + t418 * t324;
t319 = m(5) * t343 - qJDD(4) * mrSges(5,2) + t399 * mrSges(5,3) - qJD(4) * t402 + t396 * t436 + t428;
t342 = -t416 * t347 + t419 * t381;
t403 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t436;
t337 = -qJDD(4) * pkin(4) - t421 * pkin(9) + t397 * t437 - t342;
t335 = -t369 * pkin(5) - t391 * qJ(6) + t395 * t378 + qJDD(6) + t337;
t426 = -m(7) * t335 + t369 * mrSges(7,1) + t394 * t376;
t423 = -m(6) * t337 + t369 * mrSges(6,1) + t370 * t449 + t394 * t377 + t438 * t395 + t426;
t326 = m(5) * t342 + qJDD(4) * mrSges(5,1) - t398 * mrSges(5,3) + qJD(4) * t403 - t396 * t437 + t423;
t429 = t419 * t319 - t416 * t326;
t313 = m(4) * t349 - t422 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t429;
t321 = t418 * t322 + t415 * t324;
t424 = -m(5) * t346 + t399 * mrSges(5,1) - t398 * mrSges(5,2) - t402 * t437 + t403 * t436 - t321;
t316 = m(4) * t348 + qJDD(2) * mrSges(4,1) - t422 * mrSges(4,2) + t424;
t442 = t409 * t313 + t412 * t316;
t441 = -t446 * t394 - t395 * t447 - t452 * t406;
t432 = m(4) * t381 + t416 * t319 + t419 * t326;
t387 = Ifges(5,5) * qJD(4) + (t416 * Ifges(5,1) + t419 * Ifges(5,4)) * qJD(2);
t386 = Ifges(5,6) * qJD(4) + (t416 * Ifges(5,4) + t419 * Ifges(5,2)) * qJD(2);
t332 = t370 * mrSges(7,2) + t395 * t379 - t426;
t320 = mrSges(6,2) * t337 + mrSges(7,2) * t335 - mrSges(6,3) * t333 - mrSges(7,3) * t329 - qJ(6) * t327 + t448 * t369 + t456 * t370 + t447 * t392 - t441 * t394 + t440 * t406;
t314 = -mrSges(6,1) * t337 + mrSges(6,3) * t334 - mrSges(7,1) * t335 + mrSges(7,3) * t331 - pkin(5) * t332 + qJ(6) * t433 + (-qJ(6) * t379 + t450) * t406 + t441 * t395 + (-qJ(6) * mrSges(7,2) + t446) * t392 + t448 * t370 + t455 * t369;
t1 = [m(2) * t408 + (m(3) * t357 - t422 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t412 * t313 - t409 * t316) * t444 + (m(3) * t356 + qJDD(2) * mrSges(3,1) - t422 * mrSges(3,2) + t442) * t443 + t414 * (m(3) * t427 + t432); mrSges(3,1) * t356 - mrSges(3,2) * t357 + mrSges(4,1) * t348 - mrSges(4,2) * t349 + t416 * (mrSges(5,2) * t346 - mrSges(5,3) * t342 + Ifges(5,1) * t398 + Ifges(5,4) * t399 + Ifges(5,5) * qJDD(4) - pkin(9) * t321 - qJD(4) * t386 - t415 * t314 + t418 * t320) + t419 * (-mrSges(5,1) * t346 + mrSges(5,3) * t343 + Ifges(5,4) * t398 + Ifges(5,2) * t399 + Ifges(5,6) * qJDD(4) - pkin(4) * t321 + qJD(4) * t387 - t454) + pkin(3) * t424 + pkin(8) * t429 + pkin(2) * t442 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t432; Ifges(5,5) * t398 + Ifges(5,6) * t399 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t342 - mrSges(5,2) * t343 + t415 * t320 + t418 * t314 + pkin(4) * t423 + pkin(9) * t428 + (t416 * t386 - t419 * t387) * qJD(2); t454; t332;];
tauJ  = t1;
