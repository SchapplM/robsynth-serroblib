% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRP2
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
% Datum: 2019-05-04 23:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:35:18
% EndTime: 2019-05-04 23:35:20
% DurationCPUTime: 1.33s
% Computational Cost: add. (7254->211), mult. (13196->260), div. (0->0), fcn. (8939->12), ass. (0->98)
t455 = Ifges(6,1) + Ifges(7,1);
t445 = Ifges(6,4) - Ifges(7,5);
t451 = -Ifges(6,5) - Ifges(7,4);
t454 = Ifges(6,2) + Ifges(7,3);
t443 = Ifges(6,6) - Ifges(7,6);
t414 = sin(qJ(5));
t415 = sin(qJ(4));
t434 = t415 * qJD(2);
t447 = cos(qJ(5));
t390 = -qJD(4) * t447 + t414 * t434;
t417 = cos(qJ(4));
t432 = qJD(2) * qJD(4);
t428 = t417 * t432;
t394 = qJDD(2) * t415 + t428;
t365 = -qJD(5) * t390 + qJDD(4) * t414 + t394 * t447;
t391 = qJD(4) * t414 + t434 * t447;
t369 = mrSges(7,1) * t390 - mrSges(7,3) * t391;
t409 = sin(pkin(10));
t412 = cos(pkin(10));
t397 = -g(1) * t412 - g(2) * t409;
t407 = -g(3) + qJDD(1);
t416 = sin(qJ(2));
t418 = cos(qJ(2));
t410 = sin(pkin(6));
t440 = t410 * t418;
t396 = g(1) * t409 - g(2) * t412;
t413 = cos(pkin(6));
t442 = t396 * t413;
t352 = -t397 * t416 + t407 * t440 + t418 * t442;
t350 = qJDD(2) * pkin(2) + t352;
t441 = t410 * t416;
t353 = t397 * t418 + t407 * t441 + t416 * t442;
t420 = qJD(2) ^ 2;
t351 = -pkin(2) * t420 + t353;
t408 = sin(pkin(11));
t411 = cos(pkin(11));
t346 = t350 * t408 + t351 * t411;
t344 = -pkin(3) * t420 + qJDD(2) * pkin(8) + t346;
t425 = -t396 * t410 + t407 * t413;
t377 = qJDD(3) + t425;
t340 = t344 * t417 + t377 * t415;
t393 = (-pkin(4) * t417 - pkin(9) * t415) * qJD(2);
t419 = qJD(4) ^ 2;
t433 = t417 * qJD(2);
t336 = -pkin(4) * t419 + qJDD(4) * pkin(9) + t393 * t433 + t340;
t345 = t411 * t350 - t351 * t408;
t343 = -qJDD(2) * pkin(3) - t420 * pkin(8) - t345;
t429 = t415 * t432;
t395 = qJDD(2) * t417 - t429;
t338 = (-t394 - t428) * pkin(9) + (-t395 + t429) * pkin(4) + t343;
t331 = -t336 * t414 + t338 * t447;
t368 = pkin(5) * t390 - qJ(6) * t391;
t388 = qJDD(5) - t395;
t403 = qJD(5) - t433;
t402 = t403 ^ 2;
t330 = -pkin(5) * t388 - qJ(6) * t402 + t368 * t391 + qJDD(6) - t331;
t376 = -mrSges(7,2) * t390 + mrSges(7,3) * t403;
t424 = -m(7) * t330 + mrSges(7,1) * t388 + t376 * t403;
t326 = t365 * mrSges(7,2) + t391 * t369 - t424;
t332 = t336 * t447 + t338 * t414;
t329 = -t402 * pkin(5) + t388 * qJ(6) + 0.2e1 * qJD(6) * t403 - t390 * t368 + t332;
t364 = qJD(5) * t391 - qJDD(4) * t447 + t394 * t414;
t375 = -mrSges(7,1) * t403 + mrSges(7,2) * t391;
t431 = m(7) * t329 + mrSges(7,3) * t388 + t375 * t403;
t437 = t390 * t454 - t391 * t445 - t403 * t443;
t448 = t390 * t445 - t391 * t455 + t403 * t451;
t450 = -Ifges(6,3) - Ifges(7,2);
t453 = -t451 * t365 - t448 * t390 - t443 * t364 - t450 * t388 + mrSges(6,1) * t331 - mrSges(7,1) * t330 - mrSges(6,2) * t332 + mrSges(7,3) * t329 - pkin(5) * t326 + qJ(6) * (-t364 * mrSges(7,2) - t390 * t369 + t431) - t437 * t391;
t446 = -mrSges(6,3) - mrSges(7,2);
t392 = (-mrSges(5,1) * t417 + mrSges(5,2) * t415) * qJD(2);
t398 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t434;
t374 = mrSges(6,1) * t403 - mrSges(6,3) * t391;
t435 = -mrSges(6,1) * t390 - mrSges(6,2) * t391 - t369;
t321 = m(6) * t332 - t388 * mrSges(6,2) + t364 * t446 - t403 * t374 + t390 * t435 + t431;
t373 = -mrSges(6,2) * t403 - mrSges(6,3) * t390;
t322 = m(6) * t331 + t388 * mrSges(6,1) + t365 * t446 + t403 * t373 + t391 * t435 + t424;
t426 = t321 * t447 - t322 * t414;
t316 = m(5) * t340 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t395 - qJD(4) * t398 + t392 * t433 + t426;
t339 = -t344 * t415 + t417 * t377;
t399 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t433;
t335 = -qJDD(4) * pkin(4) - t419 * pkin(9) + t393 * t434 - t339;
t333 = -0.2e1 * qJD(6) * t391 + (t390 * t403 - t365) * qJ(6) + (t391 * t403 + t364) * pkin(5) + t335;
t327 = m(7) * t333 + mrSges(7,1) * t364 - mrSges(7,3) * t365 - t375 * t391 + t376 * t390;
t421 = -m(6) * t335 - mrSges(6,1) * t364 - mrSges(6,2) * t365 - t373 * t390 - t374 * t391 - t327;
t324 = m(5) * t339 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t394 + qJD(4) * t399 - t392 * t434 + t421;
t427 = t316 * t417 - t324 * t415;
t311 = m(4) * t346 - mrSges(4,1) * t420 - qJDD(2) * mrSges(4,2) + t427;
t319 = t321 * t414 + t322 * t447;
t422 = -m(5) * t343 + mrSges(5,1) * t395 - mrSges(5,2) * t394 - t398 * t434 + t399 * t433 - t319;
t313 = m(4) * t345 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t420 + t422;
t439 = t311 * t408 + t313 * t411;
t438 = t390 * t443 + t391 * t451 + t403 * t450;
t430 = m(4) * t377 + t316 * t415 + t324 * t417;
t383 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t415 + Ifges(5,4) * t417) * qJD(2);
t382 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t415 + Ifges(5,2) * t417) * qJD(2);
t318 = mrSges(6,2) * t335 + mrSges(7,2) * t330 - mrSges(6,3) * t331 - mrSges(7,3) * t333 - qJ(6) * t327 - t364 * t445 + t365 * t455 - t388 * t451 + t390 * t438 + t403 * t437;
t317 = -mrSges(6,1) * t335 - mrSges(7,1) * t333 + mrSges(7,2) * t329 + mrSges(6,3) * t332 - pkin(5) * t327 - t364 * t454 + t365 * t445 + t388 * t443 + t391 * t438 - t403 * t448;
t1 = [m(2) * t407 + (m(3) * t353 - mrSges(3,1) * t420 - qJDD(2) * mrSges(3,2) + t311 * t411 - t313 * t408) * t441 + (m(3) * t352 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t420 + t439) * t440 + t413 * (m(3) * t425 + t430); mrSges(3,1) * t352 - mrSges(3,2) * t353 + mrSges(4,1) * t345 - mrSges(4,2) * t346 + t415 * (mrSges(5,2) * t343 - mrSges(5,3) * t339 + Ifges(5,1) * t394 + Ifges(5,4) * t395 + Ifges(5,5) * qJDD(4) - pkin(9) * t319 - qJD(4) * t382 - t317 * t414 + t318 * t447) + t417 * (-mrSges(5,1) * t343 + mrSges(5,3) * t340 + Ifges(5,4) * t394 + Ifges(5,2) * t395 + Ifges(5,6) * qJDD(4) - pkin(4) * t319 + qJD(4) * t383 - t453) + pkin(3) * t422 + pkin(8) * t427 + pkin(2) * t439 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t430; Ifges(5,5) * t394 + Ifges(5,6) * t395 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t339 - mrSges(5,2) * t340 + t414 * t318 + t447 * t317 + pkin(4) * t421 + pkin(9) * t426 + (t382 * t415 - t383 * t417) * qJD(2); t453; t326;];
tauJ  = t1;
