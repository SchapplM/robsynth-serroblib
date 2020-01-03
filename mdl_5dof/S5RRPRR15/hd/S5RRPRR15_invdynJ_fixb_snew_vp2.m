% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR15_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:27
% EndTime: 2019-12-31 20:41:30
% DurationCPUTime: 1.93s
% Computational Cost: add. (10314->261), mult. (21203->317), div. (0->0), fcn. (12267->8), ass. (0->107)
t461 = -2 * qJD(3);
t460 = Ifges(3,1) + Ifges(4,2);
t455 = Ifges(3,4) + Ifges(4,6);
t454 = Ifges(3,5) - Ifges(4,4);
t459 = Ifges(3,2) + Ifges(4,3);
t453 = Ifges(3,6) - Ifges(4,5);
t458 = (Ifges(3,3) + Ifges(4,1));
t429 = qJD(1) ^ 2;
t423 = sin(qJ(1));
t427 = cos(qJ(1));
t440 = -g(1) * t427 - g(2) * t423;
t388 = -pkin(1) * t429 + qJDD(1) * pkin(6) + t440;
t422 = sin(qJ(2));
t426 = cos(qJ(2));
t373 = -g(3) * t422 + t426 * t388;
t396 = (-t426 * pkin(2) - t422 * qJ(3)) * qJD(1);
t428 = qJD(2) ^ 2;
t448 = qJD(1) * t426;
t354 = pkin(2) * t428 - qJDD(2) * qJ(3) + (qJD(2) * t461) - t396 * t448 - t373;
t457 = pkin(6) * t429;
t456 = mrSges(3,1) - mrSges(4,2);
t447 = qJD(1) * qJD(2);
t444 = t422 * t447;
t400 = qJDD(1) * t426 - t444;
t413 = t422 * qJD(1);
t407 = pkin(3) * t413 - (qJD(2) * pkin(7));
t419 = t426 ^ 2;
t445 = t426 * t447;
t399 = qJDD(1) * t422 + t445;
t443 = g(1) * t423 - t427 * g(2);
t438 = -qJDD(1) * pkin(1) - t443;
t433 = pkin(2) * t444 + t413 * t461 + (-t399 - t445) * qJ(3) + t438;
t341 = -t407 * t413 + (-pkin(3) * t419 - pkin(6)) * t429 + (-pkin(2) - pkin(7)) * t400 + t433;
t372 = -t426 * g(3) - t422 * t388;
t355 = -qJDD(2) * pkin(2) - qJ(3) * t428 + t396 * t413 + qJDD(3) - t372;
t346 = (-t422 * t426 * t429 - qJDD(2)) * pkin(7) + (t399 - t445) * pkin(3) + t355;
t421 = sin(qJ(4));
t425 = cos(qJ(4));
t330 = -t341 * t421 + t425 * t346;
t394 = -qJD(2) * t421 - t425 * t448;
t366 = qJD(4) * t394 + qJDD(2) * t425 - t400 * t421;
t393 = qJDD(4) + t399;
t395 = qJD(2) * t425 - t421 * t448;
t410 = t413 + qJD(4);
t327 = (t394 * t410 - t366) * pkin(8) + (t394 * t395 + t393) * pkin(4) + t330;
t331 = t425 * t341 + t421 * t346;
t365 = -qJD(4) * t395 - qJDD(2) * t421 - t400 * t425;
t374 = pkin(4) * t410 - pkin(8) * t395;
t392 = t394 ^ 2;
t328 = -pkin(4) * t392 + pkin(8) * t365 - t374 * t410 + t331;
t420 = sin(qJ(5));
t424 = cos(qJ(5));
t325 = t327 * t424 - t328 * t420;
t367 = t394 * t424 - t395 * t420;
t339 = qJD(5) * t367 + t365 * t420 + t366 * t424;
t368 = t394 * t420 + t395 * t424;
t351 = -mrSges(6,1) * t367 + mrSges(6,2) * t368;
t408 = qJD(5) + t410;
t356 = -mrSges(6,2) * t408 + mrSges(6,3) * t367;
t390 = qJDD(5) + t393;
t322 = m(6) * t325 + mrSges(6,1) * t390 - mrSges(6,3) * t339 - t351 * t368 + t356 * t408;
t326 = t327 * t420 + t328 * t424;
t338 = -qJD(5) * t368 + t365 * t424 - t366 * t420;
t357 = mrSges(6,1) * t408 - mrSges(6,3) * t368;
t323 = m(6) * t326 - mrSges(6,2) * t390 + mrSges(6,3) * t338 + t351 * t367 - t357 * t408;
t316 = t424 * t322 + t420 * t323;
t452 = (t458 * qJD(2)) + (t422 * t454 + t426 * t453) * qJD(1);
t451 = t454 * qJD(2) + (t422 * t460 + t426 * t455) * qJD(1);
t450 = -t453 * qJD(2) + (-t422 * t455 - t459 * t426) * qJD(1);
t405 = -mrSges(4,1) * t448 - (qJD(2) * mrSges(4,3));
t449 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t448 - t405;
t369 = -mrSges(5,1) * t394 + mrSges(5,2) * t395;
t370 = -mrSges(5,2) * t410 + mrSges(5,3) * t394;
t313 = m(5) * t330 + mrSges(5,1) * t393 - mrSges(5,3) * t366 - t369 * t395 + t370 * t410 + t316;
t371 = mrSges(5,1) * t410 - mrSges(5,3) * t395;
t441 = -t322 * t420 + t424 * t323;
t314 = m(5) * t331 - mrSges(5,2) * t393 + mrSges(5,3) * t365 + t369 * t394 - t371 * t410 + t441;
t442 = -t313 * t421 + t425 * t314;
t311 = t425 * t313 + t421 * t314;
t352 = -pkin(2) * t400 + t433 - t457;
t439 = m(4) * t352 + t442;
t345 = -pkin(7) * t419 * t429 + pkin(3) * t400 + qJD(2) * t407 - t354;
t333 = -pkin(4) * t365 - pkin(8) * t392 + t374 * t395 + t345;
t436 = m(6) * t333 - t338 * mrSges(6,1) + t339 * mrSges(6,2) - t367 * t356 + t368 * t357;
t435 = m(4) * t355 + t399 * mrSges(4,1) + t311;
t348 = Ifges(6,4) * t368 + Ifges(6,2) * t367 + Ifges(6,6) * t408;
t349 = Ifges(6,1) * t368 + Ifges(6,4) * t367 + Ifges(6,5) * t408;
t434 = mrSges(6,1) * t325 - mrSges(6,2) * t326 + Ifges(6,5) * t339 + Ifges(6,6) * t338 + Ifges(6,3) * t390 + t368 * t348 - t367 * t349;
t432 = -m(5) * t345 + t365 * mrSges(5,1) - t366 * mrSges(5,2) + t394 * t370 - t395 * t371 - t436;
t359 = Ifges(5,4) * t395 + Ifges(5,2) * t394 + Ifges(5,6) * t410;
t360 = Ifges(5,1) * t395 + Ifges(5,4) * t394 + Ifges(5,5) * t410;
t431 = mrSges(5,1) * t330 - mrSges(5,2) * t331 + Ifges(5,5) * t366 + Ifges(5,6) * t365 + Ifges(5,3) * t393 + pkin(4) * t316 + t395 * t359 - t394 * t360 + t434;
t397 = (t426 * mrSges(4,2) - t422 * mrSges(4,3)) * qJD(1);
t406 = mrSges(4,1) * t413 + qJD(2) * mrSges(4,2);
t430 = -m(4) * t354 + qJDD(2) * mrSges(4,3) + qJD(2) * t406 + t397 * t448 - t432;
t403 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t413;
t398 = (-t426 * mrSges(3,1) + t422 * mrSges(3,2)) * qJD(1);
t387 = t438 - t457;
t358 = Ifges(5,5) * t395 + Ifges(5,6) * t394 + Ifges(5,3) * t410;
t347 = Ifges(6,5) * t368 + Ifges(6,6) * t367 + Ifges(6,3) * t408;
t318 = mrSges(6,2) * t333 - mrSges(6,3) * t325 + Ifges(6,1) * t339 + Ifges(6,4) * t338 + Ifges(6,5) * t390 + t347 * t367 - t348 * t408;
t317 = -mrSges(6,1) * t333 + mrSges(6,3) * t326 + Ifges(6,4) * t339 + Ifges(6,2) * t338 + Ifges(6,6) * t390 - t347 * t368 + t349 * t408;
t310 = qJDD(2) * mrSges(4,2) + qJD(2) * t405 + t397 * t413 + t435;
t309 = mrSges(4,2) * t400 - mrSges(4,3) * t399 + (t405 * t426 - t406 * t422) * qJD(1) + t439;
t308 = mrSges(5,2) * t345 - mrSges(5,3) * t330 + Ifges(5,1) * t366 + Ifges(5,4) * t365 + Ifges(5,5) * t393 - pkin(8) * t316 - t317 * t420 + t318 * t424 + t358 * t394 - t359 * t410;
t307 = -mrSges(5,1) * t345 + mrSges(5,3) * t331 + Ifges(5,4) * t366 + Ifges(5,2) * t365 + Ifges(5,6) * t393 - pkin(4) * t436 + pkin(8) * t441 + t424 * t317 + t420 * t318 - t395 * t358 + t410 * t360;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t443 - mrSges(2,2) * t440 + t422 * (mrSges(4,1) * t355 + mrSges(3,2) * t387 - mrSges(3,3) * t372 - mrSges(4,3) * t352 + pkin(3) * t311 - qJ(3) * t309 + t450 * qJD(2) + t454 * qJDD(2) + t460 * t399 + t455 * t400 + t452 * t448 + t431) + t426 * (-mrSges(3,1) * t387 - mrSges(4,1) * t354 + mrSges(4,2) * t352 + mrSges(3,3) * t373 - pkin(2) * t309 - pkin(3) * t432 - pkin(7) * t442 + t451 * qJD(2) + t453 * qJDD(2) - t425 * t307 - t421 * t308 + t455 * t399 + t459 * t400 - t452 * t413) + pkin(1) * (-m(3) * t387 + t456 * t400 + (-mrSges(3,2) + mrSges(4,3)) * t399 + (t449 * t426 + (-t403 + t406) * t422) * qJD(1) - t439) + pkin(6) * (t426 * (t430 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t400 - qJD(2) * t403 + m(3) * t373 + t398 * t448) + (-m(3) * t372 + t399 * mrSges(3,3) - t456 * qJDD(2) - t449 * qJD(2) + (t397 + t398) * t413 + t435) * t422); mrSges(3,1) * t372 - mrSges(3,2) * t373 + mrSges(4,2) * t355 - mrSges(4,3) * t354 + t425 * t308 - t421 * t307 - pkin(7) * t311 - pkin(2) * t310 + qJ(3) * t430 + (qJ(3) * mrSges(4,1) + t453) * t400 + t454 * t399 + t458 * qJDD(2) + (-t450 * t422 - t451 * t426) * qJD(1); t310; t431; t434;];
tauJ = t1;
