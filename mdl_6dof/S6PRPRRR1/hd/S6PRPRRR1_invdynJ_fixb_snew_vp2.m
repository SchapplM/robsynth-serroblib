% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:05:36
% EndTime: 2019-05-05 00:05:38
% DurationCPUTime: 1.64s
% Computational Cost: add. (15270->236), mult. (28754->307), div. (0->0), fcn. (20527->14), ass. (0->108)
t441 = sin(qJ(5));
t442 = sin(qJ(4));
t445 = cos(qJ(5));
t446 = cos(qJ(4));
t411 = (t441 * t442 - t445 * t446) * qJD(2);
t435 = sin(pkin(11));
t438 = cos(pkin(11));
t420 = g(1) * t435 - g(2) * t438;
t439 = cos(pkin(6));
t468 = t420 * t439;
t436 = sin(pkin(6));
t443 = sin(qJ(2));
t467 = t436 * t443;
t447 = cos(qJ(2));
t466 = t436 * t447;
t421 = -g(1) * t438 - g(2) * t435;
t433 = -g(3) + qJDD(1);
t392 = -t421 * t443 + t433 * t466 + t447 * t468;
t387 = qJDD(2) * pkin(2) + t392;
t393 = t447 * t421 + t433 * t467 + t443 * t468;
t448 = qJD(2) ^ 2;
t388 = -pkin(2) * t448 + t393;
t434 = sin(pkin(12));
t437 = cos(pkin(12));
t366 = t434 * t387 + t437 * t388;
t364 = -pkin(3) * t448 + qJDD(2) * pkin(8) + t366;
t456 = -t420 * t436 + t439 * t433;
t404 = qJDD(3) + t456;
t360 = -t442 * t364 + t446 * t404;
t462 = qJD(2) * qJD(4);
t460 = t446 * t462;
t418 = qJDD(2) * t442 + t460;
t357 = (-t418 + t460) * pkin(9) + (t442 * t446 * t448 + qJDD(4)) * pkin(4) + t360;
t361 = t446 * t364 + t442 * t404;
t419 = qJDD(2) * t446 - t442 * t462;
t464 = qJD(2) * t442;
t426 = qJD(4) * pkin(4) - pkin(9) * t464;
t432 = t446 ^ 2;
t358 = -pkin(4) * t432 * t448 + pkin(9) * t419 - qJD(4) * t426 + t361;
t353 = t441 * t357 + t445 * t358;
t412 = (t441 * t446 + t442 * t445) * qJD(2);
t380 = -qJD(5) * t412 - t418 * t441 + t419 * t445;
t395 = mrSges(6,1) * t411 + mrSges(6,2) * t412;
t431 = qJD(4) + qJD(5);
t403 = mrSges(6,1) * t431 - mrSges(6,3) * t412;
t430 = qJDD(4) + qJDD(5);
t396 = pkin(5) * t411 - pkin(10) * t412;
t429 = t431 ^ 2;
t350 = -pkin(5) * t429 + pkin(10) * t430 - t396 * t411 + t353;
t365 = t437 * t387 - t434 * t388;
t454 = -qJDD(2) * pkin(3) - t365;
t359 = -t419 * pkin(4) + t426 * t464 + (-pkin(9) * t432 - pkin(8)) * t448 + t454;
t381 = -qJD(5) * t411 + t418 * t445 + t419 * t441;
t354 = (t411 * t431 - t381) * pkin(10) + (t412 * t431 - t380) * pkin(5) + t359;
t440 = sin(qJ(6));
t444 = cos(qJ(6));
t347 = -t350 * t440 + t354 * t444;
t397 = -t412 * t440 + t431 * t444;
t369 = qJD(6) * t397 + t381 * t444 + t430 * t440;
t398 = t412 * t444 + t431 * t440;
t374 = -mrSges(7,1) * t397 + mrSges(7,2) * t398;
t379 = qJDD(6) - t380;
t405 = qJD(6) + t411;
t382 = -mrSges(7,2) * t405 + mrSges(7,3) * t397;
t344 = m(7) * t347 + mrSges(7,1) * t379 - mrSges(7,3) * t369 - t374 * t398 + t382 * t405;
t348 = t350 * t444 + t354 * t440;
t368 = -qJD(6) * t398 - t381 * t440 + t430 * t444;
t383 = mrSges(7,1) * t405 - mrSges(7,3) * t398;
t345 = m(7) * t348 - mrSges(7,2) * t379 + mrSges(7,3) * t368 + t374 * t397 - t383 * t405;
t457 = -t344 * t440 + t444 * t345;
t332 = m(6) * t353 - mrSges(6,2) * t430 + mrSges(6,3) * t380 - t395 * t411 - t403 * t431 + t457;
t352 = t357 * t445 - t358 * t441;
t402 = -mrSges(6,2) * t431 - mrSges(6,3) * t411;
t349 = -pkin(5) * t430 - pkin(10) * t429 + t396 * t412 - t352;
t453 = -m(7) * t349 + t368 * mrSges(7,1) - mrSges(7,2) * t369 + t397 * t382 - t383 * t398;
t340 = m(6) * t352 + mrSges(6,1) * t430 - mrSges(6,3) * t381 - t395 * t412 + t402 * t431 + t453;
t327 = t441 * t332 + t445 * t340;
t417 = (-mrSges(5,1) * t446 + mrSges(5,2) * t442) * qJD(2);
t463 = qJD(2) * t446;
t423 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t463;
t325 = m(5) * t360 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t418 + qJD(4) * t423 - t417 * t464 + t327;
t422 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t464;
t458 = t445 * t332 - t340 * t441;
t326 = m(5) * t361 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t419 - qJD(4) * t422 + t417 * t463 + t458;
t459 = -t325 * t442 + t446 * t326;
t319 = m(4) * t366 - mrSges(4,1) * t448 - qJDD(2) * mrSges(4,2) + t459;
t363 = -t448 * pkin(8) + t454;
t334 = t444 * t344 + t440 * t345;
t452 = m(6) * t359 - t380 * mrSges(6,1) + mrSges(6,2) * t381 + t411 * t402 + t403 * t412 + t334;
t449 = -m(5) * t363 + t419 * mrSges(5,1) - mrSges(5,2) * t418 - t422 * t464 + t423 * t463 - t452;
t329 = m(4) * t365 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t448 + t449;
t465 = t434 * t319 + t437 * t329;
t461 = m(4) * t404 + t446 * t325 + t442 * t326;
t370 = Ifges(7,5) * t398 + Ifges(7,6) * t397 + Ifges(7,3) * t405;
t372 = Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t405;
t337 = -mrSges(7,1) * t349 + mrSges(7,3) * t348 + Ifges(7,4) * t369 + Ifges(7,2) * t368 + Ifges(7,6) * t379 - t370 * t398 + t372 * t405;
t371 = Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t405;
t338 = mrSges(7,2) * t349 - mrSges(7,3) * t347 + Ifges(7,1) * t369 + Ifges(7,4) * t368 + Ifges(7,5) * t379 + t370 * t397 - t371 * t405;
t390 = Ifges(6,4) * t412 - Ifges(6,2) * t411 + Ifges(6,6) * t431;
t391 = Ifges(6,1) * t412 - Ifges(6,4) * t411 + Ifges(6,5) * t431;
t451 = mrSges(6,1) * t352 - mrSges(6,2) * t353 + Ifges(6,5) * t381 + Ifges(6,6) * t380 + Ifges(6,3) * t430 + pkin(5) * t453 + pkin(10) * t457 + t444 * t337 + t440 * t338 + t412 * t390 + t411 * t391;
t450 = mrSges(7,1) * t347 - mrSges(7,2) * t348 + Ifges(7,5) * t369 + Ifges(7,6) * t368 + Ifges(7,3) * t379 + t371 * t398 - t372 * t397;
t410 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t442 + Ifges(5,4) * t446) * qJD(2);
t409 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t442 + Ifges(5,2) * t446) * qJD(2);
t389 = Ifges(6,5) * t412 - Ifges(6,6) * t411 + Ifges(6,3) * t431;
t321 = -mrSges(6,1) * t359 + mrSges(6,3) * t353 + Ifges(6,4) * t381 + Ifges(6,2) * t380 + Ifges(6,6) * t430 - pkin(5) * t334 - t389 * t412 + t391 * t431 - t450;
t320 = mrSges(6,2) * t359 - mrSges(6,3) * t352 + Ifges(6,1) * t381 + Ifges(6,4) * t380 + Ifges(6,5) * t430 - pkin(10) * t334 - t337 * t440 + t338 * t444 - t389 * t411 - t390 * t431;
t1 = [m(2) * t433 + (m(3) * t393 - mrSges(3,1) * t448 - qJDD(2) * mrSges(3,2) + t319 * t437 - t329 * t434) * t467 + (m(3) * t392 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t448 + t465) * t466 + t439 * (m(3) * t456 + t461); mrSges(3,1) * t392 - mrSges(3,2) * t393 + mrSges(4,1) * t365 - mrSges(4,2) * t366 + t442 * (mrSges(5,2) * t363 - mrSges(5,3) * t360 + Ifges(5,1) * t418 + Ifges(5,4) * t419 + Ifges(5,5) * qJDD(4) - pkin(9) * t327 - qJD(4) * t409 + t320 * t445 - t321 * t441) + t446 * (-mrSges(5,1) * t363 + mrSges(5,3) * t361 + Ifges(5,4) * t418 + Ifges(5,2) * t419 + Ifges(5,6) * qJDD(4) - pkin(4) * t452 + pkin(9) * t458 + qJD(4) * t410 + t441 * t320 + t445 * t321) + pkin(3) * t449 + pkin(8) * t459 + pkin(2) * t465 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t461; Ifges(5,3) * qJDD(4) + (t409 * t442 - t410 * t446) * qJD(2) + t451 + Ifges(5,5) * t418 + Ifges(5,6) * t419 + mrSges(5,1) * t360 - mrSges(5,2) * t361 + pkin(4) * t327; t451; t450;];
tauJ  = t1;
