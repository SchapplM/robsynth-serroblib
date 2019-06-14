% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-04 23:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:54:30
% EndTime: 2019-05-04 23:54:33
% DurationCPUTime: 1.15s
% Computational Cost: add. (5117->211), mult. (9416->251), div. (0->0), fcn. (5881->10), ass. (0->93)
t466 = Ifges(6,1) + Ifges(7,1);
t455 = Ifges(6,4) + Ifges(7,4);
t454 = Ifges(6,5) + Ifges(7,5);
t465 = Ifges(6,2) + Ifges(7,2);
t453 = Ifges(6,6) + Ifges(7,6);
t421 = sin(qJ(5));
t424 = cos(qJ(5));
t425 = cos(qJ(4));
t445 = qJD(2) * t425;
t399 = qJD(4) * t424 - t421 * t445;
t422 = sin(qJ(4));
t444 = qJD(2) * qJD(4);
t440 = t422 * t444;
t404 = qJDD(2) * t425 - t440;
t374 = qJD(5) * t399 + qJDD(4) * t421 + t404 * t424;
t400 = qJD(4) * t421 + t424 * t445;
t376 = -mrSges(7,1) * t399 + mrSges(7,2) * t400;
t417 = sin(pkin(10));
t419 = cos(pkin(10));
t406 = -g(1) * t419 - g(2) * t417;
t423 = sin(qJ(2));
t426 = cos(qJ(2));
t405 = g(1) * t417 - g(2) * t419;
t414 = -g(3) + qJDD(1);
t418 = sin(pkin(6));
t420 = cos(pkin(6));
t463 = t405 * t420 + t414 * t418;
t360 = -t423 * t406 + t426 * t463;
t428 = qJD(2) ^ 2;
t431 = -t428 * qJ(3) + qJDD(3) - t360;
t457 = -pkin(2) - pkin(8);
t356 = qJDD(2) * t457 + t431;
t385 = -t405 * t418 + t414 * t420;
t352 = t422 * t356 + t425 * t385;
t402 = (pkin(4) * t422 - pkin(9) * t425) * qJD(2);
t427 = qJD(4) ^ 2;
t446 = qJD(2) * t422;
t347 = -pkin(4) * t427 + qJDD(4) * pkin(9) - t402 * t446 + t352;
t361 = t426 * t406 + t463 * t423;
t458 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t361;
t355 = t428 * t457 - t458;
t439 = t425 * t444;
t403 = -qJDD(2) * t422 - t439;
t350 = (-t404 + t440) * pkin(9) + (-t403 + t439) * pkin(4) + t355;
t342 = -t421 * t347 + t424 * t350;
t396 = qJDD(5) - t403;
t410 = qJD(5) + t446;
t339 = -0.2e1 * qJD(6) * t400 + (t399 * t410 - t374) * qJ(6) + (t399 * t400 + t396) * pkin(5) + t342;
t380 = -mrSges(7,2) * t410 + mrSges(7,3) * t399;
t442 = m(7) * t339 + t396 * mrSges(7,1) + t410 * t380;
t336 = -t374 * mrSges(7,3) - t400 * t376 + t442;
t343 = t424 * t347 + t421 * t350;
t373 = -qJD(5) * t400 + qJDD(4) * t424 - t404 * t421;
t382 = pkin(5) * t410 - qJ(6) * t400;
t395 = t399 ^ 2;
t341 = -pkin(5) * t395 + qJ(6) * t373 + 0.2e1 * qJD(6) * t399 - t382 * t410 + t343;
t449 = -t465 * t399 - t455 * t400 - t453 * t410;
t459 = t455 * t399 + t466 * t400 + t454 * t410;
t461 = Ifges(6,3) + Ifges(7,3);
t464 = mrSges(6,1) * t342 + mrSges(7,1) * t339 - mrSges(6,2) * t343 - mrSges(7,2) * t341 + pkin(5) * t336 + t453 * t373 + t454 * t374 + t461 * t396 - t459 * t399 - t449 * t400;
t456 = -mrSges(6,2) - mrSges(7,2);
t377 = -mrSges(6,1) * t399 + mrSges(6,2) * t400;
t381 = -mrSges(6,2) * t410 + mrSges(6,3) * t399;
t331 = m(6) * t342 + t396 * mrSges(6,1) + t410 * t381 + (-t376 - t377) * t400 + (-mrSges(6,3) - mrSges(7,3)) * t374 + t442;
t441 = m(7) * t341 + t373 * mrSges(7,3) + t399 * t376;
t383 = mrSges(7,1) * t410 - mrSges(7,3) * t400;
t447 = -mrSges(6,1) * t410 + mrSges(6,3) * t400 - t383;
t334 = m(6) * t343 + t373 * mrSges(6,3) + t399 * t377 + t396 * t456 + t410 * t447 + t441;
t329 = t424 * t331 + t421 * t334;
t450 = -t453 * t399 - t454 * t400 - t461 * t410;
t438 = -t331 * t421 + t424 * t334;
t351 = t356 * t425 - t422 * t385;
t346 = -qJDD(4) * pkin(4) - pkin(9) * t427 + t402 * t445 - t351;
t344 = -pkin(5) * t373 - qJ(6) * t395 + t382 * t400 + qJDD(6) + t346;
t436 = -m(7) * t344 + t373 * mrSges(7,1) + t399 * t380;
t401 = (mrSges(5,1) * t422 + mrSges(5,2) * t425) * qJD(2);
t408 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t445;
t327 = m(5) * t352 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t403 - qJD(4) * t408 - t401 * t446 + t438;
t407 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t446;
t429 = -m(6) * t346 + t373 * mrSges(6,1) + t374 * t456 + t399 * t381 + t400 * t447 + t436;
t335 = m(5) * t351 + qJDD(4) * mrSges(5,1) - t404 * mrSges(5,3) + qJD(4) * t407 - t401 * t445 + t429;
t435 = t422 * t327 + t425 * t335;
t358 = -qJDD(2) * pkin(2) + t431;
t433 = -m(4) * t358 + t428 * mrSges(4,3) - t435;
t357 = t428 * pkin(2) + t458;
t430 = -m(4) * t357 + m(5) * t355 - mrSges(5,1) * t403 + t428 * mrSges(4,2) + t404 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t407 * t446 + t408 * t445 + t329;
t390 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t425 - Ifges(5,4) * t422) * qJD(2);
t389 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t425 - Ifges(5,2) * t422) * qJD(2);
t337 = t374 * mrSges(7,2) + t400 * t383 - t436;
t328 = mrSges(6,2) * t346 + mrSges(7,2) * t344 - mrSges(6,3) * t342 - mrSges(7,3) * t339 - qJ(6) * t336 + t455 * t373 + t466 * t374 + t454 * t396 - t450 * t399 + t449 * t410;
t326 = -mrSges(6,1) * t346 + mrSges(6,3) * t343 - mrSges(7,1) * t344 + mrSges(7,3) * t341 - pkin(5) * t337 + qJ(6) * t441 + (-qJ(6) * t383 + t459) * t410 + t450 * t400 + (-mrSges(7,2) * qJ(6) + t453) * t396 + t455 * t374 + t465 * t373;
t325 = qJDD(2) * mrSges(4,2) - t433;
t1 = [m(2) * t414 + t420 * (t425 * t327 - t422 * t335 + (m(3) + m(4)) * t385) + (t423 * (m(3) * t361 - mrSges(3,1) * t428 - qJDD(2) * mrSges(3,2) + t430) + t426 * (m(3) * t360 - t428 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t433)) * t418; mrSges(3,1) * t360 - mrSges(3,2) * t361 + mrSges(4,2) * t358 - mrSges(4,3) * t357 + t425 * (mrSges(5,2) * t355 - mrSges(5,3) * t351 + Ifges(5,1) * t404 + Ifges(5,4) * t403 + Ifges(5,5) * qJDD(4) - pkin(9) * t329 - qJD(4) * t389 - t326 * t421 + t328 * t424) - t422 * (-mrSges(5,1) * t355 + mrSges(5,3) * t352 + Ifges(5,4) * t404 + Ifges(5,2) * t403 + Ifges(5,6) * qJDD(4) - pkin(4) * t329 + qJD(4) * t390 - t464) - pkin(8) * t435 - pkin(2) * t325 + qJ(3) * t430 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t325; Ifges(5,5) * t404 + Ifges(5,6) * t403 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t351 - mrSges(5,2) * t352 + t421 * t328 + t424 * t326 + pkin(4) * t429 + pkin(9) * t438 + (t389 * t425 + t390 * t422) * qJD(2); t464; t337;];
tauJ  = t1;
