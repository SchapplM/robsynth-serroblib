% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRP6
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
% Datum: 2019-05-05 00:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:59:33
% EndTime: 2019-05-04 23:59:35
% DurationCPUTime: 1.16s
% Computational Cost: add. (5053->209), mult. (9216->251), div. (0->0), fcn. (5726->10), ass. (0->92)
t465 = Ifges(6,1) + Ifges(7,1);
t452 = Ifges(6,4) - Ifges(7,5);
t460 = -Ifges(6,5) - Ifges(7,4);
t464 = Ifges(6,2) + Ifges(7,3);
t450 = Ifges(6,6) - Ifges(7,6);
t420 = sin(qJ(5));
t423 = cos(qJ(4));
t442 = qJD(2) * t423;
t454 = cos(qJ(5));
t395 = -qJD(4) * t454 + t420 * t442;
t421 = sin(qJ(4));
t441 = qJD(2) * qJD(4);
t438 = t421 * t441;
t400 = qJDD(2) * t423 - t438;
t369 = -t395 * qJD(5) + t420 * qJDD(4) + t400 * t454;
t396 = t420 * qJD(4) + t442 * t454;
t373 = mrSges(7,1) * t395 - mrSges(7,3) * t396;
t416 = sin(pkin(10));
t418 = cos(pkin(10));
t402 = -g(1) * t418 - g(2) * t416;
t422 = sin(qJ(2));
t424 = cos(qJ(2));
t401 = g(1) * t416 - g(2) * t418;
t413 = -g(3) + qJDD(1);
t417 = sin(pkin(6));
t419 = cos(pkin(6));
t462 = t401 * t419 + t413 * t417;
t356 = -t422 * t402 + t424 * t462;
t426 = qJD(2) ^ 2;
t430 = -qJ(3) * t426 + qJDD(3) - t356;
t455 = -pkin(2) - pkin(8);
t353 = qJDD(2) * t455 + t430;
t381 = -t401 * t417 + t413 * t419;
t349 = t421 * t353 + t423 * t381;
t398 = (pkin(4) * t421 - pkin(9) * t423) * qJD(2);
t425 = qJD(4) ^ 2;
t443 = qJD(2) * t421;
t345 = -pkin(4) * t425 + qJDD(4) * pkin(9) - t398 * t443 + t349;
t357 = t424 * t402 + t462 * t422;
t456 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t357;
t352 = t426 * t455 - t456;
t437 = t423 * t441;
t399 = -qJDD(2) * t421 - t437;
t347 = (-t400 + t438) * pkin(9) + (-t399 + t437) * pkin(4) + t352;
t341 = -t420 * t345 + t347 * t454;
t372 = pkin(5) * t395 - qJ(6) * t396;
t392 = qJDD(5) - t399;
t407 = qJD(5) + t443;
t406 = t407 ^ 2;
t339 = -t392 * pkin(5) - t406 * qJ(6) + t396 * t372 + qJDD(6) - t341;
t380 = -mrSges(7,2) * t395 + mrSges(7,3) * t407;
t434 = -m(7) * t339 + t392 * mrSges(7,1) + t407 * t380;
t335 = mrSges(7,2) * t369 + t373 * t396 - t434;
t342 = t454 * t345 + t420 * t347;
t338 = -pkin(5) * t406 + qJ(6) * t392 + 0.2e1 * qJD(6) * t407 - t372 * t395 + t342;
t368 = qJD(5) * t396 - qJDD(4) * t454 + t400 * t420;
t379 = -mrSges(7,1) * t407 + mrSges(7,2) * t396;
t439 = m(7) * t338 + t392 * mrSges(7,3) + t407 * t379;
t446 = t464 * t395 - t452 * t396 - t450 * t407;
t457 = t452 * t395 - t465 * t396 + t460 * t407;
t459 = -Ifges(6,3) - Ifges(7,2);
t463 = -t369 * t460 - t457 * t395 - t450 * t368 - t459 * t392 + mrSges(6,1) * t341 - mrSges(7,1) * t339 - mrSges(6,2) * t342 + mrSges(7,3) * t338 - pkin(5) * t335 + qJ(6) * (-mrSges(7,2) * t368 - t373 * t395 + t439) - t446 * t396;
t453 = -mrSges(6,3) - mrSges(7,2);
t378 = mrSges(6,1) * t407 - mrSges(6,3) * t396;
t444 = -mrSges(6,1) * t395 - mrSges(6,2) * t396 - t373;
t331 = m(6) * t342 - mrSges(6,2) * t392 + t368 * t453 - t378 * t407 + t395 * t444 + t439;
t377 = -mrSges(6,2) * t407 - mrSges(6,3) * t395;
t333 = m(6) * t341 + mrSges(6,1) * t392 + t369 * t453 + t377 * t407 + t396 * t444 + t434;
t327 = t420 * t331 + t454 * t333;
t447 = t450 * t395 + t460 * t396 + t459 * t407;
t436 = t454 * t331 - t333 * t420;
t348 = t353 * t423 - t421 * t381;
t397 = (mrSges(5,1) * t421 + mrSges(5,2) * t423) * qJD(2);
t404 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t442;
t324 = m(5) * t349 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t399 - qJD(4) * t404 - t397 * t443 + t436;
t403 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t443;
t344 = -qJDD(4) * pkin(4) - pkin(9) * t425 + t398 * t442 - t348;
t340 = -0.2e1 * qJD(6) * t396 + (t395 * t407 - t369) * qJ(6) + (t396 * t407 + t368) * pkin(5) + t344;
t336 = m(7) * t340 + mrSges(7,1) * t368 - t369 * mrSges(7,3) - t396 * t379 + t380 * t395;
t427 = -m(6) * t344 - t368 * mrSges(6,1) - mrSges(6,2) * t369 - t395 * t377 - t378 * t396 - t336;
t328 = m(5) * t348 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t400 + qJD(4) * t403 - t397 * t442 + t427;
t433 = t324 * t421 + t328 * t423;
t355 = -qJDD(2) * pkin(2) + t430;
t431 = -m(4) * t355 + t426 * mrSges(4,3) - t433;
t354 = pkin(2) * t426 + t456;
t429 = -m(4) * t354 + m(5) * t352 - mrSges(5,1) * t399 + t426 * mrSges(4,2) + t400 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t403 * t443 + t404 * t442 + t327;
t386 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t423 - Ifges(5,4) * t421) * qJD(2);
t385 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t423 - Ifges(5,2) * t421) * qJD(2);
t326 = mrSges(6,2) * t344 + mrSges(7,2) * t339 - mrSges(6,3) * t341 - mrSges(7,3) * t340 - qJ(6) * t336 - t452 * t368 + t465 * t369 - t460 * t392 + t447 * t395 + t446 * t407;
t325 = -mrSges(6,1) * t344 - mrSges(7,1) * t340 + mrSges(7,2) * t338 + mrSges(6,3) * t342 - pkin(5) * t336 - t464 * t368 + t452 * t369 + t450 * t392 + t447 * t396 - t457 * t407;
t323 = qJDD(2) * mrSges(4,2) - t431;
t1 = [m(2) * t413 + t419 * (t324 * t423 - t328 * t421 + (m(3) + m(4)) * t381) + (t422 * (m(3) * t357 - mrSges(3,1) * t426 - qJDD(2) * mrSges(3,2) + t429) + t424 * (m(3) * t356 - mrSges(3,2) * t426 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t431)) * t417; mrSges(3,1) * t356 - mrSges(3,2) * t357 + mrSges(4,2) * t355 - mrSges(4,3) * t354 + t423 * (mrSges(5,2) * t352 - mrSges(5,3) * t348 + Ifges(5,1) * t400 + Ifges(5,4) * t399 + Ifges(5,5) * qJDD(4) - pkin(9) * t327 - qJD(4) * t385 - t420 * t325 + t454 * t326) - t421 * (-mrSges(5,1) * t352 + mrSges(5,3) * t349 + Ifges(5,4) * t400 + Ifges(5,2) * t399 + Ifges(5,6) * qJDD(4) - pkin(4) * t327 + qJD(4) * t386 - t463) - pkin(8) * t433 - pkin(2) * t323 + qJ(3) * t429 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t323; Ifges(5,5) * t400 + Ifges(5,6) * t399 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t348 - mrSges(5,2) * t349 + t420 * t326 + t454 * t325 + pkin(4) * t427 + pkin(9) * t436 + (t385 * t423 + t386 * t421) * qJD(2); t463; t335;];
tauJ  = t1;
