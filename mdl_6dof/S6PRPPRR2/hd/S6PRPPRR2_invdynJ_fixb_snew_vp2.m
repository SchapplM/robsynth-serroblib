% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-04 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:50:26
% EndTime: 2019-05-04 21:50:27
% DurationCPUTime: 0.71s
% Computational Cost: add. (4296->179), mult. (7560->226), div. (0->0), fcn. (4988->12), ass. (0->87)
t372 = sin(pkin(10));
t375 = cos(pkin(10));
t358 = -t375 * g(1) - t372 * g(2);
t368 = -g(3) + qJDD(1);
t379 = sin(qJ(2));
t382 = cos(qJ(2));
t373 = sin(pkin(6));
t403 = t373 * t382;
t357 = t372 * g(1) - t375 * g(2);
t376 = cos(pkin(6));
t405 = t357 * t376;
t325 = -t379 * t358 + t368 * t403 + t382 * t405;
t323 = qJDD(2) * pkin(2) + t325;
t404 = t373 * t379;
t326 = t382 * t358 + t368 * t404 + t379 * t405;
t384 = qJD(2) ^ 2;
t324 = -t384 * pkin(2) + t326;
t371 = sin(pkin(11));
t374 = cos(pkin(11));
t319 = t371 * t323 + t374 * t324;
t407 = -qJDD(2) * qJ(4) - (2 * qJD(4) * qJD(2)) - t319;
t406 = -pkin(3) - pkin(8);
t392 = -t373 * t357 + t376 * t368;
t340 = qJDD(3) + t392;
t378 = sin(qJ(5));
t402 = t378 * t340;
t318 = t374 * t323 - t371 * t324;
t389 = -t384 * qJ(4) + qJDD(4) - t318;
t317 = -qJDD(2) * pkin(3) + t389;
t315 = t406 * qJDD(2) + t389;
t381 = cos(qJ(5));
t311 = t378 * t315 + t381 * t340;
t353 = (t378 * mrSges(6,1) + t381 * mrSges(6,2)) * qJD(2);
t398 = qJD(2) * qJD(5);
t394 = t381 * t398;
t355 = -t378 * qJDD(2) - t394;
t400 = qJD(2) * t381;
t360 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t400;
t354 = (t378 * pkin(5) - t381 * pkin(9)) * qJD(2);
t383 = qJD(5) ^ 2;
t399 = t378 * qJD(2);
t308 = -t383 * pkin(5) + qJDD(5) * pkin(9) - t354 * t399 + t311;
t314 = t406 * t384 - t407;
t395 = t378 * t398;
t356 = t381 * qJDD(2) - t395;
t309 = (-t356 + t395) * pkin(9) + (-t355 + t394) * pkin(5) + t314;
t377 = sin(qJ(6));
t380 = cos(qJ(6));
t305 = -t377 * t308 + t380 * t309;
t351 = t380 * qJD(5) - t377 * t400;
t333 = t351 * qJD(6) + t377 * qJDD(5) + t380 * t356;
t352 = t377 * qJD(5) + t380 * t400;
t334 = -t351 * mrSges(7,1) + t352 * mrSges(7,2);
t363 = qJD(6) + t399;
t338 = -t363 * mrSges(7,2) + t351 * mrSges(7,3);
t349 = qJDD(6) - t355;
t303 = m(7) * t305 + t349 * mrSges(7,1) - t333 * mrSges(7,3) - t352 * t334 + t363 * t338;
t306 = t380 * t308 + t377 * t309;
t332 = -t352 * qJD(6) + t380 * qJDD(5) - t377 * t356;
t339 = t363 * mrSges(7,1) - t352 * mrSges(7,3);
t304 = m(7) * t306 - t349 * mrSges(7,2) + t332 * mrSges(7,3) + t351 * t334 - t363 * t339;
t393 = -t377 * t303 + t380 * t304;
t295 = m(6) * t311 - qJDD(5) * mrSges(6,2) + t355 * mrSges(6,3) - qJD(5) * t360 - t353 * t399 + t393;
t310 = t381 * t315 - t402;
t359 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t399;
t307 = -qJDD(5) * pkin(5) - t383 * pkin(9) + t402 + (qJD(2) * t354 - t315) * t381;
t387 = -m(7) * t307 + t332 * mrSges(7,1) - t333 * mrSges(7,2) + t351 * t338 - t352 * t339;
t299 = m(6) * t310 + qJDD(5) * mrSges(6,1) - t356 * mrSges(6,3) + qJD(5) * t359 - t353 * t400 + t387;
t391 = t378 * t295 + t381 * t299;
t388 = -m(5) * t317 + t384 * mrSges(5,3) - t391;
t290 = m(4) * t318 - t384 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + t388;
t296 = t380 * t303 + t377 * t304;
t316 = t384 * pkin(3) + t407;
t386 = -m(5) * t316 + m(6) * t314 - t355 * mrSges(6,1) + t384 * mrSges(5,2) + t356 * mrSges(6,2) + qJDD(2) * mrSges(5,3) + t359 * t399 + t360 * t400 + t296;
t293 = m(4) * t319 - t384 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t386;
t401 = t374 * t290 + t371 * t293;
t390 = t381 * t295 - t378 * t299 + (m(4) + m(5)) * t340;
t328 = Ifges(7,4) * t352 + Ifges(7,2) * t351 + Ifges(7,6) * t363;
t329 = Ifges(7,1) * t352 + Ifges(7,4) * t351 + Ifges(7,5) * t363;
t385 = mrSges(7,1) * t305 - mrSges(7,2) * t306 + Ifges(7,5) * t333 + Ifges(7,6) * t332 + Ifges(7,3) * t349 + t352 * t328 - t351 * t329;
t345 = (Ifges(6,5) * qJD(5)) + (t381 * Ifges(6,1) - t378 * Ifges(6,4)) * qJD(2);
t344 = (Ifges(6,6) * qJD(5)) + (t381 * Ifges(6,4) - t378 * Ifges(6,2)) * qJD(2);
t327 = Ifges(7,5) * t352 + Ifges(7,6) * t351 + Ifges(7,3) * t363;
t298 = mrSges(7,2) * t307 - mrSges(7,3) * t305 + Ifges(7,1) * t333 + Ifges(7,4) * t332 + Ifges(7,5) * t349 + t351 * t327 - t363 * t328;
t297 = -mrSges(7,1) * t307 + mrSges(7,3) * t306 + Ifges(7,4) * t333 + Ifges(7,2) * t332 + Ifges(7,6) * t349 - t352 * t327 + t363 * t329;
t291 = qJDD(2) * mrSges(5,2) - t388;
t1 = [m(2) * t368 + (m(3) * t326 - t384 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t371 * t290 + t374 * t293) * t404 + (m(3) * t325 + qJDD(2) * mrSges(3,1) - t384 * mrSges(3,2) + t401) * t403 + t376 * (m(3) * t392 + t390); pkin(2) * t401 + mrSges(3,1) * t325 - mrSges(3,2) * t326 - pkin(3) * t291 + qJ(4) * t386 + t381 * (mrSges(6,2) * t314 - mrSges(6,3) * t310 + Ifges(6,1) * t356 + Ifges(6,4) * t355 + Ifges(6,5) * qJDD(5) - pkin(9) * t296 - qJD(5) * t344 - t377 * t297 + t380 * t298) - t378 * (-mrSges(6,1) * t314 + mrSges(6,3) * t311 + Ifges(6,4) * t356 + Ifges(6,2) * t355 + Ifges(6,6) * qJDD(5) - pkin(5) * t296 + qJD(5) * t345 - t385) - pkin(8) * t391 + mrSges(4,1) * t318 - mrSges(4,2) * t319 + mrSges(5,2) * t317 - mrSges(5,3) * t316 + (Ifges(3,3) + Ifges(4,3) + Ifges(5,1)) * qJDD(2); t390; t291; Ifges(6,5) * t356 + Ifges(6,6) * t355 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t310 - mrSges(6,2) * t311 + t377 * t298 + t380 * t297 + pkin(5) * t387 + pkin(9) * t393 + (t381 * t344 + t378 * t345) * qJD(2); t385;];
tauJ  = t1;
