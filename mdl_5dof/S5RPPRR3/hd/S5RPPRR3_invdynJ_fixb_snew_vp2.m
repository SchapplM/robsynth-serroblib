% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:21
% EndTime: 2022-01-23 09:14:22
% DurationCPUTime: 1.07s
% Computational Cost: add. (7485->192), mult. (16642->248), div. (0->0), fcn. (11244->10), ass. (0->90)
t371 = qJD(1) ^ 2;
t363 = cos(pkin(9));
t393 = t363 * pkin(3);
t361 = sin(pkin(9));
t392 = mrSges(4,2) * t361;
t358 = t363 ^ 2;
t391 = t358 * t371;
t367 = sin(qJ(1));
t370 = cos(qJ(1));
t384 = t367 * g(1) - g(2) * t370;
t345 = qJDD(1) * pkin(1) + t384;
t380 = -g(1) * t370 - g(2) * t367;
t346 = -pkin(1) * t371 + t380;
t362 = sin(pkin(8));
t364 = cos(pkin(8));
t333 = t362 * t345 + t364 * t346;
t326 = -pkin(2) * t371 + qJDD(1) * qJ(3) + t333;
t360 = -g(3) + qJDD(2);
t386 = qJD(1) * qJD(3);
t389 = t363 * t360 - 0.2e1 * t361 * t386;
t312 = (-pkin(6) * qJDD(1) + t371 * t393 - t326) * t361 + t389;
t316 = t361 * t360 + (t326 + 0.2e1 * t386) * t363;
t385 = t363 * qJDD(1);
t313 = -pkin(3) * t391 + pkin(6) * t385 + t316;
t366 = sin(qJ(4));
t369 = cos(qJ(4));
t294 = t369 * t312 - t313 * t366;
t377 = t361 * t369 + t363 * t366;
t376 = -t361 * t366 + t363 * t369;
t338 = t376 * qJD(1);
t387 = qJD(4) * t338;
t331 = t377 * qJDD(1) + t387;
t339 = t377 * qJD(1);
t290 = (-t331 + t387) * pkin(7) + (t338 * t339 + qJDD(4)) * pkin(4) + t294;
t295 = t366 * t312 + t369 * t313;
t330 = -qJD(4) * t339 + t376 * qJDD(1);
t336 = qJD(4) * pkin(4) - pkin(7) * t339;
t337 = t338 ^ 2;
t291 = -pkin(4) * t337 + pkin(7) * t330 - qJD(4) * t336 + t295;
t365 = sin(qJ(5));
t368 = cos(qJ(5));
t288 = t290 * t368 - t291 * t365;
t324 = t338 * t368 - t339 * t365;
t302 = qJD(5) * t324 + t330 * t365 + t331 * t368;
t325 = t338 * t365 + t339 * t368;
t308 = -mrSges(6,1) * t324 + mrSges(6,2) * t325;
t359 = qJD(4) + qJD(5);
t317 = -mrSges(6,2) * t359 + mrSges(6,3) * t324;
t356 = qJDD(4) + qJDD(5);
t285 = m(6) * t288 + mrSges(6,1) * t356 - mrSges(6,3) * t302 - t308 * t325 + t317 * t359;
t289 = t290 * t365 + t291 * t368;
t301 = -qJD(5) * t325 + t330 * t368 - t331 * t365;
t318 = mrSges(6,1) * t359 - mrSges(6,3) * t325;
t286 = m(6) * t289 - mrSges(6,2) * t356 + mrSges(6,3) * t301 + t308 * t324 - t318 * t359;
t278 = t368 * t285 + t365 * t286;
t328 = -mrSges(5,1) * t338 + mrSges(5,2) * t339;
t334 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t338;
t276 = m(5) * t294 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t331 + qJD(4) * t334 - t328 * t339 + t278;
t335 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t339;
t381 = -t285 * t365 + t368 * t286;
t277 = m(5) * t295 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t330 - qJD(4) * t335 + t328 * t338 + t381;
t390 = t369 * t276 + t366 * t277;
t315 = -t326 * t361 + t389;
t375 = mrSges(4,3) * qJDD(1) + t371 * (-t363 * mrSges(4,1) + t392);
t271 = m(4) * t315 - t361 * t375 + t390;
t382 = -t276 * t366 + t369 * t277;
t272 = m(4) * t316 + t363 * t375 + t382;
t383 = -t271 * t361 + t363 * t272;
t332 = t345 * t364 - t362 * t346;
t379 = qJDD(3) - t332;
t357 = t361 ^ 2;
t314 = (-pkin(2) - t393) * qJDD(1) + (-qJ(3) + (-t357 - t358) * pkin(6)) * t371 + t379;
t293 = -pkin(4) * t330 - pkin(7) * t337 + t336 * t339 + t314;
t378 = m(6) * t293 - t301 * mrSges(6,1) + t302 * mrSges(6,2) - t324 * t317 + t325 * t318;
t304 = Ifges(6,4) * t325 + Ifges(6,2) * t324 + Ifges(6,6) * t359;
t305 = Ifges(6,1) * t325 + Ifges(6,4) * t324 + Ifges(6,5) * t359;
t374 = mrSges(6,1) * t288 - mrSges(6,2) * t289 + Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t356 + t325 * t304 - t305 * t324;
t373 = m(5) * t314 - t330 * mrSges(5,1) + mrSges(5,2) * t331 - t338 * t334 + t335 * t339 + t378;
t320 = -qJDD(1) * pkin(2) - qJ(3) * t371 + t379;
t372 = -m(4) * t320 + mrSges(4,1) * t385 - t373 + (t357 * t371 + t391) * mrSges(4,3);
t323 = Ifges(5,1) * t339 + Ifges(5,4) * t338 + Ifges(5,5) * qJD(4);
t322 = Ifges(5,4) * t339 + Ifges(5,2) * t338 + Ifges(5,6) * qJD(4);
t321 = Ifges(5,5) * t339 + Ifges(5,6) * t338 + Ifges(5,3) * qJD(4);
t303 = Ifges(6,5) * t325 + Ifges(6,6) * t324 + Ifges(6,3) * t359;
t281 = qJDD(1) * t392 - t372;
t280 = mrSges(6,2) * t293 - mrSges(6,3) * t288 + Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t356 + t303 * t324 - t304 * t359;
t279 = -mrSges(6,1) * t293 + mrSges(6,3) * t289 + Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t356 - t303 * t325 + t305 * t359;
t269 = mrSges(5,2) * t314 - mrSges(5,3) * t294 + Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * qJDD(4) - pkin(7) * t278 - qJD(4) * t322 - t279 * t365 + t280 * t368 + t321 * t338;
t268 = -mrSges(5,1) * t314 + mrSges(5,3) * t295 + Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * qJDD(4) - pkin(4) * t378 + pkin(7) * t381 + qJD(4) * t323 + t368 * t279 + t365 * t280 - t339 * t321;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t384 - mrSges(2,2) * t380 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t332 - mrSges(3,2) * t333 + t361 * (mrSges(4,2) * t320 - mrSges(4,3) * t315 + t369 * t269 - t366 * t268 - pkin(6) * t390 + (Ifges(4,1) * t361 + Ifges(4,4) * t363) * qJDD(1)) + t363 * (-mrSges(4,1) * t320 + mrSges(4,3) * t316 + t366 * t269 + t369 * t268 - pkin(3) * t373 + pkin(6) * t382 + (Ifges(4,4) * t361 + Ifges(4,2) * t363) * qJDD(1)) - pkin(2) * t281 + qJ(3) * t383 + pkin(1) * (t362 * (m(3) * t333 - mrSges(3,1) * t371 - qJDD(1) * mrSges(3,2) + t383) + t364 * (t372 + (mrSges(3,1) - t392) * qJDD(1) + m(3) * t332 - mrSges(3,2) * t371)); m(3) * t360 + t271 * t363 + t272 * t361; t281; mrSges(5,1) * t294 - mrSges(5,2) * t295 + Ifges(5,5) * t331 + Ifges(5,6) * t330 + Ifges(5,3) * qJDD(4) + pkin(4) * t278 + t322 * t339 - t323 * t338 + t374; t374;];
tauJ = t1;
