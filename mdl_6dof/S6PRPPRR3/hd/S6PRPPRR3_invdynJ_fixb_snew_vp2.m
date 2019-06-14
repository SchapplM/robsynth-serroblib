% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-04 22:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:00:12
% EndTime: 2019-05-04 22:00:13
% DurationCPUTime: 0.83s
% Computational Cost: add. (5623->181), mult. (9730->227), div. (0->0), fcn. (5808->12), ass. (0->85)
t340 = sin(pkin(10));
t343 = cos(pkin(10));
t326 = g(1) * t340 - g(2) * t343;
t336 = -g(3) + qJDD(1);
t341 = sin(pkin(6));
t344 = cos(pkin(6));
t374 = t326 * t344 + t336 * t341;
t327 = -g(1) * t343 - g(2) * t340;
t347 = sin(qJ(2));
t350 = cos(qJ(2));
t296 = -t347 * t327 + t374 * t350;
t373 = -pkin(2) - pkin(3);
t359 = -t326 * t341 + t336 * t344;
t309 = qJDD(4) - t359;
t349 = cos(qJ(5));
t370 = t349 * t309;
t352 = qJD(2) ^ 2;
t297 = t350 * t327 + t374 * t347;
t362 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t297;
t291 = t373 * t352 + t362;
t355 = -t352 * qJ(3) + qJDD(3) - t296;
t293 = t373 * qJDD(2) + t355;
t339 = sin(pkin(11));
t342 = cos(pkin(11));
t287 = t342 * t291 + t339 * t293;
t285 = -pkin(4) * t352 - qJDD(2) * pkin(8) + t287;
t346 = sin(qJ(5));
t282 = t349 * t285 + t346 * t309;
t369 = qJD(2) * t346;
t368 = qJD(2) * t349;
t367 = qJD(2) * qJD(5);
t366 = t346 * t367;
t365 = t349 * t367;
t322 = (mrSges(6,1) * t349 - mrSges(6,2) * t346) * qJD(2);
t325 = -qJDD(2) * t349 + t366;
t328 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t369;
t323 = (pkin(5) * t349 + pkin(9) * t346) * qJD(2);
t351 = qJD(5) ^ 2;
t279 = -pkin(5) * t351 + qJDD(5) * pkin(9) - t323 * t368 + t282;
t286 = -t339 * t291 + t342 * t293;
t284 = qJDD(2) * pkin(4) - t352 * pkin(8) - t286;
t324 = -qJDD(2) * t346 - t365;
t280 = (-t324 + t365) * pkin(9) + (-t325 - t366) * pkin(5) + t284;
t345 = sin(qJ(6));
t348 = cos(qJ(6));
t276 = -t279 * t345 + t280 * t348;
t320 = qJD(5) * t348 + t345 * t369;
t304 = qJD(6) * t320 + qJDD(5) * t345 + t324 * t348;
t321 = qJD(5) * t345 - t348 * t369;
t305 = -mrSges(7,1) * t320 + mrSges(7,2) * t321;
t331 = qJD(6) + t368;
t307 = -mrSges(7,2) * t331 + mrSges(7,3) * t320;
t317 = qJDD(6) - t325;
t274 = m(7) * t276 + mrSges(7,1) * t317 - mrSges(7,3) * t304 - t305 * t321 + t307 * t331;
t277 = t279 * t348 + t280 * t345;
t303 = -qJD(6) * t321 + qJDD(5) * t348 - t324 * t345;
t308 = mrSges(7,1) * t331 - mrSges(7,3) * t321;
t275 = m(7) * t277 - mrSges(7,2) * t317 + mrSges(7,3) * t303 + t305 * t320 - t308 * t331;
t363 = -t274 * t345 + t348 * t275;
t268 = m(6) * t282 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t325 - qJD(5) * t328 - t322 * t368 + t363;
t281 = -t285 * t346 + t370;
t329 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t368;
t278 = -qJDD(5) * pkin(5) - t351 * pkin(9) - t370 + (-qJD(2) * t323 + t285) * t346;
t356 = -m(7) * t278 + t303 * mrSges(7,1) - mrSges(7,2) * t304 + t320 * t307 - t308 * t321;
t272 = m(6) * t281 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t324 + qJD(5) * t329 + t322 * t369 + t356;
t364 = t349 * t268 - t272 * t346;
t265 = m(5) * t287 - mrSges(5,1) * t352 + qJDD(2) * mrSges(5,2) + t364;
t269 = t274 * t348 + t275 * t345;
t354 = -m(6) * t284 + t325 * mrSges(6,1) - mrSges(6,2) * t324 + t328 * t369 - t329 * t368 - t269;
t266 = m(5) * t286 - qJDD(2) * mrSges(5,1) - mrSges(5,2) * t352 + t354;
t361 = t265 * t339 + t266 * t342;
t294 = -pkin(2) * t352 + t362;
t358 = m(4) * t294 + qJDD(2) * mrSges(4,3) + t342 * t265 - t339 * t266;
t357 = m(5) * t309 + t346 * t268 + t349 * t272;
t295 = -qJDD(2) * pkin(2) + t355;
t263 = m(4) * t295 - qJDD(2) * mrSges(4,1) - t352 * mrSges(4,3) + t361;
t299 = Ifges(7,4) * t321 + Ifges(7,2) * t320 + Ifges(7,6) * t331;
t300 = Ifges(7,1) * t321 + Ifges(7,4) * t320 + Ifges(7,5) * t331;
t353 = mrSges(7,1) * t276 - mrSges(7,2) * t277 + Ifges(7,5) * t304 + Ifges(7,6) * t303 + Ifges(7,3) * t317 + t299 * t321 - t300 * t320;
t314 = (Ifges(6,5) * qJD(5)) + (-Ifges(6,1) * t346 - Ifges(6,4) * t349) * qJD(2);
t313 = (Ifges(6,6) * qJD(5)) + (-Ifges(6,4) * t346 - Ifges(6,2) * t349) * qJD(2);
t298 = Ifges(7,5) * t321 + Ifges(7,6) * t320 + Ifges(7,3) * t331;
t271 = mrSges(7,2) * t278 - mrSges(7,3) * t276 + Ifges(7,1) * t304 + Ifges(7,4) * t303 + Ifges(7,5) * t317 + t298 * t320 - t299 * t331;
t270 = -mrSges(7,1) * t278 + mrSges(7,3) * t277 + Ifges(7,4) * t304 + Ifges(7,2) * t303 + Ifges(7,6) * t317 - t298 * t321 + t300 * t331;
t1 = [m(2) * t336 + t344 * ((m(3) + m(4)) * t359 - t357) + (t347 * (m(3) * t297 - qJDD(2) * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t352 + t358) + t350 * (m(3) * t296 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t352 - t263)) * t341; -pkin(2) * t263 + qJ(3) * (-mrSges(4,1) * t352 + t358) + mrSges(3,1) * t296 - mrSges(3,2) * t297 - pkin(3) * t361 - mrSges(4,1) * t295 + mrSges(4,3) * t294 - pkin(8) * t364 - mrSges(5,1) * t286 + mrSges(5,2) * t287 - t346 * (mrSges(6,2) * t284 - mrSges(6,3) * t281 + Ifges(6,1) * t324 + Ifges(6,4) * t325 + Ifges(6,5) * qJDD(5) - pkin(9) * t269 - qJD(5) * t313 - t270 * t345 + t271 * t348) - t349 * (-mrSges(6,1) * t284 + mrSges(6,3) * t282 + Ifges(6,4) * t324 + Ifges(6,2) * t325 + Ifges(6,6) * qJDD(5) - pkin(5) * t269 + qJD(5) * t314 - t353) - pkin(4) * t354 + (Ifges(3,3) + Ifges(4,2) + Ifges(5,3)) * qJDD(2); t263; t357; Ifges(6,5) * t324 + Ifges(6,6) * t325 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t281 - mrSges(6,2) * t282 + t345 * t271 + t348 * t270 + pkin(5) * t356 + pkin(9) * t363 + (-t313 * t346 + t314 * t349) * qJD(2); t353;];
tauJ  = t1;
