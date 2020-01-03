% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:27
% EndTime: 2019-12-31 17:42:29
% DurationCPUTime: 1.53s
% Computational Cost: add. (22135->175), mult. (27758->222), div. (0->0), fcn. (16760->10), ass. (0->77)
t369 = m(3) + m(4);
t368 = cos(pkin(8));
t345 = qJD(2) + qJD(3);
t350 = sin(qJ(5));
t367 = t345 * t350;
t353 = cos(qJ(5));
t366 = t345 * t353;
t348 = sin(pkin(8));
t339 = -t368 * g(1) - t348 * g(2);
t346 = -g(3) + qJDD(1);
t352 = sin(qJ(2));
t355 = cos(qJ(2));
t323 = -t352 * t339 + t355 * t346;
t321 = qJDD(2) * pkin(2) + t323;
t324 = t355 * t339 + t352 * t346;
t356 = qJD(2) ^ 2;
t322 = -t356 * pkin(2) + t324;
t351 = sin(qJ(3));
t354 = cos(qJ(3));
t316 = t354 * t321 - t351 * t322;
t344 = qJDD(2) + qJDD(3);
t314 = t344 * pkin(3) + t316;
t317 = t351 * t321 + t354 * t322;
t343 = t345 ^ 2;
t315 = -t343 * pkin(3) + t317;
t347 = sin(pkin(9));
t349 = cos(pkin(9));
t311 = t347 * t314 + t349 * t315;
t309 = -t343 * pkin(4) + t344 * pkin(7) + t311;
t338 = t348 * g(1) - t368 * g(2);
t337 = qJDD(4) - t338;
t306 = -t350 * t309 + t353 * t337;
t330 = (-mrSges(6,1) * t353 + mrSges(6,2) * t350) * t345;
t364 = qJD(5) * t345;
t331 = t350 * t344 + t353 * t364;
t336 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t366;
t304 = m(6) * t306 + qJDD(5) * mrSges(6,1) - t331 * mrSges(6,3) + qJD(5) * t336 - t330 * t367;
t307 = t353 * t309 + t350 * t337;
t332 = t353 * t344 - t350 * t364;
t335 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t367;
t305 = m(6) * t307 - qJDD(5) * mrSges(6,2) + t332 * mrSges(6,3) - qJD(5) * t335 + t330 * t366;
t358 = -t350 * t304 + t353 * t305;
t293 = m(5) * t311 - t343 * mrSges(5,1) - t344 * mrSges(5,2) + t358;
t310 = t349 * t314 - t347 * t315;
t308 = -t344 * pkin(4) - t343 * pkin(7) - t310;
t357 = -m(6) * t308 + t332 * mrSges(6,1) - t331 * mrSges(6,2) - t335 * t367 + t336 * t366;
t300 = m(5) * t310 + t344 * mrSges(5,1) - t343 * mrSges(5,2) + t357;
t290 = t347 * t293 + t349 * t300;
t287 = m(4) * t316 + t344 * mrSges(4,1) - t343 * mrSges(4,2) + t290;
t359 = t349 * t293 - t347 * t300;
t288 = m(4) * t317 - t343 * mrSges(4,1) - t344 * mrSges(4,2) + t359;
t282 = t354 * t287 + t351 * t288;
t280 = m(3) * t323 + qJDD(2) * mrSges(3,1) - t356 * mrSges(3,2) + t282;
t360 = -t351 * t287 + t354 * t288;
t281 = m(3) * t324 - t356 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t360;
t361 = -t352 * t280 + t355 * t281;
t273 = m(2) * t339 + t361;
t296 = t353 * t304 + t350 * t305;
t363 = m(5) * t337 + t296;
t295 = (m(2) + t369) * t338 - t363;
t365 = t348 * t273 + t368 * t295;
t274 = t355 * t280 + t352 * t281;
t362 = t368 * t273 - t348 * t295;
t327 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t350 + Ifges(6,4) * t353) * t345;
t326 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t350 + Ifges(6,2) * t353) * t345;
t325 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t350 + Ifges(6,6) * t353) * t345;
t298 = mrSges(6,2) * t308 - mrSges(6,3) * t306 + Ifges(6,1) * t331 + Ifges(6,4) * t332 + Ifges(6,5) * qJDD(5) - qJD(5) * t326 + t325 * t366;
t297 = -mrSges(6,1) * t308 + mrSges(6,3) * t307 + Ifges(6,4) * t331 + Ifges(6,2) * t332 + Ifges(6,6) * qJDD(5) + qJD(5) * t327 - t325 * t367;
t289 = -mrSges(5,1) * t337 - mrSges(6,1) * t306 + mrSges(6,2) * t307 + mrSges(5,3) * t311 + t343 * Ifges(5,5) - Ifges(6,5) * t331 + Ifges(5,6) * t344 - Ifges(6,6) * t332 - Ifges(6,3) * qJDD(5) - pkin(4) * t296 + (-t326 * t350 + t327 * t353) * t345;
t283 = mrSges(5,2) * t337 - mrSges(5,3) * t310 + Ifges(5,5) * t344 - t343 * Ifges(5,6) - pkin(7) * t296 - t350 * t297 + t353 * t298;
t276 = -mrSges(4,2) * t338 - mrSges(4,3) * t316 + Ifges(4,5) * t344 - t343 * Ifges(4,6) - qJ(4) * t290 + t349 * t283 - t347 * t289;
t275 = mrSges(4,1) * t338 + mrSges(4,3) * t317 + t343 * Ifges(4,5) + Ifges(4,6) * t344 - pkin(3) * t363 + qJ(4) * t359 + t347 * t283 + t349 * t289;
t270 = -pkin(1) * t274 - pkin(2) * t282 + mrSges(3,2) * t324 - mrSges(3,1) * t323 - pkin(3) * t290 - mrSges(4,1) * t316 + mrSges(4,2) * t317 - pkin(4) * t357 - pkin(7) * t358 - mrSges(5,1) * t310 + mrSges(5,2) * t311 - t350 * t298 - t353 * t297 + mrSges(2,3) * t339 - Ifges(3,3) * qJDD(2) - mrSges(2,1) * t346 + (-Ifges(4,3) - Ifges(5,3)) * t344;
t269 = -mrSges(3,2) * t338 - mrSges(3,3) * t323 + Ifges(3,5) * qJDD(2) - t356 * Ifges(3,6) - pkin(6) * t282 - t351 * t275 + t354 * t276;
t268 = Ifges(3,6) * qJDD(2) + t356 * Ifges(3,5) + mrSges(3,1) * t338 + mrSges(3,3) * t324 + t351 * t276 + t354 * t275 - pkin(2) * (-m(4) * t338 + t363) + pkin(6) * t360;
t267 = mrSges(2,2) * t346 - mrSges(2,3) * t338 - pkin(5) * t274 - t352 * t268 + t355 * t269;
t1 = [-m(1) * g(1) + t362; -m(1) * g(2) + t365; -m(1) * g(3) + m(2) * t346 + t274; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t365 + t368 * t267 - t348 * t270; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t362 + t348 * t267 + t368 * t270; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t339 + t352 * t269 + t355 * t268 - pkin(1) * t363 + pkin(5) * t361 + (pkin(1) * t369 + mrSges(2,1)) * t338;];
tauB = t1;
