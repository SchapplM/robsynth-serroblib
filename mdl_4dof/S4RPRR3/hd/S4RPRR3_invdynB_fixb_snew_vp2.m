% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:10
% EndTime: 2019-12-31 16:49:11
% DurationCPUTime: 1.06s
% Computational Cost: add. (10437->196), mult. (20169->253), div. (0->0), fcn. (11357->8), ass. (0->80)
t356 = sin(qJ(1));
t359 = cos(qJ(1));
t342 = t356 * g(1) - t359 * g(2);
t334 = qJDD(1) * pkin(1) + t342;
t343 = -t359 * g(1) - t356 * g(2);
t360 = qJD(1) ^ 2;
t336 = -t360 * pkin(1) + t343;
t352 = sin(pkin(7));
t353 = cos(pkin(7));
t321 = t352 * t334 + t353 * t336;
t318 = -t360 * pkin(2) + qJDD(1) * pkin(5) + t321;
t351 = -g(3) + qJDD(2);
t355 = sin(qJ(3));
t358 = cos(qJ(3));
t308 = -t355 * t318 + t358 * t351;
t370 = qJD(1) * qJD(3);
t368 = t358 * t370;
t337 = t355 * qJDD(1) + t368;
t304 = (-t337 + t368) * pkin(6) + (t355 * t358 * t360 + qJDD(3)) * pkin(3) + t308;
t309 = t358 * t318 + t355 * t351;
t338 = t358 * qJDD(1) - t355 * t370;
t372 = qJD(1) * t355;
t341 = qJD(3) * pkin(3) - pkin(6) * t372;
t350 = t358 ^ 2;
t305 = -t350 * t360 * pkin(3) + t338 * pkin(6) - qJD(3) * t341 + t309;
t354 = sin(qJ(4));
t357 = cos(qJ(4));
t302 = t357 * t304 - t354 * t305;
t330 = (-t354 * t355 + t357 * t358) * qJD(1);
t311 = t330 * qJD(4) + t357 * t337 + t354 * t338;
t331 = (t354 * t358 + t355 * t357) * qJD(1);
t319 = -t330 * mrSges(5,1) + t331 * mrSges(5,2);
t349 = qJD(3) + qJD(4);
t322 = -t349 * mrSges(5,2) + t330 * mrSges(5,3);
t348 = qJDD(3) + qJDD(4);
t300 = m(5) * t302 + t348 * mrSges(5,1) - t311 * mrSges(5,3) - t331 * t319 + t349 * t322;
t303 = t354 * t304 + t357 * t305;
t310 = -t331 * qJD(4) - t354 * t337 + t357 * t338;
t323 = t349 * mrSges(5,1) - t331 * mrSges(5,3);
t301 = m(5) * t303 - t348 * mrSges(5,2) + t310 * mrSges(5,3) + t330 * t319 - t349 * t323;
t292 = t357 * t300 + t354 * t301;
t335 = (-mrSges(4,1) * t358 + mrSges(4,2) * t355) * qJD(1);
t371 = qJD(1) * t358;
t340 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t371;
t290 = m(4) * t308 + qJDD(3) * mrSges(4,1) - t337 * mrSges(4,3) + qJD(3) * t340 - t335 * t372 + t292;
t339 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t372;
t364 = -t354 * t300 + t357 * t301;
t291 = m(4) * t309 - qJDD(3) * mrSges(4,2) + t338 * mrSges(4,3) - qJD(3) * t339 + t335 * t371 + t364;
t365 = -t355 * t290 + t358 * t291;
t285 = m(3) * t321 - t360 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t365;
t320 = t353 * t334 - t352 * t336;
t363 = -qJDD(1) * pkin(2) - t320;
t317 = -t360 * pkin(5) + t363;
t306 = t341 * t372 - t338 * pkin(3) + (-pkin(6) * t350 - pkin(5)) * t360 + t363;
t362 = m(5) * t306 - t310 * mrSges(5,1) + t311 * mrSges(5,2) - t330 * t322 + t331 * t323;
t361 = -m(4) * t317 + t338 * mrSges(4,1) - t337 * mrSges(4,2) - t339 * t372 + t340 * t371 - t362;
t296 = m(3) * t320 + qJDD(1) * mrSges(3,1) - t360 * mrSges(3,2) + t361;
t281 = t352 * t285 + t353 * t296;
t279 = m(2) * t342 + qJDD(1) * mrSges(2,1) - t360 * mrSges(2,2) + t281;
t366 = t353 * t285 - t352 * t296;
t280 = m(2) * t343 - t360 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t366;
t373 = t359 * t279 + t356 * t280;
t286 = t358 * t290 + t355 * t291;
t369 = m(3) * t351 + t286;
t367 = -t356 * t279 + t359 * t280;
t329 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t355 + Ifges(4,4) * t358) * qJD(1);
t328 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t355 + Ifges(4,2) * t358) * qJD(1);
t327 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t355 + Ifges(4,6) * t358) * qJD(1);
t315 = Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * t349;
t314 = Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * t349;
t313 = Ifges(5,5) * t331 + Ifges(5,6) * t330 + Ifges(5,3) * t349;
t294 = mrSges(5,2) * t306 - mrSges(5,3) * t302 + Ifges(5,1) * t311 + Ifges(5,4) * t310 + Ifges(5,5) * t348 + t330 * t313 - t349 * t314;
t293 = -mrSges(5,1) * t306 + mrSges(5,3) * t303 + Ifges(5,4) * t311 + Ifges(5,2) * t310 + Ifges(5,6) * t348 - t331 * t313 + t349 * t315;
t282 = mrSges(4,2) * t317 - mrSges(4,3) * t308 + Ifges(4,1) * t337 + Ifges(4,4) * t338 + Ifges(4,5) * qJDD(3) - pkin(6) * t292 - qJD(3) * t328 - t354 * t293 + t357 * t294 + t327 * t371;
t275 = -mrSges(4,1) * t317 + mrSges(4,3) * t309 + Ifges(4,4) * t337 + Ifges(4,2) * t338 + Ifges(4,6) * qJDD(3) - pkin(3) * t362 + pkin(6) * t364 + qJD(3) * t329 + t357 * t293 + t354 * t294 - t327 * t372;
t274 = Ifges(3,6) * qJDD(1) + t360 * Ifges(3,5) - mrSges(3,1) * t351 + mrSges(3,3) * t321 - Ifges(4,5) * t337 - Ifges(4,6) * t338 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t308 + mrSges(4,2) * t309 - Ifges(5,5) * t311 - Ifges(5,6) * t310 - Ifges(5,3) * t348 - t331 * t314 + t330 * t315 - mrSges(5,1) * t302 + mrSges(5,2) * t303 - pkin(3) * t292 - pkin(2) * t286 + (-t355 * t328 + t358 * t329) * qJD(1);
t273 = mrSges(3,2) * t351 - mrSges(3,3) * t320 + Ifges(3,5) * qJDD(1) - t360 * Ifges(3,6) - pkin(5) * t286 - t355 * t275 + t358 * t282;
t272 = -mrSges(2,2) * g(3) - mrSges(2,3) * t342 + Ifges(2,5) * qJDD(1) - t360 * Ifges(2,6) - qJ(2) * t281 + t353 * t273 - t352 * t274;
t271 = mrSges(2,1) * g(3) + mrSges(2,3) * t343 + t360 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t369 + qJ(2) * t366 + t352 * t273 + t353 * t274;
t1 = [-m(1) * g(1) + t367; -m(1) * g(2) + t373; (-m(1) - m(2)) * g(3) + t369; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t373 - t356 * t271 + t359 * t272; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t367 + t359 * t271 + t356 * t272; pkin(1) * t281 + mrSges(2,1) * t342 - mrSges(2,2) * t343 + t355 * t282 + t358 * t275 + pkin(2) * t361 + pkin(5) * t365 + mrSges(3,1) * t320 - mrSges(3,2) * t321 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
