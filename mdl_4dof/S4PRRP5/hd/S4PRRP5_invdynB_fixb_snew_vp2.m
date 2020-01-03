% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:55
% EndTime: 2019-12-31 16:28:56
% DurationCPUTime: 0.75s
% Computational Cost: add. (3492->167), mult. (6636->204), div. (0->0), fcn. (3207->6), ass. (0->72)
t363 = Ifges(4,1) + Ifges(5,1);
t356 = Ifges(4,4) + Ifges(5,4);
t355 = Ifges(4,5) + Ifges(5,5);
t362 = Ifges(4,2) + Ifges(5,2);
t361 = Ifges(4,6) + Ifges(5,6);
t360 = Ifges(4,3) + Ifges(5,3);
t334 = qJD(2) ^ 2;
t329 = sin(pkin(6));
t353 = cos(pkin(6));
t318 = -t353 * g(1) - t329 * g(2);
t328 = -g(3) + qJDD(1);
t331 = sin(qJ(2));
t333 = cos(qJ(2));
t294 = -t331 * t318 + t333 * t328;
t336 = -qJDD(2) * pkin(2) - t294;
t292 = -t334 * pkin(5) + t336;
t330 = sin(qJ(3));
t332 = cos(qJ(3));
t345 = qJD(2) * qJD(3);
t341 = t332 * t345;
t314 = t330 * qJDD(2) + t341;
t315 = t332 * qJDD(2) - t330 * t345;
t346 = qJD(2) * t332;
t323 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t346;
t347 = qJD(2) * t330;
t319 = qJD(3) * pkin(3) - qJ(4) * t347;
t327 = t332 ^ 2;
t288 = t319 * t347 - t315 * pkin(3) + qJDD(4) + (-qJ(4) * t327 - pkin(5)) * t334 + t336;
t322 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t346;
t337 = m(5) * t288 - t315 * mrSges(5,1) - t322 * t346;
t357 = -mrSges(4,2) - mrSges(5,2);
t359 = -m(4) * t292 + t315 * mrSges(4,1) + t357 * t314 + t323 * t346 - t337;
t358 = pkin(3) * t334;
t295 = t333 * t318 + t331 * t328;
t293 = -t334 * pkin(2) + qJDD(2) * pkin(5) + t295;
t317 = t329 * g(1) - t353 * g(2);
t290 = t332 * t293 - t330 * t317;
t313 = (-mrSges(4,1) * t332 + mrSges(4,2) * t330) * qJD(2);
t344 = qJD(2) * qJD(4);
t287 = t315 * qJ(4) - qJD(3) * t319 - t327 * t358 + 0.2e1 * t332 * t344 + t290;
t312 = (-mrSges(5,1) * t332 + mrSges(5,2) * t330) * qJD(2);
t342 = m(5) * t287 + t315 * mrSges(5,3) + t312 * t346;
t320 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t347;
t348 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t347 - t320;
t282 = m(4) * t290 + t315 * mrSges(4,3) + t348 * qJD(3) + t357 * qJDD(3) + t313 * t346 + t342;
t280 = t332 * t282;
t310 = t332 * t317;
t289 = -t330 * t293 - t310;
t286 = qJDD(3) * pkin(3) - t310 + (-t314 + t341) * qJ(4) + (t332 * t358 - t293 - 0.2e1 * t344) * t330;
t343 = m(5) * t286 + qJDD(3) * mrSges(5,1) + qJD(3) * t322;
t281 = m(4) * t289 + qJDD(3) * mrSges(4,1) + qJD(3) * t323 + (-mrSges(4,3) - mrSges(5,3)) * t314 + (-t312 - t313) * t347 + t343;
t274 = m(3) * t295 - t334 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t330 * t281 + t280;
t338 = qJD(2) * t348;
t279 = m(3) * t294 + qJDD(2) * mrSges(3,1) - t334 * mrSges(3,2) + t330 * t338 + t359;
t339 = t333 * t274 - t331 * t279;
t268 = m(2) * t318 + t339;
t277 = t332 * t281 + t330 * t282;
t276 = (m(2) + m(3)) * t317 - t277;
t352 = t329 * t268 + t353 * t276;
t269 = t331 * t274 + t333 * t279;
t351 = t360 * qJD(3) + (t355 * t330 + t361 * t332) * qJD(2);
t350 = -t361 * qJD(3) + (-t356 * t330 - t362 * t332) * qJD(2);
t349 = t355 * qJD(3) + (t363 * t330 + t356 * t332) * qJD(2);
t340 = t353 * t268 - t329 * t276;
t283 = -t314 * mrSges(5,3) - t312 * t347 + t343;
t271 = mrSges(4,2) * t292 + mrSges(5,2) * t288 - mrSges(4,3) * t289 - mrSges(5,3) * t286 - qJ(4) * t283 + t350 * qJD(3) + t355 * qJDD(3) + t363 * t314 + t356 * t315 + t351 * t346;
t270 = -mrSges(4,1) * t292 + mrSges(4,3) * t290 - mrSges(5,1) * t288 + mrSges(5,3) * t287 - pkin(3) * t337 + qJ(4) * t342 + t362 * t315 + (-pkin(3) * mrSges(5,2) + t356) * t314 + (-qJ(4) * mrSges(5,2) + t361) * qJDD(3) + (-qJ(4) * t320 + t349) * qJD(3) + (-pkin(3) * t320 - t351) * t347;
t265 = mrSges(3,1) * t317 - mrSges(4,1) * t289 - mrSges(5,1) * t286 + mrSges(4,2) * t290 + mrSges(5,2) * t287 + mrSges(3,3) * t295 + t334 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t277 - pkin(3) * t283 - t361 * t315 - t355 * t314 - t360 * qJDD(3) + (t350 * t330 + t349 * t332) * qJD(2);
t264 = -mrSges(3,2) * t317 - mrSges(3,3) * t294 + Ifges(3,5) * qJDD(2) - t334 * Ifges(3,6) - pkin(5) * t277 - t330 * t270 + t332 * t271;
t263 = -mrSges(2,1) * t328 + mrSges(2,3) * t318 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t294 + mrSges(3,2) * t295 - t332 * t270 - pkin(2) * t359 - pkin(5) * t280 - pkin(1) * t269 + (-pkin(2) * t338 + pkin(5) * t281 - t271) * t330;
t262 = mrSges(2,2) * t328 - mrSges(2,3) * t317 - pkin(4) * t269 + t333 * t264 - t331 * t265;
t1 = [-m(1) * g(1) + t340; -m(1) * g(2) + t352; -m(1) * g(3) + m(2) * t328 + t269; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t352 + t353 * t262 - t329 * t263; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t340 + t329 * t262 + t353 * t263; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t317 - mrSges(2,2) * t318 + t331 * t264 + t333 * t265 + pkin(1) * (m(3) * t317 - t277) + pkin(4) * t339;];
tauB = t1;
