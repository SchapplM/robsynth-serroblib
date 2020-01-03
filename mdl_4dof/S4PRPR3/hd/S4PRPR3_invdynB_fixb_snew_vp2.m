% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRPR3
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:52
% EndTime: 2019-12-31 16:20:52
% DurationCPUTime: 0.90s
% Computational Cost: add. (7600->163), mult. (16408->213), div. (0->0), fcn. (10403->8), ass. (0->77)
t337 = qJD(2) ^ 2;
t331 = cos(pkin(7));
t359 = pkin(3) * t331;
t329 = sin(pkin(7));
t358 = mrSges(4,2) * t329;
t327 = t331 ^ 2;
t357 = t327 * t337;
t330 = sin(pkin(6));
t332 = cos(pkin(6));
t316 = g(1) * t330 - g(2) * t332;
t317 = -g(1) * t332 - g(2) * t330;
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t306 = t334 * t316 + t336 * t317;
t304 = -pkin(2) * t337 + qJDD(2) * qJ(3) + t306;
t328 = -g(3) + qJDD(1);
t353 = qJD(2) * qJD(3);
t355 = t331 * t328 - 0.2e1 * t329 * t353;
t289 = (-pkin(5) * qJDD(2) + t337 * t359 - t304) * t329 + t355;
t293 = t329 * t328 + (t304 + 0.2e1 * t353) * t331;
t352 = qJDD(2) * t331;
t290 = -pkin(3) * t357 + pkin(5) * t352 + t293;
t333 = sin(qJ(4));
t335 = cos(qJ(4));
t287 = t289 * t335 - t290 * t333;
t341 = -t329 * t333 + t331 * t335;
t309 = t341 * qJD(2);
t342 = t329 * t335 + t331 * t333;
t310 = t342 * qJD(2);
t299 = -mrSges(5,1) * t309 + mrSges(5,2) * t310;
t303 = qJD(4) * t309 + qJDD(2) * t342;
t307 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t309;
t285 = m(5) * t287 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t303 + qJD(4) * t307 - t299 * t310;
t288 = t289 * t333 + t290 * t335;
t302 = -qJD(4) * t310 + qJDD(2) * t341;
t308 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t310;
t286 = m(5) * t288 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t302 - qJD(4) * t308 + t299 * t309;
t277 = t335 * t285 + t333 * t286;
t292 = -t304 * t329 + t355;
t340 = mrSges(4,3) * qJDD(2) + t337 * (-mrSges(4,1) * t331 + t358);
t275 = m(4) * t292 - t329 * t340 + t277;
t347 = -t285 * t333 + t335 * t286;
t276 = m(4) * t293 + t331 * t340 + t347;
t348 = -t275 * t329 + t331 * t276;
t270 = m(3) * t306 - mrSges(3,1) * t337 - qJDD(2) * mrSges(3,2) + t348;
t305 = t316 * t336 - t334 * t317;
t343 = qJDD(3) - t305;
t301 = -qJDD(2) * pkin(2) - qJ(3) * t337 + t343;
t326 = t329 ^ 2;
t291 = (-pkin(2) - t359) * qJDD(2) + (-qJ(3) + (-t326 - t327) * pkin(5)) * t337 + t343;
t339 = m(5) * t291 - t302 * mrSges(5,1) + mrSges(5,2) * t303 - t309 * t307 + t308 * t310;
t338 = -m(4) * t301 + mrSges(4,1) * t352 - t339 + (t326 * t337 + t357) * mrSges(4,3);
t281 = m(3) * t305 - mrSges(3,2) * t337 + (mrSges(3,1) - t358) * qJDD(2) + t338;
t266 = t334 * t270 + t336 * t281;
t264 = m(2) * t316 + t266;
t349 = t336 * t270 - t281 * t334;
t265 = m(2) * t317 + t349;
t356 = t332 * t264 + t330 * t265;
t271 = t331 * t275 + t329 * t276;
t344 = Ifges(4,5) * t329 + Ifges(4,6) * t331;
t354 = t337 * t344;
t351 = m(3) * t328 + t271;
t350 = -t264 * t330 + t332 * t265;
t346 = Ifges(4,1) * t329 + Ifges(4,4) * t331;
t345 = Ifges(4,4) * t329 + Ifges(4,2) * t331;
t296 = Ifges(5,1) * t310 + Ifges(5,4) * t309 + Ifges(5,5) * qJD(4);
t295 = Ifges(5,4) * t310 + Ifges(5,2) * t309 + Ifges(5,6) * qJD(4);
t294 = Ifges(5,5) * t310 + Ifges(5,6) * t309 + Ifges(5,3) * qJD(4);
t279 = mrSges(5,2) * t291 - mrSges(5,3) * t287 + Ifges(5,1) * t303 + Ifges(5,4) * t302 + Ifges(5,5) * qJDD(4) - qJD(4) * t295 + t294 * t309;
t278 = -mrSges(5,1) * t291 + mrSges(5,3) * t288 + Ifges(5,4) * t303 + Ifges(5,2) * t302 + Ifges(5,6) * qJDD(4) + qJD(4) * t296 - t294 * t310;
t267 = mrSges(4,2) * t301 - mrSges(4,3) * t292 - pkin(5) * t277 + qJDD(2) * t346 - t278 * t333 + t279 * t335 + t331 * t354;
t260 = -mrSges(4,1) * t301 + mrSges(4,3) * t293 - pkin(3) * t339 + pkin(5) * t347 + qJDD(2) * t345 + t335 * t278 + t333 * t279 - t329 * t354;
t259 = -mrSges(3,1) * t328 - mrSges(4,1) * t292 - mrSges(5,1) * t287 + mrSges(4,2) * t293 + mrSges(5,2) * t288 + mrSges(3,3) * t306 - Ifges(5,5) * t303 - Ifges(5,6) * t302 - Ifges(5,3) * qJDD(4) - pkin(2) * t271 - pkin(3) * t277 - t310 * t295 + t309 * t296 + (Ifges(3,6) - t344) * qJDD(2) + (-t329 * t345 + t331 * t346 + Ifges(3,5)) * t337;
t258 = mrSges(3,2) * t328 - mrSges(3,3) * t305 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t337 - qJ(3) * t271 - t260 * t329 + t267 * t331;
t257 = mrSges(2,2) * t328 - mrSges(2,3) * t316 - pkin(4) * t266 + t258 * t336 - t259 * t334;
t256 = -mrSges(2,1) * t328 + mrSges(2,3) * t317 - pkin(1) * t351 + pkin(4) * t349 + t334 * t258 + t336 * t259;
t1 = [-m(1) * g(1) + t350; -m(1) * g(2) + t356; -m(1) * g(3) + m(2) * t328 + t351; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t356 - t330 * t256 + t332 * t257; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t350 + t332 * t256 + t330 * t257; pkin(1) * t266 + mrSges(2,1) * t316 - mrSges(2,2) * t317 + qJ(3) * t348 + mrSges(3,1) * t305 - mrSges(3,2) * t306 + t329 * t267 + t331 * t260 + pkin(2) * (-qJDD(2) * t358 + t338) - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2);];
tauB = t1;
