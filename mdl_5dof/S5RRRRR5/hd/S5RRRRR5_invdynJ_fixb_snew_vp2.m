% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:15
% EndTime: 2020-01-03 12:13:15
% DurationCPUTime: 0.67s
% Computational Cost: add. (10318->167), mult. (10575->220), div. (0->0), fcn. (5990->10), ass. (0->79)
t311 = qJD(1) + qJD(2);
t305 = qJD(3) + t311;
t303 = t305 ^ 2;
t338 = pkin(4) * t303;
t314 = sin(qJ(4));
t337 = t305 * t314;
t319 = cos(qJ(4));
t336 = t305 * t319;
t317 = sin(qJ(1));
t322 = cos(qJ(1));
t329 = -t322 * g(2) - t317 * g(3);
t298 = qJDD(1) * pkin(1) + t329;
t332 = -t317 * g(2) + t322 * g(3);
t299 = -qJD(1) ^ 2 * pkin(1) + t332;
t316 = sin(qJ(2));
t321 = cos(qJ(2));
t278 = t321 * t298 - t316 * t299;
t309 = qJDD(1) + qJDD(2);
t275 = t309 * pkin(2) + t278;
t279 = t316 * t298 + t321 * t299;
t307 = t311 ^ 2;
t276 = -t307 * pkin(2) + t279;
t315 = sin(qJ(3));
t320 = cos(qJ(3));
t260 = t315 * t275 + t320 * t276;
t304 = qJDD(3) + t309;
t257 = -t303 * pkin(3) + t304 * pkin(8) + t260;
t335 = t314 * t257;
t333 = qJD(4) * t305;
t290 = t314 * t304 + t319 * t333;
t250 = qJDD(4) * pkin(4) - t290 * pkin(9) - t335 + (pkin(9) * t333 + t314 * t338 - g(1)) * t319;
t254 = -t314 * g(1) + t319 * t257;
t291 = t319 * t304 - t314 * t333;
t297 = qJD(4) * pkin(4) - pkin(9) * t337;
t312 = t319 ^ 2;
t251 = t291 * pkin(9) - qJD(4) * t297 - t312 * t338 + t254;
t313 = sin(qJ(5));
t318 = cos(qJ(5));
t248 = t318 * t250 - t313 * t251;
t285 = (-t313 * t314 + t318 * t319) * t305;
t266 = t285 * qJD(5) + t318 * t290 + t313 * t291;
t286 = (t313 * t319 + t314 * t318) * t305;
t271 = -t285 * mrSges(6,1) + t286 * mrSges(6,2);
t310 = qJD(4) + qJD(5);
t280 = -t310 * mrSges(6,2) + t285 * mrSges(6,3);
t308 = qJDD(4) + qJDD(5);
t245 = m(6) * t248 + t308 * mrSges(6,1) - t266 * mrSges(6,3) - t286 * t271 + t310 * t280;
t249 = t313 * t250 + t318 * t251;
t265 = -t286 * qJD(5) - t313 * t290 + t318 * t291;
t281 = t310 * mrSges(6,1) - t286 * mrSges(6,3);
t246 = m(6) * t249 - t308 * mrSges(6,2) + t265 * mrSges(6,3) + t285 * t271 - t310 * t281;
t236 = t318 * t245 + t313 * t246;
t253 = -t319 * g(1) - t335;
t289 = (-mrSges(5,1) * t319 + mrSges(5,2) * t314) * t305;
t295 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t337;
t296 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t336;
t330 = -t313 * t245 + t318 * t246;
t331 = -t314 * (m(5) * t253 + qJDD(4) * mrSges(5,1) - t290 * mrSges(5,3) + qJD(4) * t296 - t289 * t337 + t236) + t319 * (m(5) * t254 - qJDD(4) * mrSges(5,2) + t291 * mrSges(5,3) - qJD(4) * t295 + t289 * t336 + t330);
t232 = m(4) * t260 - t303 * mrSges(4,1) - t304 * mrSges(4,2) + t331;
t259 = t320 * t275 - t315 * t276;
t328 = -t304 * pkin(3) - t259;
t256 = -t303 * pkin(8) + t328;
t252 = t297 * t337 - t291 * pkin(4) + (-pkin(9) * t312 - pkin(8)) * t303 + t328;
t326 = m(6) * t252 - t265 * mrSges(6,1) + t266 * mrSges(6,2) - t285 * t280 + t286 * t281;
t323 = -m(5) * t256 + t291 * mrSges(5,1) - t290 * mrSges(5,2) - t295 * t337 + t296 * t336 - t326;
t240 = m(4) * t259 + t304 * mrSges(4,1) - t303 * mrSges(4,2) + t323;
t334 = t315 * t232 + t320 * t240;
t267 = Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * t310;
t269 = Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * t310;
t237 = -mrSges(6,1) * t252 + mrSges(6,3) * t249 + Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t308 - t286 * t267 + t310 * t269;
t268 = Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * t310;
t238 = mrSges(6,2) * t252 - mrSges(6,3) * t248 + Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t308 + t285 * t267 - t310 * t268;
t282 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t314 + Ifges(5,6) * t319) * t305;
t283 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t314 + Ifges(5,2) * t319) * t305;
t284 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t314 + Ifges(5,4) * t319) * t305;
t327 = -mrSges(4,2) * t260 + t319 * (-mrSges(5,1) * t256 + mrSges(5,3) * t254 + Ifges(5,4) * t290 + Ifges(5,2) * t291 + Ifges(5,6) * qJDD(4) - pkin(4) * t326 + pkin(9) * t330 + qJD(4) * t284 + t318 * t237 + t313 * t238 - t282 * t337) + t314 * (mrSges(5,2) * t256 - mrSges(5,3) * t253 + Ifges(5,1) * t290 + Ifges(5,4) * t291 + Ifges(5,5) * qJDD(4) - pkin(9) * t236 - qJD(4) * t283 - t313 * t237 + t318 * t238 + t282 * t336) + pkin(8) * t331 + pkin(3) * t323 + mrSges(4,1) * t259 + Ifges(4,3) * t304;
t325 = mrSges(6,1) * t248 - mrSges(6,2) * t249 + Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t308 + t286 * t268 - t285 * t269;
t324 = mrSges(3,1) * t278 - mrSges(3,2) * t279 + Ifges(3,3) * t309 + pkin(2) * t334 + t327;
t1 = [t324 + pkin(1) * (t316 * (m(3) * t279 - t307 * mrSges(3,1) - t309 * mrSges(3,2) + t320 * t232 - t315 * t240) + t321 * (m(3) * t278 + t309 * mrSges(3,1) - t307 * mrSges(3,2) + t334)) + mrSges(2,1) * t329 - mrSges(2,2) * t332 + Ifges(2,3) * qJDD(1); t324; t327; mrSges(5,1) * t253 - mrSges(5,2) * t254 + Ifges(5,5) * t290 + Ifges(5,6) * t291 + Ifges(5,3) * qJDD(4) + pkin(4) * t236 + (t283 * t314 - t284 * t319) * t305 + t325; t325;];
tauJ = t1;
