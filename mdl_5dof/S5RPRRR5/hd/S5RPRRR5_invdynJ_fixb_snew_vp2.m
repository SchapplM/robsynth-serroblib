% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:59
% EndTime: 2020-01-03 11:54:00
% DurationCPUTime: 0.68s
% Computational Cost: add. (6034->168), mult. (8145->221), div. (0->0), fcn. (4613->10), ass. (0->78)
t310 = qJD(1) + qJD(3);
t316 = sin(qJ(4));
t337 = t310 * t316;
t320 = cos(qJ(4));
t336 = t310 * t320;
t318 = sin(qJ(1));
t322 = cos(qJ(1));
t329 = -g(2) * t322 - g(3) * t318;
t296 = qJDD(1) * pkin(1) + t329;
t323 = qJD(1) ^ 2;
t332 = -g(2) * t318 + t322 * g(3);
t297 = -pkin(1) * t323 + t332;
t313 = sin(pkin(9));
t314 = cos(pkin(9));
t279 = t314 * t296 - t297 * t313;
t277 = qJDD(1) * pkin(2) + t279;
t280 = t313 * t296 + t314 * t297;
t278 = -pkin(2) * t323 + t280;
t317 = sin(qJ(3));
t321 = cos(qJ(3));
t262 = t317 * t277 + t321 * t278;
t306 = t310 ^ 2;
t308 = qJDD(1) + qJDD(3);
t259 = -pkin(3) * t306 + pkin(7) * t308 + t262;
t312 = -g(1) + qJDD(2);
t255 = -t316 * t259 + t320 * t312;
t334 = qJD(4) * t310;
t333 = t320 * t334;
t291 = t308 * t316 + t333;
t252 = (-t291 + t333) * pkin(8) + (t306 * t316 * t320 + qJDD(4)) * pkin(4) + t255;
t256 = t320 * t259 + t316 * t312;
t292 = t308 * t320 - t316 * t334;
t300 = qJD(4) * pkin(4) - pkin(8) * t337;
t311 = t320 ^ 2;
t253 = -pkin(4) * t306 * t311 + pkin(8) * t292 - qJD(4) * t300 + t256;
t315 = sin(qJ(5));
t319 = cos(qJ(5));
t250 = t252 * t319 - t253 * t315;
t286 = (-t315 * t316 + t319 * t320) * t310;
t268 = qJD(5) * t286 + t291 * t319 + t292 * t315;
t287 = (t315 * t320 + t316 * t319) * t310;
t273 = -mrSges(6,1) * t286 + mrSges(6,2) * t287;
t309 = qJD(4) + qJD(5);
t281 = -mrSges(6,2) * t309 + mrSges(6,3) * t286;
t307 = qJDD(4) + qJDD(5);
t247 = m(6) * t250 + mrSges(6,1) * t307 - mrSges(6,3) * t268 - t273 * t287 + t281 * t309;
t251 = t252 * t315 + t253 * t319;
t267 = -qJD(5) * t287 - t291 * t315 + t292 * t319;
t282 = mrSges(6,1) * t309 - mrSges(6,3) * t287;
t248 = m(6) * t251 - mrSges(6,2) * t307 + mrSges(6,3) * t267 + t273 * t286 - t282 * t309;
t238 = t319 * t247 + t315 * t248;
t290 = (-mrSges(5,1) * t320 + mrSges(5,2) * t316) * t310;
t299 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t336;
t236 = m(5) * t255 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t291 + qJD(4) * t299 - t290 * t337 + t238;
t298 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t337;
t330 = -t247 * t315 + t319 * t248;
t237 = m(5) * t256 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t292 - qJD(4) * t298 + t290 * t336 + t330;
t331 = -t236 * t316 + t320 * t237;
t233 = m(4) * t262 - mrSges(4,1) * t306 - mrSges(4,2) * t308 + t331;
t261 = t321 * t277 - t317 * t278;
t328 = -t308 * pkin(3) - t261;
t258 = -pkin(7) * t306 + t328;
t254 = t300 * t337 - t292 * pkin(4) + (-pkin(8) * t311 - pkin(7)) * t306 + t328;
t326 = m(6) * t254 - t267 * mrSges(6,1) + mrSges(6,2) * t268 - t286 * t281 + t282 * t287;
t324 = -m(5) * t258 + t292 * mrSges(5,1) - mrSges(5,2) * t291 - t298 * t337 + t299 * t336 - t326;
t242 = m(4) * t261 + mrSges(4,1) * t308 - mrSges(4,2) * t306 + t324;
t335 = t317 * t233 + t321 * t242;
t269 = Ifges(6,5) * t287 + Ifges(6,6) * t286 + Ifges(6,3) * t309;
t271 = Ifges(6,1) * t287 + Ifges(6,4) * t286 + Ifges(6,5) * t309;
t239 = -mrSges(6,1) * t254 + mrSges(6,3) * t251 + Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * t307 - t269 * t287 + t271 * t309;
t270 = Ifges(6,4) * t287 + Ifges(6,2) * t286 + Ifges(6,6) * t309;
t240 = mrSges(6,2) * t254 - mrSges(6,3) * t250 + Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t307 + t269 * t286 - t270 * t309;
t283 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t316 + Ifges(5,6) * t320) * t310;
t284 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t316 + Ifges(5,2) * t320) * t310;
t285 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t316 + Ifges(5,4) * t320) * t310;
t327 = -mrSges(4,2) * t262 + t320 * (-mrSges(5,1) * t258 + mrSges(5,3) * t256 + Ifges(5,4) * t291 + Ifges(5,2) * t292 + Ifges(5,6) * qJDD(4) - pkin(4) * t326 + pkin(8) * t330 + qJD(4) * t285 + t319 * t239 + t315 * t240 - t283 * t337) + t316 * (mrSges(5,2) * t258 - mrSges(5,3) * t255 + Ifges(5,1) * t291 + Ifges(5,4) * t292 + Ifges(5,5) * qJDD(4) - pkin(8) * t238 - qJD(4) * t284 - t239 * t315 + t240 * t319 + t283 * t336) + pkin(7) * t331 + pkin(3) * t324 + mrSges(4,1) * t261 + Ifges(4,3) * t308;
t325 = mrSges(6,1) * t250 - mrSges(6,2) * t251 + Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t307 + t287 * t270 - t286 * t271;
t1 = [mrSges(2,1) * t329 - mrSges(2,2) * t332 + pkin(1) * (t313 * (m(3) * t280 - mrSges(3,1) * t323 - qJDD(1) * mrSges(3,2) + t233 * t321 - t242 * t317) + t314 * (m(3) * t279 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t323 + t335)) + pkin(2) * t335 + mrSges(3,1) * t279 - mrSges(3,2) * t280 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t327; t320 * t236 + t316 * t237 + (m(3) + m(4)) * t312; t327; mrSges(5,1) * t255 - mrSges(5,2) * t256 + Ifges(5,5) * t291 + Ifges(5,6) * t292 + Ifges(5,3) * qJDD(4) + pkin(4) * t238 + (t284 * t316 - t285 * t320) * t310 + t325; t325;];
tauJ = t1;
