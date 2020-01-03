% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:13
% EndTime: 2020-01-03 11:36:14
% DurationCPUTime: 0.71s
% Computational Cost: add. (4202->132), mult. (5880->179), div. (0->0), fcn. (3273->10), ass. (0->74)
t326 = 2 * qJD(4);
t299 = sin(qJ(1));
t302 = cos(qJ(1));
t309 = -t302 * g(2) - t299 * g(3);
t278 = qJDD(1) * pkin(1) + t309;
t303 = qJD(1) ^ 2;
t312 = -t299 * g(2) + t302 * g(3);
t279 = -t303 * pkin(1) + t312;
t294 = sin(pkin(8));
t296 = cos(pkin(8));
t265 = t296 * t278 - t294 * t279;
t260 = qJDD(1) * pkin(2) + t265;
t266 = t294 * t278 + t296 * t279;
t261 = -t303 * pkin(2) + t266;
t298 = sin(qJ(3));
t301 = cos(qJ(3));
t256 = t298 * t260 + t301 * t261;
t291 = (qJD(1) + qJD(3));
t289 = t291 ^ 2;
t290 = qJDD(1) + qJDD(3);
t253 = -t289 * pkin(3) + t290 * qJ(4) + t256;
t325 = (t291 * t326) + t253;
t293 = sin(pkin(9));
t324 = mrSges(5,2) * t293;
t322 = mrSges(5,3) * t290;
t321 = t291 * t293;
t295 = cos(pkin(9));
t320 = t295 * t290;
t319 = t295 * t291;
t292 = -g(1) + qJDD(2);
t318 = t295 * t292;
t249 = t293 * t292 + t325 * t295;
t310 = -pkin(4) * t295 - pkin(7) * t293;
t274 = t310 * t291;
t247 = t274 * t319 + t249;
t255 = t301 * t260 - t298 * t261;
t307 = -t289 * qJ(4) + qJDD(4) - t255;
t250 = (-pkin(3) + t310) * t290 + t307;
t297 = sin(qJ(5));
t300 = cos(qJ(5));
t244 = -t297 * t247 + t300 * t250;
t281 = qJD(5) - t319;
t314 = t297 * t321;
t267 = -t281 * mrSges(6,2) - mrSges(6,3) * t314;
t269 = (t297 * mrSges(6,1) + t300 * mrSges(6,2)) * t321;
t315 = qJD(5) * t291;
t271 = (t290 * t300 - t297 * t315) * t293;
t280 = qJDD(5) - t320;
t313 = t300 * t321;
t242 = m(6) * t244 + t280 * mrSges(6,1) - t271 * mrSges(6,3) + t281 * t267 - t269 * t313;
t245 = t300 * t247 + t297 * t250;
t268 = t281 * mrSges(6,1) - mrSges(6,3) * t313;
t270 = (-t290 * t297 - t300 * t315) * t293;
t243 = m(6) * t245 - t280 * mrSges(6,2) + t270 * mrSges(6,3) - t281 * t268 - t269 * t314;
t272 = (-mrSges(5,1) * t295 + t324) * t291;
t237 = m(5) * t249 - t297 * t242 + t300 * t243 + (t272 * t291 + t322) * t295;
t246 = -t318 + (t253 + (t326 + t274) * t291) * t293;
t248 = -t325 * t293 + t318;
t241 = m(5) * t248 - m(6) * t246 + t270 * mrSges(6,1) - t271 * mrSges(6,2) + (-t322 + (-t267 * t297 - t268 * t300 - t272) * t291) * t293;
t311 = t295 * t237 - t293 * t241;
t232 = m(4) * t256 - t289 * mrSges(4,1) - t290 * mrSges(4,2) + t311;
t240 = t300 * t242 + t297 * t243;
t252 = -t290 * pkin(3) + t307;
t305 = -m(5) * t252 + mrSges(5,1) * t320 - t240 + (t293 ^ 2 + t295 ^ 2) * mrSges(5,3) * t289;
t235 = m(4) * t255 - t289 * mrSges(4,2) + (mrSges(4,1) - t324) * t290 + t305;
t317 = t298 * t232 + t301 * t235;
t263 = Ifges(6,6) * t281 + (t300 * Ifges(6,4) - t297 * Ifges(6,2)) * t321;
t264 = Ifges(6,5) * t281 + (t300 * Ifges(6,1) - t297 * Ifges(6,4)) * t321;
t308 = t300 * t263 + t297 * t264;
t239 = t290 * t324 - t305;
t273 = (Ifges(5,5) * t293 + Ifges(5,6) * t295) * t291;
t304 = mrSges(6,1) * t244 - mrSges(6,2) * t245 + Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t280;
t306 = -mrSges(4,2) * t256 + t293 * (t273 * t319 + mrSges(5,2) * t252 - mrSges(5,3) * t248 + t300 * (mrSges(6,2) * t246 - mrSges(6,3) * t244 + Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t280 - t281 * t263) - t297 * (-mrSges(6,1) * t246 + mrSges(6,3) * t245 + Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t280 + t281 * t264) - pkin(7) * t240 + (Ifges(5,1) * t293 + Ifges(5,4) * t295) * t290) + t295 * (Ifges(5,2) * t320 - mrSges(5,1) * t252 + mrSges(5,3) * t249 - pkin(4) * t240 + (Ifges(5,4) * t290 + (-t273 - t308) * t291) * t293 - t304) + qJ(4) * t311 - pkin(3) * t239 + mrSges(4,1) * t255 + Ifges(4,3) * t290;
t1 = [mrSges(2,1) * t309 - mrSges(2,2) * t312 + pkin(1) * (t294 * (m(3) * t266 - t303 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t301 * t232 - t298 * t235) + t296 * (m(3) * t265 + qJDD(1) * mrSges(3,1) - t303 * mrSges(3,2) + t317)) + pkin(2) * t317 + mrSges(3,1) * t265 - mrSges(3,2) * t266 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t306; t293 * t237 + t295 * t241 + (m(3) + m(4)) * t292; t306; t239; t308 * t321 + t304;];
tauJ = t1;
