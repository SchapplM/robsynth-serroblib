% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRPR3
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:57
% EndTime: 2019-12-05 15:05:00
% DurationCPUTime: 1.23s
% Computational Cost: add. (13073->162), mult. (20018->210), div. (0->0), fcn. (13408->10), ass. (0->73)
t300 = sin(pkin(7));
t303 = cos(pkin(7));
t293 = -t303 * g(1) - t300 * g(2);
t297 = -g(3) + qJDD(1);
t299 = sin(pkin(8));
t302 = cos(pkin(8));
t280 = t302 * t293 + t299 * t297;
t292 = t300 * g(1) - t303 * g(2);
t291 = qJDD(2) - t292;
t305 = sin(qJ(3));
t307 = cos(qJ(3));
t275 = -t305 * t280 + t307 * t291;
t273 = qJDD(3) * pkin(3) + t275;
t276 = t307 * t280 + t305 * t291;
t308 = qJD(3) ^ 2;
t274 = -t308 * pkin(3) + t276;
t298 = sin(pkin(9));
t301 = cos(pkin(9));
t270 = t298 * t273 + t301 * t274;
t268 = -t308 * pkin(4) + qJDD(3) * pkin(6) + t270;
t279 = t299 * t293 - t302 * t297;
t278 = qJDD(4) + t279;
t304 = sin(qJ(5));
t306 = cos(qJ(5));
t265 = -t304 * t268 + t306 * t278;
t288 = (-mrSges(6,1) * t306 + mrSges(6,2) * t304) * qJD(3);
t317 = qJD(3) * qJD(5);
t289 = t304 * qJDD(3) + t306 * t317;
t318 = qJD(3) * t306;
t295 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t318;
t319 = qJD(3) * t304;
t263 = m(6) * t265 + qJDD(5) * mrSges(6,1) - t289 * mrSges(6,3) + qJD(5) * t295 - t288 * t319;
t266 = t306 * t268 + t304 * t278;
t290 = t306 * qJDD(3) - t304 * t317;
t294 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t319;
t264 = m(6) * t266 - qJDD(5) * mrSges(6,2) + t290 * mrSges(6,3) - qJD(5) * t294 + t288 * t318;
t312 = -t304 * t263 + t306 * t264;
t252 = m(5) * t270 - t308 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t312;
t269 = t301 * t273 - t298 * t274;
t267 = -qJDD(3) * pkin(4) - t308 * pkin(6) - t269;
t309 = -m(6) * t267 + t290 * mrSges(6,1) - t289 * mrSges(6,2) - t294 * t319 + t295 * t318;
t259 = m(5) * t269 + qJDD(3) * mrSges(5,1) - t308 * mrSges(5,2) + t309;
t249 = t298 * t252 + t301 * t259;
t247 = m(4) * t275 + qJDD(3) * mrSges(4,1) - t308 * mrSges(4,2) + t249;
t313 = t301 * t252 - t298 * t259;
t248 = m(4) * t276 - t308 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t313;
t314 = -t305 * t247 + t307 * t248;
t242 = m(3) * t280 + t314;
t255 = t306 * t263 + t304 * t264;
t311 = m(5) * t278 + t255;
t254 = (-m(3) - m(4)) * t279 - t311;
t315 = t302 * t242 - t299 * t254;
t236 = m(2) * t293 + t315;
t243 = t307 * t247 + t305 * t248;
t310 = -m(3) * t291 - t243;
t241 = m(2) * t292 + t310;
t320 = t300 * t236 + t303 * t241;
t237 = t299 * t242 + t302 * t254;
t316 = t303 * t236 - t300 * t241;
t283 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t304 + Ifges(6,4) * t306) * qJD(3);
t282 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t304 + Ifges(6,2) * t306) * qJD(3);
t281 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t304 + Ifges(6,6) * t306) * qJD(3);
t257 = mrSges(6,2) * t267 - mrSges(6,3) * t265 + Ifges(6,1) * t289 + Ifges(6,4) * t290 + Ifges(6,5) * qJDD(5) - qJD(5) * t282 + t281 * t318;
t256 = -mrSges(6,1) * t267 + mrSges(6,3) * t266 + Ifges(6,4) * t289 + Ifges(6,2) * t290 + Ifges(6,6) * qJDD(5) + qJD(5) * t283 - t281 * t319;
t245 = -mrSges(5,1) * t278 - mrSges(6,1) * t265 + mrSges(6,2) * t266 + mrSges(5,3) * t270 + t308 * Ifges(5,5) - Ifges(6,5) * t289 + Ifges(5,6) * qJDD(3) - Ifges(6,6) * t290 - Ifges(6,3) * qJDD(5) - pkin(4) * t255 + (-t282 * t304 + t283 * t306) * qJD(3);
t244 = mrSges(5,2) * t278 - mrSges(5,3) * t269 + Ifges(5,5) * qJDD(3) - t308 * Ifges(5,6) - pkin(6) * t255 - t304 * t256 + t306 * t257;
t233 = mrSges(4,2) * t279 - mrSges(4,3) * t275 + Ifges(4,5) * qJDD(3) - t308 * Ifges(4,6) - qJ(4) * t249 + t301 * t244 - t298 * t245;
t232 = -mrSges(4,1) * t279 + mrSges(4,3) * t276 + t308 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t311 + qJ(4) * t313 + t298 * t244 + t301 * t245;
t231 = -pkin(2) * t243 + mrSges(3,3) * t280 - mrSges(3,1) * t291 - pkin(3) * t249 - mrSges(4,1) * t275 + mrSges(4,2) * t276 - t304 * t257 - t306 * t256 - pkin(4) * t309 - pkin(6) * t312 - mrSges(5,1) * t269 + mrSges(5,2) * t270 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3);
t230 = mrSges(3,2) * t291 + mrSges(3,3) * t279 - pkin(5) * t243 - t305 * t232 + t307 * t233;
t229 = -mrSges(2,1) * t297 + mrSges(2,3) * t293 + mrSges(3,1) * t279 + mrSges(3,2) * t280 - t305 * t233 - t307 * t232 - pkin(2) * (-m(4) * t279 - t311) - pkin(5) * t314 - pkin(1) * t237;
t228 = mrSges(2,2) * t297 - mrSges(2,3) * t292 - qJ(2) * t237 + t302 * t230 - t299 * t231;
t1 = [-m(1) * g(1) + t316; -m(1) * g(2) + t320; -m(1) * g(3) + m(2) * t297 + t237; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t320 + t303 * t228 - t300 * t229; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t316 + t300 * t228 + t303 * t229; -mrSges(1,1) * g(2) + mrSges(2,1) * t292 + mrSges(1,2) * g(1) - mrSges(2,2) * t293 + pkin(1) * t310 + qJ(2) * t315 + t299 * t230 + t302 * t231;];
tauB = t1;
