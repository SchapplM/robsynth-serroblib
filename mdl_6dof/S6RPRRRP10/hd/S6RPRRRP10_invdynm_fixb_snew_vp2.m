% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:55:38
% EndTime: 2019-05-06 01:55:52
% DurationCPUTime: 6.34s
% Computational Cost: add. (117886->337), mult. (227978->397), div. (0->0), fcn. (146192->8), ass. (0->128)
t298 = sin(qJ(1));
t301 = cos(qJ(1));
t279 = -t301 * g(1) - t298 * g(2);
t336 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t279;
t335 = (-pkin(1) - pkin(7));
t334 = cos(qJ(5));
t333 = mrSges(2,1) - mrSges(3,2);
t332 = -mrSges(6,3) - mrSges(7,2);
t331 = -Ifges(3,4) + Ifges(2,5);
t330 = (Ifges(3,5) - Ifges(2,6));
t303 = qJD(1) ^ 2;
t244 = (t303 * t335) - t336;
t300 = cos(qJ(3));
t326 = qJD(1) * qJD(3);
t283 = t300 * t326;
t297 = sin(qJ(3));
t273 = -qJDD(1) * t297 - t283;
t323 = t297 * t326;
t274 = qJDD(1) * t300 - t323;
t214 = (-t274 + t323) * pkin(8) + (-t273 + t283) * pkin(3) + t244;
t278 = t298 * g(1) - g(2) * t301;
t317 = -t303 * qJ(2) + qJDD(2) - t278;
t245 = qJDD(1) * t335 + t317;
t238 = -g(3) * t300 + t245 * t297;
t272 = (pkin(3) * t297 - pkin(8) * t300) * qJD(1);
t285 = t297 * qJD(1);
t302 = qJD(3) ^ 2;
t219 = -pkin(3) * t302 + qJDD(3) * pkin(8) - t272 * t285 + t238;
t296 = sin(qJ(4));
t299 = cos(qJ(4));
t184 = t214 * t299 - t296 * t219;
t327 = qJD(1) * t300;
t269 = qJD(3) * t299 - t296 * t327;
t232 = qJD(4) * t269 + qJDD(3) * t296 + t274 * t299;
t268 = qJDD(4) - t273;
t270 = qJD(3) * t296 + t299 * t327;
t282 = t285 + qJD(4);
t180 = (t269 * t282 - t232) * pkin(9) + (t269 * t270 + t268) * pkin(4) + t184;
t185 = t214 * t296 + t219 * t299;
t231 = -qJD(4) * t270 + qJDD(3) * t299 - t274 * t296;
t243 = pkin(4) * t282 - pkin(9) * t270;
t267 = t269 ^ 2;
t182 = -pkin(4) * t267 + pkin(9) * t231 - t243 * t282 + t185;
t295 = sin(qJ(5));
t178 = t180 * t295 + t182 * t334;
t234 = t269 * t295 + t270 * t334;
t197 = qJD(5) * t234 - t231 * t334 + t232 * t295;
t281 = qJD(5) + t282;
t222 = mrSges(6,1) * t281 - mrSges(6,3) * t234;
t233 = -t269 * t334 + t270 * t295;
t261 = qJDD(5) + t268;
t209 = pkin(5) * t233 - qJ(6) * t234;
t280 = t281 ^ 2;
t171 = -pkin(5) * t280 + qJ(6) * t261 + 0.2e1 * qJD(6) * t281 - t209 * t233 + t178;
t223 = -mrSges(7,1) * t281 + mrSges(7,2) * t234;
t324 = m(7) * t171 + mrSges(7,3) * t261 + t223 * t281;
t210 = mrSges(7,1) * t233 - mrSges(7,3) * t234;
t328 = -mrSges(6,1) * t233 - mrSges(6,2) * t234 - t210;
t161 = m(6) * t178 - t261 * mrSges(6,2) + t197 * t332 - t281 * t222 + t233 * t328 + t324;
t177 = t180 * t334 - t182 * t295;
t198 = -qJD(5) * t233 + t231 * t295 + t232 * t334;
t221 = -mrSges(6,2) * t281 - mrSges(6,3) * t233;
t173 = -pkin(5) * t261 - qJ(6) * t280 + t209 * t234 + qJDD(6) - t177;
t220 = -mrSges(7,2) * t233 + mrSges(7,3) * t281;
t320 = -m(7) * t173 + mrSges(7,1) * t261 + t220 * t281;
t163 = m(6) * t177 + t261 * mrSges(6,1) + t198 * t332 + t281 * t221 + t234 * t328 + t320;
t155 = t161 * t295 + t163 * t334;
t236 = -mrSges(5,1) * t269 + mrSges(5,2) * t270;
t239 = -mrSges(5,2) * t282 + mrSges(5,3) * t269;
t151 = m(5) * t184 + mrSges(5,1) * t268 - mrSges(5,3) * t232 - t236 * t270 + t239 * t282 + t155;
t240 = mrSges(5,1) * t282 - mrSges(5,3) * t270;
t321 = t161 * t334 - t163 * t295;
t152 = m(5) * t185 - mrSges(5,2) * t268 + mrSges(5,3) * t231 + t236 * t269 - t240 * t282 + t321;
t146 = t151 * t299 + t152 * t296;
t202 = Ifges(7,4) * t234 + Ifges(7,2) * t281 + Ifges(7,6) * t233;
t329 = -Ifges(6,5) * t234 + Ifges(6,6) * t233 - Ifges(6,3) * t281 - t202;
t271 = (mrSges(4,1) * t297 + mrSges(4,2) * t300) * qJD(1);
t277 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t327;
t322 = -t151 * t296 + t152 * t299;
t144 = m(4) * t238 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t273 - qJD(3) * t277 - t271 * t285 + t322;
t237 = t297 * g(3) + t300 * t245;
t276 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t285;
t218 = -qJDD(3) * pkin(3) - t302 * pkin(8) + t272 * t327 - t237;
t183 = -t231 * pkin(4) - t267 * pkin(9) + t243 * t270 + t218;
t175 = -0.2e1 * qJD(6) * t234 + (t233 * t281 - t198) * qJ(6) + (t234 * t281 + t197) * pkin(5) + t183;
t164 = m(7) * t175 + mrSges(7,1) * t197 - mrSges(7,3) * t198 + t220 * t233 - t223 * t234;
t309 = m(6) * t183 + mrSges(6,1) * t197 + mrSges(6,2) * t198 + t221 * t233 + t222 * t234 + t164;
t305 = -m(5) * t218 + mrSges(5,1) * t231 - mrSges(5,2) * t232 + t239 * t269 - t240 * t270 - t309;
t156 = m(4) * t237 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t274 + qJD(3) * t276 - t271 * t327 + t305;
t139 = t144 * t300 - t156 * t297;
t319 = -mrSges(7,1) * t175 + mrSges(7,2) * t171;
t138 = t297 * t144 + t300 * t156;
t200 = Ifges(7,5) * t234 + Ifges(7,6) * t281 + Ifges(7,3) * t233;
t316 = mrSges(7,2) * t173 - mrSges(7,3) * t175 + Ifges(7,1) * t198 + Ifges(7,4) * t261 + Ifges(7,5) * t197 + t200 * t281;
t250 = -qJDD(1) * pkin(1) + t317;
t315 = -m(3) * t250 + (mrSges(3,3) * t303) - t138;
t142 = -m(4) * t244 + mrSges(4,1) * t273 - mrSges(4,2) * t274 - t276 * t285 - t277 * t327 - t146;
t204 = Ifges(7,1) * t234 + Ifges(7,4) * t281 + Ifges(7,5) * t233;
t314 = mrSges(7,1) * t173 - mrSges(7,3) * t171 - Ifges(7,4) * t198 - Ifges(7,2) * t261 - Ifges(7,6) * t197 + t200 * t234 - t204 * t233;
t205 = Ifges(6,1) * t234 - Ifges(6,4) * t233 + Ifges(6,5) * t281;
t153 = -mrSges(6,1) * t183 + mrSges(6,3) * t178 - pkin(5) * t164 + (t204 + t205) * t281 + (Ifges(6,6) - Ifges(7,6)) * t261 + t329 * t234 + (Ifges(6,4) - Ifges(7,5)) * t198 + (-Ifges(6,2) - Ifges(7,3)) * t197 + t319;
t203 = Ifges(6,4) * t234 - Ifges(6,2) * t233 + Ifges(6,6) * t281;
t154 = mrSges(6,2) * t183 - mrSges(6,3) * t177 + Ifges(6,1) * t198 - Ifges(6,4) * t197 + Ifges(6,5) * t261 - qJ(6) * t164 - t281 * t203 + t233 * t329 + t316;
t225 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * t282;
t227 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * t282;
t132 = -mrSges(5,1) * t218 + mrSges(5,3) * t185 + Ifges(5,4) * t232 + Ifges(5,2) * t231 + Ifges(5,6) * t268 - pkin(4) * t309 + pkin(9) * t321 + t153 * t334 + t295 * t154 - t270 * t225 + t282 * t227;
t226 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * t282;
t134 = mrSges(5,2) * t218 - mrSges(5,3) * t184 + Ifges(5,1) * t232 + Ifges(5,4) * t231 + Ifges(5,5) * t268 - pkin(9) * t155 - t153 * t295 + t154 * t334 + t225 * t269 - t226 * t282;
t258 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t300 - Ifges(4,6) * t297) * qJD(1);
t259 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t300 - Ifges(4,2) * t297) * qJD(1);
t129 = mrSges(4,2) * t244 - mrSges(4,3) * t237 + Ifges(4,1) * t274 + Ifges(4,4) * t273 + Ifges(4,5) * qJDD(3) - pkin(8) * t146 - qJD(3) * t259 - t132 * t296 + t134 * t299 - t258 * t285;
t260 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t300 - Ifges(4,4) * t297) * qJD(1);
t306 = mrSges(6,2) * t178 - t233 * t205 - qJ(6) * (-t197 * mrSges(7,2) - t233 * t210 + t324) - pkin(5) * (-t198 * mrSges(7,2) - t234 * t210 + t320) - mrSges(6,1) * t177 + Ifges(6,6) * t197 - Ifges(6,5) * t198 - t234 * t203 - Ifges(6,3) * t261 + t314;
t304 = mrSges(5,1) * t184 - mrSges(5,2) * t185 + Ifges(5,5) * t232 + Ifges(5,6) * t231 + Ifges(5,3) * t268 + pkin(4) * t155 + t226 * t270 - t227 * t269 - t306;
t130 = -mrSges(4,1) * t244 + mrSges(4,3) * t238 + Ifges(4,4) * t274 + Ifges(4,2) * t273 + Ifges(4,6) * qJDD(3) - pkin(3) * t146 + qJD(3) * t260 - t258 * t327 - t304;
t248 = t303 * pkin(1) + t336;
t313 = mrSges(3,2) * t250 - mrSges(3,3) * t248 + Ifges(3,1) * qJDD(1) - pkin(7) * t138 + t129 * t300 - t130 * t297;
t312 = -mrSges(3,1) * t248 - pkin(2) * t142 - pkin(7) * t139 - t297 * t129 - t300 * t130;
t311 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t274 + Ifges(4,6) * t273 + Ifges(4,3) * qJDD(3) + pkin(3) * t305 + pkin(8) * t322 + t132 * t299 + t134 * t296 + t259 * t327 + t260 * t285;
t310 = -m(3) * t248 + mrSges(3,2) * t303 + qJDD(1) * mrSges(3,3) - t142;
t308 = -mrSges(2,2) * t279 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t315) + qJ(2) * t310 + mrSges(2,1) * t278 + Ifges(2,3) * qJDD(1) + t313;
t307 = mrSges(3,1) * t250 + pkin(2) * t138 + t311;
t140 = m(2) * t279 - mrSges(2,1) * t303 - qJDD(1) * mrSges(2,2) + t310;
t137 = -m(3) * g(3) + t139;
t135 = m(2) * t278 - t303 * mrSges(2,2) + qJDD(1) * t333 + t315;
t127 = -mrSges(2,3) * t278 - qJ(2) * t137 + (t330 * t303) + t331 * qJDD(1) + t307 + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t126 = mrSges(2,3) * t279 - pkin(1) * t137 + g(3) * t333 - qJDD(1) * t330 + t303 * t331 + t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t301 * t127 - t298 * t126 - pkin(6) * (t135 * t301 + t140 * t298), t127, t313, t129, t134, t154, -t202 * t233 + t316; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t127 + t301 * t126 + pkin(6) * (-t135 * t298 + t140 * t301), t126, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t303 * Ifges(3,5)) - t307, t130, t132, t153, -t314; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t308, t308, mrSges(3,2) * g(3) + t303 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t312, t311, t304, -t306, Ifges(7,5) * t198 + Ifges(7,6) * t261 + Ifges(7,3) * t197 + t234 * t202 - t281 * t204 - t319;];
m_new  = t1;
