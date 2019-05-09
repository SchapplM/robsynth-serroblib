% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-05-05 17:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:15:10
% EndTime: 2019-05-05 17:15:23
% DurationCPUTime: 6.67s
% Computational Cost: add. (75600->345), mult. (166670->403), div. (0->0), fcn. (102783->8), ass. (0->134)
t301 = sin(qJ(1));
t304 = cos(qJ(1));
t281 = t301 * g(1) - t304 * g(2);
t306 = qJD(1) ^ 2;
t324 = -t306 * qJ(2) + qJDD(2) - t281;
t350 = -pkin(1) - pkin(7);
t250 = t350 * qJDD(1) + t324;
t300 = sin(qJ(3));
t303 = cos(qJ(3));
t236 = t300 * g(3) + t303 * t250;
t334 = qJD(1) * qJD(3);
t332 = t300 * t334;
t276 = qJDD(1) * t303 - t332;
t202 = (-t276 - t332) * qJ(4) + (-t300 * t303 * t306 + qJDD(3)) * pkin(3) + t236;
t237 = -g(3) * t303 + t300 * t250;
t275 = -qJDD(1) * t300 - t303 * t334;
t337 = qJD(1) * t303;
t279 = qJD(3) * pkin(3) - qJ(4) * t337;
t295 = t300 ^ 2;
t203 = -pkin(3) * t295 * t306 + qJ(4) * t275 - qJD(3) * t279 + t237;
t298 = cos(pkin(9));
t338 = qJD(1) * t300;
t343 = sin(pkin(9));
t263 = t298 * t337 - t343 * t338;
t355 = -2 * qJD(4);
t187 = t202 * t298 - t343 * t203 + t263 * t355;
t262 = (t298 * t300 + t343 * t303) * qJD(1);
t256 = t262 * t355;
t342 = t343 * t202 + t298 * t203;
t188 = t256 + t342;
t216 = Ifges(5,4) * t263 - Ifges(5,2) * t262 + Ifges(5,6) * qJD(3);
t223 = -mrSges(6,2) * t262 - mrSges(6,3) * t263;
t233 = -t298 * t275 + t343 * t276;
t234 = t343 * t275 + t298 * t276;
t247 = mrSges(6,1) * t262 - qJD(3) * mrSges(6,3);
t221 = pkin(4) * t262 - qJ(5) * t263;
t305 = qJD(3) ^ 2;
t183 = -qJDD(3) * pkin(4) - qJ(5) * t305 + t263 * t221 + qJDD(5) - t187;
t336 = qJD(3) * t262;
t177 = (t262 * t263 - qJDD(3)) * pkin(8) + (t234 + t336) * pkin(5) + t183;
t249 = pkin(5) * t263 - qJD(3) * pkin(8);
t261 = t262 ^ 2;
t282 = -t304 * g(1) - t301 * g(2);
t327 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t282;
t205 = -pkin(3) * t275 + qJDD(4) + t279 * t337 + (-qJ(4) * t295 + t350) * t306 + t327;
t351 = -2 * qJD(5);
t310 = (-t234 + t336) * qJ(5) + t205 + (qJD(3) * pkin(4) + t351) * t263;
t180 = (pkin(4) + pkin(8)) * t233 - pkin(5) * t261 - t249 * t263 + t310;
t299 = sin(qJ(6));
t302 = cos(qJ(6));
t174 = t177 * t302 - t180 * t299;
t238 = -qJD(3) * t299 + t262 * t302;
t198 = qJD(6) * t238 + qJDD(3) * t302 + t233 * t299;
t239 = qJD(3) * t302 + t262 * t299;
t208 = -mrSges(7,1) * t238 + mrSges(7,2) * t239;
t259 = qJD(6) + t263;
t211 = -mrSges(7,2) * t259 + mrSges(7,3) * t238;
t232 = qJDD(6) + t234;
t171 = m(7) * t174 + mrSges(7,1) * t232 - t198 * mrSges(7,3) - t208 * t239 + t211 * t259;
t175 = t177 * t299 + t180 * t302;
t197 = -qJD(6) * t239 - qJDD(3) * t299 + t233 * t302;
t212 = mrSges(7,1) * t259 - mrSges(7,3) * t239;
t172 = m(7) * t175 - mrSges(7,2) * t232 + t197 * mrSges(7,3) + t208 * t238 - t212 * t259;
t160 = t171 * t302 + t172 * t299;
t325 = t305 * pkin(4) - qJDD(3) * qJ(5) - t342;
t179 = -pkin(5) * t233 - pkin(8) * t261 - t262 * t221 + t256 + ((2 * qJD(5)) + t249) * qJD(3) - t325;
t190 = Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t259;
t192 = Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t259;
t163 = -mrSges(7,1) * t179 + mrSges(7,3) * t175 + Ifges(7,4) * t198 + Ifges(7,2) * t197 + Ifges(7,6) * t232 - t190 * t239 + t192 * t259;
t191 = Ifges(7,4) * t239 + Ifges(7,2) * t238 + Ifges(7,6) * t259;
t164 = mrSges(7,2) * t179 - mrSges(7,3) * t174 + Ifges(7,1) * t198 + Ifges(7,4) * t197 + Ifges(7,5) * t232 + t190 * t238 - t191 * t259;
t181 = qJD(3) * t351 + ((2 * qJD(4)) + t221) * t262 + t325;
t213 = Ifges(6,5) * qJD(3) - Ifges(6,6) * t263 + Ifges(6,3) * t262;
t316 = -mrSges(6,2) * t183 + mrSges(6,3) * t181 - Ifges(6,1) * qJDD(3) + Ifges(6,4) * t234 - Ifges(6,5) * t233 + pkin(8) * t160 + t299 * t163 - t302 * t164 + t263 * t213;
t176 = -m(7) * t179 + t197 * mrSges(7,1) - t198 * mrSges(7,2) + t211 * t238 - t239 * t212;
t248 = mrSges(6,1) * t263 + qJD(3) * mrSges(6,2);
t317 = -m(6) * t181 + qJDD(3) * mrSges(6,3) + qJD(3) * t248 - t176;
t322 = -m(6) * t183 - t234 * mrSges(6,1) - t263 * t223 - t160;
t215 = Ifges(6,4) * qJD(3) - Ifges(6,2) * t263 + Ifges(6,6) * t262;
t340 = Ifges(5,1) * t263 - Ifges(5,4) * t262 + Ifges(5,5) * qJD(3) - t215;
t357 = -mrSges(5,2) * t188 + pkin(4) * (-qJDD(3) * mrSges(6,2) - qJD(3) * t247 + t322) + qJ(5) * (-mrSges(6,1) * t233 - t223 * t262 + t317) + mrSges(5,1) * t187 + t263 * t216 - Ifges(5,6) * t233 + Ifges(5,5) * t234 + Ifges(5,3) * qJDD(3) - t316 + t340 * t262;
t222 = mrSges(5,1) * t262 + mrSges(5,2) * t263;
t339 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t262 - t247;
t348 = mrSges(5,1) - mrSges(6,2);
t156 = m(5) * t187 - mrSges(5,3) * t234 + t339 * qJD(3) + t348 * qJDD(3) - t222 * t263 + t322;
t246 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t263;
t167 = m(5) * t188 - qJDD(3) * mrSges(5,2) - qJD(3) * t246 + (-t222 - t223) * t262 + (-mrSges(5,3) - mrSges(6,1)) * t233 + t317;
t151 = t298 * t156 + t343 * t167;
t265 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t303 - Ifges(4,2) * t300) * qJD(1);
t266 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t303 - Ifges(4,4) * t300) * qJD(1);
t356 = mrSges(4,1) * t236 - mrSges(4,2) * t237 + Ifges(4,5) * t276 + Ifges(4,6) * t275 + Ifges(4,3) * qJDD(3) + pkin(3) * t151 + t265 * t337 + t266 * t338 + t357;
t244 = t350 * t306 + t327;
t278 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t338;
t280 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t337;
t161 = -t171 * t299 + t302 * t172;
t185 = pkin(4) * t233 + t310;
t326 = -m(6) * t185 + t234 * mrSges(6,3) + t263 * t248 - t161;
t312 = m(5) * t205 + t234 * mrSges(5,2) + t348 * t233 + t263 * t246 + t339 * t262 - t326;
t354 = -m(4) * t244 + mrSges(4,1) * t275 - t276 * mrSges(4,2) - t278 * t338 - t280 * t337 - t312;
t274 = (mrSges(4,1) * t300 + mrSges(4,2) * t303) * qJD(1);
t148 = m(4) * t236 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t276 + qJD(3) * t278 - t274 * t337 + t151;
t329 = -t343 * t156 + t298 * t167;
t149 = m(4) * t237 - qJDD(3) * mrSges(4,2) + t275 * mrSges(4,3) - qJD(3) * t280 - t274 * t338 + t329;
t144 = t148 * t303 + t149 * t300;
t260 = -qJDD(1) * pkin(1) + t324;
t353 = mrSges(3,1) * t260 + pkin(2) * t144 + t356;
t349 = mrSges(2,1) - mrSges(3,2);
t347 = Ifges(5,4) + Ifges(6,6);
t346 = Ifges(2,5) - Ifges(3,4);
t345 = -Ifges(2,6) + Ifges(3,5);
t217 = Ifges(6,1) * qJD(3) - Ifges(6,4) * t263 + Ifges(6,5) * t262;
t341 = -Ifges(5,5) * t263 + Ifges(5,6) * t262 - Ifges(5,3) * qJD(3) - t217;
t145 = -t148 * t300 + t303 * t149;
t323 = -m(3) * t260 + t306 * mrSges(3,3) - t144;
t320 = mrSges(7,1) * t174 - mrSges(7,2) * t175 + Ifges(7,5) * t198 + Ifges(7,6) * t197 + Ifges(7,3) * t232 + t239 * t191 - t238 * t192;
t157 = -mrSges(6,2) * t233 - t247 * t262 - t326;
t315 = -mrSges(6,1) * t181 + mrSges(6,2) * t185 - pkin(5) * t176 - pkin(8) * t161 - t302 * t163 - t299 * t164;
t140 = -mrSges(5,1) * t205 + mrSges(5,3) * t188 - pkin(4) * t157 + t341 * t263 + t347 * t234 + (-Ifges(5,2) - Ifges(6,3)) * t233 + (Ifges(5,6) - Ifges(6,5)) * qJDD(3) + t340 * qJD(3) + t315;
t313 = mrSges(6,1) * t183 - mrSges(6,3) * t185 + pkin(5) * t160 + t320;
t146 = t313 + t341 * t262 + (Ifges(5,1) + Ifges(6,2)) * t234 - t347 * t233 + (Ifges(5,5) - Ifges(6,4)) * qJDD(3) + (-t216 + t213) * qJD(3) + mrSges(5,2) * t205 - mrSges(5,3) * t187 - qJ(5) * t157;
t264 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t303 - Ifges(4,6) * t300) * qJD(1);
t137 = -mrSges(4,1) * t244 + mrSges(4,3) * t237 + Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * qJDD(3) - pkin(3) * t312 + qJ(4) * t329 + qJD(3) * t266 + t298 * t140 + t343 * t146 - t264 * t337;
t139 = mrSges(4,2) * t244 - mrSges(4,3) * t236 + Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * qJDD(3) - qJ(4) * t151 - qJD(3) * t265 - t343 * t140 + t298 * t146 - t264 * t338;
t253 = pkin(1) * t306 - t327;
t319 = mrSges(3,2) * t260 - mrSges(3,3) * t253 + Ifges(3,1) * qJDD(1) - pkin(7) * t144 - t137 * t300 + t303 * t139;
t318 = -mrSges(3,1) * t253 - pkin(2) * t354 - pkin(7) * t145 - t137 * t303 - t139 * t300;
t309 = -m(3) * t253 + t306 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t354;
t314 = -mrSges(2,2) * t282 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t323) + qJ(2) * t309 + mrSges(2,1) * t281 + Ifges(2,3) * qJDD(1) + t319;
t152 = m(2) * t282 - mrSges(2,1) * t306 - qJDD(1) * mrSges(2,2) + t309;
t143 = -m(3) * g(3) + t145;
t141 = m(2) * t281 - mrSges(2,2) * t306 + t349 * qJDD(1) + t323;
t136 = -qJ(2) * t143 + t345 * t306 + t346 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t281 + t353;
t135 = mrSges(2,3) * t282 - pkin(1) * t143 + t349 * g(3) - t345 * qJDD(1) + t346 * t306 + t318;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t136 - t301 * t135 - pkin(6) * (t141 * t304 + t152 * t301), t136, t319, t139, t146, -t262 * t215 - t316, t164; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t136 + t304 * t135 + pkin(6) * (-t141 * t301 + t152 * t304), t135, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t306 * Ifges(3,5) - t353, t137, t140, Ifges(6,4) * qJDD(3) - Ifges(6,2) * t234 + Ifges(6,6) * t233 - qJD(3) * t213 + t262 * t217 - t313, t163; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t314, t314, mrSges(3,2) * g(3) + Ifges(3,4) * t306 + Ifges(3,5) * qJDD(1) - t318, t356, t357, Ifges(6,5) * qJDD(3) - Ifges(6,6) * t234 + Ifges(6,3) * t233 + qJD(3) * t215 + t263 * t217 - t315, t320;];
m_new  = t1;
