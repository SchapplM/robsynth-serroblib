% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:12:38
% EndTime: 2019-05-06 01:13:02
% DurationCPUTime: 10.88s
% Computational Cost: add. (196614->335), mult. (381830->408), div. (0->0), fcn. (252157->10), ass. (0->128)
t298 = sin(qJ(4));
t302 = cos(qJ(4));
t299 = sin(qJ(3));
t328 = qJD(1) * t299;
t275 = qJD(3) * t298 + t302 * t328;
t303 = cos(qJ(3));
t326 = qJD(1) * qJD(3);
t323 = t303 * t326;
t280 = qJDD(1) * t299 + t323;
t245 = -qJD(4) * t275 + qJDD(3) * t302 - t280 * t298;
t274 = qJD(3) * t302 - t298 * t328;
t246 = qJD(4) * t274 + qJDD(3) * t298 + t280 * t302;
t297 = sin(qJ(5));
t301 = cos(qJ(5));
t249 = t274 * t301 - t275 * t297;
t207 = qJD(5) * t249 + t245 * t297 + t246 * t301;
t250 = t274 * t297 + t275 * t301;
t225 = -mrSges(7,1) * t249 + mrSges(7,2) * t250;
t300 = sin(qJ(1));
t304 = cos(qJ(1));
t285 = t300 * g(1) - g(2) * t304;
t276 = qJDD(1) * pkin(1) + t285;
t286 = -g(1) * t304 - g(2) * t300;
t306 = qJD(1) ^ 2;
t278 = -pkin(1) * t306 + t286;
t295 = sin(pkin(10));
t296 = cos(pkin(10));
t251 = t296 * t276 - t295 * t278;
t237 = -qJDD(1) * pkin(2) - t306 * pkin(7) - t251;
t290 = t299 * t326;
t281 = qJDD(1) * t303 - t290;
t218 = (-t280 - t323) * pkin(8) + (-t281 + t290) * pkin(3) + t237;
t252 = t295 * t276 + t296 * t278;
t238 = -pkin(2) * t306 + qJDD(1) * pkin(7) + t252;
t294 = -g(3) + qJDD(2);
t229 = t303 * t238 + t299 * t294;
t279 = (-pkin(3) * t303 - pkin(8) * t299) * qJD(1);
t305 = qJD(3) ^ 2;
t327 = qJD(1) * t303;
t224 = -pkin(3) * t305 + qJDD(3) * pkin(8) + t279 * t327 + t229;
t189 = t302 * t218 - t298 * t224;
t273 = qJDD(4) - t281;
t288 = qJD(4) - t327;
t185 = (t274 * t288 - t246) * pkin(9) + (t274 * t275 + t273) * pkin(4) + t189;
t190 = t298 * t218 + t302 * t224;
t256 = pkin(4) * t288 - pkin(9) * t275;
t272 = t274 ^ 2;
t187 = -pkin(4) * t272 + pkin(9) * t245 - t256 * t288 + t190;
t178 = t301 * t185 - t297 * t187;
t269 = qJDD(5) + t273;
t287 = qJD(5) + t288;
t173 = -0.2e1 * qJD(6) * t250 + (t249 * t287 - t207) * qJ(6) + (t249 * t250 + t269) * pkin(5) + t178;
t230 = -mrSges(7,2) * t287 + mrSges(7,3) * t249;
t325 = m(7) * t173 + t269 * mrSges(7,1) + t287 * t230;
t170 = -t207 * mrSges(7,3) - t250 * t225 + t325;
t179 = t297 * t185 + t301 * t187;
t206 = -qJD(5) * t250 + t245 * t301 - t246 * t297;
t215 = Ifges(6,4) * t250 + Ifges(6,2) * t249 + Ifges(6,6) * t287;
t216 = Ifges(7,1) * t250 + Ifges(7,4) * t249 + Ifges(7,5) * t287;
t217 = Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t287;
t232 = pkin(5) * t287 - qJ(6) * t250;
t248 = t249 ^ 2;
t176 = -pkin(5) * t248 + qJ(6) * t206 + 0.2e1 * qJD(6) * t249 - t232 * t287 + t179;
t214 = Ifges(7,4) * t250 + Ifges(7,2) * t249 + Ifges(7,6) * t287;
t315 = -mrSges(7,1) * t173 + mrSges(7,2) * t176 - Ifges(7,5) * t207 - Ifges(7,6) * t206 - Ifges(7,3) * t269 - t250 * t214;
t332 = mrSges(6,1) * t178 - mrSges(6,2) * t179 + Ifges(6,5) * t207 + Ifges(6,6) * t206 + Ifges(6,3) * t269 + pkin(5) * t170 + t250 * t215 - t315 + (-t217 - t216) * t249;
t226 = -mrSges(6,1) * t249 + mrSges(6,2) * t250;
t231 = -mrSges(6,2) * t287 + mrSges(6,3) * t249;
t162 = m(6) * t178 + t269 * mrSges(6,1) + t287 * t231 + (-t225 - t226) * t250 + (-mrSges(6,3) - mrSges(7,3)) * t207 + t325;
t233 = mrSges(7,1) * t287 - mrSges(7,3) * t250;
t234 = mrSges(6,1) * t287 - mrSges(6,3) * t250;
t324 = m(7) * t176 + t206 * mrSges(7,3) + t249 * t225;
t165 = m(6) * t179 + t206 * mrSges(6,3) + t249 * t226 + (-t233 - t234) * t287 + (-mrSges(6,2) - mrSges(7,2)) * t269 + t324;
t160 = t301 * t162 + t297 * t165;
t240 = Ifges(5,4) * t275 + Ifges(5,2) * t274 + Ifges(5,6) * t288;
t241 = Ifges(5,1) * t275 + Ifges(5,4) * t274 + Ifges(5,5) * t288;
t331 = mrSges(5,1) * t189 - mrSges(5,2) * t190 + Ifges(5,5) * t246 + Ifges(5,6) * t245 + Ifges(5,3) * t273 + pkin(4) * t160 + t275 * t240 - t274 * t241 + t332;
t228 = -t299 * t238 + t294 * t303;
t223 = -qJDD(3) * pkin(3) - pkin(8) * t305 + t279 * t328 - t228;
t188 = -pkin(4) * t245 - pkin(9) * t272 + t275 * t256 + t223;
t212 = Ifges(7,5) * t250 + Ifges(7,6) * t249 + Ifges(7,3) * t287;
t213 = Ifges(6,5) * t250 + Ifges(6,6) * t249 + Ifges(6,3) * t287;
t182 = -pkin(5) * t206 - qJ(6) * t248 + t232 * t250 + qJDD(6) + t188;
t316 = -mrSges(7,1) * t182 + mrSges(7,3) * t176 + Ifges(7,4) * t207 + Ifges(7,2) * t206 + Ifges(7,6) * t269 + t287 * t216;
t318 = m(7) * t182 - t206 * mrSges(7,1) + t207 * mrSges(7,2) - t249 * t230 + t250 * t233;
t155 = Ifges(6,4) * t207 + Ifges(6,2) * t206 + Ifges(6,6) * t269 + t287 * t217 - mrSges(6,1) * t188 + mrSges(6,3) * t179 - pkin(5) * t318 + qJ(6) * (-t269 * mrSges(7,2) - t287 * t233 + t324) + (-t213 - t212) * t250 + t316;
t314 = mrSges(7,2) * t182 - mrSges(7,3) * t173 + Ifges(7,1) * t207 + Ifges(7,4) * t206 + Ifges(7,5) * t269 + t249 * t212;
t159 = mrSges(6,2) * t188 - mrSges(6,3) * t178 + Ifges(6,1) * t207 + Ifges(6,4) * t206 + Ifges(6,5) * t269 - qJ(6) * t170 + t249 * t213 + (-t214 - t215) * t287 + t314;
t239 = Ifges(5,5) * t275 + Ifges(5,6) * t274 + Ifges(5,3) * t288;
t311 = m(6) * t188 - t206 * mrSges(6,1) + t207 * mrSges(6,2) - t249 * t231 + t250 * t234 + t318;
t319 = -t162 * t297 + t301 * t165;
t140 = -mrSges(5,1) * t223 + mrSges(5,3) * t190 + Ifges(5,4) * t246 + Ifges(5,2) * t245 + Ifges(5,6) * t273 - pkin(4) * t311 + pkin(9) * t319 + t301 * t155 + t297 * t159 - t275 * t239 + t288 * t241;
t141 = mrSges(5,2) * t223 - mrSges(5,3) * t189 + Ifges(5,1) * t246 + Ifges(5,4) * t245 + Ifges(5,5) * t273 - pkin(9) * t160 - t155 * t297 + t159 * t301 + t239 * t274 - t240 * t288;
t253 = -mrSges(5,1) * t274 + mrSges(5,2) * t275;
t254 = -mrSges(5,2) * t288 + mrSges(5,3) * t274;
t157 = m(5) * t189 + mrSges(5,1) * t273 - mrSges(5,3) * t246 - t253 * t275 + t254 * t288 + t160;
t255 = mrSges(5,1) * t288 - mrSges(5,3) * t275;
t158 = m(5) * t190 - mrSges(5,2) * t273 + mrSges(5,3) * t245 + t253 * t274 - t255 * t288 + t319;
t154 = -t157 * t298 + t302 * t158;
t168 = -m(5) * t223 + t245 * mrSges(5,1) - t246 * mrSges(5,2) + t274 * t254 - t275 * t255 - t311;
t267 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t299 + Ifges(4,2) * t303) * qJD(1);
t268 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t299 + Ifges(4,4) * t303) * qJD(1);
t330 = mrSges(4,1) * t228 - mrSges(4,2) * t229 + Ifges(4,5) * t280 + Ifges(4,6) * t281 + Ifges(4,3) * qJDD(3) + pkin(3) * t168 + pkin(8) * t154 + t302 * t140 + t298 * t141 + (t267 * t299 - t268 * t303) * qJD(1);
t277 = (-mrSges(4,1) * t303 + mrSges(4,2) * t299) * qJD(1);
t283 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t328;
t152 = m(4) * t229 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t281 - qJD(3) * t283 + t277 * t327 + t154;
t284 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t327;
t167 = m(4) * t228 + qJDD(3) * mrSges(4,1) - t280 * mrSges(4,3) + qJD(3) * t284 - t277 * t328 + t168;
t320 = t303 * t152 - t167 * t299;
t144 = m(3) * t252 - mrSges(3,1) * t306 - qJDD(1) * mrSges(3,2) + t320;
t153 = t157 * t302 + t158 * t298;
t310 = -m(4) * t237 + t281 * mrSges(4,1) - mrSges(4,2) * t280 - t283 * t328 + t284 * t327 - t153;
t148 = m(3) * t251 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t306 + t310;
t137 = t295 * t144 + t296 * t148;
t146 = t299 * t152 + t303 * t167;
t321 = t296 * t144 - t148 * t295;
t266 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t299 + Ifges(4,6) * t303) * qJD(1);
t133 = mrSges(4,2) * t237 - mrSges(4,3) * t228 + Ifges(4,1) * t280 + Ifges(4,4) * t281 + Ifges(4,5) * qJDD(3) - pkin(8) * t153 - qJD(3) * t267 - t140 * t298 + t141 * t302 + t266 * t327;
t139 = -mrSges(4,1) * t237 + mrSges(4,3) * t229 + Ifges(4,4) * t280 + Ifges(4,2) * t281 + Ifges(4,6) * qJDD(3) - pkin(3) * t153 + qJD(3) * t268 - t266 * t328 - t331;
t313 = mrSges(3,1) * t251 - mrSges(3,2) * t252 + Ifges(3,3) * qJDD(1) + pkin(2) * t310 + pkin(7) * t320 + t299 * t133 + t303 * t139;
t312 = mrSges(2,1) * t285 - mrSges(2,2) * t286 + Ifges(2,3) * qJDD(1) + pkin(1) * t137 + t313;
t135 = m(2) * t286 - mrSges(2,1) * t306 - qJDD(1) * mrSges(2,2) + t321;
t134 = m(2) * t285 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t306 + t137;
t131 = -mrSges(3,1) * t294 + mrSges(3,3) * t252 + t306 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t146 - t330;
t130 = mrSges(3,2) * t294 - mrSges(3,3) * t251 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t306 - pkin(7) * t146 + t133 * t303 - t139 * t299;
t129 = -mrSges(2,2) * g(3) - mrSges(2,3) * t285 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t306 - qJ(2) * t137 + t130 * t296 - t131 * t295;
t128 = Ifges(2,6) * qJDD(1) + t306 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t286 + t295 * t130 + t296 * t131 - pkin(1) * (m(3) * t294 + t146) + qJ(2) * t321;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t129 - t300 * t128 - pkin(6) * (t134 * t304 + t135 * t300), t129, t130, t133, t141, t159, -t214 * t287 + t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t300 * t129 + t304 * t128 + pkin(6) * (-t134 * t300 + t135 * t304), t128, t131, t139, t140, t155, -t250 * t212 + t316; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t312, t312, t313, t330, t331, t332, -t249 * t216 - t315;];
m_new  = t1;
