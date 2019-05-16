% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 02:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:31:40
% EndTime: 2019-05-06 02:32:10
% DurationCPUTime: 21.78s
% Computational Cost: add. (417690->342), mult. (871682->431), div. (0->0), fcn. (617753->12), ass. (0->139)
t304 = sin(qJ(1));
t309 = cos(qJ(1));
t282 = t304 * g(1) - g(2) * t309;
t273 = qJDD(1) * pkin(1) + t282;
t283 = -g(1) * t309 - g(2) * t304;
t310 = qJD(1) ^ 2;
t275 = -pkin(1) * t310 + t283;
t298 = sin(pkin(11));
t299 = cos(pkin(11));
t257 = t298 * t273 + t299 * t275;
t251 = -pkin(2) * t310 + qJDD(1) * pkin(7) + t257;
t297 = -g(3) + qJDD(2);
t303 = sin(qJ(3));
t308 = cos(qJ(3));
t237 = -t303 * t251 + t308 * t297;
t329 = qJD(1) * qJD(3);
t328 = t308 * t329;
t276 = qJDD(1) * t303 + t328;
t226 = (-t276 + t328) * pkin(8) + (t303 * t308 * t310 + qJDD(3)) * pkin(3) + t237;
t238 = t308 * t251 + t303 * t297;
t277 = qJDD(1) * t308 - t303 * t329;
t331 = qJD(1) * t303;
t281 = qJD(3) * pkin(3) - pkin(8) * t331;
t296 = t308 ^ 2;
t229 = -pkin(3) * t296 * t310 + pkin(8) * t277 - qJD(3) * t281 + t238;
t302 = sin(qJ(4));
t307 = cos(qJ(4));
t201 = t307 * t226 - t302 * t229;
t268 = (-t302 * t303 + t307 * t308) * qJD(1);
t240 = qJD(4) * t268 + t276 * t307 + t277 * t302;
t269 = (t302 * t308 + t303 * t307) * qJD(1);
t292 = qJDD(3) + qJDD(4);
t293 = qJD(3) + qJD(4);
t192 = (t268 * t293 - t240) * pkin(9) + (t268 * t269 + t292) * pkin(4) + t201;
t202 = t302 * t226 + t307 * t229;
t239 = -qJD(4) * t269 - t276 * t302 + t277 * t307;
t260 = pkin(4) * t293 - pkin(9) * t269;
t261 = t268 ^ 2;
t194 = -pkin(4) * t261 + pkin(9) * t239 - t260 * t293 + t202;
t301 = sin(qJ(5));
t306 = cos(qJ(5));
t190 = t301 * t192 + t306 * t194;
t253 = t268 * t301 + t269 * t306;
t212 = -qJD(5) * t253 + t239 * t306 - t240 * t301;
t252 = t268 * t306 - t269 * t301;
t227 = -mrSges(6,1) * t252 + mrSges(6,2) * t253;
t290 = qJD(5) + t293;
t242 = mrSges(6,1) * t290 - mrSges(6,3) * t253;
t289 = qJDD(5) + t292;
t228 = -pkin(5) * t252 - pkin(10) * t253;
t288 = t290 ^ 2;
t186 = -pkin(5) * t288 + pkin(10) * t289 + t228 * t252 + t190;
t256 = t299 * t273 - t298 * t275;
t321 = -qJDD(1) * pkin(2) - t256;
t230 = -t277 * pkin(3) + t281 * t331 + (-pkin(8) * t296 - pkin(7)) * t310 + t321;
t199 = -t239 * pkin(4) - t261 * pkin(9) + t269 * t260 + t230;
t213 = qJD(5) * t252 + t239 * t301 + t240 * t306;
t187 = (-t252 * t290 - t213) * pkin(10) + (t253 * t290 - t212) * pkin(5) + t199;
t300 = sin(qJ(6));
t305 = cos(qJ(6));
t183 = -t186 * t300 + t187 * t305;
t232 = -t253 * t300 + t290 * t305;
t197 = qJD(6) * t232 + t213 * t305 + t289 * t300;
t210 = qJDD(6) - t212;
t233 = t253 * t305 + t290 * t300;
t215 = -mrSges(7,1) * t232 + mrSges(7,2) * t233;
t247 = qJD(6) - t252;
t216 = -mrSges(7,2) * t247 + mrSges(7,3) * t232;
t179 = m(7) * t183 + mrSges(7,1) * t210 - mrSges(7,3) * t197 - t215 * t233 + t216 * t247;
t184 = t186 * t305 + t187 * t300;
t196 = -qJD(6) * t233 - t213 * t300 + t289 * t305;
t217 = mrSges(7,1) * t247 - mrSges(7,3) * t233;
t180 = m(7) * t184 - mrSges(7,2) * t210 + mrSges(7,3) * t196 + t215 * t232 - t217 * t247;
t323 = -t179 * t300 + t305 * t180;
t166 = m(6) * t190 - mrSges(6,2) * t289 + mrSges(6,3) * t212 + t227 * t252 - t242 * t290 + t323;
t189 = t192 * t306 - t194 * t301;
t241 = -mrSges(6,2) * t290 + mrSges(6,3) * t252;
t185 = -pkin(5) * t289 - pkin(10) * t288 + t228 * t253 - t189;
t318 = -m(7) * t185 + t196 * mrSges(7,1) - mrSges(7,2) * t197 + t232 * t216 - t217 * t233;
t175 = m(6) * t189 + mrSges(6,1) * t289 - mrSges(6,3) * t213 - t227 * t253 + t241 * t290 + t318;
t160 = t301 * t166 + t306 * t175;
t254 = -mrSges(5,1) * t268 + mrSges(5,2) * t269;
t258 = -mrSges(5,2) * t293 + mrSges(5,3) * t268;
t157 = m(5) * t201 + mrSges(5,1) * t292 - mrSges(5,3) * t240 - t254 * t269 + t258 * t293 + t160;
t259 = mrSges(5,1) * t293 - mrSges(5,3) * t269;
t324 = t306 * t166 - t175 * t301;
t158 = m(5) * t202 - mrSges(5,2) * t292 + mrSges(5,3) * t239 + t254 * t268 - t259 * t293 + t324;
t151 = t307 * t157 + t302 * t158;
t266 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t303 + Ifges(4,2) * t308) * qJD(1);
t267 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t303 + Ifges(4,4) * t308) * qJD(1);
t245 = Ifges(5,4) * t269 + Ifges(5,2) * t268 + Ifges(5,6) * t293;
t246 = Ifges(5,1) * t269 + Ifges(5,4) * t268 + Ifges(5,5) * t293;
t203 = Ifges(7,5) * t233 + Ifges(7,6) * t232 + Ifges(7,3) * t247;
t205 = Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t247;
t172 = -mrSges(7,1) * t185 + mrSges(7,3) * t184 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t210 - t203 * t233 + t205 * t247;
t204 = Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t247;
t173 = mrSges(7,2) * t185 - mrSges(7,3) * t183 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t210 + t203 * t232 - t204 * t247;
t219 = Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t290;
t220 = Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t290;
t316 = -mrSges(6,1) * t189 + mrSges(6,2) * t190 - Ifges(6,5) * t213 - Ifges(6,6) * t212 - Ifges(6,3) * t289 - pkin(5) * t318 - pkin(10) * t323 - t305 * t172 - t300 * t173 - t253 * t219 + t252 * t220;
t313 = -mrSges(5,1) * t201 + mrSges(5,2) * t202 - Ifges(5,5) * t240 - Ifges(5,6) * t239 - Ifges(5,3) * t292 - pkin(4) * t160 - t269 * t245 + t268 * t246 + t316;
t332 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t276 + Ifges(4,6) * t277 + Ifges(4,3) * qJDD(3) + pkin(3) * t151 + (t266 * t303 - t267 * t308) * qJD(1) - t313;
t274 = (-mrSges(4,1) * t308 + mrSges(4,2) * t303) * qJD(1);
t330 = qJD(1) * t308;
t280 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t330;
t149 = m(4) * t237 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t276 + qJD(3) * t280 - t274 * t331 + t151;
t279 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t331;
t325 = -t157 * t302 + t307 * t158;
t150 = m(4) * t238 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t277 - qJD(3) * t279 + t274 * t330 + t325;
t326 = -t149 * t303 + t308 * t150;
t142 = m(3) * t257 - mrSges(3,1) * t310 - qJDD(1) * mrSges(3,2) + t326;
t250 = -pkin(7) * t310 + t321;
t168 = t305 * t179 + t300 * t180;
t320 = m(6) * t199 - t212 * mrSges(6,1) + t213 * mrSges(6,2) - t252 * t241 + t253 * t242 + t168;
t315 = m(5) * t230 - t239 * mrSges(5,1) + mrSges(5,2) * t240 - t268 * t258 + t259 * t269 + t320;
t312 = -m(4) * t250 + t277 * mrSges(4,1) - mrSges(4,2) * t276 - t279 * t331 + t280 * t330 - t315;
t162 = m(3) * t256 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t310 + t312;
t138 = t298 * t142 + t299 * t162;
t144 = t308 * t149 + t303 * t150;
t327 = t299 * t142 - t162 * t298;
t218 = Ifges(6,5) * t253 + Ifges(6,6) * t252 + Ifges(6,3) * t290;
t152 = mrSges(6,2) * t199 - mrSges(6,3) * t189 + Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t289 - pkin(10) * t168 - t172 * t300 + t173 * t305 + t218 * t252 - t219 * t290;
t314 = mrSges(7,1) * t183 - mrSges(7,2) * t184 + Ifges(7,5) * t197 + Ifges(7,6) * t196 + Ifges(7,3) * t210 + t204 * t233 - t205 * t232;
t153 = -mrSges(6,1) * t199 + mrSges(6,3) * t190 + Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t289 - pkin(5) * t168 - t218 * t253 + t220 * t290 - t314;
t244 = Ifges(5,5) * t269 + Ifges(5,6) * t268 + Ifges(5,3) * t293;
t139 = -mrSges(5,1) * t230 + mrSges(5,3) * t202 + Ifges(5,4) * t240 + Ifges(5,2) * t239 + Ifges(5,6) * t292 - pkin(4) * t320 + pkin(9) * t324 + t301 * t152 + t306 * t153 - t269 * t244 + t293 * t246;
t145 = mrSges(5,2) * t230 - mrSges(5,3) * t201 + Ifges(5,1) * t240 + Ifges(5,4) * t239 + Ifges(5,5) * t292 - pkin(9) * t160 + t152 * t306 - t153 * t301 + t244 * t268 - t245 * t293;
t265 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t303 + Ifges(4,6) * t308) * qJD(1);
t131 = -mrSges(4,1) * t250 + mrSges(4,3) * t238 + Ifges(4,4) * t276 + Ifges(4,2) * t277 + Ifges(4,6) * qJDD(3) - pkin(3) * t315 + pkin(8) * t325 + qJD(3) * t267 + t307 * t139 + t302 * t145 - t265 * t331;
t133 = mrSges(4,2) * t250 - mrSges(4,3) * t237 + Ifges(4,1) * t276 + Ifges(4,4) * t277 + Ifges(4,5) * qJDD(3) - pkin(8) * t151 - qJD(3) * t266 - t139 * t302 + t145 * t307 + t265 * t330;
t319 = mrSges(3,1) * t256 - mrSges(3,2) * t257 + Ifges(3,3) * qJDD(1) + pkin(2) * t312 + pkin(7) * t326 + t308 * t131 + t303 * t133;
t317 = mrSges(2,1) * t282 - mrSges(2,2) * t283 + Ifges(2,3) * qJDD(1) + pkin(1) * t138 + t319;
t136 = m(2) * t283 - mrSges(2,1) * t310 - qJDD(1) * mrSges(2,2) + t327;
t135 = m(2) * t282 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t310 + t138;
t134 = -mrSges(3,1) * t297 + mrSges(3,3) * t257 + t310 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t144 - t332;
t129 = mrSges(3,2) * t297 - mrSges(3,3) * t256 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t310 - pkin(7) * t144 - t131 * t303 + t133 * t308;
t128 = -mrSges(2,2) * g(3) - mrSges(2,3) * t282 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t310 - qJ(2) * t138 + t129 * t299 - t134 * t298;
t127 = Ifges(2,6) * qJDD(1) + t310 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t283 + t298 * t129 + t299 * t134 - pkin(1) * (m(3) * t297 + t144) + qJ(2) * t327;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t128 - t304 * t127 - pkin(6) * (t135 * t309 + t136 * t304), t128, t129, t133, t145, t152, t173; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t128 + t309 * t127 + pkin(6) * (-t135 * t304 + t136 * t309), t127, t134, t131, t139, t153, t172; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t317, t317, t319, t332, -t313, -t316, t314;];
m_new  = t1;
