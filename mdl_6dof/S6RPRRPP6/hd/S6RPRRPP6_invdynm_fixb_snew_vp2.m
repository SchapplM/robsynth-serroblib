% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-05-05 21:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:41:27
% EndTime: 2019-05-05 21:41:41
% DurationCPUTime: 6.11s
% Computational Cost: add. (111531->337), mult. (222316->399), div. (0->0), fcn. (140512->8), ass. (0->128)
t293 = sin(qJ(1));
t296 = cos(qJ(1));
t278 = -t296 * g(1) - t293 * g(2);
t335 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t278;
t334 = -2 * qJD(5);
t333 = (-pkin(1) - pkin(7));
t332 = mrSges(2,1) - mrSges(3,2);
t331 = -mrSges(6,3) - mrSges(7,2);
t330 = -Ifges(3,4) + Ifges(2,5);
t329 = (Ifges(3,5) - Ifges(2,6));
t328 = cos(pkin(9));
t298 = qJD(1) ^ 2;
t244 = (t298 * t333) - t335;
t292 = sin(qJ(3));
t295 = cos(qJ(3));
t323 = qJD(1) * qJD(3);
t319 = t295 * t323;
t272 = -qJDD(1) * t292 - t319;
t320 = t292 * t323;
t273 = qJDD(1) * t295 - t320;
t213 = (-t273 + t320) * pkin(8) + (-t272 + t319) * pkin(3) + t244;
t277 = t293 * g(1) - t296 * g(2);
t312 = -t298 * qJ(2) + qJDD(2) - t277;
t245 = qJDD(1) * t333 + t312;
t238 = -g(3) * t295 + t292 * t245;
t271 = (pkin(3) * t292 - pkin(8) * t295) * qJD(1);
t297 = qJD(3) ^ 2;
t325 = qJD(1) * t292;
t218 = -pkin(3) * t297 + qJDD(3) * pkin(8) - t271 * t325 + t238;
t291 = sin(qJ(4));
t294 = cos(qJ(4));
t183 = t294 * t213 - t291 * t218;
t324 = qJD(1) * t295;
t268 = qJD(3) * t294 - t291 * t324;
t232 = qJD(4) * t268 + qJDD(3) * t291 + t273 * t294;
t267 = qJDD(4) - t272;
t269 = qJD(3) * t291 + t294 * t324;
t280 = qJD(4) + t325;
t179 = (t268 * t280 - t232) * qJ(5) + (t268 * t269 + t267) * pkin(4) + t183;
t184 = t291 * t213 + t294 * t218;
t231 = -qJD(4) * t269 + qJDD(3) * t294 - t273 * t291;
t240 = pkin(4) * t280 - qJ(5) * t269;
t266 = t268 ^ 2;
t181 = -pkin(4) * t266 + qJ(5) * t231 - t240 * t280 + t184;
t290 = sin(pkin(9));
t233 = -t268 * t328 + t290 * t269;
t175 = t290 * t179 + t328 * t181 + t233 * t334;
t203 = -t231 * t328 + t290 * t232;
t234 = t290 * t268 + t269 * t328;
t221 = mrSges(6,1) * t280 - mrSges(6,3) * t234;
t208 = pkin(5) * t233 - qJ(6) * t234;
t279 = t280 ^ 2;
t170 = -pkin(5) * t279 + qJ(6) * t267 + 0.2e1 * qJD(6) * t280 - t208 * t233 + t175;
t222 = -mrSges(7,1) * t280 + mrSges(7,2) * t234;
t321 = m(7) * t170 + t267 * mrSges(7,3) + t280 * t222;
t209 = mrSges(7,1) * t233 - mrSges(7,3) * t234;
t326 = -mrSges(6,1) * t233 - mrSges(6,2) * t234 - t209;
t157 = m(6) * t175 - t267 * mrSges(6,2) + t203 * t331 - t280 * t221 + t233 * t326 + t321;
t313 = t179 * t328 - t290 * t181;
t174 = t234 * t334 + t313;
t204 = t290 * t231 + t232 * t328;
t219 = -mrSges(6,2) * t280 - mrSges(6,3) * t233;
t172 = -t267 * pkin(5) - t279 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t208) * t234 - t313;
t220 = -mrSges(7,2) * t233 + mrSges(7,3) * t280;
t316 = -m(7) * t172 + t267 * mrSges(7,1) + t280 * t220;
t159 = m(6) * t174 + t267 * mrSges(6,1) + t204 * t331 + t280 * t219 + t234 * t326 + t316;
t152 = t290 * t157 + t328 * t159;
t236 = -mrSges(5,1) * t268 + mrSges(5,2) * t269;
t239 = -mrSges(5,2) * t280 + mrSges(5,3) * t268;
t150 = m(5) * t183 + mrSges(5,1) * t267 - mrSges(5,3) * t232 - t236 * t269 + t239 * t280 + t152;
t241 = mrSges(5,1) * t280 - mrSges(5,3) * t269;
t317 = t328 * t157 - t159 * t290;
t151 = m(5) * t184 - mrSges(5,2) * t267 + mrSges(5,3) * t231 + t236 * t268 - t241 * t280 + t317;
t145 = t294 * t150 + t291 * t151;
t199 = Ifges(7,4) * t234 + Ifges(7,2) * t280 + Ifges(7,6) * t233;
t327 = -Ifges(6,5) * t234 + Ifges(6,6) * t233 - Ifges(6,3) * t280 - t199;
t270 = (mrSges(4,1) * t292 + mrSges(4,2) * t295) * qJD(1);
t276 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t324;
t318 = -t150 * t291 + t294 * t151;
t143 = m(4) * t238 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t272 - qJD(3) * t276 - t270 * t325 + t318;
t237 = t292 * g(3) + t295 * t245;
t275 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t325;
t217 = -qJDD(3) * pkin(3) - t297 * pkin(8) + t271 * t324 - t237;
t182 = -t231 * pkin(4) - t266 * qJ(5) + t269 * t240 + qJDD(5) + t217;
t177 = -0.2e1 * qJD(6) * t234 + (t233 * t280 - t204) * qJ(6) + (t234 * t280 + t203) * pkin(5) + t182;
t167 = m(7) * t177 + t203 * mrSges(7,1) - t204 * mrSges(7,3) + t233 * t220 - t234 * t222;
t304 = m(6) * t182 + t203 * mrSges(6,1) + t204 * mrSges(6,2) + t233 * t219 + t234 * t221 + t167;
t300 = -m(5) * t217 + t231 * mrSges(5,1) - t232 * mrSges(5,2) + t268 * t239 - t269 * t241 - t304;
t160 = m(4) * t237 + qJDD(3) * mrSges(4,1) - t273 * mrSges(4,3) + qJD(3) * t275 - t270 * t324 + t300;
t138 = t295 * t143 - t160 * t292;
t315 = -mrSges(7,1) * t177 + mrSges(7,2) * t170;
t137 = t292 * t143 + t295 * t160;
t197 = Ifges(7,5) * t234 + Ifges(7,6) * t280 + Ifges(7,3) * t233;
t311 = mrSges(7,2) * t172 - mrSges(7,3) * t177 + Ifges(7,1) * t204 + Ifges(7,4) * t267 + Ifges(7,5) * t203 + t280 * t197;
t250 = -qJDD(1) * pkin(1) + t312;
t310 = -m(3) * t250 + (t298 * mrSges(3,3)) - t137;
t141 = -m(4) * t244 + mrSges(4,1) * t272 - t273 * mrSges(4,2) - t275 * t325 - t276 * t324 - t145;
t201 = Ifges(7,1) * t234 + Ifges(7,4) * t280 + Ifges(7,5) * t233;
t309 = mrSges(7,1) * t172 - mrSges(7,3) * t170 - Ifges(7,4) * t204 - Ifges(7,2) * t267 - Ifges(7,6) * t203 + t234 * t197 - t233 * t201;
t202 = Ifges(6,1) * t234 - Ifges(6,4) * t233 + Ifges(6,5) * t280;
t153 = -mrSges(6,1) * t182 + mrSges(6,3) * t175 - pkin(5) * t167 + (t201 + t202) * t280 + (Ifges(6,6) - Ifges(7,6)) * t267 + t327 * t234 + (Ifges(6,4) - Ifges(7,5)) * t204 + (-Ifges(6,2) - Ifges(7,3)) * t203 + t315;
t200 = Ifges(6,4) * t234 - Ifges(6,2) * t233 + Ifges(6,6) * t280;
t154 = mrSges(6,2) * t182 - mrSges(6,3) * t174 + Ifges(6,1) * t204 - Ifges(6,4) * t203 + Ifges(6,5) * t267 - qJ(6) * t167 - t280 * t200 + t233 * t327 + t311;
t224 = Ifges(5,5) * t269 + Ifges(5,6) * t268 + Ifges(5,3) * t280;
t226 = Ifges(5,1) * t269 + Ifges(5,4) * t268 + Ifges(5,5) * t280;
t131 = -mrSges(5,1) * t217 + mrSges(5,3) * t184 + Ifges(5,4) * t232 + Ifges(5,2) * t231 + Ifges(5,6) * t267 - pkin(4) * t304 + qJ(5) * t317 + t153 * t328 + t290 * t154 - t269 * t224 + t280 * t226;
t225 = Ifges(5,4) * t269 + Ifges(5,2) * t268 + Ifges(5,6) * t280;
t133 = mrSges(5,2) * t217 - mrSges(5,3) * t183 + Ifges(5,1) * t232 + Ifges(5,4) * t231 + Ifges(5,5) * t267 - qJ(5) * t152 - t290 * t153 + t154 * t328 + t268 * t224 - t280 * t225;
t253 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t295 - Ifges(4,6) * t292) * qJD(1);
t254 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t295 - Ifges(4,2) * t292) * qJD(1);
t128 = mrSges(4,2) * t244 - mrSges(4,3) * t237 + Ifges(4,1) * t273 + Ifges(4,4) * t272 + Ifges(4,5) * qJDD(3) - pkin(8) * t145 - qJD(3) * t254 - t131 * t291 + t133 * t294 - t253 * t325;
t255 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t295 - Ifges(4,4) * t292) * qJD(1);
t301 = mrSges(6,2) * t175 - t233 * t202 - qJ(6) * (-t203 * mrSges(7,2) - t233 * t209 + t321) - pkin(5) * (-t204 * mrSges(7,2) - t234 * t209 + t316) - mrSges(6,1) * t174 - t234 * t200 + Ifges(6,6) * t203 - Ifges(6,5) * t204 - Ifges(6,3) * t267 + t309;
t299 = mrSges(5,1) * t183 - mrSges(5,2) * t184 + Ifges(5,5) * t232 + Ifges(5,6) * t231 + Ifges(5,3) * t267 + pkin(4) * t152 + t269 * t225 - t268 * t226 - t301;
t129 = -mrSges(4,1) * t244 + mrSges(4,3) * t238 + Ifges(4,4) * t273 + Ifges(4,2) * t272 + Ifges(4,6) * qJDD(3) - pkin(3) * t145 + qJD(3) * t255 - t253 * t324 - t299;
t248 = t298 * pkin(1) + t335;
t308 = mrSges(3,2) * t250 - mrSges(3,3) * t248 + Ifges(3,1) * qJDD(1) - pkin(7) * t137 + t295 * t128 - t129 * t292;
t307 = -mrSges(3,1) * t248 - pkin(2) * t141 - pkin(7) * t138 - t292 * t128 - t295 * t129;
t306 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t273 + Ifges(4,6) * t272 + Ifges(4,3) * qJDD(3) + pkin(3) * t300 + pkin(8) * t318 + t294 * t131 + t291 * t133 + t254 * t324 + t255 * t325;
t305 = -m(3) * t248 + t298 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t141;
t303 = -mrSges(2,2) * t278 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t310) + qJ(2) * t305 + mrSges(2,1) * t277 + Ifges(2,3) * qJDD(1) + t308;
t302 = mrSges(3,1) * t250 + pkin(2) * t137 + t306;
t139 = m(2) * t278 - mrSges(2,1) * t298 - qJDD(1) * mrSges(2,2) + t305;
t136 = -m(3) * g(3) + t138;
t134 = m(2) * t277 - t298 * mrSges(2,2) + qJDD(1) * t332 + t310;
t126 = -mrSges(2,3) * t277 - qJ(2) * t136 + t302 + (t329 * t298) + t330 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t125 = mrSges(2,3) * t278 - pkin(1) * t136 + g(3) * t332 - qJDD(1) * t329 + t298 * t330 + t307;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t296 * t126 - t293 * t125 - pkin(6) * (t134 * t296 + t139 * t293), t126, t308, t128, t133, t154, -t199 * t233 + t311; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t293 * t126 + t296 * t125 + pkin(6) * (-t134 * t293 + t139 * t296), t125, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t298 * Ifges(3,5)) - t302, t129, t131, t153, -t309; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t303, t303, mrSges(3,2) * g(3) + t298 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t307, t306, t299, -t301, Ifges(7,5) * t204 + Ifges(7,6) * t267 + Ifges(7,3) * t203 + t234 * t199 - t280 * t201 - t315;];
m_new  = t1;
