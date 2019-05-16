% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:14:24
% EndTime: 2019-05-06 04:14:52
% DurationCPUTime: 13.95s
% Computational Cost: add. (267009->343), mult. (522646->419), div. (0->0), fcn. (358604->10), ass. (0->138)
t308 = sin(qJ(1));
t313 = cos(qJ(1));
t286 = -t313 * g(1) - t308 * g(2);
t329 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t286;
t342 = -pkin(1) - pkin(7);
t341 = mrSges(2,1) - mrSges(3,2);
t340 = Ifges(2,5) - Ifges(3,4);
t339 = (-Ifges(2,6) + Ifges(3,5));
t285 = t308 * g(1) - t313 * g(2);
t314 = qJD(1) ^ 2;
t328 = -t314 * qJ(2) + qJDD(2) - t285;
t257 = t342 * qJDD(1) + t328;
t307 = sin(qJ(3));
t312 = cos(qJ(3));
t246 = t307 * g(3) + t312 * t257;
t336 = qJD(1) * qJD(3);
t334 = t307 * t336;
t280 = qJDD(1) * t312 - t334;
t219 = (-t280 - t334) * pkin(8) + (-t307 * t312 * t314 + qJDD(3)) * pkin(3) + t246;
t247 = -g(3) * t312 + t307 * t257;
t279 = -qJDD(1) * t307 - t312 * t336;
t337 = qJD(1) * t312;
t284 = qJD(3) * pkin(3) - pkin(8) * t337;
t301 = t307 ^ 2;
t222 = -pkin(3) * t301 * t314 + pkin(8) * t279 - qJD(3) * t284 + t247;
t306 = sin(qJ(4));
t311 = cos(qJ(4));
t205 = t306 * t219 + t311 * t222;
t270 = (-t306 * t307 + t311 * t312) * qJD(1);
t234 = -t270 * qJD(4) + t279 * t311 - t306 * t280;
t338 = qJD(1) * t307;
t269 = -t306 * t337 - t311 * t338;
t243 = -mrSges(5,1) * t269 + mrSges(5,2) * t270;
t295 = qJD(3) + qJD(4);
t255 = mrSges(5,1) * t295 - mrSges(5,3) * t270;
t294 = qJDD(3) + qJDD(4);
t226 = -t279 * pkin(3) + t284 * t337 + (-pkin(8) * t301 + t342) * t314 + t329;
t235 = qJD(4) * t269 + t279 * t306 + t280 * t311;
t195 = (-t269 * t295 - t235) * pkin(9) + (t270 * t295 - t234) * pkin(4) + t226;
t244 = -pkin(4) * t269 - pkin(9) * t270;
t293 = t295 ^ 2;
t198 = -pkin(4) * t293 + pkin(9) * t294 + t244 * t269 + t205;
t305 = sin(qJ(5));
t310 = cos(qJ(5));
t184 = t310 * t195 - t305 * t198;
t249 = -t270 * t305 + t295 * t310;
t209 = qJD(5) * t249 + t235 * t310 + t294 * t305;
t233 = qJDD(5) - t234;
t250 = t270 * t310 + t295 * t305;
t265 = qJD(5) - t269;
t182 = (t249 * t265 - t209) * pkin(10) + (t249 * t250 + t233) * pkin(5) + t184;
t185 = t305 * t195 + t310 * t198;
t208 = -qJD(5) * t250 - t235 * t305 + t294 * t310;
t238 = pkin(5) * t265 - pkin(10) * t250;
t248 = t249 ^ 2;
t183 = -pkin(5) * t248 + pkin(10) * t208 - t238 * t265 + t185;
t304 = sin(qJ(6));
t309 = cos(qJ(6));
t180 = t182 * t309 - t183 * t304;
t220 = t249 * t309 - t250 * t304;
t191 = qJD(6) * t220 + t208 * t304 + t209 * t309;
t221 = t249 * t304 + t250 * t309;
t206 = -mrSges(7,1) * t220 + mrSges(7,2) * t221;
t262 = qJD(6) + t265;
t210 = -mrSges(7,2) * t262 + mrSges(7,3) * t220;
t228 = qJDD(6) + t233;
t175 = m(7) * t180 + mrSges(7,1) * t228 - t191 * mrSges(7,3) - t206 * t221 + t210 * t262;
t181 = t182 * t304 + t183 * t309;
t190 = -qJD(6) * t221 + t208 * t309 - t209 * t304;
t211 = mrSges(7,1) * t262 - mrSges(7,3) * t221;
t176 = m(7) * t181 - mrSges(7,2) * t228 + t190 * mrSges(7,3) + t206 * t220 - t211 * t262;
t167 = t309 * t175 + t304 * t176;
t223 = -mrSges(6,1) * t249 + mrSges(6,2) * t250;
t236 = -mrSges(6,2) * t265 + mrSges(6,3) * t249;
t165 = m(6) * t184 + mrSges(6,1) * t233 - mrSges(6,3) * t209 - t223 * t250 + t236 * t265 + t167;
t237 = mrSges(6,1) * t265 - mrSges(6,3) * t250;
t331 = -t175 * t304 + t309 * t176;
t166 = m(6) * t185 - mrSges(6,2) * t233 + mrSges(6,3) * t208 + t223 * t249 - t237 * t265 + t331;
t332 = -t165 * t305 + t310 * t166;
t158 = m(5) * t205 - mrSges(5,2) * t294 + mrSges(5,3) * t234 + t243 * t269 - t255 * t295 + t332;
t204 = t219 * t311 - t306 * t222;
t254 = -mrSges(5,2) * t295 + mrSges(5,3) * t269;
t197 = -pkin(4) * t294 - pkin(9) * t293 + t270 * t244 - t204;
t186 = -pkin(5) * t208 - pkin(10) * t248 + t238 * t250 + t197;
t326 = m(7) * t186 - t190 * mrSges(7,1) + t191 * mrSges(7,2) - t220 * t210 + t211 * t221;
t319 = -m(6) * t197 + t208 * mrSges(6,1) - mrSges(6,2) * t209 + t249 * t236 - t237 * t250 - t326;
t171 = m(5) * t204 + mrSges(5,1) * t294 - mrSges(5,3) * t235 - t243 * t270 + t254 * t295 + t319;
t150 = t306 * t158 + t311 * t171;
t160 = t310 * t165 + t305 * t166;
t278 = (mrSges(4,1) * t307 + mrSges(4,2) * t312) * qJD(1);
t282 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t338;
t145 = m(4) * t246 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t280 + qJD(3) * t282 - t278 * t337 + t150;
t283 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t337;
t333 = t311 * t158 - t171 * t306;
t146 = m(4) * t247 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t279 - qJD(3) * t283 - t278 * t338 + t333;
t142 = -t145 * t307 + t312 * t146;
t141 = t312 * t145 + t307 * t146;
t264 = -qJDD(1) * pkin(1) + t328;
t327 = -m(3) * t264 + (t314 * mrSges(3,3)) - t141;
t325 = m(5) * t226 - mrSges(5,1) * t234 + t235 * mrSges(5,2) - t254 * t269 + t270 * t255 + t160;
t200 = Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t262;
t201 = Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t262;
t324 = -mrSges(7,1) * t180 + mrSges(7,2) * t181 - Ifges(7,5) * t191 - Ifges(7,6) * t190 - Ifges(7,3) * t228 - t221 * t200 + t220 * t201;
t199 = Ifges(7,5) * t221 + Ifges(7,6) * t220 + Ifges(7,3) * t262;
t168 = -mrSges(7,1) * t186 + mrSges(7,3) * t181 + Ifges(7,4) * t191 + Ifges(7,2) * t190 + Ifges(7,6) * t228 - t199 * t221 + t201 * t262;
t169 = mrSges(7,2) * t186 - mrSges(7,3) * t180 + Ifges(7,1) * t191 + Ifges(7,4) * t190 + Ifges(7,5) * t228 + t199 * t220 - t200 * t262;
t212 = Ifges(6,5) * t250 + Ifges(6,6) * t249 + Ifges(6,3) * t265;
t214 = Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t265;
t148 = -mrSges(6,1) * t197 + mrSges(6,3) * t185 + Ifges(6,4) * t209 + Ifges(6,2) * t208 + Ifges(6,6) * t233 - pkin(5) * t326 + pkin(10) * t331 + t309 * t168 + t304 * t169 - t250 * t212 + t265 * t214;
t213 = Ifges(6,4) * t250 + Ifges(6,2) * t249 + Ifges(6,6) * t265;
t152 = mrSges(6,2) * t197 - mrSges(6,3) * t184 + Ifges(6,1) * t209 + Ifges(6,4) * t208 + Ifges(6,5) * t233 - pkin(10) * t167 - t168 * t304 + t169 * t309 + t212 * t249 - t213 * t265;
t239 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * t295;
t240 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * t295;
t137 = mrSges(5,2) * t226 - mrSges(5,3) * t204 + Ifges(5,1) * t235 + Ifges(5,4) * t234 + Ifges(5,5) * t294 - pkin(9) * t160 - t148 * t305 + t152 * t310 + t239 * t269 - t240 * t295;
t241 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * t295;
t315 = mrSges(6,1) * t184 - mrSges(6,2) * t185 + Ifges(6,5) * t209 + Ifges(6,6) * t208 + Ifges(6,3) * t233 + pkin(5) * t167 + t250 * t213 - t249 * t214 - t324;
t143 = -mrSges(5,1) * t226 + mrSges(5,3) * t205 + Ifges(5,4) * t235 + Ifges(5,2) * t234 + Ifges(5,6) * t294 - pkin(4) * t160 - t270 * t239 + t295 * t241 - t315;
t256 = t342 * t314 + t329;
t266 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t312 - Ifges(4,6) * t307) * qJD(1);
t268 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 - Ifges(4,4) * t307) * qJD(1);
t134 = -mrSges(4,1) * t256 + mrSges(4,3) * t247 + Ifges(4,4) * t280 + Ifges(4,2) * t279 + Ifges(4,6) * qJDD(3) - pkin(3) * t325 + pkin(8) * t333 + qJD(3) * t268 + t306 * t137 + t311 * t143 - t266 * t337;
t267 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 - Ifges(4,2) * t307) * qJD(1);
t136 = mrSges(4,2) * t256 - mrSges(4,3) * t246 + Ifges(4,1) * t280 + Ifges(4,4) * t279 + Ifges(4,5) * qJDD(3) - pkin(8) * t150 - qJD(3) * t267 + t137 * t311 - t143 * t306 - t266 * t338;
t260 = t314 * pkin(1) - t329;
t323 = mrSges(3,2) * t264 - mrSges(3,3) * t260 + Ifges(3,1) * qJDD(1) - pkin(7) * t141 - t134 * t307 + t312 * t136;
t155 = -m(4) * t256 + mrSges(4,1) * t279 - t280 * mrSges(4,2) - t282 * t338 - t283 * t337 - t325;
t322 = -mrSges(3,1) * t260 - pkin(2) * t155 - pkin(7) * t142 - t312 * t134 - t307 * t136;
t321 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t235 - Ifges(5,6) * t234 - Ifges(5,3) * t294 - pkin(4) * t319 - pkin(9) * t332 - t310 * t148 - t305 * t152 - t270 * t240 + t269 * t241;
t318 = -m(3) * t260 + t314 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t155;
t320 = -mrSges(2,2) * t286 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t327) + qJ(2) * t318 + mrSges(2,1) * t285 + Ifges(2,3) * qJDD(1) + t323;
t317 = -mrSges(4,1) * t246 + mrSges(4,2) * t247 - Ifges(4,5) * t280 - Ifges(4,6) * t279 - Ifges(4,3) * qJDD(3) - pkin(3) * t150 - t267 * t337 - t268 * t338 + t321;
t316 = -mrSges(3,1) * t264 - pkin(2) * t141 + t317;
t153 = m(2) * t286 - mrSges(2,1) * t314 - qJDD(1) * mrSges(2,2) + t318;
t140 = -m(3) * g(3) + t142;
t138 = m(2) * t285 - t314 * mrSges(2,2) + t341 * qJDD(1) + t327;
t133 = -t316 + (t339 * t314) + t340 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t285 - qJ(2) * t140;
t132 = mrSges(2,3) * t286 - pkin(1) * t140 + t341 * g(3) - t339 * qJDD(1) + t340 * t314 + t322;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t313 * t133 - t308 * t132 - pkin(6) * (t138 * t313 + t153 * t308), t133, t323, t136, t137, t152, t169; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t133 + t313 * t132 + pkin(6) * (-t138 * t308 + t153 * t313), t132, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t314 * Ifges(3,5)) + t316, t134, t143, t148, t168; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t320, t320, mrSges(3,2) * g(3) + t314 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t322, -t317, -t321, t315, -t324;];
m_new  = t1;
