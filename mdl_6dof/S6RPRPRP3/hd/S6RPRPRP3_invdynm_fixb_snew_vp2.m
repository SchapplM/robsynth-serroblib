% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:37:31
% EndTime: 2019-05-05 17:37:48
% DurationCPUTime: 9.24s
% Computational Cost: add. (166686->337), mult. (342771->410), div. (0->0), fcn. (221567->10), ass. (0->126)
t293 = sin(qJ(1));
t295 = cos(qJ(1));
t277 = g(1) * t293 - g(2) * t295;
t268 = qJDD(1) * pkin(1) + t277;
t278 = -g(1) * t295 - g(2) * t293;
t297 = qJD(1) ^ 2;
t271 = -pkin(1) * t297 + t278;
t288 = sin(pkin(9));
t290 = cos(pkin(9));
t238 = t268 * t288 + t271 * t290;
t228 = -pkin(2) * t297 + qJDD(1) * pkin(7) + t238;
t286 = -g(3) + qJDD(2);
t292 = sin(qJ(3));
t294 = cos(qJ(3));
t219 = -t228 * t292 + t294 * t286;
t269 = (-pkin(3) * t294 - qJ(4) * t292) * qJD(1);
t296 = qJD(3) ^ 2;
t317 = qJD(1) * t292;
t212 = -qJDD(3) * pkin(3) - t296 * qJ(4) + t269 * t317 + qJDD(4) - t219;
t315 = qJD(1) * qJD(3);
t313 = t294 * t315;
t272 = qJDD(1) * t292 + t313;
t287 = sin(pkin(10));
t289 = cos(pkin(10));
t244 = qJDD(3) * t289 - t272 * t287;
t265 = qJD(3) * t287 + t289 * t317;
t316 = t294 * qJD(1);
t246 = -pkin(4) * t316 - pkin(8) * t265;
t264 = qJD(3) * t289 - t287 * t317;
t263 = t264 ^ 2;
t184 = -t244 * pkin(4) - t263 * pkin(8) + t246 * t265 + t212;
t291 = sin(qJ(5));
t321 = cos(qJ(5));
t235 = t264 * t291 + t265 * t321;
t245 = qJDD(3) * t287 + t272 * t289;
t197 = qJD(5) * t235 - t244 * t321 + t245 * t291;
t234 = -t264 * t321 + t265 * t291;
t198 = -qJD(5) * t234 + t244 * t291 + t245 * t321;
t280 = qJD(5) - t316;
t177 = -0.2e1 * qJD(6) * t235 + (t234 * t280 - t198) * qJ(6) + (t235 * t280 + t197) * pkin(5) + t184;
t223 = -mrSges(7,1) * t280 + mrSges(7,2) * t235;
t224 = -mrSges(7,2) * t234 + mrSges(7,3) * t280;
t167 = m(7) * t177 + mrSges(7,1) * t197 - mrSges(7,3) * t198 - t223 * t235 + t224 * t234;
t237 = t290 * t268 - t271 * t288;
t227 = -qJDD(1) * pkin(2) - t297 * pkin(7) - t237;
t282 = t292 * t315;
t273 = qJDD(1) * t294 - t282;
t210 = (-t272 - t313) * qJ(4) + (-t273 + t282) * pkin(3) + t227;
t220 = t228 * t294 + t286 * t292;
t216 = -pkin(3) * t296 + qJDD(3) * qJ(4) + t269 * t316 + t220;
t182 = -0.2e1 * qJD(4) * t265 + t210 * t289 - t287 * t216;
t179 = (-t264 * t316 - t245) * pkin(8) + (t264 * t265 - t273) * pkin(4) + t182;
t183 = 0.2e1 * qJD(4) * t264 + t210 * t287 + t216 * t289;
t181 = -pkin(4) * t263 + pkin(8) * t244 + t246 * t316 + t183;
t175 = t179 * t291 + t181 * t321;
t205 = Ifges(7,1) * t235 + Ifges(7,4) * t280 + Ifges(7,5) * t234;
t206 = Ifges(6,1) * t235 - Ifges(6,4) * t234 + Ifges(6,5) * t280;
t267 = qJDD(5) - t273;
t213 = pkin(5) * t234 - qJ(6) * t235;
t279 = t280 ^ 2;
t170 = -pkin(5) * t279 + qJ(6) * t267 + 0.2e1 * qJD(6) * t280 - t213 * t234 + t175;
t308 = -mrSges(7,1) * t177 + mrSges(7,2) * t170;
t203 = Ifges(7,4) * t235 + Ifges(7,2) * t280 + Ifges(7,6) * t234;
t319 = -Ifges(6,5) * t235 + Ifges(6,6) * t234 - Ifges(6,3) * t280 - t203;
t153 = -mrSges(6,1) * t184 + mrSges(6,3) * t175 - pkin(5) * t167 + (t205 + t206) * t280 + (Ifges(6,6) - Ifges(7,6)) * t267 + t319 * t235 + (Ifges(6,4) - Ifges(7,5)) * t198 + (-Ifges(6,2) - Ifges(7,3)) * t197 + t308;
t174 = t179 * t321 - t181 * t291;
t204 = Ifges(6,4) * t235 - Ifges(6,2) * t234 + Ifges(6,6) * t280;
t172 = -pkin(5) * t267 - qJ(6) * t279 + t213 * t235 + qJDD(6) - t174;
t201 = Ifges(7,5) * t235 + Ifges(7,6) * t280 + Ifges(7,3) * t234;
t306 = mrSges(7,2) * t172 - mrSges(7,3) * t177 + Ifges(7,1) * t198 + Ifges(7,4) * t267 + Ifges(7,5) * t197 + t201 * t280;
t154 = mrSges(6,2) * t184 - mrSges(6,3) * t174 + Ifges(6,1) * t198 - Ifges(6,4) * t197 + Ifges(6,5) * t267 - qJ(6) * t167 - t280 * t204 + t234 * t319 + t306;
t229 = Ifges(5,5) * t265 + Ifges(5,6) * t264 - Ifges(5,3) * t316;
t231 = Ifges(5,1) * t265 + Ifges(5,4) * t264 - Ifges(5,5) * t316;
t221 = -mrSges(6,2) * t280 - mrSges(6,3) * t234;
t222 = mrSges(6,1) * t280 - mrSges(6,3) * t235;
t302 = m(6) * t184 + mrSges(6,1) * t197 + mrSges(6,2) * t198 + t221 * t234 + t222 * t235 + t167;
t314 = m(7) * t170 + mrSges(7,3) * t267 + t223 * t280;
t214 = mrSges(7,1) * t234 - mrSges(7,3) * t235;
t318 = -mrSges(6,1) * t234 - mrSges(6,2) * t235 - t214;
t320 = -mrSges(6,3) - mrSges(7,2);
t157 = m(6) * t175 - t267 * mrSges(6,2) + t197 * t320 - t280 * t222 + t234 * t318 + t314;
t309 = -m(7) * t172 + mrSges(7,1) * t267 + t224 * t280;
t159 = m(6) * t174 + t267 * mrSges(6,1) + t198 * t320 + t280 * t221 + t235 * t318 + t309;
t310 = t157 * t321 - t159 * t291;
t134 = -mrSges(5,1) * t212 + mrSges(5,3) * t183 + Ifges(5,4) * t245 + Ifges(5,2) * t244 - Ifges(5,6) * t273 - pkin(4) * t302 + pkin(8) * t310 + t153 * t321 + t291 * t154 - t265 * t229 - t231 * t316;
t152 = t157 * t291 + t159 * t321;
t230 = Ifges(5,4) * t265 + Ifges(5,2) * t264 - Ifges(5,6) * t316;
t135 = mrSges(5,2) * t212 - mrSges(5,3) * t182 + Ifges(5,1) * t245 + Ifges(5,4) * t244 - Ifges(5,5) * t273 - pkin(8) * t152 - t153 * t291 + t154 * t321 + t229 * t264 + t230 * t316;
t239 = -mrSges(5,1) * t264 + mrSges(5,2) * t265;
t242 = mrSges(5,2) * t316 + mrSges(5,3) * t264;
t150 = m(5) * t182 - mrSges(5,1) * t273 - mrSges(5,3) * t245 - t239 * t265 - t242 * t316 + t152;
t243 = -mrSges(5,1) * t316 - mrSges(5,3) * t265;
t151 = m(5) * t183 + mrSges(5,2) * t273 + mrSges(5,3) * t244 + t239 * t264 + t243 * t316 + t310;
t148 = -t150 * t287 + t151 * t289;
t162 = -m(5) * t212 + mrSges(5,1) * t244 - mrSges(5,2) * t245 + t242 * t264 - t243 * t265 - t302;
t255 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t292 + Ifges(4,2) * t294) * qJD(1);
t256 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t292 + Ifges(4,4) * t294) * qJD(1);
t322 = mrSges(4,1) * t219 - mrSges(4,2) * t220 + Ifges(4,5) * t272 + Ifges(4,6) * t273 + Ifges(4,3) * qJDD(3) + pkin(3) * t162 + qJ(4) * t148 + t289 * t134 + t287 * t135 + (t255 * t292 - t256 * t294) * qJD(1);
t270 = (-mrSges(4,1) * t294 + mrSges(4,2) * t292) * qJD(1);
t275 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t317;
t146 = m(4) * t220 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t273 - qJD(3) * t275 + t270 * t316 + t148;
t276 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t316;
t161 = m(4) * t219 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t272 + qJD(3) * t276 - t270 * t317 + t162;
t311 = t146 * t294 - t161 * t292;
t138 = m(3) * t238 - mrSges(3,1) * t297 - qJDD(1) * mrSges(3,2) + t311;
t147 = t150 * t289 + t151 * t287;
t301 = -m(4) * t227 + mrSges(4,1) * t273 - mrSges(4,2) * t272 - t275 * t317 + t276 * t316 - t147;
t142 = m(3) * t237 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t297 + t301;
t131 = t138 * t288 + t142 * t290;
t140 = t146 * t292 + t161 * t294;
t312 = t138 * t290 - t142 * t288;
t254 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t292 + Ifges(4,6) * t294) * qJD(1);
t127 = mrSges(4,2) * t227 - mrSges(4,3) * t219 + Ifges(4,1) * t272 + Ifges(4,4) * t273 + Ifges(4,5) * qJDD(3) - qJ(4) * t147 - qJD(3) * t255 - t134 * t287 + t135 * t289 + t254 * t316;
t304 = mrSges(7,1) * t172 - mrSges(7,3) * t170 - Ifges(7,4) * t198 - Ifges(7,2) * t267 - Ifges(7,6) * t197 + t235 * t201 - t205 * t234;
t300 = mrSges(6,2) * t175 - t234 * t206 - qJ(6) * (-t197 * mrSges(7,2) - t234 * t214 + t314) - pkin(5) * (-t198 * mrSges(7,2) - t235 * t214 + t309) - mrSges(6,1) * t174 - t235 * t204 + Ifges(6,6) * t197 - Ifges(6,5) * t198 - Ifges(6,3) * t267 + t304;
t298 = mrSges(5,1) * t182 - mrSges(5,2) * t183 + Ifges(5,5) * t245 + Ifges(5,6) * t244 + pkin(4) * t152 + t265 * t230 - t264 * t231 - t300;
t133 = Ifges(4,4) * t272 + qJD(3) * t256 - mrSges(4,1) * t227 + mrSges(4,3) * t220 - t298 - pkin(3) * t147 + Ifges(4,6) * qJDD(3) - t254 * t317 + (Ifges(4,2) + Ifges(5,3)) * t273;
t305 = mrSges(3,1) * t237 - mrSges(3,2) * t238 + Ifges(3,3) * qJDD(1) + pkin(2) * t301 + pkin(7) * t311 + t127 * t292 + t133 * t294;
t303 = mrSges(2,1) * t277 - mrSges(2,2) * t278 + Ifges(2,3) * qJDD(1) + pkin(1) * t131 + t305;
t129 = m(2) * t278 - mrSges(2,1) * t297 - qJDD(1) * mrSges(2,2) + t312;
t128 = m(2) * t277 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t297 + t131;
t125 = -mrSges(3,1) * t286 + mrSges(3,3) * t238 + t297 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t140 - t322;
t124 = mrSges(3,2) * t286 - mrSges(3,3) * t237 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t297 - pkin(7) * t140 + t127 * t294 - t133 * t292;
t123 = -mrSges(2,2) * g(3) - mrSges(2,3) * t277 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t297 - qJ(2) * t131 + t124 * t290 - t125 * t288;
t122 = Ifges(2,6) * qJDD(1) + t297 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t278 + t288 * t124 + t290 * t125 - pkin(1) * (m(3) * t286 + t140) + qJ(2) * t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t295 * t123 - t293 * t122 - pkin(6) * (t128 * t295 + t129 * t293), t123, t124, t127, t135, t154, -t203 * t234 + t306; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t293 * t123 + t295 * t122 + pkin(6) * (-t128 * t293 + t129 * t295), t122, t125, t133, t134, t153, -t304; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t303, t303, t305, t322, -Ifges(5,3) * t273 + t298, -t300, Ifges(7,5) * t198 + Ifges(7,6) * t267 + Ifges(7,3) * t197 + t235 * t203 - t280 * t205 - t308;];
m_new  = t1;
