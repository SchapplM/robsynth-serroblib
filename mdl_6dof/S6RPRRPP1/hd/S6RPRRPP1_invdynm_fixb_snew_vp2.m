% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:15:22
% EndTime: 2019-05-05 21:15:42
% DurationCPUTime: 9.64s
% Computational Cost: add. (177248->337), mult. (350973->410), div. (0->0), fcn. (227407->10), ass. (0->128)
t290 = sin(qJ(1));
t293 = cos(qJ(1));
t276 = g(1) * t290 - g(2) * t293;
t267 = qJDD(1) * pkin(1) + t276;
t277 = -g(1) * t293 - g(2) * t290;
t295 = qJD(1) ^ 2;
t269 = -pkin(1) * t295 + t277;
t286 = sin(pkin(9));
t287 = cos(pkin(9));
t243 = t267 * t286 + t269 * t287;
t229 = -pkin(2) * t295 + qJDD(1) * pkin(7) + t243;
t284 = -g(3) + qJDD(2);
t289 = sin(qJ(3));
t292 = cos(qJ(3));
t219 = -t229 * t289 + t292 * t284;
t270 = (-pkin(3) * t292 - pkin(8) * t289) * qJD(1);
t294 = qJD(3) ^ 2;
t317 = qJD(1) * t289;
t212 = -qJDD(3) * pkin(3) - t294 * pkin(8) + t270 * t317 - t219;
t288 = sin(qJ(4));
t291 = cos(qJ(4));
t266 = qJD(3) * t288 + t291 * t317;
t315 = qJD(1) * qJD(3);
t312 = t292 * t315;
t271 = qJDD(1) * t289 + t312;
t237 = -qJD(4) * t266 + qJDD(3) * t291 - t271 * t288;
t316 = t292 * qJD(1);
t279 = qJD(4) - t316;
t246 = pkin(4) * t279 - qJ(5) * t266;
t265 = qJD(3) * t291 - t288 * t317;
t263 = t265 ^ 2;
t182 = -t237 * pkin(4) - t263 * qJ(5) + t246 * t266 + qJDD(5) + t212;
t238 = qJD(4) * t265 + qJDD(3) * t288 + t271 * t291;
t285 = sin(pkin(10));
t320 = cos(pkin(10));
t206 = -t237 * t320 + t238 * t285;
t207 = t237 * t285 + t238 * t320;
t240 = -t265 * t320 + t266 * t285;
t241 = t265 * t285 + t266 * t320;
t177 = -0.2e1 * qJD(6) * t241 + (t240 * t279 - t207) * qJ(6) + (t241 * t279 + t206) * pkin(5) + t182;
t223 = -mrSges(7,1) * t279 + mrSges(7,2) * t241;
t224 = -mrSges(7,2) * t240 + mrSges(7,3) * t279;
t167 = m(7) * t177 + mrSges(7,1) * t206 - mrSges(7,3) * t207 - t223 * t241 + t224 * t240;
t242 = t287 * t267 - t269 * t286;
t228 = -qJDD(1) * pkin(2) - t295 * pkin(7) - t242;
t313 = t289 * t315;
t272 = qJDD(1) * t292 - t313;
t205 = (-t271 - t312) * pkin(8) + (-t272 + t313) * pkin(3) + t228;
t220 = t229 * t292 + t284 * t289;
t213 = -pkin(3) * t294 + qJDD(3) * pkin(8) + t270 * t316 + t220;
t183 = t205 * t291 - t288 * t213;
t264 = qJDD(4) - t272;
t179 = (t265 * t279 - t238) * qJ(5) + (t265 * t266 + t264) * pkin(4) + t183;
t184 = t205 * t288 + t213 * t291;
t181 = -pkin(4) * t263 + qJ(5) * t237 - t246 * t279 + t184;
t322 = -2 * qJD(5);
t175 = t179 * t285 + t181 * t320 + t240 * t322;
t203 = Ifges(7,1) * t241 + Ifges(7,4) * t279 + Ifges(7,5) * t240;
t204 = Ifges(6,1) * t241 - Ifges(6,4) * t240 + Ifges(6,5) * t279;
t214 = pkin(5) * t240 - qJ(6) * t241;
t278 = t279 ^ 2;
t170 = -pkin(5) * t278 + qJ(6) * t264 + 0.2e1 * qJD(6) * t279 - t214 * t240 + t175;
t307 = -mrSges(7,1) * t177 + mrSges(7,2) * t170;
t201 = Ifges(7,4) * t241 + Ifges(7,2) * t279 + Ifges(7,6) * t240;
t319 = -Ifges(6,5) * t241 + Ifges(6,6) * t240 - Ifges(6,3) * t279 - t201;
t153 = -mrSges(6,1) * t182 + mrSges(6,3) * t175 - pkin(5) * t167 + (t203 + t204) * t279 + (Ifges(6,6) - Ifges(7,6)) * t264 + t319 * t241 + (Ifges(6,4) - Ifges(7,5)) * t207 + (-Ifges(6,2) - Ifges(7,3)) * t206 + t307;
t305 = t179 * t320 - t285 * t181;
t174 = t241 * t322 + t305;
t202 = Ifges(6,4) * t241 - Ifges(6,2) * t240 + Ifges(6,6) * t279;
t172 = -t264 * pkin(5) - t278 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t214) * t241 - t305;
t199 = Ifges(7,5) * t241 + Ifges(7,6) * t279 + Ifges(7,3) * t240;
t304 = mrSges(7,2) * t172 - mrSges(7,3) * t177 + Ifges(7,1) * t207 + Ifges(7,4) * t264 + Ifges(7,5) * t206 + t199 * t279;
t154 = mrSges(6,2) * t182 - mrSges(6,3) * t174 + Ifges(6,1) * t207 - Ifges(6,4) * t206 + Ifges(6,5) * t264 - qJ(6) * t167 - t279 * t202 + t240 * t319 + t304;
t230 = Ifges(5,5) * t266 + Ifges(5,6) * t265 + Ifges(5,3) * t279;
t232 = Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * t279;
t221 = -mrSges(6,2) * t279 - mrSges(6,3) * t240;
t222 = mrSges(6,1) * t279 - mrSges(6,3) * t241;
t300 = m(6) * t182 + mrSges(6,1) * t206 + mrSges(6,2) * t207 + t221 * t240 + t222 * t241 + t167;
t314 = m(7) * t170 + mrSges(7,3) * t264 + t223 * t279;
t215 = mrSges(7,1) * t240 - mrSges(7,3) * t241;
t318 = -mrSges(6,1) * t240 - mrSges(6,2) * t241 - t215;
t321 = -mrSges(6,3) - mrSges(7,2);
t157 = m(6) * t175 - t264 * mrSges(6,2) + t206 * t321 - t279 * t222 + t240 * t318 + t314;
t308 = -m(7) * t172 + mrSges(7,1) * t264 + t224 * t279;
t159 = m(6) * t174 + t264 * mrSges(6,1) + t207 * t321 + t279 * t221 + t241 * t318 + t308;
t309 = t157 * t320 - t159 * t285;
t134 = -mrSges(5,1) * t212 + mrSges(5,3) * t184 + Ifges(5,4) * t238 + Ifges(5,2) * t237 + Ifges(5,6) * t264 - pkin(4) * t300 + qJ(5) * t309 + t153 * t320 + t285 * t154 - t266 * t230 + t279 * t232;
t152 = t157 * t285 + t159 * t320;
t231 = Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * t279;
t135 = mrSges(5,2) * t212 - mrSges(5,3) * t183 + Ifges(5,1) * t238 + Ifges(5,4) * t237 + Ifges(5,5) * t264 - qJ(5) * t152 - t153 * t285 + t154 * t320 + t230 * t265 - t231 * t279;
t244 = -mrSges(5,1) * t265 + mrSges(5,2) * t266;
t245 = -mrSges(5,2) * t279 + mrSges(5,3) * t265;
t150 = m(5) * t183 + mrSges(5,1) * t264 - mrSges(5,3) * t238 - t244 * t266 + t245 * t279 + t152;
t247 = mrSges(5,1) * t279 - mrSges(5,3) * t266;
t151 = m(5) * t184 - mrSges(5,2) * t264 + mrSges(5,3) * t237 + t244 * t265 - t247 * t279 + t309;
t148 = -t150 * t288 + t151 * t291;
t162 = -m(5) * t212 + mrSges(5,1) * t237 - mrSges(5,2) * t238 + t245 * t265 - t247 * t266 - t300;
t254 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t289 + Ifges(4,2) * t292) * qJD(1);
t255 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t289 + Ifges(4,4) * t292) * qJD(1);
t323 = mrSges(4,1) * t219 - mrSges(4,2) * t220 + Ifges(4,5) * t271 + Ifges(4,6) * t272 + Ifges(4,3) * qJDD(3) + pkin(3) * t162 + pkin(8) * t148 + t291 * t134 + t288 * t135 + (t254 * t289 - t255 * t292) * qJD(1);
t268 = (-mrSges(4,1) * t292 + mrSges(4,2) * t289) * qJD(1);
t274 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t317;
t146 = m(4) * t220 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t272 - qJD(3) * t274 + t268 * t316 + t148;
t275 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t316;
t161 = m(4) * t219 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t271 + qJD(3) * t275 - t268 * t317 + t162;
t310 = t146 * t292 - t161 * t289;
t138 = m(3) * t243 - mrSges(3,1) * t295 - qJDD(1) * mrSges(3,2) + t310;
t147 = t150 * t291 + t151 * t288;
t299 = -m(4) * t228 + mrSges(4,1) * t272 - mrSges(4,2) * t271 - t274 * t317 + t275 * t316 - t147;
t142 = m(3) * t242 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t295 + t299;
t131 = t138 * t286 + t142 * t287;
t140 = t146 * t289 + t161 * t292;
t311 = t138 * t287 - t142 * t286;
t253 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t289 + Ifges(4,6) * t292) * qJD(1);
t127 = mrSges(4,2) * t228 - mrSges(4,3) * t219 + Ifges(4,1) * t271 + Ifges(4,4) * t272 + Ifges(4,5) * qJDD(3) - pkin(8) * t147 - qJD(3) * t254 - t134 * t288 + t135 * t291 + t253 * t316;
t302 = mrSges(7,1) * t172 - mrSges(7,3) * t170 - Ifges(7,4) * t207 - Ifges(7,2) * t264 - Ifges(7,6) * t206 + t199 * t241 - t203 * t240;
t298 = mrSges(6,2) * t175 - t240 * t204 - qJ(6) * (-t206 * mrSges(7,2) - t240 * t215 + t314) - pkin(5) * (-t207 * mrSges(7,2) - t241 * t215 + t308) - mrSges(6,1) * t174 - t241 * t202 + Ifges(6,6) * t206 - Ifges(6,5) * t207 - Ifges(6,3) * t264 + t302;
t296 = mrSges(5,1) * t183 - mrSges(5,2) * t184 + Ifges(5,5) * t238 + Ifges(5,6) * t237 + Ifges(5,3) * t264 + pkin(4) * t152 + t231 * t266 - t232 * t265 - t298;
t133 = -mrSges(4,1) * t228 + mrSges(4,3) * t220 + Ifges(4,4) * t271 + Ifges(4,2) * t272 + Ifges(4,6) * qJDD(3) - pkin(3) * t147 + qJD(3) * t255 - t253 * t317 - t296;
t303 = mrSges(3,1) * t242 - mrSges(3,2) * t243 + Ifges(3,3) * qJDD(1) + pkin(2) * t299 + pkin(7) * t310 + t127 * t289 + t133 * t292;
t301 = mrSges(2,1) * t276 - mrSges(2,2) * t277 + Ifges(2,3) * qJDD(1) + pkin(1) * t131 + t303;
t129 = m(2) * t277 - mrSges(2,1) * t295 - qJDD(1) * mrSges(2,2) + t311;
t128 = m(2) * t276 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t295 + t131;
t125 = -mrSges(3,1) * t284 + mrSges(3,3) * t243 + t295 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t140 - t323;
t124 = mrSges(3,2) * t284 - mrSges(3,3) * t242 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t295 - pkin(7) * t140 + t127 * t292 - t133 * t289;
t123 = -mrSges(2,2) * g(3) - mrSges(2,3) * t276 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t295 - qJ(2) * t131 + t124 * t287 - t125 * t286;
t122 = Ifges(2,6) * qJDD(1) + t295 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t277 + t286 * t124 + t287 * t125 - pkin(1) * (m(3) * t284 + t140) + qJ(2) * t311;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t293 * t123 - t290 * t122 - pkin(6) * (t128 * t293 + t129 * t290), t123, t124, t127, t135, t154, -t201 * t240 + t304; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t290 * t123 + t293 * t122 + pkin(6) * (-t128 * t290 + t129 * t293), t122, t125, t133, t134, t153, -t302; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t301, t301, t303, t323, t296, -t298, Ifges(7,5) * t207 + Ifges(7,6) * t264 + Ifges(7,3) * t206 + t241 * t201 - t279 * t203 - t307;];
m_new  = t1;
