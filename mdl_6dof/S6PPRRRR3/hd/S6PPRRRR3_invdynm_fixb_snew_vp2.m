% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 21:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:12:34
% EndTime: 2019-05-04 21:13:44
% DurationCPUTime: 69.65s
% Computational Cost: add. (1419786->302), mult. (2699951->412), div. (0->0), fcn. (2288019->18), ass. (0->148)
t272 = sin(pkin(13));
t277 = cos(pkin(13));
t266 = -g(1) * t277 - g(2) * t272;
t271 = sin(pkin(14));
t276 = cos(pkin(14));
t265 = g(1) * t272 - g(2) * t277;
t270 = -g(3) + qJDD(1);
t275 = sin(pkin(6));
t280 = cos(pkin(6));
t297 = t265 * t280 + t270 * t275;
t239 = -t271 * t266 + t297 * t276;
t240 = t276 * t266 + t297 * t271;
t252 = -t265 * t275 + t270 * t280 + qJDD(2);
t284 = sin(qJ(3));
t279 = cos(pkin(7));
t288 = cos(qJ(3));
t304 = t279 * t288;
t274 = sin(pkin(7));
t307 = t274 * t288;
t214 = t239 * t304 - t240 * t284 + t252 * t307;
t289 = qJD(3) ^ 2;
t273 = sin(pkin(8));
t316 = pkin(10) * t273;
t209 = qJDD(3) * pkin(3) + t289 * t316 + t214;
t305 = t279 * t284;
t308 = t274 * t284;
t215 = t239 * t305 + t288 * t240 + t252 * t308;
t210 = -pkin(3) * t289 + qJDD(3) * t316 + t215;
t223 = -t239 * t274 + t252 * t279;
t287 = cos(qJ(4));
t278 = cos(pkin(8));
t283 = sin(qJ(4));
t306 = t278 * t283;
t309 = t273 * t283;
t202 = t209 * t306 + t287 * t210 + t223 * t309;
t269 = qJD(3) * t278 + qJD(4);
t303 = qJD(3) * t273;
t301 = t283 * t303;
t256 = mrSges(5,1) * t269 - mrSges(5,3) * t301;
t258 = (-mrSges(5,1) * t287 + mrSges(5,2) * t283) * t303;
t302 = qJD(3) * qJD(4);
t261 = (-qJDD(3) * t287 + t283 * t302) * t273;
t268 = qJDD(3) * t278 + qJDD(4);
t259 = (-pkin(4) * t287 - pkin(11) * t283) * t303;
t267 = t269 ^ 2;
t300 = t287 * t303;
t200 = -pkin(4) * t267 + pkin(11) * t268 + t259 * t300 + t202;
t222 = t278 * t223;
t260 = (qJDD(3) * t283 + t287 * t302) * t273;
t204 = t261 * pkin(4) - t260 * pkin(11) + t222 + (-t209 + (pkin(4) * t283 - pkin(11) * t287) * t269 * qJD(3)) * t273;
t282 = sin(qJ(5));
t286 = cos(qJ(5));
t196 = t286 * t200 + t282 * t204;
t253 = t269 * t286 - t282 * t301;
t254 = t269 * t282 + t286 * t301;
t235 = -pkin(5) * t253 - pkin(12) * t254;
t255 = qJDD(5) + t261;
t264 = qJD(5) - t300;
t263 = t264 ^ 2;
t194 = -pkin(5) * t263 + pkin(12) * t255 + t235 * t253 + t196;
t201 = -t283 * t210 + (t209 * t278 + t223 * t273) * t287;
t199 = -t268 * pkin(4) - t267 * pkin(11) + t259 * t301 - t201;
t232 = -qJD(5) * t254 - t260 * t282 + t268 * t286;
t233 = qJD(5) * t253 + t260 * t286 + t268 * t282;
t197 = (-t253 * t264 - t233) * pkin(12) + (t254 * t264 - t232) * pkin(5) + t199;
t281 = sin(qJ(6));
t285 = cos(qJ(6));
t190 = -t194 * t281 + t197 * t285;
t241 = -t254 * t281 + t264 * t285;
t213 = qJD(6) * t241 + t233 * t285 + t255 * t281;
t242 = t254 * t285 + t264 * t281;
t220 = -mrSges(7,1) * t241 + mrSges(7,2) * t242;
t251 = qJD(6) - t253;
t224 = -mrSges(7,2) * t251 + mrSges(7,3) * t241;
t230 = qJDD(6) - t232;
t188 = m(7) * t190 + mrSges(7,1) * t230 - mrSges(7,3) * t213 - t220 * t242 + t224 * t251;
t191 = t194 * t285 + t197 * t281;
t212 = -qJD(6) * t242 - t233 * t281 + t255 * t285;
t225 = mrSges(7,1) * t251 - mrSges(7,3) * t242;
t189 = m(7) * t191 - mrSges(7,2) * t230 + mrSges(7,3) * t212 + t220 * t241 - t225 * t251;
t182 = -t188 * t281 + t285 * t189;
t234 = -mrSges(6,1) * t253 + mrSges(6,2) * t254;
t244 = mrSges(6,1) * t264 - mrSges(6,3) * t254;
t180 = m(6) * t196 - mrSges(6,2) * t255 + mrSges(6,3) * t232 + t234 * t253 - t244 * t264 + t182;
t195 = -t200 * t282 + t204 * t286;
t193 = -pkin(5) * t255 - pkin(12) * t263 + t235 * t254 - t195;
t192 = -m(7) * t193 + t212 * mrSges(7,1) - mrSges(7,2) * t213 + t241 * t224 - t225 * t242;
t243 = -mrSges(6,2) * t264 + mrSges(6,3) * t253;
t186 = m(6) * t195 + mrSges(6,1) * t255 - mrSges(6,3) * t233 - t234 * t254 + t243 * t264 + t192;
t299 = t286 * t180 - t186 * t282;
t171 = m(5) * t202 - mrSges(5,2) * t268 - mrSges(5,3) * t261 - t256 * t269 + t258 * t300 + t299;
t174 = t282 * t180 + t286 * t186;
t205 = -t273 * t209 + t222;
t257 = -mrSges(5,2) * t269 + mrSges(5,3) * t300;
t173 = m(5) * t205 + t261 * mrSges(5,1) + t260 * mrSges(5,2) + (t256 * t283 - t257 * t287) * t303 + t174;
t181 = t188 * t285 + t189 * t281;
t292 = -m(6) * t199 + t232 * mrSges(6,1) - mrSges(6,2) * t233 + t253 * t243 - t244 * t254 - t181;
t177 = m(5) * t201 + mrSges(5,1) * t268 - mrSges(5,3) * t260 + t257 * t269 - t258 * t301 + t292;
t310 = t177 * t287;
t160 = t171 * t306 - t173 * t273 + t278 * t310;
t156 = m(4) * t214 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t289 + t160;
t159 = t171 * t309 + t278 * t173 + t273 * t310;
t158 = m(4) * t223 + t159;
t165 = t287 * t171 - t177 * t283;
t164 = m(4) * t215 - mrSges(4,1) * t289 - qJDD(3) * mrSges(4,2) + t165;
t146 = t156 * t304 - t158 * t274 + t164 * t305;
t143 = m(3) * t239 + t146;
t150 = -t156 * t284 + t288 * t164;
t149 = m(3) * t240 + t150;
t318 = t143 * t276 + t149 * t271;
t216 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t251;
t218 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t251;
t183 = -mrSges(7,1) * t193 + mrSges(7,3) * t191 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t230 - t216 * t242 + t218 * t251;
t217 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t251;
t184 = mrSges(7,2) * t193 - mrSges(7,3) * t190 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t230 + t216 * t241 - t217 * t251;
t226 = Ifges(6,5) * t254 + Ifges(6,6) * t253 + Ifges(6,3) * t264;
t227 = Ifges(6,4) * t254 + Ifges(6,2) * t253 + Ifges(6,6) * t264;
t166 = mrSges(6,2) * t199 - mrSges(6,3) * t195 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t255 - pkin(12) * t181 - t183 * t281 + t184 * t285 + t226 * t253 - t227 * t264;
t228 = Ifges(6,1) * t254 + Ifges(6,4) * t253 + Ifges(6,5) * t264;
t291 = mrSges(7,1) * t190 - mrSges(7,2) * t191 + Ifges(7,5) * t213 + Ifges(7,6) * t212 + Ifges(7,3) * t230 + t217 * t242 - t218 * t241;
t167 = -mrSges(6,1) * t199 + mrSges(6,3) * t196 + Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t255 - pkin(5) * t181 - t226 * t254 + t228 * t264 - t291;
t248 = Ifges(5,6) * t269 + (Ifges(5,4) * t283 + Ifges(5,2) * t287) * t303;
t249 = Ifges(5,5) * t269 + (Ifges(5,1) * t283 + Ifges(5,4) * t287) * t303;
t151 = Ifges(5,5) * t260 - Ifges(5,6) * t261 + Ifges(5,3) * t268 + mrSges(5,1) * t201 - mrSges(5,2) * t202 + t282 * t166 + t286 * t167 + pkin(4) * t292 + pkin(11) * t299 + (t248 * t283 - t249 * t287) * t303;
t247 = Ifges(5,3) * t269 + (Ifges(5,5) * t283 + Ifges(5,6) * t287) * t303;
t152 = mrSges(5,2) * t205 - mrSges(5,3) * t201 + Ifges(5,1) * t260 - Ifges(5,4) * t261 + Ifges(5,5) * t268 - pkin(11) * t174 + t166 * t286 - t167 * t282 + t247 * t300 - t248 * t269;
t290 = mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t233 + Ifges(6,6) * t232 + Ifges(6,3) * t255 + pkin(5) * t192 + pkin(12) * t182 + t285 * t183 + t281 * t184 + t254 * t227 - t253 * t228;
t153 = -mrSges(5,1) * t205 + mrSges(5,3) * t202 + Ifges(5,4) * t260 - Ifges(5,2) * t261 + Ifges(5,6) * t268 - pkin(4) * t174 - t247 * t301 + t269 * t249 - t290;
t294 = pkin(10) * t165 + t152 * t283 + t153 * t287;
t136 = mrSges(4,1) * t214 - mrSges(4,2) * t215 + Ifges(4,3) * qJDD(3) + pkin(3) * t160 + t278 * t151 + t294 * t273;
t145 = t156 * t307 + t279 * t158 + t164 * t308;
t137 = -mrSges(4,1) * t223 + mrSges(4,3) * t215 + t289 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t159 - t273 * t151 + t294 * t278;
t138 = mrSges(4,2) * t223 - mrSges(4,3) * t214 + Ifges(4,5) * qJDD(3) - t289 * Ifges(4,6) + t287 * t152 - t283 * t153 + (-t159 * t273 - t160 * t278) * pkin(10);
t295 = pkin(9) * t150 + t137 * t288 + t138 * t284;
t129 = -mrSges(3,1) * t252 + mrSges(3,3) * t240 - pkin(2) * t145 - t274 * t136 + t295 * t279;
t131 = mrSges(3,2) * t252 - mrSges(3,3) * t239 - t284 * t137 + t288 * t138 + (-t145 * t274 - t146 * t279) * pkin(9);
t141 = -t143 * t271 + t276 * t149;
t317 = qJ(2) * t141 + t129 * t276 + t131 * t271;
t144 = m(3) * t252 + t145;
t135 = -t144 * t275 + t280 * t318;
t127 = mrSges(3,1) * t239 - mrSges(3,2) * t240 + pkin(2) * t146 + t279 * t136 + t295 * t274;
t293 = mrSges(2,1) * t265 - mrSges(2,2) * t266 + pkin(1) * t135 + t280 * t127 + t275 * t317;
t139 = m(2) * t266 + t141;
t134 = t280 * t144 + t275 * t318;
t132 = m(2) * t265 + t135;
t125 = mrSges(2,2) * t270 - mrSges(2,3) * t265 - t271 * t129 + t276 * t131 + (-t134 * t275 - t135 * t280) * qJ(2);
t124 = -mrSges(2,1) * t270 + mrSges(2,3) * t266 - pkin(1) * t134 - t275 * t127 + t280 * t317;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t277 * t125 - t272 * t124 - qJ(1) * (t132 * t277 + t139 * t272), t125, t131, t138, t152, t166, t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t272 * t125 + t277 * t124 + qJ(1) * (-t132 * t272 + t139 * t277), t124, t129, t137, t153, t167, t183; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t293, t293, t127, t136, t151, t290, t291;];
m_new  = t1;
