% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-05-08 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR10V2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 20:31:52
% EndTime: 2019-05-08 20:32:51
% DurationCPUTime: 14.92s
% Computational Cost: add. (272158->351), mult. (537157->446), div. (0->0), fcn. (409983->12), ass. (0->137)
t263 = sin(qJ(3));
t264 = sin(qJ(2));
t269 = cos(qJ(3));
t270 = cos(qJ(2));
t239 = (t263 * t264 - t269 * t270) * qJD(1);
t265 = sin(qJ(1));
t271 = cos(qJ(1));
t253 = -t271 * g(1) - t265 * g(2);
t272 = qJD(1) ^ 2;
t246 = -t272 * pkin(1) + t253;
t233 = -t270 * g(3) - t264 * t246;
t226 = (t264 * t270 * t272 + qJDD(2)) * pkin(2) + t233;
t234 = -t264 * g(3) + t270 * t246;
t227 = (-t270 ^ 2 * t272 - qJD(2) ^ 2) * pkin(2) + t234;
t203 = t263 * t226 + t269 * t227;
t240 = (t263 * t270 + t264 * t269) * qJD(1);
t288 = qJD(1) * qJD(2);
t247 = t264 * qJDD(1) + t270 * t288;
t287 = t264 * t288;
t248 = t270 * qJDD(1) - t287;
t213 = -t240 * qJD(3) - t263 * t247 + t269 * t248;
t221 = t239 * mrSges(4,1) + t240 * mrSges(4,2);
t258 = qJD(2) + qJD(3);
t232 = t258 * mrSges(4,1) - t240 * mrSges(4,3);
t257 = qJDD(2) + qJDD(3);
t223 = t239 * pkin(3) - t240 * pkin(5);
t256 = t258 ^ 2;
t187 = -t256 * pkin(3) + t257 * pkin(5) - t239 * t223 + t203;
t262 = sin(qJ(4));
t268 = cos(qJ(4));
t214 = -t239 * qJD(3) + t269 * t247 + t263 * t248;
t252 = t265 * g(1) - t271 * g(2);
t244 = -qJDD(1) * pkin(1) - t252;
t224 = (-t248 + t287) * pkin(2) + t244;
t275 = (t239 * t258 - t214) * pkin(5) + (t240 * t258 - t213) * pkin(3) + t224;
t171 = t268 * t187 + t262 * t275;
t202 = t269 * t226 - t263 * t227;
t186 = -t257 * pkin(3) - t256 * pkin(5) + t240 * t223 - t202;
t261 = sin(qJ(5));
t267 = cos(qJ(5));
t166 = t267 * t171 + t261 * t186;
t230 = t268 * t240 + t262 * t258;
t191 = -t230 * qJD(4) - t262 * t214 + t268 * t257;
t190 = qJDD(5) - t191;
t235 = qJD(4) + t239;
t207 = -t261 * t230 + t267 * t235;
t208 = t267 * t230 + t261 * t235;
t162 = (-t207 * t208 + t190) * pkin(6) + t166;
t170 = t262 * t187 - t268 * t275;
t229 = -t262 * t240 + t268 * t258;
t192 = t229 * qJD(4) + t268 * t214 + t262 * t257;
t212 = qJDD(4) - t213;
t177 = t207 * qJD(5) + t267 * t192 + t261 * t212;
t228 = qJD(5) - t229;
t163 = (-t207 * t228 - t177) * pkin(6) + t170;
t260 = sin(qJ(6));
t266 = cos(qJ(6));
t160 = -t260 * t162 + t266 * t163;
t193 = -t260 * t208 + t266 * t228;
t168 = t193 * qJD(6) + t266 * t177 + t260 * t190;
t176 = -t208 * qJD(5) - t261 * t192 + t267 * t212;
t175 = qJDD(6) - t176;
t194 = t266 * t208 + t260 * t228;
t179 = -mrSges(7,1) * t193 + mrSges(7,2) * t194;
t206 = qJD(6) - t207;
t180 = -mrSges(7,2) * t206 + mrSges(7,3) * t193;
t158 = m(7) * t160 + mrSges(7,1) * t175 - mrSges(7,3) * t168 - t179 * t194 + t180 * t206;
t161 = t266 * t162 + t260 * t163;
t167 = -t194 * qJD(6) - t260 * t177 + t266 * t190;
t181 = mrSges(7,1) * t206 - mrSges(7,3) * t194;
t159 = m(7) * t161 - mrSges(7,2) * t175 + mrSges(7,3) * t167 + t179 * t193 - t181 * t206;
t153 = -t260 * t158 + t266 * t159;
t188 = -mrSges(6,1) * t207 + mrSges(6,2) * t208;
t196 = mrSges(6,1) * t228 - mrSges(6,3) * t208;
t152 = m(6) * t166 - t190 * mrSges(6,2) + t176 * mrSges(6,3) + t207 * t188 - t228 * t196 + t153;
t165 = -t261 * t171 + t267 * t186;
t164 = (-t208 ^ 2 - t228 ^ 2) * pkin(6) - t165;
t195 = -mrSges(6,2) * t228 + mrSges(6,3) * t207;
t156 = m(6) * t165 - m(7) * t164 + mrSges(6,1) * t190 + mrSges(7,1) * t167 - mrSges(7,2) * t168 - mrSges(6,3) * t177 + t180 * t193 - t181 * t194 - t188 * t208 + t195 * t228;
t204 = -mrSges(5,1) * t229 + mrSges(5,2) * t230;
t216 = mrSges(5,1) * t235 - mrSges(5,3) * t230;
t146 = m(5) * t171 - t212 * mrSges(5,2) + t191 * mrSges(5,3) + t267 * t152 - t261 * t156 + t229 * t204 - t235 * t216;
t215 = -mrSges(5,2) * t235 + mrSges(5,3) * t229;
t285 = t266 * t158 + t260 * t159;
t150 = t212 * mrSges(5,1) + t176 * mrSges(6,1) - t177 * mrSges(6,2) - t192 * mrSges(5,3) + t207 * t195 - t208 * t196 - t230 * t204 + t235 * t215 + (-m(5) - m(6)) * t170 - t285;
t286 = t268 * t146 - t262 * t150;
t135 = m(4) * t203 - t257 * mrSges(4,2) + t213 * mrSges(4,3) - t239 * t221 - t258 * t232 + t286;
t231 = -t258 * mrSges(4,2) - t239 * mrSges(4,3);
t279 = -m(5) * t186 + t191 * mrSges(5,1) - t192 * mrSges(5,2) - t261 * t152 - t267 * t156 + t229 * t215 - t230 * t216;
t143 = m(4) * t202 + t257 * mrSges(4,1) - t214 * mrSges(4,3) - t240 * t221 + t258 * t231 + t279;
t129 = t263 * t135 + t269 * t143;
t237 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t264 + Ifges(3,2) * t270) * qJD(1);
t238 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t264 + Ifges(3,4) * t270) * qJD(1);
t172 = Ifges(7,5) * t194 + Ifges(7,6) * t193 + Ifges(7,3) * t206;
t174 = Ifges(7,1) * t194 + Ifges(7,4) * t193 + Ifges(7,5) * t206;
t154 = -mrSges(7,1) * t164 + mrSges(7,3) * t161 + Ifges(7,4) * t168 + Ifges(7,2) * t167 + Ifges(7,6) * t175 - t172 * t194 + t174 * t206;
t173 = Ifges(7,4) * t194 + Ifges(7,2) * t193 + Ifges(7,6) * t206;
t155 = mrSges(7,2) * t164 - mrSges(7,3) * t160 + Ifges(7,1) * t168 + Ifges(7,4) * t167 + Ifges(7,5) * t175 + t172 * t193 - t173 * t206;
t182 = Ifges(6,5) * t208 + Ifges(6,6) * t207 + Ifges(6,3) * t228;
t183 = Ifges(6,4) * t208 + Ifges(6,2) * t207 + Ifges(6,6) * t228;
t141 = mrSges(6,2) * t170 - mrSges(6,3) * t165 + Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * t190 - pkin(6) * t285 - t260 * t154 + t266 * t155 + t207 * t182 - t228 * t183;
t184 = Ifges(6,1) * t208 + Ifges(6,4) * t207 + Ifges(6,5) * t228;
t278 = mrSges(7,1) * t160 - mrSges(7,2) * t161 + Ifges(7,5) * t168 + Ifges(7,6) * t167 + Ifges(7,3) * t175 + t173 * t194 - t193 * t174;
t151 = -mrSges(6,1) * t170 + mrSges(6,3) * t166 + Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * t190 - t182 * t208 + t184 * t228 - t278;
t197 = Ifges(5,5) * t230 + Ifges(5,6) * t229 + Ifges(5,3) * t235;
t198 = Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t235;
t131 = mrSges(5,2) * t186 + mrSges(5,3) * t170 + Ifges(5,1) * t192 + Ifges(5,4) * t191 + Ifges(5,5) * t212 + t267 * t141 - t261 * t151 + t229 * t197 - t235 * t198;
t199 = Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t235;
t274 = mrSges(6,1) * t165 - mrSges(6,2) * t166 + Ifges(6,5) * t177 + Ifges(6,6) * t176 + Ifges(6,3) * t190 + pkin(6) * t153 + t266 * t154 + t260 * t155 + t208 * t183 - t207 * t184;
t140 = -mrSges(5,1) * t186 + mrSges(5,3) * t171 + Ifges(5,4) * t192 + Ifges(5,2) * t191 + Ifges(5,6) * t212 - t230 * t197 + t235 * t199 - t274;
t218 = Ifges(4,4) * t240 - Ifges(4,2) * t239 + Ifges(4,6) * t258;
t219 = Ifges(4,1) * t240 - Ifges(4,4) * t239 + Ifges(4,5) * t258;
t280 = -mrSges(4,1) * t202 + mrSges(4,2) * t203 - Ifges(4,5) * t214 - Ifges(4,6) * t213 - Ifges(4,3) * t257 - pkin(3) * t279 - pkin(5) * t286 - t262 * t131 - t268 * t140 - t240 * t218 - t239 * t219;
t291 = mrSges(3,1) * t233 - mrSges(3,2) * t234 + Ifges(3,5) * t247 + Ifges(3,6) * t248 + Ifges(3,3) * qJDD(2) + pkin(2) * t129 + (t264 * t237 - t270 * t238) * qJD(1) - t280;
t137 = t262 * t146 + t268 * t150;
t290 = qJD(1) * t264;
t289 = qJD(1) * t270;
t217 = Ifges(4,5) * t240 - Ifges(4,6) * t239 + Ifges(4,3) * t258;
t125 = mrSges(4,2) * t224 - mrSges(4,3) * t202 + Ifges(4,1) * t214 + Ifges(4,4) * t213 + Ifges(4,5) * t257 - pkin(5) * t137 + t268 * t131 - t262 * t140 - t239 * t217 - t258 * t218;
t276 = mrSges(5,1) * t170 + mrSges(5,2) * t171 - Ifges(5,5) * t192 - Ifges(5,6) * t191 - Ifges(5,3) * t212 - t261 * t141 - t267 * t151 - t230 * t198 + t229 * t199;
t126 = -mrSges(4,1) * t224 + mrSges(4,3) * t203 + Ifges(4,4) * t214 + Ifges(4,2) * t213 + Ifges(4,6) * t257 - pkin(3) * t137 - t240 * t217 + t258 * t219 + t276;
t236 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t264 + Ifges(3,6) * t270) * qJD(1);
t281 = m(4) * t224 - t213 * mrSges(4,1) + t214 * mrSges(4,2) + t239 * t231 + t240 * t232 + t137;
t121 = -mrSges(3,1) * t244 + mrSges(3,3) * t234 + Ifges(3,4) * t247 + Ifges(3,2) * t248 + Ifges(3,6) * qJDD(2) - pkin(2) * t281 + qJD(2) * t238 + t263 * t125 + t269 * t126 - t236 * t290;
t123 = mrSges(3,2) * t244 - mrSges(3,3) * t233 + Ifges(3,1) * t247 + Ifges(3,4) * t248 + Ifges(3,5) * qJDD(2) - qJD(2) * t237 + t269 * t125 - t263 * t126 + t236 * t289;
t250 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t290;
t251 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t289;
t277 = -m(3) * t244 + t248 * mrSges(3,1) - t247 * mrSges(3,2) - t250 * t290 + t251 * t289 - t281;
t282 = mrSges(2,1) * t252 - mrSges(2,2) * t253 + Ifges(2,3) * qJDD(1) + pkin(1) * t277 + t270 * t121 + t264 * t123;
t245 = (-mrSges(3,1) * t270 + mrSges(3,2) * t264) * qJD(1);
t132 = m(2) * t252 + qJDD(1) * mrSges(2,1) - t272 * mrSges(2,2) + t277;
t128 = m(3) * t234 - qJDD(2) * mrSges(3,2) + t248 * mrSges(3,3) - qJD(2) * t250 + t269 * t135 - t263 * t143 + t245 * t289;
t127 = m(3) * t233 + qJDD(2) * mrSges(3,1) - t247 * mrSges(3,3) + qJD(2) * t251 - t245 * t290 + t129;
t124 = m(2) * t253 - t272 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t264 * t127 + t270 * t128;
t119 = -pkin(1) * (t270 * t127 + t264 * t128) + t272 * Ifges(2,5) + mrSges(2,3) * t253 + mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) - t291;
t118 = -mrSges(2,2) * g(3) - mrSges(2,3) * t252 + Ifges(2,5) * qJDD(1) - t272 * Ifges(2,6) - t264 * t121 + t270 * t123;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t271 * t118 - t265 * t119 - pkin(4) * (t265 * t124 + t271 * t132), t118, t123, t125, t131, t141, t155; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t265 * t118 + t271 * t119 + pkin(4) * (t271 * t124 - t265 * t132), t119, t121, t126, t140, t151, t154; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t282, t282, t291, -t280, -t276, t274, t278;];
m_new  = t1;
