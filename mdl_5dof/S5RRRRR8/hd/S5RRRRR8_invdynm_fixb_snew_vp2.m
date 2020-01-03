% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:33
% EndTime: 2019-12-31 22:24:52
% DurationCPUTime: 9.97s
% Computational Cost: add. (176107->312), mult. (354017->395), div. (0->0), fcn. (248857->10), ass. (0->126)
t266 = sin(qJ(2));
t271 = cos(qJ(2));
t289 = qJD(1) * qJD(2);
t247 = qJDD(1) * t266 + t271 * t289;
t267 = sin(qJ(1));
t272 = cos(qJ(1));
t254 = -g(1) * t272 - g(2) * t267;
t273 = qJD(1) ^ 2;
t241 = -pkin(1) * t273 + qJDD(1) * pkin(6) + t254;
t292 = t266 * t241;
t293 = pkin(2) * t273;
t203 = qJDD(2) * pkin(2) - t247 * pkin(7) - t292 + (pkin(7) * t289 + t266 * t293 - g(3)) * t271;
t228 = -g(3) * t266 + t271 * t241;
t248 = qJDD(1) * t271 - t266 * t289;
t291 = qJD(1) * t266;
t252 = qJD(2) * pkin(2) - pkin(7) * t291;
t262 = t271 ^ 2;
t204 = pkin(7) * t248 - qJD(2) * t252 - t262 * t293 + t228;
t265 = sin(qJ(3));
t270 = cos(qJ(3));
t186 = t265 * t203 + t270 * t204;
t239 = (t265 * t271 + t266 * t270) * qJD(1);
t212 = -t239 * qJD(3) - t265 * t247 + t248 * t270;
t290 = qJD(1) * t271;
t238 = -t265 * t291 + t270 * t290;
t222 = -mrSges(4,1) * t238 + mrSges(4,2) * t239;
t260 = qJD(2) + qJD(3);
t230 = mrSges(4,1) * t260 - mrSges(4,3) * t239;
t259 = qJDD(2) + qJDD(3);
t213 = qJD(3) * t238 + t247 * t270 + t248 * t265;
t253 = t267 * g(1) - t272 * g(2);
t283 = -qJDD(1) * pkin(1) - t253;
t214 = -t248 * pkin(2) + t252 * t291 + (-pkin(7) * t262 - pkin(6)) * t273 + t283;
t175 = (-t238 * t260 - t213) * pkin(8) + (t239 * t260 - t212) * pkin(3) + t214;
t223 = -pkin(3) * t238 - pkin(8) * t239;
t258 = t260 ^ 2;
t178 = -pkin(3) * t258 + pkin(8) * t259 + t223 * t238 + t186;
t264 = sin(qJ(4));
t269 = cos(qJ(4));
t164 = t175 * t269 - t264 * t178;
t225 = -t239 * t264 + t260 * t269;
t189 = qJD(4) * t225 + t213 * t269 + t259 * t264;
t211 = qJDD(4) - t212;
t226 = t239 * t269 + t260 * t264;
t234 = qJD(4) - t238;
t162 = (t225 * t234 - t189) * pkin(9) + (t225 * t226 + t211) * pkin(4) + t164;
t165 = t175 * t264 + t178 * t269;
t188 = -qJD(4) * t226 - t213 * t264 + t259 * t269;
t217 = pkin(4) * t234 - pkin(9) * t226;
t224 = t225 ^ 2;
t163 = -pkin(4) * t224 + pkin(9) * t188 - t217 * t234 + t165;
t263 = sin(qJ(5));
t268 = cos(qJ(5));
t160 = t162 * t268 - t163 * t263;
t196 = t225 * t268 - t226 * t263;
t171 = qJD(5) * t196 + t188 * t263 + t189 * t268;
t197 = t225 * t263 + t226 * t268;
t183 = -mrSges(6,1) * t196 + mrSges(6,2) * t197;
t232 = qJD(5) + t234;
t190 = -mrSges(6,2) * t232 + mrSges(6,3) * t196;
t206 = qJDD(5) + t211;
t155 = m(6) * t160 + mrSges(6,1) * t206 - mrSges(6,3) * t171 - t183 * t197 + t190 * t232;
t161 = t162 * t263 + t163 * t268;
t170 = -qJD(5) * t197 + t188 * t268 - t189 * t263;
t191 = mrSges(6,1) * t232 - mrSges(6,3) * t197;
t156 = m(6) * t161 - mrSges(6,2) * t206 + mrSges(6,3) * t170 + t183 * t196 - t191 * t232;
t147 = t155 * t268 + t156 * t263;
t201 = -mrSges(5,1) * t225 + mrSges(5,2) * t226;
t215 = -mrSges(5,2) * t234 + mrSges(5,3) * t225;
t145 = m(5) * t164 + mrSges(5,1) * t211 - mrSges(5,3) * t189 - t201 * t226 + t215 * t234 + t147;
t216 = mrSges(5,1) * t234 - mrSges(5,3) * t226;
t285 = -t155 * t263 + t156 * t268;
t146 = m(5) * t165 - mrSges(5,2) * t211 + mrSges(5,3) * t188 + t201 * t225 - t216 * t234 + t285;
t286 = -t145 * t264 + t146 * t269;
t138 = m(4) * t186 - mrSges(4,2) * t259 + mrSges(4,3) * t212 + t222 * t238 - t230 * t260 + t286;
t185 = t203 * t270 - t265 * t204;
t229 = -mrSges(4,2) * t260 + mrSges(4,3) * t238;
t177 = -pkin(3) * t259 - pkin(8) * t258 + t239 * t223 - t185;
t166 = -pkin(4) * t188 - pkin(9) * t224 + t217 * t226 + t177;
t281 = m(6) * t166 - mrSges(6,1) * t170 + mrSges(6,2) * t171 - t190 * t196 + t191 * t197;
t277 = -m(5) * t177 + mrSges(5,1) * t188 - mrSges(5,2) * t189 + t225 * t215 - t216 * t226 - t281;
t151 = m(4) * t185 + mrSges(4,1) * t259 - mrSges(4,3) * t213 - t222 * t239 + t229 * t260 + t277;
t131 = t138 * t265 + t151 * t270;
t227 = -t271 * g(3) - t292;
t236 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t266 + Ifges(3,2) * t271) * qJD(1);
t237 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t266 + Ifges(3,4) * t271) * qJD(1);
t179 = Ifges(6,5) * t197 + Ifges(6,6) * t196 + Ifges(6,3) * t232;
t181 = Ifges(6,1) * t197 + Ifges(6,4) * t196 + Ifges(6,5) * t232;
t148 = -mrSges(6,1) * t166 + mrSges(6,3) * t161 + Ifges(6,4) * t171 + Ifges(6,2) * t170 + Ifges(6,6) * t206 - t179 * t197 + t181 * t232;
t180 = Ifges(6,4) * t197 + Ifges(6,2) * t196 + Ifges(6,6) * t232;
t149 = mrSges(6,2) * t166 - mrSges(6,3) * t160 + Ifges(6,1) * t171 + Ifges(6,4) * t170 + Ifges(6,5) * t206 + t179 * t196 - t180 * t232;
t192 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * t234;
t194 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * t234;
t130 = -mrSges(5,1) * t177 + mrSges(5,3) * t165 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * t211 - pkin(4) * t281 + pkin(9) * t285 + t268 * t148 + t263 * t149 - t226 * t192 + t234 * t194;
t193 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * t234;
t133 = mrSges(5,2) * t177 - mrSges(5,3) * t164 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * t211 - pkin(9) * t147 - t148 * t263 + t149 * t268 + t192 * t225 - t234 * t193;
t219 = Ifges(4,4) * t239 + Ifges(4,2) * t238 + Ifges(4,6) * t260;
t220 = Ifges(4,1) * t239 + Ifges(4,4) * t238 + Ifges(4,5) * t260;
t278 = -mrSges(4,1) * t185 + mrSges(4,2) * t186 - Ifges(4,5) * t213 - Ifges(4,6) * t212 - Ifges(4,3) * t259 - pkin(3) * t277 - pkin(8) * t286 - t130 * t269 - t133 * t264 - t239 * t219 + t238 * t220;
t294 = mrSges(3,1) * t227 - mrSges(3,2) * t228 + Ifges(3,5) * t247 + Ifges(3,6) * t248 + Ifges(3,3) * qJDD(2) + pkin(2) * t131 + (t236 * t266 - t237 * t271) * qJD(1) - t278;
t140 = t145 * t269 + t146 * t264;
t246 = (-mrSges(3,1) * t271 + mrSges(3,2) * t266) * qJD(1);
t251 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t290;
t127 = m(3) * t227 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t247 + qJD(2) * t251 - t246 * t291 + t131;
t250 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t291;
t287 = t138 * t270 - t151 * t265;
t128 = m(3) * t228 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t248 - qJD(2) * t250 + t246 * t290 + t287;
t288 = -t127 * t266 + t128 * t271;
t218 = Ifges(4,5) * t239 + Ifges(4,6) * t238 + Ifges(4,3) * t260;
t121 = mrSges(4,2) * t214 - mrSges(4,3) * t185 + Ifges(4,1) * t213 + Ifges(4,4) * t212 + Ifges(4,5) * t259 - pkin(8) * t140 - t130 * t264 + t133 * t269 + t218 * t238 - t219 * t260;
t280 = -mrSges(6,1) * t160 + mrSges(6,2) * t161 - Ifges(6,5) * t171 - Ifges(6,6) * t170 - Ifges(6,3) * t206 - t180 * t197 + t181 * t196;
t274 = mrSges(5,1) * t164 - mrSges(5,2) * t165 + Ifges(5,5) * t189 + Ifges(5,6) * t188 + Ifges(5,3) * t211 + pkin(4) * t147 + t193 * t226 - t194 * t225 - t280;
t125 = -mrSges(4,1) * t214 + mrSges(4,3) * t186 + Ifges(4,4) * t213 + Ifges(4,2) * t212 + Ifges(4,6) * t259 - pkin(3) * t140 - t218 * t239 + t220 * t260 - t274;
t235 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t266 + Ifges(3,6) * t271) * qJD(1);
t240 = -t273 * pkin(6) + t283;
t279 = m(4) * t214 - t212 * mrSges(4,1) + mrSges(4,2) * t213 - t238 * t229 + t230 * t239 + t140;
t117 = -mrSges(3,1) * t240 + mrSges(3,3) * t228 + Ifges(3,4) * t247 + Ifges(3,2) * t248 + Ifges(3,6) * qJDD(2) - pkin(2) * t279 + pkin(7) * t287 + qJD(2) * t237 + t265 * t121 + t270 * t125 - t235 * t291;
t120 = mrSges(3,2) * t240 - mrSges(3,3) * t227 + Ifges(3,1) * t247 + Ifges(3,4) * t248 + Ifges(3,5) * qJDD(2) - pkin(7) * t131 - qJD(2) * t236 + t121 * t270 - t125 * t265 + t235 * t290;
t276 = -m(3) * t240 + t248 * mrSges(3,1) - mrSges(3,2) * t247 - t250 * t291 + t251 * t290 - t279;
t282 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + Ifges(2,3) * qJDD(1) + pkin(1) * t276 + pkin(6) * t288 + t117 * t271 + t120 * t266;
t134 = m(2) * t253 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t273 + t276;
t124 = t127 * t271 + t128 * t266;
t122 = m(2) * t254 - mrSges(2,1) * t273 - qJDD(1) * mrSges(2,2) + t288;
t118 = mrSges(2,1) * g(3) + mrSges(2,3) * t254 + t273 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t124 - t294;
t115 = -mrSges(2,2) * g(3) - mrSges(2,3) * t253 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t273 - pkin(6) * t124 - t117 * t266 + t120 * t271;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t272 * t115 - t267 * t118 - pkin(5) * (t122 * t267 + t134 * t272), t115, t120, t121, t133, t149; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t267 * t115 + t272 * t118 + pkin(5) * (t122 * t272 - t134 * t267), t118, t117, t125, t130, t148; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t282, t282, t294, -t278, t274, -t280;];
m_new = t1;
