% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:27
% EndTime: 2019-12-31 21:43:43
% DurationCPUTime: 7.32s
% Computational Cost: add. (116963->328), mult. (252760->407), div. (0->0), fcn. (188108->10), ass. (0->132)
t269 = sin(pkin(5));
t273 = sin(qJ(2));
t276 = cos(qJ(2));
t295 = qJD(1) * qJD(2);
t254 = (-qJDD(1) * t276 + t273 * t295) * t269;
t297 = qJD(1) * t269;
t252 = (-pkin(2) * t276 - pkin(8) * t273) * t297;
t270 = cos(pkin(5));
t265 = t270 * qJD(1) + qJD(2);
t263 = t265 ^ 2;
t264 = t270 * qJDD(1) + qJDD(2);
t296 = qJD(1) * t276;
t274 = sin(qJ(1));
t277 = cos(qJ(1));
t261 = t274 * g(1) - t277 * g(2);
t278 = qJD(1) ^ 2;
t308 = pkin(7) * t269;
t249 = qJDD(1) * pkin(1) + t278 * t308 + t261;
t262 = -t277 * g(1) - t274 * g(2);
t250 = -t278 * pkin(1) + qJDD(1) * t308 + t262;
t302 = t270 * t273;
t298 = t249 * t302 + t276 * t250;
t192 = -t263 * pkin(2) + t264 * pkin(8) + (-g(3) * t273 + t252 * t296) * t269 + t298;
t253 = (qJDD(1) * t273 + t276 * t295) * t269;
t307 = t270 * g(3);
t193 = t254 * pkin(2) - t253 * pkin(8) - t307 + (-t249 + (pkin(2) * t273 - pkin(8) * t276) * t265 * qJD(1)) * t269;
t272 = sin(qJ(3));
t309 = cos(qJ(3));
t176 = -t272 * t192 + t309 * t193;
t177 = t309 * t192 + t272 * t193;
t294 = t273 * t297;
t240 = -t309 * t265 + t272 * t294;
t241 = t272 * t265 + t309 * t294;
t293 = t269 * t296;
t259 = -qJD(3) + t293;
t204 = Ifges(4,4) * t241 - Ifges(4,2) * t240 - Ifges(4,6) * t259;
t215 = t241 * qJD(3) + t272 * t253 - t309 * t264;
t216 = -t240 * qJD(3) + t309 * t253 + t272 * t264;
t221 = -t240 * mrSges(5,2) - t241 * mrSges(5,3);
t227 = t240 * mrSges(5,1) + t259 * mrSges(5,3);
t246 = qJDD(3) + t254;
t219 = t240 * pkin(3) - t241 * qJ(4);
t258 = t259 ^ 2;
t174 = -t246 * pkin(3) - t258 * qJ(4) + t241 * t219 + qJDD(4) - t176;
t305 = t240 * t259;
t168 = (t240 * t241 - t246) * pkin(9) + (t216 - t305) * pkin(4) + t174;
t229 = t241 * pkin(4) + t259 * pkin(9);
t239 = t240 ^ 2;
t301 = t270 * t276;
t303 = t269 * t276;
t217 = -g(3) * t303 + t249 * t301 - t273 * t250;
t191 = -t264 * pkin(2) - t263 * pkin(8) + t252 * t294 - t217;
t310 = -2 * qJD(4);
t280 = (-t216 - t305) * qJ(4) + t191 + (-t259 * pkin(3) + t310) * t241;
t171 = -t239 * pkin(4) - t241 * t229 + (pkin(3) + pkin(9)) * t215 + t280;
t271 = sin(qJ(5));
t275 = cos(qJ(5));
t165 = t275 * t168 - t271 * t171;
t223 = t275 * t240 + t271 * t259;
t183 = t223 * qJD(5) + t271 * t215 + t275 * t246;
t224 = t271 * t240 - t275 * t259;
t194 = -t223 * mrSges(6,1) + t224 * mrSges(6,2);
t238 = qJD(5) + t241;
t198 = -t238 * mrSges(6,2) + t223 * mrSges(6,3);
t212 = qJDD(5) + t216;
t162 = m(6) * t165 + t212 * mrSges(6,1) - t183 * mrSges(6,3) - t224 * t194 + t238 * t198;
t166 = t271 * t168 + t275 * t171;
t182 = -t224 * qJD(5) + t275 * t215 - t271 * t246;
t199 = t238 * mrSges(6,1) - t224 * mrSges(6,3);
t163 = m(6) * t166 - t212 * mrSges(6,2) + t182 * mrSges(6,3) + t223 * t194 - t238 * t199;
t151 = t275 * t162 + t271 * t163;
t288 = -t258 * pkin(3) + t246 * qJ(4) - t240 * t219 + t177;
t170 = -t215 * pkin(4) - t239 * pkin(9) + (t310 - t229) * t259 + t288;
t184 = Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * t238;
t186 = Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t238;
t154 = -mrSges(6,1) * t170 + mrSges(6,3) * t166 + Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t212 - t224 * t184 + t238 * t186;
t185 = Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t238;
t155 = mrSges(6,2) * t170 - mrSges(6,3) * t165 + Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t212 + t223 * t184 - t238 * t185;
t172 = 0.2e1 * qJD(4) * t259 - t288;
t201 = -Ifges(5,5) * t259 - Ifges(5,6) * t241 + Ifges(5,3) * t240;
t284 = -mrSges(5,2) * t174 + mrSges(5,3) * t172 - Ifges(5,1) * t246 + Ifges(5,4) * t216 - Ifges(5,5) * t215 + pkin(9) * t151 + t271 * t154 - t275 * t155 + t241 * t201;
t167 = -m(6) * t170 + t182 * mrSges(6,1) - t183 * mrSges(6,2) + t223 * t198 - t224 * t199;
t228 = t241 * mrSges(5,1) - t259 * mrSges(5,2);
t285 = -m(5) * t172 + t246 * mrSges(5,3) - t259 * t228 - t167;
t289 = -m(5) * t174 - t216 * mrSges(5,1) - t241 * t221 - t151;
t203 = -Ifges(5,4) * t259 - Ifges(5,2) * t241 + Ifges(5,6) * t240;
t299 = -Ifges(4,1) * t241 + Ifges(4,4) * t240 + Ifges(4,5) * t259 + t203;
t311 = -t299 * t240 + mrSges(4,1) * t176 - mrSges(4,2) * t177 + Ifges(4,5) * t216 - Ifges(4,6) * t215 + Ifges(4,3) * t246 + pkin(3) * (-t246 * mrSges(5,2) + t259 * t227 + t289) + qJ(4) * (-t215 * mrSges(5,1) - t240 * t221 + t285) + t241 * t204 - t284;
t306 = Ifges(4,4) + Ifges(5,6);
t304 = t269 * t273;
t220 = t240 * mrSges(4,1) + t241 * mrSges(4,2);
t225 = t259 * mrSges(4,2) - t240 * mrSges(4,3);
t148 = m(4) * t176 - t216 * mrSges(4,3) - t241 * t220 + (-t225 + t227) * t259 + (mrSges(4,1) - mrSges(5,2)) * t246 + t289;
t226 = -t259 * mrSges(4,1) - t241 * mrSges(4,3);
t158 = m(4) * t177 - t246 * mrSges(4,2) + t259 * t226 + (-t220 - t221) * t240 + (-mrSges(4,3) - mrSges(5,1)) * t215 + t285;
t144 = t309 * t148 + t272 * t158;
t152 = -t271 * t162 + t275 * t163;
t205 = -Ifges(5,1) * t259 - Ifges(5,4) * t241 + Ifges(5,5) * t240;
t300 = -Ifges(4,5) * t241 + Ifges(4,6) * t240 + Ifges(4,3) * t259 - t205;
t218 = -g(3) * t304 + t298;
t247 = t265 * mrSges(3,1) - mrSges(3,3) * t294;
t251 = (-mrSges(3,1) * t276 + mrSges(3,2) * t273) * t297;
t292 = -t272 * t148 + t309 * t158;
t142 = m(3) * t218 - t264 * mrSges(3,2) - t254 * mrSges(3,3) - t265 * t247 + t251 * t293 + t292;
t248 = -t265 * mrSges(3,2) + mrSges(3,3) * t293;
t175 = t215 * pkin(3) + t280;
t291 = -m(5) * t175 + t215 * mrSges(5,2) + t240 * t227 - t152;
t281 = -m(4) * t191 - t215 * mrSges(4,1) - t240 * t225 + (-t226 + t228) * t241 + (-mrSges(4,2) + mrSges(5,3)) * t216 + t291;
t146 = m(3) * t217 + t264 * mrSges(3,1) - t253 * mrSges(3,3) + t265 * t248 - t251 * t294 + t281;
t139 = t276 * t142 - t273 * t146;
t234 = -t269 * t249 - t307;
t143 = m(3) * t234 + t254 * mrSges(3,1) + t253 * mrSges(3,2) + (t247 * t273 - t248 * t276) * t297 + t144;
t134 = t142 * t302 - t269 * t143 + t146 * t301;
t149 = -t216 * mrSges(5,3) - t241 * t228 - t291;
t283 = -mrSges(5,1) * t172 + mrSges(5,2) * t175 - pkin(4) * t167 - pkin(9) * t152 - t275 * t154 - t271 * t155;
t135 = -mrSges(4,1) * t191 + mrSges(4,3) * t177 - pkin(3) * t149 + t299 * t259 + (Ifges(4,6) - Ifges(5,5)) * t246 + t300 * t241 + t306 * t216 + (-Ifges(4,2) - Ifges(5,3)) * t215 + t283;
t286 = mrSges(6,1) * t165 - mrSges(6,2) * t166 + Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t212 + t224 * t185 - t223 * t186;
t282 = mrSges(5,1) * t174 - mrSges(5,3) * t175 + pkin(4) * t151 + t286;
t136 = t282 + mrSges(4,2) * t191 - mrSges(4,3) * t176 - qJ(4) * t149 - t306 * t215 + (Ifges(4,1) + Ifges(5,2)) * t216 + t300 * t240 + (Ifges(4,5) - Ifges(5,4)) * t246 + (t204 - t201) * t259;
t232 = Ifges(3,6) * t265 + (Ifges(3,4) * t273 + Ifges(3,2) * t276) * t297;
t233 = Ifges(3,5) * t265 + (Ifges(3,1) * t273 + Ifges(3,4) * t276) * t297;
t126 = Ifges(3,5) * t253 - Ifges(3,6) * t254 + Ifges(3,3) * t264 + mrSges(3,1) * t217 - mrSges(3,2) * t218 + t272 * t136 + t309 * t135 + pkin(2) * t281 + pkin(8) * t292 + (t232 * t273 - t233 * t276) * t297;
t231 = Ifges(3,3) * t265 + (Ifges(3,5) * t273 + Ifges(3,6) * t276) * t297;
t128 = mrSges(3,2) * t234 - mrSges(3,3) * t217 + Ifges(3,1) * t253 - Ifges(3,4) * t254 + Ifges(3,5) * t264 - pkin(8) * t144 - t272 * t135 + t309 * t136 + t231 * t293 - t265 * t232;
t130 = -mrSges(3,1) * t234 + mrSges(3,3) * t218 + Ifges(3,4) * t253 - Ifges(3,2) * t254 + Ifges(3,6) * t264 - pkin(2) * t144 - t231 * t294 + t265 * t233 - t311;
t287 = mrSges(2,1) * t261 - mrSges(2,2) * t262 + Ifges(2,3) * qJDD(1) + pkin(1) * t134 + t270 * t126 + t128 * t304 + t130 * t303 + t139 * t308;
t137 = m(2) * t262 - t278 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t139;
t133 = t270 * t143 + (t142 * t273 + t146 * t276) * t269;
t131 = m(2) * t261 + qJDD(1) * mrSges(2,1) - t278 * mrSges(2,2) + t134;
t124 = -mrSges(2,2) * g(3) - mrSges(2,3) * t261 + Ifges(2,5) * qJDD(1) - t278 * Ifges(2,6) + t276 * t128 - t273 * t130 + (-t133 * t269 - t134 * t270) * pkin(7);
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t262 + t278 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t133 - t269 * t126 + (pkin(7) * t139 + t128 * t273 + t130 * t276) * t270;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t277 * t124 - t274 * t123 - pkin(6) * (t277 * t131 + t274 * t137), t124, t128, t136, -t240 * t203 - t284, t155; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t274 * t124 + t277 * t123 + pkin(6) * (-t274 * t131 + t277 * t137), t123, t130, t135, Ifges(5,4) * t246 - Ifges(5,2) * t216 + Ifges(5,6) * t215 + t259 * t201 + t240 * t205 - t282, t154; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t287, t287, t126, t311, Ifges(5,5) * t246 - Ifges(5,6) * t216 + Ifges(5,3) * t215 - t259 * t203 + t241 * t205 - t283, t286;];
m_new = t1;
