% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:31
% EndTime: 2019-12-31 21:32:41
% DurationCPUTime: 4.52s
% Computational Cost: add. (54992->308), mult. (107923->372), div. (0->0), fcn. (68259->8), ass. (0->117)
t267 = sin(qJ(3));
t268 = sin(qJ(2));
t293 = qJD(1) * t268;
t299 = cos(qJ(3));
t245 = -t299 * qJD(2) + t267 * t293;
t271 = cos(qJ(2));
t291 = qJD(1) * qJD(2);
t289 = t271 * t291;
t249 = t268 * qJDD(1) + t289;
t211 = -t245 * qJD(3) + t267 * qJDD(2) + t299 * t249;
t269 = sin(qJ(1));
t272 = cos(qJ(1));
t256 = -t272 * g(1) - t269 * g(2);
t274 = qJD(1) ^ 2;
t234 = -t274 * pkin(1) + qJDD(1) * pkin(6) + t256;
t223 = -t271 * g(3) - t268 * t234;
t248 = (-pkin(2) * t271 - pkin(7) * t268) * qJD(1);
t273 = qJD(2) ^ 2;
t283 = qJDD(2) * pkin(2) + t273 * pkin(7) - t248 * t293 + t223;
t292 = t271 * qJD(1);
t259 = qJD(3) - t292;
t297 = t245 * t259;
t303 = (-t211 + t297) * qJ(4) - t283;
t255 = t269 * g(1) - t272 * g(2);
t233 = -qJDD(1) * pkin(1) - t274 * pkin(6) - t255;
t290 = t268 * t291;
t250 = t271 * qJDD(1) - t290;
t185 = (-t249 - t289) * pkin(7) + (-t250 + t290) * pkin(2) + t233;
t224 = -t268 * g(3) + t271 * t234;
t190 = -t273 * pkin(2) + qJDD(2) * pkin(7) + t248 * t292 + t224;
t176 = t299 * t185 - t267 * t190;
t177 = t267 * t185 + t299 * t190;
t246 = t267 * qJD(2) + t299 * t293;
t196 = Ifges(5,5) * t246 + Ifges(5,6) * t259 + Ifges(5,3) * t245;
t199 = Ifges(4,4) * t246 - Ifges(4,2) * t245 + Ifges(4,6) * t259;
t201 = Ifges(4,1) * t246 - Ifges(4,4) * t245 + Ifges(4,5) * t259;
t210 = t246 * qJD(3) - t299 * qJDD(2) + t267 * t249;
t217 = t245 * mrSges(5,1) - t246 * mrSges(5,3);
t244 = qJDD(3) - t250;
t216 = t245 * pkin(3) - t246 * qJ(4);
t258 = t259 ^ 2;
t167 = -t244 * pkin(3) - t258 * qJ(4) + t246 * t216 + qJDD(4) - t176;
t158 = (-t211 - t297) * pkin(8) + (t245 * t246 - t244) * pkin(4) + t167;
t300 = 2 * qJD(4);
t164 = -t258 * pkin(3) + t244 * qJ(4) - t245 * t216 + t259 * t300 + t177;
t225 = -t259 * pkin(4) - t246 * pkin(8);
t243 = t245 ^ 2;
t160 = -t243 * pkin(4) + t210 * pkin(8) + t259 * t225 + t164;
t266 = sin(qJ(5));
t270 = cos(qJ(5));
t156 = t270 * t158 - t266 * t160;
t212 = t270 * t245 - t266 * t246;
t175 = t212 * qJD(5) + t266 * t210 + t270 * t211;
t213 = t266 * t245 + t270 * t246;
t183 = -t212 * mrSges(6,1) + t213 * mrSges(6,2);
t257 = qJD(5) - t259;
t192 = -t257 * mrSges(6,2) + t212 * mrSges(6,3);
t240 = qJDD(5) - t244;
t151 = m(6) * t156 + t240 * mrSges(6,1) - t175 * mrSges(6,3) - t213 * t183 + t257 * t192;
t157 = t266 * t158 + t270 * t160;
t174 = -t213 * qJD(5) + t270 * t210 - t266 * t211;
t193 = t257 * mrSges(6,1) - t213 * mrSges(6,3);
t152 = m(6) * t157 - t240 * mrSges(6,2) + t174 * mrSges(6,3) + t212 * t183 - t257 * t193;
t142 = t270 * t151 + t266 * t152;
t200 = Ifges(5,1) * t246 + Ifges(5,4) * t259 + Ifges(5,5) * t245;
t179 = Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t257;
t180 = Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t257;
t284 = mrSges(6,1) * t156 - mrSges(6,2) * t157 + Ifges(6,5) * t175 + Ifges(6,6) * t174 + Ifges(6,3) * t240 + t213 * t179 - t212 * t180;
t277 = mrSges(5,1) * t167 - mrSges(5,3) * t164 - Ifges(5,4) * t211 - Ifges(5,2) * t244 - Ifges(5,6) * t210 + pkin(4) * t142 - t245 * t200 + t284;
t222 = -t245 * mrSges(5,2) + t259 * mrSges(5,3);
t281 = -m(5) * t167 + t244 * mrSges(5,1) + t259 * t222 - t142;
t143 = -t266 * t151 + t270 * t152;
t221 = -t259 * mrSges(5,1) + t246 * mrSges(5,2);
t285 = m(5) * t164 + t244 * mrSges(5,3) + t259 * t221 + t143;
t302 = (t199 - t196) * t246 + mrSges(4,1) * t176 - mrSges(4,2) * t177 + Ifges(4,5) * t211 - Ifges(4,6) * t210 + Ifges(4,3) * t244 + pkin(3) * (-t211 * mrSges(5,2) - t246 * t217 + t281) + qJ(4) * (-t210 * mrSges(5,2) - t245 * t217 + t285) + t245 * t201 - t277;
t161 = -t243 * pkin(8) + (-pkin(3) - pkin(4)) * t210 + (-pkin(3) * t259 + t225 + t300) * t246 - t303;
t153 = -m(6) * t161 + t174 * mrSges(6,1) - t175 * mrSges(6,2) + t212 * t192 - t213 * t193;
t166 = -0.2e1 * qJD(4) * t246 + (t246 * t259 + t210) * pkin(3) + t303;
t149 = m(5) * t166 + t210 * mrSges(5,1) - t211 * mrSges(5,3) - t246 * t221 + t245 * t222 + t153;
t178 = Ifges(6,5) * t213 + Ifges(6,6) * t212 + Ifges(6,3) * t257;
t145 = -mrSges(6,1) * t161 + mrSges(6,3) * t157 + Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t240 - t213 * t178 + t257 * t180;
t146 = mrSges(6,2) * t161 - mrSges(6,3) * t156 + Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t240 + t212 * t178 - t257 * t179;
t279 = -mrSges(5,1) * t166 + mrSges(5,2) * t164 - pkin(4) * t153 - pkin(8) * t143 - t270 * t145 - t266 * t146;
t198 = Ifges(5,4) * t246 + Ifges(5,2) * t259 + Ifges(5,6) * t245;
t296 = -Ifges(4,5) * t246 + Ifges(4,6) * t245 - Ifges(4,3) * t259 - t198;
t126 = mrSges(4,1) * t283 + mrSges(4,3) * t177 - pkin(3) * t149 + (t201 + t200) * t259 + t296 * t246 + (Ifges(4,6) - Ifges(5,6)) * t244 + (Ifges(4,4) - Ifges(5,5)) * t211 + (-Ifges(4,2) - Ifges(5,3)) * t210 + t279;
t280 = mrSges(5,2) * t167 - mrSges(5,3) * t166 + Ifges(5,1) * t211 + Ifges(5,4) * t244 + Ifges(5,5) * t210 - pkin(8) * t142 - t266 * t145 + t270 * t146 + t259 * t196;
t127 = -mrSges(4,2) * t283 - mrSges(4,3) * t176 + Ifges(4,1) * t211 - Ifges(4,4) * t210 + Ifges(4,5) * t244 - qJ(4) * t149 - t259 * t199 + t296 * t245 + t280;
t220 = t259 * mrSges(4,1) - t246 * mrSges(4,3);
t294 = -t245 * mrSges(4,1) - t246 * mrSges(4,2) - t217;
t298 = -mrSges(4,3) - mrSges(5,2);
t138 = m(4) * t177 - t244 * mrSges(4,2) + t298 * t210 - t259 * t220 + t294 * t245 + t285;
t219 = -t259 * mrSges(4,2) - t245 * mrSges(4,3);
t139 = m(4) * t176 + t244 * mrSges(4,1) + t298 * t211 + t259 * t219 + t294 * t246 + t281;
t136 = t299 * t138 - t267 * t139;
t148 = m(4) * t283 - t210 * mrSges(4,1) - t211 * mrSges(4,2) - t245 * t219 - t246 * t220 - t149;
t231 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t268 + Ifges(3,2) * t271) * qJD(1);
t232 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t268 + Ifges(3,4) * t271) * qJD(1);
t301 = mrSges(3,1) * t223 - mrSges(3,2) * t224 + Ifges(3,5) * t249 + Ifges(3,6) * t250 + Ifges(3,3) * qJDD(2) + pkin(2) * t148 + pkin(7) * t136 + (t231 * t268 - t232 * t271) * qJD(1) + t299 * t126 + t267 * t127;
t247 = (-mrSges(3,1) * t271 + mrSges(3,2) * t268) * qJD(1);
t252 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t293;
t134 = m(3) * t224 - qJDD(2) * mrSges(3,2) + t250 * mrSges(3,3) - qJD(2) * t252 + t247 * t292 + t136;
t253 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t292;
t147 = m(3) * t223 + qJDD(2) * mrSges(3,1) - t249 * mrSges(3,3) + qJD(2) * t253 - t247 * t293 + t148;
t288 = t271 * t134 - t268 * t147;
t135 = t267 * t138 + t299 * t139;
t230 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t268 + Ifges(3,6) * t271) * qJD(1);
t123 = mrSges(3,2) * t233 - mrSges(3,3) * t223 + Ifges(3,1) * t249 + Ifges(3,4) * t250 + Ifges(3,5) * qJDD(2) - pkin(7) * t135 - qJD(2) * t231 - t267 * t126 + t299 * t127 + t230 * t292;
t125 = -mrSges(3,1) * t233 + mrSges(3,3) * t224 + Ifges(3,4) * t249 + Ifges(3,2) * t250 + Ifges(3,6) * qJDD(2) - pkin(2) * t135 + qJD(2) * t232 - t230 * t293 - t302;
t278 = -m(3) * t233 + t250 * mrSges(3,1) - t249 * mrSges(3,2) - t252 * t293 + t253 * t292 - t135;
t282 = mrSges(2,1) * t255 - mrSges(2,2) * t256 + Ifges(2,3) * qJDD(1) + pkin(1) * t278 + pkin(6) * t288 + t268 * t123 + t271 * t125;
t131 = m(2) * t255 + qJDD(1) * mrSges(2,1) - t274 * mrSges(2,2) + t278;
t130 = t268 * t134 + t271 * t147;
t128 = m(2) * t256 - t274 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t288;
t121 = mrSges(2,1) * g(3) + mrSges(2,3) * t256 + t274 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t301;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t255 + Ifges(2,5) * qJDD(1) - t274 * Ifges(2,6) - pkin(6) * t130 + t271 * t123 - t268 * t125;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t272 * t120 - t269 * t121 - pkin(5) * (t269 * t128 + t272 * t131), t120, t123, t127, -t245 * t198 + t280, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t269 * t120 + t272 * t121 + pkin(5) * (t272 * t128 - t269 * t131), t121, t125, t126, -t246 * t196 - t277, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t282, t282, t301, t302, Ifges(5,5) * t211 + Ifges(5,6) * t244 + Ifges(5,3) * t210 + t246 * t198 - t259 * t200 - t279, t284;];
m_new = t1;
