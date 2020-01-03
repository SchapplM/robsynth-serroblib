% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:40
% EndTime: 2019-12-31 17:47:44
% DurationCPUTime: 3.01s
% Computational Cost: add. (26418->266), mult. (67823->347), div. (0->0), fcn. (38671->8), ass. (0->131)
t315 = -2 * qJD(4);
t245 = sin(qJ(1));
t247 = cos(qJ(1));
t222 = -g(1) * t247 - g(2) * t245;
t248 = qJD(1) ^ 2;
t314 = -pkin(1) * t248 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t222;
t241 = sin(pkin(7));
t290 = qJD(1) * t241;
t223 = -0.2e1 * qJD(3) * t290;
t237 = t241 ^ 2;
t243 = cos(pkin(7));
t238 = t243 ^ 2;
t221 = g(1) * t245 - t247 * g(2);
t276 = qJDD(2) - t221;
t301 = qJ(3) * t241;
t171 = t223 + (-qJ(2) + (-t237 - t238) * pkin(3)) * t248 + (-t301 - pkin(1) + (-pkin(2) - qJ(4)) * t243) * qJDD(1) + t276;
t289 = qJD(1) * t243;
t313 = t289 * t315 + t171;
t188 = -g(3) * t241 + t314 * t243;
t271 = -pkin(2) * t243 - t301;
t207 = t271 * qJD(1);
t178 = -t207 * t289 - t188;
t298 = t238 * t248;
t299 = t237 * t248;
t312 = t298 + t299;
t187 = -t243 * g(3) - t314 * t241;
t303 = mrSges(4,2) * t243;
t208 = (-mrSges(4,3) * t241 + t303) * qJD(1);
t272 = -Ifges(4,6) * t241 - Ifges(4,3) * t243;
t210 = t272 * qJD(1);
t177 = t207 * t290 + qJDD(3) - t187;
t310 = t243 * t248;
t167 = (pkin(3) * qJDD(1) - qJ(4) * t310) * t241 + t177;
t240 = sin(pkin(8));
t242 = cos(pkin(8));
t164 = t240 * t167 + t313 * t242;
t198 = (pkin(4) * t242 + pkin(6) * t240) * t289;
t285 = t242 * t289;
t287 = qJDD(1) * t241;
t161 = -pkin(4) * t299 + pkin(6) * t287 - t198 * t285 + t164;
t286 = qJDD(1) * t243;
t169 = pkin(3) * t286 - qJ(4) * t298 + qJDD(4) - t178;
t296 = t241 * t248;
t165 = ((qJDD(1) * t240 + t242 * t296) * pkin(6) + (qJDD(1) * t242 - t240 * t296) * pkin(4)) * t243 + t169;
t244 = sin(qJ(5));
t246 = cos(qJ(5));
t158 = -t161 * t244 + t165 * t246;
t297 = t240 * t243;
t264 = t241 * t246 + t244 * t297;
t195 = t264 * qJD(1);
t265 = t241 * t244 - t246 * t297;
t196 = t265 * qJD(1);
t179 = -mrSges(6,1) * t195 + mrSges(6,2) * t196;
t182 = t195 * qJD(5) + t265 * qJDD(1);
t218 = qJD(5) + t285;
t184 = -mrSges(6,2) * t218 + mrSges(6,3) * t195;
t281 = t242 * t286;
t217 = qJDD(5) + t281;
t154 = m(6) * t158 + mrSges(6,1) * t217 - t182 * mrSges(6,3) - t179 * t196 + t184 * t218;
t159 = t161 * t246 + t165 * t244;
t181 = -t196 * qJD(5) + t264 * qJDD(1);
t185 = mrSges(6,1) * t218 - mrSges(6,3) * t196;
t155 = m(6) * t159 - mrSges(6,2) * t217 + t181 * mrSges(6,3) + t179 * t195 - t185 * t218;
t143 = t246 * t154 + t244 * t155;
t295 = t242 * t167;
t160 = -pkin(4) * t287 - pkin(6) * t299 - t295 + (t171 + (t315 - t198) * t289) * t240;
t172 = Ifges(6,5) * t196 + Ifges(6,6) * t195 + Ifges(6,3) * t218;
t174 = Ifges(6,1) * t196 + Ifges(6,4) * t195 + Ifges(6,5) * t218;
t147 = -mrSges(6,1) * t160 + mrSges(6,3) * t159 + Ifges(6,4) * t182 + Ifges(6,2) * t181 + Ifges(6,6) * t217 - t172 * t196 + t174 * t218;
t173 = Ifges(6,4) * t196 + Ifges(6,2) * t195 + Ifges(6,6) * t218;
t148 = mrSges(6,2) * t160 - mrSges(6,3) * t158 + Ifges(6,1) * t182 + Ifges(6,4) * t181 + Ifges(6,5) * t217 + t172 * t195 - t173 * t218;
t163 = -t313 * t240 + t295;
t273 = -Ifges(5,5) * t240 - Ifges(5,6) * t242;
t189 = (Ifges(5,3) * t241 + t273 * t243) * qJD(1);
t256 = Ifges(5,6) * t241 + (-Ifges(5,4) * t240 - Ifges(5,2) * t242) * t243;
t190 = t256 * qJD(1);
t257 = Ifges(5,5) * t241 + (-Ifges(5,1) * t240 - Ifges(5,4) * t242) * t243;
t294 = t242 * t243;
t127 = mrSges(5,2) * t169 - mrSges(5,3) * t163 - pkin(6) * t143 - t147 * t244 + t148 * t246 + (-t189 * t294 - t190 * t241) * qJD(1) + t257 * qJDD(1);
t191 = t257 * qJD(1);
t250 = mrSges(6,1) * t158 - mrSges(6,2) * t159 + Ifges(6,5) * t182 + Ifges(6,6) * t181 + Ifges(6,3) * t217 + t196 * t173 - t195 * t174;
t128 = -mrSges(5,1) * t169 + mrSges(5,3) * t164 - pkin(4) * t143 + (t189 * t297 + t191 * t241) * qJD(1) + t256 * qJDD(1) - t250;
t302 = mrSges(5,2) * t240;
t197 = (mrSges(5,1) * t242 - t302) * t289;
t267 = mrSges(5,1) * t241 + mrSges(5,3) * t297;
t204 = t267 * qJD(1);
t266 = -mrSges(5,2) * t241 - mrSges(5,3) * t294;
t278 = -t244 * t154 + t246 * t155;
t140 = m(5) * t164 + t266 * qJDD(1) + (-t197 * t294 - t204 * t241) * qJD(1) + t278;
t205 = t266 * qJD(1);
t259 = -m(6) * t160 + t181 * mrSges(6,1) - t182 * mrSges(6,2) + t195 * t184 - t185 * t196;
t150 = m(5) * t163 + t267 * qJDD(1) + (t197 * t297 + t205 * t241) * qJD(1) + t259;
t134 = t240 * t140 + t242 * t150;
t254 = mrSges(4,2) * t177 - mrSges(4,3) * t178 - qJ(4) * t134 + t242 * t127 - t240 * t128 + (-Ifges(4,2) * t241 - Ifges(4,6) * t243) * t310;
t270 = -m(5) * t169 - mrSges(5,1) * t281 - t205 * t285 - t143;
t258 = -m(4) * t178 + mrSges(4,1) * t286 + t208 * t289 - t270;
t261 = -m(4) * t177 - t134;
t309 = (mrSges(5,2) * qJDD(1) + qJD(1) * t204) * t297;
t311 = ((t210 - (Ifges(3,4) * t241 + Ifges(3,2) * t243) * qJD(1)) * t241 + (Ifges(3,1) * t241 + Ifges(3,4) * t243) * t289) * qJD(1) - mrSges(3,1) * t187 + mrSges(3,2) * t188 - pkin(2) * ((-qJDD(1) * mrSges(4,1) - qJD(1) * t208) * t241 + t261) - qJ(3) * (t258 - t309) - t254;
t307 = -(Ifges(4,4) - Ifges(3,5)) * t241 - (Ifges(4,5) - Ifges(3,6)) * t243;
t304 = mrSges(3,2) * t241;
t300 = t190 * t240;
t135 = t242 * t140 - t240 * t150;
t209 = (-mrSges(3,1) * t243 + t304) * qJD(1);
t129 = m(3) * t187 + ((-mrSges(4,1) - mrSges(3,3)) * qJDD(1) + (-t208 - t209) * qJD(1)) * t241 + t261;
t137 = m(3) * t188 + ((mrSges(3,3) - t302) * qJDD(1) + (-t204 * t240 + t209) * qJD(1)) * t243 + t258;
t279 = -t129 * t241 + t243 * t137;
t274 = -Ifges(4,4) * t241 - Ifges(4,5) * t243;
t213 = t274 * qJD(1);
t277 = -t213 + t300;
t275 = -t303 - t304;
t268 = -Ifges(4,6) - t273;
t263 = -qJ(2) * t248 + t276;
t183 = t223 + (-pkin(1) + t271) * qJDD(1) + t263;
t262 = -m(4) * t183 + t312 * mrSges(4,1) + mrSges(4,3) * t287 - t135;
t133 = mrSges(4,2) * t286 - t262;
t203 = -qJDD(1) * pkin(1) + t263;
t212 = (Ifges(3,5) * t241 + Ifges(3,6) * t243) * qJD(1);
t252 = mrSges(4,1) * t178 - mrSges(4,2) * t183 + pkin(3) * (t270 + t309) + qJ(4) * t135 + t127 * t240 + t128 * t242;
t120 = -mrSges(3,1) * t203 + mrSges(3,3) * t188 - pkin(2) * t133 + (-t212 - t213) * t290 + ((Ifges(3,2) + Ifges(4,3)) * t243 + (Ifges(3,4) + Ifges(4,6)) * t241) * qJDD(1) - t252;
t255 = mrSges(5,1) * t163 - mrSges(5,2) * t164 + Ifges(5,3) * t287 + pkin(4) * t259 + pkin(6) * t278 + t246 * t147 + t244 * t148 + t191 * t285;
t251 = -mrSges(4,1) * t177 + mrSges(4,3) * t183 - pkin(3) * t134 - t255;
t122 = -t251 + ((t212 - t277) * qJD(1) + (Ifges(3,4) - t268) * qJDD(1)) * t243 + mrSges(3,2) * t203 - mrSges(3,3) * t187 - qJ(3) * t133 + (Ifges(3,1) + Ifges(4,2)) * t287;
t253 = -m(3) * t203 + mrSges(3,1) * t286 + t312 * mrSges(3,3) + t262;
t260 = -mrSges(2,2) * t222 + qJ(2) * t279 + t243 * t120 + t241 * t122 + pkin(1) * (t275 * qJDD(1) + t253) + mrSges(2,1) * t221 + Ifges(2,3) * qJDD(1);
t130 = m(2) * t221 - mrSges(2,2) * t248 + (mrSges(2,1) + t275) * qJDD(1) + t253;
t125 = t129 * t243 + t137 * t241;
t123 = m(2) * t222 - mrSges(2,1) * t248 - qJDD(1) * mrSges(2,2) + t279;
t118 = (Ifges(2,6) - t307) * qJDD(1) + t248 * Ifges(2,5) + mrSges(2,3) * t222 + mrSges(2,1) * g(3) - pkin(1) * t125 + t311;
t117 = -mrSges(2,2) * g(3) - mrSges(2,3) * t221 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t248 - qJ(2) * t125 - t120 * t241 + t122 * t243;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247 * t117 - t245 * t118 - pkin(5) * (t123 * t245 + t130 * t247), t117, t122, t274 * qJDD(1) - t210 * t290 + t254, t127, t148; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t245 * t117 + t247 * t118 + pkin(5) * (t123 * t247 - t130 * t245), t118, t120, t251 + (t277 * qJD(1) + t268 * qJDD(1)) * t243 - Ifges(4,2) * t287, t128, t147; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t260, t260, t307 * qJDD(1) - t311, t272 * qJDD(1) + t213 * t290 + t252, (-qJD(1) * t300 + t273 * qJDD(1)) * t243 + t255, t250;];
m_new = t1;
