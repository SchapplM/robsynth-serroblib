% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:53
% EndTime: 2019-12-05 18:45:09
% DurationCPUTime: 6.50s
% Computational Cost: add. (95016->306), mult. (209749->378), div. (0->0), fcn. (145497->8), ass. (0->115)
t272 = sin(qJ(3));
t273 = sin(qJ(2));
t276 = cos(qJ(3));
t277 = cos(qJ(2));
t243 = (t272 * t277 + t273 * t276) * qJD(1);
t299 = qJD(1) * qJD(2);
t250 = qJDD(1) * t273 + t277 * t299;
t251 = qJDD(1) * t277 - t273 * t299;
t215 = -qJD(3) * t243 - t250 * t272 + t251 * t276;
t242 = (-t272 * t273 + t276 * t277) * qJD(1);
t216 = qJD(3) * t242 + t250 * t276 + t251 * t272;
t271 = sin(qJ(4));
t275 = cos(qJ(4));
t229 = t242 * t275 - t243 * t271;
t184 = qJD(4) * t229 + t215 * t271 + t216 * t275;
t230 = t242 * t271 + t243 * t275;
t203 = -mrSges(6,1) * t229 + mrSges(6,2) * t230;
t274 = sin(qJ(1));
t278 = cos(qJ(1));
t257 = -g(1) * t278 - g(2) * t274;
t279 = qJD(1) ^ 2;
t245 = -pkin(1) * t279 + qJDD(1) * pkin(6) + t257;
t303 = t273 * t245;
t304 = pkin(2) * t279;
t209 = qJDD(2) * pkin(2) - t250 * pkin(7) - t303 + (pkin(7) * t299 + t273 * t304 - g(3)) * t277;
t233 = -g(3) * t273 + t245 * t277;
t301 = qJD(1) * t273;
t255 = qJD(2) * pkin(2) - pkin(7) * t301;
t270 = t277 ^ 2;
t210 = pkin(7) * t251 - qJD(2) * t255 - t270 * t304 + t233;
t190 = t209 * t276 - t272 * t210;
t267 = qJDD(2) + qJDD(3);
t268 = qJD(2) + qJD(3);
t165 = (t242 * t268 - t216) * pkin(8) + (t242 * t243 + t267) * pkin(3) + t190;
t191 = t209 * t272 + t210 * t276;
t236 = pkin(3) * t268 - pkin(8) * t243;
t238 = t242 ^ 2;
t167 = -pkin(3) * t238 + pkin(8) * t215 - t236 * t268 + t191;
t158 = t165 * t275 - t271 * t167;
t264 = qJDD(4) + t267;
t265 = qJD(4) + t268;
t153 = -0.2e1 * qJD(5) * t230 + (t229 * t265 - t184) * qJ(5) + (t229 * t230 + t264) * pkin(4) + t158;
t218 = -mrSges(6,2) * t265 + mrSges(6,3) * t229;
t298 = m(6) * t153 + t264 * mrSges(6,1) + t218 * t265;
t150 = -t184 * mrSges(6,3) - t230 * t203 + t298;
t159 = t165 * t271 + t167 * t275;
t183 = -qJD(4) * t230 + t215 * t275 - t216 * t271;
t196 = Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t265;
t197 = Ifges(6,1) * t230 + Ifges(6,4) * t229 + Ifges(6,5) * t265;
t198 = Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t265;
t220 = pkin(4) * t265 - qJ(5) * t230;
t228 = t229 ^ 2;
t156 = -pkin(4) * t228 + qJ(5) * t183 + 0.2e1 * qJD(5) * t229 - t220 * t265 + t159;
t195 = Ifges(6,4) * t230 + Ifges(6,2) * t229 + Ifges(6,6) * t265;
t288 = -mrSges(6,1) * t153 + mrSges(6,2) * t156 - Ifges(6,5) * t184 - Ifges(6,6) * t183 - Ifges(6,3) * t264 - t195 * t230;
t308 = mrSges(5,1) * t158 - mrSges(5,2) * t159 + Ifges(5,5) * t184 + Ifges(5,6) * t183 + Ifges(5,3) * t264 + pkin(4) * t150 + t230 * t196 - t288 - (t197 + t198) * t229;
t204 = -mrSges(5,1) * t229 + mrSges(5,2) * t230;
t219 = -mrSges(5,2) * t265 + mrSges(5,3) * t229;
t143 = m(5) * t158 + t264 * mrSges(5,1) + t265 * t219 + (-t203 - t204) * t230 + (-mrSges(5,3) - mrSges(6,3)) * t184 + t298;
t221 = mrSges(6,1) * t265 - mrSges(6,3) * t230;
t222 = mrSges(5,1) * t265 - mrSges(5,3) * t230;
t297 = m(6) * t156 + mrSges(6,3) * t183 + t203 * t229;
t147 = m(5) * t159 + t183 * mrSges(5,3) + t229 * t204 + (-t221 - t222) * t265 + (-mrSges(5,2) - mrSges(6,2)) * t264 + t297;
t141 = t143 * t275 + t147 * t271;
t224 = Ifges(4,4) * t243 + Ifges(4,2) * t242 + Ifges(4,6) * t268;
t225 = Ifges(4,1) * t243 + Ifges(4,4) * t242 + Ifges(4,5) * t268;
t307 = mrSges(4,1) * t190 - mrSges(4,2) * t191 + Ifges(4,5) * t216 + Ifges(4,6) * t215 + Ifges(4,3) * t267 + pkin(3) * t141 + t224 * t243 - t242 * t225 + t308;
t231 = -mrSges(4,1) * t242 + mrSges(4,2) * t243;
t234 = -mrSges(4,2) * t268 + mrSges(4,3) * t242;
t137 = m(4) * t190 + mrSges(4,1) * t267 - mrSges(4,3) * t216 - t231 * t243 + t234 * t268 + t141;
t235 = mrSges(4,1) * t268 - mrSges(4,3) * t243;
t293 = -t143 * t271 + t147 * t275;
t138 = m(4) * t191 - mrSges(4,2) * t267 + mrSges(4,3) * t215 + t231 * t242 - t235 * t268 + t293;
t132 = t137 * t276 + t138 * t272;
t232 = -t277 * g(3) - t303;
t240 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t273 + Ifges(3,2) * t277) * qJD(1);
t241 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t273 + Ifges(3,4) * t277) * qJD(1);
t306 = mrSges(3,1) * t232 - mrSges(3,2) * t233 + Ifges(3,5) * t250 + Ifges(3,6) * t251 + Ifges(3,3) * qJDD(2) + pkin(2) * t132 + (t240 * t273 - t241 * t277) * qJD(1) + t307;
t300 = qJD(1) * t277;
t256 = t274 * g(1) - t278 * g(2);
t249 = (-mrSges(3,1) * t277 + mrSges(3,2) * t273) * qJD(1);
t254 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t300;
t130 = m(3) * t232 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t250 + qJD(2) * t254 - t249 * t301 + t132;
t253 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t301;
t294 = -t137 * t272 + t138 * t276;
t131 = m(3) * t233 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t251 - qJD(2) * t253 + t249 * t300 + t294;
t295 = -t130 * t273 + t131 * t277;
t290 = -qJDD(1) * pkin(1) - t256;
t217 = -t251 * pkin(2) + t255 * t301 + (-pkin(7) * t270 - pkin(6)) * t279 + t290;
t169 = -t215 * pkin(3) - t238 * pkin(8) + t236 * t243 + t217;
t162 = -t183 * pkin(4) - t228 * qJ(5) + t230 * t220 + qJDD(5) + t169;
t292 = m(6) * t162 - mrSges(6,1) * t183 + mrSges(6,2) * t184 - t218 * t229 + t221 * t230;
t289 = -mrSges(6,1) * t162 + mrSges(6,3) * t156 + Ifges(6,4) * t184 + Ifges(6,2) * t183 + Ifges(6,6) * t264 + t197 * t265;
t193 = Ifges(6,5) * t230 + Ifges(6,6) * t229 + Ifges(6,3) * t265;
t287 = mrSges(6,2) * t162 - mrSges(6,3) * t153 + Ifges(6,1) * t184 + Ifges(6,4) * t183 + Ifges(6,5) * t264 + t193 * t229;
t194 = Ifges(5,5) * t230 + Ifges(5,6) * t229 + Ifges(5,3) * t265;
t133 = Ifges(5,4) * t184 + Ifges(5,2) * t183 + Ifges(5,6) * t264 + t265 * t198 - mrSges(5,1) * t169 + mrSges(5,3) * t159 - pkin(4) * t292 + qJ(5) * (-t264 * mrSges(6,2) - t265 * t221 + t297) + (-t194 - t193) * t230 + t289;
t139 = mrSges(5,2) * t169 - mrSges(5,3) * t158 + Ifges(5,1) * t184 + Ifges(5,4) * t183 + Ifges(5,5) * t264 - qJ(5) * t150 + t229 * t194 + (-t195 - t196) * t265 + t287;
t223 = Ifges(4,5) * t243 + Ifges(4,6) * t242 + Ifges(4,3) * t268;
t285 = m(5) * t169 - mrSges(5,1) * t183 + mrSges(5,2) * t184 - t219 * t229 + t222 * t230 + t292;
t127 = -mrSges(4,1) * t217 + mrSges(4,3) * t191 + Ifges(4,4) * t216 + Ifges(4,2) * t215 + Ifges(4,6) * t267 - pkin(3) * t285 + pkin(8) * t293 + t275 * t133 + t271 * t139 - t243 * t223 + t268 * t225;
t128 = mrSges(4,2) * t217 - mrSges(4,3) * t190 + Ifges(4,1) * t216 + Ifges(4,4) * t215 + Ifges(4,5) * t267 - pkin(8) * t141 - t133 * t271 + t139 * t275 + t223 * t242 - t224 * t268;
t239 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t273 + Ifges(3,6) * t277) * qJD(1);
t244 = -pkin(6) * t279 + t290;
t283 = m(4) * t217 - mrSges(4,1) * t215 + mrSges(4,2) * t216 - t234 * t242 + t235 * t243 + t285;
t120 = -mrSges(3,1) * t244 + mrSges(3,3) * t233 + Ifges(3,4) * t250 + Ifges(3,2) * t251 + Ifges(3,6) * qJDD(2) - pkin(2) * t283 + pkin(7) * t294 + qJD(2) * t241 + t276 * t127 + t272 * t128 - t239 * t301;
t122 = mrSges(3,2) * t244 - mrSges(3,3) * t232 + Ifges(3,1) * t250 + Ifges(3,4) * t251 + Ifges(3,5) * qJDD(2) - pkin(7) * t132 - qJD(2) * t240 - t127 * t272 + t128 * t276 + t239 * t300;
t281 = -m(3) * t244 + mrSges(3,1) * t251 - mrSges(3,2) * t250 - t253 * t301 + t254 * t300 - t283;
t286 = mrSges(2,1) * t256 - mrSges(2,2) * t257 + Ifges(2,3) * qJDD(1) + pkin(1) * t281 + pkin(6) * t295 + t120 * t277 + t122 * t273;
t144 = m(2) * t256 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t279 + t281;
t126 = t130 * t277 + t131 * t273;
t124 = m(2) * t257 - mrSges(2,1) * t279 - qJDD(1) * mrSges(2,2) + t295;
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t257 + t279 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t126 - t306;
t118 = -mrSges(2,2) * g(3) - mrSges(2,3) * t256 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t279 - pkin(6) * t126 - t120 * t273 + t122 * t277;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t278 * t118 - t274 * t123 - pkin(5) * (t124 * t274 + t144 * t278), t118, t122, t128, t139, -t195 * t265 + t287; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t274 * t118 + t278 * t123 + pkin(5) * (t124 * t278 - t144 * t274), t123, t120, t127, t133, -t230 * t193 + t289; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t286, t286, t306, t307, t308, -t229 * t197 - t288;];
m_new = t1;
