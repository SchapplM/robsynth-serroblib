% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:46
% EndTime: 2019-12-31 22:00:01
% DurationCPUTime: 6.61s
% Computational Cost: add. (81170->304), mult. (162509->372), div. (0->0), fcn. (109193->8), ass. (0->114)
t262 = sin(qJ(3));
t266 = cos(qJ(3));
t263 = sin(qJ(2));
t290 = qJD(1) * t263;
t245 = qJD(2) * t262 + t266 * t290;
t267 = cos(qJ(2));
t288 = qJD(1) * qJD(2);
t285 = t267 * t288;
t248 = qJDD(1) * t263 + t285;
t215 = -qJD(3) * t245 + qJDD(2) * t266 - t248 * t262;
t244 = qJD(2) * t266 - t262 * t290;
t216 = qJD(3) * t244 + qJDD(2) * t262 + t248 * t266;
t261 = sin(qJ(4));
t265 = cos(qJ(4));
t218 = t244 * t265 - t245 * t261;
t181 = qJD(4) * t218 + t215 * t261 + t216 * t265;
t219 = t244 * t261 + t245 * t265;
t195 = -mrSges(6,1) * t218 + mrSges(6,2) * t219;
t264 = sin(qJ(1));
t268 = cos(qJ(1));
t253 = t264 * g(1) - t268 * g(2);
t270 = qJD(1) ^ 2;
t237 = -qJDD(1) * pkin(1) - t270 * pkin(6) - t253;
t257 = t263 * t288;
t249 = qJDD(1) * t267 - t257;
t199 = (-t248 - t285) * pkin(7) + (-t249 + t257) * pkin(2) + t237;
t254 = -g(1) * t268 - g(2) * t264;
t238 = -pkin(1) * t270 + qJDD(1) * pkin(6) + t254;
t224 = -g(3) * t263 + t267 * t238;
t247 = (-pkin(2) * t267 - pkin(7) * t263) * qJD(1);
t269 = qJD(2) ^ 2;
t289 = t267 * qJD(1);
t203 = -pkin(2) * t269 + qJDD(2) * pkin(7) + t247 * t289 + t224;
t182 = t199 * t266 - t262 * t203;
t243 = qJDD(3) - t249;
t256 = qJD(3) - t289;
t161 = (t244 * t256 - t216) * pkin(8) + (t244 * t245 + t243) * pkin(3) + t182;
t183 = t199 * t262 + t266 * t203;
t225 = pkin(3) * t256 - pkin(8) * t245;
t242 = t244 ^ 2;
t163 = -pkin(3) * t242 + pkin(8) * t215 - t225 * t256 + t183;
t154 = t161 * t265 - t261 * t163;
t239 = qJDD(4) + t243;
t255 = qJD(4) + t256;
t149 = -0.2e1 * qJD(5) * t219 + (t218 * t255 - t181) * qJ(5) + (t218 * t219 + t239) * pkin(4) + t154;
t204 = -mrSges(6,2) * t255 + mrSges(6,3) * t218;
t287 = m(6) * t149 + t239 * mrSges(6,1) + t255 * t204;
t146 = -t181 * mrSges(6,3) - t219 * t195 + t287;
t155 = t161 * t261 + t163 * t265;
t180 = -qJD(4) * t219 + t215 * t265 - t216 * t261;
t189 = Ifges(5,4) * t219 + Ifges(5,2) * t218 + Ifges(5,6) * t255;
t190 = Ifges(6,1) * t219 + Ifges(6,4) * t218 + Ifges(6,5) * t255;
t191 = Ifges(5,1) * t219 + Ifges(5,4) * t218 + Ifges(5,5) * t255;
t206 = pkin(4) * t255 - qJ(5) * t219;
t217 = t218 ^ 2;
t152 = -pkin(4) * t217 + t180 * qJ(5) + 0.2e1 * qJD(5) * t218 - t206 * t255 + t155;
t188 = Ifges(6,4) * t219 + Ifges(6,2) * t218 + Ifges(6,6) * t255;
t278 = -mrSges(6,1) * t149 + mrSges(6,2) * t152 - Ifges(6,5) * t181 - Ifges(6,6) * t180 - Ifges(6,3) * t239 - t188 * t219;
t294 = mrSges(5,1) * t154 - mrSges(5,2) * t155 + Ifges(5,5) * t181 + Ifges(5,6) * t180 + Ifges(5,3) * t239 + pkin(4) * t146 + t219 * t189 - t278 + (-t191 - t190) * t218;
t196 = -mrSges(5,1) * t218 + mrSges(5,2) * t219;
t205 = -mrSges(5,2) * t255 + mrSges(5,3) * t218;
t139 = m(5) * t154 + t239 * mrSges(5,1) + t255 * t205 + (-t195 - t196) * t219 + (-mrSges(5,3) - mrSges(6,3)) * t181 + t287;
t207 = mrSges(6,1) * t255 - mrSges(6,3) * t219;
t208 = mrSges(5,1) * t255 - mrSges(5,3) * t219;
t286 = m(6) * t152 + t180 * mrSges(6,3) + t195 * t218;
t142 = m(5) * t155 + t180 * mrSges(5,3) + t218 * t196 + (-t207 - t208) * t255 + (-mrSges(5,2) - mrSges(6,2)) * t239 + t286;
t137 = t139 * t265 + t142 * t261;
t210 = Ifges(4,4) * t245 + Ifges(4,2) * t244 + Ifges(4,6) * t256;
t211 = Ifges(4,1) * t245 + Ifges(4,4) * t244 + Ifges(4,5) * t256;
t293 = mrSges(4,1) * t182 - mrSges(4,2) * t183 + Ifges(4,5) * t216 + Ifges(4,6) * t215 + Ifges(4,3) * t243 + pkin(3) * t137 + t245 * t210 - t244 * t211 + t294;
t223 = -t267 * g(3) - t263 * t238;
t202 = -qJDD(2) * pkin(2) - pkin(7) * t269 + t247 * t290 - t223;
t164 = -pkin(3) * t215 - pkin(8) * t242 + t245 * t225 + t202;
t186 = Ifges(6,5) * t219 + Ifges(6,6) * t218 + Ifges(6,3) * t255;
t187 = Ifges(5,5) * t219 + Ifges(5,6) * t218 + Ifges(5,3) * t255;
t158 = -t180 * pkin(4) - qJ(5) * t217 + t206 * t219 + qJDD(5) + t164;
t279 = -mrSges(6,1) * t158 + mrSges(6,3) * t152 + Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t239 + t190 * t255;
t281 = m(6) * t158 - t180 * mrSges(6,1) + t181 * mrSges(6,2) - t204 * t218 + t207 * t219;
t132 = Ifges(5,4) * t181 + Ifges(5,2) * t180 + Ifges(5,6) * t239 + t255 * t191 - mrSges(5,1) * t164 + mrSges(5,3) * t155 - pkin(4) * t281 + qJ(5) * (-t239 * mrSges(6,2) - t255 * t207 + t286) + (-t187 - t186) * t219 + t279;
t277 = mrSges(6,2) * t158 - mrSges(6,3) * t149 + Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t239 + t186 * t218;
t136 = mrSges(5,2) * t164 - mrSges(5,3) * t154 + Ifges(5,1) * t181 + Ifges(5,4) * t180 + Ifges(5,5) * t239 - qJ(5) * t146 + t218 * t187 + (-t188 - t189) * t255 + t277;
t209 = Ifges(4,5) * t245 + Ifges(4,6) * t244 + Ifges(4,3) * t256;
t275 = m(5) * t164 - t180 * mrSges(5,1) + t181 * mrSges(5,2) - t205 * t218 + t208 * t219 + t281;
t282 = -t139 * t261 + t142 * t265;
t121 = -mrSges(4,1) * t202 + mrSges(4,3) * t183 + Ifges(4,4) * t216 + Ifges(4,2) * t215 + Ifges(4,6) * t243 - pkin(3) * t275 + pkin(8) * t282 + t265 * t132 + t261 * t136 - t245 * t209 + t256 * t211;
t122 = mrSges(4,2) * t202 - mrSges(4,3) * t182 + Ifges(4,1) * t216 + Ifges(4,4) * t215 + Ifges(4,5) * t243 - pkin(8) * t137 - t132 * t261 + t136 * t265 + t209 * t244 - t210 * t256;
t220 = -mrSges(4,1) * t244 + mrSges(4,2) * t245;
t221 = -mrSges(4,2) * t256 + mrSges(4,3) * t244;
t134 = m(4) * t182 + mrSges(4,1) * t243 - mrSges(4,3) * t216 - t220 * t245 + t221 * t256 + t137;
t222 = mrSges(4,1) * t256 - mrSges(4,3) * t245;
t135 = m(4) * t183 - mrSges(4,2) * t243 + mrSges(4,3) * t215 + t220 * t244 - t222 * t256 + t282;
t131 = -t134 * t262 + t135 * t266;
t144 = -m(4) * t202 + t215 * mrSges(4,1) - mrSges(4,2) * t216 + t244 * t221 - t222 * t245 - t275;
t235 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t263 + Ifges(3,2) * t267) * qJD(1);
t236 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t263 + Ifges(3,4) * t267) * qJD(1);
t292 = mrSges(3,1) * t223 - mrSges(3,2) * t224 + Ifges(3,5) * t248 + Ifges(3,6) * t249 + Ifges(3,3) * qJDD(2) + pkin(2) * t144 + pkin(7) * t131 + t266 * t121 + t262 * t122 + (t235 * t263 - t236 * t267) * qJD(1);
t246 = (-mrSges(3,1) * t267 + mrSges(3,2) * t263) * qJD(1);
t251 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t290;
t129 = m(3) * t224 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t249 - qJD(2) * t251 + t246 * t289 + t131;
t252 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t289;
t143 = m(3) * t223 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t248 + qJD(2) * t252 - t246 * t290 + t144;
t283 = t129 * t267 - t143 * t263;
t130 = t134 * t266 + t135 * t262;
t234 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t263 + Ifges(3,6) * t267) * qJD(1);
t118 = mrSges(3,2) * t237 - mrSges(3,3) * t223 + Ifges(3,1) * t248 + Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) - pkin(7) * t130 - qJD(2) * t235 - t121 * t262 + t122 * t266 + t234 * t289;
t120 = -mrSges(3,1) * t237 + mrSges(3,3) * t224 + Ifges(3,4) * t248 + Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2) - pkin(2) * t130 + qJD(2) * t236 - t234 * t290 - t293;
t274 = -m(3) * t237 + t249 * mrSges(3,1) - mrSges(3,2) * t248 - t251 * t290 + t252 * t289 - t130;
t276 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + Ifges(2,3) * qJDD(1) + pkin(1) * t274 + pkin(6) * t283 + t118 * t263 + t120 * t267;
t126 = m(2) * t253 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t270 + t274;
t125 = t129 * t263 + t143 * t267;
t123 = m(2) * t254 - mrSges(2,1) * t270 - qJDD(1) * mrSges(2,2) + t283;
t116 = mrSges(2,1) * g(3) + mrSges(2,3) * t254 + t270 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t125 - t292;
t115 = -mrSges(2,2) * g(3) - mrSges(2,3) * t253 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t270 - pkin(6) * t125 + t118 * t267 - t120 * t263;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t268 * t115 - t264 * t116 - pkin(5) * (t123 * t264 + t126 * t268), t115, t118, t122, t136, -t188 * t255 + t277; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t115 + t268 * t116 + pkin(5) * (t123 * t268 - t126 * t264), t116, t120, t121, t132, -t219 * t186 + t279; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t276, t276, t292, t293, t294, -t218 * t190 - t278;];
m_new = t1;
