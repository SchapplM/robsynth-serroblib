% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:56
% EndTime: 2019-12-31 21:53:08
% DurationCPUTime: 5.88s
% Computational Cost: add. (74929->307), mult. (151286->378), div. (0->0), fcn. (101743->8), ass. (0->116)
t259 = sin(qJ(3));
t260 = sin(qJ(2));
t263 = cos(qJ(3));
t264 = cos(qJ(2));
t235 = (t259 * t260 - t263 * t264) * qJD(1);
t286 = qJD(1) * qJD(2);
t243 = qJDD(1) * t260 + t264 * t286;
t244 = qJDD(1) * t264 - t260 * t286;
t210 = -qJD(3) * t235 + t243 * t263 + t244 * t259;
t236 = (t259 * t264 + t260 * t263) * qJD(1);
t255 = qJD(2) + qJD(3);
t258 = sin(qJ(4));
t262 = cos(qJ(4));
t224 = -t236 * t258 + t255 * t262;
t254 = qJDD(2) + qJDD(3);
t181 = qJD(4) * t224 + t210 * t262 + t254 * t258;
t225 = t236 * t262 + t255 * t258;
t195 = -mrSges(6,1) * t224 + mrSges(6,2) * t225;
t209 = -qJD(3) * t236 - t243 * t259 + t244 * t263;
t288 = qJD(1) * t260;
t248 = qJD(2) * pkin(2) - pkin(7) * t288;
t257 = t264 ^ 2;
t266 = qJD(1) ^ 2;
t261 = sin(qJ(1));
t265 = cos(qJ(1));
t249 = g(1) * t261 - t265 * g(2);
t277 = -qJDD(1) * pkin(1) - t249;
t211 = -pkin(2) * t244 + t248 * t288 + (-pkin(7) * t257 - pkin(6)) * t266 + t277;
t161 = (t235 * t255 - t210) * pkin(8) + (t236 * t255 - t209) * pkin(3) + t211;
t250 = -g(1) * t265 - g(2) * t261;
t238 = -pkin(1) * t266 + qJDD(1) * pkin(6) + t250;
t291 = t238 * t260;
t293 = pkin(2) * t266;
t198 = qJDD(2) * pkin(2) - pkin(7) * t243 - t291 + (pkin(7) * t286 + t260 * t293 - g(3)) * t264;
t227 = -g(3) * t260 + t264 * t238;
t199 = pkin(7) * t244 - qJD(2) * t248 - t257 * t293 + t227;
t167 = t259 * t198 + t263 * t199;
t222 = pkin(3) * t235 - pkin(8) * t236;
t253 = t255 ^ 2;
t164 = -pkin(3) * t253 + pkin(8) * t254 - t222 * t235 + t167;
t155 = t262 * t161 - t164 * t258;
t208 = qJDD(4) - t209;
t231 = qJD(4) + t235;
t151 = -0.2e1 * qJD(5) * t225 + (t224 * t231 - t181) * qJ(5) + (t224 * t225 + t208) * pkin(4) + t155;
t212 = -mrSges(6,2) * t231 + mrSges(6,3) * t224;
t285 = m(6) * t151 + t208 * mrSges(6,1) + t231 * t212;
t148 = -mrSges(6,3) * t181 - t195 * t225 + t285;
t156 = t258 * t161 + t262 * t164;
t180 = -qJD(4) * t225 - t210 * t258 + t254 * t262;
t186 = Ifges(5,4) * t225 + Ifges(5,2) * t224 + Ifges(5,6) * t231;
t187 = Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t231;
t188 = Ifges(5,1) * t225 + Ifges(5,4) * t224 + Ifges(5,5) * t231;
t214 = pkin(4) * t231 - qJ(5) * t225;
t223 = t224 ^ 2;
t154 = -pkin(4) * t223 + qJ(5) * t180 + 0.2e1 * qJD(5) * t224 - t214 * t231 + t156;
t185 = Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t231;
t275 = -mrSges(6,1) * t151 + mrSges(6,2) * t154 - Ifges(6,5) * t181 - Ifges(6,6) * t180 - Ifges(6,3) * t208 - t225 * t185;
t295 = mrSges(5,1) * t155 - mrSges(5,2) * t156 + Ifges(5,5) * t181 + Ifges(5,6) * t180 + Ifges(5,3) * t208 + pkin(4) * t148 + t225 * t186 - (t188 + t187) * t224 - t275;
t221 = mrSges(4,1) * t235 + mrSges(4,2) * t236;
t229 = mrSges(4,1) * t255 - mrSges(4,3) * t236;
t196 = -mrSges(5,1) * t224 + mrSges(5,2) * t225;
t213 = -mrSges(5,2) * t231 + mrSges(5,3) * t224;
t140 = m(5) * t155 + mrSges(5,1) * t208 + t213 * t231 + (-t195 - t196) * t225 + (-mrSges(5,3) - mrSges(6,3)) * t181 + t285;
t284 = m(6) * t154 + t180 * mrSges(6,3) + t224 * t195;
t215 = mrSges(6,1) * t231 - mrSges(6,3) * t225;
t289 = -mrSges(5,1) * t231 + mrSges(5,3) * t225 - t215;
t292 = -mrSges(5,2) - mrSges(6,2);
t143 = m(5) * t156 + mrSges(5,3) * t180 + t196 * t224 + t208 * t292 + t289 * t231 + t284;
t281 = -t140 * t258 + t262 * t143;
t133 = m(4) * t167 - mrSges(4,2) * t254 + mrSges(4,3) * t209 - t221 * t235 - t229 * t255 + t281;
t166 = t198 * t263 - t259 * t199;
t228 = -mrSges(4,2) * t255 - mrSges(4,3) * t235;
t163 = -pkin(3) * t254 - pkin(8) * t253 + t236 * t222 - t166;
t158 = -pkin(4) * t180 - qJ(5) * t223 + t214 * t225 + qJDD(5) + t163;
t280 = -m(6) * t158 + t180 * mrSges(6,1) + t224 * t212;
t270 = -m(5) * t163 + t180 * mrSges(5,1) + t181 * t292 + t224 * t213 + t289 * t225 + t280;
t145 = m(4) * t166 + mrSges(4,1) * t254 - mrSges(4,3) * t210 - t221 * t236 + t228 * t255 + t270;
t126 = t259 * t133 + t263 * t145;
t226 = -g(3) * t264 - t291;
t233 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t260 + Ifges(3,2) * t264) * qJD(1);
t234 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t260 + Ifges(3,4) * t264) * qJD(1);
t183 = Ifges(6,5) * t225 + Ifges(6,6) * t224 + Ifges(6,3) * t231;
t184 = Ifges(5,5) * t225 + Ifges(5,6) * t224 + Ifges(5,3) * t231;
t276 = -mrSges(6,1) * t158 + mrSges(6,3) * t154 + Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t208 + t231 * t187;
t128 = Ifges(5,4) * t181 + Ifges(5,2) * t180 + Ifges(5,6) * t208 + t231 * t188 - mrSges(5,1) * t163 + mrSges(5,3) * t156 - pkin(4) * (mrSges(6,2) * t181 - t280) + qJ(5) * (-mrSges(6,2) * t208 - t215 * t231 + t284) + (-pkin(4) * t215 - t183 - t184) * t225 + t276;
t274 = mrSges(6,2) * t158 - mrSges(6,3) * t151 + Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t208 + t224 * t183;
t135 = mrSges(5,2) * t163 - mrSges(5,3) * t155 + Ifges(5,1) * t181 + Ifges(5,4) * t180 + Ifges(5,5) * t208 - qJ(5) * t148 + t184 * t224 + (-t185 - t186) * t231 + t274;
t218 = Ifges(4,4) * t236 - Ifges(4,2) * t235 + Ifges(4,6) * t255;
t219 = Ifges(4,1) * t236 - Ifges(4,4) * t235 + Ifges(4,5) * t255;
t271 = -mrSges(4,1) * t166 + mrSges(4,2) * t167 - Ifges(4,5) * t210 - Ifges(4,6) * t209 - Ifges(4,3) * t254 - pkin(3) * t270 - pkin(8) * t281 - t262 * t128 - t258 * t135 - t236 * t218 - t235 * t219;
t294 = mrSges(3,1) * t226 - mrSges(3,2) * t227 + Ifges(3,5) * t243 + Ifges(3,6) * t244 + Ifges(3,3) * qJDD(2) + pkin(2) * t126 + (t260 * t233 - t264 * t234) * qJD(1) - t271;
t137 = t262 * t140 + t258 * t143;
t287 = qJD(1) * t264;
t242 = (-mrSges(3,1) * t264 + mrSges(3,2) * t260) * qJD(1);
t247 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t287;
t124 = m(3) * t226 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t243 + qJD(2) * t247 - t242 * t288 + t126;
t246 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t288;
t282 = t263 * t133 - t259 * t145;
t125 = m(3) * t227 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t244 - qJD(2) * t246 + t242 * t287 + t282;
t283 = -t124 * t260 + t264 * t125;
t217 = Ifges(4,5) * t236 - Ifges(4,6) * t235 + Ifges(4,3) * t255;
t118 = mrSges(4,2) * t211 - mrSges(4,3) * t166 + Ifges(4,1) * t210 + Ifges(4,4) * t209 + Ifges(4,5) * t254 - pkin(8) * t137 - t128 * t258 + t135 * t262 - t217 * t235 - t218 * t255;
t122 = -mrSges(4,1) * t211 + mrSges(4,3) * t167 + Ifges(4,4) * t210 + Ifges(4,2) * t209 + Ifges(4,6) * t254 - pkin(3) * t137 - t236 * t217 + t255 * t219 - t295;
t232 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t260 + Ifges(3,6) * t264) * qJD(1);
t237 = -pkin(6) * t266 + t277;
t272 = m(4) * t211 - t209 * mrSges(4,1) + mrSges(4,2) * t210 + t235 * t228 + t229 * t236 + t137;
t115 = -mrSges(3,1) * t237 + mrSges(3,3) * t227 + Ifges(3,4) * t243 + Ifges(3,2) * t244 + Ifges(3,6) * qJDD(2) - pkin(2) * t272 + pkin(7) * t282 + qJD(2) * t234 + t259 * t118 + t263 * t122 - t232 * t288;
t117 = mrSges(3,2) * t237 - mrSges(3,3) * t226 + Ifges(3,1) * t243 + Ifges(3,4) * t244 + Ifges(3,5) * qJDD(2) - pkin(7) * t126 - qJD(2) * t233 + t118 * t263 - t122 * t259 + t232 * t287;
t269 = -m(3) * t237 + t244 * mrSges(3,1) - mrSges(3,2) * t243 - t246 * t288 + t247 * t287 - t272;
t273 = mrSges(2,1) * t249 - mrSges(2,2) * t250 + Ifges(2,3) * qJDD(1) + pkin(1) * t269 + pkin(6) * t283 + t264 * t115 + t260 * t117;
t129 = m(2) * t249 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t266 + t269;
t121 = t124 * t264 + t125 * t260;
t119 = m(2) * t250 - mrSges(2,1) * t266 - qJDD(1) * mrSges(2,2) + t283;
t113 = mrSges(2,1) * g(3) + mrSges(2,3) * t250 + t266 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t121 - t294;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t266 - pkin(6) * t121 - t115 * t260 + t117 * t264;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t265 * t112 - t261 * t113 - pkin(5) * (t119 * t261 + t129 * t265), t112, t117, t118, t135, -t185 * t231 + t274; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t261 * t112 + t265 * t113 + pkin(5) * (t119 * t265 - t261 * t129), t113, t115, t122, t128, -t225 * t183 + t276; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t273, t273, t294, -t271, t295, -t224 * t187 - t275;];
m_new = t1;
