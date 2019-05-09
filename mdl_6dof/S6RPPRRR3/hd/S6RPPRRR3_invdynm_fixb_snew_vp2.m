% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-05-05 15:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:30:53
% EndTime: 2019-05-05 15:31:05
% DurationCPUTime: 6.28s
% Computational Cost: add. (120851->297), mult. (226126->361), div. (0->0), fcn. (136215->10), ass. (0->124)
t271 = sin(qJ(1));
t275 = cos(qJ(1));
t245 = t271 * g(1) - t275 * g(2);
t236 = qJDD(1) * pkin(1) + t245;
t246 = -t275 * g(1) - t271 * g(2);
t277 = qJD(1) ^ 2;
t238 = -t277 * pkin(1) + t246;
t266 = sin(pkin(10));
t267 = cos(pkin(10));
t211 = t266 * t236 + t267 * t238;
t303 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t211;
t302 = -pkin(2) - pkin(7);
t301 = mrSges(3,1) - mrSges(4,2);
t300 = -Ifges(4,4) + Ifges(3,5);
t299 = Ifges(4,5) - Ifges(3,6);
t210 = t267 * t236 - t266 * t238;
t290 = -t277 * qJ(3) + qJDD(3) - t210;
t196 = t302 * qJDD(1) + t290;
t263 = -g(3) + qJDD(2);
t270 = sin(qJ(4));
t274 = cos(qJ(4));
t190 = t270 * t196 + t274 * t263;
t237 = (mrSges(5,1) * t270 + mrSges(5,2) * t274) * qJD(1);
t297 = qJD(1) * qJD(4);
t251 = t274 * t297;
t240 = -t270 * qJDD(1) - t251;
t298 = qJD(1) * t274;
t244 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t298;
t253 = t270 * qJD(1);
t193 = t302 * t277 - t303;
t294 = t270 * t297;
t241 = t274 * qJDD(1) - t294;
t179 = (-t241 + t294) * pkin(8) + (-t240 + t251) * pkin(4) + t193;
t239 = (pkin(4) * t270 - pkin(8) * t274) * qJD(1);
t276 = qJD(4) ^ 2;
t186 = -t276 * pkin(4) + qJDD(4) * pkin(8) - t239 * t253 + t190;
t269 = sin(qJ(5));
t273 = cos(qJ(5));
t168 = t273 * t179 - t269 * t186;
t234 = t273 * qJD(4) - t269 * t298;
t206 = t234 * qJD(5) + t269 * qJDD(4) + t273 * t241;
t233 = qJDD(5) - t240;
t235 = t269 * qJD(4) + t273 * t298;
t248 = t253 + qJD(5);
t166 = (t234 * t248 - t206) * pkin(9) + (t234 * t235 + t233) * pkin(5) + t168;
t169 = t269 * t179 + t273 * t186;
t205 = -t235 * qJD(5) + t273 * qJDD(4) - t269 * t241;
t215 = t248 * pkin(5) - t235 * pkin(9);
t232 = t234 ^ 2;
t167 = -t232 * pkin(5) + t205 * pkin(9) - t248 * t215 + t169;
t268 = sin(qJ(6));
t272 = cos(qJ(6));
t164 = t272 * t166 - t268 * t167;
t208 = t272 * t234 - t268 * t235;
t176 = t208 * qJD(6) + t268 * t205 + t272 * t206;
t209 = t268 * t234 + t272 * t235;
t187 = -t208 * mrSges(7,1) + t209 * mrSges(7,2);
t247 = qJD(6) + t248;
t194 = -t247 * mrSges(7,2) + t208 * mrSges(7,3);
t226 = qJDD(6) + t233;
t159 = m(7) * t164 + t226 * mrSges(7,1) - t176 * mrSges(7,3) - t209 * t187 + t247 * t194;
t165 = t268 * t166 + t272 * t167;
t175 = -t209 * qJD(6) + t272 * t205 - t268 * t206;
t195 = t247 * mrSges(7,1) - t209 * mrSges(7,3);
t160 = m(7) * t165 - t226 * mrSges(7,2) + t175 * mrSges(7,3) + t208 * t187 - t247 * t195;
t152 = t272 * t159 + t268 * t160;
t212 = -t234 * mrSges(6,1) + t235 * mrSges(6,2);
t213 = -t248 * mrSges(6,2) + t234 * mrSges(6,3);
t150 = m(6) * t168 + t233 * mrSges(6,1) - t206 * mrSges(6,3) - t235 * t212 + t248 * t213 + t152;
t214 = t248 * mrSges(6,1) - t235 * mrSges(6,3);
t291 = -t268 * t159 + t272 * t160;
t151 = m(6) * t169 - t233 * mrSges(6,2) + t205 * mrSges(6,3) + t234 * t212 - t248 * t214 + t291;
t292 = -t269 * t150 + t273 * t151;
t143 = m(5) * t190 - qJDD(4) * mrSges(5,2) + t240 * mrSges(5,3) - qJD(4) * t244 - t237 * t253 + t292;
t189 = t274 * t196 - t270 * t263;
t243 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t253;
t185 = -qJDD(4) * pkin(4) - t276 * pkin(8) + t239 * t298 - t189;
t170 = -t205 * pkin(5) - t232 * pkin(9) + t235 * t215 + t185;
t288 = m(7) * t170 - t175 * mrSges(7,1) + t176 * mrSges(7,2) - t208 * t194 + t209 * t195;
t280 = -m(6) * t185 + t205 * mrSges(6,1) - t206 * mrSges(6,2) + t234 * t213 - t235 * t214 - t288;
t155 = m(5) * t189 + qJDD(4) * mrSges(5,1) - t241 * mrSges(5,3) + qJD(4) * t243 - t237 * t298 + t280;
t133 = t270 * t143 + t274 * t155;
t199 = -qJDD(1) * pkin(2) + t290;
t289 = -m(4) * t199 + t277 * mrSges(4,3) - t133;
t128 = m(3) * t210 - t277 * mrSges(3,2) + t301 * qJDD(1) + t289;
t145 = t273 * t150 + t269 * t151;
t142 = -m(5) * t193 + t240 * mrSges(5,1) - t241 * mrSges(5,2) - t243 * t253 - t244 * t298 - t145;
t197 = t277 * pkin(2) + t303;
t283 = -m(4) * t197 + t277 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t142;
t139 = m(3) * t211 - t277 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t283;
t125 = t267 * t128 + t266 * t139;
t293 = -t266 * t128 + t267 * t139;
t134 = t274 * t143 - t270 * t155;
t132 = m(4) * t263 + t134;
t182 = Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t247;
t183 = Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t247;
t287 = -mrSges(7,1) * t164 + mrSges(7,2) * t165 - Ifges(7,5) * t176 - Ifges(7,6) * t175 - Ifges(7,3) * t226 - t209 * t182 + t208 * t183;
t181 = Ifges(7,5) * t209 + Ifges(7,6) * t208 + Ifges(7,3) * t247;
t153 = -mrSges(7,1) * t170 + mrSges(7,3) * t165 + Ifges(7,4) * t176 + Ifges(7,2) * t175 + Ifges(7,6) * t226 - t209 * t181 + t247 * t183;
t154 = mrSges(7,2) * t170 - mrSges(7,3) * t164 + Ifges(7,1) * t176 + Ifges(7,4) * t175 + Ifges(7,5) * t226 + t208 * t181 - t247 * t182;
t200 = Ifges(6,5) * t235 + Ifges(6,6) * t234 + Ifges(6,3) * t248;
t202 = Ifges(6,1) * t235 + Ifges(6,4) * t234 + Ifges(6,5) * t248;
t130 = -mrSges(6,1) * t185 + mrSges(6,3) * t169 + Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t233 - pkin(5) * t288 + pkin(9) * t291 + t272 * t153 + t268 * t154 - t235 * t200 + t248 * t202;
t201 = Ifges(6,4) * t235 + Ifges(6,2) * t234 + Ifges(6,6) * t248;
t136 = mrSges(6,2) * t185 - mrSges(6,3) * t168 + Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t233 - pkin(9) * t152 - t268 * t153 + t272 * t154 + t234 * t200 - t248 * t201;
t223 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t274 - Ifges(5,6) * t270) * qJD(1);
t224 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t274 - Ifges(5,2) * t270) * qJD(1);
t121 = mrSges(5,2) * t193 - mrSges(5,3) * t189 + Ifges(5,1) * t241 + Ifges(5,4) * t240 + Ifges(5,5) * qJDD(4) - pkin(8) * t145 - qJD(4) * t224 - t269 * t130 + t273 * t136 - t223 * t253;
t225 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t274 - Ifges(5,4) * t270) * qJD(1);
t278 = mrSges(6,1) * t168 - mrSges(6,2) * t169 + Ifges(6,5) * t206 + Ifges(6,6) * t205 + Ifges(6,3) * t233 + pkin(5) * t152 + t235 * t201 - t234 * t202 - t287;
t126 = -mrSges(5,1) * t193 + mrSges(5,3) * t190 + Ifges(5,4) * t241 + Ifges(5,2) * t240 + Ifges(5,6) * qJDD(4) - pkin(4) * t145 + qJD(4) * t225 - t223 * t298 - t278;
t286 = mrSges(4,2) * t199 - mrSges(4,3) * t197 + Ifges(4,1) * qJDD(1) - pkin(7) * t133 + t274 * t121 - t270 * t126;
t285 = -mrSges(4,1) * t197 - pkin(3) * t142 - pkin(7) * t134 - t270 * t121 - t274 * t126;
t284 = mrSges(5,1) * t189 - mrSges(5,2) * t190 + Ifges(5,5) * t241 + Ifges(5,6) * t240 + Ifges(5,3) * qJDD(4) + pkin(4) * t280 + pkin(8) * t292 + t273 * t130 + t269 * t136 + t224 * t298 + t225 * t253;
t282 = -mrSges(3,2) * t211 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t289) + qJ(3) * t283 + mrSges(3,1) * t210 + Ifges(3,3) * qJDD(1) + t286;
t281 = mrSges(4,1) * t199 + pkin(3) * t133 + t284;
t279 = mrSges(2,1) * t245 - mrSges(2,2) * t246 + Ifges(2,3) * qJDD(1) + pkin(1) * t125 + t282;
t123 = m(2) * t246 - t277 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t293;
t122 = m(2) * t245 + qJDD(1) * mrSges(2,1) - t277 * mrSges(2,2) + t125;
t119 = (mrSges(3,2) - mrSges(4,3)) * t263 + t299 * t277 + t300 * qJDD(1) + t281 - mrSges(3,3) * t210 - qJ(3) * t132;
t118 = mrSges(3,3) * t211 - pkin(2) * t132 - t299 * qJDD(1) - t301 * t263 + t300 * t277 + t285;
t117 = -mrSges(2,2) * g(3) - mrSges(2,3) * t245 + Ifges(2,5) * qJDD(1) - t277 * Ifges(2,6) - qJ(2) * t125 - t266 * t118 + t267 * t119;
t116 = Ifges(2,6) * qJDD(1) + t277 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t246 + t266 * t119 + t267 * t118 - pkin(1) * (m(3) * t263 + t132) + qJ(2) * t293;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t117 - t271 * t116 - pkin(6) * (t275 * t122 + t271 * t123), t117, t119, t286, t121, t136, t154; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t117 + t275 * t116 + pkin(6) * (-t271 * t122 + t275 * t123), t116, t118, mrSges(4,3) * t263 + Ifges(4,4) * qJDD(1) - t277 * Ifges(4,5) - t281, t126, t130, t153; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t279, t279, t282, -mrSges(4,2) * t263 + t277 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t285, t284, t278, -t287;];
m_new  = t1;
