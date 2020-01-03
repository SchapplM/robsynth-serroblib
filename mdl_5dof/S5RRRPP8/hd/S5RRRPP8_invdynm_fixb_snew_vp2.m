% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:34
% EndTime: 2019-12-31 21:07:40
% DurationCPUTime: 2.67s
% Computational Cost: add. (26565->313), mult. (51618->352), div. (0->0), fcn. (30295->6), ass. (0->113)
t256 = sin(qJ(1));
t258 = cos(qJ(1));
t244 = t256 * g(1) - t258 * g(2);
t260 = qJD(1) ^ 2;
t222 = -qJDD(1) * pkin(1) - t260 * pkin(6) - t244;
t255 = sin(qJ(2));
t257 = cos(qJ(2));
t281 = qJD(1) * qJD(2);
t278 = t257 * t281;
t239 = qJDD(1) * t255 + t278;
t279 = t255 * t281;
t240 = qJDD(1) * t257 - t279;
t159 = (-t239 - t278) * pkin(7) + (-t240 + t279) * pkin(2) + t222;
t245 = -g(1) * t258 - g(2) * t256;
t223 = -pkin(1) * t260 + qJDD(1) * pkin(6) + t245;
t213 = -g(3) * t255 + t257 * t223;
t238 = (-pkin(2) * t257 - pkin(7) * t255) * qJD(1);
t259 = qJD(2) ^ 2;
t282 = qJD(1) * t257;
t166 = -pkin(2) * t259 + qJDD(2) * pkin(7) + t238 * t282 + t213;
t254 = sin(qJ(3));
t293 = cos(qJ(3));
t156 = t293 * t159 - t254 * t166;
t283 = qJD(1) * t255;
t235 = -t293 * qJD(2) + t254 * t283;
t236 = t254 * qJD(2) + t293 * t283;
t202 = pkin(3) * t235 - qJ(4) * t236;
t234 = qJDD(3) - t240;
t247 = -qJD(3) + t282;
t246 = t247 ^ 2;
t155 = -t234 * pkin(3) - t246 * qJ(4) + t236 * t202 + qJDD(4) - t156;
t197 = -t235 * qJD(3) + t254 * qJDD(2) + t293 * t239;
t204 = -mrSges(5,2) * t235 - mrSges(5,3) * t236;
t298 = -m(5) * t155 - t197 * mrSges(5,1) - t236 * t204;
t211 = mrSges(5,1) * t236 - mrSges(5,2) * t247;
t196 = qJD(3) * t236 - t293 * qJDD(2) + t239 * t254;
t207 = pkin(4) * t236 + qJ(5) * t247;
t233 = t235 ^ 2;
t212 = -t257 * g(3) - t255 * t223;
t165 = -qJDD(2) * pkin(2) - t259 * pkin(7) + t238 * t283 - t212;
t289 = t235 * t247;
t295 = -2 * qJD(4);
t266 = (-t197 - t289) * qJ(4) + t165 + (-t247 * pkin(3) + t295) * t236;
t294 = 2 * qJD(5);
t148 = -t233 * pkin(4) + t235 * t294 - t236 * t207 + (pkin(3) + qJ(5)) * t196 + t266;
t208 = mrSges(6,1) * t236 + mrSges(6,3) * t247;
t210 = -mrSges(6,1) * t235 - mrSges(6,2) * t247;
t141 = m(6) * t148 - t197 * mrSges(6,2) + t196 * mrSges(6,3) - t236 * t208 + t235 * t210;
t154 = t196 * pkin(3) + t266;
t209 = mrSges(5,1) * t235 + mrSges(5,3) * t247;
t269 = -m(5) * t154 + t196 * mrSges(5,2) + t235 * t209 - t141;
t138 = -t197 * mrSges(5,3) - t236 * t211 - t269;
t157 = t254 * t159 + t293 * t166;
t170 = -Ifges(6,5) * t247 + Ifges(6,6) * t235 + Ifges(6,3) * t236;
t172 = Ifges(4,5) * t236 - Ifges(4,6) * t235 - Ifges(4,3) * t247;
t177 = -Ifges(5,1) * t247 - Ifges(5,4) * t236 + Ifges(5,5) * t235;
t268 = -t246 * pkin(3) + t234 * qJ(4) - t235 * t202 + t157;
t152 = 0.2e1 * qJD(4) * t247 - t268;
t201 = -mrSges(6,2) * t236 + mrSges(6,3) * t235;
t151 = -t196 * pkin(4) - t233 * qJ(5) + qJDD(5) + (t295 - t207) * t247 + t268;
t176 = -Ifges(6,1) * t247 + Ifges(6,4) * t235 + Ifges(6,5) * t236;
t273 = mrSges(6,1) * t151 - mrSges(6,3) * t148 - Ifges(6,4) * t234 - Ifges(6,2) * t196 - Ifges(6,6) * t197 - t236 * t176;
t280 = m(6) * t151 + t234 * mrSges(6,2) - t247 * t208;
t264 = mrSges(5,1) * t152 - mrSges(5,2) * t154 + pkin(4) * (t196 * mrSges(6,1) + t235 * t201 - t280) + qJ(5) * t141 - t273;
t174 = -Ifges(5,4) * t247 - Ifges(5,2) * t236 + Ifges(5,6) * t235;
t285 = Ifges(4,1) * t236 - Ifges(4,4) * t235 - Ifges(4,5) * t247 - t174;
t290 = Ifges(4,4) + Ifges(5,6);
t122 = -t264 + (-t170 - t285) * t247 + (-t172 - t177) * t236 + (Ifges(4,6) - Ifges(5,5)) * t234 + t290 * t197 + (-Ifges(4,2) - Ifges(5,3)) * t196 - mrSges(4,1) * t165 + mrSges(4,3) * t157 - pkin(3) * t138;
t171 = -Ifges(5,5) * t247 - Ifges(5,6) * t236 + Ifges(5,3) * t235;
t175 = Ifges(4,4) * t236 - Ifges(4,2) * t235 - Ifges(4,6) * t247;
t145 = t247 * t294 + (t235 * t236 - t234) * qJ(5) + (t197 - t289) * pkin(4) + t155;
t276 = -m(6) * t145 + t234 * mrSges(6,3) - t247 * t210;
t142 = t197 * mrSges(6,1) + t236 * t201 - t276;
t173 = -Ifges(6,4) * t247 + Ifges(6,2) * t235 + Ifges(6,6) * t236;
t272 = -mrSges(6,1) * t145 + mrSges(6,2) * t148 - Ifges(6,5) * t234 - Ifges(6,6) * t196 - Ifges(6,3) * t197 + t247 * t173;
t267 = -mrSges(5,1) * t155 + mrSges(5,3) * t154 - pkin(4) * t142 + t272;
t286 = t176 + t177;
t126 = -t267 + (t175 - t171) * t247 + (-t172 - t286) * t235 + (Ifges(4,5) - Ifges(5,4)) * t234 + (Ifges(4,1) + Ifges(5,2)) * t197 - t290 * t196 + mrSges(4,2) * t165 - mrSges(4,3) * t156 - qJ(4) * t138;
t203 = mrSges(4,1) * t235 + mrSges(4,2) * t236;
t205 = mrSges(4,2) * t247 - mrSges(4,3) * t235;
t291 = -mrSges(6,1) - mrSges(4,3);
t134 = m(4) * t156 + (-t205 + t209) * t247 + (-t201 - t203) * t236 + (mrSges(4,1) - mrSges(5,2)) * t234 + t291 * t197 + t276 + t298;
t206 = -mrSges(4,1) * t247 - mrSges(4,3) * t236;
t274 = -m(5) * t152 + t234 * mrSges(5,3) - t247 * t211 + t280;
t284 = -t201 - t204;
t137 = m(4) * t157 - t234 * mrSges(4,2) + t247 * t206 + (-t203 + t284) * t235 + (-mrSges(5,1) + t291) * t196 + t274;
t132 = -t134 * t254 + t293 * t137;
t135 = -m(4) * t165 - t196 * mrSges(4,1) - t235 * t205 + (-t206 + t211) * t236 + (-mrSges(4,2) + mrSges(5,3)) * t197 + t269;
t220 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t255 + Ifges(3,2) * t257) * qJD(1);
t221 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t255 + Ifges(3,4) * t257) * qJD(1);
t297 = mrSges(3,1) * t212 - mrSges(3,2) * t213 + Ifges(3,5) * t239 + Ifges(3,6) * t240 + Ifges(3,3) * qJDD(2) + pkin(2) * t135 + pkin(7) * t132 + (t220 * t255 - t221 * t257) * qJD(1) + t293 * t122 + t254 * t126;
t271 = -mrSges(6,2) * t151 + mrSges(6,3) * t145 - Ifges(6,1) * t234 - Ifges(6,4) * t196 - Ifges(6,5) * t197 - t235 * t170;
t263 = -mrSges(5,2) * t155 + mrSges(5,3) * t152 - Ifges(5,1) * t234 + Ifges(5,4) * t197 - Ifges(5,5) * t196 + qJ(5) * t142 + t236 * t171 + t271;
t296 = t285 * t235 + (t175 - t173) * t236 + mrSges(4,1) * t156 - mrSges(4,2) * t157 + Ifges(4,5) * t197 - Ifges(4,6) * t196 + Ifges(4,3) * t234 + pkin(3) * (-t234 * mrSges(5,2) + t247 * t209 - t142 + t298) + qJ(4) * (t284 * t235 + (-mrSges(5,1) - mrSges(6,1)) * t196 + t274) - t263;
t288 = t236 * t173;
t237 = (-mrSges(3,1) * t257 + mrSges(3,2) * t255) * qJD(1);
t242 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t283;
t130 = m(3) * t213 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t240 - qJD(2) * t242 + t237 * t282 + t132;
t243 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t282;
t133 = m(3) * t212 + qJDD(2) * mrSges(3,1) - t239 * mrSges(3,3) + qJD(2) * t243 - t237 * t283 + t135;
t277 = t257 * t130 - t133 * t255;
t131 = t293 * t134 + t254 * t137;
t219 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t255 + Ifges(3,6) * t257) * qJD(1);
t119 = mrSges(3,2) * t222 - mrSges(3,3) * t212 + Ifges(3,1) * t239 + Ifges(3,4) * t240 + Ifges(3,5) * qJDD(2) - pkin(7) * t131 - qJD(2) * t220 - t254 * t122 + t293 * t126 + t219 * t282;
t121 = -mrSges(3,1) * t222 + mrSges(3,3) * t213 + Ifges(3,4) * t239 + Ifges(3,2) * t240 + Ifges(3,6) * qJDD(2) - pkin(2) * t131 + qJD(2) * t221 - t219 * t283 - t296;
t265 = -m(3) * t222 + t240 * mrSges(3,1) - t239 * mrSges(3,2) - t242 * t283 + t243 * t282 - t131;
t270 = mrSges(2,1) * t244 - mrSges(2,2) * t245 + Ifges(2,3) * qJDD(1) + pkin(1) * t265 + pkin(6) * t277 + t255 * t119 + t257 * t121;
t127 = m(2) * t244 + qJDD(1) * mrSges(2,1) - t260 * mrSges(2,2) + t265;
t125 = t130 * t255 + t133 * t257;
t123 = m(2) * t245 - mrSges(2,1) * t260 - qJDD(1) * mrSges(2,2) + t277;
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t245 + t260 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t125 - t297;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t244 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t260 - pkin(6) * t125 + t119 * t257 - t121 * t255;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t258 * t116 - t256 * t117 - pkin(5) * (t123 * t256 + t127 * t258), t116, t119, t126, -t235 * t174 - t263 - t288, -t271 - t288; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t256 * t116 + t258 * t117 + pkin(5) * (t123 * t258 - t127 * t256), t117, t121, t122, Ifges(5,4) * t234 - Ifges(5,2) * t197 + Ifges(5,6) * t196 + t247 * t171 + t286 * t235 + t267, t247 * t170 - t273; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t270, t270, t297, t296, t264 + t236 * t177 + Ifges(5,5) * t234 + Ifges(5,3) * t196 - Ifges(5,6) * t197 + (t170 - t174) * t247, -t235 * t176 - t272;];
m_new = t1;
