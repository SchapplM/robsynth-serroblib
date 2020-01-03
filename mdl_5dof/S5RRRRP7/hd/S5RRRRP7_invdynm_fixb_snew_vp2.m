% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP7
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:14
% EndTime: 2019-12-31 21:56:25
% DurationCPUTime: 5.72s
% Computational Cost: add. (72873->308), mult. (147081->378), div. (0->0), fcn. (98614->8), ass. (0->115)
t257 = sin(qJ(3));
t258 = sin(qJ(2));
t260 = cos(qJ(3));
t261 = cos(qJ(2));
t231 = (t257 * t258 - t260 * t261) * qJD(1);
t232 = (t257 * t261 + t258 * t260) * qJD(1);
t282 = qJD(1) * qJD(2);
t239 = t258 * qJDD(1) + t261 * t282;
t240 = t261 * qJDD(1) - t258 * t282;
t206 = -t232 * qJD(3) - t257 * t239 + t260 * t240;
t207 = -t231 * qJD(3) + t260 * t239 + t257 * t240;
t284 = qJD(1) * t258;
t244 = qJD(2) * pkin(2) - pkin(7) * t284;
t255 = t261 ^ 2;
t263 = qJD(1) ^ 2;
t259 = sin(qJ(1));
t262 = cos(qJ(1));
t245 = t259 * g(1) - t262 * g(2);
t273 = -qJDD(1) * pkin(1) - t245;
t208 = -t240 * pkin(2) + t244 * t284 + (-pkin(7) * t255 - pkin(6)) * t263 + t273;
t253 = qJD(2) + qJD(3);
t161 = (t231 * t253 - t207) * pkin(8) + (t232 * t253 - t206) * pkin(3) + t208;
t246 = -t262 * g(1) - t259 * g(2);
t234 = -t263 * pkin(1) + qJDD(1) * pkin(6) + t246;
t288 = t258 * t234;
t290 = pkin(2) * t263;
t195 = qJDD(2) * pkin(2) - t239 * pkin(7) - t288 + (pkin(7) * t282 + t258 * t290 - g(3)) * t261;
t222 = -t258 * g(3) + t261 * t234;
t196 = t240 * pkin(7) - qJD(2) * t244 - t255 * t290 + t222;
t167 = t257 * t195 + t260 * t196;
t218 = t231 * pkin(3) - t232 * pkin(8);
t251 = t253 ^ 2;
t252 = qJDD(2) + qJDD(3);
t164 = -t251 * pkin(3) + t252 * pkin(8) - t231 * t218 + t167;
t256 = sin(qJ(4));
t291 = cos(qJ(4));
t158 = t291 * t161 - t256 * t164;
t159 = t256 * t161 + t291 * t164;
t220 = t291 * t232 + t256 * t253;
t176 = t220 * qJD(4) + t256 * t207 - t291 * t252;
t219 = t256 * t232 - t291 * t253;
t177 = -t219 * qJD(4) + t291 * t207 + t256 * t252;
t227 = qJD(4) + t231;
t178 = Ifges(6,5) * t220 + Ifges(6,6) * t227 + Ifges(6,3) * t219;
t181 = Ifges(5,4) * t220 - Ifges(5,2) * t219 + Ifges(5,6) * t227;
t183 = Ifges(5,1) * t220 - Ifges(5,4) * t219 + Ifges(5,5) * t227;
t192 = t219 * mrSges(6,1) - t220 * mrSges(6,3);
t205 = qJDD(4) - t206;
t191 = t219 * pkin(4) - t220 * qJ(5);
t226 = t227 ^ 2;
t154 = -t226 * pkin(4) + t205 * qJ(5) + 0.2e1 * qJD(5) * t227 - t219 * t191 + t159;
t156 = -t205 * pkin(4) - t226 * qJ(5) + t220 * t191 + qJDD(5) - t158;
t182 = Ifges(6,1) * t220 + Ifges(6,4) * t227 + Ifges(6,5) * t219;
t272 = mrSges(6,1) * t156 - mrSges(6,3) * t154 - Ifges(6,4) * t177 - Ifges(6,2) * t205 - Ifges(6,6) * t176 - t219 * t182;
t209 = -t219 * mrSges(6,2) + t227 * mrSges(6,3);
t277 = -m(6) * t156 + t205 * mrSges(6,1) + t227 * t209;
t212 = -t227 * mrSges(6,1) + t220 * mrSges(6,2);
t281 = m(6) * t154 + t205 * mrSges(6,3) + t227 * t212;
t293 = -(-t181 + t178) * t220 + mrSges(5,1) * t158 - mrSges(5,2) * t159 + Ifges(5,5) * t177 - Ifges(5,6) * t176 + Ifges(5,3) * t205 + pkin(4) * (-t177 * mrSges(6,2) - t220 * t192 + t277) + qJ(5) * (-t176 * mrSges(6,2) - t219 * t192 + t281) + t219 * t183 - t272;
t217 = t231 * mrSges(4,1) + t232 * mrSges(4,2);
t224 = t253 * mrSges(4,1) - t232 * mrSges(4,3);
t211 = t227 * mrSges(5,1) - t220 * mrSges(5,3);
t285 = -t219 * mrSges(5,1) - t220 * mrSges(5,2) - t192;
t289 = -mrSges(5,3) - mrSges(6,2);
t144 = m(5) * t159 - t205 * mrSges(5,2) + t289 * t176 - t227 * t211 + t285 * t219 + t281;
t210 = -t227 * mrSges(5,2) - t219 * mrSges(5,3);
t146 = m(5) * t158 + t205 * mrSges(5,1) + t289 * t177 + t227 * t210 + t285 * t220 + t277;
t278 = t291 * t144 - t256 * t146;
t132 = m(4) * t167 - t252 * mrSges(4,2) + t206 * mrSges(4,3) - t231 * t217 - t253 * t224 + t278;
t166 = t260 * t195 - t257 * t196;
t223 = -t253 * mrSges(4,2) - t231 * mrSges(4,3);
t163 = -t252 * pkin(3) - t251 * pkin(8) + t232 * t218 - t166;
t157 = -0.2e1 * qJD(5) * t220 + (t219 * t227 - t177) * qJ(5) + (t220 * t227 + t176) * pkin(4) + t163;
t151 = m(6) * t157 + t176 * mrSges(6,1) - t177 * mrSges(6,3) + t219 * t209 - t220 * t212;
t267 = -m(5) * t163 - t176 * mrSges(5,1) - t177 * mrSges(5,2) - t219 * t210 - t220 * t211 - t151;
t141 = m(4) * t166 + t252 * mrSges(4,1) - t207 * mrSges(4,3) - t232 * t217 + t253 * t223 + t267;
t127 = t257 * t132 + t260 * t141;
t221 = -t261 * g(3) - t288;
t229 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t258 + Ifges(3,2) * t261) * qJD(1);
t230 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t258 + Ifges(3,4) * t261) * qJD(1);
t276 = -mrSges(6,1) * t157 + mrSges(6,2) * t154;
t180 = Ifges(6,4) * t220 + Ifges(6,2) * t227 + Ifges(6,6) * t219;
t287 = -Ifges(5,5) * t220 + Ifges(5,6) * t219 - Ifges(5,3) * t227 - t180;
t134 = -mrSges(5,1) * t163 + mrSges(5,3) * t159 - pkin(4) * t151 + (t182 + t183) * t227 + t287 * t220 + (Ifges(5,6) - Ifges(6,6)) * t205 + (Ifges(5,4) - Ifges(6,5)) * t177 + (-Ifges(5,2) - Ifges(6,3)) * t176 + t276;
t271 = mrSges(6,2) * t156 - mrSges(6,3) * t157 + Ifges(6,1) * t177 + Ifges(6,4) * t205 + Ifges(6,5) * t176 + t227 * t178;
t136 = mrSges(5,2) * t163 - mrSges(5,3) * t158 + Ifges(5,1) * t177 - Ifges(5,4) * t176 + Ifges(5,5) * t205 - qJ(5) * t151 - t227 * t181 + t287 * t219 + t271;
t214 = Ifges(4,4) * t232 - Ifges(4,2) * t231 + Ifges(4,6) * t253;
t215 = Ifges(4,1) * t232 - Ifges(4,4) * t231 + Ifges(4,5) * t253;
t268 = -mrSges(4,1) * t166 + mrSges(4,2) * t167 - Ifges(4,5) * t207 - Ifges(4,6) * t206 - Ifges(4,3) * t252 - pkin(3) * t267 - pkin(8) * t278 - t291 * t134 - t256 * t136 - t232 * t214 - t231 * t215;
t292 = mrSges(3,1) * t221 - mrSges(3,2) * t222 + Ifges(3,5) * t239 + Ifges(3,6) * t240 + Ifges(3,3) * qJDD(2) + pkin(2) * t127 + (t258 * t229 - t261 * t230) * qJD(1) - t268;
t138 = t256 * t144 + t291 * t146;
t283 = qJD(1) * t261;
t238 = (-mrSges(3,1) * t261 + mrSges(3,2) * t258) * qJD(1);
t243 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t283;
t125 = m(3) * t221 + qJDD(2) * mrSges(3,1) - t239 * mrSges(3,3) + qJD(2) * t243 - t238 * t284 + t127;
t242 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t284;
t279 = t260 * t132 - t257 * t141;
t126 = m(3) * t222 - qJDD(2) * mrSges(3,2) + t240 * mrSges(3,3) - qJD(2) * t242 + t238 * t283 + t279;
t280 = -t258 * t125 + t261 * t126;
t213 = Ifges(4,5) * t232 - Ifges(4,6) * t231 + Ifges(4,3) * t253;
t119 = mrSges(4,2) * t208 - mrSges(4,3) * t166 + Ifges(4,1) * t207 + Ifges(4,4) * t206 + Ifges(4,5) * t252 - pkin(8) * t138 - t256 * t134 + t291 * t136 - t231 * t213 - t253 * t214;
t123 = -mrSges(4,1) * t208 + mrSges(4,3) * t167 + Ifges(4,4) * t207 + Ifges(4,2) * t206 + Ifges(4,6) * t252 - pkin(3) * t138 - t232 * t213 + t253 * t215 - t293;
t228 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t258 + Ifges(3,6) * t261) * qJD(1);
t233 = -t263 * pkin(6) + t273;
t269 = m(4) * t208 - t206 * mrSges(4,1) + t207 * mrSges(4,2) + t231 * t223 + t232 * t224 + t138;
t116 = -mrSges(3,1) * t233 + mrSges(3,3) * t222 + Ifges(3,4) * t239 + Ifges(3,2) * t240 + Ifges(3,6) * qJDD(2) - pkin(2) * t269 + pkin(7) * t279 + qJD(2) * t230 + t257 * t119 + t260 * t123 - t228 * t284;
t118 = mrSges(3,2) * t233 - mrSges(3,3) * t221 + Ifges(3,1) * t239 + Ifges(3,4) * t240 + Ifges(3,5) * qJDD(2) - pkin(7) * t127 - qJD(2) * t229 + t260 * t119 - t257 * t123 + t228 * t283;
t266 = -m(3) * t233 + t240 * mrSges(3,1) - t239 * mrSges(3,2) - t242 * t284 + t243 * t283 - t269;
t270 = mrSges(2,1) * t245 - mrSges(2,2) * t246 + Ifges(2,3) * qJDD(1) + pkin(1) * t266 + pkin(6) * t280 + t261 * t116 + t258 * t118;
t128 = m(2) * t245 + qJDD(1) * mrSges(2,1) - t263 * mrSges(2,2) + t266;
t122 = t261 * t125 + t258 * t126;
t120 = m(2) * t246 - t263 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t280;
t114 = mrSges(2,1) * g(3) + mrSges(2,3) * t246 + t263 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t122 - t292;
t113 = -mrSges(2,2) * g(3) - mrSges(2,3) * t245 + Ifges(2,5) * qJDD(1) - t263 * Ifges(2,6) - pkin(6) * t122 - t258 * t116 + t261 * t118;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t262 * t113 - t259 * t114 - pkin(5) * (t259 * t120 + t262 * t128), t113, t118, t119, t136, -t219 * t180 + t271; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t259 * t113 + t262 * t114 + pkin(5) * (t262 * t120 - t259 * t128), t114, t116, t123, t134, -t220 * t178 - t272; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t270, t270, t292, -t268, t293, Ifges(6,5) * t177 + Ifges(6,6) * t205 + Ifges(6,3) * t176 + t220 * t180 - t227 * t182 - t276;];
m_new = t1;
