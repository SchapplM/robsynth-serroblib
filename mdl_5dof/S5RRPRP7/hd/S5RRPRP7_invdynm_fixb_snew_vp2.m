% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:00:01
% EndTime: 2019-12-31 20:00:12
% DurationCPUTime: 5.45s
% Computational Cost: add. (61099->307), mult. (137644->379), div. (0->0), fcn. (90442->8), ass. (0->114)
t299 = -2 * qJD(3);
t261 = sin(qJ(2));
t263 = cos(qJ(2));
t286 = qJD(1) * qJD(2);
t244 = t261 * qJDD(1) + t263 * t286;
t262 = sin(qJ(1));
t264 = cos(qJ(1));
t251 = -t264 * g(1) - t262 * g(2);
t266 = qJD(1) ^ 2;
t239 = -t266 * pkin(1) + qJDD(1) * pkin(6) + t251;
t293 = t261 * t239;
t295 = pkin(2) * t266;
t195 = qJDD(2) * pkin(2) - t244 * qJ(3) - t293 + (qJ(3) * t286 + t261 * t295 - g(3)) * t263;
t225 = -t261 * g(3) + t263 * t239;
t245 = t263 * qJDD(1) - t261 * t286;
t289 = qJD(1) * t261;
t247 = qJD(2) * pkin(2) - qJ(3) * t289;
t257 = t263 ^ 2;
t196 = t245 * qJ(3) - qJD(2) * t247 - t257 * t295 + t225;
t258 = sin(pkin(8));
t259 = cos(pkin(8));
t234 = (t258 * t263 + t259 * t261) * qJD(1);
t169 = t259 * t195 - t258 * t196 + t234 * t299;
t233 = (t258 * t261 - t259 * t263) * qJD(1);
t170 = t258 * t195 + t259 * t196 + t233 * t299;
t211 = t233 * pkin(3) - t234 * pkin(7);
t265 = qJD(2) ^ 2;
t165 = -t265 * pkin(3) + qJDD(2) * pkin(7) - t233 * t211 + t170;
t250 = t262 * g(1) - t264 * g(2);
t276 = -qJDD(1) * pkin(1) - t250;
t201 = -t245 * pkin(2) + qJDD(3) + t247 * t289 + (-qJ(3) * t257 - pkin(6)) * t266 + t276;
t220 = -t258 * t244 + t259 * t245;
t221 = t259 * t244 + t258 * t245;
t167 = (qJD(2) * t233 - t221) * pkin(7) + (qJD(2) * t234 - t220) * pkin(3) + t201;
t260 = sin(qJ(4));
t296 = cos(qJ(4));
t161 = -t260 * t165 + t296 * t167;
t162 = t296 * t165 + t260 * t167;
t222 = -t296 * qJD(2) + t260 * t234;
t223 = t260 * qJD(2) + t296 * t234;
t232 = qJD(4) + t233;
t173 = Ifges(6,5) * t223 + Ifges(6,6) * t232 + Ifges(6,3) * t222;
t176 = Ifges(5,4) * t223 - Ifges(5,2) * t222 + Ifges(5,6) * t232;
t178 = Ifges(5,1) * t223 - Ifges(5,4) * t222 + Ifges(5,5) * t232;
t187 = t223 * qJD(4) - t296 * qJDD(2) + t260 * t221;
t188 = -t222 * qJD(4) + t260 * qJDD(2) + t296 * t221;
t198 = t222 * mrSges(6,1) - t223 * mrSges(6,3);
t219 = qJDD(4) - t220;
t197 = t222 * pkin(4) - t223 * qJ(5);
t231 = t232 ^ 2;
t157 = -t231 * pkin(4) + t219 * qJ(5) + 0.2e1 * qJD(5) * t232 - t222 * t197 + t162;
t159 = -t219 * pkin(4) - t231 * qJ(5) + t223 * t197 + qJDD(5) - t161;
t177 = Ifges(6,1) * t223 + Ifges(6,4) * t232 + Ifges(6,5) * t222;
t275 = mrSges(6,1) * t159 - mrSges(6,3) * t157 - Ifges(6,4) * t188 - Ifges(6,2) * t219 - Ifges(6,6) * t187 - t222 * t177;
t202 = -t222 * mrSges(6,2) + t232 * mrSges(6,3);
t280 = -m(6) * t159 + t219 * mrSges(6,1) + t232 * t202;
t205 = -t232 * mrSges(6,1) + t223 * mrSges(6,2);
t285 = m(6) * t157 + t219 * mrSges(6,3) + t232 * t205;
t298 = -(-t176 + t173) * t223 + mrSges(5,1) * t161 - mrSges(5,2) * t162 + Ifges(5,5) * t188 - Ifges(5,6) * t187 + Ifges(5,3) * t219 + pkin(4) * (-t188 * mrSges(6,2) - t223 * t198 + t280) + qJ(5) * (-t187 * mrSges(6,2) - t222 * t198 + t285) + t222 * t178 - t275;
t210 = t233 * mrSges(4,1) + t234 * mrSges(4,2);
t227 = qJD(2) * mrSges(4,1) - t234 * mrSges(4,3);
t204 = t232 * mrSges(5,1) - t223 * mrSges(5,3);
t290 = -t222 * mrSges(5,1) - t223 * mrSges(5,2) - t198;
t294 = -mrSges(5,3) - mrSges(6,2);
t147 = m(5) * t162 - t219 * mrSges(5,2) + t294 * t187 - t232 * t204 + t290 * t222 + t285;
t203 = -t232 * mrSges(5,2) - t222 * mrSges(5,3);
t149 = m(5) * t161 + t219 * mrSges(5,1) + t294 * t188 + t232 * t203 + t290 * t223 + t280;
t282 = t296 * t147 - t260 * t149;
t137 = m(4) * t170 - qJDD(2) * mrSges(4,2) + t220 * mrSges(4,3) - qJD(2) * t227 - t233 * t210 + t282;
t226 = -qJD(2) * mrSges(4,2) - t233 * mrSges(4,3);
t164 = -qJDD(2) * pkin(3) - t265 * pkin(7) + t234 * t211 - t169;
t160 = -0.2e1 * qJD(5) * t223 + (t222 * t232 - t188) * qJ(5) + (t223 * t232 + t187) * pkin(4) + t164;
t154 = m(6) * t160 + t187 * mrSges(6,1) - t188 * mrSges(6,3) + t222 * t202 - t223 * t205;
t270 = -m(5) * t164 - t187 * mrSges(5,1) - t188 * mrSges(5,2) - t222 * t203 - t223 * t204 - t154;
t144 = m(4) * t169 + qJDD(2) * mrSges(4,1) - t221 * mrSges(4,3) + qJD(2) * t226 - t234 * t210 + t270;
t130 = t258 * t137 + t259 * t144;
t224 = -t263 * g(3) - t293;
t236 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t261 + Ifges(3,2) * t263) * qJD(1);
t237 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t261 + Ifges(3,4) * t263) * qJD(1);
t279 = -mrSges(6,1) * t160 + mrSges(6,2) * t157;
t175 = Ifges(6,4) * t223 + Ifges(6,2) * t232 + Ifges(6,6) * t222;
t292 = -Ifges(5,5) * t223 + Ifges(5,6) * t222 - Ifges(5,3) * t232 - t175;
t134 = -mrSges(5,1) * t164 + mrSges(5,3) * t162 - pkin(4) * t154 + (t177 + t178) * t232 + t292 * t223 + (Ifges(5,6) - Ifges(6,6)) * t219 + (Ifges(5,4) - Ifges(6,5)) * t188 + (-Ifges(5,2) - Ifges(6,3)) * t187 + t279;
t274 = mrSges(6,2) * t159 - mrSges(6,3) * t160 + Ifges(6,1) * t188 + Ifges(6,4) * t219 + Ifges(6,5) * t187 + t232 * t173;
t139 = mrSges(5,2) * t164 - mrSges(5,3) * t161 + Ifges(5,1) * t188 - Ifges(5,4) * t187 + Ifges(5,5) * t219 - qJ(5) * t154 - t232 * t176 + t292 * t222 + t274;
t207 = Ifges(4,4) * t234 - Ifges(4,2) * t233 + Ifges(4,6) * qJD(2);
t208 = Ifges(4,1) * t234 - Ifges(4,4) * t233 + Ifges(4,5) * qJD(2);
t271 = -mrSges(4,1) * t169 + mrSges(4,2) * t170 - Ifges(4,5) * t221 - Ifges(4,6) * t220 - Ifges(4,3) * qJDD(2) - pkin(3) * t270 - pkin(7) * t282 - t296 * t134 - t260 * t139 - t234 * t207 - t233 * t208;
t297 = mrSges(3,1) * t224 - mrSges(3,2) * t225 + Ifges(3,5) * t244 + Ifges(3,6) * t245 + Ifges(3,3) * qJDD(2) + pkin(2) * t130 + (t261 * t236 - t263 * t237) * qJD(1) - t271;
t141 = t260 * t147 + t296 * t149;
t288 = qJD(1) * t263;
t243 = (-mrSges(3,1) * t263 + mrSges(3,2) * t261) * qJD(1);
t249 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t288;
t128 = m(3) * t224 + qJDD(2) * mrSges(3,1) - t244 * mrSges(3,3) + qJD(2) * t249 - t243 * t289 + t130;
t248 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t289;
t283 = t259 * t137 - t258 * t144;
t129 = m(3) * t225 - qJDD(2) * mrSges(3,2) + t245 * mrSges(3,3) - qJD(2) * t248 + t243 * t288 + t283;
t284 = -t261 * t128 + t263 * t129;
t206 = Ifges(4,5) * t234 - Ifges(4,6) * t233 + Ifges(4,3) * qJD(2);
t125 = mrSges(4,2) * t201 - mrSges(4,3) * t169 + Ifges(4,1) * t221 + Ifges(4,4) * t220 + Ifges(4,5) * qJDD(2) - pkin(7) * t141 - qJD(2) * t207 - t260 * t134 + t296 * t139 - t233 * t206;
t126 = -mrSges(4,1) * t201 + mrSges(4,3) * t170 + Ifges(4,4) * t221 + Ifges(4,2) * t220 + Ifges(4,6) * qJDD(2) - pkin(3) * t141 + qJD(2) * t208 - t234 * t206 - t298;
t235 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t261 + Ifges(3,6) * t263) * qJD(1);
t238 = -t266 * pkin(6) + t276;
t272 = m(4) * t201 - t220 * mrSges(4,1) + t221 * mrSges(4,2) + t233 * t226 + t234 * t227 + t141;
t119 = -mrSges(3,1) * t238 + mrSges(3,3) * t225 + Ifges(3,4) * t244 + Ifges(3,2) * t245 + Ifges(3,6) * qJDD(2) - pkin(2) * t272 + qJ(3) * t283 + qJD(2) * t237 + t258 * t125 + t259 * t126 - t235 * t289;
t121 = mrSges(3,2) * t238 - mrSges(3,3) * t224 + Ifges(3,1) * t244 + Ifges(3,4) * t245 + Ifges(3,5) * qJDD(2) - qJ(3) * t130 - qJD(2) * t236 + t259 * t125 - t258 * t126 + t235 * t288;
t269 = -m(3) * t238 + t245 * mrSges(3,1) - t244 * mrSges(3,2) - t248 * t289 + t249 * t288 - t272;
t273 = mrSges(2,1) * t250 - mrSges(2,2) * t251 + Ifges(2,3) * qJDD(1) + pkin(1) * t269 + pkin(6) * t284 + t263 * t119 + t261 * t121;
t131 = m(2) * t250 + qJDD(1) * mrSges(2,1) - t266 * mrSges(2,2) + t269;
t124 = t263 * t128 + t261 * t129;
t122 = m(2) * t251 - t266 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t284;
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t251 + t266 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t124 - t297;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t250 + Ifges(2,5) * qJDD(1) - t266 * Ifges(2,6) - pkin(6) * t124 - t261 * t119 + t263 * t121;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t264 * t116 - t262 * t117 - pkin(5) * (t262 * t122 + t264 * t131), t116, t121, t125, t139, -t222 * t175 + t274; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t262 * t116 + t264 * t117 + pkin(5) * (t264 * t122 - t262 * t131), t117, t119, t126, t134, -t223 * t173 - t275; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t273, t273, t297, -t271, t298, Ifges(6,5) * t188 + Ifges(6,6) * t219 + Ifges(6,3) * t187 + t223 * t175 - t232 * t177 - t279;];
m_new = t1;
