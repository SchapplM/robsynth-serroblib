% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:14:58
% EndTime: 2019-05-05 14:15:10
% DurationCPUTime: 6.14s
% Computational Cost: add. (117354->299), mult. (232357->371), div. (0->0), fcn. (127280->10), ass. (0->123)
t248 = sin(pkin(10));
t250 = cos(pkin(10));
t253 = sin(qJ(4));
t256 = cos(qJ(4));
t215 = (t248 * t253 - t250 * t256) * qJD(1);
t260 = qJD(1) ^ 2;
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t230 = -t257 * g(1) - t254 * g(2);
t273 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t230;
t286 = -pkin(1) - pkin(2);
t204 = t286 * t260 + t273;
t229 = t254 * g(1) - t257 * g(2);
t272 = -t260 * qJ(2) + qJDD(2) - t229;
t209 = t286 * qJDD(1) + t272;
t249 = sin(pkin(9));
t251 = cos(pkin(9));
t184 = t251 * t204 + t249 * t209;
t179 = -t260 * pkin(3) - qJDD(1) * pkin(7) + t184;
t244 = g(3) + qJDD(3);
t175 = -t253 * t179 + t256 * t244;
t282 = qJD(1) * qJD(4);
t281 = t256 * t282;
t223 = -t253 * qJDD(1) - t281;
t165 = (-t223 - t281) * qJ(5) + (t253 * t256 * t260 + qJDD(4)) * pkin(4) + t175;
t176 = t256 * t179 + t253 * t244;
t224 = -t256 * qJDD(1) + t253 * t282;
t284 = qJD(1) * t253;
t226 = qJD(4) * pkin(4) + qJ(5) * t284;
t243 = t256 ^ 2;
t166 = -t243 * t260 * pkin(4) + t224 * qJ(5) - qJD(4) * t226 + t176;
t287 = 2 * qJD(5);
t161 = t248 * t165 + t250 * t166 + t215 * t287;
t216 = (t248 * t256 + t250 * t253) * qJD(1);
t191 = -t215 * mrSges(6,1) - t216 * mrSges(6,2);
t197 = -t248 * t223 + t250 * t224;
t206 = qJD(4) * mrSges(6,1) + t216 * mrSges(6,3);
t192 = -t215 * pkin(5) + t216 * pkin(8);
t259 = qJD(4) ^ 2;
t158 = -t259 * pkin(5) + qJDD(4) * pkin(8) + t215 * t192 + t161;
t183 = -t249 * t204 + t251 * t209;
t278 = qJDD(1) * pkin(3) - t183;
t167 = -t226 * t284 - t224 * pkin(4) + qJDD(5) + (-qJ(5) * t243 - pkin(7)) * t260 + t278;
t198 = t250 * t223 + t248 * t224;
t162 = (-qJD(4) * t215 - t198) * pkin(8) + (-qJD(4) * t216 - t197) * pkin(5) + t167;
t252 = sin(qJ(6));
t255 = cos(qJ(6));
t155 = -t252 * t158 + t255 * t162;
t201 = t255 * qJD(4) + t252 * t216;
t174 = t201 * qJD(6) + t252 * qJDD(4) + t255 * t198;
t202 = t252 * qJD(4) - t255 * t216;
t181 = -t201 * mrSges(7,1) + t202 * mrSges(7,2);
t213 = qJD(6) - t215;
t185 = -t213 * mrSges(7,2) + t201 * mrSges(7,3);
t196 = qJDD(6) - t197;
t151 = m(7) * t155 + t196 * mrSges(7,1) - t174 * mrSges(7,3) - t202 * t181 + t213 * t185;
t156 = t255 * t158 + t252 * t162;
t173 = -t202 * qJD(6) + t255 * qJDD(4) - t252 * t198;
t186 = t213 * mrSges(7,1) - t202 * mrSges(7,3);
t152 = m(7) * t156 - t196 * mrSges(7,2) + t173 * mrSges(7,3) + t201 * t181 - t213 * t186;
t279 = -t252 * t151 + t255 * t152;
t138 = m(6) * t161 - qJDD(4) * mrSges(6,2) + t197 * mrSges(6,3) - qJD(4) * t206 + t215 * t191 + t279;
t277 = -t250 * t165 + t248 * t166;
t160 = t216 * t287 - t277;
t205 = -qJD(4) * mrSges(6,2) + t215 * mrSges(6,3);
t157 = -qJDD(4) * pkin(5) - t259 * pkin(8) + (-(2 * qJD(5)) - t192) * t216 + t277;
t271 = -m(7) * t157 + t173 * mrSges(7,1) - t174 * mrSges(7,2) + t201 * t185 - t202 * t186;
t147 = m(6) * t160 + qJDD(4) * mrSges(6,1) - t198 * mrSges(6,3) + qJD(4) * t205 + t216 * t191 + t271;
t133 = t248 * t138 + t250 * t147;
t218 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t253 - Ifges(5,2) * t256) * qJD(1);
t219 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t253 - Ifges(5,4) * t256) * qJD(1);
t168 = Ifges(7,5) * t202 + Ifges(7,6) * t201 + Ifges(7,3) * t213;
t170 = Ifges(7,1) * t202 + Ifges(7,4) * t201 + Ifges(7,5) * t213;
t144 = -mrSges(7,1) * t157 + mrSges(7,3) * t156 + Ifges(7,4) * t174 + Ifges(7,2) * t173 + Ifges(7,6) * t196 - t202 * t168 + t213 * t170;
t169 = Ifges(7,4) * t202 + Ifges(7,2) * t201 + Ifges(7,6) * t213;
t145 = mrSges(7,2) * t157 - mrSges(7,3) * t155 + Ifges(7,1) * t174 + Ifges(7,4) * t173 + Ifges(7,5) * t196 + t201 * t168 - t213 * t169;
t188 = -Ifges(6,4) * t216 + Ifges(6,2) * t215 + Ifges(6,6) * qJD(4);
t189 = -Ifges(6,1) * t216 + Ifges(6,4) * t215 + Ifges(6,5) * qJD(4);
t266 = -mrSges(6,1) * t160 + mrSges(6,2) * t161 - Ifges(6,5) * t198 - Ifges(6,6) * t197 - Ifges(6,3) * qJDD(4) - pkin(5) * t271 - pkin(8) * t279 - t255 * t144 - t252 * t145 + t216 * t188 + t215 * t189;
t288 = mrSges(5,1) * t175 - mrSges(5,2) * t176 + Ifges(5,5) * t223 + Ifges(5,6) * t224 + Ifges(5,3) * qJDD(4) + pkin(4) * t133 - (t253 * t218 - t256 * t219) * qJD(1) - t266;
t285 = mrSges(2,1) + mrSges(3,1);
t140 = t255 * t151 + t252 * t152;
t283 = qJD(1) * t256;
t222 = (mrSges(5,1) * t256 - mrSges(5,2) * t253) * qJD(1);
t228 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t283;
t131 = m(5) * t175 + qJDD(4) * mrSges(5,1) - t223 * mrSges(5,3) + qJD(4) * t228 + t222 * t284 + t133;
t227 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t284;
t280 = t250 * t138 - t248 * t147;
t132 = m(5) * t176 - qJDD(4) * mrSges(5,2) + t224 * mrSges(5,3) - qJD(4) * t227 - t222 * t283 + t280;
t127 = -t253 * t131 + t256 * t132;
t123 = m(4) * t184 - t260 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t127;
t178 = -t260 * pkin(7) + t278;
t268 = m(6) * t167 - t197 * mrSges(6,1) + t198 * mrSges(6,2) - t215 * t205 - t216 * t206 + t140;
t135 = -m(5) * t178 + t224 * mrSges(5,1) - t223 * mrSges(5,2) + t227 * t284 - t228 * t283 - t268;
t134 = m(4) * t183 - qJDD(1) * mrSges(4,1) - t260 * mrSges(4,2) + t135;
t121 = t251 * t123 - t249 * t134;
t120 = t249 * t123 + t251 * t134;
t126 = t256 * t131 + t253 * t132;
t210 = -t260 * pkin(1) + t273;
t274 = m(3) * t210 + qJDD(1) * mrSges(3,3) + t121;
t125 = -m(4) * t244 - t126;
t214 = -qJDD(1) * pkin(1) + t272;
t270 = -m(3) * t214 + qJDD(1) * mrSges(3,1) + t260 * mrSges(3,3) - t120;
t187 = -Ifges(6,5) * t216 + Ifges(6,6) * t215 + Ifges(6,3) * qJD(4);
t128 = mrSges(6,2) * t167 - mrSges(6,3) * t160 + Ifges(6,1) * t198 + Ifges(6,4) * t197 + Ifges(6,5) * qJDD(4) - pkin(8) * t140 - qJD(4) * t188 - t252 * t144 + t255 * t145 + t215 * t187;
t264 = mrSges(7,1) * t155 - mrSges(7,2) * t156 + Ifges(7,5) * t174 + Ifges(7,6) * t173 + Ifges(7,3) * t196 + t202 * t169 - t201 * t170;
t129 = -mrSges(6,1) * t167 + mrSges(6,3) * t161 + Ifges(6,4) * t198 + Ifges(6,2) * t197 + Ifges(6,6) * qJDD(4) - pkin(5) * t140 + qJD(4) * t189 + t216 * t187 - t264;
t217 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t253 - Ifges(5,6) * t256) * qJD(1);
t114 = -mrSges(5,1) * t178 + mrSges(5,3) * t176 + Ifges(5,4) * t223 + Ifges(5,2) * t224 + Ifges(5,6) * qJDD(4) - pkin(4) * t268 + qJ(5) * t280 + qJD(4) * t219 + t248 * t128 + t250 * t129 + t217 * t284;
t115 = mrSges(5,2) * t178 - mrSges(5,3) * t175 + Ifges(5,1) * t223 + Ifges(5,4) * t224 + Ifges(5,5) * qJDD(4) - qJ(5) * t133 - qJD(4) * t218 + t250 * t128 - t248 * t129 - t217 * t283;
t112 = mrSges(4,2) * t244 - mrSges(4,3) * t183 - Ifges(4,5) * qJDD(1) - t260 * Ifges(4,6) - pkin(7) * t126 - t253 * t114 + t256 * t115;
t113 = -mrSges(4,1) * t244 + mrSges(4,3) * t184 + t260 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t126 - t288;
t269 = mrSges(3,2) * t214 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t260 * Ifges(3,6) - qJ(3) * t120 + t251 * t112 - t249 * t113;
t267 = mrSges(3,2) * t210 - pkin(2) * t125 - qJ(3) * t121 - t249 * t112 - t251 * t113;
t265 = mrSges(4,1) * t183 - mrSges(4,2) * t184 - Ifges(4,3) * qJDD(1) + pkin(3) * t135 + pkin(7) * t127 + t256 * t114 + t253 * t115;
t263 = -mrSges(3,1) * t214 + mrSges(3,3) * t210 + Ifges(3,2) * qJDD(1) - pkin(2) * t120 - t265;
t262 = -mrSges(2,2) * t230 + qJ(2) * (-t260 * mrSges(3,1) + t274) + pkin(1) * t270 + mrSges(2,1) * t229 + Ifges(2,3) * qJDD(1) + t263;
t124 = -m(3) * g(3) + t125;
t117 = m(2) * t229 + qJDD(1) * mrSges(2,1) - t260 * mrSges(2,2) + t270;
t116 = m(2) * t230 - qJDD(1) * mrSges(2,2) - t285 * t260 + t274;
t110 = -mrSges(2,2) * g(3) - mrSges(2,3) * t229 + Ifges(2,5) * qJDD(1) - t260 * Ifges(2,6) - qJ(2) * t124 + t269;
t109 = mrSges(2,3) * t230 - pkin(1) * t124 + (Ifges(3,4) + Ifges(2,5)) * t260 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t285 * g(3) + t267;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t257 * t110 - t254 * t109 - pkin(6) * (t254 * t116 + t257 * t117), t110, t269, t112, t115, t128, t145; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t254 * t110 + t257 * t109 + pkin(6) * (t257 * t116 - t254 * t117), t109, t263, t113, t114, t129, t144; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t262, t262, -mrSges(3,1) * g(3) - t260 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t267, t265, t288, -t266, t264;];
m_new  = t1;
