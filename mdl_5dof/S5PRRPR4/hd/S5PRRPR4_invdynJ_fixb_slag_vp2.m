% Calculate vector of inverse dynamics joint torques for
% S5PRRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:29
% EndTime: 2019-12-05 16:21:51
% DurationCPUTime: 8.48s
% Computational Cost: add. (3291->442), mult. (7455->612), div. (0->0), fcn. (5281->14), ass. (0->203)
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t171 = -qJ(4) - pkin(6);
t204 = qJD(3) * t171;
t116 = qJD(4) * t176 + t173 * t204;
t117 = -qJD(4) * t173 + t176 * t204;
t167 = sin(pkin(9));
t169 = cos(pkin(9));
t128 = t167 * t176 + t169 * t173;
t177 = cos(qJ(2));
t186 = t128 * t177;
t275 = qJD(1) * t186 - t116 * t167 + t169 * t117;
t192 = t167 * t173 - t169 * t176;
t220 = qJD(1) * t177;
t274 = t169 * t116 + t167 * t117 + t192 * t220;
t121 = t192 * qJD(3);
t285 = pkin(7) * t121 + t275;
t120 = t128 * qJD(3);
t284 = -pkin(7) * t120 + t274;
t164 = qJ(3) + pkin(9);
t158 = cos(t164);
t250 = pkin(3) * t176;
t139 = pkin(4) * t158 + t250;
t160 = qJ(5) + t164;
t149 = sin(t160);
t150 = cos(t160);
t154 = pkin(2) + t250;
t157 = sin(t164);
t198 = -mrSges(4,1) * t176 + mrSges(4,2) * t173;
t283 = mrSges(3,1) + m(6) * (pkin(2) + t139) + mrSges(6,1) * t150 - mrSges(6,2) * t149 + m(5) * t154 + mrSges(5,1) * t158 - mrSges(5,2) * t157 + m(4) * pkin(2) - t198;
t282 = mrSges(3,2) + m(6) * (-pkin(7) + t171) - mrSges(6,3) + m(5) * t171 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t168 = sin(pkin(8));
t170 = cos(pkin(8));
t281 = g(1) * t170 + g(2) * t168;
t280 = -t120 / 0.2e1;
t279 = -t121 / 0.2e1;
t213 = qJD(2) * qJD(3);
t134 = qJDD(2) * t176 - t173 * t213;
t135 = qJDD(2) * t173 + t176 * t213;
t76 = t134 * t169 - t135 * t167;
t77 = t134 * t167 + t135 * t169;
t32 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t118 = t192 * qJD(2);
t217 = qJD(2) * t176;
t219 = qJD(2) * t173;
t119 = -t167 * t217 - t169 * t219;
t202 = -t175 * t118 + t119 * t172;
t22 = qJD(5) * t202 + t172 * t76 + t175 * t77;
t66 = -t118 * t172 - t119 * t175;
t23 = -qJD(5) * t66 - t172 * t77 + t175 * t76;
t4 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t278 = -t32 - t4;
t142 = t171 * t173;
t143 = t171 * t176;
t80 = t169 * t142 + t143 * t167;
t58 = -pkin(7) * t128 + t80;
t81 = t167 * t142 - t169 * t143;
t59 = -pkin(7) * t192 + t81;
t25 = t172 * t58 + t175 * t59;
t277 = -qJD(5) * t25 - t284 * t172 + t285 * t175;
t24 = -t172 * t59 + t175 * t58;
t276 = qJD(5) * t24 + t285 * t172 + t284 * t175;
t234 = qJDD(3) / 0.2e1;
t263 = m(5) * pkin(3);
t273 = -mrSges(4,1) - t263;
t252 = pkin(3) * t169;
t151 = pkin(4) + t252;
t253 = pkin(3) * t167;
t114 = t151 * t175 - t172 * t253;
t249 = pkin(7) * t118;
t174 = sin(qJ(2));
t221 = qJD(1) * t174;
t144 = qJD(2) * pkin(6) + t221;
t201 = qJ(4) * qJD(2) + t144;
t108 = t201 * t173;
t109 = t201 * t176;
t227 = t169 * t109;
t50 = t108 * t167 - t227;
t35 = t50 + t249;
t248 = pkin(7) * t119;
t90 = t167 * t109;
t51 = -t169 * t108 - t90;
t36 = t51 + t248;
t272 = t114 * qJD(5) - t172 * t35 - t175 * t36;
t115 = t151 * t172 + t175 * t253;
t271 = -t115 * qJD(5) + t172 * t36 - t175 * t35;
t270 = t174 * t281;
t269 = t176 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t134) - t173 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t135);
t140 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t219;
t141 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t217;
t268 = t140 * t176 + t141 * t173;
t267 = t140 * t173 - t141 * t176;
t214 = qJD(1) * qJD(2);
t153 = t177 * t214;
t137 = t174 * qJDD(1) + t153;
t126 = qJDD(2) * pkin(6) + t137;
t216 = qJD(3) * t173;
t74 = t176 * t126 - t144 * t216;
t215 = qJD(3) * t176;
t208 = t144 * t215;
t75 = -t126 * t173 - t208;
t193 = -t173 * t75 + t176 * t74;
t28 = -mrSges(6,1) * t202 + mrSges(6,2) * t66;
t69 = mrSges(5,1) * t118 - mrSges(5,2) * t119;
t266 = t198 * qJD(2) + t28 + t69;
t265 = 0.2e1 * t234;
t178 = qJD(2) ^ 2;
t262 = -t202 / 0.2e1;
t261 = -t66 / 0.2e1;
t260 = t66 / 0.2e1;
t259 = m(5) + m(6);
t96 = qJD(3) * pkin(3) - t108;
t46 = t169 * t96 - t90;
t31 = qJD(3) * pkin(4) + t248 + t46;
t47 = t167 * t96 + t227;
t33 = t47 - t249;
t9 = -t172 * t33 + t175 * t31;
t258 = t9 * mrSges(6,3);
t256 = -t119 / 0.2e1;
t163 = qJD(3) + qJD(5);
t255 = -t163 / 0.2e1;
t254 = Ifges(6,4) * t66;
t251 = pkin(3) * t173;
t245 = g(3) * t174;
t10 = t172 * t31 + t175 * t33;
t244 = t10 * mrSges(6,3);
t212 = qJD(2) * qJD(4);
t42 = -t208 + qJDD(3) * pkin(3) - qJ(4) * t135 + (-t126 - t212) * t173;
t43 = qJ(4) * t134 + t176 * t212 + t74;
t18 = t167 * t42 + t169 * t43;
t228 = t168 * t177;
t243 = (-t149 * t228 - t150 * t170) * mrSges(6,1) + (t149 * t170 - t150 * t228) * mrSges(6,2);
t226 = t170 * t177;
t242 = (-t149 * t226 + t150 * t168) * mrSges(6,1) + (-t149 * t168 - t150 * t226) * mrSges(6,2);
t241 = mrSges(5,3) * t118;
t240 = mrSges(5,3) * t119;
t239 = Ifges(4,4) * t173;
t238 = Ifges(4,4) * t176;
t237 = Ifges(5,4) * t119;
t224 = t173 * t177;
t222 = t176 * t177;
t218 = qJD(2) * t174;
t211 = pkin(3) * t219;
t210 = pkin(3) * t216;
t152 = t174 * t214;
t17 = -t167 * t43 + t169 * t42;
t136 = qJDD(1) * t177 - t152;
t197 = mrSges(4,1) * t173 + mrSges(4,2) * t176;
t196 = -mrSges(6,1) * t149 - mrSges(6,2) * t150;
t195 = t176 * Ifges(4,2) + t239;
t194 = Ifges(4,5) * t176 - Ifges(4,6) * t173;
t106 = t128 * t174;
t107 = t192 * t174;
t52 = -t106 * t175 + t107 * t172;
t53 = -t106 * t172 - t107 * t175;
t70 = -t128 * t172 - t175 * t192;
t71 = t128 * t175 - t172 * t192;
t161 = qJDD(3) + qJDD(5);
t7 = qJDD(3) * pkin(4) - pkin(7) * t77 + t17;
t8 = pkin(7) * t76 + t18;
t2 = qJD(5) * t9 + t172 * t7 + t175 * t8;
t3 = -qJD(5) * t10 - t172 * t8 + t175 * t7;
t191 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t22 + Ifges(6,6) * t23 + Ifges(6,3) * t161;
t145 = -qJD(2) * pkin(2) - t220;
t188 = t145 * t197;
t187 = t173 * (Ifges(4,1) * t176 - t239);
t125 = -qJDD(2) * pkin(2) - t136;
t182 = t145 * t174 + (t173 ^ 2 + t176 ^ 2) * t177 * t144;
t122 = -qJD(2) * t154 + qJD(4) - t220;
t78 = -pkin(3) * t134 + qJDD(4) + t125;
t155 = Ifges(4,4) * t217;
t138 = -pkin(4) * t157 - t251;
t124 = Ifges(4,1) * t219 + Ifges(4,5) * qJD(3) + t155;
t123 = Ifges(4,6) * qJD(3) + qJD(2) * t195;
t111 = Ifges(5,4) * t118;
t99 = pkin(4) * t192 - t154;
t95 = qJD(3) * mrSges(5,1) + t240;
t94 = -qJD(3) * mrSges(5,2) - t241;
t89 = pkin(4) * t120 + t210;
t88 = -pkin(4) * t119 + t211;
t79 = -mrSges(4,1) * t134 + mrSges(4,2) * t135;
t73 = pkin(4) * t118 + t122;
t68 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t77;
t67 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t76;
t62 = -t119 * Ifges(5,1) + Ifges(5,5) * qJD(3) - t111;
t61 = -t118 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t237;
t60 = Ifges(6,4) * t202;
t55 = -t118 * t177 - t120 * t174;
t54 = -qJD(2) * t186 + t121 * t174;
t49 = mrSges(6,1) * t163 - mrSges(6,3) * t66;
t48 = -mrSges(6,2) * t163 + mrSges(6,3) * t202;
t34 = -pkin(4) * t76 + t78;
t30 = -qJD(5) * t71 - t120 * t175 + t121 * t172;
t29 = qJD(5) * t70 - t120 * t172 - t121 * t175;
t27 = t66 * Ifges(6,1) + t163 * Ifges(6,5) + t60;
t26 = Ifges(6,2) * t202 + t163 * Ifges(6,6) + t254;
t16 = -mrSges(6,2) * t161 + mrSges(6,3) * t23;
t15 = mrSges(6,1) * t161 - mrSges(6,3) * t22;
t12 = -qJD(5) * t53 - t172 * t55 + t175 * t54;
t11 = qJD(5) * t52 + t172 * t54 + t175 * t55;
t1 = [m(2) * qJDD(1) - t106 * t68 - t107 * t67 + t11 * t48 + t12 * t49 + t52 * t15 + t53 * t16 + t54 * t95 + t55 * t94 + (-m(2) - m(3) - m(4) - t259) * g(3) + (qJDD(2) * mrSges(3,1) - t178 * mrSges(3,2) - qJD(2) * t267 + t278 - t79) * t177 + (-t178 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + qJD(2) * t266 - qJD(3) * t268 + t269) * t174 + m(4) * (qJD(2) * t182 - t125 * t177 + t174 * t193) + m(3) * (t136 * t177 + t137 * t174) + m(5) * (-t106 * t17 - t107 * t18 + t122 * t218 - t177 * t78 + t46 * t54 + t47 * t55) + m(6) * (t10 * t11 + t12 * t9 - t177 * t34 + t2 * t53 + t218 * t73 + t3 * t52); t176 * (Ifges(4,4) * t135 + Ifges(4,2) * t134 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t163 * (Ifges(6,5) * t29 + Ifges(6,6) * t30) / 0.2e1 - t154 * t32 + (mrSges(6,2) * t34 - mrSges(6,3) * t3 + Ifges(6,1) * t22 + Ifges(6,4) * t23 + Ifges(6,5) * t161) * t71 + (Ifges(6,1) * t29 + Ifges(6,4) * t30) * t260 + (Ifges(4,5) * t173 + Ifges(4,6) * t176) * t234 + t30 * t244 + (t176 * (-Ifges(4,2) * t173 + t238) + t187) * t213 / 0.2e1 + t202 * (Ifges(6,4) * t29 + Ifges(6,2) * t30) / 0.2e1 + (-t120 * t47 + t121 * t46) * mrSges(5,3) + (-Ifges(5,1) * t121 - Ifges(5,4) * t120) * t256 - t118 * (-Ifges(5,4) * t121 - Ifges(5,2) * t120) / 0.2e1 + t122 * (mrSges(5,1) * t120 - mrSges(5,2) * t121) + t135 * t173 * Ifges(4,1) + t89 * t28 + t99 * t4 + (-mrSges(6,1) * t34 + mrSges(6,3) * t2 + Ifges(6,4) * t22 + Ifges(6,2) * t23 + Ifges(6,6) * t161) * t70 - pkin(2) * t79 + t80 * t68 + t81 * t67 + t73 * (-mrSges(6,1) * t30 + mrSges(6,2) * t29) + t29 * t27 / 0.2e1 + t30 * t26 / 0.2e1 + t24 * t15 + t25 * t16 + t274 * t94 + (t122 * t210 - t154 * t78 + t17 * t80 + t18 * t81 + t274 * t47 + t275 * t46) * m(5) + t275 * t95 + t69 * t210 + (t78 * mrSges(5,2) - t17 * mrSges(5,3) + Ifges(5,1) * t77 + Ifges(5,4) * t76 + Ifges(5,5) * t265) * t128 + (-m(5) * t122 - m(6) * t73 - t266) * t221 + t193 * mrSges(4,3) + t267 * t220 + (m(4) * t193 - t140 * t215 - t141 * t216 + t269) * pkin(6) - (-t78 * mrSges(5,1) + t18 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,2) * t76 + Ifges(5,6) * t265) * t192 + (-t137 + t153) * mrSges(3,2) + (-pkin(2) * t125 - qJD(1) * t182) * m(4) + (t282 * g(3) + t281 * t283) * t174 + (-t283 * g(3) + t281 * t282) * t177 + (t136 + t152) * mrSges(3,1) + t173 * (Ifges(4,4) * t134 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t135 * t238 / 0.2e1 + t62 * t279 + t61 * t280 + t124 * t215 / 0.2e1 - t123 * t216 / 0.2e1 + t276 * t48 + t277 * t49 + (t10 * t276 + t2 * t25 + t24 * t3 + t277 * t9 + t34 * t99 + t73 * t89) * m(6) + (Ifges(5,5) * t279 + Ifges(5,6) * t280 + t188 + t194 * qJD(3) / 0.2e1) * qJD(3) + Ifges(3,3) * qJDD(2) - t29 * t258 + t134 * t195 / 0.2e1 + t125 * t198; Ifges(4,6) * t134 + Ifges(4,5) * t135 + t61 * t256 + (t167 * t18 + t169 * t17) * t263 - t46 * t241 + t68 * t252 + t67 * t253 - (-Ifges(4,2) * t219 + t124 + t155) * t217 / 0.2e1 + (Ifges(6,5) * t255 + t258 + Ifges(6,1) * t261 + Ifges(6,4) * t262 - t73 * mrSges(6,2) - t27 / 0.2e1) * t202 - qJD(3) * (-Ifges(5,5) * t118 + Ifges(5,6) * t119) / 0.2e1 - t122 * (-mrSges(5,1) * t119 - mrSges(5,2) * t118) + t119 * (-Ifges(5,1) * t118 + t237) / 0.2e1 + t114 * t15 + t115 * t16 - t51 * t94 - t50 * t95 + Ifges(5,6) * t76 + Ifges(5,5) * t77 - t88 * t28 - t74 * mrSges(4,2) + t75 * mrSges(4,1) - t18 * mrSges(5,2) + t17 * mrSges(5,1) + (-(-t168 * t222 + t170 * t173) * mrSges(4,2) - m(6) * (t138 * t228 - t139 * t170) - t243 - (-t157 * t228 - t158 * t170) * mrSges(5,1) - (t157 * t170 - t158 * t228) * mrSges(5,2) + t273 * (-t168 * t224 - t170 * t176)) * g(2) + (-(-t168 * t173 - t170 * t222) * mrSges(4,2) - m(6) * (t138 * t226 + t139 * t168) - t242 - (-t157 * t226 + t158 * t168) * mrSges(5,1) - (-t157 * t168 - t158 * t226) * mrSges(5,2) + t273 * (t168 * t176 - t170 * t224)) * g(1) + t191 + t268 * t144 - (Ifges(6,6) * t255 + Ifges(6,4) * t261 + Ifges(6,2) * t262 + t73 * mrSges(6,1) - t26 / 0.2e1 - t244) * t66 + (Ifges(5,2) * t119 - t111 + t62) * t118 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) - qJD(2) * t188 + (m(5) * t251 - m(6) * t138 + mrSges(5,1) * t157 + mrSges(5,2) * t158 - t196 + t197) * t245 - t194 * t213 / 0.2e1 + t123 * t219 / 0.2e1 + t271 * t49 + t272 * t48 + (t272 * t10 + t114 * t3 + t115 * t2 + t271 * t9 - t73 * t88) * m(6) - t47 * t240 - t178 * t187 / 0.2e1 - t69 * t211 - m(5) * (t122 * t211 + t46 * t50 + t47 * t51); t259 * t177 * g(3) + t118 * t94 - t119 * t95 - t202 * t48 + t66 * t49 + (-t10 * t202 + t66 * t9 - t270 + t34) * m(6) + (t118 * t47 - t119 * t46 - t270 + t78) * m(5) - t278; -t73 * (mrSges(6,1) * t66 + mrSges(6,2) * t202) + (Ifges(6,1) * t202 - t254) * t261 + t26 * t260 + (Ifges(6,5) * t202 - Ifges(6,6) * t66) * t255 - t9 * t48 + t10 * t49 - g(1) * t242 - g(2) * t243 - t196 * t245 + (t10 * t66 + t202 * t9) * mrSges(6,3) + t191 + (-Ifges(6,2) * t66 + t27 + t60) * t262;];
tau = t1;
