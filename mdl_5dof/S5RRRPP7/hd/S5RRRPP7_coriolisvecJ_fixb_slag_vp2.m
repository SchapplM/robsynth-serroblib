% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:57
% EndTime: 2019-12-31 21:04:13
% DurationCPUTime: 6.56s
% Computational Cost: add. (2392->452), mult. (6160->562), div. (0->0), fcn. (3360->4), ass. (0->222)
t298 = qJD(2) / 0.2e1;
t296 = Ifges(6,4) + Ifges(5,5);
t274 = Ifges(5,4) + Ifges(4,5);
t297 = Ifges(6,5) - t274;
t195 = Ifges(3,5) * t298;
t292 = Ifges(6,2) + Ifges(5,3);
t295 = Ifges(4,6) - Ifges(5,6);
t294 = Ifges(5,6) - Ifges(6,6);
t142 = cos(qJ(2));
t213 = qJD(1) * t142;
t293 = qJD(3) - t213;
t287 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t136 = pkin(6) * t213;
t139 = sin(qJ(3));
t291 = qJD(4) * t139 + t136;
t141 = cos(qJ(3));
t290 = t296 * t141;
t289 = t296 * t139;
t134 = Ifges(3,4) * t213;
t140 = sin(qJ(2));
t214 = qJD(1) * t140;
t121 = -qJD(2) * pkin(2) + pkin(6) * t214;
t122 = qJD(2) * pkin(7) + t136;
t116 = -pkin(2) * t142 - pkin(7) * t140 - pkin(1);
t99 = t116 * qJD(1);
t54 = -t139 * t122 + t141 * t99;
t55 = t141 * t122 + t139 * t99;
t155 = t55 * t139 + t54 * t141;
t270 = qJD(4) - t54;
t33 = -pkin(3) * t293 + t270;
t124 = t293 * qJ(4);
t34 = t124 + t55;
t156 = t34 * t139 - t33 * t141;
t234 = Ifges(4,4) * t141;
t168 = -Ifges(4,2) * t139 + t234;
t175 = -mrSges(6,1) * t139 + mrSges(6,2) * t141;
t177 = mrSges(5,1) * t139 - mrSges(5,3) * t141;
t179 = mrSges(4,1) * t139 + mrSges(4,2) * t141;
t197 = t139 * t214;
t204 = t141 * qJD(2);
t109 = t197 - t204;
t30 = qJ(5) * t109 + t55;
t20 = t124 + t30;
t196 = t141 * t214;
t110 = qJD(2) * t139 + t196;
t151 = qJ(4) * t110 - t121;
t259 = pkin(3) + pkin(4);
t21 = -t109 * t259 + qJD(5) + t151;
t247 = t141 / 0.2e1;
t249 = t139 / 0.2e1;
t250 = -t139 / 0.2e1;
t251 = -t293 / 0.2e1;
t252 = t293 / 0.2e1;
t253 = t110 / 0.2e1;
t255 = t109 / 0.2e1;
t256 = -t109 / 0.2e1;
t235 = Ifges(4,4) * t139;
t266 = t141 * t287 - t235 + t289;
t106 = Ifges(4,4) * t109;
t277 = t296 * t109;
t267 = t287 * t110 - t293 * t297 - t106 + t277;
t268 = t139 * t292 + t290;
t280 = t296 * t110;
t272 = t292 * t109 + t293 * t294 + t280;
t36 = pkin(3) * t109 - t151;
t236 = Ifges(4,4) * t110;
t42 = -Ifges(4,2) * t109 + Ifges(4,6) * t293 + t236;
t29 = qJ(5) * t110 + t54;
t271 = qJD(4) - t29;
t9 = -t259 * t293 + t271;
t264 = t156 * mrSges(5,2) + t155 * mrSges(4,3) - (t20 * t139 - t9 * t141) * mrSges(6,3) - t42 * t250 - t36 * t177 - t21 * t175 - (Ifges(6,5) * t141 + Ifges(6,6) * t139) * t251 - t121 * t179 - t168 * t256 - t268 * t255 - (-t139 * t295 + t141 * t274) * t252 - t272 * t249 - t266 * t253 - t267 * t247;
t288 = -t264 + Ifges(3,1) * t214 / 0.2e1 + t134 / 0.2e1 + t195;
t286 = -qJD(2) / 0.2e1;
t203 = qJD(1) * qJD(2);
t193 = t140 * t203;
t202 = qJD(2) * qJD(3);
t210 = qJD(3) * t139;
t72 = t141 * t202 + (-t140 * t210 + t142 * t204) * qJD(1);
t209 = qJD(3) * t141;
t211 = qJD(2) * t142;
t73 = t139 * t202 + (t139 * t211 + t140 * t209) * qJD(1);
t285 = t294 * t193 + t292 * t73 + t296 * t72;
t222 = qJ(4) * t141;
t150 = -t139 * t259 + t222;
t284 = t150 * t293 + t291;
t132 = qJ(4) * t214;
t205 = qJD(5) * t141;
t217 = t140 * t141;
t218 = t139 * t142;
t243 = pkin(7) - qJ(5);
t184 = pkin(2) * t140 - pkin(7) * t142;
t113 = t184 * qJD(1);
t97 = t139 * t113;
t283 = -t210 * t243 - t205 - t132 - t97 - (-pkin(6) * t217 + qJ(5) * t218) * qJD(1);
t158 = pkin(3) * t139 - t222;
t282 = t158 * t293 - t291;
t120 = t243 * t141;
t246 = pkin(6) * t139;
t198 = -pkin(3) - t246;
t192 = -pkin(4) + t198;
t216 = t141 * t142;
t220 = t113 * t141;
t281 = qJD(3) * t120 - qJD(5) * t139 + t220 - (-qJ(5) * t216 + t140 * t192) * qJD(1);
t279 = -t292 * t141 + t289;
t278 = t274 * t139 + t295 * t141;
t276 = (-Ifges(4,4) + t296) * t73 + t287 * t72 - t297 * t193;
t275 = t139 * t287 + t234 - t290;
t254 = -t110 / 0.2e1;
t194 = Ifges(3,6) * t286;
t212 = qJD(2) * t140;
t269 = qJ(4) * t212 - qJD(4) * t142;
t223 = qJ(4) * t139;
t265 = -t141 * t259 - t223;
t199 = Ifges(4,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t188 = Ifges(6,6) / 0.2e1 + t199;
t201 = -Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t189 = -Ifges(6,5) / 0.2e1 - t201;
t200 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1;
t237 = Ifges(3,4) * t140;
t263 = t188 * t109 - t189 * t110 - (Ifges(6,3) / 0.2e1 + t200) * t293 - t20 * mrSges(6,2) - t34 * mrSges(5,3) - t54 * mrSges(4,1) - Ifges(6,5) * t254 - Ifges(5,6) * t255 - t194 + (t142 * Ifges(3,2) + t237) * qJD(1) / 0.2e1 + t33 * mrSges(5,1) + t55 * mrSges(4,2) + t9 * mrSges(6,1) - (Ifges(6,6) + Ifges(4,6)) * t256 - t274 * t253 - (Ifges(6,3) + Ifges(4,3) + Ifges(5,2)) * t252;
t262 = t72 / 0.2e1;
t261 = -t73 / 0.2e1;
t260 = t73 / 0.2e1;
t258 = pkin(1) * mrSges(3,1);
t257 = pkin(1) * mrSges(3,2);
t248 = -t141 / 0.2e1;
t245 = t72 * mrSges(6,3);
t76 = mrSges(6,2) * t293 + mrSges(6,3) * t109;
t81 = -mrSges(5,2) * t109 + mrSges(5,3) * t293;
t242 = t76 + t81;
t239 = mrSges(4,3) * t109;
t77 = -mrSges(4,2) * t293 - t239;
t241 = -t77 - t81;
t238 = mrSges(4,3) * t110;
t79 = mrSges(4,1) * t293 - t238;
t80 = -mrSges(5,1) * t293 + mrSges(5,2) * t110;
t240 = -t79 + t80;
t114 = t184 * qJD(2);
t227 = t139 * t114 + t116 * t209;
t224 = qJ(4) * t109;
t221 = qJD(2) * mrSges(3,2);
t219 = t139 * t140;
t131 = pkin(6) * t216;
t215 = qJD(3) * t131 + t116 * t210;
t85 = t139 * t116 + t131;
t207 = qJD(4) * t141;
t23 = -t73 * mrSges(6,1) + t72 * mrSges(6,2);
t130 = pkin(6) * t218;
t84 = t116 * t141 - t130;
t191 = pkin(6) * t193;
t190 = m(4) * t121 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t109 + mrSges(4,2) * t110 + mrSges(3,3) * t214;
t100 = qJD(1) * t114;
t186 = -t100 * t141 + t122 * t209 + t99 * t210;
t185 = t198 * t140;
t11 = t139 * t100 - t122 * t210 - t141 * t191 + t99 * t209;
t6 = qJ(4) * t193 + qJD(4) * t293 + t11;
t152 = qJD(2) * t185;
t8 = qJD(1) * t152 + t186;
t183 = t139 * t8 + t141 * t6;
t70 = -qJ(4) * t142 + t85;
t182 = -t114 * t141 + t215;
t180 = mrSges(4,1) * t141 - mrSges(4,2) * t139;
t178 = mrSges(5,1) * t141 + mrSges(5,3) * t139;
t176 = mrSges(6,1) * t141 + mrSges(6,2) * t139;
t167 = Ifges(4,2) * t141 + t235;
t160 = Ifges(6,5) * t139 - Ifges(6,6) * t141;
t159 = pkin(3) * t141 + t223;
t12 = t139 * t191 - t186;
t157 = t11 * t141 - t12 * t139;
t75 = -pkin(6) * t196 + t97;
t154 = t192 * qJD(2);
t153 = pkin(6) + t158;
t149 = -pkin(6) + t150;
t148 = qJ(4) * t72 - qJD(2) * t136 + qJD(4) * t110;
t27 = (-t140 * t204 - t142 * t210) * pkin(6) + t227;
t1 = qJ(5) * t73 + qJD(5) * t109 + t6;
t2 = -qJ(5) * t72 - qJD(5) * t110 + t154 * t214 + t186;
t147 = -t12 * mrSges(4,1) + t8 * mrSges(5,1) + t2 * mrSges(6,1) + t11 * mrSges(4,2) - t1 * mrSges(6,2) - t6 * mrSges(5,3);
t138 = t142 * pkin(3);
t128 = Ifges(5,2) * t193;
t127 = Ifges(4,3) * t193;
t119 = t243 * t139;
t118 = mrSges(3,3) * t213 - t221;
t115 = -pkin(2) - t159;
t103 = pkin(2) - t265;
t88 = t153 * t140;
t78 = -mrSges(6,1) * t293 - mrSges(6,3) * t110;
t74 = pkin(6) * t197 + t220;
t71 = t138 - t84;
t69 = t149 * t140;
t68 = Ifges(5,4) * t72;
t67 = Ifges(4,5) * t72;
t66 = Ifges(4,6) * t73;
t65 = Ifges(5,6) * t73;
t64 = t72 * mrSges(5,2);
t61 = -mrSges(6,1) * t109 + mrSges(6,2) * t110;
t60 = mrSges(5,1) * t109 - mrSges(5,3) * t110;
t59 = pkin(3) * t110 + t224;
t58 = qJD(1) * t185 - t220;
t57 = t132 + t75;
t53 = qJ(5) * t219 + t70;
t52 = -mrSges(4,2) * t193 - mrSges(4,3) * t73;
t51 = mrSges(6,2) * t193 + mrSges(6,3) * t73;
t50 = -mrSges(5,1) * t193 + t64;
t49 = mrSges(4,1) * t193 - mrSges(4,3) * t72;
t48 = -mrSges(6,1) * t193 - t245;
t47 = -mrSges(5,2) * t73 + mrSges(5,3) * t193;
t46 = pkin(4) * t142 + t130 + t138 + (-qJ(5) * t140 - t116) * t141;
t31 = -t110 * t259 - t224;
t28 = t212 * t246 - t182;
t26 = (qJD(3) * t159 - t207) * t140 + t153 * t211;
t25 = t152 + t182;
t24 = mrSges(4,1) * t73 + mrSges(4,2) * t72;
t22 = mrSges(5,1) * t73 - mrSges(5,3) * t72;
t19 = t27 + t269;
t15 = t72 * Ifges(4,4) - t73 * Ifges(4,2) + Ifges(4,6) * t193;
t10 = (qJD(3) * t265 + t207) * t140 + t149 * t211;
t7 = pkin(3) * t73 - t148;
t5 = (-qJD(2) * pkin(6) + qJ(5) * qJD(3)) * t217 + (qJD(5) * t140 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t142) * t139 + t227 + t269;
t4 = (-qJ(5) * t211 - t114) * t141 + (qJ(5) * t210 + t154 - t205) * t140 + t215;
t3 = -t259 * t73 + t148;
t13 = [t84 * t49 + t85 * t52 + t88 * t22 + t5 * t76 + t27 * t77 + t4 * t78 + t28 * t79 + t25 * t80 + t19 * t81 + t69 * t23 + t70 * t47 + t71 * t50 + t53 * t51 + t26 * t60 + t10 * t61 + t46 * t48 + m(5) * (t19 * t34 + t25 * t33 + t26 * t36 + t6 * t70 + t7 * t88 + t71 * t8) + m(4) * (t11 * t85 + t12 * t84 + t27 * t55 + t28 * t54) + m(6) * (t1 * t53 + t10 * t21 + t2 * t46 + t20 * t5 + t3 * t69 + t4 * t9) + (-t127 / 0.2e1 - t128 / 0.2e1 - t68 / 0.2e1 + t66 / 0.2e1 - t67 / 0.2e1 - t65 / 0.2e1 + (Ifges(6,6) + t199) * t73 + (Ifges(6,5) + t201) * t72 + (t190 * pkin(6) + (0.3e1 / 0.2e1 * Ifges(3,4) * t142 - 0.2e1 * t257) * qJD(1) + t195 + t288) * qJD(2) + t147) * t142 + (pkin(6) * t24 + t7 * t177 + t3 * t175 + t168 * t261 + (t1 * t139 - t141 * t2) * mrSges(6,3) + (-t11 * t139 - t12 * t141) * mrSges(4,3) + (-t139 * t6 + t141 * t8) * mrSges(5,2) + (-pkin(6) * t118 + t194 - t263) * qJD(2) + (-0.3e1 / 0.2e1 * t237 - 0.2e1 * t258 + t189 * t217 - t188 * t219 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(6,3) + (m(4) * pkin(6) + t179) * pkin(6) - t200) * t142) * t203 + t268 * t260 + t285 * t249 + t266 * t262 + t15 * t250 + t276 * t247 + (-t21 * t176 + t160 * t252 + t121 * t180 + t167 * t255 + t36 * t178 + t42 * t248 + (t139 * t9 + t141 * t20) * mrSges(6,3) + (t139 * t54 - t141 * t55) * mrSges(4,3) + (-t139 * t33 - t141 * t34) * mrSges(5,2) + t279 * t256 + t278 * t251 + t275 * t254 + t267 * t250 + t272 * t247) * qJD(3)) * t140; (-t1 * t141 - t139 * t2) * mrSges(6,3) - t264 * qJD(3) + t275 * t262 + t276 * t249 + t183 * mrSges(5,2) + t3 * t176 - t7 * t178 + t157 * mrSges(4,3) + t281 * t78 + t282 * t60 + (t115 * t7 + t282 * t36 - t33 * t58 - t34 * t57) * m(5) + t283 * t76 + (t1 * t120 + t103 * t3 + t119 * t2 + t283 * t20 + t284 * t21 + t281 * t9) * m(6) + t284 * t61 + t285 * t248 - m(4) * (t54 * t74 + t55 * t75) + t115 * t22 + t119 * t48 + t120 * t51 + t103 * t23 - t75 * t77 - t74 * t79 - t58 * t80 - t57 * t81 - pkin(2) * t24 + t279 * t260 + (((t118 + t221) * pkin(6) + (t258 + t237 / 0.2e1) * qJD(1) + t160 * t286 + t194 + t278 * t298 + t263) * t140 + ((t257 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t140) * qJD(1) - t134 / 0.2e1 + ((-m(4) * pkin(2) - mrSges(3,1) - t180) * qJD(2) - t190) * pkin(6) + t195 - t288) * t142) * qJD(1) + ((-m(4) * t155 - m(5) * t156 + t139 * t241 + t141 * t240) * qJD(3) + m(4) * t157 + m(5) * t183 + (t47 + t52) * t141 + (-t49 + t50) * t139) * pkin(7) + t167 * t261 + t15 * t247; (-t239 + t241) * t54 + t242 * qJD(4) + (t238 - t240) * t55 + (-t109 * t9 - t110 * t20) * mrSges(6,3) + (t109 * t33 + t110 * t34) * mrSges(5,2) - t147 + t127 + t128 + t68 - t66 + t67 + t65 - t121 * (mrSges(4,1) * t110 - mrSges(4,2) * t109) - t36 * (mrSges(5,1) * t110 + mrSges(5,3) * t109) - t21 * (-mrSges(6,1) * t110 - mrSges(6,2) * t109) - Ifges(6,6) * t73 - t29 * t76 - t30 * t78 - Ifges(6,5) * t72 - t59 * t60 - t31 * t61 - pkin(3) * t50 + (t110 * t292 - t277) * t256 + (-pkin(3) * t8 + qJ(4) * t6 + t270 * t34 - t33 * t55 - t36 * t59) * m(5) + (-Ifges(4,2) * t110 - t106 + t267) * t255 + (-t109 * t274 - t110 * t295) * t251 + (-t109 * t287 - t236 + t272 + t280) * t254 + (-Ifges(6,5) * t109 + Ifges(6,6) * t110) * t252 + t42 * t253 - t259 * t48 + (t1 * qJ(4) - t2 * t259 + t20 * t271 - t21 * t31 - t30 * t9) * m(6) + (t47 + t51) * qJ(4) + Ifges(6,3) * t193; -t245 + t64 - t242 * t293 + (t60 - t61) * t110 + (-mrSges(5,1) - mrSges(6,1)) * t193 + (-t110 * t21 - t20 * t293 + t2) * m(6) + (t110 * t36 - t293 * t34 + t8) * m(5); -t109 * t76 + t110 * t78 + 0.2e1 * (t3 / 0.2e1 + t20 * t256 + t9 * t253) * m(6) + t23;];
tauc = t13(:);
