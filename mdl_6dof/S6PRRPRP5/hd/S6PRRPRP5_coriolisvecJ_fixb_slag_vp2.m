% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:45
% EndTime: 2018-11-23 15:13:51
% DurationCPUTime: 5.71s
% Computational Cost: add. (3218->493), mult. (7913->648), div. (0->0), fcn. (4745->8), ass. (0->234)
t296 = Ifges(6,1) + Ifges(7,1);
t292 = Ifges(7,4) + Ifges(6,5);
t152 = sin(qJ(5));
t155 = cos(qJ(5));
t156 = cos(qJ(3));
t226 = qJD(2) * t156;
t108 = qJD(3) * t152 + t155 * t226;
t205 = t152 * t226;
t109 = qJD(3) * t155 - t205;
t149 = qJD(3) * qJ(4);
t154 = sin(qJ(2));
t150 = sin(pkin(6));
t231 = qJD(1) * t150;
t209 = t154 * t231;
t117 = qJD(2) * pkin(8) + t209;
t151 = cos(pkin(6));
t153 = sin(qJ(3));
t230 = qJD(1) * t153;
t133 = t151 * t230;
t84 = t156 * t117 + t133;
t64 = pkin(4) * t226 + t84;
t51 = t149 + t64;
t16 = pkin(5) * t108 - qJ(6) * t109 + t51;
t158 = -pkin(3) - pkin(9);
t240 = t151 * t156;
t134 = qJD(1) * t240;
t200 = pkin(4) * qJD(2) + t117;
t174 = t200 * t153;
t63 = t134 - t174;
t44 = qJD(3) * t158 + qJD(4) - t63;
t201 = -qJ(4) * t153 - pkin(2);
t103 = t156 * t158 + t201;
t157 = cos(qJ(2));
t208 = t157 * t231;
t65 = qJD(2) * t103 - t208;
t14 = -t152 * t65 + t155 * t44;
t15 = t152 * t44 + t155 * t65;
t177 = t14 * t152 - t15 * t155;
t247 = Ifges(7,5) * t152;
t181 = -Ifges(7,3) * t155 + t247;
t250 = Ifges(6,4) * t152;
t185 = Ifges(6,2) * t155 + t250;
t192 = mrSges(7,1) * t155 + mrSges(7,3) * t152;
t194 = mrSges(6,1) * t155 - mrSges(6,2) * t152;
t228 = qJD(2) * t153;
t140 = qJD(5) + t228;
t277 = qJD(6) - t14;
t7 = -pkin(5) * t140 + t277;
t8 = qJ(6) * t140 + t15;
t195 = t152 * t7 + t155 * t8;
t259 = t155 / 0.2e1;
t260 = -t155 / 0.2e1;
t262 = -t152 / 0.2e1;
t264 = -t140 / 0.2e1;
t266 = -t109 / 0.2e1;
t267 = t108 / 0.2e1;
t268 = -t108 / 0.2e1;
t105 = Ifges(6,4) * t108;
t248 = Ifges(7,5) * t108;
t278 = t109 * t296 + t292 * t140 - t105 + t248;
t279 = Ifges(6,6) - Ifges(7,6);
t246 = Ifges(7,5) * t155;
t249 = Ifges(6,4) * t155;
t288 = t152 * t296 - t246 + t249;
t104 = Ifges(7,5) * t109;
t35 = t140 * Ifges(7,6) + t108 * Ifges(7,3) + t104;
t251 = Ifges(6,4) * t109;
t38 = -t108 * Ifges(6,2) + t140 * Ifges(6,6) + t251;
t297 = t259 * t35 + t260 * t38 + t262 * t278 + t16 * t192 + t181 * t268 + t185 * t267 + t51 * t194 + (t152 * t292 + t155 * t279) * t264 + t288 * t266 - t195 * mrSges(7,2) + t177 * mrSges(6,3);
t245 = qJD(2) * pkin(2);
t118 = -t208 - t245;
t282 = qJD(3) / 0.2e1;
t283 = -qJD(3) / 0.2e1;
t284 = qJD(2) / 0.2e1;
t67 = -t149 - t84;
t120 = -pkin(3) * t156 + t201;
t85 = qJD(2) * t120 - t208;
t295 = t84 * mrSges(4,3) + t85 * mrSges(5,2) + Ifges(4,6) * t282 + (Ifges(4,4) * t153 + t156 * Ifges(4,2)) * t284 + Ifges(5,5) * t283 - (-Ifges(5,6) * t153 - t156 * Ifges(5,3)) * qJD(2) / 0.2e1 - t118 * mrSges(4,1) - t67 * mrSges(5,1) + t297;
t269 = pkin(4) + pkin(8);
t294 = -t226 / 0.2e1;
t281 = mrSges(5,1) + mrSges(4,3);
t293 = mrSges(5,2) - mrSges(4,1);
t220 = qJD(2) * qJD(3);
t203 = t156 * t220;
t219 = qJD(3) * qJD(5);
t221 = qJD(5) * t155;
t225 = qJD(3) * t153;
t74 = -t152 * t219 + (t152 * t225 - t156 * t221) * qJD(2);
t202 = t153 * t220;
t75 = -qJD(5) * t205 + (-t202 + t219) * t155;
t291 = (-Ifges(6,4) + Ifges(7,5)) * t75 + t296 * t74 + t292 * t203;
t83 = t117 * t153 - t134;
t290 = -qJD(4) - t83;
t289 = t155 * t296 + t247 - t250;
t287 = -Ifges(4,1) / 0.2e1;
t285 = Ifges(4,4) * t294;
t229 = qJD(2) * t150;
t204 = qJD(1) * t229;
t199 = t157 * t204;
t31 = t153 * t199 + (t156 * t200 + t133) * qJD(3);
t178 = pkin(9) * t153 - qJ(4) * t156;
t223 = qJD(4) * t153;
t163 = qJD(3) * t178 - t223;
t232 = pkin(3) * t202 + t154 * t204;
t43 = qJD(2) * t163 + t232;
t4 = -qJD(5) * t15 - t152 * t43 + t155 * t31;
t224 = qJD(3) * t156;
t115 = t269 * t224;
t127 = t269 * t153;
t243 = t155 * t103 + t152 * t127;
t145 = pkin(3) * t225;
t86 = t145 + t163;
t11 = -qJD(5) * t243 + t115 * t155 - t152 * t86;
t211 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t212 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t213 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t66 = -qJD(3) * pkin(3) - t290;
t275 = t211 * t108 - t213 * t109 - t212 * t140 + t15 * mrSges(6,2) + t7 * mrSges(7,1) + t85 * mrSges(5,3) + Ifges(6,6) * t267 + Ifges(7,6) * t268 + t228 * t287 + Ifges(4,5) * t283 + t285 + Ifges(5,4) * t282 + (-t153 * Ifges(5,2) - Ifges(5,6) * t156) * t284 - t118 * mrSges(4,2) - t14 * mrSges(6,1) - t66 * mrSges(5,1) - t8 * mrSges(7,3) - t83 * mrSges(4,3) + t292 * t266 + (Ifges(6,3) + Ifges(7,2)) * t264;
t273 = -m(5) / 0.2e1;
t271 = -t75 / 0.2e1;
t270 = t75 / 0.2e1;
t265 = t109 / 0.2e1;
t261 = t152 / 0.2e1;
t206 = t157 * t229;
t42 = t117 * t224 + (qJD(3) * t151 + t206) * t230;
t242 = t150 * t154;
t210 = t153 * t242;
t93 = t210 - t240;
t258 = t42 * t93;
t47 = -mrSges(7,2) * t75 + mrSges(7,3) * t203;
t50 = -mrSges(6,2) * t203 - mrSges(6,3) * t75;
t257 = t47 + t50;
t48 = mrSges(6,1) * t203 - mrSges(6,3) * t74;
t49 = -mrSges(7,1) * t203 + t74 * mrSges(7,2);
t256 = t48 - t49;
t253 = mrSges(6,3) * t108;
t78 = -mrSges(6,2) * t140 - t253;
t79 = -mrSges(7,2) * t108 + mrSges(7,3) * t140;
t255 = t78 + t79;
t252 = mrSges(6,3) * t109;
t80 = mrSges(6,1) * t140 - t252;
t81 = -mrSges(7,1) * t140 + mrSges(7,2) * t109;
t254 = -t81 + t80;
t143 = pkin(3) * t228;
t90 = qJD(2) * t178 + t143;
t24 = t152 * t64 + t155 * t90;
t125 = -mrSges(5,1) * t226 - qJD(3) * mrSges(5,3);
t58 = mrSges(6,1) * t108 + mrSges(6,2) * t109;
t244 = -t125 + t58;
t241 = t150 * t157;
t239 = t152 * t157;
t238 = t152 * t158;
t237 = t155 * t158;
t111 = (mrSges(5,2) * t156 - mrSges(5,3) * t153) * qJD(2);
t236 = t111 + (-mrSges(4,1) * t156 + mrSges(4,2) * t153) * qJD(2);
t235 = qJD(3) * t134 + t156 * t199;
t234 = -qJD(3) * t293 - t228 * t281;
t124 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t226;
t233 = t125 - t124;
t128 = t269 * t156;
t227 = qJD(2) * t154;
t222 = qJD(5) * t152;
t218 = pkin(8) * t42 / 0.2e1;
t57 = mrSges(7,1) * t108 - mrSges(7,3) * t109;
t217 = -t57 - t244;
t216 = -0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(5,6);
t215 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t214 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t207 = t150 * t227;
t198 = -t208 / 0.2e1;
t3 = t152 * t31 + t155 * t43 + t44 * t221 - t222 * t65;
t1 = qJ(6) * t203 + qJD(6) * t140 + t3;
t2 = -pkin(5) * t203 - t4;
t197 = -t1 * t152 + t155 * t2;
t196 = -t152 * t3 - t155 * t4;
t193 = mrSges(6,1) * t152 + mrSges(6,2) * t155;
t191 = mrSges(7,1) * t152 - mrSges(7,3) * t155;
t186 = -Ifges(6,2) * t152 + t249;
t182 = Ifges(7,3) * t152 + t246;
t180 = pkin(5) * t155 + qJ(6) * t152;
t179 = -pkin(5) * t152 + qJ(6) * t155;
t23 = -t152 * t90 + t155 * t64;
t53 = -t103 * t152 + t127 * t155;
t171 = -pkin(4) - t180;
t41 = -t117 * t225 + t235;
t61 = t150 * t239 + t155 * t93;
t94 = t151 * t153 + t156 * t242;
t168 = -m(4) * t84 + m(5) * t67 + t233;
t167 = -qJ(4) * t224 - t223;
t10 = -t103 * t222 + t152 * t115 + t127 * t221 + t155 * t86;
t166 = t152 * t257 + t155 * t256;
t165 = -t152 * t254 + t155 * t255;
t164 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t148 = qJD(3) * qJD(4);
t28 = -qJD(3) * t174 + t148 + t235;
t159 = qJD(2) ^ 2;
t138 = Ifges(7,2) * t203;
t137 = Ifges(6,3) * t203;
t119 = qJ(4) - t179;
t114 = t269 * t225;
t113 = -qJ(4) * t226 + t143;
t101 = (mrSges(4,1) * t153 + mrSges(4,2) * t156) * t220;
t100 = (-mrSges(5,2) * t153 - mrSges(5,3) * t156) * t220;
t92 = t145 + t167;
t89 = qJD(5) * t180 - qJD(6) * t155 + qJD(4);
t82 = t156 * t180 + t128;
t77 = (t153 * t239 + t154 * t155) * t231;
t76 = -t153 * t155 * t208 + t152 * t209;
t72 = Ifges(7,4) * t74;
t71 = Ifges(6,5) * t74;
t70 = Ifges(6,6) * t75;
t69 = Ifges(7,6) * t75;
t62 = t152 * t93 - t155 * t241;
t60 = qJD(3) * t94 + t153 * t206;
t59 = qJD(3) * t210 - t151 * t224 - t156 * t206;
t56 = pkin(5) * t109 + qJ(6) * t108;
t55 = qJD(2) * t167 + t232;
t46 = -pkin(5) * t153 - t53;
t45 = qJ(6) * t153 + t243;
t32 = -t148 - t41;
t30 = t134 + (qJD(2) * t171 - t117) * t153;
t27 = mrSges(6,1) * t75 + mrSges(6,2) * t74;
t26 = mrSges(7,1) * t75 - mrSges(7,3) * t74;
t25 = (qJD(5) * t179 + qJD(6) * t152) * t156 + (-pkin(8) + t171) * t225;
t22 = -pkin(5) * t226 - t23;
t21 = qJ(6) * t226 + t24;
t18 = t74 * Ifges(6,4) - t75 * Ifges(6,2) + Ifges(6,6) * t203;
t17 = t74 * Ifges(7,5) + Ifges(7,6) * t203 + t75 * Ifges(7,3);
t13 = -t152 * t207 - t93 * t222 + (qJD(5) * t241 + t60) * t155;
t12 = qJD(5) * t61 + t152 * t60 + t155 * t207;
t9 = -pkin(5) * t224 - t11;
t6 = qJ(6) * t224 + qJD(6) * t153 + t10;
t5 = pkin(5) * t75 - qJ(6) * t74 - qJD(6) * t109 + t28;
t19 = [(t26 + t27) * t94 + t257 * t62 + t256 * t61 - t234 * t60 + t254 * t13 + t255 * t12 + (-t124 + t217) * t59 + (-t159 * t154 * mrSges(3,1) + (-mrSges(3,2) * t159 - t100 - t101) * t157) * t150 + (t236 * t242 + t281 * qJD(3) * (-t153 * t94 + t156 * t93)) * qJD(2) + m(4) * (t41 * t94 + t258 - t59 * t84 + t60 * t83 + (t118 - t208) * t207) + m(5) * (-t32 * t94 + t258 + t59 * t67 + t60 * t66 + (-t157 * t55 + t227 * t85) * t150) + m(6) * (t12 * t15 + t13 * t14 + t28 * t94 + t3 * t62 + t4 * t61 - t51 * t59) + m(7) * (t1 * t62 + t12 * t8 - t13 * t7 - t16 * t59 - t2 * t61 + t5 * t94); -t114 * t58 + t120 * t100 + t128 * t27 + t92 * t111 - pkin(2) * t101 + t9 * t81 + t82 * t26 + t10 * t78 + t6 * t79 + t11 * t80 + t53 * t48 + t25 * t57 + t45 * t47 + t46 * t49 + t243 * t50 + m(6) * (t10 * t15 + t11 * t14 - t114 * t51 + t128 * t28 + t243 * t3 + t4 * t53) - m(7) * (t7 * t76 + t77 * t8) - m(6) * (-t14 * t76 + t15 * t77) + m(7) * (t1 * t45 + t16 * t25 + t2 * t46 + t5 * t82 + t6 * t8 + t7 * t9) + (t137 / 0.2e1 + t138 / 0.2e1 + t72 / 0.2e1 + t69 / 0.2e1 - t70 / 0.2e1 + t71 / 0.2e1 - t55 * mrSges(5,3) - t211 * t75 + t213 * t74 + t281 * t42 + (mrSges(4,2) * t227 + t157 * t234) * t231 + 0.2e1 * (t198 * t66 + t218) * m(5) + 0.2e1 * (t198 * t83 + t218) * m(4) + (pkin(8) * t168 + t214 * qJD(3) + t216 * t228 - t295) * qJD(3) + t164) * t153 + (t185 * t270 + t28 * t194 + t5 * t192 + t55 * mrSges(5,2) + t41 * mrSges(4,3) - t32 * mrSges(5,1) + t17 * t259 + t18 * t260 + t181 * t271 + (t4 * t152 - t3 * t155) * mrSges(6,3) + (-t1 * t155 - t2 * t152) * mrSges(7,2) + (-mrSges(4,1) * t227 + (-m(6) * t51 - m(7) * t16 + t168 - t57 - t58) * t157) * t231 + (t38 * t261 - t51 * t193 - t16 * t191 + t186 * t267 + t182 * t268 + (t14 * t155 + t15 * t152) * mrSges(6,3) + (t8 * t152 - t7 * t155) * mrSges(7,2) + t289 * t266 + (t152 * t279 - t155 * t292) * t140 / 0.2e1 + t278 * t260) * qJD(5) + (((-t152 * t213 - t155 * t211 - t216) * t156 + (0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + t212) * t153) * qJD(2) + t215 * qJD(3) - t275) * qJD(3) - t288 * t74 / 0.2e1 + (qJD(5) * t35 + t291) * t262 + (0.2e1 * t32 * t273 + m(4) * t41 + (m(4) * t83 + m(5) * t66 - t234) * qJD(3)) * pkin(8)) * t156 + 0.2e1 * (t85 * t273 + (-t118 / 0.2e1 - t245 / 0.2e1) * m(4)) * t209 + t254 * t76 - t255 * t77 - t236 * t209 + m(5) * (t120 * t55 + t85 * t92); t119 * t26 - t113 * t111 - t22 * t81 - t24 * t78 - t21 * t79 - t23 * t80 - t63 * t58 - t41 * mrSges(4,2) - t32 * mrSges(5,3) + ((-m(6) * t177 + m(7) * t195 + t165) * t158 + t297) * qJD(5) - m(6) * (t14 * t23 + t15 * t24 + t51 * t63) - m(7) * (t16 * t30 + t21 * t8 + t22 * t7) + qJ(4) * t27 + t28 * t193 + t196 * mrSges(6,3) + t197 * mrSges(7,2) + t5 * t191 + t291 * t259 + t293 * t42 + t289 * t74 / 0.2e1 + (-pkin(3) * t42 - qJ(4) * t32 - t113 * t85 + t290 * t67 - t66 * t84) * m(5) + t166 * t158 + (t89 - t30) * t57 + ((t285 + Ifges(5,6) * t294 + (-pkin(3) * mrSges(5,1) - t152 * t211 + t155 * t213 + t215) * qJD(3) + t275) * t156 + ((-qJ(4) * mrSges(5,1) + t214) * qJD(3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t287 - Ifges(5,2) / 0.2e1) * t226 + (Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t228 + t295) * t153) * qJD(2) + t182 * t270 + t186 * t271 + t17 * t261 + t18 * t262 + t244 * qJD(4) + m(6) * (qJ(4) * t28 + qJD(4) * t51 + t237 * t4 + t238 * t3) + m(7) * (t1 * t238 + t119 * t5 + t16 * t89 - t2 * t237) - t233 * t83 + t234 * t84; t217 * qJD(3) + t165 * qJD(5) + (mrSges(5,1) * t224 + (t111 + t165) * t153) * qJD(2) + t166 + (-qJD(3) * t16 + t140 * t195 - t197) * m(7) + (-qJD(3) * t51 - t140 * t177 - t196) * m(6) + (qJD(3) * t67 + t228 * t85 + t42) * m(5); t137 + t138 + t72 + t69 - t70 + t71 - t16 * (mrSges(7,1) * t109 + mrSges(7,3) * t108) - t51 * (mrSges(6,1) * t109 - mrSges(6,2) * t108) + qJD(6) * t79 - t56 * t57 + qJ(6) * t47 - pkin(5) * t49 + (t108 * t7 + t109 * t8) * mrSges(7,2) + t164 + t38 * t265 + (Ifges(7,3) * t109 - t248) * t268 + (t252 + t254) * t15 + (-t253 - t255) * t14 + (-t108 * t292 - t109 * t279) * t264 + (-pkin(5) * t2 + qJ(6) * t1 - t15 * t7 - t16 * t56 + t277 * t8) * m(7) + (-Ifges(6,2) * t109 - t105 + t278) * t267 + (-t108 * t296 + t104 - t251 + t35) * t266; t109 * t57 - t140 * t79 + 0.2e1 * (t2 / 0.2e1 + t16 * t265 + t8 * t264) * m(7) + t49;];
tauc  = t19(:);
