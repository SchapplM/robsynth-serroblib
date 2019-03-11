% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:09
% EndTime: 2019-03-08 21:04:22
% DurationCPUTime: 7.12s
% Computational Cost: add. (4172->489), mult. (10987->645), div. (0->0), fcn. (7814->10), ass. (0->249)
t158 = sin(qJ(3));
t253 = -qJ(4) - pkin(8);
t195 = qJD(3) * t253;
t161 = cos(qJ(3));
t214 = qJD(4) * t161;
t115 = t158 * t195 + t214;
t116 = -qJD(4) * t158 + t161 * t195;
t154 = sin(pkin(11));
t230 = cos(pkin(11));
t128 = t154 * t161 + t158 * t230;
t162 = cos(qJ(2));
t155 = sin(pkin(6));
t221 = qJD(1) * t155;
t202 = t162 * t221;
t291 = t115 * t154 - t230 * t116 - t128 * t202;
t160 = cos(qJ(6));
t157 = sin(qJ(6));
t185 = mrSges(7,1) * t160 - mrSges(7,2) * t157;
t194 = t230 * t161;
t189 = qJD(2) * t194;
t219 = qJD(2) * t158;
t118 = t154 * t219 - t189;
t263 = pkin(5) * t118;
t159 = sin(qJ(2));
t203 = t159 * t221;
t133 = qJD(2) * pkin(8) + t203;
t192 = qJ(4) * qJD(2) + t133;
t156 = cos(pkin(6));
t220 = qJD(1) * t158;
t201 = t156 * t220;
t90 = t161 * t192 + t201;
t198 = t230 * t90;
t224 = t156 * t161;
t144 = qJD(1) * t224;
t89 = -t158 * t192 + t144;
t86 = qJD(3) * pkin(3) + t89;
t37 = t154 * t86 + t198;
t34 = -qJD(3) * qJ(5) - t37;
t23 = -t34 - t263;
t267 = -t157 / 0.2e1;
t120 = t128 * qJD(2);
t114 = qJD(6) + t120;
t97 = qJD(3) * t160 + t118 * t157;
t265 = Ifges(7,4) * t97;
t96 = -qJD(3) * t157 + t118 * t160;
t30 = Ifges(7,2) * t96 + Ifges(7,6) * t114 + t265;
t95 = Ifges(7,4) * t96;
t31 = Ifges(7,1) * t97 + Ifges(7,5) * t114 + t95;
t300 = t23 * t185 - t160 * t30 / 0.2e1 + t31 * t267;
t294 = mrSges(6,2) - mrSges(5,1);
t216 = qJD(3) * t158;
t121 = qJD(3) * t194 - t154 * t216;
t299 = pkin(5) * t121 + t291;
t119 = t128 * qJD(3);
t153 = pkin(3) * t216;
t169 = -qJ(5) * t121 - qJD(5) * t128 + t153;
t274 = pkin(4) + pkin(9);
t298 = -t119 * t274 - t169 + t203;
t150 = -pkin(3) * t161 - pkin(2);
t111 = qJD(2) * t150 + qJD(4) - t202;
t165 = -qJ(5) * t120 + t111;
t49 = pkin(4) * t118 + t165;
t297 = t111 * mrSges(5,1) - t49 * mrSges(6,2) - t300;
t173 = -qJ(5) * t128 + t150;
t286 = -t154 * t158 + t194;
t56 = -t274 * t286 + t173;
t137 = t253 * t158;
t138 = t253 * t161;
t93 = -t230 * t137 - t138 * t154;
t66 = pkin(5) * t128 + t93;
t19 = -t157 * t56 + t160 * t66;
t296 = qJD(6) * t19 + t157 * t299 - t160 * t298;
t20 = t157 * t66 + t160 * t56;
t295 = -qJD(6) * t20 + t157 * t298 + t160 * t299;
t293 = Ifges(6,4) - Ifges(5,5);
t292 = Ifges(6,5) - Ifges(5,6);
t63 = t115 * t230 + t154 * t116;
t92 = t286 * t202;
t290 = t63 - t92;
t78 = mrSges(5,1) * t118 + mrSges(5,2) * t120;
t79 = -mrSges(6,2) * t118 - mrSges(6,3) * t120;
t289 = t78 + t79;
t251 = mrSges(6,1) * t118;
t102 = -qJD(3) * mrSges(6,3) + t251;
t50 = -mrSges(7,1) * t96 + mrSges(7,2) * t97;
t288 = t50 - t102;
t211 = qJD(2) * qJD(3);
t197 = t158 * t211;
t249 = mrSges(5,3) * t118;
t223 = -qJD(3) * mrSges(5,2) - t102 - t249;
t248 = mrSges(5,3) * t120;
t250 = mrSges(6,1) * t120;
t287 = -qJD(3) * t294 - t248 - t250;
t110 = qJD(3) * t189 - t154 * t197;
t109 = qJD(2) * t119;
t57 = qJD(6) * t96 + t109 * t157;
t32 = mrSges(7,1) * t110 - mrSges(7,3) * t57;
t58 = -qJD(6) * t97 + t109 * t160;
t33 = -mrSges(7,2) * t110 + mrSges(7,3) * t58;
t284 = t157 * t33 + t160 * t32;
t215 = qJD(3) * t161;
t99 = t133 * t161 + t201;
t164 = -t99 * qJD(3) + (-qJ(4) * t215 + (-qJD(4) - t202) * t158) * qJD(2);
t225 = t155 * t162;
t199 = qJD(2) * t225;
t191 = t161 * t199;
t76 = qJD(1) * t191 + qJD(3) * t144 - t133 * t216;
t55 = (-qJ(4) * t216 + t214) * qJD(2) + t76;
t16 = t154 * t55 - t230 * t164;
t12 = pkin(5) * t110 + t16;
t218 = qJD(2) * t159;
t200 = t155 * t218;
t117 = pkin(3) * t197 + qJD(1) * t200;
t167 = -qJ(5) * t110 - qJD(5) * t120 + t117;
t21 = t109 * t274 + t167;
t82 = t154 * t90;
t36 = t230 * t86 - t82;
t172 = qJD(5) - t36;
t262 = t120 * pkin(5);
t22 = -qJD(3) * t274 + t172 + t262;
t35 = t118 * t274 + t165;
t5 = -t157 * t35 + t160 * t22;
t1 = qJD(6) * t5 + t12 * t157 + t160 * t21;
t6 = t157 * t22 + t160 * t35;
t2 = -qJD(6) * t6 + t12 * t160 - t157 * t21;
t283 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t57 + Ifges(7,6) * t58;
t281 = t5 * mrSges(7,1) + t111 * mrSges(5,2) - t6 * mrSges(7,2) - t49 * mrSges(6,3);
t163 = qJD(2) ^ 2;
t280 = t57 / 0.2e1;
t279 = t58 / 0.2e1;
t278 = -t96 / 0.2e1;
t277 = t96 / 0.2e1;
t276 = -t97 / 0.2e1;
t275 = t97 / 0.2e1;
t273 = -t114 / 0.2e1;
t272 = t114 / 0.2e1;
t271 = -t118 / 0.2e1;
t270 = t118 / 0.2e1;
t269 = -t120 / 0.2e1;
t268 = t120 / 0.2e1;
t266 = t160 / 0.2e1;
t264 = pkin(3) * t154;
t226 = t155 * t159;
t122 = -t158 * t226 + t224;
t123 = t156 * t158 + t161 * t226;
t70 = -t122 * t230 + t123 * t154;
t261 = t16 * t70;
t260 = t16 * t93;
t259 = t96 * Ifges(7,6);
t258 = t97 * Ifges(7,5);
t257 = -qJD(3) / 0.2e1;
t255 = mrSges(6,1) + mrSges(5,3);
t254 = -Ifges(5,4) - Ifges(6,6);
t17 = t154 * t164 + t230 * t55;
t247 = mrSges(7,3) * t157;
t246 = mrSges(7,3) * t160;
t245 = Ifges(4,4) * t158;
t244 = Ifges(7,4) * t157;
t243 = Ifges(7,4) * t160;
t241 = Ifges(7,5) * t157;
t240 = Ifges(7,6) * t160;
t239 = qJD(2) * pkin(2);
t238 = t110 * mrSges(6,1);
t237 = t114 * Ifges(7,3);
t236 = t120 * Ifges(5,4);
t235 = t120 * Ifges(6,6);
t229 = Ifges(4,5) * qJD(3);
t228 = Ifges(4,6) * qJD(3);
t217 = qJD(2) * t161;
t213 = qJD(6) * t157;
t212 = qJD(6) * t160;
t152 = pkin(3) * t219;
t210 = -Ifges(5,4) / 0.2e1 - Ifges(6,6) / 0.2e1;
t209 = -t50 - t223;
t208 = mrSges(4,3) * t219;
t207 = mrSges(4,3) * t217;
t204 = t230 * pkin(3);
t39 = t154 * t89 + t198;
t193 = qJ(5) * t118 + t152;
t149 = -t204 - pkin(4);
t188 = t1 * t160 - t2 * t157;
t187 = t6 * t157 + t5 * t160;
t186 = t157 * t5 - t160 * t6;
t15 = -qJD(3) * qJD(5) - t17;
t184 = mrSges(7,1) * t157 + mrSges(7,2) * t160;
t183 = Ifges(7,1) * t160 - t244;
t182 = Ifges(7,1) * t157 + t243;
t181 = -Ifges(7,2) * t157 + t243;
t180 = Ifges(7,2) * t160 + t244;
t179 = Ifges(7,5) * t160 - Ifges(7,6) * t157;
t178 = t240 + t241;
t94 = t154 * t137 - t138 * t230;
t177 = -t109 * t94 + t110 * t93;
t59 = -mrSges(7,2) * t114 + mrSges(7,3) * t96;
t60 = mrSges(7,1) * t114 - mrSges(7,3) * t97;
t176 = -t157 * t60 + t160 * t59;
t175 = -t157 * t59 - t160 * t60;
t77 = -t133 * t215 + (-qJD(3) * t156 - t199) * t220;
t174 = -t158 * t77 + t161 * t76;
t171 = -t157 * t70 + t160 * t225;
t51 = t157 * t225 + t160 * t70;
t41 = t230 * t89 - t82;
t168 = (mrSges(4,1) * t158 + mrSges(4,2) * t161) * qJD(2);
t166 = -qJD(6) * t186 + t1 * t157 + t160 * t2;
t151 = Ifges(4,4) * t217;
t147 = qJ(5) + t264;
t136 = -qJD(3) * mrSges(4,2) + t207;
t135 = qJD(3) * mrSges(4,1) - t208;
t134 = -t202 - t239;
t126 = qJD(3) * t168;
t125 = Ifges(4,1) * t219 + t151 + t229;
t124 = t228 + (t161 * Ifges(4,2) + t245) * qJD(2);
t113 = Ifges(5,4) * t118;
t112 = Ifges(6,6) * t118;
t107 = Ifges(7,3) * t110;
t106 = t110 * mrSges(6,3);
t105 = t110 * mrSges(5,2);
t98 = -t133 * t158 + t144;
t88 = qJD(3) * t122 + t191;
t87 = -qJD(3) * t123 - t158 * t199;
t80 = -pkin(4) * t286 + t173;
t75 = t120 * Ifges(5,1) + Ifges(5,5) * qJD(3) - t113;
t74 = -t118 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t236;
t73 = Ifges(6,4) * qJD(3) - t120 * Ifges(6,2) + t112;
t72 = Ifges(6,5) * qJD(3) + t118 * Ifges(6,3) - t235;
t71 = t154 * t122 + t123 * t230;
t67 = pkin(5) * t286 + t94;
t65 = -t109 * mrSges(6,2) - t106;
t64 = t109 * mrSges(5,1) + t105;
t61 = pkin(4) * t120 + t193;
t45 = pkin(4) * t119 + t169;
t44 = -t119 * pkin(5) + t63;
t42 = t120 * t274 + t193;
t40 = t154 * t87 + t230 * t88;
t38 = t154 * t88 - t230 * t87;
t29 = t237 + t258 + t259;
t28 = -qJD(3) * pkin(4) + t172;
t26 = pkin(4) * t109 + t167;
t25 = t41 - t262;
t24 = t39 - t263;
t18 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t14 = t57 * Ifges(7,1) + t58 * Ifges(7,4) + t110 * Ifges(7,5);
t13 = t57 * Ifges(7,4) + t58 * Ifges(7,2) + t110 * Ifges(7,6);
t11 = -pkin(5) * t109 - t15;
t10 = qJD(6) * t51 + t157 * t38 + t160 * t200;
t9 = qJD(6) * t171 - t157 * t200 + t160 * t38;
t8 = t157 * t24 + t160 * t42;
t7 = -t157 * t42 + t160 * t24;
t3 = [t10 * t59 + t87 * t135 + t88 * t136 + t51 * t32 - t171 * t33 + t9 * t60 - t287 * t38 + t255 * t70 * t110 + (-t122 * t161 - t123 * t158) * mrSges(4,3) * t211 + (-t109 * t255 + t18) * t71 - t209 * t40 + ((-mrSges(3,2) * t163 - t126 - t64 - t65) * t162 + (-mrSges(3,1) * t163 + (qJD(2) * (-mrSges(4,1) * t161 + mrSges(4,2) * t158) + t289) * qJD(2)) * t159) * t155 + m(7) * (-t1 * t171 + t10 * t6 + t11 * t71 + t2 * t51 + t23 * t40 + t5 * t9) + m(6) * (-t15 * t71 + t261 + t28 * t38 - t34 * t40 + (-t162 * t26 + t218 * t49) * t155) + m(5) * (t261 + t17 * t71 - t36 * t38 + t37 * t40 + (t111 * t218 - t117 * t162) * t155) + m(4) * (t122 * t77 + t123 * t76 + t87 * t98 + t88 * t99 + (t134 - t202) * t200); t295 * t60 + (t1 * t20 + t11 * t67 + t19 * t2 + t296 * t6 + t295 * t5 + (t44 - t92) * t23) * m(7) + t296 * t59 + (t259 / 0.2e1 + t258 / 0.2e1 + t237 / 0.2e1 + t29 / 0.2e1 + t75 / 0.2e1 - t73 / 0.2e1 + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t120 + t210 * t118 + t281) * t121 + (t180 * t277 + t182 * t275 + t178 * t272 + t72 / 0.2e1 - t74 / 0.2e1 + t210 * t120 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t118 - t186 * mrSges(7,3) + t297) * t119 - t287 * t291 + (-t15 * t94 + t26 * t80 + t28 * t291 - t290 * t34 + t45 * t49 + t260) * m(6) + (t111 * t153 + t117 * t150 + t17 * t94 + t290 * t37 - t291 * t36 + t260) * m(5) + (-t119 * t37 - t121 * t36 + t177) * mrSges(5,3) + (t119 * t34 + t121 * t28 + t177) * mrSges(6,1) + t174 * mrSges(4,3) + ((t135 * t158 - t136 * t161) * t162 + (-m(5) * t111 - m(6) * t49 - t289) * t159 + ((t158 * t98 - t161 * t99) * t162 + (-t239 - t134) * t159) * m(4)) * t221 + (t107 / 0.2e1 + t117 * mrSges(5,2) - t26 * mrSges(6,3) + t255 * t16 + t254 * t109 + (Ifges(5,1) + Ifges(6,2) + Ifges(7,3) / 0.2e1) * t110 + t283) * t128 + m(4) * ((-t158 * t99 - t161 * t98) * qJD(3) + t174) * pkin(8) + ((-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t121 + (-Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1) * t119 + (t125 / 0.2e1 - t98 * mrSges(4,3) - pkin(8) * t135 + t134 * mrSges(4,2) + t229 / 0.2e1) * t161 + (-t124 / 0.2e1 + pkin(3) * t78 - t99 * mrSges(4,3) - pkin(8) * t136 + t134 * mrSges(4,1) - t228 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t217) * t158) * qJD(3) + t223 * t63 + (-0.3e1 / 0.2e1 * t158 ^ 2 + 0.3e1 / 0.2e1 * t161 ^ 2) * Ifges(4,4) * t211 + t19 * t32 + t20 * t33 + t44 * t50 - (t182 * t280 + t180 * t279 + t15 * mrSges(6,1) - t17 * mrSges(5,3) + t157 * t14 / 0.2e1 + t13 * t266 + t117 * mrSges(5,1) - t26 * mrSges(6,2) - t11 * t185 + (t241 / 0.2e1 + t240 / 0.2e1 + t254) * t110 + (Ifges(6,3) + Ifges(5,2)) * t109 + t188 * mrSges(7,3) + (-mrSges(7,3) * t187 + t179 * t272 + t181 * t277 + t183 * t275 + t184 * t23 + t266 * t31 + t267 * t30) * qJD(6)) * t286 + t67 * t18 + t45 * t79 + t80 * t65 - pkin(2) * t126 + t150 * t64 + t209 * t92; (m(7) * t166 + t212 * t59 - t213 * t60 + t284) * (-pkin(9) + t149) + t287 * t39 - t223 * t41 + (t208 + t135) * t99 + (-t147 * t15 + t149 * t16 - t28 * t39 - t49 * t61 + (-qJD(5) + t41) * t34) * m(6) + (-Ifges(5,2) * t270 + Ifges(6,3) * t271 + t178 * t273 + t180 * t278 + t182 * t276 - t246 * t6 + t247 * t5 + t257 * t292 - t297) * t120 + (-t212 * t6 + t213 * t5) * mrSges(7,3) + ((t154 * t17 - t16 * t230) * pkin(3) - t111 * t152 + t36 * t39 - t37 * t41) * m(5) + t211 * Ifges(4,5) * t161 / 0.2e1 + (-t236 + t72) * t269 + (t11 * t147 - t5 * t7 - t6 * t8 + (qJD(5) - t25) * t23) * m(7) - (-Ifges(4,2) * t219 + t125 + t151) * t217 / 0.2e1 - (t114 * t178 + t180 * t96 + t182 * t97) * qJD(6) / 0.2e1 - t134 * t168 + (-mrSges(6,1) * t147 - mrSges(5,3) * t264 + t292) * t109 + (t179 / 0.2e1 - mrSges(5,3) * t204 - t293) * t110 + (-Ifges(5,1) * t269 - Ifges(7,5) * t276 + Ifges(6,2) * t268 - Ifges(7,6) * t278 - Ifges(7,3) * t273 + t257 * t293 + t281) * t118 + t294 * t16 + (t207 - t136) * t98 + t11 * t184 + (t235 + t74) * t268 - t36 * t249 - t34 * t250 + t300 * qJD(6) - t2 * t246 - t1 * t247 + t124 * t219 / 0.2e1 - t158 * t163 * (Ifges(4,1) * t161 - t245) / 0.2e1 + t288 * qJD(5) - t15 * mrSges(6,3) - t17 * mrSges(5,2) - t25 * t50 + (t112 + t73) * t271 - t8 * t59 - t7 * t60 + (-t113 + t75 + t29) * t270 - t76 * mrSges(4,2) + t77 * mrSges(4,1) - t61 * t79 + t147 * t18 - Ifges(4,6) * t197 / 0.2e1 - t78 * t152 + t149 * t238 + t37 * t248 + t28 * t251 + t14 * t266 + t13 * t267 + t181 * t279 + t183 * t280; -t157 * t32 + t160 * t33 + t105 - t106 - t294 * t109 + t175 * qJD(6) - t209 * t118 + (t175 + t287) * t120 + (-t114 * t187 + t118 * t23 + t188) * m(7) + (-t118 * t34 - t120 * t28 + t26) * m(6) + (t118 * t37 + t120 * t36 + t117) * m(5); t238 + t176 * qJD(6) - t288 * qJD(3) + (t176 + t79) * t120 + (-qJD(3) * t23 - t120 * t186 + t166) * m(7) + (qJD(3) * t34 + t120 * t49 + t16) * m(6) + t284; t107 - t23 * (mrSges(7,1) * t97 + mrSges(7,2) * t96) + (Ifges(7,1) * t96 - t265) * t276 + t30 * t275 + (Ifges(7,5) * t96 - Ifges(7,6) * t97) * t273 - t5 * t59 + t6 * t60 + (t5 * t96 + t6 * t97) * mrSges(7,3) + (-Ifges(7,2) * t97 + t31 + t95) * t278 + t283;];
tauc  = t3(:);
