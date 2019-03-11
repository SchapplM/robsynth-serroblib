% Calculate time derivative of joint inertia matrix for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:30
% EndTime: 2019-03-09 00:42:39
% DurationCPUTime: 4.31s
% Computational Cost: add. (7499->481), mult. (17927->710), div. (0->0), fcn. (17424->12), ass. (0->213)
t192 = sin(qJ(5));
t197 = cos(qJ(5));
t198 = cos(qJ(4));
t298 = (t192 ^ 2 + t197 ^ 2) * t198;
t172 = -mrSges(6,1) * t197 + mrSges(6,2) * t192;
t306 = -mrSges(5,1) + t172;
t245 = qJD(5) * t197;
t191 = sin(qJ(6));
t196 = cos(qJ(6));
t213 = t191 * t192 - t196 * t197;
t296 = qJD(5) + qJD(6);
t124 = t296 * t213;
t160 = t191 * t197 + t192 * t196;
t125 = t296 * t160;
t249 = -Ifges(7,5) * t124 - Ifges(7,6) * t125;
t305 = Ifges(6,5) * t245 + t249;
t190 = cos(pkin(6));
t194 = sin(qJ(3));
t199 = cos(qJ(3));
t189 = sin(pkin(6));
t195 = sin(qJ(2));
t254 = t189 * t195;
t149 = t190 * t199 - t194 * t254;
t150 = t190 * t194 + t199 * t254;
t193 = sin(qJ(4));
t102 = t149 * t193 + t150 * t198;
t200 = cos(qJ(2));
t253 = t189 * t200;
t210 = -t197 * t102 + t192 * t253;
t248 = qJD(2) * t195;
t231 = t189 * t248;
t247 = qJD(2) * t200;
t230 = t189 * t247;
t134 = -t150 * qJD(3) - t194 * t230;
t135 = t149 * qJD(3) + t199 * t230;
t214 = t198 * t149 - t150 * t193;
t49 = t214 * qJD(4) + t134 * t193 + t135 * t198;
t92 = -t192 * t102 - t197 * t253;
t27 = t92 * qJD(5) + t192 * t231 + t197 * t49;
t28 = t210 * qJD(5) - t192 * t49 + t197 * t231;
t304 = -t28 * t192 + (t192 * t210 - t197 * t92) * qJD(5) + t197 * t27;
t159 = t193 * t194 - t198 * t199;
t297 = qJD(3) + qJD(4);
t126 = t297 * t159;
t161 = t193 * t199 + t194 * t198;
t209 = -t126 * t192 + t161 * t245;
t183 = -pkin(3) * t199 - pkin(2);
t114 = pkin(4) * t159 - pkin(10) * t161 + t183;
t286 = -pkin(9) - pkin(8);
t176 = t286 * t194;
t178 = t286 * t199;
t143 = t176 * t193 - t178 * t198;
t246 = qJD(5) * t192;
t127 = t297 * t161;
t238 = pkin(3) * qJD(3) * t194;
t70 = pkin(4) * t127 + pkin(10) * t126 + t238;
t233 = qJD(3) * t286;
t170 = t194 * t233;
t223 = t199 * t233;
t299 = t198 * t176 + t178 * t193;
t86 = t299 * qJD(4) + t198 * t170 + t193 * t223;
t18 = t114 * t245 - t143 * t246 + t192 * t70 + t197 * t86;
t225 = -t192 * t86 + t197 * t70;
t136 = t197 * t143;
t73 = t192 * t114 + t136;
t19 = -t73 * qJD(5) + t225;
t303 = t18 * t197 - t19 * t192;
t180 = pkin(3) * t193 + pkin(10);
t277 = -pkin(11) - t180;
t156 = t277 * t192;
t185 = t197 * pkin(11);
t157 = t180 * t197 + t185;
t109 = t156 * t196 - t157 * t191;
t224 = qJD(5) * t277;
t273 = pkin(3) * qJD(4);
t237 = t198 * t273;
t144 = t192 * t224 + t197 * t237;
t145 = -t192 * t237 + t197 * t224;
t57 = t109 * qJD(6) + t144 * t196 + t145 * t191;
t110 = t156 * t191 + t157 * t196;
t58 = -t110 * qJD(6) - t144 * t191 + t145 * t196;
t302 = t58 * mrSges(7,1) - t57 * mrSges(7,2);
t285 = -pkin(11) - pkin(10);
t175 = t285 * t192;
t177 = pkin(10) * t197 + t185;
t140 = t175 * t196 - t177 * t191;
t232 = qJD(5) * t285;
t167 = t192 * t232;
t168 = t197 * t232;
t84 = t140 * qJD(6) + t167 * t196 + t168 * t191;
t142 = t175 * t191 + t177 * t196;
t85 = -t142 * qJD(6) - t167 * t191 + t168 * t196;
t301 = t85 * mrSges(7,1) - t84 * mrSges(7,2);
t106 = t213 * t161;
t261 = t126 * t197;
t300 = -Ifges(6,5) * t261 + Ifges(6,3) * t127;
t295 = 2 * m(6);
t294 = 2 * m(7);
t293 = -2 * mrSges(5,3);
t74 = mrSges(7,1) * t125 - mrSges(7,2) * t124;
t292 = 0.2e1 * t74;
t87 = t143 * qJD(4) + t170 * t193 - t198 * t223;
t291 = 0.2e1 * t87;
t130 = mrSges(7,1) * t213 + mrSges(7,2) * t160;
t290 = 0.2e1 * t130;
t289 = -0.2e1 * t299;
t288 = m(5) / 0.2e1;
t287 = m(6) * pkin(4);
t280 = pkin(3) * t198;
t181 = -pkin(4) - t280;
t281 = m(6) * t181;
t276 = Ifges(6,4) * t192;
t275 = Ifges(6,4) * t197;
t274 = Ifges(6,6) * t192;
t272 = pkin(5) * qJD(6);
t50 = t102 * qJD(4) - t198 * t134 + t135 * t193;
t33 = t214 * t50;
t271 = t125 * mrSges(7,3);
t270 = t299 * t87;
t269 = t213 * mrSges(7,3);
t266 = t193 * mrSges(5,1);
t264 = t198 * mrSges(5,2);
t263 = t214 * t193;
t260 = t299 * t193;
t259 = t161 * t192;
t258 = t161 * t197;
t252 = t192 * t198;
t251 = t193 * t172;
t250 = t197 * t198;
t244 = qJD(6) * t191;
t243 = qJD(6) * t196;
t242 = 0.2e1 * mrSges(7,3);
t241 = 0.2e1 * t194;
t240 = mrSges(7,3) * t272;
t31 = -t125 * t161 + t213 * t126;
t32 = t296 * t106 + t160 * t126;
t239 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t127;
t236 = pkin(5) * t246;
t235 = t196 * t124 * mrSges(7,3);
t44 = t191 * t210 + t196 * t92;
t7 = t44 * qJD(6) + t191 * t28 + t196 * t27;
t45 = t191 * t92 - t196 * t210;
t8 = -t45 * qJD(6) - t191 * t27 + t196 * t28;
t234 = t8 * mrSges(7,1) - t7 * mrSges(7,2);
t182 = -pkin(5) * t197 - pkin(4);
t229 = t161 * t246;
t227 = -t246 / 0.2e1;
t226 = -(2 * Ifges(5,4)) - t274;
t72 = t197 * t114 - t143 * t192;
t222 = t189 ^ 2 * t195 * t247;
t221 = mrSges(6,3) * t298;
t220 = -mrSges(4,1) * t199 + mrSges(4,2) * t194;
t219 = mrSges(6,1) * t192 + mrSges(6,2) * t197;
t218 = Ifges(6,1) * t197 - t276;
t217 = -Ifges(6,2) * t192 + t275;
t216 = Ifges(6,5) * t192 + Ifges(6,6) * t197;
t215 = -t214 * t87 - t299 * t50;
t55 = pkin(5) * t159 - pkin(11) * t258 + t72;
t63 = -pkin(11) * t259 + t73;
t20 = -t191 * t63 + t196 * t55;
t21 = t191 * t55 + t196 * t63;
t14 = pkin(11) * t261 + pkin(5) * t127 + (-t136 + (pkin(11) * t161 - t114) * t192) * qJD(5) + t225;
t15 = -t209 * pkin(11) + t18;
t3 = t20 * qJD(6) + t14 * t191 + t15 * t196;
t4 = -t21 * qJD(6) + t14 * t196 - t15 * t191;
t212 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t239;
t211 = -t196 * t213 * t240 + (-pkin(5) * t271 + t160 * t240) * t191 + t305;
t208 = t229 + t261;
t132 = Ifges(7,4) * t160 - Ifges(7,2) * t213;
t133 = Ifges(7,1) * t160 - Ifges(7,4) * t213;
t165 = t217 * qJD(5);
t166 = t218 * qJD(5);
t173 = Ifges(6,2) * t197 + t276;
t174 = Ifges(6,1) * t192 + t275;
t76 = -Ifges(7,4) * t124 - Ifges(7,2) * t125;
t77 = -Ifges(7,1) * t124 - Ifges(7,4) * t125;
t206 = -t124 * t133 - t125 * t132 + t160 * t77 + t197 * t165 + t192 * t166 - t173 * t246 + t174 * t245 - t213 * t76;
t205 = -t134 * t194 + t135 * t199 + (-t149 * t199 - t150 * t194) * qJD(3);
t204 = m(6) * t304;
t116 = -mrSges(6,2) * t159 - mrSges(6,3) * t259;
t117 = mrSges(6,1) * t159 - mrSges(6,3) * t258;
t61 = mrSges(6,1) * t127 + t208 * mrSges(6,3);
t62 = -mrSges(6,2) * t127 - t209 * mrSges(6,3);
t203 = -t117 * t245 - t116 * t246 + m(6) * (-t72 * t245 - t73 * t246 + t303) + t197 * t62 - t192 * t61;
t163 = t219 * qJD(5);
t202 = -t45 * t271 - t7 * t269 + (t124 * t44 - t160 * t8) * mrSges(7,3) - t49 * mrSges(5,2) + (-t74 - t163) * t214 + t304 * mrSges(6,3) + (t130 + t306) * t50;
t105 = t160 * t161;
t11 = Ifges(7,4) * t31 + Ifges(7,2) * t32 + Ifges(7,6) * t127;
t12 = Ifges(7,1) * t31 + Ifges(7,4) * t32 + Ifges(7,5) * t127;
t38 = -t208 * Ifges(6,4) - t209 * Ifges(6,2) + Ifges(6,6) * t127;
t39 = -t208 * Ifges(6,1) - t209 * Ifges(6,4) + Ifges(6,5) * t127;
t46 = t209 * pkin(5) + t87;
t59 = -Ifges(7,4) * t106 - Ifges(7,2) * t105 + Ifges(7,6) * t159;
t60 = -Ifges(7,1) * t106 - Ifges(7,4) * t105 + Ifges(7,5) * t159;
t95 = t159 * Ifges(6,6) + t217 * t161;
t96 = t159 * Ifges(6,5) + t218 * t161;
t99 = pkin(5) * t259 - t299;
t201 = (Ifges(7,5) * t160 - Ifges(7,6) * t213 + t216) * t127 / 0.2e1 - t213 * t11 / 0.2e1 + (t161 * t227 - t261 / 0.2e1) * t174 + (t124 * t20 - t160 * t4) * mrSges(7,3) - t299 * t163 - t3 * t269 - t21 * t271 + t166 * t258 / 0.2e1 - t165 * t259 / 0.2e1 + t96 * t245 / 0.2e1 + t306 * t87 - t209 * t173 / 0.2e1 + (-Ifges(6,6) * t246 + t305) * t159 / 0.2e1 + ((-t192 * t73 - t197 * t72) * qJD(5) + t303) * mrSges(6,3) + t197 * t38 / 0.2e1 + t192 * t39 / 0.2e1 + t95 * t227 - t86 * mrSges(5,2) + t99 * t74 - t105 * t76 / 0.2e1 - t106 * t77 / 0.2e1 - t124 * t60 / 0.2e1 - t125 * t59 / 0.2e1 - Ifges(5,5) * t126 - Ifges(5,6) * t127 + t46 * t130 + t32 * t132 / 0.2e1 + t31 * t133 / 0.2e1 + t160 * t12 / 0.2e1;
t171 = t182 - t280;
t169 = t193 * t273 + t236;
t164 = (mrSges(4,1) * t194 + mrSges(4,2) * t199) * qJD(3);
t155 = (-mrSges(7,1) * t191 - mrSges(7,2) * t196) * t272;
t131 = mrSges(5,1) * t159 + mrSges(5,2) * t161;
t112 = t219 * t161;
t91 = mrSges(7,1) * t159 + mrSges(7,3) * t106;
t90 = -mrSges(7,2) * t159 - mrSges(7,3) * t105;
t75 = mrSges(5,1) * t127 - mrSges(5,2) * t126;
t64 = mrSges(7,1) * t105 - mrSges(7,2) * t106;
t51 = t209 * mrSges(6,1) - t208 * mrSges(6,2);
t23 = -mrSges(7,2) * t127 + mrSges(7,3) * t32;
t22 = mrSges(7,1) * t127 - mrSges(7,3) * t31;
t13 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t1 = [0.2e1 * m(7) * (t44 * t8 + t45 * t7 - t33) + 0.2e1 * m(6) * (-t210 * t27 + t28 * t92 - t33) + 0.2e1 * m(5) * (t102 * t49 - t222 - t33) + 0.2e1 * m(4) * (t149 * t134 + t150 * t135 - t222); t27 * t116 + t28 * t117 + t44 * t22 + t45 * t23 + t92 * t61 - t210 * t62 + t7 * t90 + t8 * t91 + (t112 + t64) * t50 - (t13 + t51) * t214 + (-t102 * t127 + t126 * t214 - t159 * t49 + t161 * t50) * mrSges(5,3) + t205 * mrSges(4,3) + ((-t164 - t75) * t200 + (-t200 * mrSges(3,2) + (-mrSges(3,1) + t131 + t220) * t195) * qJD(2)) * t189 + m(5) * (t86 * t102 + t143 * t49 + (t183 * t248 - t200 * t238) * t189 + t215) + m(6) * (-t18 * t210 + t19 * t92 + t27 * t73 + t28 * t72 + t215) + m(7) * (t20 * t8 + t21 * t7 - t214 * t46 + t3 * t45 + t4 * t44 + t50 * t99) + (-pkin(2) * t231 + t205 * pkin(8)) * m(4); (t86 * t293 - t226 * t126 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t127 + t239 + t300) * t159 + (mrSges(5,3) * t291 - 0.2e1 * Ifges(5,1) * t126 - t192 * t38 + t197 * t39 + (Ifges(6,5) * t197 + t226) * t127 + (-t159 * t216 - t192 * t96 - t197 * t95) * qJD(5)) * t161 + (t18 * t73 + t19 * t72 - t270) * t295 + 0.2e1 * m(5) * (t143 * t86 + t183 * t238 - t270) + ((0.2e1 * Ifges(4,4) * t199 + (Ifges(4,1) - Ifges(4,2)) * t241) * t199 + (-Ifges(4,4) * t194 + pkin(3) * t131) * t241) * qJD(3) + t51 * t289 + t112 * t291 + (t20 * t4 + t21 * t3 + t46 * t99) * t294 - (mrSges(5,3) * t289 - t192 * t95 + t197 * t96) * t126 + (-Ifges(7,5) * t106 - Ifges(7,6) * t105 + t143 * t293) * t127 + 0.2e1 * t20 * t22 + 0.2e1 * t21 * t23 + t32 * t59 + t31 * t60 + 0.2e1 * t46 * t64 + 0.2e1 * t72 * t61 + 0.2e1 * t73 * t62 + 0.2e1 * t3 * t90 + 0.2e1 * t4 * t91 + 0.2e1 * t99 * t13 - t105 * t11 - t106 * t12 + 0.2e1 * t18 * t116 + 0.2e1 * t19 * t117 - 0.2e1 * pkin(2) * t164 + 0.2e1 * t183 * t75; t202 + t180 * t204 + t50 * t281 + 0.2e1 * ((t193 * t49 - t198 * t50) * t288 + (m(6) * (-t210 * t250 - t92 * t252 - t263) / 0.2e1 + (t102 * t198 - t263) * t288) * qJD(4)) * pkin(3) + m(7) * (t109 * t8 + t110 * t7 - t169 * t214 + t171 * t50 + t44 * t58 + t45 * t57) + t134 * mrSges(4,1) - t135 * mrSges(4,2); (Ifges(4,5) * t199 - Ifges(4,6) * t194 + t220 * pkin(8)) * qJD(3) + m(7) * (t109 * t4 + t110 * t3 + t169 * t99 + t171 * t46 + t20 * t58 + t21 * t57) + t201 + (m(5) * (t193 * t86 - t198 * t87) + (t126 * t198 - t127 * t193) * mrSges(5,3) + ((t161 * mrSges(5,3) + t112) * t193 + (-t159 * mrSges(5,3) + t197 * t116 - t192 * t117) * t198 + m(6) * (t73 * t250 - t72 * t252 - t260) + m(5) * (t143 * t198 - t260)) * qJD(4)) * pkin(3) + t203 * t180 + t87 * t281 + t57 * t90 + t58 * t91 + t109 * t22 + t110 * t23 + t169 * t64 + t171 * t13 + t181 * t51; (t109 * t58 + t110 * t57 + t169 * t171) * t294 + t169 * t290 + t171 * t292 + 0.2e1 * t181 * t163 + (t109 * t124 - t110 * t125 - t58 * t160 - t213 * t57) * t242 + (-0.2e1 * t264 - 0.2e1 * t266 + 0.2e1 * t251 + (t298 * t180 + t181 * t193) * t295 + 0.2e1 * t221) * t273 + t206; t202 + pkin(10) * t204 - t50 * t287 + m(7) * (t140 * t8 + t142 * t7 + t182 * t50 - t214 * t236 + t44 * t85 + t45 * t84); t64 * t236 + t201 + m(7) * (t140 * t4 + t142 * t3 + t182 * t46 + t20 * t85 + t21 * t84 + t99 * t236) - t87 * t287 + t203 * pkin(10) - pkin(4) * t51 + t84 * t90 + t85 * t91 + t140 * t22 + t142 * t23 + t182 * t13; m(7) * (t109 * t85 + t110 * t84 + t140 * t58 + t142 * t57 + t169 * t182 + t171 * t236) + (t171 + t182) * t74 + (-pkin(4) + t181) * t163 + (t169 + t236) * t130 + (-t264 + t251 - t266 + m(6) * (-pkin(4) * t193 + t298 * pkin(10)) + t221) * t273 + ((-t58 - t85) * t160 - (t57 + t84) * t213 - (t110 + t142) * t125 - (-t109 - t140) * t124) * mrSges(7,3) + t206; -0.2e1 * pkin(4) * t163 + t236 * t290 + t182 * t292 + (t140 * t85 + t142 * t84 + t182 * t236) * t294 + (t140 * t124 - t142 * t125 - t85 * t160 - t213 * t84) * t242 + t206; t28 * mrSges(6,1) - t27 * mrSges(6,2) + m(7) * (t191 * t7 + t196 * t8 + (-t191 * t44 + t196 * t45) * qJD(6)) * pkin(5) + t234; -Ifges(6,5) * t229 + t19 * mrSges(6,1) - t18 * mrSges(6,2) - t209 * Ifges(6,6) + (-t91 * t244 + t196 * t22 + m(7) * (t191 * t3 + t196 * t4 - t20 * t244 + t21 * t243) + t90 * t243 + t191 * t23) * pkin(5) + t212 + t300; -t219 * t237 + (t172 * t180 - t274) * qJD(5) + (m(7) * (-t109 * t244 + t110 * t243 + t191 * t57 + t196 * t58) + t235) * pkin(5) + t211 + t302; (pkin(10) * t172 - t274) * qJD(5) + (m(7) * (-t140 * t244 + t142 * t243 + t191 * t84 + t196 * t85) + t235) * pkin(5) + t211 + t301; 0.2e1 * t155; t234; t212; t249 + t302; t249 + t301; t155; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
