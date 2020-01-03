% Calculate time derivative of joint inertia matrix for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:02
% EndTime: 2019-12-31 22:32:15
% DurationCPUTime: 5.16s
% Computational Cost: add. (7332->482), mult. (19127->728), div. (0->0), fcn. (18346->10), ass. (0->217)
t199 = sin(qJ(5));
t203 = cos(qJ(5));
t204 = cos(qJ(4));
t300 = (t199 ^ 2 + t203 ^ 2) * t204;
t200 = sin(qJ(4));
t201 = sin(qJ(3));
t205 = cos(qJ(3));
t165 = t200 * t201 - t204 * t205;
t299 = qJD(3) + qJD(4);
t127 = t299 * t165;
t166 = t200 * t205 + t201 * t204;
t239 = qJD(5) * t203;
t213 = -t199 * t127 + t166 * t239;
t240 = qJD(5) * t199;
t248 = t203 * t127;
t212 = t166 * t240 + t248;
t128 = t299 * t166;
t243 = qJD(3) * t201;
t235 = pkin(3) * t243;
t77 = pkin(4) * t128 + pkin(10) * t127 + t235;
t192 = -pkin(3) * t205 - pkin(2);
t118 = pkin(4) * t165 - pkin(10) * t166 + t192;
t288 = -pkin(9) - pkin(8);
t182 = t288 * t205;
t234 = t201 * t288;
t142 = -t204 * t182 + t200 * t234;
t80 = t118 * t203 - t142 * t199;
t141 = -t200 * t182 - t204 * t234;
t231 = qJD(3) * t288;
t175 = t201 * t231;
t224 = t205 * t231;
t90 = -qJD(4) * t141 + t204 * t175 + t200 * t224;
t25 = qJD(5) * t80 + t199 * t77 + t203 * t90;
t81 = t118 * t199 + t142 * t203;
t26 = -qJD(5) * t81 - t199 * t90 + t203 * t77;
t302 = -t26 * t199 + t203 * t25;
t241 = qJD(4) * t204;
t198 = cos(pkin(5));
t202 = sin(qJ(2));
t197 = sin(pkin(5));
t206 = cos(qJ(2));
t252 = t197 * t206;
t163 = t198 * t202 * pkin(1) + pkin(7) * t252;
t151 = pkin(8) * t198 + t163;
t152 = (-pkin(2) * t206 - pkin(8) * t202 - pkin(1)) * t197;
t104 = t205 * t151 + t201 * t152;
t253 = t197 * t202;
t157 = t198 * t205 - t201 * t253;
t87 = pkin(9) * t157 + t104;
t263 = t200 * t87;
t245 = qJD(2) * t197;
t229 = t206 * t245;
t135 = qJD(3) * t157 + t205 * t229;
t244 = qJD(2) * t202;
t230 = t197 * t244;
t154 = (pkin(2) * t202 - pkin(8) * t206) * t245;
t186 = pkin(7) * t253;
t276 = pkin(1) * t206;
t162 = t198 * t276 - t186;
t155 = t162 * qJD(2);
t66 = -qJD(3) * t104 + t205 * t154 - t155 * t201;
t49 = pkin(3) * t230 - pkin(9) * t135 + t66;
t158 = t198 * t201 + t205 * t253;
t134 = -qJD(3) * t158 - t201 * t229;
t242 = qJD(3) * t205;
t65 = -t151 * t243 + t152 * t242 + t201 * t154 + t205 * t155;
t51 = pkin(9) * t134 + t65;
t103 = -t201 * t151 + t205 * t152;
t74 = -pkin(3) * t252 - t158 * pkin(9) + t103;
t13 = -qJD(4) * t263 + t200 * t49 + t204 * t51 + t74 * t241;
t10 = pkin(10) * t230 + t13;
t273 = t200 * t74 + t204 * t87;
t37 = -pkin(10) * t252 + t273;
t113 = t157 * t200 + t158 * t204;
t150 = t186 + (-pkin(2) - t276) * t198;
t115 = -t157 * pkin(3) + t150;
t215 = t204 * t157 - t158 * t200;
t52 = -pkin(4) * t215 - t113 * pkin(10) + t115;
t18 = -t199 * t37 + t203 * t52;
t156 = t163 * qJD(2);
t102 = -t134 * pkin(3) + t156;
t60 = qJD(4) * t215 + t134 * t200 + t135 * t204;
t61 = qJD(4) * t113 - t204 * t134 + t135 * t200;
t20 = t61 * pkin(4) - t60 * pkin(10) + t102;
t2 = qJD(5) * t18 + t10 * t203 + t199 * t20;
t19 = t199 * t52 + t203 * t37;
t3 = -qJD(5) * t19 - t10 * t199 + t20 * t203;
t301 = -t3 * t199 + t2 * t203;
t281 = t199 / 0.2e1;
t279 = t203 / 0.2e1;
t226 = -t240 / 0.2e1;
t298 = Ifges(4,5) * t135 + Ifges(4,6) * t134 + Ifges(4,3) * t230;
t14 = -qJD(4) * t273 - t200 * t51 + t204 * t49;
t297 = 2 * m(5);
t296 = 2 * m(6);
t295 = -2 * mrSges(3,3);
t294 = -2 * mrSges(5,3);
t91 = qJD(4) * t142 + t200 * t175 - t204 * t224;
t293 = 0.2e1 * t91;
t292 = 0.2e1 * t141;
t291 = 0.2e1 * t156;
t214 = -t203 * t113 + t199 * t252;
t32 = qJD(5) * t214 - t199 * t60 + t203 * t230;
t290 = t32 / 0.2e1;
t95 = -t199 * t113 - t203 * t252;
t289 = t95 / 0.2e1;
t193 = Ifges(6,5) * t239;
t287 = Ifges(6,6) * t226 + t193 / 0.2e1;
t270 = Ifges(6,4) * t199;
t220 = Ifges(6,1) * t203 - t270;
t173 = t220 * qJD(5);
t286 = t173 / 0.2e1;
t285 = Ifges(6,5) * t281 + Ifges(6,6) * t279;
t269 = Ifges(6,4) * t203;
t180 = Ifges(6,1) * t199 + t269;
t283 = t180 / 0.2e1;
t282 = -t199 / 0.2e1;
t280 = t201 / 0.2e1;
t278 = t205 / 0.2e1;
t272 = Ifges(4,4) * t201;
t271 = Ifges(4,4) * t205;
t268 = Ifges(6,6) * t199;
t267 = pkin(3) * qJD(4);
t266 = t141 * t91;
t265 = t155 * mrSges(3,2);
t264 = t200 * mrSges(5,1);
t261 = t204 * mrSges(5,2);
t101 = -mrSges(5,1) * t252 - t113 * mrSges(5,3);
t56 = -mrSges(6,1) * t95 - mrSges(6,2) * t214;
t259 = -t101 + t56;
t258 = t141 * t200;
t257 = t166 * t199;
t256 = t166 * t203;
t250 = t199 * t204;
t176 = -mrSges(6,1) * t203 + mrSges(6,2) * t199;
t249 = t200 * t176;
t247 = t203 * t204;
t246 = -Ifges(5,5) * t127 - Ifges(5,6) * t128;
t238 = 0.2e1 * t197;
t31 = qJD(5) * t95 + t199 * t230 + t203 * t60;
t6 = Ifges(6,5) * t31 + Ifges(6,6) * t32 + Ifges(6,3) * t61;
t236 = Ifges(5,5) * t60 - Ifges(5,6) * t61 + Ifges(5,3) * t230;
t11 = -pkin(4) * t230 - t14;
t15 = -mrSges(6,1) * t32 + mrSges(6,2) * t31;
t233 = m(6) * t11 + t15;
t62 = mrSges(6,1) * t213 - mrSges(6,2) * t212;
t232 = m(6) * t91 + t62;
t225 = t239 / 0.2e1;
t223 = mrSges(6,3) * t300;
t222 = -mrSges(4,1) * t205 + mrSges(4,2) * t201;
t221 = mrSges(6,1) * t199 + mrSges(6,2) * t203;
t219 = -Ifges(6,2) * t199 + t269;
t43 = t204 * t74 - t263;
t216 = -t201 * t66 + t205 * t65;
t171 = t219 * qJD(5);
t178 = Ifges(6,2) * t203 + t270;
t211 = t203 * t171 + t199 * t173 - t178 * t240 + t180 * t239;
t45 = -Ifges(6,5) * t212 - Ifges(6,6) * t213 + Ifges(6,3) * t128;
t16 = mrSges(6,1) * t61 - mrSges(6,3) * t31;
t17 = -mrSges(6,2) * t61 + mrSges(6,3) * t32;
t63 = mrSges(6,2) * t215 + mrSges(6,3) * t95;
t64 = -mrSges(6,1) * t215 + mrSges(6,3) * t214;
t210 = -t63 * t240 + m(6) * (-t18 * t239 - t19 * t240 + t301) + t203 * t17 - t199 * t16 - t64 * t239;
t120 = -mrSges(6,2) * t165 - mrSges(6,3) * t257;
t121 = mrSges(6,1) * t165 - mrSges(6,3) * t256;
t69 = mrSges(6,1) * t128 + mrSges(6,3) * t212;
t70 = -mrSges(6,2) * t128 - mrSges(6,3) * t213;
t209 = -t120 * t240 + m(6) * (-t239 * t80 - t240 * t81 + t302) + t203 * t70 - t199 * t69 - t121 * t239;
t168 = t221 * qJD(5);
t36 = pkin(4) * t252 - t43;
t39 = -Ifges(6,4) * t214 + Ifges(6,2) * t95 - Ifges(6,6) * t215;
t40 = -Ifges(6,1) * t214 + Ifges(6,4) * t95 - Ifges(6,5) * t215;
t7 = Ifges(6,4) * t31 + Ifges(6,2) * t32 + Ifges(6,6) * t61;
t8 = Ifges(6,1) * t31 + Ifges(6,4) * t32 + Ifges(6,5) * t61;
t208 = t14 * mrSges(5,1) - t13 * mrSges(5,2) + t11 * t176 - t215 * t287 + t36 * t168 + t171 * t289 + t178 * t290 + t40 * t225 + t39 * t226 + t7 * t279 + t8 * t281 + t31 * t283 + t61 * t285 - t214 * t286 + t236 + ((-t18 * t203 - t19 * t199) * qJD(5) + t301) * mrSges(6,3);
t46 = -Ifges(6,4) * t212 - Ifges(6,2) * t213 + Ifges(6,6) * t128;
t47 = -Ifges(6,1) * t212 - Ifges(6,4) * t213 + Ifges(6,5) * t128;
t98 = Ifges(6,6) * t165 + t166 * t219;
t99 = Ifges(6,5) * t165 + t166 * t220;
t207 = -t248 * t283 + t128 * t285 + t141 * t168 - t171 * t257 / 0.2e1 + t256 * t286 + t165 * t287 + t47 * t281 + t46 * t279 + t99 * t225 - t90 * mrSges(5,2) + t246 + (t176 - mrSges(5,1)) * t91 - t213 * t178 / 0.2e1 + (t166 * t180 + t98) * t226 + ((-t199 * t81 - t203 * t80) * qJD(5) + t302) * mrSges(6,3);
t194 = Ifges(4,5) * t242;
t191 = -pkin(3) * t204 - pkin(4);
t190 = pkin(3) * t200 + pkin(10);
t185 = Ifges(3,5) * t229;
t181 = Ifges(4,1) * t201 + t271;
t179 = Ifges(4,2) * t205 + t272;
t174 = (Ifges(4,1) * t205 - t272) * qJD(3);
t172 = (-Ifges(4,2) * t201 + t271) * qJD(3);
t169 = (mrSges(4,1) * t201 + mrSges(4,2) * t205) * qJD(3);
t140 = -mrSges(4,1) * t252 - t158 * mrSges(4,3);
t139 = mrSges(4,2) * t252 + t157 * mrSges(4,3);
t133 = Ifges(5,1) * t166 - Ifges(5,4) * t165;
t132 = Ifges(5,4) * t166 - Ifges(5,2) * t165;
t131 = mrSges(5,1) * t165 + mrSges(5,2) * t166;
t116 = t221 * t166;
t111 = mrSges(4,1) * t230 - mrSges(4,3) * t135;
t110 = -mrSges(4,2) * t230 + mrSges(4,3) * t134;
t108 = Ifges(4,1) * t158 + Ifges(4,4) * t157 - Ifges(4,5) * t252;
t107 = Ifges(4,4) * t158 + Ifges(4,2) * t157 - Ifges(4,6) * t252;
t100 = mrSges(5,2) * t252 + mrSges(5,3) * t215;
t97 = Ifges(6,3) * t165 + (Ifges(6,5) * t203 - t268) * t166;
t88 = -mrSges(4,1) * t134 + mrSges(4,2) * t135;
t86 = -Ifges(5,1) * t127 - Ifges(5,4) * t128;
t85 = -Ifges(5,4) * t127 - Ifges(5,2) * t128;
t84 = mrSges(5,1) * t128 - mrSges(5,2) * t127;
t76 = Ifges(4,1) * t135 + Ifges(4,4) * t134 + Ifges(4,5) * t230;
t75 = Ifges(4,4) * t135 + Ifges(4,2) * t134 + Ifges(4,6) * t230;
t71 = -mrSges(5,1) * t215 + mrSges(5,2) * t113;
t68 = Ifges(5,1) * t113 + Ifges(5,4) * t215 - Ifges(5,5) * t252;
t67 = Ifges(5,4) * t113 + Ifges(5,2) * t215 - Ifges(5,6) * t252;
t55 = -mrSges(5,2) * t230 - mrSges(5,3) * t61;
t54 = mrSges(5,1) * t230 - mrSges(5,3) * t60;
t38 = -Ifges(6,5) * t214 + Ifges(6,6) * t95 - Ifges(6,3) * t215;
t23 = mrSges(5,1) * t61 + mrSges(5,2) * t60;
t22 = Ifges(5,1) * t60 - Ifges(5,4) * t61 + Ifges(5,5) * t230;
t21 = Ifges(5,4) * t60 - Ifges(5,2) * t61 + Ifges(5,6) * t230;
t1 = [(t102 * t115 + t13 * t273 + t14 * t43) * t297 + 0.2e1 * t273 * t55 - t214 * t8 + (mrSges(3,3) * t202 * t291 + (0.2e1 * t155 * mrSges(3,3) - t236 - t298) * t206 + ((t162 * t295 + Ifges(3,5) * t198 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t206) * t238) * t206 + (t163 * t295 + Ifges(4,5) * t158 + Ifges(5,5) * t113 - 0.2e1 * Ifges(3,6) * t198 + Ifges(4,6) * t157 + Ifges(5,6) * t215 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t202) * t238 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3)) * t252) * t202) * qJD(2)) * t197 - (t6 - t21) * t215 + (-mrSges(4,1) * t157 + mrSges(4,2) * t158) * t291 + (t11 * t36 + t18 * t3 + t19 * t2) * t296 + 0.2e1 * m(3) * (t155 * t163 - t156 * t162) + 0.2e1 * m(4) * (t103 * t66 + t104 * t65 + t150 * t156) + (t38 - t67) * t61 + (-0.2e1 * t156 * mrSges(3,1) + t185 - 0.2e1 * t265) * t198 + 0.2e1 * t18 * t16 + 0.2e1 * t19 * t17 + 0.2e1 * t36 * t15 + t32 * t39 + t31 * t40 + 0.2e1 * t43 * t54 + 0.2e1 * t11 * t56 + 0.2e1 * t2 * t63 + 0.2e1 * t3 * t64 + t60 * t68 + t95 * t7 + 0.2e1 * t13 * t100 + 0.2e1 * t14 * t101 + 0.2e1 * t102 * t71 + 0.2e1 * t104 * t110 + 0.2e1 * t103 * t111 + t113 * t22 + 0.2e1 * t115 * t23 + t134 * t107 + t135 * t108 + 0.2e1 * t65 * t139 + 0.2e1 * t66 * t140 + 0.2e1 * t150 * t88 + t157 * t75 + t158 * t76; t185 + ((Ifges(5,5) * t166 / 0.2e1 - Ifges(5,6) * t165 / 0.2e1 - Ifges(3,6) + Ifges(4,5) * t280 + Ifges(4,6) * t278) * t244 - (-Ifges(4,6) * t243 + t194 + t246) * t206 / 0.2e1) * t197 + m(5) * (t102 * t192 + t115 * t235 + t13 * t142 - t14 * t141 + t273 * t90 - t43 * t91) + (-t273 * mrSges(5,3) + t38 / 0.2e1 - t67 / 0.2e1) * t128 - t214 * t47 / 0.2e1 - (t45 / 0.2e1 - t85 / 0.2e1) * t215 - t265 + (-mrSges(3,1) + t222) * t156 + ((-t103 * t205 - t104 * t201) * qJD(3) + t216) * mrSges(4,3) + t46 * t289 + t98 * t290 + t75 * t278 + t76 * t280 + (-m(4) * t156 - t88) * pkin(2) + (-t13 * mrSges(5,3) + t6 / 0.2e1 - t21 / 0.2e1) * t165 + (t15 - t54) * t141 + m(6) * (t11 * t141 + t18 * t26 + t19 * t25 + t2 * t81 + t3 * t80 + t36 * t91) + (m(4) * (-t103 * t242 - t104 * t243 + t216) + t205 * t110 - t201 * t111 - t139 * t243 - t140 * t242) * pkin(8) + t259 * t91 + (t108 * t278 + (-t107 / 0.2e1 + pkin(3) * t71) * t201) * qJD(3) - (t40 * t279 + t39 * t282 - t43 * mrSges(5,3) + t68 / 0.2e1) * t127 + (t8 * t279 + t7 * t282 - t14 * mrSges(5,3) + t22 / 0.2e1 + (t40 * t282 - t203 * t39 / 0.2e1) * qJD(5)) * t166 + t36 * t62 + t25 * t63 + t26 * t64 + t18 * t69 + t19 * t70 + t80 * t16 + t81 * t17 + t31 * t99 / 0.2e1 + t90 * t100 + t113 * t86 / 0.2e1 + t115 * t84 + t11 * t116 + t2 * t120 + t3 * t121 + t102 * t131 + t60 * t133 / 0.2e1 + t142 * t55 + t150 * t169 + t157 * t172 / 0.2e1 + t158 * t174 / 0.2e1 + t134 * t179 / 0.2e1 + t135 * t181 / 0.2e1 + t192 * t23 + (t97 / 0.2e1 - t132 / 0.2e1) * t61; -0.2e1 * pkin(2) * t169 + t116 * t293 + 0.2e1 * t25 * t120 + 0.2e1 * t26 * t121 + t62 * t292 + t205 * t172 + t201 * t174 + 0.2e1 * t192 * t84 + 0.2e1 * t80 * t69 + 0.2e1 * t81 * t70 + (t205 * t181 + (0.2e1 * pkin(3) * t131 - t179) * t201) * qJD(3) + (t142 * t90 + t192 * t235 + t266) * t297 + (t25 * t81 + t26 * t80 + t266) * t296 + (t294 * t90 + t45 - t85) * t165 + (t142 * t294 - t132 + t97) * t128 - (mrSges(5,3) * t292 - t199 * t98 + t203 * t99 + t133) * t127 + (mrSges(5,3) * t293 - t199 * t46 + t203 * t47 + t86 + (-t199 * t99 - t203 * t98) * qJD(5)) * t166; t208 + (m(5) * (t13 * t200 + t14 * t204) + t204 * t54 + t200 * t55 + (t259 * t200 + (-t199 * t64 + t203 * t63 + t100) * t204 + m(6) * (-t18 * t250 + t19 * t247 + t200 * t36) + m(5) * (-t200 * t43 + t204 * t273)) * qJD(4)) * pkin(3) + t233 * t191 + t210 * t190 - t65 * mrSges(4,2) + t66 * mrSges(4,1) + t298; t194 + t207 + t232 * t191 + t209 * t190 + (m(5) * (t200 * t90 - t204 * t91) + (t127 * t204 - t128 * t200) * mrSges(5,3) + ((t166 * mrSges(5,3) + t116) * t200 + (-t165 * mrSges(5,3) + t203 * t120 - t199 * t121) * t204 + m(6) * (t247 * t81 - t250 * t80 + t258) + m(5) * (t142 * t204 + t258)) * qJD(4)) * pkin(3) + (-Ifges(4,6) * t201 + pkin(8) * t222) * qJD(3); 0.2e1 * t191 * t168 + (-0.2e1 * t261 - 0.2e1 * t264 + 0.2e1 * t249 + (t190 * t300 + t191 * t200) * t296 + 0.2e1 * t223) * t267 + t211; -pkin(4) * t233 + pkin(10) * t210 + t208; -pkin(4) * t232 + pkin(10) * t209 + t207; (-pkin(4) + t191) * t168 + (-t261 - t264 + m(6) * (-pkin(4) * t200 + pkin(10) * t300) + t249 + t223) * t267 + t211; -0.2e1 * pkin(4) * t168 + t211; mrSges(6,1) * t3 - mrSges(6,2) * t2 + t6; mrSges(6,1) * t26 - mrSges(6,2) * t25 + t45; t193 - t221 * pkin(3) * t241 + (t176 * t190 - t268) * qJD(5); t193 + (pkin(10) * t176 - t268) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
