% Calculate time derivative of joint inertia matrix for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:11
% EndTime: 2019-03-09 18:25:27
% DurationCPUTime: 6.43s
% Computational Cost: add. (13911->564), mult. (32908->844), div. (0->0), fcn. (31448->10), ass. (0->235)
t297 = -Ifges(4,3) - Ifges(5,3);
t236 = sin(qJ(3));
t237 = sin(qJ(2));
t240 = cos(qJ(3));
t269 = qJD(3) * t240;
t241 = cos(qJ(2));
t272 = qJD(2) * t241;
t245 = t236 * t272 + t237 * t269;
t284 = -qJ(4) - pkin(8);
t211 = t284 * t236;
t212 = t284 * t240;
t232 = sin(pkin(11));
t233 = cos(pkin(11));
t162 = t233 * t211 + t212 * t232;
t197 = t232 * t240 + t233 * t236;
t142 = -pkin(9) * t197 + t162;
t163 = t232 * t211 - t233 * t212;
t248 = t232 * t236 - t233 * t240;
t143 = -pkin(9) * t248 + t163;
t235 = sin(qJ(5));
t239 = cos(qJ(5));
t91 = t235 * t142 + t239 * t143;
t222 = pkin(3) * t233 + pkin(4);
t287 = pkin(3) * t232;
t186 = t239 * t222 - t235 * t287;
t210 = -pkin(2) * t241 - t237 * pkin(8) - pkin(1);
t276 = t240 * t241;
t221 = pkin(7) * t276;
t171 = t236 * t210 + t221;
t273 = qJD(2) * t237;
t180 = t197 * t237;
t181 = t248 * t237;
t131 = -t180 * t239 + t181 * t235;
t270 = qJD(3) * t237;
t138 = -t197 * t272 + t248 * t270;
t190 = t197 * qJD(3);
t139 = -t190 * t237 - t248 * t272;
t68 = qJD(5) * t131 + t138 * t235 + t139 * t239;
t132 = -t180 * t235 - t181 * t239;
t69 = -qJD(5) * t132 + t138 * t239 - t139 * t235;
t296 = -Ifges(6,5) * t68 - Ifges(6,6) * t69 - Ifges(6,3) * t273;
t260 = t240 * t272;
t295 = -Ifges(4,5) * t260 - Ifges(5,5) * t139 - Ifges(5,6) * t138 + t297 * t273;
t294 = 2 * m(4);
t293 = 2 * m(5);
t292 = 2 * m(6);
t291 = 2 * m(7);
t290 = -0.2e1 * pkin(1);
t289 = 0.2e1 * pkin(7);
t288 = -t236 / 0.2e1;
t286 = pkin(7) * t236;
t199 = t240 * t210;
t278 = t237 * t240;
t155 = -qJ(4) * t278 + t199 + (-pkin(3) - t286) * t241;
t279 = t236 * t237;
t161 = -qJ(4) * t279 + t171;
t110 = t233 * t155 - t232 * t161;
t88 = -pkin(4) * t241 + t181 * pkin(9) + t110;
t111 = t232 * t155 + t233 * t161;
t92 = -pkin(9) * t180 + t111;
t54 = t235 * t88 + t239 * t92;
t283 = Ifges(4,4) * t236;
t282 = Ifges(4,4) * t240;
t281 = Ifges(4,6) * t236;
t280 = t241 * Ifges(4,6);
t268 = qJD(4) * t240;
t208 = (pkin(2) * t237 - pkin(8) * t241) * qJD(2);
t274 = t240 * t208 + t273 * t286;
t100 = -t237 * t268 + (pkin(3) * t237 - qJ(4) * t276) * qJD(2) + (-t221 + (qJ(4) * t237 - t210) * t236) * qJD(3) + t274;
t275 = t236 * t208 + t210 * t269;
t109 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t278 + (-qJD(4) * t237 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t241) * t236 + t275;
t56 = t232 * t100 + t233 * t109;
t214 = Ifges(4,1) * t236 + t282;
t277 = t240 * t214;
t255 = qJD(3) * t284;
t188 = t236 * t255 + t268;
t189 = -qJD(4) * t236 + t240 * t255;
t141 = t233 * t188 + t232 * t189;
t209 = pkin(3) * t279 + t237 * pkin(7);
t271 = qJD(3) * t236;
t267 = qJD(5) * t235;
t266 = qJD(5) * t239;
t234 = sin(qJ(6));
t265 = qJD(6) * t234;
t238 = cos(qJ(6));
t264 = qJD(6) * t238;
t78 = t131 * t238 - t132 * t234;
t25 = qJD(6) * t78 + t234 * t69 + t238 * t68;
t79 = t131 * t234 + t132 * t238;
t26 = -qJD(6) * t79 - t234 * t68 + t238 * t69;
t262 = -Ifges(7,5) * t25 - Ifges(7,6) * t26 - Ifges(7,3) * t273;
t230 = pkin(7) * t272;
t229 = pkin(3) * t271;
t169 = pkin(3) * t245 + t230;
t223 = -pkin(3) * t240 - pkin(2);
t259 = t236 * t270;
t31 = -t69 * mrSges(6,1) + t68 * mrSges(6,2);
t9 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t153 = -t197 * t235 - t239 * t248;
t191 = t248 * qJD(3);
t101 = qJD(5) * t153 - t190 * t235 - t191 * t239;
t154 = t197 * t239 - t235 * t248;
t102 = -qJD(5) * t154 - t190 * t239 + t191 * t235;
t103 = t153 * t238 - t154 * t234;
t38 = qJD(6) * t103 + t101 * t238 + t102 * t234;
t104 = t153 * t234 + t154 * t238;
t39 = -qJD(6) * t104 - t101 * t234 + t102 * t238;
t17 = -t39 * mrSges(7,1) + t38 * mrSges(7,2);
t182 = pkin(5) + t186;
t187 = t222 * t235 + t239 * t287;
t136 = t182 * t238 - t187 * t234;
t176 = t186 * qJD(5);
t177 = t187 * qJD(5);
t86 = qJD(6) * t136 + t176 * t238 - t177 * t234;
t137 = t182 * t234 + t187 * t238;
t87 = -qJD(6) * t137 - t176 * t234 - t177 * t238;
t257 = t87 * mrSges(7,1) - t86 * mrSges(7,2);
t57 = -t102 * mrSges(6,1) + t101 * mrSges(6,2);
t256 = (2 * Ifges(3,4)) + t281;
t166 = pkin(4) * t190 + t229;
t53 = -t235 * t92 + t239 * t88;
t93 = -t138 * mrSges(5,1) + t139 * mrSges(5,2);
t147 = t190 * mrSges(5,1) - t191 * mrSges(5,2);
t55 = t233 * t100 - t109 * t232;
t90 = t239 * t142 - t143 * t235;
t140 = -t188 * t232 + t233 * t189;
t156 = pkin(4) * t180 + t209;
t36 = Ifges(7,6) * t39;
t37 = Ifges(7,5) * t38;
t116 = pkin(9) * t191 + t140;
t117 = -pkin(9) * t190 + t141;
t42 = t235 * t116 + t239 * t117 + t142 * t266 - t143 * t267;
t27 = pkin(10) * t102 + t42;
t43 = -qJD(5) * t91 + t239 * t116 - t117 * t235;
t28 = -pkin(10) * t101 + t43;
t70 = -pkin(10) * t154 + t90;
t71 = pkin(10) * t153 + t91;
t32 = -t234 * t71 + t238 * t70;
t5 = qJD(6) * t32 + t234 * t28 + t238 * t27;
t33 = t234 * t70 + t238 * t71;
t6 = -qJD(6) * t33 - t234 * t27 + t238 * t28;
t254 = t6 * mrSges(7,1) - t5 * mrSges(7,2) + t36 + t37;
t253 = -mrSges(4,1) * t240 + mrSges(4,2) * t236;
t252 = mrSges(4,1) * t236 + mrSges(4,2) * t240;
t251 = Ifges(4,1) * t240 - t283;
t250 = -Ifges(4,2) * t236 + t282;
t213 = Ifges(4,2) * t240 + t283;
t249 = Ifges(4,5) * t236 + Ifges(4,6) * t240;
t34 = -pkin(5) * t241 - t132 * pkin(10) + t53;
t40 = pkin(10) * t131 + t54;
t15 = -t234 * t40 + t238 * t34;
t16 = t234 * t34 + t238 * t40;
t172 = pkin(4) * t248 + t223;
t112 = -pkin(4) * t138 + t169;
t46 = pkin(4) * t273 - pkin(9) * t139 + t55;
t51 = pkin(9) * t138 + t56;
t14 = -qJD(5) * t54 - t235 * t51 + t239 * t46;
t10 = pkin(5) * t273 - pkin(10) * t68 + t14;
t13 = t235 * t46 + t239 * t51 + t88 * t266 - t267 * t92;
t11 = pkin(10) * t69 + t13;
t2 = qJD(6) * t15 + t10 * t234 + t11 * t238;
t3 = -qJD(6) * t16 + t10 * t238 - t11 * t234;
t247 = t3 * mrSges(7,1) - t2 * mrSges(7,2) - t262;
t246 = -t259 + t260;
t244 = -t177 * mrSges(6,1) - t176 * mrSges(6,2) + t257;
t98 = Ifges(6,6) * t102;
t99 = Ifges(6,5) * t101;
t243 = t43 * mrSges(6,1) - t42 * mrSges(6,2) + t254 + t98 + t99;
t242 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t247 - t296;
t228 = Ifges(4,5) * t269;
t207 = -mrSges(4,1) * t241 - mrSges(4,3) * t278;
t206 = mrSges(4,2) * t241 - mrSges(4,3) * t279;
t205 = t251 * qJD(3);
t204 = t250 * qJD(3);
t203 = t252 * qJD(3);
t195 = (-mrSges(7,1) * t234 - mrSges(7,2) * t238) * qJD(6) * pkin(5);
t185 = Ifges(5,5) * t191;
t184 = Ifges(5,6) * t190;
t179 = -Ifges(4,5) * t241 + t237 * t251;
t178 = t237 * t250 - t280;
t170 = -t241 * t286 + t199;
t168 = -mrSges(4,2) * t273 - mrSges(4,3) * t245;
t167 = mrSges(4,1) * t273 - mrSges(4,3) * t246;
t165 = -mrSges(5,1) * t241 + t181 * mrSges(5,3);
t164 = mrSges(5,2) * t241 - t180 * mrSges(5,3);
t160 = Ifges(5,1) * t197 - Ifges(5,4) * t248;
t159 = Ifges(5,4) * t197 - Ifges(5,2) * t248;
t158 = mrSges(5,1) * t248 + mrSges(5,2) * t197;
t152 = mrSges(4,1) * t245 + mrSges(4,2) * t246;
t149 = -Ifges(5,1) * t191 - Ifges(5,4) * t190;
t148 = -Ifges(5,4) * t191 - Ifges(5,2) * t190;
t146 = -t214 * t270 + (Ifges(4,5) * t237 + t241 * t251) * qJD(2);
t145 = -t213 * t270 + (Ifges(4,6) * t237 + t241 * t250) * qJD(2);
t144 = mrSges(5,1) * t180 - mrSges(5,2) * t181;
t127 = -Ifges(5,1) * t181 - Ifges(5,4) * t180 - Ifges(5,5) * t241;
t126 = -Ifges(5,4) * t181 - Ifges(5,2) * t180 - Ifges(5,6) * t241;
t124 = -qJD(3) * t171 + t274;
t123 = (-t240 * t273 - t241 * t271) * pkin(7) + t275;
t122 = mrSges(5,1) * t273 - mrSges(5,3) * t139;
t121 = -mrSges(5,2) * t273 + mrSges(5,3) * t138;
t120 = -mrSges(6,1) * t241 - t132 * mrSges(6,3);
t119 = mrSges(6,2) * t241 + t131 * mrSges(6,3);
t118 = -pkin(5) * t153 + t172;
t108 = Ifges(6,1) * t154 + Ifges(6,4) * t153;
t107 = Ifges(6,4) * t154 + Ifges(6,2) * t153;
t106 = -mrSges(6,1) * t153 + mrSges(6,2) * t154;
t96 = -pkin(5) * t131 + t156;
t85 = -mrSges(6,1) * t131 + mrSges(6,2) * t132;
t84 = Ifges(5,1) * t139 + Ifges(5,4) * t138 + Ifges(5,5) * t273;
t83 = Ifges(5,4) * t139 + Ifges(5,2) * t138 + Ifges(5,6) * t273;
t76 = Ifges(6,1) * t132 + Ifges(6,4) * t131 - Ifges(6,5) * t241;
t75 = Ifges(6,4) * t132 + Ifges(6,2) * t131 - Ifges(6,6) * t241;
t74 = -pkin(5) * t102 + t166;
t73 = -mrSges(7,1) * t241 - t79 * mrSges(7,3);
t72 = mrSges(7,2) * t241 + t78 * mrSges(7,3);
t64 = -mrSges(6,2) * t273 + mrSges(6,3) * t69;
t63 = mrSges(6,1) * t273 - mrSges(6,3) * t68;
t62 = Ifges(7,1) * t104 + Ifges(7,4) * t103;
t61 = Ifges(7,4) * t104 + Ifges(7,2) * t103;
t60 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t59 = Ifges(6,1) * t101 + Ifges(6,4) * t102;
t58 = Ifges(6,4) * t101 + Ifges(6,2) * t102;
t52 = -mrSges(7,1) * t78 + mrSges(7,2) * t79;
t50 = Ifges(7,1) * t79 + Ifges(7,4) * t78 - Ifges(7,5) * t241;
t49 = Ifges(7,4) * t79 + Ifges(7,2) * t78 - Ifges(7,6) * t241;
t47 = -pkin(5) * t69 + t112;
t30 = Ifges(6,1) * t68 + Ifges(6,4) * t69 + Ifges(6,5) * t273;
t29 = Ifges(6,4) * t68 + Ifges(6,2) * t69 + Ifges(6,6) * t273;
t21 = -mrSges(7,2) * t273 + mrSges(7,3) * t26;
t20 = mrSges(7,1) * t273 - mrSges(7,3) * t25;
t19 = Ifges(7,1) * t38 + Ifges(7,4) * t39;
t18 = Ifges(7,4) * t38 + Ifges(7,2) * t39;
t8 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t273;
t7 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t273;
t1 = [0.2e1 * t123 * t206 + 0.2e1 * t124 * t207 + 0.2e1 * t209 * t93 - t180 * t83 - t181 * t84 + 0.2e1 * t56 * t164 + 0.2e1 * t55 * t165 + 0.2e1 * t169 * t144 + 0.2e1 * t170 * t167 + 0.2e1 * t171 * t168 + 0.2e1 * t156 * t31 + t138 * t126 + t139 * t127 + t131 * t29 + t132 * t30 + 0.2e1 * t13 * t119 + 0.2e1 * t14 * t120 + 0.2e1 * t111 * t121 + 0.2e1 * t110 * t122 + 0.2e1 * t112 * t85 + 0.2e1 * t96 * t9 + t69 * t75 + t68 * t76 + t78 * t7 + t79 * t8 + 0.2e1 * t2 * t72 + 0.2e1 * t3 * t73 + 0.2e1 * t53 * t63 + 0.2e1 * t54 * t64 + 0.2e1 * t47 * t52 + t26 * t49 + t25 * t50 + (t152 * t289 - t236 * t145 + t240 * t146 + (-t240 * t178 - t236 * t179 + t241 * t249) * qJD(3) + (mrSges(3,1) * t290 + Ifges(7,5) * t79 + Ifges(7,6) * t78 + Ifges(6,5) * t132 + Ifges(6,6) * t131 - Ifges(5,5) * t181 - Ifges(5,6) * t180 + (Ifges(4,5) * t240 - t256) * t237 + (pkin(7) ^ 2 * t294 + t252 * t289 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(6,3) - Ifges(7,3) + t297) * t241) * qJD(2)) * t237 + 0.2e1 * t16 * t21 + 0.2e1 * t15 * t20 + (t15 * t3 + t16 * t2 + t47 * t96) * t291 + (t112 * t156 + t13 * t54 + t14 * t53) * t292 + (t110 * t55 + t111 * t56 + t169 * t209) * t293 + (t171 * t123 + t170 * t124) * t294 + ((mrSges(3,2) * t290 - t236 * t178 + t240 * t179 + t241 * t256) * qJD(2) + t262 + t295 + t296) * t241; ((t123 * t240 - t124 * t236 + (-t170 * t240 - t171 * t236) * qJD(3)) * pkin(8) - pkin(2) * t230) * m(4) + (-t101 * t53 + t102 * t54 + t13 * t153 - t14 * t154) * mrSges(6,3) + m(7) * (t118 * t47 + t15 * t6 + t16 * t5 + t2 * t33 + t3 * t32 + t74 * t96) + m(6) * (t112 * t172 + t13 * t91 + t14 * t90 + t156 * t166 + t42 * t54 + t43 * t53) + (t145 / 0.2e1 + t123 * mrSges(4,3) + pkin(8) * t168) * t240 + (t146 / 0.2e1 - t124 * mrSges(4,3) - pkin(8) * t167) * t236 + (t110 * t191 - t111 * t190 - t197 * t55 - t248 * t56) * mrSges(5,3) + (-Ifges(3,6) * qJD(2) + t204 * t288 + t240 * t205 / 0.2e1 + (t214 * t288 - t240 * t213 / 0.2e1) * qJD(3) + (mrSges(3,2) * qJD(2) + t203) * pkin(7) + (Ifges(5,5) * t197 + Ifges(6,5) * t154 + Ifges(7,5) * t104 - Ifges(5,6) * t248 + Ifges(6,6) * t153 + Ifges(7,6) * t103 + t249) * qJD(2) / 0.2e1) * t237 - t248 * t83 / 0.2e1 + (t185 / 0.2e1 + t184 / 0.2e1 - t99 / 0.2e1 - t98 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 - t228 / 0.2e1 + (Ifges(3,5) + t213 * t288 + t277 / 0.2e1 + (-mrSges(3,1) + t253) * pkin(7)) * qJD(2)) * t241 + ((-pkin(8) * t207 - t170 * mrSges(4,3) + t179 / 0.2e1) * t240 + (t280 / 0.2e1 - pkin(8) * t206 - t171 * mrSges(4,3) + pkin(3) * t144 - t178 / 0.2e1) * t236) * qJD(3) + m(5) * (t110 * t140 + t111 * t141 + t162 * t55 + t163 * t56 + t169 * t223 + t209 * t229) + t223 * t93 + t209 * t147 + t197 * t84 / 0.2e1 - t190 * t126 / 0.2e1 - t191 * t127 / 0.2e1 - t180 * t148 / 0.2e1 - t181 * t149 / 0.2e1 + t162 * t122 + t163 * t121 + t141 * t164 + t140 * t165 + t166 * t85 + t169 * t158 + t172 * t31 - pkin(2) * t152 + t153 * t29 / 0.2e1 + t154 * t30 / 0.2e1 + t156 * t57 + t138 * t159 / 0.2e1 + t139 * t160 / 0.2e1 + t131 * t58 / 0.2e1 + t132 * t59 / 0.2e1 + t118 * t9 + t42 * t119 + t43 * t120 + t112 * t106 + t101 * t76 / 0.2e1 + t102 * t75 / 0.2e1 + t103 * t7 / 0.2e1 + t104 * t8 / 0.2e1 + t69 * t107 / 0.2e1 + t68 * t108 / 0.2e1 + t96 * t17 + t90 * t63 + t91 * t64 + t78 * t18 / 0.2e1 + t79 * t19 / 0.2e1 + t5 * t72 + t6 * t73 + t74 * t52 + t47 * t60 + t26 * t61 / 0.2e1 + t25 * t62 / 0.2e1 + t39 * t49 / 0.2e1 + t38 * t50 / 0.2e1 + t32 * t20 + t33 * t21 + (t103 * t2 - t104 * t3 - t15 * t38 + t16 * t39) * mrSges(7,3); 0.2e1 * (-t101 * t90 + t102 * t91 + t153 * t42 - t154 * t43) * mrSges(6,3) + 0.2e1 * (t103 * t5 - t104 * t6 - t32 * t38 + t33 * t39) * mrSges(7,3) + 0.2e1 * (-t140 * t197 - t141 * t248 + t162 * t191 - t163 * t190) * mrSges(5,3) - t248 * t148 + (t277 + (0.2e1 * pkin(3) * t158 - t213) * t236) * qJD(3) + t240 * t204 + t236 * t205 + 0.2e1 * t223 * t147 + t197 * t149 - 0.2e1 * pkin(2) * t203 - t190 * t159 - t191 * t160 + 0.2e1 * t166 * t106 + 0.2e1 * t172 * t57 + t153 * t58 + t154 * t59 + 0.2e1 * t118 * t17 + t103 * t18 + t104 * t19 + t102 * t107 + t101 * t108 + 0.2e1 * t74 * t60 + t39 * t61 + t38 * t62 + (t118 * t74 + t32 * t6 + t33 * t5) * t291 + (t166 * t172 + t42 * t91 + t43 * t90) * t292 + (t140 * t162 + t141 * t163 + t223 * t229) * t293; m(7) * (t136 * t3 + t137 * t2 + t15 * t87 + t16 * t86) + m(6) * (t13 * t187 + t14 * t186 + t176 * t54 - t177 * t53) + (t232 * t121 + t233 * t122 + m(5) * (t232 * t56 + t233 * t55)) * pkin(3) - Ifges(4,5) * t259 - t245 * Ifges(4,6) + t186 * t63 + t187 * t64 + t176 * t119 - t177 * t120 + t136 * t20 + t137 * t21 - t123 * mrSges(4,2) + t124 * mrSges(4,1) + t86 * t72 + t87 * t73 + t55 * mrSges(5,1) - t56 * mrSges(5,2) - t295 + t242; t228 + t243 + (pkin(8) * t253 - t281) * qJD(3) + (m(5) * (t140 * t233 + t141 * t232) + (-t190 * t232 + t191 * t233) * mrSges(5,3)) * pkin(3) + m(7) * (t136 * t6 + t137 * t5 + t32 * t87 + t33 * t86) + m(6) * (t176 * t91 - t177 * t90 + t186 * t43 + t187 * t42) + (t103 * t86 - t104 * t87 - t136 * t38 + t137 * t39) * mrSges(7,3) + (-t101 * t186 + t102 * t187 + t153 * t176 + t154 * t177) * mrSges(6,3) - t184 - t185 + t140 * mrSges(5,1) - t141 * mrSges(5,2); 0.2e1 * m(7) * (t136 * t87 + t137 * t86) + 0.2e1 * m(6) * (t176 * t187 - t177 * t186) + 0.2e1 * t244; m(5) * t169 + m(6) * t112 + m(7) * t47 + t31 + t9 + t93; m(5) * t229 + m(6) * t166 + m(7) * t74 + t147 + t17 + t57; 0; 0; (-t73 * t265 + t238 * t20 + t72 * t264 + t234 * t21 + m(7) * (-t15 * t265 + t16 * t264 + t2 * t234 + t238 * t3)) * pkin(5) + t242; (m(7) * (t234 * t5 + t238 * t6 + (-t234 * t32 + t238 * t33) * qJD(6)) + (t234 * t39 - t238 * t38 + (t103 * t238 + t104 * t234) * qJD(6)) * mrSges(7,3)) * pkin(5) + t243; (-mrSges(7,1) * t265 + m(7) * (-t136 * t265 + t137 * t264 + t234 * t86 + t238 * t87) - mrSges(7,2) * t264) * pkin(5) + t244; 0; 0.2e1 * t195; t247; t254; t257; 0; t195; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
