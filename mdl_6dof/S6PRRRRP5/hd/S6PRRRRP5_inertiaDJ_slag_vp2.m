% Calculate time derivative of joint inertia matrix for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:19
% EndTime: 2019-03-09 00:20:36
% DurationCPUTime: 7.15s
% Computational Cost: add. (5772->654), mult. (16801->908), div. (0->0), fcn. (15969->12), ass. (0->257)
t214 = cos(qJ(5));
t325 = Ifges(6,6) + Ifges(7,6);
t328 = t325 * t214;
t210 = sin(qJ(5));
t327 = (Ifges(6,5) + Ifges(7,5)) * t210;
t215 = cos(qJ(4));
t264 = qJD(4) * t215;
t241 = t214 * t264;
t211 = sin(qJ(4));
t262 = qJD(5) * t211;
t244 = t210 * t262;
t220 = t241 - t244;
t242 = t210 * t264;
t261 = qJD(5) * t214;
t219 = t211 * t261 + t242;
t206 = sin(pkin(7));
t212 = sin(qJ(3));
t266 = qJD(3) * t212;
t246 = t206 * t266;
t208 = cos(pkin(7));
t283 = t206 * t212;
t137 = -t215 * t208 + t211 * t283;
t216 = cos(qJ(3));
t267 = qJD(3) * t206;
t245 = t216 * t267;
t97 = -qJD(4) * t137 + t215 * t245;
t138 = t208 * t211 + t215 * t283;
t282 = t206 * t216;
t99 = -t138 * t210 - t214 * t282;
t50 = qJD(5) * t99 + t210 * t246 + t214 * t97;
t222 = -t138 * t214 + t210 * t282;
t51 = qJD(5) * t222 - t210 * t97 + t214 * t246;
t98 = qJD(4) * t138 + t211 * t245;
t8 = Ifges(7,5) * t50 + Ifges(7,6) * t51 + Ifges(7,3) * t98;
t9 = Ifges(6,5) * t50 + Ifges(6,6) * t51 + Ifges(6,3) * t98;
t324 = t8 + t9;
t320 = (-m(6) * pkin(11) - mrSges(6,3));
t192 = pkin(9) * t283;
t304 = pkin(2) * t216;
t141 = t208 * t304 - t192;
t169 = -pkin(4) * t215 - pkin(11) * t211 - pkin(3);
t276 = t214 * t215;
t196 = pkin(10) * t276;
t119 = t210 * t169 + t196;
t217 = cos(qJ(2));
t275 = t216 * t217;
t213 = sin(qJ(2));
t279 = t212 * t213;
t319 = t208 * t275 - t279;
t172 = -mrSges(6,1) * t214 + mrSges(6,2) * t210;
t318 = -m(6) * pkin(4) - mrSges(5,1) + t172;
t288 = t211 * mrSges(5,2);
t317 = -m(5) * pkin(3) - t215 * mrSges(5,1) - mrSges(4,1) + t288;
t171 = -mrSges(7,1) * t214 + mrSges(7,2) * t210;
t198 = -pkin(5) * t214 - pkin(4);
t316 = m(7) * t198 + t171;
t315 = 0.2e1 * m(6);
t314 = 0.2e1 * m(7);
t313 = 0.2e1 * pkin(10);
t312 = -2 * mrSges(4,3);
t311 = -2 * mrSges(7,3);
t310 = m(6) / 0.2e1;
t306 = m(7) * pkin(5);
t142 = t208 * t212 * pkin(2) + pkin(9) * t282;
t122 = pkin(10) * t208 + t142;
t123 = (-pkin(3) * t216 - pkin(10) * t212 - pkin(2)) * t206;
t131 = (pkin(3) * t212 - pkin(10) * t216) * t267;
t132 = t141 * qJD(3);
t265 = qJD(4) * t211;
t33 = -t122 * t264 - t123 * t265 + t131 * t215 - t211 * t132;
t29 = -pkin(4) * t246 - t33;
t305 = m(6) * t29;
t303 = pkin(10) * t210;
t207 = sin(pkin(6));
t268 = qJD(2) * t207;
t247 = t213 * t268;
t231 = t206 * t247;
t209 = cos(pkin(6));
t60 = t209 * t245 + (t319 * qJD(3) + (-t208 * t279 + t275) * qJD(2)) * t207;
t136 = -t206 * t207 * t217 + t209 * t208;
t277 = t213 * t216;
t278 = t212 * t217;
t221 = t208 * t278 + t277;
t92 = t207 * t221 + t209 * t283;
t65 = t136 * t211 + t215 * t92;
t19 = qJD(4) * t65 + t211 * t60 - t215 * t231;
t224 = t136 * t215 - t211 * t92;
t302 = t19 * t224;
t59 = t209 * t246 + (t221 * qJD(3) + (t208 * t277 + t278) * qJD(2)) * t207;
t91 = -t207 * t319 - t209 * t282;
t301 = t59 * t91;
t299 = -qJ(6) - pkin(11);
t17 = -mrSges(6,1) * t51 + mrSges(6,2) * t50;
t79 = mrSges(5,1) * t246 - mrSges(5,3) * t97;
t298 = t17 - t79;
t121 = t192 + (-pkin(3) - t304) * t208;
t66 = pkin(4) * t137 - pkin(11) * t138 + t121;
t76 = t215 * t122 + t211 * t123;
t68 = -pkin(11) * t282 + t76;
t22 = t210 * t66 + t214 * t68;
t297 = Ifges(5,4) * t211;
t296 = Ifges(5,4) * t215;
t295 = Ifges(6,4) * t210;
t294 = Ifges(6,4) * t214;
t293 = Ifges(7,4) * t210;
t292 = Ifges(7,4) * t214;
t133 = t142 * qJD(3);
t291 = t133 * t91;
t290 = t19 * t211;
t20 = qJD(4) * t224 + t211 * t231 + t215 * t60;
t289 = t20 * t215;
t286 = -mrSges(4,1) * t208 + mrSges(5,1) * t137 + mrSges(5,2) * t138 + mrSges(4,3) * t283;
t103 = -mrSges(5,1) * t282 - mrSges(5,3) * t138;
t57 = -mrSges(6,1) * t99 - mrSges(6,2) * t222;
t284 = t57 - t103;
t281 = t210 * t211;
t280 = t211 * t214;
t225 = -Ifges(7,2) * t210 + t292;
t126 = -Ifges(7,6) * t215 + t211 * t225;
t226 = -Ifges(6,2) * t210 + t294;
t127 = -Ifges(6,6) * t215 + t211 * t226;
t274 = -t126 - t127;
t227 = Ifges(7,1) * t214 - t293;
t128 = -Ifges(7,5) * t215 + t211 * t227;
t228 = Ifges(6,1) * t214 - t295;
t129 = -Ifges(6,5) * t215 + t211 * t228;
t273 = t128 + t129;
t166 = (pkin(4) * t211 - pkin(11) * t215) * qJD(4);
t272 = t210 * t166 + t169 * t261;
t271 = t214 * t166 + t265 * t303;
t270 = Ifges(7,5) * t241 + Ifges(7,3) * t265;
t269 = Ifges(6,5) * t241 + Ifges(6,3) * t265;
t263 = qJD(5) * t210;
t151 = mrSges(7,1) * t263 + mrSges(7,2) * t261;
t260 = qJD(6) * t214;
t258 = Ifges(5,5) * t97 - Ifges(5,6) * t98 + Ifges(5,3) * t246;
t257 = pkin(5) * t263;
t256 = Ifges(5,6) * t282;
t10 = Ifges(7,4) * t50 + Ifges(7,2) * t51 + Ifges(7,6) * t98;
t11 = Ifges(6,4) * t50 + Ifges(6,2) * t51 + Ifges(6,6) * t98;
t255 = -t10 / 0.2e1 - t11 / 0.2e1;
t12 = Ifges(7,1) * t50 + Ifges(7,4) * t51 + Ifges(7,5) * t98;
t13 = Ifges(6,1) * t50 + Ifges(6,4) * t51 + Ifges(6,5) * t98;
t254 = t12 / 0.2e1 + t13 / 0.2e1;
t39 = -Ifges(7,4) * t222 + Ifges(7,2) * t99 + Ifges(7,6) * t137;
t40 = -Ifges(6,4) * t222 + Ifges(6,2) * t99 + Ifges(6,6) * t137;
t253 = -t39 / 0.2e1 - t40 / 0.2e1;
t41 = -Ifges(7,1) * t222 + Ifges(7,4) * t99 + Ifges(7,5) * t137;
t42 = -Ifges(6,1) * t222 + Ifges(6,4) * t99 + Ifges(6,5) * t137;
t252 = t41 / 0.2e1 + t42 / 0.2e1;
t177 = Ifges(7,2) * t214 + t293;
t83 = -t177 * t262 + (Ifges(7,6) * t211 + t215 * t225) * qJD(4);
t178 = Ifges(6,2) * t214 + t295;
t84 = -t178 * t262 + (Ifges(6,6) * t211 + t215 * t226) * qJD(4);
t251 = t83 / 0.2e1 + t84 / 0.2e1;
t180 = Ifges(7,1) * t210 + t292;
t85 = -t180 * t262 + (Ifges(7,5) * t211 + t215 * t227) * qJD(4);
t181 = Ifges(6,1) * t210 + t294;
t86 = -t181 * t262 + (Ifges(6,5) * t211 + t215 * t228) * qJD(4);
t250 = t85 / 0.2e1 + t86 / 0.2e1;
t249 = mrSges(7,1) + t306;
t240 = t126 / 0.2e1 + t127 / 0.2e1;
t239 = t128 / 0.2e1 + t129 / 0.2e1;
t203 = Ifges(7,5) * t261;
t204 = Ifges(6,5) * t261;
t238 = t203 / 0.2e1 + t204 / 0.2e1 - t325 * t263 / 0.2e1;
t156 = t225 * qJD(5);
t157 = t226 * qJD(5);
t237 = t156 / 0.2e1 + t157 / 0.2e1;
t159 = t227 * qJD(5);
t160 = t228 * qJD(5);
t236 = t160 / 0.2e1 + t159 / 0.2e1;
t235 = t327 / 0.2e1 + t328 / 0.2e1;
t234 = t177 / 0.2e1 + t178 / 0.2e1;
t233 = t180 / 0.2e1 + t181 / 0.2e1;
t16 = -t51 * mrSges(7,1) + t50 * mrSges(7,2);
t21 = -t210 * t68 + t214 * t66;
t232 = qJD(5) * t299;
t75 = -t211 * t122 + t123 * t215;
t230 = t246 / 0.2e1;
t67 = pkin(4) * t282 - t75;
t229 = mrSges(6,1) * t210 + mrSges(6,2) * t214;
t31 = t210 * t91 + t214 * t65;
t30 = -t210 * t65 + t214 * t91;
t223 = -t224 * t264 + t290;
t32 = -t122 * t265 + t123 * t264 + t211 * t131 + t215 * t132;
t28 = pkin(11) * t246 + t32;
t43 = pkin(4) * t98 - pkin(11) * t97 + t133;
t4 = t210 * t43 + t214 * t28 + t66 * t261 - t263 * t68;
t88 = mrSges(7,1) * t219 + mrSges(7,2) * t220;
t5 = -qJD(5) * t22 - t210 * t28 + t214 * t43;
t205 = Ifges(5,5) * t264;
t185 = Ifges(4,5) * t245;
t182 = Ifges(5,1) * t211 + t296;
t179 = Ifges(5,2) * t215 + t297;
t174 = t299 * t214;
t170 = t299 * t210;
t168 = (pkin(5) * t210 + pkin(10)) * t211;
t165 = -mrSges(6,1) * t215 - mrSges(6,3) * t280;
t164 = -mrSges(7,1) * t215 - mrSges(7,3) * t280;
t163 = mrSges(6,2) * t215 - mrSges(6,3) * t281;
t162 = mrSges(7,2) * t215 - mrSges(7,3) * t281;
t161 = (Ifges(5,1) * t215 - t297) * qJD(4);
t158 = (-Ifges(5,2) * t211 + t296) * qJD(4);
t153 = (mrSges(5,1) * t211 + mrSges(5,2) * t215) * qJD(4);
t152 = t229 * qJD(5);
t150 = -mrSges(4,2) * t208 + mrSges(4,3) * t282;
t148 = t214 * t169;
t144 = t229 * t211;
t143 = (mrSges(7,1) * t210 + mrSges(7,2) * t214) * t211;
t135 = -qJD(6) * t210 + t214 * t232;
t134 = t210 * t232 + t260;
t130 = (mrSges(4,1) * t212 + mrSges(4,2) * t216) * t267;
t125 = -Ifges(6,3) * t215 + (Ifges(6,5) * t214 - Ifges(6,6) * t210) * t211;
t124 = -Ifges(7,3) * t215 + (Ifges(7,5) * t214 - Ifges(7,6) * t210) * t211;
t118 = -t215 * t303 + t148;
t114 = pkin(5) * t219 + pkin(10) * t264;
t110 = -mrSges(6,2) * t265 - mrSges(6,3) * t219;
t109 = -mrSges(7,2) * t265 - mrSges(7,3) * t219;
t108 = mrSges(6,1) * t265 - mrSges(6,3) * t220;
t107 = mrSges(7,1) * t265 - mrSges(7,3) * t220;
t102 = mrSges(5,2) * t282 - mrSges(5,3) * t137;
t101 = -qJ(6) * t281 + t119;
t90 = -qJ(6) * t280 + t148 + (-pkin(5) - t303) * t215;
t89 = mrSges(6,1) * t219 + mrSges(6,2) * t220;
t82 = -Ifges(6,5) * t244 - Ifges(6,6) * t219 + t269;
t81 = -Ifges(7,5) * t244 - Ifges(7,6) * t219 + t270;
t80 = -mrSges(5,2) * t246 - mrSges(5,3) * t98;
t78 = Ifges(5,1) * t138 - Ifges(5,4) * t137 - Ifges(5,5) * t282;
t77 = Ifges(5,4) * t138 - Ifges(5,2) * t137 - t256;
t74 = -qJD(5) * t119 + t271;
t73 = (-t214 * t265 - t215 * t263) * pkin(10) + t272;
t72 = mrSges(6,1) * t137 + mrSges(6,3) * t222;
t71 = mrSges(7,1) * t137 + mrSges(7,3) * t222;
t70 = -mrSges(6,2) * t137 + mrSges(6,3) * t99;
t69 = -mrSges(7,2) * t137 + mrSges(7,3) * t99;
t56 = -mrSges(7,1) * t99 - mrSges(7,2) * t222;
t55 = mrSges(5,1) * t98 + mrSges(5,2) * t97;
t54 = Ifges(5,1) * t97 - Ifges(5,4) * t98 + Ifges(5,5) * t246;
t53 = Ifges(5,4) * t97 - Ifges(5,2) * t98 + Ifges(5,6) * t246;
t52 = (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t280 + (-qJD(6) * t211 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t215) * t210 + t272;
t44 = -t211 * t260 + (pkin(5) * t211 - qJ(6) * t276) * qJD(4) + (-t196 + (qJ(6) * t211 - t169) * t210) * qJD(5) + t271;
t38 = -Ifges(6,5) * t222 + Ifges(6,6) * t99 + Ifges(6,3) * t137;
t37 = -Ifges(7,5) * t222 + Ifges(7,6) * t99 + Ifges(7,3) * t137;
t34 = -pkin(5) * t99 + t67;
t26 = -mrSges(6,2) * t98 + mrSges(6,3) * t51;
t25 = -mrSges(7,2) * t98 + mrSges(7,3) * t51;
t24 = mrSges(6,1) * t98 - mrSges(6,3) * t50;
t23 = mrSges(7,1) * t98 - mrSges(7,3) * t50;
t18 = qJ(6) * t99 + t22;
t15 = pkin(5) * t137 + qJ(6) * t222 + t21;
t14 = -pkin(5) * t51 + t29;
t7 = qJD(5) * t30 + t20 * t214 + t210 * t59;
t6 = -qJD(5) * t31 - t20 * t210 + t214 * t59;
t3 = qJ(6) * t51 + qJD(6) * t99 + t4;
t2 = pkin(5) * t98 - qJ(6) * t50 + qJD(6) * t222 + t5;
t1 = [0.2e1 * m(5) * (t20 * t65 + t301 - t302) + 0.2e1 * m(4) * (t136 * t231 + t60 * t92 + t301) + 0.4e1 * (m(7) / 0.2e1 + t310) * (t30 * t6 + t31 * t7 - t302); t20 * t102 + t136 * t130 + t60 * t150 + t91 * t55 + t65 * t80 + (t69 + t70) * t7 + (t71 + t72) * t6 + t286 * t59 + (t25 + t26) * t31 + (t23 + t24) * t30 + (-mrSges(3,1) * t213 - mrSges(3,2) * t217) * t268 - (t16 + t298) * t224 + (t56 + t284) * t19 + ((-mrSges(4,1) * t216 + mrSges(4,2) * t212) * t231 + (-t212 * t92 + t216 * t91) * qJD(3) * mrSges(4,3)) * t206 + m(4) * (-pkin(2) * t206 ^ 2 * t247 + t132 * t92 - t141 * t59 + t142 * t60 + t291) + m(5) * (t121 * t59 - t19 * t75 + t20 * t76 + t224 * t33 + t32 * t65 + t291) + m(7) * (-t14 * t224 + t15 * t6 + t18 * t7 + t19 * t34 + t2 * t30 + t3 * t31) + m(6) * (t19 * t67 + t21 * t6 + t22 * t7 - t224 * t29 + t30 * t5 + t31 * t4); (t14 * t34 + t15 * t2 + t18 * t3) * t314 + (t21 * t5 + t22 * t4 + t29 * t67) * t315 + 0.2e1 * t132 * t150 + t138 * t54 + 0.2e1 * t121 * t55 + 0.2e1 * t32 * t102 + 0.2e1 * t33 * t103 + t97 * t78 + 0.2e1 * t75 * t79 + 0.2e1 * t76 * t80 + 0.2e1 * t67 * t17 + 0.2e1 * t3 * t69 + 0.2e1 * t4 * t70 + 0.2e1 * t2 * t71 + 0.2e1 * t5 * t72 + 0.2e1 * t14 * t56 + 0.2e1 * t29 * t57 + 0.2e1 * t34 * t16 + 0.2e1 * t15 * t23 + 0.2e1 * t21 * t24 + 0.2e1 * t18 * t25 + 0.2e1 * t22 * t26 + 0.2e1 * m(4) * (t132 * t142 - t133 * t141) + 0.2e1 * m(5) * (t121 * t133 + t32 * t76 + t33 * t75) + (t10 + t11) * t99 + (t37 + t38 - t77) * t98 + (t39 + t40) * t51 + (t41 + t42) * t50 + (-t216 * t258 - 0.2e1 * pkin(2) * t130 + ((0.2e1 * Ifges(4,4) * t282 + Ifges(4,5) * t208 + t141 * t312) * t216 + (-0.2e1 * Ifges(4,4) * t283 + t142 * t312 + Ifges(5,5) * t138 - 0.2e1 * t208 * Ifges(4,6) - Ifges(5,6) * t137 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t282) * t212) * qJD(3)) * t206 + 0.2e1 * t286 * t133 + t208 * t185 + (-t53 + t324) * t137 - (t12 + t13) * t222; -t60 * mrSges(4,2) + t91 * t153 + (t162 + t163) * t7 - (t88 + t89) * t224 + (t164 + t165) * t6 + (t109 + t110) * t31 + (t107 + t108) * t30 + (t143 + t144) * t19 + m(7) * (t101 * t7 - t114 * t224 + t168 * t19 + t30 * t44 + t31 * t52 + t6 * t90) + m(6) * (t118 * t6 + t119 * t7 + t30 * t74 + t31 * t73) + (t223 * t310 + m(5) * (-t265 * t65 + t223 + t289) / 0.2e1) * t313 + (t290 + t289 + (-t211 * t65 - t215 * t224) * qJD(4)) * mrSges(5,3) + t317 * t59; t185 + t97 * t182 / 0.2e1 + t138 * t161 / 0.2e1 + t3 * t162 + t4 * t163 + t2 * t164 + t5 * t165 + t168 * t16 + t14 * t143 + t29 * t144 + t121 * t153 - t132 * mrSges(4,2) + t118 * t24 + t119 * t26 + t15 * t107 + t21 * t108 + t18 * t109 + t22 * t110 + t114 * t56 + t101 * t25 + t34 * t88 + t67 * t89 + t90 * t23 + t73 * t70 + t74 * t72 + t52 * t69 + t44 * t71 - pkin(3) * t55 + m(7) * (t101 * t3 + t114 * t34 + t14 * t168 + t15 * t44 + t18 * t52 + t2 * t90) + (-t179 / 0.2e1 + t124 / 0.2e1 + t125 / 0.2e1) * t98 + (t54 / 0.2e1 - t33 * mrSges(5,3) + Ifges(5,5) * t230 + t254 * t214 + t255 * t210 + (-t210 * t252 + t214 * t253) * qJD(5) + (t256 / 0.2e1 - t77 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1 - t76 * mrSges(5,3)) * qJD(4) + (-qJD(4) * t102 + m(5) * (-qJD(4) * t76 - t33) + t305 + t298) * pkin(10)) * t211 + (t53 / 0.2e1 - t8 / 0.2e1 - t9 / 0.2e1 + t32 * mrSges(5,3) + Ifges(5,6) * t230 + (m(5) * t32 + t80) * pkin(10) + (t78 / 0.2e1 - t75 * mrSges(5,3) + t252 * t214 + t253 * t210 + (-m(5) * t75 + m(6) * t67 + t284) * pkin(10)) * qJD(4)) * t215 + (-t216 * t205 / 0.2e1 - Ifges(4,6) * t266) * t206 + t251 * t99 + t239 * t50 + t240 * t51 + (-t158 / 0.2e1 + t81 / 0.2e1 + t82 / 0.2e1) * t137 + t317 * t133 + m(6) * (t118 * t5 + t119 * t4 + t21 * t74 + t22 * t73) - t250 * t222; -0.2e1 * pkin(3) * t153 + 0.2e1 * t101 * t109 + 0.2e1 * t90 * t107 + 0.2e1 * t118 * t108 + 0.2e1 * t119 * t110 + 0.2e1 * t114 * t143 + 0.2e1 * t52 * t162 + 0.2e1 * t73 * t163 + 0.2e1 * t44 * t164 + 0.2e1 * t74 * t165 + 0.2e1 * t168 * t88 + (t118 * t74 + t119 * t73) * t315 + (t101 * t52 + t114 * t168 + t44 * t90) * t314 + (t158 - t81 - t82 + (t144 * t313 + t210 * t274 + t214 * t273 + t182) * qJD(4)) * t215 + (t89 * t313 + t161 + (t85 + t86) * t214 + (-t83 - t84) * t210 + (-t210 * t273 + t214 * t274) * qJD(5) + (pkin(10) ^ 2 * t215 * t315 + t124 + t125 - t179) * qJD(4)) * t211; -t20 * mrSges(5,2) - (t151 + t152) * t224 + m(7) * (t134 * t31 + t135 * t30 + t170 * t6 - t174 * t7 - t224 * t257) + (t316 + t318) * t19 + (mrSges(7,3) - t320) * (-t6 * t210 + t7 * t214 + (-t210 * t31 - t214 * t30) * qJD(5)); t198 * t16 + t170 * t23 + t14 * t171 + t29 * t172 - t174 * t25 + t34 * t151 + t67 * t152 + t135 * t71 + t134 * t69 - t32 * mrSges(5,2) + t33 * mrSges(5,1) + (-t5 * mrSges(6,3) - t2 * mrSges(7,3) + (-m(6) * t5 - t24) * pkin(11) + (-t18 * mrSges(7,3) + pkin(5) * t56 - pkin(11) * t70 + t22 * t320 + t306 * t34 + t253) * qJD(5) + t254) * t210 + (t3 * mrSges(7,3) + t4 * mrSges(6,3) + (-t21 * mrSges(6,3) - t15 * mrSges(7,3) + t252) * qJD(5) + (m(6) * (-qJD(5) * t21 + t4) + t26 - qJD(5) * t72) * pkin(11) - t255) * t214 + t238 * t137 + t233 * t50 + t234 * t51 + t235 * t98 - t236 * t222 + t237 * t99 + m(7) * (t134 * t18 + t135 * t15 + t14 * t198 + t170 * t2 - t174 * t3) + t258 + (-t17 - t305) * pkin(4); t205 + t198 * t88 + t134 * t162 + t135 * t164 + t168 * t151 + t170 * t107 + t114 * t171 - t174 * t109 - pkin(4) * t89 + t211 * pkin(10) * t152 + m(7) * (t101 * t134 + t114 * t198 + t135 * t90 + t170 * t44 - t174 * t52) - t238 * t215 + ((-Ifges(5,6) + t235) * t211 + (t215 * t318 + t288) * pkin(10)) * qJD(4) + (-t74 * mrSges(6,3) - t44 * mrSges(7,3) - t237 * t211 - t234 * t264 + (-m(6) * t74 - t108) * pkin(11) + (-t101 * mrSges(7,3) + pkin(5) * t143 - pkin(11) * t163 + t119 * t320 + t168 * t306 - t211 * t233 - t240) * qJD(5) + t250) * t210 + (t73 * mrSges(6,3) + t52 * mrSges(7,3) + t236 * t211 + t233 * t264 + (m(6) * t73 + t110) * pkin(11) + (-t118 * mrSges(6,3) - t90 * mrSges(7,3) - t234 * t211 + (-m(6) * t118 - t165) * pkin(11) + t239) * qJD(5) + t251) * t214; -0.2e1 * pkin(4) * t152 + 0.2e1 * t198 * t151 + (-t134 * t174 + t135 * t170) * t314 + (t135 * t311 + t159 + t160 + (0.2e1 * pkin(5) * t316 - t174 * t311 - t177 - t178) * qJD(5)) * t210 + (0.2e1 * t134 * mrSges(7,3) + t156 + t157 + (t170 * t311 + t180 + t181) * qJD(5)) * t214; (-mrSges(6,2) - mrSges(7,2)) * t7 + (mrSges(6,1) + t249) * t6; mrSges(6,1) * t5 + mrSges(7,1) * t2 - mrSges(6,2) * t4 - mrSges(7,2) * t3 + (m(7) * t2 + t23) * pkin(5) + t324; mrSges(6,1) * t74 + mrSges(7,1) * t44 - mrSges(6,2) * t73 - mrSges(7,2) * t52 - t325 * t242 + (m(7) * t44 + t107) * pkin(5) + (-t327 - t328) * t262 + t269 + t270; -mrSges(7,2) * t134 + t203 + t204 + t249 * t135 + ((-mrSges(6,1) * pkin(11) - mrSges(7,3) * pkin(5)) * t214 + (mrSges(6,2) * pkin(11) - t325) * t210) * qJD(5); 0; m(7) * t19; m(7) * t14 + t16; m(7) * t114 + t88; m(7) * t257 + t151; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
