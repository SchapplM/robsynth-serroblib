% Calculate time derivative of joint inertia matrix for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:45
% EndTime: 2019-03-09 22:18:59
% DurationCPUTime: 5.91s
% Computational Cost: add. (14530->586), mult. (34924->870), div. (0->0), fcn. (32441->10), ass. (0->244)
t305 = -Ifges(5,3) - Ifges(6,3);
t241 = sin(qJ(3));
t242 = sin(qJ(2));
t245 = cos(qJ(3));
t270 = qJD(3) * t245;
t246 = cos(qJ(2));
t273 = qJD(2) * t246;
t250 = t241 * t273 + t242 * t270;
t237 = sin(pkin(11));
t238 = cos(pkin(11));
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t201 = t240 * t245 + t241 * t244;
t302 = qJD(3) + qJD(4);
t162 = t302 * t201;
t253 = t240 * t241 - t244 * t245;
t125 = -t162 * t242 - t253 * t273;
t186 = t253 * t242;
t274 = qJD(2) * t242;
t214 = -pkin(2) * t246 - t242 * pkin(8) - pkin(1);
t276 = t245 * t246;
t225 = pkin(7) * t276;
t212 = (pkin(2) * t242 - pkin(8) * t246) * qJD(2);
t291 = pkin(7) * t241;
t275 = t245 * t212 + t274 * t291;
t111 = (pkin(3) * t242 - pkin(9) * t276) * qJD(2) + (-t225 + (pkin(9) * t242 - t214) * t241) * qJD(3) + t275;
t199 = t245 * t214;
t278 = t242 * t245;
t156 = -pkin(9) * t278 + t199 + (-pkin(3) - t291) * t246;
t176 = t241 * t214 + t225;
t279 = t241 * t242;
t166 = -pkin(9) * t279 + t176;
t117 = t240 * t156 + t244 * t166;
t272 = qJD(3) * t241;
t133 = t241 * t212 + t214 * t270 + (-t245 * t274 - t246 * t272) * pkin(7);
t122 = -pkin(9) * t250 + t133;
t52 = -qJD(4) * t117 + t244 * t111 - t122 * t240;
t29 = pkin(4) * t274 - qJ(5) * t125 + qJD(5) * t186 + t52;
t126 = t186 * t302 - t201 * t273;
t185 = t201 * t242;
t268 = qJD(4) * t244;
t269 = qJD(4) * t240;
t51 = t240 * t111 + t244 * t122 + t156 * t268 - t166 * t269;
t31 = qJ(5) * t126 - qJD(5) * t185 + t51;
t14 = t237 * t29 + t238 * t31;
t71 = -t125 * t237 + t126 * t238;
t11 = pkin(10) * t71 + t14;
t239 = sin(qJ(6));
t243 = cos(qJ(6));
t138 = -t185 * t237 - t186 * t238;
t116 = t244 * t156 - t240 * t166;
t91 = -pkin(4) * t246 + t186 * qJ(5) + t116;
t98 = -qJ(5) * t185 + t117;
t56 = -t237 * t98 + t238 * t91;
t39 = -pkin(5) * t246 - t138 * pkin(10) + t56;
t137 = -t185 * t238 + t186 * t237;
t57 = t237 * t91 + t238 * t98;
t40 = pkin(10) * t137 + t57;
t15 = -t239 * t40 + t243 * t39;
t13 = -t237 * t31 + t238 * t29;
t72 = t125 * t238 + t126 * t237;
t9 = pkin(5) * t274 - pkin(10) * t72 + t13;
t2 = qJD(6) * t15 + t11 * t243 + t239 * t9;
t16 = t239 * t39 + t243 * t40;
t3 = -qJD(6) * t16 - t11 * t239 + t243 * t9;
t304 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t295 = -pkin(9) - pkin(8);
t218 = t295 * t241;
t219 = t295 * t245;
t168 = t240 * t218 - t244 * t219;
t262 = t245 * t273;
t303 = -Ifges(4,5) * t262 - Ifges(4,3) * t274;
t301 = 2 * m(4);
t300 = 2 * m(5);
t299 = 2 * m(6);
t298 = 2 * m(7);
t297 = -0.2e1 * pkin(1);
t296 = 0.2e1 * pkin(7);
t293 = -t241 / 0.2e1;
t292 = pkin(4) * t237;
t227 = pkin(3) * t244 + pkin(4);
t281 = t237 * t240;
t187 = -pkin(3) * t281 + t238 * t227;
t182 = pkin(5) + t187;
t280 = t238 * t240;
t189 = pkin(3) * t280 + t227 * t237;
t142 = t182 * t243 - t189 * t239;
t284 = pkin(3) * qJD(4);
t183 = (-t237 * t244 - t280) * t284;
t184 = (t238 * t244 - t281) * t284;
t96 = qJD(6) * t142 + t183 * t239 + t184 * t243;
t290 = t96 * mrSges(7,2);
t266 = qJD(3) * t295;
t210 = t241 * t266;
t211 = t245 * t266;
t128 = t244 * t210 + t240 * t211 + t218 * t268 + t219 * t269;
t80 = -qJ(5) * t162 - qJD(5) * t253 + t128;
t129 = -qJD(4) * t168 - t210 * t240 + t244 * t211;
t161 = t302 * t253;
t81 = qJ(5) * t161 - qJD(5) * t201 + t129;
t48 = t237 * t81 + t238 * t80;
t288 = Ifges(4,4) * t241;
t287 = Ifges(4,4) * t245;
t286 = Ifges(4,6) * t241;
t285 = Ifges(4,6) * t246;
t283 = t183 * mrSges(6,1);
t282 = t184 * mrSges(6,2);
t216 = Ifges(4,1) * t241 + t287;
t277 = t245 * t216;
t167 = t244 * t218 + t219 * t240;
t146 = -qJ(5) * t201 + t167;
t147 = -qJ(5) * t253 + t168;
t95 = t237 * t146 + t238 * t147;
t213 = pkin(3) * t279 + t242 * pkin(7);
t271 = qJD(3) * t242;
t86 = t137 * t243 - t138 * t239;
t25 = qJD(6) * t86 + t239 * t71 + t243 * t72;
t87 = t137 * t239 + t138 * t243;
t26 = -qJD(6) * t87 - t239 * t72 + t243 * t71;
t267 = -Ifges(7,5) * t25 - Ifges(7,6) * t26 - Ifges(7,3) * t274;
t235 = pkin(7) * t273;
t234 = pkin(3) * t272;
t173 = pkin(3) * t250 + t235;
t228 = -pkin(3) * t245 - pkin(2);
t265 = t241 * t271;
t36 = -t71 * mrSges(6,1) + t72 * mrSges(6,2);
t10 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t154 = -t201 * t237 - t238 * t253;
t155 = t201 * t238 - t237 * t253;
t101 = t154 * t243 - t155 * t239;
t114 = t161 * t237 - t162 * t238;
t115 = -t161 * t238 - t162 * t237;
t45 = qJD(6) * t101 + t114 * t239 + t115 * t243;
t102 = t154 * t239 + t155 * t243;
t46 = -qJD(6) * t102 + t114 * t243 - t115 * t239;
t17 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t261 = (2 * Ifges(3,4)) + t286;
t149 = pkin(4) * t162 + t234;
t47 = -t237 * t80 + t238 * t81;
t61 = -t114 * mrSges(6,1) + t115 * mrSges(6,2);
t226 = pkin(4) * t238 + pkin(5);
t190 = t226 * t239 + t243 * t292;
t178 = t190 * qJD(6);
t174 = t178 * mrSges(7,1);
t188 = t226 * t243 - t239 * t292;
t177 = t188 * qJD(6);
t260 = -t177 * mrSges(7,2) - t174;
t94 = t238 * t146 - t147 * t237;
t159 = pkin(4) * t185 + t213;
t42 = Ifges(7,6) * t46;
t43 = Ifges(7,5) * t45;
t32 = -pkin(10) * t115 + t47;
t33 = pkin(10) * t114 + t48;
t74 = -pkin(10) * t155 + t94;
t75 = pkin(10) * t154 + t95;
t37 = -t239 * t75 + t243 * t74;
t5 = qJD(6) * t37 + t239 * t32 + t243 * t33;
t38 = t239 * t74 + t243 * t75;
t6 = -qJD(6) * t38 - t239 * t33 + t243 * t32;
t259 = t6 * mrSges(7,1) - t5 * mrSges(7,2) + t42 + t43;
t258 = -mrSges(4,1) * t245 + mrSges(4,2) * t241;
t257 = mrSges(4,1) * t241 + mrSges(4,2) * t245;
t256 = Ifges(4,1) * t245 - t288;
t255 = -Ifges(4,2) * t241 + t287;
t215 = Ifges(4,2) * t245 + t288;
t254 = Ifges(4,5) * t241 + Ifges(4,6) * t245;
t143 = t182 * t239 + t189 * t243;
t179 = pkin(4) * t253 + t228;
t99 = -pkin(4) * t126 + t173;
t252 = (-mrSges(5,1) * t240 - mrSges(5,2) * t244) * t284;
t251 = t262 - t265;
t249 = -Ifges(5,5) * t125 - Ifges(6,5) * t72 - Ifges(5,6) * t126 - Ifges(6,6) * t71 + t274 * t305 + t267;
t109 = Ifges(6,6) * t114;
t110 = Ifges(6,5) * t115;
t157 = Ifges(5,6) * t162;
t158 = Ifges(5,5) * t161;
t248 = t129 * mrSges(5,1) + t47 * mrSges(6,1) - t128 * mrSges(5,2) - t48 * mrSges(6,2) + t109 + t110 - t157 - t158 + t259;
t247 = t52 * mrSges(5,1) + t13 * mrSges(6,1) - t51 * mrSges(5,2) - t14 * mrSges(6,2) - t249 + t304;
t233 = Ifges(4,5) * t270;
t209 = -mrSges(4,1) * t246 - mrSges(4,3) * t278;
t208 = mrSges(4,2) * t246 - mrSges(4,3) * t279;
t207 = t256 * qJD(3);
t206 = t255 * qJD(3);
t205 = t257 * qJD(3);
t181 = -Ifges(4,5) * t246 + t242 * t256;
t180 = t242 * t255 - t285;
t175 = -t246 * t291 + t199;
t172 = -mrSges(4,2) * t274 - mrSges(4,3) * t250;
t171 = mrSges(4,1) * t274 - mrSges(4,3) * t251;
t170 = -mrSges(5,1) * t246 + t186 * mrSges(5,3);
t169 = mrSges(5,2) * t246 - t185 * mrSges(5,3);
t165 = Ifges(5,1) * t201 - Ifges(5,4) * t253;
t164 = Ifges(5,4) * t201 - Ifges(5,2) * t253;
t163 = mrSges(5,1) * t253 + mrSges(5,2) * t201;
t153 = mrSges(4,1) * t250 + mrSges(4,2) * t251;
t148 = mrSges(5,1) * t185 - mrSges(5,2) * t186;
t145 = -t216 * t271 + (Ifges(4,5) * t242 + t246 * t256) * qJD(2);
t144 = -t215 * t271 + (Ifges(4,6) * t242 + t246 * t255) * qJD(2);
t136 = -Ifges(5,1) * t186 - Ifges(5,4) * t185 - Ifges(5,5) * t246;
t135 = -Ifges(5,4) * t186 - Ifges(5,2) * t185 - Ifges(5,6) * t246;
t134 = -qJD(3) * t176 + t275;
t132 = -mrSges(6,1) * t246 - t138 * mrSges(6,3);
t131 = mrSges(6,2) * t246 + t137 * mrSges(6,3);
t130 = -pkin(5) * t154 + t179;
t121 = -Ifges(5,1) * t161 - Ifges(5,4) * t162;
t120 = -Ifges(5,4) * t161 - Ifges(5,2) * t162;
t119 = mrSges(5,1) * t162 - mrSges(5,2) * t161;
t113 = -mrSges(5,2) * t274 + mrSges(5,3) * t126;
t112 = mrSges(5,1) * t274 - mrSges(5,3) * t125;
t107 = Ifges(6,1) * t155 + Ifges(6,4) * t154;
t106 = Ifges(6,4) * t155 + Ifges(6,2) * t154;
t105 = -mrSges(6,1) * t154 + mrSges(6,2) * t155;
t100 = -pkin(5) * t137 + t159;
t97 = -qJD(6) * t143 + t183 * t243 - t184 * t239;
t93 = t97 * mrSges(7,1);
t90 = -mrSges(6,1) * t137 + mrSges(6,2) * t138;
t85 = Ifges(6,1) * t138 + Ifges(6,4) * t137 - Ifges(6,5) * t246;
t84 = Ifges(6,4) * t138 + Ifges(6,2) * t137 - Ifges(6,6) * t246;
t83 = -mrSges(7,1) * t246 - t87 * mrSges(7,3);
t82 = mrSges(7,2) * t246 + t86 * mrSges(7,3);
t79 = -pkin(5) * t114 + t149;
t73 = -mrSges(5,1) * t126 + mrSges(5,2) * t125;
t67 = Ifges(5,1) * t125 + Ifges(5,4) * t126 + Ifges(5,5) * t274;
t66 = Ifges(5,4) * t125 + Ifges(5,2) * t126 + Ifges(5,6) * t274;
t65 = mrSges(6,1) * t274 - mrSges(6,3) * t72;
t64 = -mrSges(6,2) * t274 + mrSges(6,3) * t71;
t63 = Ifges(6,1) * t115 + Ifges(6,4) * t114;
t62 = Ifges(6,4) * t115 + Ifges(6,2) * t114;
t60 = Ifges(7,1) * t102 + Ifges(7,4) * t101;
t59 = Ifges(7,4) * t102 + Ifges(7,2) * t101;
t58 = -mrSges(7,1) * t101 + mrSges(7,2) * t102;
t55 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t54 = Ifges(7,1) * t87 + Ifges(7,4) * t86 - Ifges(7,5) * t246;
t53 = Ifges(7,4) * t87 + Ifges(7,2) * t86 - Ifges(7,6) * t246;
t49 = -pkin(5) * t71 + t99;
t35 = Ifges(6,1) * t72 + Ifges(6,4) * t71 + Ifges(6,5) * t274;
t34 = Ifges(6,4) * t72 + Ifges(6,2) * t71 + Ifges(6,6) * t274;
t21 = -mrSges(7,2) * t274 + mrSges(7,3) * t26;
t20 = mrSges(7,1) * t274 - mrSges(7,3) * t25;
t19 = Ifges(7,1) * t45 + Ifges(7,4) * t46;
t18 = Ifges(7,4) * t45 + Ifges(7,2) * t46;
t8 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t274;
t7 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t274;
t1 = [((mrSges(3,2) * t297 - t241 * t180 + t245 * t181 + t246 * t261) * qJD(2) + t249 + t303) * t246 + (t153 * t296 - t241 * t144 + t245 * t145 + (-t245 * t180 - t241 * t181 + t246 * t254) * qJD(3) + (mrSges(3,1) * t297 - Ifges(5,5) * t186 - Ifges(5,6) * t185 + Ifges(6,5) * t138 + Ifges(6,6) * t137 + Ifges(7,5) * t87 + Ifges(7,6) * t86 + (Ifges(4,5) * t245 - t261) * t242 + (pkin(7) ^ 2 * t301 + t257 * t296 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(7,3) + t305) * t246) * qJD(2)) * t242 + 0.2e1 * t133 * t208 + 0.2e1 * t134 * t209 + 0.2e1 * t213 * t73 - t186 * t67 - t185 * t66 + 0.2e1 * t51 * t169 + 0.2e1 * t52 * t170 + 0.2e1 * t173 * t148 + 0.2e1 * t175 * t171 + 0.2e1 * t176 * t172 + 0.2e1 * t159 * t36 + t137 * t34 + t138 * t35 + 0.2e1 * t14 * t131 + 0.2e1 * t13 * t132 + t126 * t135 + t125 * t136 + 0.2e1 * t116 * t112 + 0.2e1 * t117 * t113 + 0.2e1 * t99 * t90 + 0.2e1 * t100 * t10 + t87 * t8 + 0.2e1 * t2 * t82 + 0.2e1 * t3 * t83 + t71 * t84 + t72 * t85 + t86 * t7 + 0.2e1 * t57 * t64 + 0.2e1 * t56 * t65 + t26 * t53 + t25 * t54 + 0.2e1 * t49 * t55 + 0.2e1 * t15 * t20 + 0.2e1 * t16 * t21 + (t100 * t49 + t15 * t3 + t16 * t2) * t298 + (t13 * t56 + t14 * t57 + t159 * t99) * t299 + (t116 * t52 + t117 * t51 + t173 * t213) * t300 + (t176 * t133 + t175 * t134) * t301; (t116 * t161 - t117 * t162 - t201 * t52 - t253 * t51) * mrSges(5,3) + (t245 * t207 / 0.2e1 + t206 * t293 - Ifges(3,6) * qJD(2) + (t216 * t293 - t245 * t215 / 0.2e1) * qJD(3) + (mrSges(3,2) * qJD(2) + t205) * pkin(7) + (Ifges(5,5) * t201 + Ifges(6,5) * t155 + Ifges(7,5) * t102 - Ifges(5,6) * t253 + Ifges(6,6) * t154 + Ifges(7,6) * t101 + t254) * qJD(2) / 0.2e1) * t242 - t253 * t66 / 0.2e1 + (t144 / 0.2e1 + pkin(8) * t172 + t133 * mrSges(4,3)) * t245 + (t145 / 0.2e1 - pkin(8) * t171 - t134 * mrSges(4,3)) * t241 + ((t181 / 0.2e1 - t175 * mrSges(4,3) - pkin(8) * t209) * t245 + (t285 / 0.2e1 - t180 / 0.2e1 + pkin(3) * t148 - t176 * mrSges(4,3) - pkin(8) * t208) * t241) * qJD(3) + m(5) * (t116 * t129 + t117 * t128 + t167 * t52 + t168 * t51 + t173 * t228 + t213 * t234) + (-pkin(2) * t235 + (t133 * t245 - t134 * t241 + (-t175 * t245 - t176 * t241) * qJD(3)) * pkin(8)) * m(4) + (t101 * t2 - t102 * t3 - t15 * t45 + t16 * t46) * mrSges(7,3) + (t114 * t57 - t115 * t56 - t13 * t155 + t14 * t154) * mrSges(6,3) + t228 * t73 + t213 * t119 + t201 * t67 / 0.2e1 - t186 * t121 / 0.2e1 + t179 * t36 - t185 * t120 / 0.2e1 + t167 * t112 + t168 * t113 + t128 * t169 + t129 * t170 + t173 * t163 + t159 * t61 - t161 * t136 / 0.2e1 - t162 * t135 / 0.2e1 + t126 * t164 / 0.2e1 + t125 * t165 / 0.2e1 - pkin(2) * t153 + t154 * t34 / 0.2e1 + t155 * t35 / 0.2e1 + t137 * t62 / 0.2e1 + t138 * t63 / 0.2e1 + t149 * t90 + t130 * t10 + t48 * t131 + t47 * t132 + t115 * t85 / 0.2e1 + t114 * t84 / 0.2e1 + t101 * t7 / 0.2e1 + t102 * t8 / 0.2e1 + t99 * t105 + t71 * t106 / 0.2e1 + t72 * t107 / 0.2e1 + t94 * t65 + t95 * t64 + t100 * t17 + t87 * t19 / 0.2e1 + t5 * t82 + t6 * t83 + t86 * t18 / 0.2e1 + t79 * t55 + t49 * t58 + t26 * t59 / 0.2e1 + t25 * t60 / 0.2e1 + t46 * t53 / 0.2e1 + t45 * t54 / 0.2e1 + t37 * t20 + t38 * t21 + m(6) * (t13 * t94 + t14 * t95 + t149 * t159 + t179 * t99 + t47 * t56 + t48 * t57) + m(7) * (t100 * t79 + t130 * t49 + t15 * t6 + t16 * t5 + t2 * t38 + t3 * t37) + (-t233 / 0.2e1 + t158 / 0.2e1 + t157 / 0.2e1 - t110 / 0.2e1 - t109 / 0.2e1 - t43 / 0.2e1 - t42 / 0.2e1 + (Ifges(3,5) + t215 * t293 + t277 / 0.2e1 + (-mrSges(3,1) + t258) * pkin(7)) * qJD(2)) * t246; 0.2e1 * (-t128 * t253 - t129 * t201 + t161 * t167 - t162 * t168) * mrSges(5,3) - t253 * t120 + (t277 + (0.2e1 * pkin(3) * t163 - t215) * t241) * qJD(3) + t245 * t206 + t241 * t207 + 0.2e1 * t228 * t119 + t201 * t121 - 0.2e1 * pkin(2) * t205 + 0.2e1 * t179 * t61 - t162 * t164 - t161 * t165 + t154 * t62 + t155 * t63 + 0.2e1 * t149 * t105 + 0.2e1 * t130 * t17 + t115 * t107 + t114 * t106 + t101 * t18 + t102 * t19 + 0.2e1 * t79 * t58 + t46 * t59 + t45 * t60 + (t130 * t79 + t37 * t6 + t38 * t5) * t298 + (t149 * t179 + t47 * t94 + t48 * t95) * t299 + (t128 * t168 + t129 * t167 + t228 * t234) * t300 + 0.2e1 * (t101 * t5 - t102 * t6 - t37 * t45 + t38 * t46) * mrSges(7,3) + 0.2e1 * (t114 * t95 - t115 * t94 + t154 * t48 - t155 * t47) * mrSges(6,3); m(6) * (t13 * t187 + t14 * t189 + t183 * t56 + t184 * t57) + m(7) * (t142 * t3 + t143 * t2 + t15 * t97 + t16 * t96) - t250 * Ifges(4,6) + t247 + t187 * t65 + t189 * t64 + t183 * t132 + t184 * t131 + t142 * t20 + t143 * t21 - t133 * mrSges(4,2) + t134 * mrSges(4,1) + t96 * t82 + t97 * t83 + (m(5) * (-t116 * t269 + t117 * t268 + t240 * t51 + t244 * t52) + t244 * t112 + t240 * t113 - t170 * t269 + t169 * t268) * pkin(3) - Ifges(4,5) * t265 - t303; (pkin(8) * t258 - t286) * qJD(3) + m(6) * (t183 * t94 + t184 * t95 + t187 * t47 + t189 * t48) + m(7) * (t142 * t6 + t143 * t5 + t37 * t97 + t38 * t96) + t248 + (t101 * t96 - t102 * t97 - t142 * t45 + t143 * t46) * mrSges(7,3) + (t114 * t189 - t115 * t187 + t154 * t184 - t155 * t183) * mrSges(6,3) + t233 + (m(5) * (t128 * t240 + t129 * t244 + (-t167 * t240 + t168 * t244) * qJD(4)) + (t244 * t161 - t240 * t162 + (t201 * t240 - t244 * t253) * qJD(4)) * mrSges(5,3)) * pkin(3); 0.2e1 * t283 - 0.2e1 * t282 - 0.2e1 * t290 + 0.2e1 * t93 + 0.2e1 * t252 + (t142 * t97 + t143 * t96) * t298 + (t183 * t187 + t184 * t189) * t299; t247 + t188 * t20 + t190 * t21 - t178 * t83 + t177 * t82 + m(7) * (-t15 * t178 + t16 * t177 + t188 * t3 + t190 * t2) + (t237 * t64 + t238 * t65 + m(6) * (t13 * t238 + t14 * t237)) * pkin(4); t248 + (t101 * t177 + t102 * t178 - t188 * t45 + t190 * t46) * mrSges(7,3) + (m(6) * (t237 * t48 + t238 * t47) + (t114 * t237 - t115 * t238) * mrSges(6,3)) * pkin(4) + m(7) * (t177 * t38 - t178 * t37 + t188 * t6 + t190 * t5); t283 - t282 - t174 + t93 + t252 + (-t177 - t96) * mrSges(7,2) + m(7) * (-t142 * t178 + t143 * t177 + t188 * t97 + t190 * t96) + m(6) * (t183 * t238 + t184 * t237) * pkin(4); 0.2e1 * m(7) * (t177 * t190 - t178 * t188) + 0.2e1 * t260; m(6) * t99 + m(7) * t49 + t10 + t36; m(6) * t149 + m(7) * t79 + t17 + t61; 0; 0; 0; -t267 + t304; t259; t93 - t290; t260; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
