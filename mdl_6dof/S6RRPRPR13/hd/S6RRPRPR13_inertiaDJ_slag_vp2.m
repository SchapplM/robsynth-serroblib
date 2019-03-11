% Calculate time derivative of joint inertia matrix for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:13
% EndTime: 2019-03-09 11:24:25
% DurationCPUTime: 5.39s
% Computational Cost: add. (6796->620), mult. (17470->915), div. (0->0), fcn. (16131->10), ass. (0->254)
t223 = (-pkin(2) - pkin(9));
t306 = 2 * t223;
t213 = sin(pkin(11));
t215 = cos(pkin(11));
t217 = sin(qJ(6));
t220 = cos(qJ(6));
t225 = t213 * t217 - t215 * t220;
t279 = -t225 / 0.2e1;
t177 = t213 * t220 + t215 * t217;
t278 = t177 / 0.2e1;
t305 = t215 / 0.2e1;
t216 = cos(pkin(6));
t221 = cos(qJ(4));
t218 = sin(qJ(4));
t214 = sin(pkin(6));
t222 = cos(qJ(2));
t256 = t214 * t222;
t243 = t218 * t256;
t167 = t216 * t221 - t243;
t219 = sin(qJ(2));
t257 = t214 * t219;
t125 = -t167 * t213 + t215 * t257;
t126 = t167 * t215 + t213 * t257;
t166 = t216 * t218 + t221 * t256;
t62 = t125 * t220 - t126 * t217;
t63 = t125 * t217 + t126 * t220;
t304 = Ifges(6,5) * t126 + Ifges(7,5) * t63 + Ifges(6,6) * t125 + Ifges(7,6) * t62 + (Ifges(6,3) + Ifges(7,3)) * t166;
t248 = qJD(2) * t222;
t239 = t214 * t248;
t249 = qJD(2) * t214;
t240 = t219 * t249;
t201 = pkin(2) * t240;
t247 = qJD(3) * t219;
t108 = t201 + (-t247 + (pkin(9) * t219 - qJ(3) * t222) * qJD(2)) * t214;
t203 = pkin(8) * t257;
t275 = pkin(1) * t222;
t241 = -pkin(2) - t275;
t113 = pkin(3) * t257 + t203 + (-pkin(9) + t241) * t216;
t262 = qJ(3) * t219;
t132 = (t222 * t223 - pkin(1) - t262) * t214;
t206 = t216 * t219 * pkin(1);
t290 = pkin(3) + pkin(8);
t133 = (t256 * t290 + t206) * qJD(2);
t245 = qJD(4) * t221;
t246 = qJD(4) * t218;
t35 = -t218 * t108 - t113 * t246 - t132 * t245 + t133 * t221;
t30 = -pkin(4) * t239 - t35;
t127 = -qJD(4) * t166 + t218 * t240;
t82 = -t127 * t213 + t215 * t239;
t83 = t127 * t215 + t213 * t239;
t46 = -t82 * mrSges(6,1) + t83 * mrSges(6,2);
t303 = -m(6) * t30 - t46;
t250 = t213 ^ 2 + t215 ^ 2;
t302 = qJD(5) * m(6) * t250;
t301 = qJD(6) * t225;
t300 = 0.2e1 * m(6);
t299 = 2 * m(7);
t298 = -0.2e1 * pkin(1);
t297 = 2 * mrSges(4,1);
t296 = -2 * mrSges(3,3);
t227 = -pkin(2) * t222 - t262;
t148 = (-pkin(1) + t227) * t214;
t295 = -0.2e1 * t148;
t128 = -qJD(4) * t243 + t216 * t245 - t221 * t240;
t38 = Ifges(6,1) * t83 + Ifges(6,4) * t82 + Ifges(6,5) * t128;
t294 = t38 / 0.2e1;
t293 = t62 / 0.2e1;
t292 = t63 / 0.2e1;
t291 = t83 / 0.2e1;
t165 = t177 * qJD(6);
t100 = -Ifges(7,5) * t301 - Ifges(7,6) * t165;
t288 = t100 / 0.2e1;
t117 = Ifges(7,4) * t177 - Ifges(7,2) * t225;
t287 = t117 / 0.2e1;
t118 = Ifges(7,1) * t177 - Ifges(7,4) * t225;
t286 = t118 / 0.2e1;
t267 = Ifges(6,4) * t215;
t229 = -Ifges(6,2) * t213 + t267;
t140 = (Ifges(6,6) * t221 - t218 * t229) * qJD(4);
t285 = t140 / 0.2e1;
t268 = Ifges(6,4) * t213;
t230 = Ifges(6,1) * t215 - t268;
t141 = (Ifges(6,5) * t221 - t218 * t230) * qJD(4);
t284 = t141 / 0.2e1;
t154 = t177 * t221;
t283 = -t154 / 0.2e1;
t156 = t225 * t221;
t282 = -t156 / 0.2e1;
t281 = -t301 / 0.2e1;
t280 = -t165 / 0.2e1;
t277 = Ifges(6,2) * t305 + t268 / 0.2e1;
t274 = pkin(4) * t218;
t273 = Ifges(3,4) + Ifges(4,6);
t272 = pkin(10) + qJ(5);
t34 = t221 * t108 + t113 * t245 - t132 * t246 + t218 * t133;
t29 = (qJ(5) * t248 + qJD(5) * t219) * t214 + t34;
t244 = t216 * t275;
t202 = qJD(2) * t244;
t209 = t216 * qJD(3);
t112 = -t240 * t290 + t202 + t209;
t41 = pkin(4) * t128 - qJ(5) * t127 - qJD(5) * t167 + t112;
t12 = t213 * t41 + t215 * t29;
t67 = t218 * t113 + t221 * t132;
t58 = qJ(5) * t257 + t67;
t172 = pkin(8) * t256 + t206;
t147 = -t216 * qJ(3) - t172;
t131 = pkin(3) * t256 - t147;
t69 = pkin(4) * t166 - qJ(5) * t167 + t131;
t32 = t213 * t69 + t215 * t58;
t95 = mrSges(5,1) * t239 - mrSges(5,3) * t127;
t271 = t95 - t46;
t270 = Ifges(5,4) * t218;
t269 = Ifges(5,4) * t221;
t159 = -pkin(8) * t240 + t202;
t266 = t159 * mrSges(3,2);
t66 = t113 * t221 - t218 * t132;
t59 = -pkin(4) * t257 - t66;
t265 = t218 * t59;
t136 = mrSges(5,1) * t257 - mrSges(5,3) * t167;
t68 = -mrSges(6,1) * t125 + mrSges(6,2) * t126;
t264 = -t136 + t68;
t188 = -mrSges(6,1) * t215 + mrSges(6,2) * t213;
t263 = t188 - mrSges(5,1);
t261 = qJ(5) * t221;
t160 = t172 * qJD(2);
t260 = t160 * t219;
t259 = t213 * t218;
t258 = t213 * t221;
t255 = t215 * t218;
t254 = t215 * t221;
t253 = t218 * t223;
t252 = t221 * t223;
t157 = -qJD(5) * t221 + qJD(3) + (pkin(4) * t221 + qJ(5) * t218) * qJD(4);
t238 = t223 * t245;
t110 = t213 * t157 + t215 * t238;
t251 = Ifges(3,5) * t239 + Ifges(4,5) * t240;
t184 = qJ(3) - t261 + t274;
t143 = t213 * t184 + t215 * t253;
t23 = qJD(6) * t62 + t217 * t82 + t220 * t83;
t24 = -qJD(6) * t63 - t217 * t83 + t220 * t82;
t5 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t128;
t90 = -t165 * t221 + t225 * t246;
t92 = t177 * t246 + t221 * t301;
t43 = Ifges(7,5) * t90 + Ifges(7,6) * t92 + Ifges(7,3) * t245;
t242 = Ifges(5,5) * t127 - Ifges(5,6) * t128 + Ifges(5,3) * t239;
t237 = Ifges(7,5) * t278 + Ifges(7,6) * t279 + Ifges(6,5) * t213 / 0.2e1 + Ifges(6,6) * t305;
t9 = -t24 * mrSges(7,1) + t23 * mrSges(7,2);
t47 = -t92 * mrSges(7,1) + t90 * mrSges(7,2);
t235 = (mrSges(4,2) - mrSges(3,1)) * t160;
t234 = -t213 * t223 + pkin(5);
t233 = pkin(5) * t213 - t223;
t11 = -t213 * t29 + t215 * t41;
t31 = -t213 * t58 + t215 * t69;
t232 = t250 * mrSges(6,3);
t231 = mrSges(6,1) * t213 + mrSges(6,2) * t215;
t228 = Ifges(6,5) * t215 - Ifges(6,6) * t213;
t15 = pkin(5) * t166 - pkin(10) * t126 + t31;
t17 = pkin(10) * t125 + t32;
t3 = t15 * t220 - t17 * t217;
t4 = t15 * t217 + t17 * t220;
t226 = t66 * t218 - t67 * t221;
t174 = t215 * t184;
t106 = -pkin(10) * t254 + t218 * t234 + t174;
t114 = -pkin(10) * t258 + t143;
t56 = t106 * t220 - t114 * t217;
t57 = t106 * t217 + t114 * t220;
t187 = t272 * t213;
t189 = t272 * t215;
t129 = -t187 * t220 - t189 * t217;
t130 = -t187 * t217 + t189 * t220;
t207 = -pkin(5) * t215 - pkin(4);
t196 = Ifges(5,1) * t221 - t270;
t195 = -Ifges(5,2) * t218 + t269;
t194 = t218 * mrSges(5,1) + t221 * mrSges(5,2);
t192 = Ifges(6,1) * t213 + t267;
t183 = (-Ifges(5,1) * t218 - t269) * qJD(4);
t182 = (-Ifges(5,2) * t221 - t270) * qJD(4);
t181 = (mrSges(5,1) * t221 - mrSges(5,2) * t218) * qJD(4);
t180 = mrSges(6,1) * t218 - mrSges(6,3) * t254;
t179 = -mrSges(6,2) * t218 - mrSges(6,3) * t258;
t178 = -mrSges(4,1) * t256 - mrSges(4,3) * t216;
t175 = t233 * t221;
t171 = -t203 + t244;
t170 = (mrSges(6,1) * t221 + mrSges(6,3) * t255) * qJD(4);
t169 = (-mrSges(6,2) * t221 + mrSges(6,3) * t259) * qJD(4);
t168 = t231 * t221;
t163 = t233 * t246;
t158 = t231 * t246;
t155 = t225 * t218;
t153 = t177 * t218;
t152 = t216 * t241 + t203;
t151 = Ifges(6,5) * t218 + t221 * t230;
t150 = Ifges(6,6) * t218 + t221 * t229;
t149 = Ifges(6,3) * t218 + t221 * t228;
t146 = t215 * t157;
t144 = -t159 - t209;
t142 = -t213 * t253 + t174;
t139 = (Ifges(6,3) * t221 - t218 * t228) * qJD(4);
t138 = mrSges(7,1) * t218 + mrSges(7,3) * t156;
t137 = -mrSges(7,2) * t218 - mrSges(7,3) * t154;
t135 = -mrSges(5,2) * t257 - mrSges(5,3) * t166;
t134 = t201 + (-qJ(3) * t248 - t247) * t214;
t115 = mrSges(7,1) * t225 + mrSges(7,2) * t177;
t109 = -t213 * t238 + t146;
t103 = mrSges(5,1) * t166 + mrSges(5,2) * t167;
t102 = -Ifges(7,1) * t301 - Ifges(7,4) * t165;
t101 = -Ifges(7,4) * t301 - Ifges(7,2) * t165;
t99 = mrSges(7,1) * t165 - mrSges(7,2) * t301;
t96 = -mrSges(5,2) * t239 - mrSges(5,3) * t128;
t94 = mrSges(7,1) * t154 - mrSges(7,2) * t156;
t93 = pkin(10) * t213 * t246 + t110;
t91 = -qJD(4) * t154 + t218 * t301;
t89 = -qJD(4) * t156 - t165 * t218;
t88 = Ifges(5,1) * t167 - Ifges(5,4) * t166 + Ifges(5,5) * t257;
t87 = Ifges(5,4) * t167 - Ifges(5,2) * t166 + Ifges(5,6) * t257;
t81 = -Ifges(7,1) * t156 - Ifges(7,4) * t154 + Ifges(7,5) * t218;
t80 = -Ifges(7,4) * t156 - Ifges(7,2) * t154 + Ifges(7,6) * t218;
t79 = -Ifges(7,5) * t156 - Ifges(7,6) * t154 + Ifges(7,3) * t218;
t77 = -qJD(5) * t177 - qJD(6) * t130;
t76 = -qJD(5) * t225 + qJD(6) * t129;
t75 = t146 + (pkin(10) * t255 + t221 * t234) * qJD(4);
t74 = mrSges(6,1) * t166 - mrSges(6,3) * t126;
t73 = -mrSges(6,2) * t166 + mrSges(6,3) * t125;
t72 = -mrSges(7,2) * t245 + mrSges(7,3) * t92;
t71 = mrSges(7,1) * t245 - mrSges(7,3) * t90;
t70 = mrSges(5,1) * t128 + mrSges(5,2) * t127;
t61 = Ifges(5,1) * t127 - Ifges(5,4) * t128 + Ifges(5,5) * t239;
t60 = Ifges(5,4) * t127 - Ifges(5,2) * t128 + Ifges(5,6) * t239;
t54 = mrSges(6,1) * t128 - mrSges(6,3) * t83;
t53 = -mrSges(6,2) * t128 + mrSges(6,3) * t82;
t52 = Ifges(6,1) * t126 + Ifges(6,4) * t125 + Ifges(6,5) * t166;
t51 = Ifges(6,4) * t126 + Ifges(6,2) * t125 + Ifges(6,6) * t166;
t49 = mrSges(7,1) * t166 - mrSges(7,3) * t63;
t48 = -mrSges(7,2) * t166 + mrSges(7,3) * t62;
t45 = Ifges(7,1) * t90 + Ifges(7,4) * t92 + Ifges(7,5) * t245;
t44 = Ifges(7,4) * t90 + Ifges(7,2) * t92 + Ifges(7,6) * t245;
t42 = -pkin(5) * t125 + t59;
t37 = Ifges(6,4) * t83 + Ifges(6,2) * t82 + Ifges(6,6) * t128;
t36 = Ifges(6,5) * t83 + Ifges(6,6) * t82 + Ifges(6,3) * t128;
t33 = -mrSges(7,1) * t62 + mrSges(7,2) * t63;
t28 = Ifges(7,1) * t63 + Ifges(7,4) * t62 + Ifges(7,5) * t166;
t27 = Ifges(7,4) * t63 + Ifges(7,2) * t62 + Ifges(7,6) * t166;
t19 = -qJD(6) * t57 - t217 * t93 + t220 * t75;
t18 = qJD(6) * t56 + t217 * t75 + t220 * t93;
t16 = -pkin(5) * t82 + t30;
t14 = -mrSges(7,2) * t128 + mrSges(7,3) * t24;
t13 = mrSges(7,1) * t128 - mrSges(7,3) * t23;
t10 = pkin(10) * t82 + t12;
t8 = pkin(5) * t128 - pkin(10) * t83 + t11;
t7 = Ifges(7,1) * t23 + Ifges(7,4) * t24 + Ifges(7,5) * t128;
t6 = Ifges(7,4) * t23 + Ifges(7,2) * t24 + Ifges(7,6) * t128;
t2 = -qJD(6) * t4 - t10 * t217 + t220 * t8;
t1 = qJD(6) * t3 + t10 * t220 + t217 * t8;
t20 = [(t1 * t4 + t16 * t42 + t2 * t3) * t299 + (t11 * t31 + t12 * t32 + t30 * t59) * t300 + (t5 + t36 - t60) * t166 + (t219 * t242 + t260 * t297 + 0.2e1 * t134 * (mrSges(4,2) * t222 - mrSges(4,3) * t219) + 0.2e1 * (t159 * t222 + t260) * mrSges(3,3) + ((t147 * t297 + mrSges(4,2) * t295 + t172 * t296 + (Ifges(4,5) - (2 * Ifges(3,6))) * t216 + (mrSges(3,1) * t298 - 0.2e1 * t219 * t273) * t214) * t219 + (t152 * t297 + t171 * t296 + mrSges(4,3) * t295 + Ifges(5,5) * t167 - Ifges(5,6) * t166 + (Ifges(3,5) - (2 * Ifges(4,4))) * t216 + (mrSges(3,2) * t298 + 0.2e1 * t222 * t273) * t214 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3)) * t257) * t222) * qJD(2)) * t214 + (0.2e1 * t235 + t251 - 0.2e1 * t266) * t216 + (-t87 + t304) * t128 + 0.2e1 * t3 * t13 + 0.2e1 * t4 * t14 + t24 * t27 + t23 * t28 + 0.2e1 * t16 * t33 + 0.2e1 * t42 * t9 + 0.2e1 * t1 * t48 + 0.2e1 * t2 * t49 + 0.2e1 * t32 * t53 + 0.2e1 * t31 * t54 + 0.2e1 * t59 * t46 + t62 * t6 + t63 * t7 + 0.2e1 * t30 * t68 + 0.2e1 * m(3) * (t159 * t172 - t160 * t171) + 0.2e1 * m(4) * (t134 * t148 + t144 * t147 + t152 * t160) + 0.2e1 * m(5) * (t112 * t131 + t34 * t67 + t35 * t66) + 0.2e1 * t12 * t73 + 0.2e1 * t11 * t74 + t82 * t51 + t83 * t52 + 0.2e1 * t66 * t95 + 0.2e1 * t67 * t96 + 0.2e1 * t112 * t103 + t125 * t37 + t126 * t38 + t127 * t88 + 0.2e1 * t131 * t70 + 0.2e1 * t34 * t135 + 0.2e1 * t35 * t136 + t167 * t61 + 0.2e1 * t144 * t178; m(4) * (-pkin(2) * t160 - qJ(3) * t144 - qJD(3) * t147) + (t215 * t294 - t213 * t37 / 0.2e1 - t35 * mrSges(5,3) + t61 / 0.2e1 + t271 * t223) * t221 + (t43 / 0.2e1 + t139 / 0.2e1 - t182 / 0.2e1) * t166 + m(5) * (qJ(3) * t112 + qJD(3) * t131 + t252 * t35 + t253 * t34) + m(6) * (t109 * t31 + t142 * t11 + t110 * t32 + t143 * t12 - t252 * t30) + t45 * t292 + t44 * t293 + t7 * t282 + t6 * t283 + t126 * t284 + t125 * t285 + t151 * t291 + (t79 / 0.2e1 + t149 / 0.2e1 - t195 / 0.2e1) * t128 + (-Ifges(3,6) * t219 - Ifges(4,4) * t222 + t222 * (Ifges(5,5) * t221 - Ifges(5,6) * t218) / 0.2e1 + t227 * mrSges(4,1)) * t249 - t266 + m(7) * (t1 * t57 + t16 * t175 - t163 * t42 + t18 * t4 + t19 * t3 + t2 * t56) + t251 + t235 + (t103 - t178) * qJD(3) + ((-Ifges(5,5) * t218 - Ifges(5,6) * t221) * t257 / 0.2e1 - t221 * t87 / 0.2e1 - t218 * t88 / 0.2e1 + t51 * t259 / 0.2e1 - t52 * t255 / 0.2e1 + t226 * mrSges(5,3) + (-m(5) * t226 + m(6) * t265 + t221 * t135 + t218 * t264) * t223 + t304 * t221 / 0.2e1) * qJD(4) + t42 * t47 + t18 * t48 + t19 * t49 + t56 * t13 + t57 * t14 + qJ(3) * t70 + t3 * t71 + t4 * t72 + t24 * t80 / 0.2e1 + t23 * t81 / 0.2e1 + t90 * t28 / 0.2e1 + t92 * t27 / 0.2e1 + t16 * t94 + t109 * t74 + t110 * t73 + t1 * t137 + t2 * t138 + t142 * t54 + t143 * t53 - t144 * mrSges(4,3) + t82 * t150 / 0.2e1 - t59 * t158 - t163 * t33 + t30 * t168 + t32 * t169 + t31 * t170 + t175 * t9 + t12 * t179 + t11 * t180 + t131 * t181 + t167 * t183 / 0.2e1 + t112 * t194 + t127 * t196 / 0.2e1 + (t223 * t96 - t34 * mrSges(5,3) + t5 / 0.2e1 + t36 / 0.2e1 - t60 / 0.2e1) * t218; 0.2e1 * qJ(3) * t181 + 0.2e1 * t109 * t180 + 0.2e1 * t110 * t179 + 0.2e1 * t18 * t137 + 0.2e1 * t19 * t138 + 0.2e1 * t142 * t170 + 0.2e1 * t143 * t169 - t154 * t44 - t156 * t45 - 0.2e1 * t163 * t94 + 0.2e1 * t175 * t47 + 0.2e1 * t56 * t71 + 0.2e1 * t57 * t72 + t92 * t80 + t90 * t81 + (t142 * t109 + t143 * t110) * t300 + (-t163 * t175 + t18 * t57 + t19 * t56) * t299 + (t43 + t139 - t182) * t218 + (-t213 * t140 + t215 * t141 + t158 * t306 + t183) * t221 + 0.2e1 * (mrSges(4,3) + t194 + (m(4) + m(5)) * qJ(3)) * qJD(3) + ((t149 - t195 + t79) * t221 + (t213 * t150 - t215 * t151 - t196 + (-m(6) * t252 + t168) * t306) * t218) * qJD(4); mrSges(4,1) * t239 - t153 * t13 - t155 * t14 + t89 * t48 + t91 * t49 + (-t9 + t271) * t221 + (-t213 * t54 + t215 * t53 + t96) * t218 + ((-t213 * t74 + t215 * t73 + t135) * t221 + (t33 + t264) * t218) * qJD(4) + m(7) * (-t1 * t155 - t153 * t2 - t16 * t221 + t246 * t42 + t3 * t91 + t4 * t89) + m(6) * (-t221 * t30 + (-t11 * t213 + t12 * t215) * t218 + (t265 + (-t213 * t31 + t215 * t32) * t221) * qJD(4)) + m(5) * (-qJD(4) * t226 + t34 * t218 + t35 * t221) + m(4) * t160; t89 * t137 + t91 * t138 - t153 * t71 - t155 * t72 + (t158 - t47) * t221 + (t215 * t169 - t213 * t170) * t218 + ((t215 * t179 - t213 * t180) * t221 + (t168 + t94) * t218) * qJD(4) + m(7) * (-t153 * t19 - t155 * t18 + t163 * t221 + t175 * t246 + t56 * t91 + t57 * t89) + m(6) * ((-t109 * t213 + t110 * t215) * t218 + (-t142 * t213 + t143 * t215 - 0.2e1 * t253) * t245); (-t153 * t91 - t155 * t89) * t299 + 0.4e1 * (m(6) * (-0.1e1 + t250) / 0.2e1 - m(7) / 0.2e1) * t218 * t245; t242 + (-t1 * t225 - t165 * t4 - t177 * t2 + t3 * t301) * mrSges(7,3) + (m(6) * (-qJ(5) * t11 - qJD(5) * t31) - qJD(5) * t74 - qJ(5) * t54 - t11 * mrSges(6,3) + t294) * t213 + (m(6) * (qJ(5) * t12 + qJD(5) * t32) + qJD(5) * t73 + qJ(5) * t53 + t12 * mrSges(6,3) + t37 / 0.2e1) * t215 + t237 * t128 + t101 * t293 + t82 * t277 + t7 * t278 + t6 * t279 + t27 * t280 + t28 * t281 + t23 * t286 + t24 * t287 + t166 * t288 + t192 * t291 + t102 * t292 + m(7) * (t1 * t130 + t129 * t2 + t16 * t207 + t3 * t77 + t4 * t76) + t303 * pkin(4) - t34 * mrSges(5,2) + t35 * mrSges(5,1) + t76 * t48 + t77 * t49 + t42 * t99 + t16 * t115 + t129 * t13 + t130 * t14 + t30 * t188 + t207 * t9; m(7) * (t129 * t19 + t130 * t18 - t163 * t207 + t56 * t77 + t57 * t76) + t92 * t287 + t90 * t286 + t129 * t71 + t130 * t72 + t76 * t137 + t77 * t138 + t101 * t283 + t102 * t282 + pkin(4) * t158 - t163 * t115 + t81 * t281 + t80 * t280 + t175 * t99 + t44 * t279 + t45 * t278 + t207 * t47 + t218 * t288 + (m(6) * (qJ(5) * t110 + qJD(5) * t143) + qJ(5) * t169 + qJD(5) * t179 + t110 * mrSges(6,3) + t285) * t215 + (m(6) * (-qJ(5) * t109 - qJD(5) * t142) - qJ(5) * t170 - qJD(5) * t180 - t109 * mrSges(6,3) + t284) * t213 + (-t165 * t57 - t177 * t19 - t18 * t225 + t301 * t56) * mrSges(7,3) + ((-t223 * mrSges(5,2) - Ifges(5,6) + t237) * t221 + (-Ifges(5,5) + t213 * t277 - t215 * t192 / 0.2e1 + (-m(6) * pkin(4) + t263) * t223) * t218) * qJD(4); -t221 * t99 + m(7) * (t129 * t91 + t130 * t89 - t153 * t77 - t155 * t76) + t218 * t302 + ((-mrSges(5,2) + t232) * t221 + m(6) * (t250 * t261 - t274) + (m(7) * t207 + t115 + t263) * t218) * qJD(4) + (-t153 * t301 + t155 * t165 - t177 * t91 - t225 * t89) * mrSges(7,3); -t301 * t118 + t177 * t102 - t165 * t117 - t225 * t101 + (t129 * t77 + t130 * t76) * t299 + 0.2e1 * t207 * t99 + 0.2e1 * qJ(5) * t302 + 0.2e1 * t232 * qJD(5) + 0.2e1 * (t129 * t301 - t130 * t165 - t177 * t77 - t225 * t76) * mrSges(7,3); m(7) * t16 - t303 + t9; -m(7) * t163 + (m(6) * t223 - t231) * t246 + t47; (m(6) + m(7)) * t246; t99; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t19 - mrSges(7,2) * t18 + t43; mrSges(7,1) * t91 - mrSges(7,2) * t89; mrSges(7,1) * t77 - t76 * mrSges(7,2) + t100; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
