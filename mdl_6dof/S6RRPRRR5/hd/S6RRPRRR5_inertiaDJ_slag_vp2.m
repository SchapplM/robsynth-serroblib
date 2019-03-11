% Calculate time derivative of joint inertia matrix for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:07
% EndTime: 2019-03-09 13:38:23
% DurationCPUTime: 7.14s
% Computational Cost: add. (15030->699), mult. (40693->1017), div. (0->0), fcn. (40908->12), ass. (0->287)
t233 = sin(pkin(6));
t349 = 0.2e1 * t233;
t237 = sin(qJ(5));
t241 = cos(qJ(5));
t238 = sin(qJ(4));
t277 = qJD(5) * t238;
t242 = cos(qJ(4));
t279 = qJD(4) * t242;
t247 = -t237 * t277 + t241 * t279;
t234 = cos(pkin(12));
t224 = -pkin(2) * t234 - pkin(3);
t188 = -pkin(4) * t242 - pkin(10) * t238 + t224;
t232 = sin(pkin(12));
t223 = pkin(2) * t232 + pkin(9);
t286 = t241 * t242;
t203 = t223 * t286;
t144 = t237 * t188 + t203;
t348 = qJD(5) * t144;
t236 = sin(qJ(6));
t240 = cos(qJ(6));
t252 = t236 * t237 - t240 * t241;
t320 = -t252 / 0.2e1;
t192 = t236 * t241 + t237 * t240;
t319 = t192 / 0.2e1;
t347 = t237 / 0.2e1;
t315 = t241 / 0.2e1;
t276 = qJD(5) * t241;
t246 = t237 * t279 + t238 * t276;
t139 = t246 * mrSges(6,1) + t247 * mrSges(6,2);
t344 = qJD(5) + qJD(6);
t146 = t344 * t192;
t105 = -t146 * t238 - t252 * t279;
t179 = t252 * t238;
t106 = t179 * t344 - t192 * t279;
t67 = -t106 * mrSges(7,1) + t105 * mrSges(7,2);
t346 = -t139 - t67;
t284 = t237 ^ 2 + t241 ^ 2;
t239 = sin(qJ(2));
t235 = cos(pkin(6));
t314 = pkin(1) * t235;
t220 = t239 * t314;
t243 = cos(qJ(2));
t290 = t233 * t243;
t291 = t233 * t239;
t307 = -pkin(8) - qJ(3);
t345 = (t290 * t307 - t220) * qJD(2) - qJD(3) * t291;
t171 = (t232 * t243 + t234 * t239) * t233;
t147 = t171 * t238 - t235 * t242;
t283 = qJD(2) * t233;
t289 = t234 * t243;
t169 = (-t232 * t239 + t289) * t283;
t118 = -qJD(4) * t147 + t169 * t242;
t148 = t171 * t242 + t235 * t238;
t170 = t232 * t291 - t233 * t289;
t124 = t148 * t241 + t170 * t237;
t168 = qJD(2) * t171;
t55 = -qJD(5) * t124 - t118 * t237 + t168 * t241;
t123 = -t148 * t237 + t170 * t241;
t56 = qJD(5) * t123 + t118 * t241 + t168 * t237;
t26 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t267 = t239 * t283;
t258 = pkin(2) * t267;
t122 = pkin(3) * t168 - pkin(9) * t169 + t258;
t193 = (-pkin(2) * t243 - pkin(1)) * t233;
t125 = t170 * pkin(3) - t171 * pkin(9) + t193;
t281 = qJD(4) * t238;
t221 = t243 * t314;
t217 = qJD(2) * t221;
t261 = t307 * t239;
t138 = t217 + (qJD(2) * t261 + qJD(3) * t243) * t233;
t88 = t234 * t138 + t232 * t345;
t154 = pkin(2) * t235 + t233 * t261 + t221;
t186 = pkin(8) * t290 + t220;
t167 = qJ(3) * t290 + t186;
t120 = t232 * t154 + t234 * t167;
t99 = pkin(9) * t235 + t120;
t36 = t122 * t242 - t125 * t281 - t238 * t88 - t99 * t279;
t70 = t238 * t125 + t242 * t99;
t83 = mrSges(5,1) * t168 - mrSges(5,3) * t118;
t343 = m(5) * (qJD(4) * t70 + t36) - t26 + t83;
t342 = 2 * m(6);
t341 = 2 * m(7);
t340 = -2 * mrSges(3,3);
t339 = -2 * mrSges(4,3);
t338 = -2 * Ifges(4,4);
t337 = 0.2e1 * t193;
t336 = 0.2e1 * t223;
t335 = m(4) * pkin(2);
t334 = t55 / 0.2e1;
t71 = t123 * t240 - t124 * t236;
t333 = t71 / 0.2e1;
t72 = t123 * t236 + t124 * t240;
t332 = t72 / 0.2e1;
t331 = -pkin(11) - pkin(10);
t29 = -pkin(4) * t168 - t36;
t330 = m(6) * t29;
t303 = Ifges(6,4) * t237;
t209 = Ifges(6,2) * t241 + t303;
t302 = Ifges(6,4) * t241;
t255 = -Ifges(6,2) * t237 + t302;
t133 = -t209 * t277 + (Ifges(6,6) * t238 + t242 * t255) * qJD(4);
t329 = t133 / 0.2e1;
t211 = Ifges(6,1) * t237 + t302;
t256 = Ifges(6,1) * t241 - t303;
t134 = -t211 * t277 + (Ifges(6,5) * t238 + t242 * t256) * qJD(4);
t328 = t134 / 0.2e1;
t145 = t344 * t252;
t327 = -t145 / 0.2e1;
t326 = -t146 / 0.2e1;
t151 = Ifges(7,4) * t192 - Ifges(7,2) * t252;
t325 = t151 / 0.2e1;
t152 = Ifges(7,1) * t192 - Ifges(7,4) * t252;
t324 = t152 / 0.2e1;
t175 = -Ifges(6,5) * t242 + t238 * t256;
t323 = t175 / 0.2e1;
t178 = t192 * t238;
t322 = -t178 / 0.2e1;
t321 = -t179 / 0.2e1;
t318 = t211 / 0.2e1;
t317 = -t237 / 0.2e1;
t316 = -t241 / 0.2e1;
t313 = pkin(4) * t238;
t312 = pkin(10) * t242;
t87 = t138 * t232 - t234 * t345;
t311 = t87 * mrSges(4,1);
t310 = t87 * mrSges(5,1);
t309 = t87 * mrSges(5,2);
t308 = t88 * mrSges(4,2);
t58 = pkin(10) * t170 + t70;
t119 = t154 * t234 - t232 * t167;
t98 = -pkin(3) * t235 - t119;
t68 = pkin(4) * t147 - pkin(10) * t148 + t98;
t31 = t237 * t68 + t241 * t58;
t305 = Ifges(5,4) * t238;
t304 = Ifges(5,4) * t242;
t301 = Ifges(6,6) * t237;
t300 = t168 * Ifges(5,5);
t299 = t168 * Ifges(5,6);
t298 = t170 * Ifges(5,6);
t176 = -pkin(8) * t267 + t217;
t297 = t176 * mrSges(3,2);
t177 = t186 * qJD(2);
t296 = t177 * mrSges(3,1);
t127 = mrSges(5,1) * t170 - mrSges(5,3) * t148;
t74 = -mrSges(6,1) * t123 + mrSges(6,2) * t124;
t295 = -t127 + t74;
t207 = -mrSges(6,1) * t241 + mrSges(6,2) * t237;
t294 = t207 - mrSges(5,1);
t69 = t125 * t242 - t238 * t99;
t57 = -pkin(4) * t170 - t69;
t293 = qJD(4) * t57;
t292 = t223 * t237;
t288 = t237 * t238;
t287 = t238 * t241;
t95 = -Ifges(7,5) * t145 - Ifges(7,6) * t146;
t206 = (-t312 + t313) * qJD(4);
t285 = t241 * t206 + t281 * t292;
t282 = qJD(4) * t223;
t280 = qJD(4) * t241;
t278 = qJD(5) * t237;
t275 = qJD(6) * t236;
t274 = qJD(6) * t240;
t117 = qJD(4) * t148 + t169 * t238;
t19 = qJD(6) * t71 + t236 * t55 + t240 * t56;
t20 = -qJD(6) * t72 - t236 * t56 + t240 * t55;
t6 = Ifges(7,5) * t19 + Ifges(7,6) * t20 + Ifges(7,3) * t117;
t22 = Ifges(6,5) * t56 + Ifges(6,6) * t55 + Ifges(6,3) * t117;
t273 = pkin(5) * t278;
t51 = Ifges(6,1) * t124 + Ifges(6,4) * t123 + Ifges(6,5) * t147;
t272 = t51 * t315;
t228 = Ifges(6,5) * t276;
t271 = -Ifges(6,6) * t278 / 0.2e1 + t228 / 0.2e1 + t95 / 0.2e1;
t61 = Ifges(7,5) * t105 + Ifges(7,6) * t106 + Ifges(7,3) * t281;
t270 = Ifges(5,5) * t118 - Ifges(5,6) * t117 + Ifges(5,3) * t168;
t269 = Ifges(3,5) * t243 * t283 + Ifges(4,5) * t169 - Ifges(4,6) * t168;
t268 = qJD(5) * t331;
t265 = t242 * t278;
t262 = Ifges(6,5) * t347 + Ifges(7,5) * t319 + Ifges(6,6) * t315 + Ifges(7,6) * t320;
t30 = -t237 * t58 + t241 * t68;
t103 = t188 * t276 + t237 * t206 + (-t238 * t280 - t265) * t223;
t181 = t241 * t188;
t143 = -t242 * t292 + t181;
t259 = -qJD(5) * t143 + t103;
t257 = mrSges(6,1) * t237 + mrSges(6,2) * t241;
t35 = t238 * t122 + t125 * t279 + t242 * t88 - t281 * t99;
t28 = pkin(10) * t168 + t35;
t46 = pkin(4) * t117 - pkin(10) * t118 + t87;
t10 = t237 * t46 + t241 * t28 + t68 * t276 - t278 * t58;
t11 = -qJD(5) * t31 - t237 * t28 + t241 * t46;
t254 = t10 * t241 - t11 * t237;
t21 = pkin(5) * t147 - pkin(11) * t124 + t30;
t25 = pkin(11) * t123 + t31;
t12 = t21 * t240 - t236 * t25;
t13 = t21 * t236 + t240 * t25;
t41 = -mrSges(6,2) * t117 + mrSges(6,3) * t55;
t42 = mrSges(6,1) * t117 - mrSges(6,3) * t56;
t253 = -t237 * t42 + t241 * t41;
t131 = -pkin(11) * t287 + t181 + (-pkin(5) - t292) * t242;
t137 = -pkin(11) * t288 + t144;
t77 = t131 * t240 - t137 * t236;
t78 = t131 * t236 + t137 * t240;
t214 = t331 * t237;
t215 = t331 * t241;
t157 = t214 * t240 + t215 * t236;
t158 = t214 * t236 - t215 * t240;
t204 = t237 * t268;
t205 = t241 * t268;
t113 = qJD(6) * t157 + t204 * t240 + t205 * t236;
t114 = -qJD(6) * t158 - t204 * t236 + t205 * t240;
t251 = t114 * mrSges(7,1) - t113 * mrSges(7,2) + t95;
t4 = pkin(5) * t117 - pkin(11) * t56 + t11;
t5 = pkin(11) * t55 + t10;
t2 = qJD(6) * t12 + t236 * t4 + t240 * t5;
t3 = -qJD(6) * t13 - t236 * t5 + t240 * t4;
t250 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t82 = -mrSges(5,2) * t168 - mrSges(5,3) * t117;
t249 = t82 + m(5) * (-qJD(4) * t69 + t35);
t81 = (pkin(5) * t238 - pkin(11) * t286) * qJD(4) + (-t203 + (pkin(11) * t238 - t188) * t237) * qJD(5) + t285;
t86 = -pkin(11) * t246 + t103;
t38 = qJD(6) * t77 + t236 * t81 + t240 * t86;
t39 = -qJD(6) * t78 - t236 * t86 + t240 * t81;
t248 = t39 * mrSges(7,1) - t38 * mrSges(7,2) + t61;
t244 = -t276 * t30 - t278 * t31 + t254;
t132 = Ifges(6,5) * t247 - Ifges(6,6) * t246 + Ifges(6,3) * t281;
t229 = Ifges(5,5) * t279;
t225 = -pkin(5) * t241 - pkin(4);
t212 = Ifges(5,1) * t238 + t304;
t210 = Ifges(5,2) * t242 + t305;
t202 = -mrSges(6,1) * t242 - mrSges(6,3) * t287;
t201 = mrSges(6,2) * t242 - mrSges(6,3) * t288;
t200 = (Ifges(5,1) * t242 - t305) * qJD(4);
t199 = t256 * qJD(5);
t198 = (-Ifges(5,2) * t238 + t304) * qJD(4);
t197 = t255 * qJD(5);
t195 = (t238 * mrSges(5,1) + t242 * mrSges(5,2)) * qJD(4);
t194 = t257 * qJD(5);
t189 = (-mrSges(7,1) * t236 - mrSges(7,2) * t240) * qJD(6) * pkin(5);
t187 = t257 * t238;
t185 = -pkin(8) * t291 + t221;
t182 = (pkin(5) * t237 + t223) * t238;
t174 = -Ifges(6,6) * t242 + t238 * t255;
t173 = -Ifges(6,3) * t242 + (Ifges(6,5) * t241 - t301) * t238;
t166 = -mrSges(6,2) * t281 - mrSges(6,3) * t246;
t165 = mrSges(6,1) * t281 - mrSges(6,3) * t247;
t161 = t169 * mrSges(4,2);
t160 = -mrSges(7,1) * t242 + mrSges(7,3) * t179;
t159 = mrSges(7,2) * t242 - mrSges(7,3) * t178;
t153 = pkin(5) * t246 + t223 * t279;
t149 = mrSges(7,1) * t252 + mrSges(7,2) * t192;
t135 = mrSges(7,1) * t178 - mrSges(7,2) * t179;
t130 = -Ifges(7,1) * t179 - Ifges(7,4) * t178 - Ifges(7,5) * t242;
t129 = -Ifges(7,4) * t179 - Ifges(7,2) * t178 - Ifges(7,6) * t242;
t128 = -Ifges(7,5) * t179 - Ifges(7,6) * t178 - Ifges(7,3) * t242;
t126 = -mrSges(5,2) * t170 - mrSges(5,3) * t147;
t104 = t285 - t348;
t97 = -Ifges(7,1) * t145 - Ifges(7,4) * t146;
t96 = -Ifges(7,4) * t145 - Ifges(7,2) * t146;
t94 = mrSges(7,1) * t146 - mrSges(7,2) * t145;
t90 = -mrSges(7,2) * t281 + mrSges(7,3) * t106;
t89 = mrSges(7,1) * t281 - mrSges(7,3) * t105;
t80 = Ifges(5,1) * t148 - Ifges(5,4) * t147 + Ifges(5,5) * t170;
t79 = Ifges(5,4) * t148 - Ifges(5,2) * t147 + t298;
t76 = mrSges(6,1) * t147 - mrSges(6,3) * t124;
t75 = -mrSges(6,2) * t147 + mrSges(6,3) * t123;
t73 = mrSges(5,1) * t117 + mrSges(5,2) * t118;
t63 = Ifges(7,1) * t105 + Ifges(7,4) * t106 + Ifges(7,5) * t281;
t62 = Ifges(7,4) * t105 + Ifges(7,2) * t106 + Ifges(7,6) * t281;
t60 = Ifges(5,1) * t118 - Ifges(5,4) * t117 + t300;
t59 = Ifges(5,4) * t118 - Ifges(5,2) * t117 + t299;
t50 = Ifges(6,4) * t124 + Ifges(6,2) * t123 + Ifges(6,6) * t147;
t49 = Ifges(6,5) * t124 + Ifges(6,6) * t123 + Ifges(6,3) * t147;
t48 = mrSges(7,1) * t147 - mrSges(7,3) * t72;
t47 = -mrSges(7,2) * t147 + mrSges(7,3) * t71;
t43 = -pkin(5) * t123 + t57;
t40 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t34 = Ifges(7,1) * t72 + Ifges(7,4) * t71 + Ifges(7,5) * t147;
t33 = Ifges(7,4) * t72 + Ifges(7,2) * t71 + Ifges(7,6) * t147;
t32 = Ifges(7,5) * t72 + Ifges(7,6) * t71 + Ifges(7,3) * t147;
t24 = Ifges(6,1) * t56 + Ifges(6,4) * t55 + Ifges(6,5) * t117;
t23 = Ifges(6,4) * t56 + Ifges(6,2) * t55 + Ifges(6,6) * t117;
t16 = -pkin(5) * t55 + t29;
t15 = -mrSges(7,2) * t117 + mrSges(7,3) * t20;
t14 = mrSges(7,1) * t117 - mrSges(7,3) * t19;
t9 = -mrSges(7,1) * t20 + mrSges(7,2) * t19;
t8 = Ifges(7,1) * t19 + Ifges(7,4) * t20 + Ifges(7,5) * t117;
t7 = Ifges(7,4) * t19 + Ifges(7,2) * t20 + Ifges(7,6) * t117;
t1 = [t161 * t337 + (t12 * t3 + t13 * t2 + t16 * t43) * t341 + (t10 * t31 + t11 * t30 + t29 * t57) * t342 + (t269 - 0.2e1 * t296 - 0.2e1 * t297 - 0.2e1 * t308 - 0.2e1 * t311) * t235 + 0.2e1 * m(3) * (t176 * t186 - t177 * t185) + 0.2e1 * m(4) * (-t119 * t87 + t120 * t88) + (mrSges(4,1) * t337 + t120 * t339 + t171 * t338 + Ifges(5,5) * t148 - Ifges(4,6) * t235 - Ifges(5,6) * t147 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t170) * t168 + (0.2e1 * Ifges(4,1) * t171 + Ifges(4,5) * t235 + t119 * t339 + t170 * t338) * t169 + (t60 + 0.2e1 * t309) * t148 + (t22 - t59 + t6 + 0.2e1 * t310) * t147 + 0.2e1 * (-t170 * t88 + t171 * t87) * mrSges(4,3) + t170 * t270 + 0.2e1 * m(5) * (t35 * t70 + t36 * t69 + t87 * t98) + t118 * t80 + t123 * t23 + t124 * t24 + 0.2e1 * t35 * t126 + 0.2e1 * t36 * t127 + 0.2e1 * t98 * t73 + 0.2e1 * t70 * t82 + 0.2e1 * t69 * t83 + 0.2e1 * t10 * t75 + 0.2e1 * t11 * t76 + t71 * t7 + t72 * t8 + 0.2e1 * t29 * t74 + t56 * t51 + 0.2e1 * t57 * t26 + 0.2e1 * t3 * t48 + t55 * t50 + 0.2e1 * t2 * t47 + 0.2e1 * t16 * t40 + 0.2e1 * t31 * t41 + 0.2e1 * t30 * t42 + 0.2e1 * t43 * t9 + t20 * t33 + t19 * t34 + 0.2e1 * t12 * t14 + 0.2e1 * t13 * t15 + (t49 + t32 - t79) * t117 + (0.2e1 * (t176 * t243 + t177 * t239) * mrSges(3,3) + ((t185 * t340 + Ifges(3,5) * t235 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t243) * t349) * t243 + (-0.2e1 * Ifges(3,6) * t235 + t335 * t337 + 0.2e1 * pkin(2) * (mrSges(4,1) * t170 + mrSges(4,2) * t171) + t186 * t340 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t239 + (Ifges(3,1) - Ifges(3,2)) * t243) * t349) * t239) * qJD(2)) * t233; -t308 - t311 + t8 * t321 + t7 * t322 + t56 * t323 + t124 * t328 + t123 * t329 + t63 * t332 + t62 * t333 + t174 * t334 + (t232 * t88 - t234 * t87) * t335 - t296 - t297 + t269 + (m(5) * t87 + t73) * t224 + (t300 / 0.2e1 + t309 + t60 / 0.2e1 + t24 * t315 + t23 * t317 - t36 * mrSges(5,3) + (t316 * t50 + t317 * t51) * qJD(5) + (-t70 * mrSges(5,3) - t298 / 0.2e1 - t79 / 0.2e1 + t49 / 0.2e1 + t32 / 0.2e1) * qJD(4) + (-qJD(4) * t126 + t330 - t343) * t223) * t238 + m(6) * (t10 * t144 + t103 * t31 + t104 * t30 + t11 * t143) + (t299 / 0.2e1 - t310 + t59 / 0.2e1 - t22 / 0.2e1 - t6 / 0.2e1 + t35 * mrSges(5,3) + (t272 + t50 * t317 - t69 * mrSges(5,3) + t80 / 0.2e1) * qJD(4) + (m(6) * t293 + qJD(4) * t295 + t249) * t223) * t242 - Ifges(3,6) * t267 + t11 * t202 + t118 * t212 / 0.2e1 + t98 * t195 + t148 * t200 / 0.2e1 + t10 * t201 + t182 * t9 + t29 * t187 + t30 * t165 + t31 * t166 + t153 * t40 + t2 * t159 + t3 * t160 + t57 * t139 + t143 * t42 + t144 * t41 + t20 * t129 / 0.2e1 + t19 * t130 / 0.2e1 + t16 * t135 + t105 * t34 / 0.2e1 + t106 * t33 / 0.2e1 + t103 * t75 + t104 * t76 + t12 * t89 + t13 * t90 + t77 * t14 + t78 * t15 + t43 * t67 + t39 * t48 + t38 * t47 + (-t168 * t232 - t169 * t234) * pkin(2) * mrSges(4,3) + m(7) * (t12 * t39 + t13 * t38 + t153 * t43 + t16 * t182 + t2 * t78 + t3 * t77) + (-t198 / 0.2e1 + t132 / 0.2e1 + t61 / 0.2e1) * t147 + (-t210 / 0.2e1 + t173 / 0.2e1 + t128 / 0.2e1) * t117 + t170 * t229 / 0.2e1; 0.2e1 * t103 * t201 + 0.2e1 * t104 * t202 + t105 * t130 + t106 * t129 + 0.2e1 * t153 * t135 + 0.2e1 * t143 * t165 + 0.2e1 * t144 * t166 + 0.2e1 * t38 * t159 + 0.2e1 * t39 * t160 - t178 * t62 - t179 * t63 + 0.2e1 * t182 * t67 + 0.2e1 * t224 * t195 + 0.2e1 * t77 * t89 + 0.2e1 * t78 * t90 + (t103 * t144 + t104 * t143) * t342 + (t153 * t182 + t38 * t78 + t39 * t77) * t341 + (-t132 + t198 - t61 + (-t174 * t237 + t175 * t241 + t187 * t336 + t212) * qJD(4)) * t242 + (-t237 * t133 + t241 * t134 + t139 * t336 + t200 + (-t174 * t241 - t175 * t237) * qJD(5) + (t223 ^ 2 * t242 * t342 + t128 + t173 - t210) * qJD(4)) * t238; m(7) * (t105 * t13 + t106 * t12 - t178 * t3 - t179 * t2) + t105 * t47 - t179 * t15 + t106 * t48 - t178 * t14 + m(4) * t258 + t161 + t168 * mrSges(4,1) + (-t9 + (-t237 * t76 + t241 * t75 + t126) * qJD(4) - m(7) * t16 + m(6) * (-qJD(4) * t237 * t30 + t280 * t31 - t29) + t343) * t242 + ((-t237 * t75 - t241 * t76) * qJD(5) + m(6) * (t244 + t293) + t249 + t253 + (m(7) * t43 + t295 + t40) * qJD(4)) * t238; t105 * t159 + t106 * t160 - t178 * t89 - t179 * t90 + t346 * t242 + (t135 + t187) * t281 + m(7) * (t105 * t78 + t106 * t77 - t153 * t242 - t178 * t39 - t179 * t38 + t182 * t281) + m(6) * (t238 ^ 2 - t242 ^ 2) * t282 + (m(6) * (t103 * t238 - t143 * t277 + t144 * t279) + t201 * t279 + t238 * t166 - t202 * t277) * t241 + (m(6) * (-t104 * t238 - t143 * t279 - t144 * t277) - t201 * t277 - t202 * t279 - t238 * t165) * t237; (-t105 * t179 - t106 * t178) * t341 + 0.4e1 * (m(6) * (-0.1e1 + t284) / 0.2e1 - m(7) / 0.2e1) * t238 * t279; (t12 * t145 - t13 * t146 - t192 * t3 - t2 * t252) * mrSges(7,3) + t56 * t318 + t8 * t319 + t7 * t320 + t19 * t324 + t20 * t325 + t33 * t326 + t34 * t327 + t97 * t332 + t96 * t333 + t209 * t334 + t23 * t315 + t270 + t24 * t347 + (m(6) * t244 - t276 * t76 - t278 * t75 + t253) * pkin(10) + m(7) * (t113 * t13 + t114 * t12 + t157 * t3 + t158 * t2 + t16 * t225 + t273 * t43) + (t272 + (pkin(5) * t40 - t50 / 0.2e1) * t237) * qJD(5) + t271 * t147 + t262 * t117 + ((-t237 * t31 - t241 * t30) * qJD(5) + t254) * mrSges(6,3) + t225 * t9 + t29 * t207 + t57 * t194 + t123 * t197 / 0.2e1 + t124 * t199 / 0.2e1 + t157 * t14 + t158 * t15 + t16 * t149 + t114 * t48 + t113 * t47 + t43 * t94 - t35 * mrSges(5,2) + t36 * mrSges(5,1) + (-t330 - t26) * pkin(4); t229 + t225 * t67 + t62 * t320 + t63 * t319 + t96 * t322 + t97 * t321 + t182 * t94 + t153 * t149 + t157 * t89 + t158 * t90 + t113 * t159 + t114 * t160 + t130 * t327 + t129 * t326 + t106 * t325 + t105 * t324 - pkin(4) * t139 + m(7) * (t113 * t78 + t114 * t77 + t153 * t225 + t157 * t39 + t158 * t38) + ((-m(6) * pkin(4) + t294) * t282 - t271) * t242 + (-t209 * t279 / 0.2e1 + t328 - t104 * mrSges(6,3) + (-t144 * mrSges(6,3) - t174 / 0.2e1 + (m(7) * t182 + t135) * pkin(5)) * qJD(5) + (m(6) * (-t104 - t348) - qJD(5) * t201 - t165) * pkin(10)) * t237 + (t145 * t77 - t146 * t78 - t192 * t39 - t252 * t38) * mrSges(7,3) + (t279 * t318 + t329 + qJD(5) * t323 + t259 * mrSges(6,3) + (m(6) * t259 - qJD(5) * t202 + t166) * pkin(10)) * t241 + (t199 * t315 + t197 * t317 + t223 * t194 + (t209 * t316 + t211 * t317) * qJD(5) + (t223 * mrSges(5,2) - Ifges(5,6) + t262) * qJD(4)) * t238; m(7) * (-pkin(5) * t265 + t105 * t158 + t106 * t157 - t113 * t179 - t114 * t178) + (-t105 * t252 - t106 * t192 - t145 * t178 + t146 * t179) * mrSges(7,3) + (m(6) * (t284 * t312 - t313) + (m(7) * t225 + t149 + t294) * t238) * qJD(4) + (-t194 - t94 + (mrSges(6,3) * t284 - mrSges(5,2)) * qJD(4)) * t242; (t113 * t158 + t114 * t157 + t225 * t273) * t341 + 0.2e1 * t149 * t273 + 0.2e1 * t225 * t94 - t145 * t152 + t192 * t97 - t146 * t151 - t252 * t96 - 0.2e1 * pkin(4) * t194 + t237 * t199 - t209 * t278 + (qJD(5) * t211 + t197) * t241 + 0.2e1 * (-t113 * t252 - t114 * t192 + t145 * t157 - t146 * t158) * mrSges(7,3); t11 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (-t12 * t275 + t13 * t274 + t2 * t236 + t240 * t3) + t47 * t274 + t236 * t15 - t48 * t275 + t240 * t14) * pkin(5) + t250 + t22; t104 * mrSges(6,1) - t103 * mrSges(6,2) + (m(7) * (t236 * t38 + t240 * t39 + t274 * t78 - t275 * t77) + t159 * t274 + t236 * t90 - t160 * t275 + t240 * t89) * pkin(5) + t132 + t248; m(7) * (t105 * t236 + t106 * t240 + (t178 * t236 - t179 * t240) * qJD(6)) * pkin(5) + t346; t228 + (pkin(10) * t207 - t301) * qJD(5) + (m(7) * (t113 * t236 + t114 * t240 + (-t157 * t236 + t158 * t240) * qJD(6)) + (t240 * t145 - t236 * t146 + (t192 * t236 - t240 * t252) * qJD(6)) * mrSges(7,3)) * pkin(5) + t251; 0.2e1 * t189; t250; t248; -t67; t251; t189; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
