% Calculate time derivative of joint inertia matrix for
% S6RRPRRR4
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:59
% EndTime: 2019-03-09 13:28:17
% DurationCPUTime: 8.28s
% Computational Cost: add. (15251->589), mult. (40690->886), div. (0->0), fcn. (41525->12), ass. (0->259)
t233 = sin(qJ(5));
t234 = sin(qJ(4));
t237 = cos(qJ(5));
t238 = cos(qJ(4));
t196 = t233 * t234 - t237 * t238;
t357 = qJD(4) + qJD(5);
t158 = t357 * t196;
t232 = sin(qJ(6));
t236 = cos(qJ(6));
t287 = t232 ^ 2 + t236 ^ 2;
t265 = t287 * t158;
t358 = t287 * t237;
t228 = sin(pkin(12));
t229 = sin(pkin(6));
t230 = cos(pkin(12));
t239 = cos(qJ(2));
t294 = t230 * t239;
t235 = sin(qJ(2));
t296 = t229 * t235;
t183 = t228 * t296 - t229 * t294;
t184 = (t228 * t239 + t230 * t235) * t229;
t231 = cos(pkin(6));
t161 = t184 * t238 + t231 * t234;
t332 = pkin(1) * t231;
t218 = t239 * t332;
t328 = -pkin(8) - qJ(3);
t270 = t328 * t235;
t169 = pkin(2) * t231 + t229 * t270 + t218;
t217 = t235 * t332;
t295 = t229 * t239;
t193 = pkin(8) * t295 + t217;
t179 = qJ(3) * t295 + t193;
t122 = t228 * t169 + t230 * t179;
t112 = pkin(9) * t231 + t122;
t198 = (-pkin(2) * t239 - pkin(1)) * t229;
t128 = t183 * pkin(3) - t184 * pkin(9) + t198;
t77 = -t112 * t234 + t238 * t128;
t61 = pkin(4) * t183 - pkin(10) * t161 + t77;
t160 = -t184 * t234 + t231 * t238;
t78 = t238 * t112 + t234 * t128;
t67 = pkin(10) * t160 + t78;
t324 = t233 * t61 + t237 * t67;
t26 = pkin(11) * t183 + t324;
t106 = t160 * t233 + t161 * t237;
t251 = t237 * t160 - t161 * t233;
t121 = t169 * t230 - t228 * t179;
t111 = -pkin(3) * t231 - t121;
t83 = -pkin(4) * t160 + t111;
t49 = -pkin(5) * t251 - pkin(11) * t106 + t83;
t16 = -t232 * t26 + t236 * t49;
t286 = qJD(2) * t229;
t181 = (-t228 * t235 + t294) * t286;
t119 = -qJD(4) * t161 - t181 * t234;
t120 = qJD(4) * t160 + t181 * t238;
t56 = qJD(5) * t251 + t119 * t233 + t120 * t237;
t57 = qJD(5) * t106 - t237 * t119 + t120 * t233;
t215 = qJD(2) * t218;
t144 = t215 + (qJD(2) * t270 + qJD(3) * t239) * t229;
t359 = (t295 * t328 - t217) * qJD(2) - qJD(3) * t296;
t102 = t144 * t228 - t230 * t359;
t76 = -pkin(4) * t119 + t102;
t18 = pkin(5) * t57 - pkin(11) * t56 + t76;
t180 = qJD(2) * t184;
t103 = t230 * t144 + t359 * t228;
t275 = t235 * t286;
t264 = pkin(2) * t275;
t125 = pkin(3) * t180 - pkin(9) * t181 + t264;
t40 = -t78 * qJD(4) - t103 * t234 + t238 * t125;
t28 = pkin(4) * t180 - pkin(10) * t120 + t40;
t283 = qJD(5) * t237;
t310 = t233 * t67;
t284 = qJD(4) * t238;
t285 = qJD(4) * t234;
t39 = t238 * t103 - t112 * t285 + t234 * t125 + t128 * t284;
t34 = pkin(10) * t119 + t39;
t8 = -qJD(5) * t310 + t233 * t28 + t237 * t34 + t61 * t283;
t5 = pkin(11) * t180 + t8;
t2 = qJD(6) * t16 + t18 * t232 + t236 * t5;
t17 = t232 * t49 + t236 * t26;
t3 = -qJD(6) * t17 + t18 * t236 - t232 * t5;
t365 = t2 * t236 - t232 * t3;
t364 = 0.2e1 * t229;
t363 = (-t16 * t236 - t17 * t232) * qJD(6) + t365;
t197 = t233 * t238 + t234 * t237;
t281 = qJD(6) * t232;
t290 = t236 * t158;
t248 = t197 * t281 + t290;
t280 = qJD(6) * t236;
t362 = -t232 * t158 + t197 * t280;
t219 = pkin(2) * t228 + pkin(9);
t327 = pkin(10) + t219;
t194 = t327 * t238;
t269 = t327 * t234;
t141 = t233 * t194 + t237 * t269;
t266 = qJD(4) * t327;
t189 = t234 * t266;
t263 = t238 * t266;
t100 = -qJD(5) * t141 - t237 * t189 - t233 * t263;
t278 = pkin(4) * t285;
t159 = t357 * t197;
t331 = pkin(5) * t159;
t104 = pkin(11) * t158 + t278 + t331;
t220 = -pkin(2) * t230 - pkin(3);
t206 = -pkin(4) * t238 + t220;
t139 = pkin(5) * t196 - pkin(11) * t197 + t206;
t142 = t237 * t194 - t233 * t269;
t95 = t139 * t236 - t142 * t232;
t47 = qJD(6) * t95 + t100 * t236 + t104 * t232;
t96 = t139 * t232 + t142 * t236;
t48 = -qJD(6) * t96 - t100 * t232 + t104 * t236;
t360 = -t48 * t232 + t236 * t47;
t334 = t232 / 0.2e1;
t333 = t236 / 0.2e1;
t268 = -t281 / 0.2e1;
t200 = (mrSges(5,1) * t234 + mrSges(5,2) * t238) * qJD(4);
t356 = Ifges(5,5) * t120 + Ifges(5,6) * t119 + Ifges(5,3) * t180;
t9 = -qJD(5) * t324 - t233 * t34 + t237 * t28;
t355 = 2 * m(6);
t354 = 2 * m(7);
t353 = 0.2e1 * pkin(4);
t352 = -2 * mrSges(3,3);
t351 = -2 * mrSges(4,3);
t350 = -2 * mrSges(6,3);
t349 = -2 * Ifges(4,4);
t101 = qJD(5) * t142 - t233 * t189 + t237 * t263;
t348 = 0.2e1 * t101;
t347 = 0.2e1 * t184;
t346 = 0.2e1 * t198;
t345 = m(6) / 0.2e1;
t344 = m(4) * pkin(2);
t85 = t106 * t236 + t183 * t232;
t38 = -qJD(6) * t85 + t180 * t236 - t232 * t56;
t343 = t38 / 0.2e1;
t84 = -t106 * t232 + t183 * t236;
t342 = t84 / 0.2e1;
t224 = Ifges(7,5) * t280;
t340 = Ifges(7,6) * t268 + t224 / 0.2e1;
t321 = Ifges(7,4) * t232;
t258 = Ifges(7,1) * t236 - t321;
t204 = t258 * qJD(6);
t339 = t204 / 0.2e1;
t338 = Ifges(7,5) * t334 + Ifges(7,6) * t333;
t320 = Ifges(7,4) * t236;
t211 = Ifges(7,1) * t232 + t320;
t336 = t211 / 0.2e1;
t335 = -t232 / 0.2e1;
t37 = qJD(6) * t84 + t180 * t232 + t236 * t56;
t15 = -mrSges(7,1) * t38 + mrSges(7,2) * t37;
t50 = mrSges(6,1) * t180 - mrSges(6,3) * t56;
t326 = t15 - t50;
t58 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t87 = mrSges(6,1) * t183 - mrSges(6,3) * t106;
t325 = t58 - t87;
t323 = Ifges(5,4) * t234;
t322 = Ifges(5,4) * t238;
t319 = Ifges(7,6) * t232;
t318 = pkin(4) * qJD(5);
t317 = t103 * mrSges(4,2);
t316 = t180 * Ifges(6,5);
t315 = t180 * Ifges(6,6);
t314 = t183 * Ifges(5,6);
t186 = -pkin(8) * t275 + t215;
t313 = t186 * mrSges(3,2);
t187 = t193 * qJD(2);
t312 = t187 * mrSges(3,1);
t311 = t233 * mrSges(6,1);
t308 = t237 * mrSges(6,2);
t306 = t101 * t141;
t305 = t141 * t233;
t302 = t159 * t196;
t301 = t196 * t233;
t300 = t197 * t232;
t299 = t197 * t236;
t292 = t232 * t237;
t207 = -mrSges(7,1) * t236 + mrSges(7,2) * t232;
t291 = t233 * t207;
t289 = t236 * t237;
t288 = -Ifges(6,5) * t158 - Ifges(6,6) * t159;
t282 = qJD(6) * t197;
t12 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t57;
t279 = Ifges(6,5) * t56 - Ifges(6,6) * t57 + Ifges(6,3) * t180;
t6 = -pkin(5) * t180 - t9;
t277 = m(7) * t6 + t15;
t276 = Ifges(3,5) * t239 * t286 + Ifges(4,5) * t181 - Ifges(4,6) * t180;
t82 = -mrSges(7,1) * t362 + t248 * mrSges(7,2);
t271 = m(7) * t101 - t82;
t267 = t280 / 0.2e1;
t108 = t159 * mrSges(6,1) - t158 * mrSges(6,2);
t262 = mrSges(7,3) * t358;
t261 = -mrSges(5,1) * t238 + mrSges(5,2) * t234;
t259 = mrSges(7,1) * t232 + mrSges(7,2) * t236;
t257 = -Ifges(7,2) * t232 + t320;
t19 = mrSges(7,1) * t57 - mrSges(7,3) * t37;
t20 = -mrSges(7,2) * t57 + mrSges(7,3) * t38;
t256 = -t232 * t19 + t236 * t20;
t29 = t237 * t61 - t310;
t253 = -t40 * t234 + t39 * t238;
t252 = t101 * t196 + t141 * t159;
t62 = mrSges(7,2) * t251 + mrSges(7,3) * t84;
t63 = -mrSges(7,1) * t251 - mrSges(7,3) * t85;
t86 = -mrSges(6,2) * t183 + mrSges(6,3) * t251;
t250 = -t232 * t63 + t236 * t62 + t86;
t202 = t257 * qJD(6);
t209 = Ifges(7,2) * t236 + t321;
t247 = t236 * t202 + t232 * t204 - t209 * t281 + t211 * t280;
t199 = t259 * qJD(6);
t246 = -mrSges(7,3) * t265 + t159 * t207 + t196 * t199 - t108;
t72 = -Ifges(7,5) * t248 - Ifges(7,6) * t362 + Ifges(7,3) * t159;
t243 = -t63 * t280 - t62 * t281 + m(7) * (-t16 * t280 - t17 * t281 + t365) + t256;
t151 = -mrSges(7,2) * t196 - mrSges(7,3) * t300;
t152 = mrSges(7,1) * t196 - mrSges(7,3) * t299;
t89 = mrSges(7,1) * t159 + mrSges(7,3) * t248;
t90 = -mrSges(7,2) * t159 - mrSges(7,3) * t362;
t242 = -t152 * t280 - t151 * t281 + m(7) * (-t280 * t95 - t281 * t96 + t360) + t236 * t90 - t232 * t89;
t13 = Ifges(7,4) * t37 + Ifges(7,2) * t38 + Ifges(7,6) * t57;
t14 = Ifges(7,1) * t37 + Ifges(7,4) * t38 + Ifges(7,5) * t57;
t25 = -pkin(5) * t183 - t29;
t44 = Ifges(7,4) * t85 + Ifges(7,2) * t84 - Ifges(7,6) * t251;
t45 = Ifges(7,1) * t85 + Ifges(7,4) * t84 - Ifges(7,5) * t251;
t241 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + mrSges(7,3) * t363 + t13 * t333 + t14 * t334 + t25 * t199 + t202 * t342 + t6 * t207 + t209 * t343 - t251 * t340 + t45 * t267 + t44 * t268 + t37 * t336 + t57 * t338 + t85 * t339 + t279;
t132 = Ifges(7,6) * t196 + t197 * t257;
t133 = Ifges(7,5) * t196 + t197 * t258;
t73 = -Ifges(7,4) * t248 - Ifges(7,2) * t362 + Ifges(7,6) * t159;
t74 = -Ifges(7,1) * t248 - Ifges(7,4) * t362 + Ifges(7,5) * t159;
t240 = t133 * t267 + t141 * t199 - t290 * t336 + t159 * t338 - t202 * t300 / 0.2e1 + t299 * t339 + t196 * t340 + t74 * t334 + t73 * t333 - t100 * mrSges(6,2) + t288 - t362 * t209 / 0.2e1 + (t197 * t211 + t132) * t268 + (t207 - mrSges(6,1)) * t101 + ((-t232 * t96 - t236 * t95) * qJD(6) + t360) * mrSges(7,3);
t225 = Ifges(5,5) * t284;
t223 = -pkin(4) * t237 - pkin(5);
t222 = pkin(4) * t233 + pkin(11);
t212 = Ifges(5,1) * t234 + t322;
t210 = Ifges(5,2) * t238 + t323;
t205 = (Ifges(5,1) * t238 - t323) * qJD(4);
t203 = (-Ifges(5,2) * t234 + t322) * qJD(4);
t192 = -pkin(8) * t296 + t218;
t174 = t181 * mrSges(4,2);
t166 = Ifges(6,1) * t197 - Ifges(6,4) * t196;
t165 = Ifges(6,4) * t197 - Ifges(6,2) * t196;
t164 = mrSges(6,1) * t196 + mrSges(6,2) * t197;
t143 = t259 * t197;
t131 = Ifges(7,3) * t196 + (Ifges(7,5) * t236 - t319) * t197;
t130 = mrSges(5,1) * t183 - mrSges(5,3) * t161;
t129 = -mrSges(5,2) * t183 + mrSges(5,3) * t160;
t110 = -Ifges(6,1) * t158 - Ifges(6,4) * t159;
t109 = -Ifges(6,4) * t158 - Ifges(6,2) * t159;
t94 = mrSges(5,1) * t180 - mrSges(5,3) * t120;
t93 = -mrSges(5,2) * t180 + mrSges(5,3) * t119;
t92 = Ifges(5,1) * t161 + Ifges(5,4) * t160 + Ifges(5,5) * t183;
t91 = Ifges(5,4) * t161 + Ifges(5,2) * t160 + t314;
t79 = -mrSges(5,1) * t119 + mrSges(5,2) * t120;
t75 = -mrSges(6,1) * t251 + mrSges(6,2) * t106;
t69 = Ifges(5,1) * t120 + Ifges(5,4) * t119 + t180 * Ifges(5,5);
t68 = Ifges(5,4) * t120 + Ifges(5,2) * t119 + t180 * Ifges(5,6);
t65 = Ifges(6,1) * t106 + Ifges(6,4) * t251 + Ifges(6,5) * t183;
t64 = Ifges(6,4) * t106 + Ifges(6,2) * t251 + Ifges(6,6) * t183;
t51 = -mrSges(6,2) * t180 - mrSges(6,3) * t57;
t43 = Ifges(7,5) * t85 + Ifges(7,6) * t84 - Ifges(7,3) * t251;
t23 = mrSges(6,1) * t57 + mrSges(6,2) * t56;
t22 = Ifges(6,1) * t56 - Ifges(6,4) * t57 + t316;
t21 = Ifges(6,4) * t56 - Ifges(6,2) * t57 + t315;
t1 = [0.2e1 * t324 * t51 + (t29 * t9 + t324 * t8 + t76 * t83) * t355 + (mrSges(4,1) * t346 + t122 * t351 + t184 * t349 + Ifges(5,5) * t161 + Ifges(6,5) * t106 - Ifges(4,6) * t231 + Ifges(5,6) * t160 + Ifges(6,6) * t251 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t183) * t180 - (t12 - t21) * t251 + (Ifges(4,1) * t347 + Ifges(4,5) * t231 + t121 * t351 + t183 * t349) * t181 + t174 * t346 + t160 * t68 + 0.2e1 * t102 * (-mrSges(5,1) * t160 + mrSges(5,2) * t161) + t161 * t69 + 0.2e1 * t39 * t129 + 0.2e1 * t40 * t130 + 0.2e1 * t111 * t79 + t119 * t91 + t120 * t92 + t106 * t22 + 0.2e1 * t78 * t93 + 0.2e1 * t77 * t94 + 0.2e1 * t83 * t23 + t84 * t13 + t85 * t14 + 0.2e1 * t8 * t86 + 0.2e1 * t9 * t87 + 0.2e1 * t76 * t75 + 0.2e1 * t6 * t58 + 0.2e1 * t2 * t62 + 0.2e1 * t3 * t63 + t56 * t65 + 0.2e1 * t29 * t50 + t38 * t44 + t37 * t45 + 0.2e1 * t25 * t15 + 0.2e1 * t16 * t19 + 0.2e1 * t17 * t20 + 0.2e1 * m(3) * (t186 * t193 - t187 * t192) + (-0.2e1 * t102 * mrSges(4,1) + t276 - 0.2e1 * t312 - 0.2e1 * t313 - 0.2e1 * t317) * t231 + (t103 * t351 + t279 + t356) * t183 + (-t64 + t43) * t57 + 0.2e1 * m(4) * (-t102 * t121 + t103 * t122) + (t16 * t3 + t17 * t2 + t25 * t6) * t354 + t102 * mrSges(4,3) * t347 + 0.2e1 * m(5) * (t102 * t111 + t39 * t78 + t40 * t77) + (0.2e1 * (t186 * t239 + t187 * t235) * mrSges(3,3) + ((t192 * t352 + Ifges(3,5) * t231 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t239) * t364) * t239 + (-0.2e1 * Ifges(3,6) * t231 + t344 * t346 + 0.2e1 * pkin(2) * (mrSges(4,1) * t183 + mrSges(4,2) * t184) + t193 * t352 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t235 + (Ifges(3,1) - Ifges(3,2)) * t239) * t364) * t235) * qJD(2)) * t229; t276 + (t288 + t225) * t183 / 0.2e1 + (t43 / 0.2e1 - t64 / 0.2e1 - t324 * mrSges(6,3)) * t159 + m(6) * (t100 * t324 - t101 * t29 - t141 * t9 + t142 * t8 + t206 * t76 + t278 * t83) - (t72 / 0.2e1 - t109 / 0.2e1) * t251 - t312 - t313 + (-t165 / 0.2e1 + t131 / 0.2e1) * t57 - (t65 / 0.2e1 + t45 * t333 + t44 * t335 - t29 * mrSges(6,3)) * t158 + (t316 / 0.2e1 + t22 / 0.2e1 + t14 * t333 + t13 * t335 - t9 * mrSges(6,3) + (-t236 * t44 / 0.2e1 + t45 * t335) * qJD(6)) * t197 + t325 * t101 + t326 * t141 + ((-t219 * t130 - t77 * mrSges(5,3) + t92 / 0.2e1) * t238 + (-t219 * t129 + pkin(4) * t75 - t78 * mrSges(5,3) - t314 / 0.2e1 - t91 / 0.2e1) * t234) * qJD(4) + (-t21 / 0.2e1 - t315 / 0.2e1 + t12 / 0.2e1 - t8 * mrSges(6,3)) * t196 - Ifges(3,6) * t275 + (-mrSges(4,1) + t261) * t102 + t253 * mrSges(5,3) + t73 * t342 + t132 * t343 + (-t102 * t230 + t103 * t228) * t344 + t180 * (Ifges(5,5) * t234 + Ifges(5,6) * t238) / 0.2e1 + t238 * t68 / 0.2e1 + t234 * t69 / 0.2e1 + t206 * t23 + t119 * t210 / 0.2e1 + t120 * t212 / 0.2e1 + t220 * t79 + t111 * t200 + t160 * t203 / 0.2e1 + t161 * t205 / 0.2e1 + t56 * t166 / 0.2e1 + t76 * t164 + t2 * t151 + t3 * t152 + t142 * t51 + t6 * t143 + t37 * t133 / 0.2e1 + t83 * t108 + t106 * t110 / 0.2e1 + t17 * t90 + t95 * t19 + t96 * t20 + t100 * t86 + t85 * t74 / 0.2e1 + t16 * t89 - t25 * t82 + t47 * t62 + t48 * t63 + (-t180 * t228 - t181 * t230) * pkin(2) * mrSges(4,3) - t317 + (t102 * t220 + ((-t234 * t78 - t238 * t77) * qJD(4) + t253) * t219) * m(5) + m(7) * (t101 * t25 + t141 * t6 + t16 * t48 + t17 * t47 + t2 * t96 + t3 * t95) + (-t234 * t94 + t238 * t93) * t219; t143 * t348 + 0.2e1 * t206 * t108 - 0.2e1 * t141 * t82 + 0.2e1 * t47 * t151 + 0.2e1 * t48 * t152 + 0.2e1 * t220 * t200 + t238 * t203 + t234 * t205 + 0.2e1 * t95 * t89 + 0.2e1 * t96 * t90 + (t238 * t212 + (t164 * t353 - t210) * t234) * qJD(4) + (t100 * t142 + t206 * t278 + t306) * t355 + (t47 * t96 + t48 * t95 + t306) * t354 + (t100 * t350 - t109 + t72) * t196 + (t142 * t350 + t131 - t165) * t159 - (0.2e1 * mrSges(6,3) * t141 - t132 * t232 + t133 * t236 + t166) * t158 + (mrSges(6,3) * t348 - t232 * t73 + t236 * t74 + t110 + (-t132 * t236 - t133 * t232) * qJD(6)) * t197; m(4) * t264 + t180 * mrSges(4,1) + t234 * t93 + t238 * t94 + t174 + t326 * t196 + t325 * t159 + (t129 * t238 - t130 * t234) * qJD(4) - t250 * t158 + (t51 + (-t232 * t62 - t236 * t63) * qJD(6) + t256) * t197 + m(7) * (t159 * t25 + t196 * t6 - (-t16 * t232 + t17 * t236) * t158 + t363 * t197) + m(6) * (-t158 * t324 - t159 * t29 - t196 * t9 + t197 * t8) + m(5) * (t234 * t39 + t238 * t40 + (-t234 * t77 + t238 * t78) * qJD(4)); t159 * t143 - t196 * t82 + m(7) * t252 + m(6) * (t100 * t197 - t142 * t158 + t252) + (m(7) * (-t158 * t96 + t197 * t47 - t282 * t95) - t158 * t151 + t197 * t90 - t152 * t282) * t236 + (m(7) * (t158 * t95 - t197 * t48 - t282 * t96) - t151 * t282 + t158 * t152 - t197 * t89) * t232; 0.2e1 * m(6) * (-t158 * t197 + t302) + 0.2e1 * m(7) * (-t197 * t265 + t302); t40 * mrSges(5,1) - t39 * mrSges(5,2) + t241 + t243 * t222 + t277 * t223 + (m(6) * (t233 * t8 + t237 * t9) + t237 * t50 + t233 * t51 + (t325 * t233 + t250 * t237 + m(7) * (-t16 * t292 + t17 * t289 + t233 * t25) + m(6) * (-t233 * t29 + t237 * t324)) * qJD(5)) * pkin(4) + t356; t225 + t240 + (-Ifges(5,6) * t234 + t219 * t261) * qJD(4) + (m(6) * (t100 * t233 - t101 * t237) + (t158 * t237 - t159 * t233) * mrSges(6,3) + ((t197 * mrSges(6,3) + t143) * t233 + (-t196 * mrSges(6,3) + t236 * t151 - t232 * t152) * t237 + m(7) * (t289 * t96 - t292 * t95 + t305) + m(6) * (t142 * t237 + t305)) * qJD(5)) * pkin(4) + t271 * t223 + t242 * t222; m(7) * (t159 * t223 - t222 * t265) - t200 + ((-t158 * t233 - t159 * t237) * t345 + ((t197 * t237 + t301) * t345 + m(7) * (t358 * t197 + t301) / 0.2e1) * qJD(5)) * t353 + t246; 0.2e1 * t223 * t199 + (0.2e1 * t291 + (t358 * t222 + t223 * t233) * t354 - 0.2e1 * t308 - 0.2e1 * t311 + 0.2e1 * t262) * t318 + t247; -pkin(5) * t277 + pkin(11) * t243 + t241; -pkin(5) * t271 + pkin(11) * t242 + t240; m(7) * (-pkin(11) * t265 - t331) + t246; (-pkin(5) + t223) * t199 + (m(7) * (-pkin(5) * t233 + t358 * pkin(11)) + t291 - t308 - t311 + t262) * t318 + t247; -0.2e1 * pkin(5) * t199 + t247; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t48 - mrSges(7,2) * t47 + t72; t82; t224 - t259 * pkin(4) * t283 + (t207 * t222 - t319) * qJD(6); t224 + (pkin(11) * t207 - t319) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
