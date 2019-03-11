% Calculate time derivative of joint inertia matrix for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:57
% EndTime: 2019-03-09 16:54:17
% DurationCPUTime: 9.00s
% Computational Cost: add. (9759->683), mult. (25747->974), div. (0->0), fcn. (24951->10), ass. (0->283)
t377 = Ifges(6,5) + Ifges(7,5);
t376 = Ifges(6,6) + Ifges(7,6);
t375 = Ifges(4,3) + Ifges(5,3);
t247 = cos(qJ(5));
t244 = sin(qJ(5));
t323 = Ifges(7,4) * t244;
t325 = Ifges(6,4) * t244;
t364 = t323 + t325 + (Ifges(6,2) + Ifges(7,2)) * t247;
t322 = Ifges(7,4) * t247;
t324 = Ifges(6,4) * t247;
t363 = t322 + t324 + (Ifges(6,1) + Ifges(7,1)) * t244;
t241 = sin(pkin(11));
t248 = cos(qJ(3));
t245 = sin(qJ(3));
t316 = cos(pkin(11));
t268 = t316 * t245;
t200 = t241 * t248 + t268;
t296 = qJD(5) * t244;
t276 = t200 * t296;
t309 = t241 * t245;
t251 = t248 * t316 - t309;
t192 = t251 * qJD(3);
t306 = t247 * t192;
t252 = t276 - t306;
t295 = qJD(5) * t247;
t279 = t200 * t295;
t314 = t192 * t244;
t253 = t279 + t314;
t243 = cos(pkin(6));
t242 = sin(pkin(6));
t246 = sin(qJ(2));
t308 = t242 * t246;
t194 = t243 * t245 + t248 * t308;
t249 = cos(qJ(2));
t300 = qJD(2) * t242;
t280 = t249 * t300;
t157 = -qJD(3) * t194 - t245 * t280;
t193 = t243 * t248 - t245 * t308;
t158 = qJD(3) * t193 + t248 * t280;
t102 = -t157 * t316 + t158 * t241;
t103 = t241 * t157 + t158 * t316;
t136 = t241 * t193 + t194 * t316;
t307 = t242 * t249;
t110 = -t244 * t136 - t247 * t307;
t299 = qJD(2) * t246;
t281 = t242 * t299;
t59 = qJD(5) * t110 + t247 * t103 + t244 * t281;
t254 = -t247 * t136 + t244 * t307;
t60 = qJD(5) * t254 - t244 * t103 + t247 * t281;
t7 = Ifges(7,5) * t59 + Ifges(7,6) * t60 + Ifges(7,3) * t102;
t8 = Ifges(6,5) * t59 + Ifges(6,6) * t60 + Ifges(6,3) * t102;
t374 = t7 + t8;
t337 = t244 / 0.2e1;
t336 = t247 / 0.2e1;
t271 = -t296 / 0.2e1;
t10 = Ifges(6,4) * t59 + Ifges(6,2) * t60 + Ifges(6,6) * t102;
t9 = Ifges(7,4) * t59 + Ifges(7,2) * t60 + Ifges(7,6) * t102;
t373 = t9 + t10;
t11 = Ifges(7,1) * t59 + Ifges(7,4) * t60 + Ifges(7,5) * t102;
t12 = Ifges(6,1) * t59 + Ifges(6,4) * t60 + Ifges(6,5) * t102;
t372 = t11 + t12;
t135 = -t193 * t316 + t194 * t241;
t39 = -Ifges(7,4) * t254 + Ifges(7,2) * t110 + Ifges(7,6) * t135;
t40 = -Ifges(6,4) * t254 + Ifges(6,2) * t110 + Ifges(6,6) * t135;
t371 = t39 + t40;
t41 = -Ifges(7,1) * t254 + Ifges(7,4) * t110 + Ifges(7,5) * t135;
t42 = -Ifges(6,1) * t254 + Ifges(6,4) * t110 + Ifges(6,5) * t135;
t370 = t41 + t42;
t191 = t200 * qJD(3);
t70 = -Ifges(7,4) * t252 - Ifges(7,2) * t253 + Ifges(7,6) * t191;
t71 = -Ifges(6,4) * t252 - Ifges(6,2) * t253 + Ifges(6,6) * t191;
t369 = t71 + t70;
t72 = -Ifges(7,1) * t252 - Ifges(7,4) * t253 + Ifges(7,5) * t191;
t73 = -Ifges(6,1) * t252 - Ifges(6,4) * t253 + Ifges(6,5) * t191;
t368 = t72 + t73;
t257 = -Ifges(7,2) * t244 + t322;
t114 = -Ifges(7,6) * t251 + t200 * t257;
t258 = -Ifges(6,2) * t244 + t324;
t115 = -Ifges(6,6) * t251 + t200 * t258;
t367 = t115 + t114;
t259 = Ifges(7,1) * t247 - t323;
t116 = -Ifges(7,5) * t251 + t200 * t259;
t260 = Ifges(6,1) * t247 - t325;
t117 = -Ifges(6,5) * t251 + t200 * t260;
t303 = t116 + t117;
t235 = -pkin(3) * t248 - pkin(2);
t145 = -pkin(4) * t251 - pkin(10) * t200 + t235;
t330 = -qJ(4) - pkin(9);
t216 = t330 * t248;
t160 = -t216 * t316 + t309 * t330;
t156 = t247 * t160;
t94 = t244 * t145 + t156;
t196 = t243 * t246 * pkin(1) + pkin(8) * t307;
t178 = pkin(9) * t243 + t196;
t179 = (-pkin(2) * t249 - pkin(9) * t246 - pkin(1)) * t242;
t122 = t248 * t178 + t245 * t179;
t366 = (t257 + t258) * qJD(5);
t365 = (t259 + t260) * qJD(5);
t228 = pkin(8) * t308;
t335 = pkin(1) * t249;
t195 = t243 * t335 - t228;
t269 = qJD(3) * t330;
t190 = qJD(4) * t248 + t245 * t269;
t250 = -qJD(4) * t245 + t248 * t269;
t130 = t190 * t316 + t241 * t250;
t298 = qJD(3) * t245;
t292 = pkin(3) * t298;
t131 = pkin(4) * t191 - pkin(10) * t192 + t292;
t283 = t247 * t130 + t244 * t131 + t145 * t295;
t33 = -t160 * t296 + t283;
t267 = -t130 * t244 + t247 * t131;
t34 = -qJD(5) * t94 + t267;
t360 = -t244 * t34 + t247 * t33;
t182 = (pkin(2) * t246 - pkin(9) * t249) * t300;
t183 = t195 * qJD(2);
t79 = -t122 * qJD(3) + t248 * t182 - t183 * t245;
t45 = pkin(3) * t281 - qJ(4) * t158 - qJD(4) * t194 + t79;
t297 = qJD(3) * t248;
t78 = -t178 * t298 + t179 * t297 + t245 * t182 + t248 * t183;
t52 = qJ(4) * t157 + qJD(4) * t193 + t78;
t18 = t241 * t45 + t316 * t52;
t16 = pkin(10) * t281 + t18;
t184 = t196 * qJD(2);
t120 = -t157 * pkin(3) + t184;
t32 = t102 * pkin(4) - t103 * pkin(10) + t120;
t104 = qJ(4) * t193 + t122;
t121 = -t245 * t178 + t248 * t179;
t90 = -pkin(3) * t307 - t194 * qJ(4) + t121;
t51 = t316 * t104 + t241 * t90;
t44 = -pkin(10) * t307 + t51;
t177 = t228 + (-pkin(2) - t335) * t243;
t139 = -t193 * pkin(3) + t177;
t64 = t135 * pkin(4) - t136 * pkin(10) + t139;
t3 = t247 * t16 + t244 * t32 + t64 * t295 - t296 * t44;
t20 = t244 * t64 + t247 * t44;
t4 = -qJD(5) * t20 - t16 * t244 + t247 * t32;
t359 = -t244 * t4 + t247 * t3;
t358 = t376 * t336 + t337 * t377;
t238 = Ifges(7,5) * t295;
t239 = Ifges(6,5) * t295;
t357 = t239 / 0.2e1 + t238 / 0.2e1 + t376 * t271;
t356 = 2 * m(5);
t355 = 2 * m(6);
t354 = 2 * m(7);
t353 = -2 * mrSges(3,3);
t352 = -2 * mrSges(5,3);
t351 = -2 * mrSges(7,3);
t129 = t190 * t241 - t316 * t250;
t350 = 0.2e1 * t129;
t159 = -t216 * t241 - t330 * t268;
t349 = 0.2e1 * t159;
t348 = 0.2e1 * t184;
t347 = m(5) * pkin(3);
t346 = m(7) * pkin(5);
t334 = pkin(3) * t241;
t329 = mrSges(6,2) * t247;
t328 = mrSges(7,3) * t247;
t327 = Ifges(4,4) * t245;
t326 = Ifges(4,4) * t248;
t321 = t183 * mrSges(3,2);
t320 = t184 * mrSges(3,1);
t319 = t244 * mrSges(7,3);
t315 = t129 * t159;
t313 = t200 * t244;
t312 = t200 * t247;
t233 = pkin(10) + t334;
t311 = t233 * t244;
t310 = t233 * t247;
t305 = qJ(6) + t233;
t302 = Ifges(7,5) * t306 + Ifges(7,3) * t191;
t301 = Ifges(6,5) * t306 + Ifges(6,3) * t191;
t202 = mrSges(7,1) * t296 + mrSges(7,2) * t295;
t294 = 0.2e1 * t242;
t291 = pkin(5) * t296;
t290 = Ifges(4,6) * t307;
t289 = -t39 / 0.2e1 - t40 / 0.2e1;
t288 = t41 / 0.2e1 + t42 / 0.2e1;
t287 = mrSges(6,3) * t296;
t286 = mrSges(6,3) * t295;
t285 = mrSges(7,3) * t296;
t284 = mrSges(7,3) * t295;
t282 = t316 * pkin(3);
t278 = t233 * t296;
t277 = t233 * t295;
t21 = -t60 * mrSges(7,1) + t59 * mrSges(7,2);
t270 = t295 / 0.2e1;
t58 = t102 * mrSges(5,1) + t103 * mrSges(5,2);
t19 = -t244 * t44 + t247 * t64;
t140 = t191 * mrSges(5,1) + t192 * mrSges(5,2);
t93 = t247 * t145 - t160 * t244;
t266 = qJD(5) * t305;
t265 = Ifges(5,5) * t281;
t264 = Ifges(5,6) * t281;
t234 = -t282 - pkin(4);
t17 = -t241 * t52 + t316 * t45;
t50 = -t241 * t104 + t316 * t90;
t261 = mrSges(6,1) * t244 + t329;
t256 = -qJ(6) * t192 - qJD(6) * t200;
t43 = pkin(4) * t307 - t50;
t255 = Ifges(4,5) * t158 + Ifges(5,5) * t103 + Ifges(4,6) * t157 - Ifges(5,6) * t102 + t281 * t375;
t85 = mrSges(7,1) * t253 - mrSges(7,2) * t252;
t15 = -pkin(4) * t281 - t17;
t240 = Ifges(4,5) * t297;
t227 = Ifges(3,5) * t280;
t224 = Ifges(4,1) * t245 + t326;
t221 = Ifges(4,2) * t248 + t327;
t215 = -mrSges(6,1) * t247 + mrSges(6,2) * t244;
t214 = -mrSges(7,1) * t247 + mrSges(7,2) * t244;
t213 = -t247 * pkin(5) + t234;
t212 = (Ifges(4,1) * t248 - t327) * qJD(3);
t209 = (-Ifges(4,2) * t245 + t326) * qJD(3);
t204 = (mrSges(4,1) * t245 + mrSges(4,2) * t248) * qJD(3);
t203 = t261 * qJD(5);
t198 = t305 * t247;
t197 = t305 * t244;
t189 = Ifges(5,5) * t192;
t188 = Ifges(5,6) * t191;
t168 = -qJD(6) * t244 - t247 * t266;
t167 = qJD(6) * t247 - t244 * t266;
t162 = -mrSges(4,1) * t307 - t194 * mrSges(4,3);
t161 = mrSges(4,2) * t307 + t193 * mrSges(4,3);
t153 = Ifges(5,1) * t200 + Ifges(5,4) * t251;
t152 = Ifges(5,4) * t200 + Ifges(5,2) * t251;
t151 = -mrSges(5,1) * t251 + mrSges(5,2) * t200;
t149 = -mrSges(6,1) * t251 - mrSges(6,3) * t312;
t148 = -mrSges(7,1) * t251 - mrSges(7,3) * t312;
t147 = mrSges(6,2) * t251 - mrSges(6,3) * t313;
t146 = mrSges(7,2) * t251 - mrSges(7,3) * t313;
t144 = t261 * t200;
t143 = (mrSges(7,1) * t244 + mrSges(7,2) * t247) * t200;
t142 = Ifges(5,1) * t192 - Ifges(5,4) * t191;
t141 = Ifges(5,4) * t192 - Ifges(5,2) * t191;
t133 = mrSges(4,1) * t281 - mrSges(4,3) * t158;
t132 = -mrSges(4,2) * t281 + mrSges(4,3) * t157;
t128 = Ifges(4,1) * t194 + Ifges(4,4) * t193 - Ifges(4,5) * t307;
t127 = Ifges(4,4) * t194 + Ifges(4,2) * t193 - t290;
t126 = pkin(5) * t313 + t159;
t119 = -mrSges(5,1) * t307 - t136 * mrSges(5,3);
t118 = mrSges(5,2) * t307 - t135 * mrSges(5,3);
t113 = -Ifges(6,3) * t251 + (Ifges(6,5) * t247 - Ifges(6,6) * t244) * t200;
t112 = -Ifges(7,3) * t251 + (Ifges(7,5) * t247 - Ifges(7,6) * t244) * t200;
t109 = -mrSges(6,2) * t191 - mrSges(6,3) * t253;
t108 = -mrSges(7,2) * t191 - mrSges(7,3) * t253;
t107 = mrSges(6,1) * t191 + mrSges(6,3) * t252;
t106 = mrSges(7,1) * t191 + mrSges(7,3) * t252;
t105 = -mrSges(4,1) * t157 + mrSges(4,2) * t158;
t92 = Ifges(4,1) * t158 + Ifges(4,4) * t157 + Ifges(4,5) * t281;
t91 = Ifges(4,4) * t158 + Ifges(4,2) * t157 + Ifges(4,6) * t281;
t88 = mrSges(5,1) * t281 - mrSges(5,3) * t103;
t87 = -mrSges(5,2) * t281 - mrSges(5,3) * t102;
t86 = mrSges(6,1) * t253 - mrSges(6,2) * t252;
t84 = pkin(5) * t253 + t129;
t83 = mrSges(5,1) * t135 + mrSges(5,2) * t136;
t82 = -qJ(6) * t313 + t94;
t81 = Ifges(5,1) * t136 - Ifges(5,4) * t135 - Ifges(5,5) * t307;
t80 = Ifges(5,4) * t136 - Ifges(5,2) * t135 - Ifges(5,6) * t307;
t77 = mrSges(6,1) * t135 + mrSges(6,3) * t254;
t76 = mrSges(7,1) * t135 + mrSges(7,3) * t254;
t75 = -mrSges(6,2) * t135 + mrSges(6,3) * t110;
t74 = -mrSges(7,2) * t135 + mrSges(7,3) * t110;
t69 = -Ifges(6,5) * t276 - Ifges(6,6) * t253 + t301;
t68 = -Ifges(7,5) * t276 - Ifges(7,6) * t253 + t302;
t67 = -pkin(5) * t251 - qJ(6) * t312 + t93;
t66 = -mrSges(6,1) * t110 - mrSges(6,2) * t254;
t65 = -mrSges(7,1) * t110 - mrSges(7,2) * t254;
t47 = Ifges(5,1) * t103 - Ifges(5,4) * t102 + t265;
t46 = Ifges(5,4) * t103 - Ifges(5,2) * t102 + t264;
t38 = -Ifges(6,5) * t254 + Ifges(6,6) * t110 + Ifges(6,3) * t135;
t37 = -Ifges(7,5) * t254 + Ifges(7,6) * t110 + Ifges(7,3) * t135;
t29 = -t110 * pkin(5) + t43;
t28 = -mrSges(6,2) * t102 + mrSges(6,3) * t60;
t27 = -mrSges(7,2) * t102 + mrSges(7,3) * t60;
t26 = mrSges(6,1) * t102 - mrSges(6,3) * t59;
t25 = mrSges(7,1) * t102 - mrSges(7,3) * t59;
t24 = -qJ(6) * t279 + (-qJD(5) * t160 + t256) * t244 + t283;
t23 = pkin(5) * t191 + t256 * t247 + (-t156 + (qJ(6) * t200 - t145) * t244) * qJD(5) + t267;
t22 = -mrSges(6,1) * t60 + mrSges(6,2) * t59;
t13 = qJ(6) * t110 + t20;
t6 = pkin(5) * t135 + qJ(6) * t254 + t19;
t5 = -t60 * pkin(5) + t15;
t2 = qJ(6) * t60 + qJD(6) * t110 + t3;
t1 = pkin(5) * t102 - qJ(6) * t59 + qJD(6) * t254 + t4;
t14 = [-t372 * t254 + (mrSges(3,3) * t246 * t348 + (0.2e1 * mrSges(3,3) * t183 - t255) * t249 + ((t195 * t353 + Ifges(3,5) * t243 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t249) * t294) * t249 + (t196 * t353 + Ifges(4,5) * t194 + Ifges(5,5) * t136 - 0.2e1 * Ifges(3,6) * t243 + Ifges(4,6) * t193 - Ifges(5,6) * t135 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t246) * t294 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t375) * t307) * t246) * qJD(2)) * t242 + t193 * t91 + t194 * t92 + 0.2e1 * t177 * t105 + t158 * t128 + 0.2e1 * t78 * t161 + 0.2e1 * t79 * t162 + (t37 + t38 - t80) * t102 + t157 * t127 + t136 * t47 + 0.2e1 * t139 * t58 + 0.2e1 * t122 * t132 + 0.2e1 * t121 * t133 + 0.2e1 * m(3) * (t183 * t196 - t184 * t195) + 0.2e1 * m(4) * (t121 * t79 + t122 * t78 + t177 * t184) + 0.2e1 * t18 * t118 + 0.2e1 * t17 * t119 + 0.2e1 * t120 * t83 + t103 * t81 + 0.2e1 * t51 * t87 + 0.2e1 * t50 * t88 + 0.2e1 * t2 * t74 + 0.2e1 * t3 * t75 + 0.2e1 * t1 * t76 + 0.2e1 * t4 * t77 + 0.2e1 * t5 * t65 + 0.2e1 * t15 * t66 + 0.2e1 * t43 * t22 + 0.2e1 * t20 * t28 + 0.2e1 * t29 * t21 + 0.2e1 * t6 * t25 + 0.2e1 * t19 * t26 + 0.2e1 * t13 * t27 + (t227 - 0.2e1 * t320 - 0.2e1 * t321) * t243 + (-mrSges(4,1) * t193 + mrSges(4,2) * t194) * t348 + (t1 * t6 + t13 * t2 + t29 * t5) * t354 + (t15 * t43 + t19 * t4 + t20 * t3) * t355 + (t120 * t139 + t17 * t50 + t18 * t51) * t356 + t370 * t59 + t371 * t60 + t373 * t110 + (-t46 + t374) * t135; -t321 - t320 + t227 + (t22 - t88) * t159 + (t66 - t119) * t129 + m(6) * (t129 * t43 + t15 * t159 + t19 * t34 + t20 * t33 + t3 * t94 + t4 * t93) + m(7) * (t1 * t67 + t126 * t5 + t13 * t24 + t2 * t82 + t23 * t6 + t29 * t84) + (t114 / 0.2e1 + t115 / 0.2e1) * t60 + (t116 / 0.2e1 + t117 / 0.2e1) * t59 + ((-t79 * t245 + t78 * t248 + (-t121 * t248 - t122 * t245) * qJD(3)) * pkin(9) - pkin(2) * t184) * m(4) - (t72 / 0.2e1 + t73 / 0.2e1) * t254 - (-t18 * mrSges(5,3) - t264 / 0.2e1 + t7 / 0.2e1 + t8 / 0.2e1 - t46 / 0.2e1) * t251 + t158 * t224 / 0.2e1 + t235 * t58 + t177 * t204 + t193 * t209 / 0.2e1 + t194 * t212 / 0.2e1 + t157 * t221 / 0.2e1 + t160 * t87 + t2 * t146 + t3 * t147 + t1 * t148 + t4 * t149 + t120 * t151 + t103 * t153 / 0.2e1 + t139 * t140 + t136 * t142 / 0.2e1 + t5 * t143 + t15 * t144 + t130 * t118 + t126 * t21 - pkin(2) * t105 + t6 * t106 + t19 * t107 + t13 * t108 + t20 * t109 + t93 * t26 + t94 * t28 + t82 * t27 + t84 * t65 + t29 * t85 + t43 * t86 + t24 * t74 + t33 * t75 + t23 * t76 + t34 * t77 + t67 * t25 + (-t51 * mrSges(5,3) - t80 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1) * t191 + (t68 / 0.2e1 + t69 / 0.2e1 - t141 / 0.2e1) * t135 + ((-t189 / 0.2e1 + t188 / 0.2e1 - t240 / 0.2e1) * t249 + (Ifges(4,5) * t245 / 0.2e1 + Ifges(4,6) * t248 / 0.2e1 - Ifges(3,6)) * t299) * t242 + m(5) * (t120 * t235 - t129 * t50 + t130 * t51 + t139 * t292 - t159 * t17 + t160 * t18) + (-t50 * mrSges(5,3) + t81 / 0.2e1 + t288 * t247 + t289 * t244) * t192 + (-t17 * mrSges(5,3) + t265 / 0.2e1 + t47 / 0.2e1 + (t11 / 0.2e1 + t12 / 0.2e1) * t247 + (-t9 / 0.2e1 - t10 / 0.2e1) * t244 + (-t244 * t288 + t247 * t289) * qJD(5)) * t200 + ((t128 / 0.2e1 - pkin(9) * t162 - t121 * mrSges(4,3)) * t248 + (t290 / 0.2e1 - t127 / 0.2e1 - pkin(9) * t161 - t122 * mrSges(4,3) + pkin(3) * t83) * t245) * qJD(3) + (t70 / 0.2e1 + t71 / 0.2e1) * t110 + (-t152 / 0.2e1 + t112 / 0.2e1 + t113 / 0.2e1) * t102 + (-pkin(9) * t133 - t79 * mrSges(4,3) + t184 * mrSges(4,2) + t92 / 0.2e1) * t245 + (pkin(9) * t132 + t78 * mrSges(4,3) - t184 * mrSges(4,1) + t91 / 0.2e1) * t248; t248 * t209 + t245 * t212 + 0.2e1 * t235 * t140 - 0.2e1 * pkin(2) * t204 + t86 * t349 + 0.2e1 * t33 * t147 + 0.2e1 * t23 * t148 + 0.2e1 * t34 * t149 + 0.2e1 * t84 * t143 + t144 * t350 + 0.2e1 * t24 * t146 + 0.2e1 * t126 * t85 + 0.2e1 * t67 * t106 + 0.2e1 * t93 * t107 + 0.2e1 * t82 * t108 + 0.2e1 * t94 * t109 + (t248 * t224 + (0.2e1 * pkin(3) * t151 - t221) * t245) * qJD(3) + (t130 * t160 + t235 * t292 + t315) * t356 + (t33 * t94 + t34 * t93 + t315) * t355 + (t126 * t84 + t23 * t67 + t24 * t82) * t354 - (t130 * t352 - t141 + t68 + t69) * t251 + (t160 * t352 + t112 + t113 - t152) * t191 + (mrSges(5,3) * t349 - t244 * t367 + t247 * t303 + t153) * t192 + (mrSges(5,3) * t350 + t142 + t368 * t247 - t369 * t244 + (-t244 * t303 - t247 * t367) * qJD(5)) * t200; t255 + t65 * t291 - t365 * t254 / 0.2e1 + t234 * t22 + t213 * t21 + t5 * t214 + t15 * t215 - t197 * t25 + t198 * t27 + t29 * t202 + t43 * t203 + t168 * t76 + t167 * t74 - t78 * mrSges(4,2) + t79 * mrSges(4,1) - t18 * mrSges(5,2) + t17 * mrSges(5,1) - t1 * t319 - t26 * t311 - t20 * t287 - t6 * t284 - t13 * t285 - t19 * t286 - t77 * t277 - t75 * t278 + t88 * t282 + t28 * t310 + t2 * t328 + t87 * t334 + (t17 * t316 + t18 * t241) * t347 + t357 * t135 + t358 * t102 + m(6) * (t15 * t234 + ((-t19 * t247 - t20 * t244) * qJD(5) + t359) * t233) + t359 * mrSges(6,3) + t363 * t59 / 0.2e1 + t364 * t60 / 0.2e1 + t366 * t110 / 0.2e1 + t370 * t270 + t371 * t271 + t372 * t337 + t373 * t336 + m(7) * (-t1 * t197 + t13 * t167 + t168 * t6 + t198 * t2 + t213 * t5 + t29 * t291); t189 - t188 + t240 - t192 * mrSges(5,3) * t282 + t143 * t291 - t357 * t251 + t234 * t86 + t213 * t85 + t84 * t214 + t363 * (t200 * t271 + t306 / 0.2e1) + t364 * (-t314 / 0.2e1 - t279 / 0.2e1) - t197 * t106 + t198 * t108 + t126 * t202 + t159 * t203 + t168 * t148 + t167 * t146 + (t241 * t347 - mrSges(5,2)) * t130 + (m(6) * t234 - t316 * t347 - mrSges(5,1) + t215) * t129 + (-mrSges(4,1) * t297 + mrSges(4,2) * t298) * pkin(9) - t23 * t319 - t107 * t311 - Ifges(4,6) * t298 - t94 * t287 - t67 * t284 - t82 * t285 - t93 * t286 - t149 * t277 - t147 * t278 + m(6) * ((-t244 * t94 - t247 * t93) * qJD(5) + t360) * t233 + t109 * t310 + t24 * t328 + (-mrSges(5,3) * t334 + t358) * t191 + t360 * mrSges(6,3) + t365 * t312 / 0.2e1 - t366 * t313 / 0.2e1 + t303 * t270 + t367 * t271 + t368 * t337 + t369 * t336 + m(7) * (t126 * t291 + t167 * t82 + t168 * t67 - t197 * t23 + t198 * t24 + t213 * t84); 0.2e1 * t234 * t203 + 0.2e1 * t213 * t202 + (t167 * t198 - t168 * t197) * t354 + (t168 * t351 + (t198 * t351 + 0.2e1 * (m(7) * t213 + t214) * pkin(5) - t364) * qJD(5) + t365) * t244 + (0.2e1 * t167 * mrSges(7,3) + (-t197 * t351 + t363) * qJD(5) + t366) * t247; (t25 + t26) * t247 + (t27 + t28) * t244 + ((t74 + t75) * t247 + (-t76 - t77) * t244) * qJD(5) + m(7) * (t1 * t247 + t2 * t244 + (t13 * t247 - t244 * t6) * qJD(5)) + m(6) * (t244 * t3 + t247 * t4 + (-t19 * t244 + t20 * t247) * qJD(5)) + m(5) * t120 + t58; m(5) * t292 + (t106 + t107) * t247 + (t108 + t109) * t244 + ((t146 + t147) * t247 + (-t148 - t149) * t244) * qJD(5) + m(7) * (t23 * t247 + t24 * t244 + (-t244 * t67 + t247 * t82) * qJD(5)) + m(6) * (t244 * t33 + t247 * t34 + (-t244 * t93 + t247 * t94) * qJD(5)) + t140; m(7) * (t167 * t244 + t168 * t247 + (t197 * t244 + t198 * t247) * qJD(5)); 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t374; mrSges(6,1) * t34 + mrSges(7,1) * t23 - mrSges(6,2) * t33 - mrSges(7,2) * t24 - t376 * t314 + (m(7) * t23 + t106) * pkin(5) + (-t244 * t377 - t247 * t376) * t200 * qJD(5) + t301 + t302; -mrSges(7,2) * t167 + t238 + t239 + (mrSges(7,1) + t346) * t168 + ((-mrSges(6,1) * t233 - (mrSges(7,3) * pkin(5))) * t247 + (mrSges(6,2) * t233 - t376) * t244) * qJD(5); (-t329 + (-mrSges(6,1) - t346) * t244) * qJD(5) - t202; 0; m(7) * t5 + t21; m(7) * t84 + t85; m(7) * t291 + t202; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
