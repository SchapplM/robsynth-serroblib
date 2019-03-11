% Calculate time derivative of joint inertia matrix for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:52
% EndTime: 2019-03-09 07:28:12
% DurationCPUTime: 8.54s
% Computational Cost: add. (24317->633), mult. (69148->949), div. (0->0), fcn. (74999->14), ass. (0->271)
t244 = sin(qJ(6));
t248 = cos(qJ(6));
t242 = cos(pkin(7));
t246 = sin(qJ(4));
t250 = cos(qJ(4));
t239 = sin(pkin(7));
t247 = sin(qJ(3));
t315 = t239 * t247;
t199 = t242 * t250 - t246 * t315;
t251 = cos(qJ(3));
t300 = qJD(3) * t251;
t283 = t239 * t300;
t188 = qJD(4) * t199 + t250 * t283;
t200 = t242 * t246 + t250 * t315;
t189 = -qJD(4) * t200 - t246 * t283;
t245 = sin(qJ(5));
t249 = cos(qJ(5));
t264 = t249 * t199 - t200 * t245;
t110 = qJD(5) * t264 + t188 * t249 + t189 * t245;
t156 = t199 * t245 + t200 * t249;
t314 = t239 * t251;
t146 = -t244 * t156 - t248 * t314;
t301 = qJD(3) * t247;
t284 = t239 * t301;
t77 = qJD(6) * t146 + t248 * t110 + t244 * t284;
t263 = -t248 * t156 + t244 * t314;
t78 = qJD(6) * t263 - t244 * t110 + t248 * t284;
t382 = -t78 * t244 + t248 * t77;
t374 = (t244 ^ 2 + t248 ^ 2) * t249;
t217 = -mrSges(7,1) * t248 + mrSges(7,2) * t244;
t381 = t217 - mrSges(6,1);
t243 = cos(pkin(6));
t238 = sin(pkin(13));
t240 = sin(pkin(6));
t241 = cos(pkin(13));
t313 = t240 * t241;
t351 = pkin(1) * t243;
t303 = qJ(2) * t313 + t238 * t351;
t169 = (t239 * t243 + t242 * t313) * pkin(9) + t303;
t228 = t241 * t351;
t317 = t238 * t240;
t179 = pkin(2) * t243 + t228 + (-pkin(9) * t242 - qJ(2)) * t317;
t196 = (-pkin(9) * t238 * t239 - pkin(2) * t241 - pkin(1)) * t240;
t265 = t179 * t242 + t196 * t239;
t113 = -t247 * t169 + t265 * t251;
t380 = (-t146 * t248 + t244 * t263) * qJD(6) + t382;
t205 = t245 * t246 - t249 * t250;
t373 = qJD(4) + qJD(5);
t180 = t373 * t205;
t206 = t245 * t250 + t246 * t249;
t295 = qJD(6) * t248;
t262 = -t244 * t180 + t206 * t295;
t296 = qJD(6) * t244;
t307 = t248 * t180;
t261 = t206 * t296 + t307;
t181 = t373 * t206;
t299 = qJD(4) * t246;
t292 = pkin(4) * t299;
t128 = pkin(5) * t181 + pkin(12) * t180 + t292;
t232 = -pkin(4) * t250 - pkin(3);
t170 = pkin(5) * t205 - pkin(12) * t206 + t232;
t362 = -pkin(11) - pkin(10);
t225 = t362 * t250;
t290 = t362 * t246;
t195 = -t249 * t225 + t245 * t290;
t130 = t170 * t248 - t195 * t244;
t194 = -t245 * t225 - t249 * t290;
t286 = qJD(4) * t362;
t216 = t246 * t286;
t277 = t250 * t286;
t138 = -qJD(5) * t194 + t249 * t216 + t245 * t277;
t72 = qJD(6) * t130 + t128 * t244 + t138 * t248;
t131 = t170 * t244 + t195 * t248;
t73 = -qJD(6) * t131 + t128 * t248 - t138 * t244;
t379 = -t73 * t244 + t248 * t72;
t312 = t242 * t247;
t178 = t243 * t315 + (t238 * t251 + t241 * t312) * t240;
t198 = -t239 * t313 + t242 * t243;
t141 = t178 * t250 + t198 * t246;
t311 = t242 * t251;
t377 = (-t238 * t247 + t241 * t311) * t240 + t243 * t314;
t135 = -t179 * t239 + t242 * t196;
t101 = -pkin(3) * t377 - pkin(10) * t178 + t135;
t160 = t251 * t169;
t114 = t179 * t312 + t196 * t315 + t160;
t104 = pkin(10) * t198 + t114;
t65 = t250 * t101 - t104 * t246;
t51 = -pkin(4) * t377 - pkin(11) * t141 + t65;
t140 = -t178 * t246 + t198 * t250;
t66 = t246 * t101 + t250 * t104;
t60 = pkin(11) * t140 + t66;
t347 = t245 * t51 + t249 * t60;
t23 = -pkin(12) * t377 + t347;
t106 = t140 * t245 + t141 * t249;
t266 = t249 * t140 - t141 * t245;
t103 = -t198 * pkin(3) - t113;
t75 = -t140 * pkin(4) + t103;
t46 = -pkin(5) * t266 - t106 * pkin(12) + t75;
t15 = -t23 * t244 + t248 * t46;
t165 = t377 * qJD(3);
t118 = -qJD(4) * t141 - t165 * t246;
t119 = qJD(4) * t140 + t165 * t250;
t56 = qJD(5) * t266 + t118 * t245 + t119 * t249;
t57 = qJD(5) * t106 - t249 * t118 + t119 * t245;
t302 = qJD(2) * t240;
t96 = (t238 * t311 + t241 * t247) * t302 + (t247 * t265 + t160) * qJD(3);
t70 = -t118 * pkin(4) + t96;
t18 = t57 * pkin(5) - t56 * pkin(12) + t70;
t166 = t178 * qJD(3);
t297 = qJD(5) * t249;
t285 = t238 * t302;
t276 = t239 * t285;
t124 = pkin(3) * t166 - pkin(10) * t165 + t276;
t95 = (-t238 * t312 + t241 * t251) * t302 + t113 * qJD(3);
t40 = -qJD(4) * t66 + t250 * t124 - t246 * t95;
t30 = pkin(4) * t166 - pkin(11) * t119 + t40;
t298 = qJD(4) * t250;
t39 = t101 * t298 - t104 * t299 + t246 * t124 + t250 * t95;
t32 = pkin(11) * t118 + t39;
t333 = t245 * t60;
t8 = -qJD(5) * t333 + t245 * t30 + t249 * t32 + t51 * t297;
t5 = pkin(12) * t166 + t8;
t2 = qJD(6) * t15 + t18 * t244 + t248 * t5;
t16 = t23 * t248 + t244 * t46;
t3 = -qJD(6) * t16 + t18 * t248 - t244 * t5;
t378 = t2 * t248 - t3 * t244;
t376 = -2 * Ifges(4,4);
t354 = t244 / 0.2e1;
t353 = t248 / 0.2e1;
t280 = -t296 / 0.2e1;
t9 = -qJD(5) * t347 - t245 * t32 + t249 * t30;
t372 = 2 * m(6);
t371 = 2 * m(7);
t370 = 0.2e1 * pkin(4);
t369 = -2 * mrSges(6,3);
t139 = qJD(5) * t195 + t245 * t216 - t249 * t277;
t368 = 0.2e1 * t139;
t367 = 0.2e1 * t194;
t366 = m(6) / 0.2e1;
t365 = m(5) * pkin(3);
t82 = t106 * t248 - t244 * t377;
t38 = -qJD(6) * t82 + t166 * t248 - t244 * t56;
t364 = t38 / 0.2e1;
t81 = -t106 * t244 - t248 * t377;
t363 = t81 / 0.2e1;
t233 = Ifges(7,5) * t295;
t360 = Ifges(7,6) * t280 + t233 / 0.2e1;
t341 = Ifges(7,4) * t244;
t272 = Ifges(7,1) * t248 - t341;
t213 = t272 * qJD(6);
t359 = t213 / 0.2e1;
t358 = Ifges(7,5) * t354 + Ifges(7,6) * t353;
t340 = Ifges(7,4) * t248;
t222 = Ifges(7,1) * t244 + t340;
t356 = t222 / 0.2e1;
t355 = -t244 / 0.2e1;
t352 = t250 / 0.2e1;
t37 = qJD(6) * t81 + t166 * t244 + t248 * t56;
t17 = -mrSges(7,1) * t38 + mrSges(7,2) * t37;
t47 = mrSges(6,1) * t166 - mrSges(6,3) * t56;
t348 = t17 - t47;
t59 = -mrSges(7,1) * t81 + mrSges(7,2) * t82;
t84 = -mrSges(6,1) * t377 - mrSges(6,3) * t106;
t346 = t59 - t84;
t345 = mrSges(4,3) * t165;
t344 = mrSges(4,3) * t166;
t343 = Ifges(5,4) * t246;
t342 = Ifges(5,4) * t250;
t339 = Ifges(7,6) * t244;
t338 = pkin(4) * qJD(5);
t337 = t166 * Ifges(6,5);
t336 = t166 * Ifges(6,6);
t335 = t377 * Ifges(5,6);
t334 = t245 * mrSges(6,1);
t330 = t249 * mrSges(6,2);
t329 = t251 * t96;
t218 = -mrSges(5,1) * t250 + mrSges(5,2) * t246;
t326 = t218 - mrSges(4,1);
t325 = t139 * t194;
t111 = qJD(5) * t156 + t188 * t245 - t249 * t189;
t324 = t264 * t111;
t323 = t264 * t245;
t322 = t194 * t245;
t321 = t206 * t244;
t320 = t206 * t248;
t309 = t244 * t249;
t308 = t245 * t217;
t306 = t248 * t249;
t305 = Ifges(4,5) * t165 - Ifges(4,6) * t166;
t304 = -Ifges(6,5) * t180 - Ifges(6,6) * t181;
t12 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t57;
t293 = Ifges(6,5) * t56 - Ifges(6,6) * t57 + Ifges(6,3) * t166;
t6 = -pkin(5) * t166 - t9;
t291 = m(7) * t6 + t17;
t287 = Ifges(5,5) * t119 + Ifges(5,6) * t118 + Ifges(5,3) * t166;
t279 = t295 / 0.2e1;
t112 = mrSges(7,1) * t262 - mrSges(7,2) * t261;
t278 = m(7) * t139 + t112;
t275 = t239 ^ 2 * t247 * t300;
t274 = mrSges(7,3) * t374;
t273 = mrSges(7,1) * t244 + mrSges(7,2) * t248;
t271 = -Ifges(7,2) * t244 + t340;
t26 = t249 * t51 - t333;
t268 = -t246 * t40 + t250 * t39;
t267 = t194 * t111 - t139 * t264;
t211 = t271 * qJD(6);
t220 = Ifges(7,2) * t248 + t341;
t259 = t248 * t211 + t244 * t213 - t220 * t296 + t222 * t295;
t89 = -Ifges(7,5) * t261 - Ifges(7,6) * t262 + Ifges(7,3) * t181;
t208 = t273 * qJD(6);
t256 = -t110 * mrSges(6,2) + t380 * mrSges(7,3) + t111 * t381 - t264 * t208;
t19 = mrSges(7,1) * t57 - mrSges(7,3) * t37;
t20 = -mrSges(7,2) * t57 + mrSges(7,3) * t38;
t61 = mrSges(7,2) * t266 + mrSges(7,3) * t81;
t62 = -mrSges(7,1) * t266 - mrSges(7,3) * t82;
t255 = m(7) * (-t15 * t295 - t16 * t296 + t378) + t248 * t20 - t244 * t19 - t62 * t295 - t61 * t296;
t125 = mrSges(7,1) * t181 + mrSges(7,3) * t261;
t126 = -mrSges(7,2) * t181 - mrSges(7,3) * t262;
t172 = -mrSges(7,2) * t205 - mrSges(7,3) * t321;
t173 = mrSges(7,1) * t205 - mrSges(7,3) * t320;
t254 = m(7) * (-t130 * t295 - t131 * t296 + t379) + t248 * t126 - t244 * t125 - t173 * t295 - t172 * t296;
t13 = Ifges(7,4) * t37 + Ifges(7,2) * t38 + Ifges(7,6) * t57;
t14 = Ifges(7,1) * t37 + Ifges(7,4) * t38 + Ifges(7,5) * t57;
t22 = pkin(5) * t377 - t26;
t44 = Ifges(7,4) * t82 + Ifges(7,2) * t81 - Ifges(7,6) * t266;
t45 = Ifges(7,1) * t82 + Ifges(7,4) * t81 - Ifges(7,5) * t266;
t253 = t9 * mrSges(6,1) - t8 * mrSges(6,2) - t266 * t360 + t13 * t353 + t14 * t354 + t22 * t208 + t211 * t363 + t6 * t217 + t220 * t364 + t45 * t279 + t44 * t280 + t37 * t356 + t57 * t358 + t82 * t359 + t293 + ((-t15 * t248 - t16 * t244) * qJD(6) + t378) * mrSges(7,3);
t150 = Ifges(7,6) * t205 + t206 * t271;
t151 = Ifges(7,5) * t205 + t206 * t272;
t90 = -Ifges(7,4) * t261 - Ifges(7,2) * t262 + Ifges(7,6) * t181;
t91 = -Ifges(7,1) * t261 - Ifges(7,4) * t262 + Ifges(7,5) * t181;
t252 = t151 * t279 - t307 * t356 + t181 * t358 + t194 * t208 - t211 * t321 / 0.2e1 + t320 * t359 + t205 * t360 + t91 * t354 + t90 * t353 - t138 * mrSges(6,2) + t304 - t262 * t220 / 0.2e1 + (t206 * t222 + t150) * t280 + t381 * t139 + ((-t130 * t248 - t131 * t244) * qJD(6) + t379) * mrSges(7,3);
t234 = Ifges(5,5) * t298;
t231 = -pkin(4) * t249 - pkin(5);
t230 = pkin(4) * t245 + pkin(12);
t223 = Ifges(5,1) * t246 + t342;
t221 = Ifges(5,2) * t250 + t343;
t214 = (Ifges(5,1) * t250 - t343) * qJD(4);
t212 = (-Ifges(5,2) * t246 + t342) * qJD(4);
t209 = (mrSges(5,1) * t246 + mrSges(5,2) * t250) * qJD(4);
t187 = Ifges(6,1) * t206 - Ifges(6,4) * t205;
t186 = Ifges(6,4) * t206 - Ifges(6,2) * t205;
t185 = mrSges(6,1) * t205 + mrSges(6,2) * t206;
t167 = t273 * t206;
t149 = Ifges(7,3) * t205 + (Ifges(7,5) * t248 - t339) * t206;
t145 = mrSges(4,1) * t198 - mrSges(4,3) * t178;
t144 = -mrSges(4,2) * t198 + mrSges(4,3) * t377;
t134 = -Ifges(6,1) * t180 - Ifges(6,4) * t181;
t133 = -Ifges(6,4) * t180 - Ifges(6,2) * t181;
t132 = mrSges(6,1) * t181 - mrSges(6,2) * t180;
t127 = mrSges(4,1) * t166 + mrSges(4,2) * t165;
t121 = -mrSges(5,1) * t377 - mrSges(5,3) * t141;
t120 = mrSges(5,2) * t377 + mrSges(5,3) * t140;
t108 = -mrSges(5,1) * t140 + mrSges(5,2) * t141;
t93 = mrSges(5,1) * t166 - mrSges(5,3) * t119;
t92 = -mrSges(5,2) * t166 + mrSges(5,3) * t118;
t86 = Ifges(5,1) * t141 + Ifges(5,4) * t140 - Ifges(5,5) * t377;
t85 = Ifges(5,4) * t141 + Ifges(5,2) * t140 - t335;
t83 = mrSges(6,2) * t377 + mrSges(6,3) * t266;
t74 = -mrSges(5,1) * t118 + mrSges(5,2) * t119;
t69 = Ifges(5,1) * t119 + Ifges(5,4) * t118 + t166 * Ifges(5,5);
t68 = Ifges(5,4) * t119 + Ifges(5,2) * t118 + t166 * Ifges(5,6);
t67 = -mrSges(6,1) * t266 + mrSges(6,2) * t106;
t64 = Ifges(6,1) * t106 + Ifges(6,4) * t266 - Ifges(6,5) * t377;
t63 = Ifges(6,4) * t106 + Ifges(6,2) * t266 - Ifges(6,6) * t377;
t48 = -mrSges(6,2) * t166 - mrSges(6,3) * t57;
t43 = Ifges(7,5) * t82 + Ifges(7,6) * t81 - Ifges(7,3) * t266;
t28 = mrSges(6,1) * t57 + mrSges(6,2) * t56;
t25 = Ifges(6,1) * t56 - Ifges(6,4) * t57 + t337;
t24 = Ifges(6,4) * t56 - Ifges(6,2) * t57 + t336;
t1 = [-(t12 - t24) * t266 + (t178 * t376 + Ifges(5,5) * t141 + Ifges(6,5) * t106 - Ifges(4,6) * t198 + Ifges(5,6) * t140 + Ifges(6,6) * t266 - ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t377) * t166 + 0.2e1 * (-mrSges(4,1) * t377 + mrSges(4,2) * t178) * t276 - t377 * t287 - t377 * t293 + (0.2e1 * Ifges(4,1) * t178 + Ifges(4,5) * t198 - t376 * t377) * t165 - 0.2e1 * t96 * t145 + t140 * t68 + t141 * t69 + 0.2e1 * t95 * t144 + 0.2e1 * t135 * t127 + 0.2e1 * t40 * t121 + t118 * t85 + t119 * t86 + 0.2e1 * t39 * t120 + t106 * t25 + 0.2e1 * t96 * t108 + 0.2e1 * t103 * t74 + 0.2e1 * t66 * t92 + 0.2e1 * t65 * t93 + t81 * t13 + t82 * t14 + 0.2e1 * t8 * t83 + 0.2e1 * t9 * t84 + 0.2e1 * t75 * t28 + 0.2e1 * t70 * t67 + 0.2e1 * t2 * t61 + 0.2e1 * t3 * t62 + t56 * t64 + 0.2e1 * t6 * t59 + t37 * t45 + 0.2e1 * t26 * t47 + t38 * t44 + 0.2e1 * t22 * t17 + 0.2e1 * t16 * t20 + 0.2e1 * t15 * t19 + 0.2e1 * (t241 * (-mrSges(3,2) * t243 + mrSges(3,3) * t313) + m(3) * (t303 * t241 + (qJ(2) * t317 - t228) * t238)) * t302 + 0.2e1 * t347 * t48 + (t26 * t9 + t347 * t8 + t70 * t75) * t372 + 0.2e1 * m(5) * (t103 * t96 + t39 * t66 + t40 * t65) + (-t63 + t43) * t57 + 0.2e1 * m(4) * (-t113 * t96 + t114 * t95 + t135 * t276) + (t15 * t3 + t16 * t2 + t22 * t6) * t371 + t198 * t305 - 0.2e1 * (mrSges(3,1) * t243 - mrSges(3,3) * t317) * t285 - 0.2e1 * t114 * t344 - 0.2e1 * t113 * t345; t110 * t83 + t188 * t120 + t189 * t121 + t242 * t127 + t146 * t19 - t263 * t20 + t156 * t48 + t199 * t93 + t200 * t92 + t77 * t61 + t78 * t62 - t348 * t264 + t346 * t111 + m(5) * (t188 * t66 + t189 * t65 + t199 * t40 + t200 * t39) + m(6) * (t110 * t347 - t111 * t26 + t156 * t8 + t264 * t9) + m(7) * (t111 * t22 + t146 * t3 + t15 * t78 + t16 * t77 - t2 * t263 - t264 * t6) + (-t247 * t344 + (-t28 - t74 - t345) * t251 + (t251 * t144 + (t108 - t145 + t67) * t247) * qJD(3) + m(4) * (-t113 * t301 + t114 * t300 + t242 * t285 + t247 * t95 - t329) + m(5) * (t103 * t301 - t329) + m(6) * (-t251 * t70 + t301 * t75)) * t239; 0.2e1 * m(7) * (t146 * t78 - t263 * t77 - t324) + 0.2e1 * m(6) * (t156 * t110 - t275 - t324) + 0.2e1 * m(5) * (t200 * t188 + t199 * t189 - t275); ((-t246 * t66 - t250 * t65) * qJD(4) + t268) * mrSges(5,3) + t305 + (t86 * t352 + (t335 / 0.2e1 - t85 / 0.2e1 + pkin(4) * t67) * t246) * qJD(4) + (t326 - t365) * t96 - (-t133 / 0.2e1 + t89 / 0.2e1) * t266 + m(7) * (t130 * t3 + t131 * t2 + t139 * t22 + t15 * t73 + t16 * t72 + t194 * t6) - (t304 + t234) * t377 / 0.2e1 + t166 * (Ifges(5,5) * t246 + Ifges(5,6) * t250) / 0.2e1 + t246 * t69 / 0.2e1 + t232 * t28 + t103 * t209 + t140 * t212 / 0.2e1 + t141 * t214 / 0.2e1 + t118 * t221 / 0.2e1 + t119 * t223 / 0.2e1 + t195 * t48 + t70 * t185 + t56 * t187 / 0.2e1 + t6 * t167 + t2 * t172 + t3 * t173 + t37 * t151 / 0.2e1 + t131 * t20 + t75 * t132 + t106 * t134 / 0.2e1 + t138 * t83 + t15 * t125 + t16 * t126 + t130 * t19 + t22 * t112 + t82 * t91 / 0.2e1 - t95 * mrSges(4,2) + t72 * t61 + t73 * t62 - pkin(3) * t74 + (-t186 / 0.2e1 + t149 / 0.2e1) * t57 + (t43 / 0.2e1 - t63 / 0.2e1 - t347 * mrSges(6,3)) * t181 + m(6) * (t138 * t347 - t139 * t26 - t194 * t9 + t195 * t8 + t232 * t70 + t292 * t75) + t68 * t352 + t90 * t363 + t150 * t364 + (-t121 * t298 + m(5) * (-t298 * t65 - t299 * t66 + t268) + t250 * t92 - t246 * t93 - t120 * t299) * pkin(10) + (t12 / 0.2e1 - t24 / 0.2e1 - t336 / 0.2e1 - t8 * mrSges(6,3)) * t205 + t346 * t139 + t348 * t194 - (t64 / 0.2e1 + t45 * t353 + t44 * t355 - t26 * mrSges(6,3)) * t180 + (t337 / 0.2e1 + t25 / 0.2e1 + t14 * t353 + t13 * t355 - t9 * mrSges(6,3) + (-t248 * t44 / 0.2e1 + t45 * t355) * qJD(6)) * t206; t111 * t167 - t264 * t112 + t146 * t125 - t263 * t126 + t77 * t172 + t78 * t173 + ((-t132 - t209) * t251 + (-mrSges(4,2) * t251 + (t185 + t326) * t247) * qJD(3)) * t239 + m(7) * (t130 * t78 + t131 * t77 + t146 * t73 - t263 * t72 + t267) + m(6) * (t195 * t110 + t138 * t156 + (t232 * t301 - t251 * t292) * t239 + t267) - t284 * t365 + (-t110 * t205 + t111 * t206 - t156 * t181 + t180 * t264) * mrSges(6,3) + (m(5) * pkin(10) + mrSges(5,3)) * (t188 * t250 - t189 * t246 + (-t199 * t250 - t200 * t246) * qJD(4)); -0.2e1 * pkin(3) * t209 + t112 * t367 + 0.2e1 * t130 * t125 + 0.2e1 * t131 * t126 + 0.2e1 * t232 * t132 + t167 * t368 + 0.2e1 * t72 * t172 + 0.2e1 * t73 * t173 + t250 * t212 + t246 * t214 + (t250 * t223 + (t185 * t370 - t221) * t246) * qJD(4) + (t138 * t195 + t232 * t292 + t325) * t372 + (t130 * t73 + t131 * t72 + t325) * t371 + (t138 * t369 - t133 + t89) * t205 + (t195 * t369 + t149 - t186) * t181 - (mrSges(6,3) * t367 - t150 * t244 + t151 * t248 + t187) * t180 + (mrSges(6,3) * t368 - t244 * t90 + t248 * t91 + t134 + (-t150 * t248 - t151 * t244) * qJD(6)) * t206; t253 + (m(6) * (t245 * t8 + t249 * t9) + t249 * t47 + t245 * t48 + (t346 * t245 + (-t244 * t62 + t248 * t61 + t83) * t249 + m(7) * (-t15 * t309 + t16 * t306 + t22 * t245) + m(6) * (-t245 * t26 + t249 * t347)) * qJD(5)) * pkin(4) + t255 * t230 + t291 * t231 - t39 * mrSges(5,2) + t40 * mrSges(5,1) + t287; m(7) * t111 * t231 - t188 * mrSges(5,2) + t189 * mrSges(5,1) + ((t110 * t245 - t111 * t249) * t366 + (m(7) * (-t146 * t309 - t263 * t306 - t323) / 0.2e1 + (t156 * t249 - t323) * t366) * qJD(5)) * t370 + t256 + m(7) * (-t146 * t295 + t263 * t296 + t382) * t230; (-Ifges(5,6) * t246 + pkin(10) * t218) * qJD(4) + t252 + t234 + (m(6) * (t138 * t245 - t139 * t249) + (t180 * t249 - t181 * t245) * mrSges(6,3) + ((mrSges(6,3) * t206 + t167) * t245 + (-t205 * mrSges(6,3) + t172 * t248 - t173 * t244) * t249 + m(7) * (-t130 * t309 + t131 * t306 + t322) + m(6) * (t195 * t249 + t322)) * qJD(5)) * pkin(4) + t278 * t231 + t254 * t230; 0.2e1 * t231 * t208 + (0.2e1 * t308 + (t374 * t230 + t231 * t245) * t371 - 0.2e1 * t330 - 0.2e1 * t334 + 0.2e1 * t274) * t338 + t259; -pkin(5) * t291 + pkin(12) * t255 + t253; m(7) * (-pkin(5) * t111 + pkin(12) * t380) + t256; -pkin(5) * t278 + pkin(12) * t254 + t252; (-pkin(5) + t231) * t208 + (m(7) * (-pkin(5) * t245 + t374 * pkin(12)) + t308 - t330 - t334 + t274) * t338 + t259; -0.2e1 * pkin(5) * t208 + t259; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t78 - mrSges(7,2) * t77; mrSges(7,1) * t73 - mrSges(7,2) * t72 + t89; t233 - t273 * pkin(4) * t297 + (t217 * t230 - t339) * qJD(6); t233 + (pkin(12) * t217 - t339) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
