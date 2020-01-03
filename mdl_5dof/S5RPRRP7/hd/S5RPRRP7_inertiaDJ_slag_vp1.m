% Calculate time derivative of joint inertia matrix for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:34
% EndTime: 2019-12-31 18:44:50
% DurationCPUTime: 8.61s
% Computational Cost: add. (18890->638), mult. (26058->903), div. (0->0), fcn. (25250->8), ass. (0->320)
t236 = sin(qJ(3));
t239 = cos(qJ(3));
t235 = sin(qJ(4));
t238 = cos(qJ(4));
t275 = Icges(5,5) * t238 - Icges(5,6) * t235;
t188 = -Icges(5,3) * t239 + t236 * t275;
t277 = Icges(6,4) * t238 + Icges(6,6) * t235;
t189 = -Icges(6,2) * t239 + t236 * t277;
t419 = t188 + t189;
t234 = qJ(1) + pkin(8);
t231 = sin(t234);
t347 = t235 * t239;
t232 = cos(t234);
t350 = t232 * t238;
t183 = t231 * t347 + t350;
t345 = t238 * t239;
t184 = t231 * t345 - t232 * t235;
t412 = rSges(6,3) + qJ(5);
t413 = rSges(6,1) + pkin(4);
t418 = t412 * t183 + t413 * t184;
t330 = qJD(3) * t236;
t308 = t235 * t330;
t326 = qJD(4) * t238;
t328 = qJD(4) * t235;
t334 = qJD(1) * t232;
t335 = qJD(1) * t231;
t117 = -t231 * t308 - t232 * t328 - t238 * t335 + (t231 * t326 + t235 * t334) * t239;
t245 = -t238 * t330 + (-qJD(4) * t239 + qJD(1)) * t235;
t332 = qJD(1) * t239;
t299 = -qJD(4) + t332;
t118 = t231 * t245 + t299 * t350;
t417 = -t183 * qJD(5) - t412 * t117 - t413 * t118;
t353 = t231 * t236;
t121 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t353;
t125 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t353;
t129 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t353;
t352 = t231 * t238;
t185 = t232 * t347 - t352;
t321 = t232 * t345;
t186 = t231 * t235 + t321;
t351 = t232 * t236;
t53 = t121 * t351 - t125 * t185 + t129 * t186;
t122 = Icges(5,5) * t186 - Icges(5,6) * t185 + Icges(5,3) * t351;
t126 = Icges(5,4) * t186 - Icges(5,2) * t185 + Icges(5,6) * t351;
t130 = Icges(5,1) * t186 - Icges(5,4) * t185 + Icges(5,5) * t351;
t54 = t122 * t351 - t126 * t185 + t130 * t186;
t286 = t231 * t53 + t232 * t54;
t119 = Icges(6,5) * t184 + Icges(6,6) * t353 + Icges(6,3) * t183;
t123 = Icges(6,4) * t184 + Icges(6,2) * t353 + Icges(6,6) * t183;
t127 = Icges(6,1) * t184 + Icges(6,4) * t353 + Icges(6,5) * t183;
t51 = t119 * t185 + t123 * t351 + t127 * t186;
t120 = Icges(6,5) * t186 + Icges(6,6) * t351 + Icges(6,3) * t185;
t124 = Icges(6,4) * t186 + Icges(6,2) * t351 + Icges(6,6) * t185;
t128 = Icges(6,1) * t186 + Icges(6,4) * t351 + Icges(6,5) * t185;
t52 = t120 * t185 + t124 * t351 + t128 * t186;
t287 = t231 * t51 + t232 * t52;
t416 = t286 + t287;
t49 = t121 * t353 - t125 * t183 + t129 * t184;
t50 = t122 * t353 - t126 * t183 + t130 * t184;
t288 = t231 * t49 + t232 * t50;
t47 = t119 * t183 + t123 * t353 + t127 * t184;
t48 = t120 * t183 + t124 * t353 + t128 * t184;
t289 = t231 * t47 + t232 * t48;
t415 = t288 + t289;
t333 = qJD(1) * t236;
t414 = -t333 / 0.2e1;
t365 = Icges(6,5) * t238;
t274 = Icges(6,3) * t235 + t365;
t187 = -Icges(6,6) * t239 + t236 * t274;
t368 = Icges(5,4) * t238;
t278 = -Icges(5,2) * t235 + t368;
t190 = -Icges(5,6) * t239 + t236 * t278;
t366 = Icges(6,5) * t235;
t281 = Icges(6,1) * t238 + t366;
t191 = -Icges(6,4) * t239 + t236 * t281;
t369 = Icges(5,4) * t235;
t282 = Icges(5,1) * t238 - t369;
t192 = -Icges(5,5) * t239 + t236 * t282;
t411 = t419 * t239 + ((-t191 - t192) * t238 + (-t187 + t190) * t235) * t236;
t304 = t236 * t326;
t329 = qJD(3) * t239;
t307 = t235 * t329;
t247 = t304 + t307;
t272 = t120 * t235 + t128 * t238;
t59 = -t124 * t239 + t236 * t272;
t270 = -t126 * t235 + t130 * t238;
t61 = -t122 * t239 + t236 * t270;
t377 = t59 + t61;
t273 = t119 * t235 + t127 * t238;
t58 = -t123 * t239 + t236 * t273;
t271 = -t125 * t235 + t129 * t238;
t60 = -t121 * t239 + t236 * t271;
t378 = t58 + t60;
t410 = t231 * t378 + t232 * t377;
t115 = qJD(1) * t183 - qJD(4) * t321 - t231 * t328 + t232 * t308;
t116 = t232 * t245 - t299 * t352;
t309 = t232 * t329;
t312 = t231 * t333;
t248 = t309 - t312;
t311 = t231 * t329;
t249 = t232 * t333 + t311;
t68 = Icges(6,5) * t118 + Icges(6,6) * t249 + Icges(6,3) * t117;
t72 = Icges(6,4) * t118 + Icges(6,2) * t249 + Icges(6,6) * t117;
t76 = Icges(6,1) * t118 + Icges(6,4) * t249 + Icges(6,5) * t117;
t11 = -t115 * t119 + t116 * t127 + t123 * t248 + t185 * t68 + t186 * t76 + t351 * t72;
t67 = Icges(6,5) * t116 + Icges(6,6) * t248 - Icges(6,3) * t115;
t71 = Icges(6,4) * t116 + Icges(6,2) * t248 - Icges(6,6) * t115;
t75 = Icges(6,1) * t116 + Icges(6,4) * t248 - Icges(6,5) * t115;
t12 = -t115 * t120 + t116 * t128 + t124 * t248 + t185 * t67 + t186 * t75 + t351 * t71;
t70 = Icges(5,5) * t118 - Icges(5,6) * t117 + Icges(5,3) * t249;
t74 = Icges(5,4) * t118 - Icges(5,2) * t117 + Icges(5,6) * t249;
t78 = Icges(5,1) * t118 - Icges(5,4) * t117 + Icges(5,5) * t249;
t13 = t115 * t125 + t116 * t129 + t121 * t248 - t185 * t74 + t186 * t78 + t351 * t70;
t69 = Icges(5,5) * t116 + Icges(5,6) * t115 + Icges(5,3) * t248;
t73 = Icges(5,4) * t116 + Icges(5,2) * t115 + Icges(5,6) * t248;
t77 = Icges(5,1) * t116 + Icges(5,4) * t115 + Icges(5,5) * t248;
t14 = t115 * t126 + t116 * t130 + t122 * t248 - t185 * t73 + t186 * t77 + t351 * t69;
t409 = (-t11 - t13) * t232 + (t12 + t14) * t231 + t416 * qJD(1);
t15 = t117 * t119 + t118 * t127 + t123 * t249 + t183 * t68 + t184 * t76 + t353 * t72;
t16 = t117 * t120 + t118 * t128 + t124 * t249 + t183 * t67 + t184 * t75 + t353 * t71;
t17 = -t117 * t125 + t118 * t129 + t121 * t249 - t183 * t74 + t184 * t78 + t353 * t70;
t18 = -t117 * t126 + t118 * t130 + t122 * t249 - t183 * t73 + t184 * t77 + t353 * t69;
t408 = (-t15 - t17) * t232 + (t16 + t18) * t231 + t415 * qJD(1);
t19 = (qJD(3) * t273 - t72) * t239 + (qJD(3) * t123 + t235 * t68 + t238 * t76 + (t119 * t238 - t127 * t235) * qJD(4)) * t236;
t21 = (qJD(3) * t271 - t70) * t239 + (qJD(3) * t121 - t235 * t74 + t238 * t78 + (-t125 * t238 - t129 * t235) * qJD(4)) * t236;
t407 = -t19 - t21;
t20 = (qJD(3) * t272 - t71) * t239 + (qJD(3) * t124 + t235 * t67 + t238 * t75 + (t120 * t238 - t128 * t235) * qJD(4)) * t236;
t22 = (qJD(3) * t270 - t69) * t239 + (qJD(3) * t122 - t235 * t73 + t238 * t77 + (-t126 * t238 - t130 * t235) * qJD(4)) * t236;
t406 = t20 + t22;
t88 = t183 * t187 + t184 * t191 + t189 * t353;
t89 = -t183 * t190 + t184 * t192 + t188 * t353;
t405 = (-t88 - t89) * t239 + t415 * t236;
t90 = t185 * t187 + t186 * t191 + t189 * t351;
t91 = -t185 * t190 + t186 * t192 + t188 * t351;
t404 = (-t90 - t91) * t239 + t416 * t236;
t396 = 2 * m(4);
t374 = rSges(4,1) * t239;
t294 = -rSges(4,2) * t236 + t374;
t373 = rSges(4,3) * t232;
t171 = t231 * t294 - t373;
t349 = t232 * t239;
t399 = -rSges(4,2) * t351 + t231 * rSges(4,3);
t172 = rSges(4,1) * t349 + t399;
t217 = rSges(4,1) * t236 + rSges(4,2) * t239;
t256 = qJD(3) * t217;
t244 = rSges(4,2) * t312 + rSges(4,3) * t334 - t232 * t256;
t57 = (qJD(1) * t171 + t244) * t232 + (-t231 * t256 + (-t172 + t399) * qJD(1)) * t231;
t403 = t396 * t57;
t327 = qJD(4) * t236;
t151 = (Icges(6,3) * t238 - t366) * t327 + (Icges(6,6) * t236 + t239 * t274) * qJD(3);
t152 = (-Icges(5,5) * t235 - Icges(5,6) * t238) * t327 + (Icges(5,3) * t236 + t239 * t275) * qJD(3);
t153 = (-Icges(6,4) * t235 + Icges(6,6) * t238) * t327 + (Icges(6,2) * t236 + t239 * t277) * qJD(3);
t155 = (-Icges(6,1) * t235 + t365) * t327 + (Icges(6,4) * t236 + t239 * t281) * qJD(3);
t156 = (-Icges(5,1) * t235 - t368) * t327 + (Icges(5,5) * t236 + t239 * t282) * qJD(3);
t305 = t235 * t327;
t306 = t238 * t329;
t348 = t235 * t236;
t402 = -t190 * t304 + t192 * t306 + t151 * t348 + (-t305 + t306) * t191 + t247 * t187 + (t156 + t155) * t236 * t238 + t419 * t330 + (-t152 - t153) * t239;
t401 = rSges(6,2) * t309 + qJD(5) * t185 - t412 * t115 + t413 * t116;
t371 = Icges(4,4) * t236;
t284 = Icges(4,1) * t239 - t371;
t170 = Icges(4,5) * t231 + t232 * t284;
t354 = t170 * t239;
t370 = Icges(4,4) * t239;
t280 = -Icges(4,2) * t236 + t370;
t168 = Icges(4,6) * t231 + t232 * t280;
t359 = t168 * t236;
t264 = -t354 + t359;
t400 = t264 * t232;
t233 = cos(qJ(1)) * pkin(1);
t298 = -t231 * pkin(6) - t233;
t276 = Icges(4,5) * t239 - Icges(4,6) * t236;
t165 = -Icges(4,3) * t232 + t231 * t276;
t167 = -Icges(4,6) * t232 + t231 * t280;
t169 = -Icges(4,5) * t232 + t231 * t284;
t398 = t239 * t377 - t404;
t342 = rSges(6,2) * t351 + t412 * t185 + t413 * t186;
t343 = rSges(6,2) * t353 + t418;
t397 = -t231 * t343 - t232 * t342;
t395 = 2 * m(5);
t394 = 2 * m(6);
t393 = t231 ^ 2;
t392 = t232 ^ 2;
t389 = -t239 / 0.2e1;
t387 = -rSges(6,2) - pkin(7);
t386 = -rSges(5,3) - pkin(7);
t385 = m(4) * t217;
t384 = sin(qJ(1)) * pkin(1);
t383 = pkin(3) * t239;
t154 = (-Icges(5,2) * t238 - t369) * t327 + (Icges(5,6) * t236 + t239 * t278) * qJD(3);
t379 = (-t190 * t329 + (-qJD(4) * t192 - t154) * t236) * t235 + t402;
t376 = -rSges(6,2) * t312 + t401;
t375 = rSges(6,2) * t249 - t417;
t292 = -rSges(5,1) * t184 + rSges(5,2) * t183;
t132 = rSges(5,3) * t353 - t292;
t362 = t132 * t232;
t361 = t167 * t236;
t360 = t167 * t239;
t358 = t168 * t239;
t357 = t169 * t236;
t356 = t169 * t239;
t355 = t170 * t236;
t344 = t411 * t330;
t290 = rSges(6,1) * t238 + rSges(6,3) * t235;
t341 = (pkin(4) * t329 + qJ(5) * t327) * t238 + (qJ(5) * t329 + (-pkin(4) * qJD(4) + qJD(5)) * t236) * t235 + (-rSges(6,1) * t235 + rSges(6,3) * t238) * t327 + (rSges(6,2) * t236 + t239 * t290) * qJD(3);
t291 = rSges(5,1) * t238 - rSges(5,2) * t235;
t158 = (-rSges(5,1) * t235 - rSges(5,2) * t238) * t327 + (rSges(5,3) * t236 + t239 * t291) * qJD(3);
t297 = pkin(7) * t236 + t383;
t205 = t297 * qJD(3);
t340 = -t158 - t205;
t195 = t297 * t231;
t223 = pkin(3) * t349;
t196 = pkin(7) * t351 + t223;
t339 = t231 * t195 + t232 * t196;
t338 = -rSges(6,2) * t239 + (pkin(4) * t238 + qJ(5) * t235 + t290) * t236;
t194 = -rSges(5,3) * t239 + t236 * t291;
t221 = pkin(3) * t236 - pkin(7) * t239;
t337 = -t194 - t221;
t166 = Icges(4,3) * t231 + t232 * t276;
t336 = qJD(1) * t166;
t331 = qJD(3) * t231;
t31 = t231 * t48 - t232 * t47;
t32 = t231 * t50 - t232 * t49;
t323 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = t231 * t52 - t232 * t51;
t34 = t231 * t54 - t232 * t53;
t322 = t33 / 0.2e1 + t34 / 0.2e1;
t320 = t116 * rSges(5,1) + t115 * rSges(5,2) + rSges(5,3) * t309;
t209 = t231 * pkin(3) * t330;
t210 = pkin(7) * t309;
t310 = t232 * t330;
t318 = t231 * (pkin(7) * t249 + qJD(1) * t223 - t209) + t232 * (-pkin(7) * t312 + t210 + (-t231 * t332 - t310) * pkin(3)) + t195 * t334;
t317 = -t205 - t341;
t134 = t186 * rSges(5,1) - t185 * rSges(5,2) + rSges(5,3) * t351;
t316 = -t221 - t338;
t315 = t232 * pkin(2) - t298;
t314 = -pkin(2) - t383;
t313 = t194 * t335;
t303 = t231 * t338;
t302 = t343 * t232;
t301 = t338 * t232;
t160 = t337 * t232;
t300 = qJD(1) * t338;
t109 = t316 * t232;
t296 = t231 * t300;
t293 = t118 * rSges(5,1) - t117 * rSges(5,2);
t285 = t315 + t196;
t279 = Icges(4,2) * t239 + t371;
t269 = -t134 * t231 + t362;
t268 = -t132 * t231 - t134 * t232;
t265 = -t356 + t361;
t263 = -t239 * t378 + t405;
t262 = -pkin(2) - t294;
t37 = t117 * t187 + t118 * t191 + t151 * t183 + t153 * t353 + t155 * t184 + t189 * t249;
t38 = -t117 * t190 + t118 * t192 + t152 * t353 - t154 * t183 + t156 * t184 + t188 * t249;
t261 = t21 / 0.2e1 + t19 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1;
t35 = -t115 * t187 + t116 * t191 + t151 * t185 + t153 * t351 + t155 * t186 + t189 * t248;
t36 = t115 * t190 + t116 * t192 + t152 * t351 - t154 * t185 + t156 * t186 + t188 * t248;
t260 = t22 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1;
t259 = t60 / 0.2e1 + t58 / 0.2e1 + t88 / 0.2e1 + t89 / 0.2e1;
t258 = -t61 / 0.2e1 - t59 / 0.2e1 - t90 / 0.2e1 - t91 / 0.2e1;
t226 = pkin(6) * t334;
t257 = -pkin(3) * t310 + t210 + t226;
t255 = t236 * t387 + t314;
t254 = t236 * t386 + t314;
t252 = qJD(3) * t279;
t251 = qJD(3) * (-Icges(4,5) * t236 - Icges(4,6) * t239);
t246 = -t231 * t342 + t302;
t242 = t231 * t255 - t384;
t241 = t231 * t254 - t384;
t229 = t232 * pkin(6);
t202 = t294 * qJD(3);
t197 = t221 * t335;
t159 = t337 * t231;
t149 = t172 + t315;
t148 = t231 * t262 + t229 + t373 - t384;
t136 = t231 * t251 + t336;
t135 = -qJD(1) * t165 + t232 * t251;
t108 = t316 * t231;
t107 = t217 * t331 + (-t233 + (-rSges(4,3) - pkin(6)) * t231 + t262 * t232) * qJD(1);
t106 = t226 + (-t384 + (-pkin(2) - t374) * t231) * qJD(1) + t244;
t101 = -t134 * t239 - t194 * t351;
t100 = t132 * t239 + t194 * t353;
t99 = t285 + t134;
t98 = t229 + t241 + t292;
t97 = t166 * t231 - t400;
t96 = t165 * t231 - t265 * t232;
t95 = -t166 * t232 - t231 * t264;
t94 = -t165 * t232 - t231 * t265;
t93 = qJD(1) * t160 + t231 * t340;
t92 = t232 * t340 + t197 + t313;
t87 = t269 * t236;
t86 = t285 + t342;
t85 = t229 + t242 - t418;
t84 = rSges(5,3) * t249 + t293;
t82 = -rSges(5,3) * t312 + t320;
t66 = -t236 * t301 - t239 * t342;
t65 = t236 * t303 + t239 * t343;
t64 = -t268 + t339;
t63 = qJD(1) * t109 + t231 * t317;
t62 = t232 * t317 + t197 + t296;
t56 = t209 + t386 * t311 + (t232 * t254 + t298) * qJD(1) - t293;
t55 = qJD(1) * t241 + t257 + t320;
t46 = t246 * t236;
t45 = t339 - t397;
t42 = (t194 * t331 + t84) * t239 + (-qJD(3) * t132 + t158 * t231 + t194 * t334) * t236;
t41 = (-qJD(3) * t194 * t232 - t82) * t239 + (qJD(3) * t134 - t158 * t232 + t313) * t236;
t40 = t209 + t387 * t311 + (t232 * t255 + t298) * qJD(1) + t417;
t39 = qJD(1) * t242 + t257 + t401;
t30 = t269 * t329 + (qJD(1) * t268 - t231 * t82 + t232 * t84) * t236;
t25 = t231 * t84 + t232 * t82 + (t362 + (-t134 - t196) * t231) * qJD(1) + t318;
t24 = (qJD(3) * t303 + t375) * t239 + (-qJD(3) * t343 + t231 * t341 + t232 * t300) * t236;
t23 = (-qJD(3) * t301 - t376) * t239 + (qJD(3) * t342 - t232 * t341 + t296) * t236;
t10 = t376 * t232 + t375 * t231 + (t302 + (-t196 - t342) * t231) * qJD(1) + t318;
t9 = t246 * t329 + (qJD(1) * t397 - t376 * t231 + t375 * t232) * t236;
t4 = (qJD(3) * t288 - t38) * t239 + (-qJD(1) * t32 + qJD(3) * t89 + t17 * t231 + t18 * t232) * t236;
t3 = (qJD(3) * t289 - t37) * t239 + (-qJD(1) * t31 + qJD(3) * t88 + t15 * t231 + t16 * t232) * t236;
t2 = (qJD(3) * t286 - t36) * t239 + (-qJD(1) * t34 + qJD(3) * t91 + t13 * t231 + t14 * t232) * t236;
t1 = (qJD(3) * t287 - t35) * t239 + (-qJD(1) * t33 + qJD(3) * t90 + t11 * t231 + t12 * t232) * t236;
t5 = [(t55 * t99 + t56 * t98) * t395 + (t39 * t86 + t40 * t85) * t394 + (t106 * t149 + t107 * t148) * t396 - t190 * t307 - t192 * t305 - t154 * t348 + t402 + (t284 - t279) * t330 + (Icges(4,1) * t236 + t280 + t370) * t329; 0; 0; m(6) * (t108 * t39 + t109 * t40 + t62 * t85 + t63 * t86) + m(5) * (t159 * t55 + t160 * t56 + t92 * t98 + t93 * t99) + (t392 / 0.2e1 + t393 / 0.2e1) * t276 * qJD(3) + ((qJD(1) * t168 - t231 * t252) * t389 + t170 * t414 + (t361 / 0.2e1 - t356 / 0.2e1) * qJD(3) - t261 + m(4) * (-t107 * t217 - t148 * t202) + (-t149 * t385 + t358 / 0.2e1 + t355 / 0.2e1 - t258) * qJD(1)) * t232 + ((-qJD(1) * t167 - t232 * t252) * t239 / 0.2e1 + t169 * t414 + (-t359 / 0.2e1 + t354 / 0.2e1) * qJD(3) + t260 + m(4) * (-t106 * t217 - t149 * t202) + (t148 * t385 + t360 / 0.2e1 + t357 / 0.2e1 + t259) * qJD(1)) * t231; m(4) * t57 + m(5) * t25 + m(6) * t10; (t45 * t10 + t108 * t63 + t109 * t62) * t394 + (t159 * t93 + t160 * t92 + t64 * t25) * t395 + (t392 + t393) * t217 * t202 * t396 + (t172 * t403 + (-t95 * qJD(1) + (-t265 * qJD(1) - t136) * t232) * t232 - t408) * t232 + (t171 * t403 + (t96 * qJD(1) + (t264 * qJD(1) + t135) * t231) * t231 + ((-t136 + (-t355 - t358) * qJD(3) + t168 * t329 + t170 * t330 - t336) * t231 + (t167 * t329 + t169 * t330 + t135 - (t357 + t360) * qJD(3)) * t232 + (t97 - t94 + (t166 - t265) * t231 + t400) * qJD(1)) * t232 + t409) * t231 + (t231 * t95 - t232 * t94 + t31 + t32) * t335 + (t231 * t97 - t232 * t96 + t33 + t34) * t334; m(5) * (t100 * t56 + t101 * t55 + t41 * t99 + t42 * t98) + m(6) * (t23 * t86 + t24 * t85 + t66 * t39 + t65 * t40) + ((t231 * t259 - t232 * t258) * qJD(3) - t379) * t239 + (t260 * t232 + t261 * t231 + (t231 * t258 + t232 * t259) * qJD(1)) * t236 - t344; m(5) * t30 + m(6) * t9; m(6) * (t46 * t10 + t108 * t23 + t109 * t24 + t9 * t45 + t65 * t62 + t66 * t63) + m(5) * (t100 * t92 + t101 * t93 + t159 * t41 + t160 * t42 + t25 * t87 + t30 * t64) + (-t3 / 0.2e1 - t4 / 0.2e1 + t322 * t329) * t232 + (t2 / 0.2e1 + t1 / 0.2e1 + t323 * t329) * t231 + ((-t231 * t322 + t232 * t323) * qJD(1) + t408 * t231 / 0.2e1 + t409 * t232 / 0.2e1 + (t231 * t377 - t232 * t378) * qJD(3) / 0.2e1) * t236 + (t410 * qJD(1) + t406 * t231 + t407 * t232) * t389 + (t405 * t231 + t404 * t232) * qJD(1) / 0.2e1; (t23 * t66 + t24 * t65 + t46 * t9) * t394 + (t100 * t42 + t101 * t41 + t30 * t87) * t395 + (t379 * t239 + (t263 * t231 - t232 * t398) * qJD(3) + t344) * t239 + ((-t239 * t406 + t1 + t2) * t232 + (t239 * t407 + t3 + t4) * t231 + (t410 * t236 + t411 * t239) * qJD(3) + (t231 * t398 + t263 * t232) * qJD(1)) * t236; m(6) * (-t115 * t85 + t117 * t86 + t183 * t39 + t185 * t40); t247 * m(6); m(6) * (t45 * t304 + t108 * t117 - t109 * t115 + t183 * t63 + t185 * t62 + (t10 * t236 + t329 * t45) * t235); m(6) * (t46 * t304 - t115 * t65 + t117 * t66 + t183 * t23 + t185 * t24 + (t236 * t9 + t46 * t329) * t235); (-t115 * t185 + t117 * t183 + t247 * t348) * t394;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
