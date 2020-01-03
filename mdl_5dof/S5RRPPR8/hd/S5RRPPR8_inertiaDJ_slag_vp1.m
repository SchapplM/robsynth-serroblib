% Calculate time derivative of joint inertia matrix for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:27
% DurationCPUTime: 10.87s
% Computational Cost: add. (8484->645), mult. (15491->879), div. (0->0), fcn. (14496->8), ass. (0->312)
t238 = pkin(8) + qJ(5);
t231 = sin(t238);
t232 = cos(t238);
t244 = sin(qJ(2));
t246 = cos(qJ(2));
t284 = t231 * t246 - t232 * t244;
t412 = qJD(2) - qJD(5);
t416 = t412 * t284;
t242 = cos(pkin(8));
t225 = pkin(4) * t242 + pkin(3);
t243 = -pkin(7) - qJ(4);
t245 = sin(qJ(1));
t226 = t245 * t243;
t247 = cos(qJ(1));
t241 = sin(pkin(8));
t362 = t241 * t244;
t331 = pkin(4) * t362;
t359 = t246 * t247;
t148 = t284 * t247;
t283 = t231 * t244 + t232 * t246;
t149 = t283 * t247;
t87 = rSges(6,1) * t149 - rSges(6,2) * t148 - rSges(6,3) * t245;
t415 = t225 * t359 + t247 * t331 + t226 + t87;
t374 = Icges(3,4) * t246;
t298 = -Icges(3,2) * t244 + t374;
t160 = Icges(3,6) * t245 + t247 * t298;
t375 = Icges(3,4) * t244;
t302 = Icges(3,1) * t246 - t375;
t164 = Icges(3,5) * t245 + t247 * t302;
t285 = t160 * t244 - t164 * t246;
t269 = t285 * t245;
t159 = -Icges(3,6) * t247 + t245 * t298;
t163 = -Icges(3,5) * t247 + t245 * t302;
t286 = t159 * t244 - t163 * t246;
t270 = t286 * t247;
t370 = Icges(4,5) * t246;
t294 = Icges(4,3) * t244 + t370;
t154 = Icges(4,6) * t245 + t247 * t294;
t371 = Icges(4,5) * t244;
t300 = Icges(4,1) * t246 + t371;
t162 = Icges(4,4) * t245 + t247 * t300;
t287 = t154 * t244 + t162 * t246;
t271 = t287 * t245;
t153 = -Icges(4,6) * t247 + t245 * t294;
t161 = -Icges(4,4) * t247 + t245 * t300;
t288 = t153 * t244 + t161 * t246;
t272 = t288 * t247;
t360 = t244 * t247;
t414 = -rSges(3,2) * t360 + rSges(3,3) * t245;
t361 = t241 * t246;
t282 = -t242 * t244 + t361;
t413 = qJD(2) * t282;
t295 = Icges(3,5) * t246 - Icges(3,6) * t244;
t155 = -Icges(3,3) * t247 + t245 * t295;
t296 = Icges(4,4) * t246 + Icges(4,6) * t244;
t157 = -Icges(4,2) * t247 + t245 * t296;
t411 = 2 * m(3);
t410 = 2 * m(4);
t409 = 2 * m(5);
t408 = 2 * m(6);
t239 = t245 ^ 2;
t240 = t247 ^ 2;
t407 = m(4) / 0.2e1;
t406 = m(5) / 0.2e1;
t405 = m(6) / 0.2e1;
t404 = -pkin(2) - pkin(3);
t281 = t242 * t246 + t362;
t133 = -Icges(5,4) * t282 - Icges(5,2) * t281;
t403 = t133 / 0.2e1;
t134 = -Icges(5,1) * t282 - Icges(5,4) * t281;
t402 = t134 / 0.2e1;
t401 = -t281 / 0.2e1;
t400 = -t282 / 0.2e1;
t397 = -rSges(4,1) - pkin(2);
t210 = rSges(3,1) * t244 + rSges(3,2) * t246;
t396 = m(3) * t210;
t395 = pkin(4) * t241;
t393 = -pkin(2) - t225;
t392 = pkin(3) - t225;
t123 = t412 * t283;
t340 = qJD(1) * t245;
t72 = t123 * t247 + t284 * t340;
t73 = t247 * t416 - t283 * t340;
t391 = rSges(6,1) * t73 + rSges(6,2) * t72;
t390 = rSges(4,1) * t244;
t389 = rSges(4,2) * t247;
t388 = rSges(3,3) * t247;
t387 = rSges(5,3) * t245;
t386 = rSges(5,3) * t247;
t385 = rSges(6,3) * t247;
t234 = t245 * rSges(4,2);
t81 = Icges(6,5) * t149 - Icges(6,6) * t148 - Icges(6,3) * t245;
t384 = t245 * t81;
t382 = t247 * t81;
t380 = -rSges(4,3) - qJ(3);
t379 = -rSges(5,3) - qJ(4);
t378 = -rSges(6,3) + t243;
t260 = -t246 * t392 + t331;
t146 = t284 * t245;
t147 = t283 * t245;
t307 = -rSges(6,1) * t147 + rSges(6,2) * t146;
t86 = -t307 + t385;
t377 = (-qJ(4) - t243) * t247 + t260 * t245 + t86;
t223 = pkin(3) * t359;
t198 = -qJ(4) * t245 + t223;
t376 = -t198 + t415;
t365 = qJ(3) * t244;
t364 = qJ(4) * t247;
t363 = t225 * t244;
t181 = t281 * qJD(2);
t102 = t181 * t247 + t282 * t340;
t103 = t247 * t413 - t281 * t340;
t358 = rSges(5,1) * t103 + rSges(5,2) * t102;
t127 = -rSges(6,1) * t284 - rSges(6,2) * t283;
t330 = pkin(4) * t361;
t357 = -t244 * t392 + t127 - t330;
t172 = t282 * t247;
t173 = t281 * t247;
t356 = rSges(5,1) * t173 - rSges(5,2) * t172;
t306 = pkin(2) * t246 + t365;
t184 = t306 * t245;
t224 = pkin(2) * t359;
t185 = qJ(3) * t360 + t224;
t355 = t184 * t245 + t185 * t247;
t178 = qJD(2) * t306 - qJD(3) * t246;
t310 = rSges(4,1) * t246 + rSges(4,3) * t244;
t354 = -qJD(2) * t310 - t178;
t353 = -t185 - t198;
t208 = pkin(2) * t244 - qJ(3) * t246;
t186 = t208 * t340;
t325 = t244 * t340;
t352 = pkin(3) * t325 + t186;
t334 = qJD(2) * t247;
t322 = t246 * t334;
t338 = qJD(1) * t247;
t351 = t243 * t338 + t322 * t395;
t209 = -rSges(4,3) * t246 + t390;
t350 = -t208 - t209;
t333 = qJD(3) * t244;
t349 = qJ(3) * t322 + t247 * t333;
t348 = rSges(4,2) * t338 + rSges(4,3) * t322;
t347 = rSges(3,2) * t325 + rSges(3,3) * t338;
t336 = qJD(2) * t245;
t324 = t244 * t336;
t346 = pkin(3) * t324 + qJ(4) * t340;
t345 = pkin(1) * t247 + pkin(6) * t245;
t344 = t239 + t240;
t343 = qJD(1) * t127;
t156 = Icges(3,3) * t245 + t247 * t295;
t342 = qJD(1) * t156;
t158 = Icges(4,2) * t245 + t247 * t296;
t341 = qJD(1) * t158;
t339 = qJD(1) * t246;
t337 = qJD(2) * t244;
t335 = qJD(2) * t246;
t332 = qJD(4) * t245;
t215 = pkin(2) * t324;
t323 = t244 * t334;
t257 = -t245 * t339 - t323;
t329 = t184 * t338 + t245 * (qJD(1) * t224 + t245 * t333 - t215 + (t244 * t338 + t245 * t335) * qJ(3)) + t247 * (pkin(2) * t257 - qJ(3) * t325 + t349);
t230 = pkin(6) * t338;
t326 = t230 + t349;
t168 = rSges(4,1) * t359 + rSges(4,3) * t360 + t234;
t321 = -pkin(3) * t244 - t208;
t320 = -qJ(3) - t395;
t319 = t127 * t344;
t142 = t350 * t247;
t318 = -qJD(4) * t247 + t215;
t197 = pkin(3) * t245 * t246 + t364;
t317 = t197 * t245 + t198 * t247 + t355;
t316 = t345 + t185;
t137 = -rSges(5,1) * t282 - rSges(5,2) * t281;
t315 = -t137 + t321;
t314 = -pkin(3) * t335 - t178;
t74 = -qJD(1) * t148 + t123 * t245;
t75 = qJD(1) * t149 + t245 * t416;
t313 = -rSges(6,1) * t75 - rSges(6,2) * t74;
t83 = Icges(6,4) * t149 - Icges(6,2) * t148 - Icges(6,6) * t245;
t85 = Icges(6,1) * t149 - Icges(6,4) * t148 - Icges(6,5) * t245;
t23 = -t148 * t83 + t149 * t85 - t384;
t253 = Icges(6,5) * t75 + Icges(6,6) * t74 - Icges(6,3) * t340;
t254 = Icges(6,4) * t75 + Icges(6,2) * t74 - Icges(6,6) * t340;
t255 = Icges(6,1) * t75 + Icges(6,4) * t74 - Icges(6,5) * t340;
t80 = Icges(6,5) * t147 - Icges(6,6) * t146 + Icges(6,3) * t247;
t82 = Icges(6,4) * t147 - Icges(6,2) * t146 + Icges(6,6) * t247;
t84 = Icges(6,1) * t147 - Icges(6,4) * t146 + Icges(6,5) * t247;
t258 = -t148 * t82 + t149 * t84 - t245 * t80;
t33 = Icges(6,5) * t73 + Icges(6,6) * t72 - Icges(6,3) * t338;
t34 = Icges(6,4) * t73 + Icges(6,2) * t72 - Icges(6,6) * t338;
t35 = Icges(6,1) * t73 + Icges(6,4) * t72 - Icges(6,5) * t338;
t1 = (-t148 * t34 + t149 * t35 - t245 * t33 + t72 * t83 + t73 * t85) * t245 - (-t148 * t254 + t149 * t255 - t245 * t253 - t338 * t80 + t72 * t82 + t73 * t84) * t247 + (t23 * t247 + (t258 - t382) * t245) * qJD(1);
t22 = -t146 * t83 + t147 * t85 + t382;
t259 = -t146 * t82 + t147 * t84 + t247 * t80;
t2 = (-t146 * t34 + t147 * t35 + t247 * t33 + t74 * t83 + t75 * t85) * t245 - (-t146 * t254 + t147 * t255 + t247 * t253 - t340 * t80 + t74 * t82 + t75 * t84) * t247 + (t22 * t247 + (t259 - t384) * t245) * qJD(1);
t312 = t245 * t1 - t247 * t2;
t311 = rSges(3,1) * t246 - rSges(3,2) * t244;
t104 = -qJD(1) * t172 + t181 * t245;
t105 = qJD(1) * t173 + t245 * t413;
t309 = rSges(5,1) * t105 + rSges(5,2) * t104;
t170 = t282 * t245;
t171 = t281 * t245;
t308 = -rSges(5,1) * t171 + rSges(5,2) * t170;
t236 = t247 * pkin(6);
t251 = t244 * t320 + t246 * t393 - pkin(1);
t250 = t251 * t245;
t43 = t247 * t378 + t236 + t250 + t307;
t44 = t316 + t415;
t305 = t245 * t44 + t247 * t43;
t261 = t246 * t404 - pkin(1) - t365;
t249 = t245 * t261 + t247 * t379;
t60 = t236 + t249 + t308;
t61 = t245 * t379 + t223 + t316 + t356;
t304 = t245 * t61 + t247 * t60;
t280 = t321 - t357;
t169 = rSges(3,1) * t359 + t414;
t279 = -rSges(5,1) * t181 + rSges(5,2) * t413 + t314;
t278 = -pkin(1) - t311;
t277 = t245 * ((pkin(3) * t339 + qJD(4)) * t247 - t346) + t247 * (pkin(3) * t257 - qJ(4) * t338 - t332) + t197 * t338 + t329;
t276 = t326 - t332;
t79 = t315 * t247;
t58 = rSges(6,1) * t123 - rSges(6,2) * t416;
t274 = -qJD(2) * t260 + t314 - t58;
t273 = qJD(2) * t210;
t265 = qJD(2) * (-Icges(4,4) * t244 + Icges(4,6) * t246);
t264 = qJD(2) * (-Icges(3,5) * t244 - Icges(3,6) * t246);
t54 = t280 * t247;
t256 = t244 * t380 + t246 * t397 - pkin(1);
t252 = t256 * t245;
t196 = t311 * qJD(2);
t167 = t245 * t311 - t388;
t166 = t245 * t310 - t389;
t141 = t350 * t245;
t140 = t169 + t345;
t139 = t245 * t278 + t236 + t388;
t130 = Icges(5,1) * t181 - Icges(5,4) * t413;
t129 = Icges(5,4) * t181 - Icges(5,2) * t413;
t128 = Icges(5,5) * t181 - Icges(5,6) * t413;
t126 = -Icges(6,1) * t284 - Icges(6,4) * t283;
t125 = -Icges(6,4) * t284 - Icges(6,2) * t283;
t124 = -Icges(6,5) * t284 - Icges(6,6) * t283;
t113 = t245 * t265 + t341;
t112 = -qJD(1) * t157 + t247 * t265;
t111 = t245 * t264 + t342;
t110 = -qJD(1) * t155 + t247 * t264;
t107 = t316 + t168;
t106 = t236 + t252 + t389;
t97 = t356 - t387;
t96 = -t308 + t386;
t95 = Icges(5,1) * t173 - Icges(5,4) * t172 - Icges(5,5) * t245;
t94 = Icges(5,1) * t171 - Icges(5,4) * t170 + Icges(5,5) * t247;
t93 = Icges(5,4) * t173 - Icges(5,2) * t172 - Icges(5,6) * t245;
t92 = Icges(5,4) * t171 - Icges(5,2) * t170 + Icges(5,6) * t247;
t91 = Icges(5,5) * t173 - Icges(5,6) * t172 - Icges(5,3) * t245;
t90 = Icges(5,5) * t171 - Icges(5,6) * t170 + Icges(5,3) * t247;
t89 = t210 * t336 + ((-rSges(3,3) - pkin(6)) * t245 + t278 * t247) * qJD(1);
t88 = rSges(3,1) * t257 - rSges(3,2) * t322 - pkin(1) * t340 + t230 + t347;
t78 = t315 * t245;
t77 = qJD(1) * t142 + t245 * t354;
t76 = t209 * t340 + t247 * t354 + t186;
t71 = t156 * t245 - t247 * t285;
t70 = t155 * t245 - t270;
t69 = t158 * t245 + t247 * t287;
t68 = t157 * t245 + t272;
t67 = -t156 * t247 - t269;
t66 = -t155 * t247 - t245 * t286;
t65 = -t158 * t247 + t271;
t64 = -t157 * t247 + t245 * t288;
t59 = t166 * t245 + t168 * t247 + t355;
t57 = Icges(6,1) * t123 - Icges(6,4) * t416;
t56 = Icges(6,4) * t123 - Icges(6,2) * t416;
t55 = Icges(6,5) * t123 - Icges(6,6) * t416;
t53 = t280 * t245;
t52 = t215 + (-t333 + (t246 * t380 + t390) * qJD(2)) * t245 + ((-rSges(4,2) - pkin(6)) * t245 + t256 * t247) * qJD(1);
t51 = qJD(1) * t252 + t323 * t397 + t326 + t348;
t50 = Icges(5,1) * t105 + Icges(5,4) * t104 - Icges(5,5) * t340;
t49 = Icges(5,1) * t103 + Icges(5,4) * t102 - Icges(5,5) * t338;
t48 = Icges(5,4) * t105 + Icges(5,2) * t104 - Icges(5,6) * t340;
t47 = Icges(5,4) * t103 + Icges(5,2) * t102 - Icges(5,6) * t338;
t46 = Icges(5,5) * t105 + Icges(5,6) * t104 - Icges(5,3) * t340;
t45 = Icges(5,5) * t103 + Icges(5,6) * t102 - Icges(5,3) * t338;
t42 = -t245 * t86 - t247 * t87;
t41 = qJD(1) * t79 + t245 * t279;
t40 = t137 * t340 + t247 * t279 + t352;
t39 = -t283 * t83 - t284 * t85;
t38 = -t283 * t82 - t284 * t84;
t37 = -rSges(6,3) * t340 - t313;
t36 = -rSges(6,3) * t338 + t391;
t32 = t245 * t96 + t247 * t97 + t317;
t31 = -t124 * t245 - t125 * t148 + t126 * t149;
t30 = t124 * t247 - t125 * t146 + t126 * t147;
t29 = (-qJ(3) * t335 - t333) * t245 + ((rSges(5,3) - pkin(6)) * t245 + t261 * t247) * qJD(1) - t309 + t318 + t346;
t28 = qJD(1) * t249 + t323 * t404 + t276 + t358;
t27 = -t172 * t93 + t173 * t95 - t245 * t91;
t26 = -t172 * t92 + t173 * t94 - t245 * t90;
t25 = -t170 * t93 + t171 * t95 + t247 * t91;
t24 = -t170 * t92 + t171 * t94 + t247 * t90;
t19 = qJD(1) * t54 + t245 * t274;
t18 = t247 * t274 + t340 * t357 + t352;
t17 = t247 * t348 + (-t209 * t239 - t240 * t390) * qJD(2) + (t247 * t166 + (-t168 - t185 + t234) * t245) * qJD(1) + t329;
t16 = t245 * t377 + t247 * t376 + t317;
t15 = (-t333 + (t246 * t320 + t363) * qJD(2)) * t245 + ((-pkin(6) - t378) * t245 + t251 * t247) * qJD(1) + t313 + t318;
t14 = t393 * t323 + (t250 - t385) * qJD(1) + t276 + t351 + t391;
t13 = t23 * t245 - t247 * t258;
t12 = t22 * t245 - t247 * t259;
t11 = -t245 * t37 - t247 * t36 + (t245 * t87 - t247 * t86) * qJD(1);
t10 = t245 * t309 + t247 * t358 + ((t96 - t386) * t247 + (t353 - t97 - t387) * t245) * qJD(1) + t277;
t9 = t123 * t85 - t283 * t34 - t284 * t35 - t416 * t83;
t8 = t123 * t84 - t254 * t283 - t255 * t284 - t416 * t82;
t7 = -t124 * t340 + t125 * t74 + t126 * t75 - t146 * t56 + t147 * t57 + t247 * t55;
t6 = -t124 * t338 + t125 * t72 + t126 * t73 - t148 * t56 + t149 * t57 - t245 * t55;
t3 = (t323 * t392 + t351 + t36) * t247 + (t37 + (t330 - t363) * t336 + t346) * t245 + ((t364 + t377) * t247 + (t353 + t226 - t376) * t245) * qJD(1) + t277;
t4 = [t123 * t126 - t284 * t57 - t416 * t125 - t283 * t56 + (t14 * t44 + t15 * t43) * t408 + (t106 * t52 + t107 * t51) * t410 + (t28 * t61 + t29 * t60) * t409 + (t139 * t89 + t140 * t88) * t411 + t181 * t134 - t282 * t130 - t413 * t133 - t281 * t129 + (t300 + t302 + t371 - t375 + (-Icges(3,2) - Icges(4,3)) * t246) * t337 + (-t294 + t298 - t370 + t374 + (Icges(3,1) + Icges(4,1)) * t244) * t335; m(3) * ((-t245 * t88 - t247 * t89) * t210 + (-t139 * t247 - t140 * t245) * t196) + m(4) * (t106 * t76 + t107 * t77 + t141 * t51 + t142 * t52) + m(5) * (t28 * t78 + t29 * t79 + t40 * t60 + t41 * t61) + m(6) * (t14 * t53 + t15 * t54 + t18 * t43 + t19 * t44) + ((-t140 * t396 - t172 * t403 + t173 * t402 + t93 * t401 + t95 * t400 + t31 / 0.2e1 + t39 / 0.2e1) * t247 + (t139 * t396 - t170 * t403 + t171 * t402 + t92 * t401 + t94 * t400 + t30 / 0.2e1 + t38 / 0.2e1) * t245) * qJD(1) + (t102 * t133 + t103 * t134 - t245 * t128 - t172 * t129 + t173 * t130 + t181 * t95 - t281 * t47 - t282 * t49 - t413 * t93 + t6 + t9) * t245 / 0.2e1 - (t104 * t133 + t105 * t134 + t128 * t247 - t170 * t129 + t171 * t130 + t181 * t94 - t281 * t48 - t282 * t50 - t413 * t92 + t7 + t8) * t247 / 0.2e1 + (t271 / 0.2e1 - t269 / 0.2e1 - t272 / 0.2e1 + t270 / 0.2e1 + (t295 + t296) * (t240 / 0.2e1 + t239 / 0.2e1)) * qJD(2); (t16 * t3 + t18 * t54 + t19 * t53) * t408 + (t10 * t32 + t40 * t79 + t41 * t78) * t409 + (t141 * t77 + t142 * t76 + t17 * t59) * t410 - t247 * ((t104 * t93 + t105 * t95 - t170 * t47 + t171 * t49 + t247 * t45) * t245 - (t104 * t92 + t105 * t94 - t170 * t48 + t171 * t50 + t247 * t46) * t247 + (t24 * t245 + t25 * t247) * qJD(1)) - t247 * ((t247 * t113 + (t65 - t272) * qJD(1)) * t247 + (t64 * qJD(1) + (t154 * t335 - t162 * t337 + t341) * t245 + (-t112 + (-t153 * t246 + t161 * t244) * qJD(2) + t287 * qJD(1)) * t247) * t245) + t245 * ((t245 * t112 + (t68 - t271) * qJD(1)) * t245 + (t69 * qJD(1) + (-t153 * t335 + t161 * t337) * t247 + (-t113 + (t154 * t246 - t162 * t244) * qJD(2) + (t158 + t288) * qJD(1)) * t245) * t247) + t245 * ((t102 * t93 + t103 * t95 - t172 * t47 + t173 * t49 - t245 * t45) * t245 - (t102 * t92 + t103 * t94 - t172 * t48 + t173 * t50 - t245 * t46) * t247 + (t26 * t245 + t247 * t27) * qJD(1)) + t245 * ((t245 * t110 + (t70 + t269) * qJD(1)) * t245 + (t71 * qJD(1) + (t159 * t335 + t163 * t337) * t247 + (-t111 + (-t160 * t246 - t164 * t244) * qJD(2) + (t156 - t286) * qJD(1)) * t245) * t247) + ((t167 * t245 + t169 * t247) * ((qJD(1) * t167 - t247 * t273 + t347) * t247 + (-t245 * t273 + (-t169 + t414) * qJD(1)) * t245) + t344 * t210 * t196) * t411 - t247 * ((t247 * t111 + (t67 + t270) * qJD(1)) * t247 + (t66 * qJD(1) + (-t160 * t335 - t164 * t337 + t342) * t245 + (-t110 + (t159 * t246 + t163 * t244) * qJD(2) - t285 * qJD(1)) * t247) * t245) + t312 + (t12 + (-t24 - t64 - t66) * t247 + (t25 + t65 + t67) * t245) * t340 + (t13 + (-t26 - t68 - t70) * t247 + (t27 + t69 + t71) * t245) * t338; 0.2e1 * (t305 * t405 + (t106 * t247 + t107 * t245) * t407 + t304 * t406) * t335 + 0.2e1 * ((t14 * t245 + t15 * t247 + t338 * t44 - t340 * t43) * t405 + (-t106 * t340 + t107 * t338 + t245 * t51 + t247 * t52) * t407 + (t245 * t28 + t247 * t29 + t338 * t61 - t340 * t60) * t406) * t244; 0.2e1 * ((t334 * t54 + t336 * t53 - t3) * t405 + (t334 * t79 + t336 * t78 - t10) * t406 + (t141 * t336 + t142 * t334 - t17) * t407) * t246 + 0.2e1 * ((qJD(2) * t16 + t18 * t247 + t19 * t245 + t338 * t53 - t340 * t54) * t405 + (qJD(2) * t32 + t245 * t41 + t247 * t40 + t338 * t78 - t340 * t79) * t406 + (qJD(2) * t59 + t141 * t338 - t142 * t340 + t245 * t77 + t247 * t76) * t407) * t244; 0.4e1 * (t407 + t406 + t405) * (-0.1e1 + t344) * t244 * t335; m(6) * (-qJD(1) * t305 + t14 * t247 - t15 * t245) + m(5) * (-qJD(1) * t304 - t245 * t29 + t247 * t28); m(6) * (-t245 * t18 + t19 * t247 + (-t245 * t53 - t247 * t54) * qJD(1)) + m(5) * (-t245 * t40 + t247 * t41 + (-t245 * t78 - t247 * t79) * qJD(1)); 0; 0; (m(6) * (t127 * t15 + t343 * t44 + t43 * t58) + t7 / 0.2e1 + t8 / 0.2e1) * t247 + (m(6) * (t127 * t14 - t343 * t43 + t44 * t58) - t6 / 0.2e1 - t9 / 0.2e1) * t245 - ((t31 + t39) * t247 + (t30 + t38) * t245) * qJD(1) / 0.2e1; m(6) * (t11 * t16 + t3 * t42) + (m(6) * (t127 * t18 + t343 * t53 + t54 * t58) + t2 - qJD(1) * t13) * t247 + (m(6) * (t127 * t19 - t343 * t54 + t53 * t58) - qJD(1) * t12 - t1) * t245; m(6) * (-t11 * t246 + t344 * t58 * t244 + (t244 * t42 + t246 * t319) * qJD(2)); 0; (t11 * t42 + t319 * t58) * t408 + (t245 * t12 + t247 * t13) * qJD(1) + t312;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
