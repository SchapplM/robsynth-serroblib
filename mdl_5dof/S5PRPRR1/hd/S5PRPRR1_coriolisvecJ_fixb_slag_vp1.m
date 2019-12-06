% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:52
% DurationCPUTime: 6.79s
% Computational Cost: add. (13571->523), mult. (9624->700), div. (0->0), fcn. (7492->8), ass. (0->310)
t243 = pkin(9) + qJ(4);
t240 = qJ(5) + t243;
t234 = cos(t240);
t214 = Icges(6,4) * t234;
t233 = sin(t240);
t168 = Icges(6,1) * t233 + t214;
t288 = -Icges(6,2) * t233 + t214;
t445 = t168 + t288;
t444 = 2 * qJD(4);
t236 = sin(t243);
t443 = rSges(5,2) * t236;
t244 = pkin(8) + qJ(2);
t239 = cos(t244);
t218 = t239 * qJ(3);
t237 = sin(t244);
t248 = -pkin(6) - qJ(3);
t371 = t239 * t248;
t247 = cos(pkin(9));
t235 = t247 * pkin(3) + pkin(2);
t418 = pkin(2) - t235;
t279 = t237 * t418 - t371;
t112 = -t218 + t279;
t183 = pkin(2) * t237 - t218;
t172 = qJD(2) * t183;
t442 = qJD(2) * t112 - t172;
t227 = t237 * rSges(5,3);
t238 = cos(t243);
t372 = t238 * t239;
t375 = t236 * t239;
t130 = rSges(5,1) * t372 - rSges(5,2) * t375 + t227;
t200 = t239 * t235;
t303 = -t237 * t248 + t200;
t441 = t130 + t303;
t217 = t237 * qJ(3);
t186 = t239 * pkin(2) + t217;
t413 = rSges(4,2) * sin(pkin(9));
t417 = rSges(4,1) * t247;
t280 = t237 * rSges(4,3) + (-t413 + t417) * t239;
t440 = t186 + t280;
t225 = Icges(5,4) * t238;
t289 = -Icges(5,2) * t236 + t225;
t179 = Icges(5,1) * t236 + t225;
t245 = qJD(4) + qJD(5);
t410 = rSges(6,2) * t234;
t331 = t245 * t410;
t415 = rSges(6,1) * t233;
t439 = t245 * t415 + t331;
t176 = Icges(5,5) * t238 - Icges(5,6) * t236;
t175 = Icges(5,5) * t236 + Icges(5,6) * t238;
t269 = qJD(4) * t175;
t400 = Icges(5,4) * t236;
t180 = Icges(5,1) * t238 - t400;
t128 = Icges(5,5) * t237 + t180 * t239;
t126 = Icges(5,6) * t237 + t239 * t289;
t386 = t126 * t236;
t284 = -t128 * t238 + t386;
t392 = Icges(5,3) * t239;
t438 = -t239 * t269 + (-t176 * t237 + t284 + t392) * qJD(2);
t376 = t236 * t237;
t199 = Icges(5,4) * t376;
t374 = t237 * t238;
t398 = Icges(5,5) * t239;
t127 = Icges(5,1) * t374 - t199 - t398;
t394 = Icges(5,6) * t239;
t125 = Icges(5,4) * t374 - Icges(5,2) * t376 - t394;
t387 = t125 * t236;
t285 = -t127 * t238 + t387;
t124 = Icges(5,3) * t237 + t176 * t239;
t343 = qJD(2) * t124;
t437 = qJD(2) * t285 - t237 * t269 + t343;
t399 = Icges(6,4) * t233;
t169 = Icges(6,1) * t234 - t399;
t119 = Icges(6,5) * t237 + t169 * t239;
t165 = Icges(6,5) * t234 - Icges(6,6) * t233;
t164 = Icges(6,5) * t233 + Icges(6,6) * t234;
t273 = t164 * t245;
t117 = Icges(6,6) * t237 + t239 * t288;
t389 = t117 * t233;
t391 = Icges(6,3) * t239;
t436 = -t239 * t273 + (-t119 * t234 - t165 * t237 + t389 + t391) * qJD(2);
t378 = t234 * t237;
t380 = t233 * t237;
t393 = Icges(6,6) * t239;
t116 = Icges(6,4) * t378 - Icges(6,2) * t380 - t393;
t191 = Icges(6,4) * t380;
t397 = Icges(6,5) * t239;
t118 = Icges(6,1) * t378 - t191 - t397;
t287 = t116 * t233 - t118 * t234;
t115 = Icges(6,3) * t237 + t165 * t239;
t344 = qJD(2) * t115;
t435 = qJD(2) * t287 - t237 * t273 + t344;
t123 = Icges(5,5) * t374 - Icges(5,6) * t376 - t392;
t50 = -t123 * t239 - t237 * t285;
t166 = Icges(6,2) * t234 + t399;
t282 = t166 * t233 - t168 * t234;
t434 = qJD(2) * t282 + t165 * t245;
t177 = Icges(5,2) * t238 + t400;
t281 = t236 * t177 - t238 * t179;
t433 = t281 * qJD(2) + t176 * qJD(4);
t194 = rSges(6,2) * t380;
t120 = rSges(6,1) * t378 - t239 * rSges(6,3) - t194;
t226 = t237 * rSges(6,3);
t377 = t234 * t239;
t351 = rSges(6,1) * t377 + t226;
t379 = t233 * t239;
t121 = -rSges(6,2) * t379 + t351;
t170 = t410 + t415;
t139 = t170 * t237;
t140 = t170 * t239;
t411 = rSges(6,2) * t233;
t414 = rSges(6,1) * t234;
t171 = -t411 + t414;
t173 = t237 * t245;
t174 = t239 * t245;
t421 = pkin(4) * t238;
t196 = t235 + t421;
t242 = -pkin(7) + t248;
t358 = -t237 * t196 - t239 * t242;
t93 = t235 * t237 + t358 + t371;
t156 = t239 * t196;
t345 = -t242 + t248;
t94 = t237 * t345 + t156 - t200;
t35 = t120 * t173 + t121 * t174 + qJD(1) + (-t237 * t93 + t239 * t94) * qJD(4);
t215 = qJD(3) * t237;
t338 = qJD(4) * t239;
t325 = t236 * t338;
t301 = pkin(4) * t325;
t267 = -t174 * t170 + t215 - t301;
t365 = t112 - t183;
t404 = -t120 + t93;
t40 = (t365 + t404) * qJD(2) + t267;
t409 = pkin(4) * qJD(4);
t330 = t236 * t409;
t192 = t237 * t330;
t216 = qJD(3) * t239;
t352 = -t192 - t216;
t403 = -t121 - t94;
t41 = -t170 * t173 + (t303 - t403) * qJD(2) + t352;
t300 = qJD(2) * t245;
t157 = t237 * t300;
t158 = t239 * t300;
t341 = qJD(2) * t237;
t340 = qJD(2) * t239;
t354 = rSges(6,3) * t340 + qJD(2) * t194;
t75 = -t239 * t331 + (-t174 * t233 - t234 * t341) * rSges(6,1) + t354;
t76 = -t245 * t139 + (t171 * t239 + t226) * qJD(2);
t350 = t196 - t235;
t79 = -t301 + (-t237 * t350 + t239 * t345) * qJD(2);
t204 = t248 * t341;
t373 = t237 * t242;
t80 = -t192 + t204 + (t239 * t350 - t373) * qJD(2);
t8 = t120 * t158 - t121 * t157 + t173 * t76 + t174 * t75 + (t237 * t80 + t239 * t79 + (-t237 * t94 - t239 * t93) * qJD(2)) * qJD(4);
t432 = -t40 * (qJD(2) * t139 - t174 * t171) - t35 * (-t173 * t139 - t140 * t174) - t41 * (-qJD(2) * t140 - t171 * t173) + t8 * (t237 * t120 + t239 * t121);
t431 = t237 * (-t177 * t239 + t128) - t239 * (-Icges(5,2) * t374 + t127 - t199);
t430 = qJD(2) * t445 + t173 * (-t166 * t239 + t119) - t174 * (-Icges(6,2) * t378 + t118 - t191);
t429 = t157 / 0.2e1;
t428 = t158 / 0.2e1;
t427 = -t173 / 0.2e1;
t426 = t173 / 0.2e1;
t425 = -t174 / 0.2e1;
t424 = t174 / 0.2e1;
t423 = t237 / 0.2e1;
t422 = -t239 / 0.2e1;
t420 = -qJD(2) / 0.2e1;
t419 = qJD(2) / 0.2e1;
t416 = rSges(5,1) * t238;
t412 = rSges(5,2) * t238;
t182 = rSges(5,1) * t236 + t412;
t150 = t182 * t239;
t339 = qJD(4) * t237;
t327 = t182 * t339;
t55 = qJD(2) * t441 - t216 - t327;
t408 = t150 * t55;
t349 = rSges(5,2) * t376 + t239 * rSges(5,3);
t129 = rSges(5,1) * t374 - t349;
t326 = t182 * t338;
t296 = t215 - t326;
t54 = (-t129 + t365) * qJD(2) + t296;
t407 = t237 * t54;
t114 = Icges(6,5) * t378 - Icges(6,6) * t380 - t391;
t406 = -t237 * t114 - t118 * t377;
t405 = t237 * t115 + t119 * t377;
t390 = t114 * t239;
t385 = t164 * t237;
t384 = t164 * t239;
t383 = t166 * t245;
t382 = t175 * t237;
t381 = t175 * t239;
t59 = -t237 * t282 - t384;
t370 = t59 * qJD(2);
t73 = -t237 * t281 - t381;
t369 = t73 * qJD(2);
t368 = -t237 * t123 - t127 * t372;
t367 = t237 * t124 + t128 * t372;
t148 = qJD(2) * t186 - t216;
t366 = t204 - (-t239 * t418 - t217) * qJD(2) - t148;
t337 = qJD(2) * qJD(3);
t210 = qJ(3) * t340;
t346 = t210 + t215;
t359 = qJD(2) * (-pkin(2) * t341 + t346) + t237 * t337;
t356 = -t177 + t180;
t355 = t179 + t289;
t353 = rSges(5,3) * t340 + t341 * t443;
t205 = t237 * t413;
t348 = rSges(4,3) * t340 + qJD(2) * t205;
t347 = t239 * rSges(4,3) + t205;
t342 = qJD(2) * t176;
t336 = (qJD(4) ^ 2) * t421;
t335 = t120 * t340 + t237 * t76 + t239 * t75;
t334 = t237 * t417;
t329 = qJD(2) * (qJD(2) * t279 - t210) + t359;
t324 = -pkin(2) - t417;
t322 = t341 / 0.2e1;
t321 = t340 / 0.2e1;
t320 = -t339 / 0.2e1;
t317 = t338 / 0.2e1;
t316 = -pkin(4) * t236 - t170;
t314 = -t235 - t416;
t95 = t119 * t378;
t313 = t115 * t239 - t95;
t275 = t168 * t245;
t312 = qJD(2) * t119 - t116 * t245 - t237 * t275;
t311 = -t117 * t245 - t239 * t275 + (-t169 * t237 + t397) * qJD(2);
t310 = qJD(2) * t117 + t118 * t245 - t237 * t383;
t309 = t119 * t245 - t239 * t383 + (-t237 * t288 + t393) * qJD(2);
t100 = t128 * t374;
t308 = t124 * t239 - t100;
t307 = -t114 + t389;
t306 = -t123 + t386;
t305 = t445 * t245;
t304 = t169 * t245 - t383;
t154 = t171 * t245;
t297 = -t238 * t409 - t154;
t294 = t416 - t443;
t293 = -t237 * t41 - t239 * t40;
t292 = -t237 * t55 - t239 * t54;
t57 = t116 * t234 + t118 * t233;
t64 = t125 * t238 + t127 * t236;
t65 = t126 * t238 + t128 * t236;
t283 = t129 * t237 + t130 * t239;
t149 = t182 * t237;
t51 = -t126 * t376 - t308;
t277 = (t237 * t51 - t239 * t50) * qJD(4);
t52 = -t125 * t375 - t368;
t53 = -t126 * t375 + t367;
t276 = (t237 * t53 - t239 * t52) * qJD(4);
t272 = t287 * t237;
t271 = qJD(4) * t179;
t270 = qJD(4) * t177;
t268 = qJD(2) * t165 - t173 * t384 + t174 * t385;
t266 = t125 * t239 - t126 * t237;
t257 = -t233 * t309 + t234 * t311 + t344;
t10 = t237 * t436 + t257 * t239;
t258 = qJD(2) * t114 - t233 * t310 + t234 * t312;
t11 = t258 * t237 - t239 * t435;
t12 = t257 * t237 - t239 * t436;
t44 = -t272 - t390;
t45 = -t117 * t380 - t313;
t20 = t173 * t45 - t174 * t44 + t370;
t46 = -t116 * t379 - t406;
t47 = -t117 * t379 + t405;
t60 = -t239 * t282 + t385;
t56 = t60 * qJD(2);
t21 = t173 * t47 - t174 * t46 + t56;
t261 = (-t168 * t239 - t117) * t173 - (-t168 * t237 - t116) * t174 + (-t166 + t169) * qJD(2);
t251 = -t233 * t430 + t261 * t234;
t256 = qJD(2) * t164 - t233 * t305 + t234 * t304;
t28 = t237 * t434 + t256 * t239;
t29 = t256 * t237 - t239 * t434;
t30 = t233 * t312 + t234 * t310;
t31 = t233 * t311 + t234 * t309;
t58 = t117 * t234 + t119 * t233;
t9 = t237 * t435 + t258 * t239;
t265 = (qJD(2) * t28 + t10 * t173 + t157 * t46 + t158 * t47 - t174 * t9) * t423 + (t261 * t233 + t234 * t430) * t420 + t20 * t322 + t21 * t321 + (qJD(2) * t29 - t11 * t174 + t12 * t173 + t157 * t44 + t158 * t45) * t422 + (t237 * t45 - t239 * t44) * t429 + (t237 * t47 - t239 * t46) * t428 + (t10 * t237 - t239 * t9 + (t237 * t46 + t239 * t47) * qJD(2)) * t426 + (-t11 * t239 + t12 * t237 + (t237 * t44 + t239 * t45) * qJD(2)) * t425 + (t237 * t31 - t239 * t30 + (t237 * t57 + t239 * t58) * qJD(2)) * t419 + (t237 * t268 + t239 * t251) * t427 + (t237 * t251 - t239 * t268) * t424;
t264 = (-t236 * t355 + t238 * t356) * qJD(2);
t87 = -t338 * t412 + (-t238 * t341 - t325) * rSges(5,1) + t353;
t88 = -qJD(4) * t149 + (t239 * t294 + t227) * qJD(2);
t262 = t237 * t88 + t239 * t87 + (t129 * t239 - t130 * t237) * qJD(2);
t84 = qJD(2) * t126 - t237 * t270;
t86 = qJD(2) * t128 - t237 * t271;
t255 = qJD(2) * t123 - qJD(4) * t64 - t236 * t84 + t238 * t86;
t83 = -t239 * t270 + (-t237 * t289 + t394) * qJD(2);
t85 = -t239 * t271 + (-t180 * t237 + t398) * qJD(2);
t254 = -qJD(4) * t65 - t236 * t83 + t238 * t85 + t343;
t160 = t289 * qJD(4);
t161 = t180 * qJD(4);
t253 = qJD(2) * t175 - t160 * t236 + t161 * t238 + (-t177 * t238 - t179 * t236) * qJD(4);
t252 = -t236 * t431 + t266 * t238;
t209 = t239 * t337;
t162 = t294 * qJD(4);
t132 = t334 - t347;
t92 = qJD(2) * t440 - t216;
t91 = t215 + (-t132 - t183) * qJD(2);
t78 = t209 + (-qJD(2) * t280 - t148) * qJD(2);
t77 = qJD(2) * (-qJD(2) * t334 + t348) + t359;
t74 = -t239 * t281 + t382;
t66 = t74 * qJD(2);
t61 = qJD(4) * t283 + qJD(1);
t43 = -t162 * t338 + t209 + (-t88 + t327 + t366) * qJD(2);
t42 = -t162 * t339 + (t87 - t326) * qJD(2) + t329;
t37 = t253 * t237 - t239 * t433;
t36 = t237 * t433 + t253 * t239;
t34 = -qJD(4) * t284 + t236 * t85 + t238 * t83;
t33 = -t285 * qJD(4) + t236 * t86 + t238 * t84;
t32 = t262 * qJD(4);
t25 = -t239 * t336 - t154 * t174 + t157 * t170 + t209 + (-t76 - t80 + t192 + t366) * qJD(2);
t24 = -t237 * t336 - t154 * t173 - t158 * t170 + (t75 + t79 - t301) * qJD(2) + t329;
t23 = t66 + t276;
t22 = t277 + t369;
t1 = [m(5) * t32 + m(6) * t8; (t56 + (t45 + (t116 * t239 + t117 * t237) * t233 + t313 + t406) * t174 + (-t118 * t378 + t390 + t44 + (t116 * t237 - t117 * t239) * t233 + t405) * t173) * t424 + (t66 + ((t51 - t100 + (t124 + t387) * t239 + t368) * t239 + t367 * t237) * qJD(4)) * t317 + (t57 + t59) * t429 + (t58 + t60) * t428 + (-t370 + (t47 - t272 - t405) * t174 + (t237 * t307 + t46 - t95) * t173 + ((t115 + t287) * t173 + t307 * t174) * t239 + t20) * t427 + (t31 + t28) * t426 + (-t369 + ((t239 * t306 - t367 + t53) * t239 + (t237 * t306 + t308 + t52) * t237) * qJD(4) + t22) * t320 + (t34 + t36) * t339 / 0.2e1 + (-qJD(4) * t281 + t160 * t238 + t161 * t236 + t233 * t304 + t234 * t305) * qJD(2) + (-(qJD(2) * t404 + t267 - t40 + t442) * t41 + t25 * (-t120 + t358) + t40 * (t439 * t237 - t352) + t24 * (t156 + t351 - t373) + t41 * (t215 + t354) + (-t24 * t411 + t41 * (-t330 - t439)) * t239 + ((t40 * (-t171 - t196) - t41 * t242) * t239 + (t40 * (-rSges(6,3) + t242) + t41 * (-t196 - t414)) * t237) * qJD(2)) * m(6) + (-(-qJD(2) * t129 + t296 + t442 - t54) * t55 + t43 * (t237 * t314 + t349 - t371) + t54 * (t204 + t216) + t42 * t441 + t55 * (t215 + t353) + (t182 * t407 - t408) * qJD(4) + ((-t54 * rSges(5,3) + t314 * t55) * t237 + (t54 * (-t235 - t294) - t55 * t248) * t239) * qJD(2)) * m(5) + (-(-qJD(2) * t132 - t172 + t215 - t91) * t92 + t78 * (t237 * t324 + t218 + t347) + t91 * t216 + t77 * t440 + t92 * (t346 + t348) + (t91 * (t324 + t413) * t239 + (t91 * (-rSges(4,3) - qJ(3)) + t92 * t324) * t237) * qJD(2)) * m(4) + (t21 + t30 + t29) * t425 - (t33 + t37 + t23) * t338 / 0.2e1 + ((t64 + t73) * t237 + (t65 + t74) * t239) * qJD(4) * t419; 0.2e1 * (t24 * t422 + t25 * t423) * m(6) + 0.2e1 * (t42 * t422 + t423 * t43) * m(5) + 0.2e1 * (t422 * t77 + t423 * t78) * m(4); ((t236 * t356 + t238 * t355) * qJD(2) + (t266 * t236 + t238 * t431) * qJD(4)) * t420 + ((-t338 * t382 - t342) * t239 + (t264 + (t239 * t381 + t252) * qJD(4)) * t237) * t317 + ((-t339 * t381 + t342) * t237 + (t264 + (t237 * t382 + t252) * qJD(4)) * t239) * t320 + (t237 * t34 - t239 * t33 + (t64 * t237 + t239 * t65) * qJD(2)) * t419 + t265 + (qJD(2) * t36 + (-(t237 * t437 + t255 * t239) * t239 + (t237 * t438 + t254 * t239) * t237 + (t52 * t237 + t53 * t239) * qJD(2)) * t444) * t423 + (qJD(2) * t37 + (-(t255 * t237 - t239 * t437) * t239 + (t254 * t237 - t239 * t438) * t237 + (t50 * t237 + t51 * t239) * qJD(2)) * t444) * t422 + (t277 + t22) * t322 + (t276 + t23) * t321 + (t35 * t335 + (t25 * t316 + t40 * t297 + t8 * t94 + t35 * t79 + (t316 * t41 - t35 * t93) * qJD(2)) * t239 + (t24 * t316 + t41 * t297 - t8 * t93 + t35 * t80 + (t40 * t170 + t35 * t403) * qJD(2)) * t237 - (-t41 * t236 * t340 + (t293 * t238 + t35 * (-t237 ^ 2 - t239 ^ 2) * t236) * qJD(4)) * pkin(4) + t432) * m(6) + (-(t149 * t54 - t408) * qJD(2) - (t61 * (-t149 * t237 - t150 * t239) + t292 * t294) * qJD(4) + t32 * t283 + t61 * t262 + t292 * t162 + (-t42 * t237 - t43 * t239 + (-t239 * t55 + t407) * qJD(2)) * t182) * m(5); t265 + (t35 * (-t121 * t341 + t335) + t293 * t154 + (-t24 * t237 - t25 * t239 + (t237 * t40 - t239 * t41) * qJD(2)) * t170 + t432) * m(6);];
tauc = t1(:);
