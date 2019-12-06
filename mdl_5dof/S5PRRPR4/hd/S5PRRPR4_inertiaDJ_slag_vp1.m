% Calculate time derivative of joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:30
% EndTime: 2019-12-05 16:22:07
% DurationCPUTime: 16.12s
% Computational Cost: add. (21740->700), mult. (28332->1043), div. (0->0), fcn. (26687->10), ass. (0->349)
t279 = qJ(3) + pkin(9);
t271 = sin(t279);
t272 = cos(t279);
t280 = sin(pkin(8));
t281 = cos(pkin(8));
t286 = cos(qJ(2));
t391 = t281 * t286;
t242 = -t271 * t391 + t280 * t272;
t243 = t280 * t271 + t272 * t391;
t284 = sin(qJ(2));
t392 = t281 * t284;
t167 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t392;
t285 = cos(qJ(3));
t283 = sin(qJ(3));
t390 = t283 * t286;
t258 = t280 * t285 - t281 * t390;
t389 = t285 * t286;
t396 = t280 * t283;
t259 = t281 * t389 + t396;
t187 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t392;
t454 = -t187 - t167;
t169 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t392;
t171 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t392;
t189 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t392;
t191 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t392;
t394 = t280 * t286;
t240 = -t271 * t394 - t281 * t272;
t241 = -t281 * t271 + t272 * t394;
t256 = -t280 * t390 - t281 * t285;
t393 = t281 * t283;
t257 = t280 * t389 - t393;
t395 = t280 * t284;
t457 = t169 * t240 + t171 * t241 + t189 * t256 + t191 * t257 - t454 * t395;
t166 = Icges(5,5) * t241 + Icges(5,6) * t240 + Icges(5,3) * t395;
t186 = Icges(4,5) * t257 + Icges(4,6) * t256 + Icges(4,3) * t395;
t455 = -t186 - t166;
t320 = Icges(5,5) * t272 - Icges(5,6) * t271;
t226 = -Icges(5,3) * t286 + t284 * t320;
t321 = Icges(4,5) * t285 - Icges(4,6) * t283;
t244 = -Icges(4,3) * t286 + t284 * t321;
t456 = t226 + t244;
t453 = -t169 * t271 + t171 * t272 - t189 * t283 + t191 * t285;
t168 = Icges(5,4) * t241 + Icges(5,2) * t240 + Icges(5,6) * t395;
t170 = Icges(5,1) * t241 + Icges(5,4) * t240 + Icges(5,5) * t395;
t188 = Icges(4,4) * t257 + Icges(4,2) * t256 + Icges(4,6) * t395;
t190 = Icges(4,1) * t257 + Icges(4,4) * t256 + Icges(4,5) * t395;
t452 = -t168 * t271 + t170 * t272 - t188 * t283 + t190 * t285;
t441 = t168 * t240 + t170 * t241 + t188 * t256 + t190 * t257 - t455 * t395;
t404 = Icges(5,4) * t272;
t323 = -Icges(5,2) * t271 + t404;
t227 = -Icges(5,6) * t286 + t284 * t323;
t405 = Icges(5,4) * t271;
t327 = Icges(5,1) * t272 - t405;
t228 = -Icges(5,5) * t286 + t284 * t327;
t406 = Icges(4,4) * t285;
t324 = -Icges(4,2) * t283 + t406;
t245 = -Icges(4,6) * t286 + t284 * t324;
t407 = Icges(4,4) * t283;
t328 = Icges(4,1) * t285 - t407;
t246 = -Icges(4,5) * t286 + t284 * t328;
t451 = t227 * t240 + t228 * t241 + t245 * t256 + t246 * t257 + t456 * t395;
t367 = qJD(3) * t284;
t386 = t286 * ((-Icges(4,5) * t283 - Icges(4,6) * t285) * t367 + (Icges(4,3) * t284 + t286 * t321) * qJD(2));
t387 = t286 * ((-Icges(5,5) * t271 - Icges(5,6) * t272) * t367 + (Icges(5,3) * t284 + t286 * t320) * qJD(2));
t450 = -t387 - t386;
t399 = t244 * t286;
t400 = t226 * t286;
t449 = -t399 - t400;
t448 = t457 * t281;
t447 = t452 * t284 + t455 * t286;
t446 = t453 * t284 + t454 * t286;
t445 = t227 * t271 - t228 * t272 + t245 * t283 - t246 * t285;
t176 = (-Icges(5,2) * t272 - t405) * t367 + (Icges(5,6) * t284 + t286 * t323) * qJD(2);
t177 = (-Icges(5,1) * t271 - t404) * t367 + (Icges(5,5) * t284 + t286 * t327) * qJD(2);
t202 = (-Icges(4,2) * t285 - t407) * t367 + (Icges(4,6) * t284 + t286 * t324) * qJD(2);
t203 = (-Icges(4,1) * t283 - t406) * t367 + (Icges(4,5) * t284 + t286 * t328) * qJD(2);
t370 = qJD(2) * t284;
t353 = t280 * t370;
t208 = -qJD(3) * t241 + t271 * t353;
t209 = qJD(3) * t240 - t272 * t353;
t349 = t283 * t370;
t216 = -qJD(3) * t257 + t280 * t349;
t348 = t285 * t370;
t217 = qJD(3) * t256 - t280 * t348;
t369 = qJD(2) * t286;
t352 = t280 * t369;
t119 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t352;
t121 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t352;
t117 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t352;
t302 = t117 * t284 + t166 * t369;
t36 = t240 * t119 + t241 * t121 + t208 * t168 + t209 * t170 + t280 * t302;
t351 = t281 * t370;
t210 = -qJD(3) * t243 + t271 * t351;
t211 = qJD(3) * t242 - t272 * t351;
t350 = t281 * t369;
t120 = Icges(5,4) * t211 + Icges(5,2) * t210 + Icges(5,6) * t350;
t122 = Icges(5,1) * t211 + Icges(5,4) * t210 + Icges(5,5) * t350;
t118 = Icges(5,5) * t211 + Icges(5,6) * t210 + Icges(5,3) * t350;
t301 = t118 * t284 + t167 * t369;
t37 = t240 * t120 + t241 * t122 + t208 * t169 + t209 * t171 + t280 * t301;
t136 = Icges(4,4) * t217 + Icges(4,2) * t216 + Icges(4,6) * t352;
t138 = Icges(4,1) * t217 + Icges(4,4) * t216 + Icges(4,5) * t352;
t134 = Icges(4,5) * t217 + Icges(4,6) * t216 + Icges(4,3) * t352;
t300 = t134 * t284 + t186 * t369;
t45 = t256 * t136 + t257 * t138 + t216 * t188 + t217 * t190 + t280 * t300;
t218 = -qJD(3) * t259 + t281 * t349;
t219 = qJD(3) * t258 - t281 * t348;
t137 = Icges(4,4) * t219 + Icges(4,2) * t218 + Icges(4,6) * t350;
t139 = Icges(4,1) * t219 + Icges(4,4) * t218 + Icges(4,5) * t350;
t135 = Icges(4,5) * t219 + Icges(4,6) * t218 + Icges(4,3) * t350;
t299 = t135 * t284 + t187 * t369;
t46 = t256 * t137 + t257 * t139 + t216 * t189 + t217 * t191 + t280 * t299;
t444 = (-t240 * t176 - t241 * t177 - t256 * t202 - t257 * t203 - t208 * t227 - t209 * t228 - t216 * t245 - t217 * t246) * t286 + ((t46 + t37) * t281 + (t45 + t36 + t450) * t280) * t284 + (((t449 + t441) * t280 + t448) * t286 + t451 * t284) * qJD(2);
t443 = (-t452 * qJD(2) + t117 + t134) * t286 + (t119 * t271 - t121 * t272 + t136 * t283 - t138 * t285 + (t168 * t272 + t170 * t271 + t188 * t285 + t190 * t283) * qJD(3) + t455 * qJD(2)) * t284;
t442 = (-t453 * qJD(2) + t118 + t135) * t286 + (t120 * t271 - t122 * t272 + t137 * t283 - t139 * t285 + (t169 * t272 + t171 * t271 + t189 * t285 + t191 * t283) * qJD(3) + t454 * qJD(2)) * t284;
t276 = t280 ^ 2;
t277 = t281 ^ 2;
t436 = t276 + t277;
t70 = t167 * t392 + t169 * t242 + t171 * t243;
t83 = t187 * t392 + t189 * t258 + t191 * t259;
t440 = t70 + t83;
t437 = t445 * t284 - t449;
t435 = t447 * t280 + t446 * t281;
t434 = qJD(2) * (t284 * rSges(3,1) + rSges(3,2) * t286);
t419 = pkin(4) * t272;
t205 = -pkin(7) * t286 + t284 * t419;
t274 = t285 * pkin(3);
t225 = -qJ(4) * t286 + t274 * t284;
t431 = 2 * m(4);
t430 = 2 * m(5);
t429 = 2 * m(6);
t428 = 0.2e1 * qJD(2);
t427 = m(5) / 0.2e1;
t426 = m(6) / 0.2e1;
t425 = t280 / 0.2e1;
t424 = -t281 / 0.2e1;
t421 = -t286 / 0.2e1;
t420 = pkin(3) * t283;
t273 = qJ(5) + t279;
t268 = sin(t273);
t269 = cos(t273);
t230 = -t268 * t394 - t281 * t269;
t231 = -t281 * t268 + t269 * t394;
t158 = Icges(6,5) * t231 + Icges(6,6) * t230 + Icges(6,3) * t395;
t160 = Icges(6,4) * t231 + Icges(6,2) * t230 + Icges(6,6) * t395;
t162 = Icges(6,1) * t231 + Icges(6,4) * t230 + Icges(6,5) * t395;
t232 = -t268 * t391 + t280 * t269;
t233 = t280 * t268 + t269 * t391;
t64 = t158 * t392 + t160 * t232 + t162 * t233;
t416 = t280 * t64;
t69 = t166 * t392 + t168 * t242 + t170 * t243;
t415 = t280 * t69;
t82 = t186 * t392 + t188 * t258 + t190 * t259;
t414 = t280 * t82;
t159 = Icges(6,5) * t233 + Icges(6,6) * t232 + Icges(6,3) * t392;
t161 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t392;
t163 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t392;
t63 = t159 * t395 + t161 * t230 + t163 * t231;
t413 = t281 * t63;
t278 = qJD(3) + qJD(5);
t298 = -t278 * t280 + t351;
t361 = t278 * t391;
t194 = t268 * t298 - t269 * t361;
t195 = -t268 * t361 - t269 * t298;
t110 = t195 * rSges(6,1) + t194 * rSges(6,2) + rSges(6,3) * t350;
t265 = pkin(4) * t271 + t420;
t368 = qJD(3) * t283;
t365 = pkin(3) * t368;
t340 = -t265 * qJD(3) + t365;
t288 = -t205 * qJD(2) + t340 * t286;
t366 = qJD(3) * t285;
t364 = pkin(3) * t366;
t339 = -(t274 + t419) * qJD(3) + t364;
t410 = -t280 * t339 + t281 * t288 + t110;
t403 = Icges(6,4) * t268;
t402 = Icges(6,4) * t269;
t319 = Icges(6,5) * t269 - Icges(6,6) * t268;
t221 = -Icges(6,3) * t286 + t284 * t319;
t401 = t221 * t286;
t397 = t278 * t284;
t388 = t286 * ((-Icges(6,5) * t268 - Icges(6,6) * t269) * t397 + (Icges(6,3) * t284 + t286 * t319) * qJD(2));
t297 = t278 * t281 + t353;
t362 = t278 * t394;
t192 = t268 * t297 - t269 * t362;
t193 = -t268 * t362 - t269 * t297;
t109 = t193 * rSges(6,1) + t192 * rSges(6,2) + rSges(6,3) * t352;
t164 = rSges(6,1) * t231 + rSges(6,2) * t230 + rSges(6,3) * t395;
t385 = t109 * t392 + t164 * t350;
t124 = t211 * rSges(5,1) + t210 * rSges(5,2) + rSges(5,3) * t350;
t287 = -t225 * qJD(2) + qJD(4) * t284 - t286 * t365;
t144 = t280 * t364 + t281 * t287;
t384 = -t124 - t144;
t143 = t280 * t287 - t281 * t364;
t290 = qJ(4) * t284 + t274 * t286;
t196 = -pkin(3) * t393 + t280 * t290;
t383 = t143 * t392 + t196 * t350;
t289 = pkin(7) * t284 + t286 * t419;
t346 = -t265 + t420;
t131 = t280 * t289 + t281 * t346;
t382 = t131 + t164;
t132 = -t280 * t346 + t281 * t289;
t165 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t392;
t381 = t132 + t165;
t173 = rSges(5,1) * t243 + rSges(5,2) * t242 + rSges(5,3) * t392;
t197 = pkin(3) * t396 + t281 * t290;
t380 = -t173 - t197;
t334 = rSges(5,1) * t272 - rSges(5,2) * t271;
t179 = (-rSges(5,1) * t271 - rSges(5,2) * t272) * t367 + (rSges(5,3) * t284 + t286 * t334) * qJD(2);
t200 = qJD(2) * t290 - qJD(4) * t286 - t284 * t365;
t379 = -t179 - t200;
t378 = t286 * t196 + t225 * t395;
t335 = rSges(4,1) * t285 - rSges(4,2) * t283;
t204 = (-rSges(4,1) * t283 - rSges(4,2) * t285) * t367 + (rSges(4,3) * t284 + t286 * t335) * qJD(2);
t337 = pkin(2) * t286 + pkin(6) * t284;
t264 = t337 * qJD(2);
t377 = -t204 - t264;
t333 = rSges(6,1) * t269 - rSges(6,2) * t268;
t224 = -rSges(6,3) * t286 + t284 * t333;
t114 = t286 * t164 + t224 * t395;
t229 = -rSges(5,3) * t286 + t284 * t334;
t376 = -t225 - t229;
t267 = t284 * pkin(2) - pkin(6) * t286;
t375 = t436 * qJD(2) * t267;
t247 = -rSges(4,3) * t286 + t284 * t335;
t374 = -t247 - t267;
t373 = t436 * t337;
t363 = -t144 - t410;
t157 = (-rSges(6,1) * t268 - rSges(6,2) * t269) * t397 + (rSges(6,3) * t284 + t286 * t333) * qJD(2);
t360 = t286 * t109 + t157 * t395 + t224 * t352;
t359 = -t197 - t381;
t358 = t286 * t143 + t200 * t395 + t225 * t352;
t142 = qJD(2) * t289 + t284 * t340;
t357 = -t142 - t157 - t200;
t356 = -t264 + t379;
t355 = -t205 - t224 - t225;
t354 = -t267 + t376;
t345 = t380 * t286;
t344 = t280 * t143 + t281 * t144 - t375;
t343 = -t264 + t357;
t342 = t280 * t196 + t281 * t197 + t373;
t341 = -t267 + t355;
t338 = t359 * t286;
t318 = -t160 * t268 + t162 * t269;
t73 = -t158 * t286 + t284 * t318;
t317 = -t161 * t268 + t163 * t269;
t74 = -t159 * t286 + t284 * t317;
t332 = t73 * t280 + t74 * t281;
t326 = Icges(6,1) * t269 - t403;
t322 = -Icges(6,2) * t268 + t402;
t198 = rSges(4,1) * t257 + rSges(4,2) * t256 + rSges(4,3) * t395;
t199 = rSges(4,1) * t259 + rSges(4,2) * t258 + rSges(4,3) * t392;
t312 = t198 * t281 - t199 * t280;
t222 = -Icges(6,6) * t286 + t284 * t322;
t223 = -Icges(6,5) * t286 + t284 * t326;
t311 = t222 * t268 - t223 * t269;
t103 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t352;
t304 = t103 * t284 + t158 * t369;
t104 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t350;
t303 = t104 * t284 + t159 * t369;
t102 = -t284 * t311 - t401;
t155 = (-Icges(6,2) * t269 - t403) * t397 + (Icges(6,6) * t284 + t286 * t322) * qJD(2);
t156 = (-Icges(6,1) * t268 - t402) * t397 + (Icges(6,5) * t284 + t286 * t326) * qJD(2);
t105 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t352;
t107 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t352;
t24 = (qJD(2) * t318 - t103) * t286 + (qJD(2) * t158 + (-t160 * t278 + t107) * t269 + (-t162 * t278 - t105) * t268) * t284;
t106 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t350;
t108 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t350;
t25 = (qJD(2) * t317 - t104) * t286 + (qJD(2) * t159 + (-t161 * t278 + t108) * t269 + (-t163 * t278 - t106) * t268) * t284;
t27 = t230 * t105 + t231 * t107 + t192 * t160 + t193 * t162 + t280 * t304;
t28 = t230 * t106 + t231 * t108 + t192 * t161 + t193 * t163 + t280 * t303;
t62 = t158 * t395 + t160 * t230 + t162 * t231;
t88 = t221 * t395 + t222 * t230 + t223 * t231;
t5 = -(t230 * t155 + t231 * t156 + t192 * t222 + t193 * t223) * t286 + (t28 * t281 + (t27 - t388) * t280) * t284 + (t88 * t284 + (t413 + (t62 - t401) * t280) * t286) * qJD(2);
t29 = t232 * t105 + t233 * t107 + t194 * t160 + t195 * t162 + t281 * t304;
t30 = t232 * t106 + t233 * t108 + t194 * t161 + t195 * t163 + t281 * t303;
t65 = t159 * t392 + t161 * t232 + t163 * t233;
t89 = t221 * t392 + t222 * t232 + t223 * t233;
t6 = -(t232 * t155 + t233 * t156 + t194 * t222 + t195 * t223) * t286 + (t29 * t280 + (t30 - t388) * t281) * t284 + (t89 * t284 + (t416 + (t65 - t401) * t281) * t286) * qJD(2);
t296 = -t286 * ((t388 + (t286 * t311 + t332) * qJD(2)) * t286 + (t25 * t281 + t24 * t280 - (qJD(2) * t221 + (-t222 * t278 + t156) * t269 + (-t223 * t278 - t155) * t268) * t286 + t102 * qJD(2)) * t284) + t6 * t392 + t5 * t395 + (-t286 * t88 + (t280 * t62 + t413) * t284) * t352 + (-t286 * t89 + (t281 * t65 + t416) * t284) * t350 + (-t102 * t286 + t284 * t332) * t370;
t15 = -t27 * t281 + t28 * t280;
t16 = t280 * t30 - t281 * t29;
t295 = t5 * t424 + t6 * t425 + (-t24 * t281 + t25 * t280) * t421 + t15 * t395 / 0.2e1 + t16 * t392 / 0.2e1 + (t280 * t74 - t281 * t73) * t370 / 0.2e1 + (t280 * (t280 * t63 - t281 * t62) + t281 * (t280 * t65 - t281 * t64)) * t369 / 0.2e1;
t291 = qJD(2) * (-Icges(3,5) * t284 - Icges(3,6) * t286);
t251 = t281 * t291;
t250 = t280 * t291;
t207 = t374 * t281;
t206 = t374 * t280;
t185 = t436 * t434;
t181 = t197 * t370;
t178 = t196 * t392;
t172 = rSges(5,1) * t241 + rSges(5,2) * t240 + rSges(5,3) * t395;
t153 = t377 * t281;
t152 = t377 * t280;
t150 = t165 * t370;
t149 = t164 * t392;
t146 = t354 * t281;
t145 = t354 * t280;
t141 = t219 * rSges(4,1) + t218 * rSges(4,2) + rSges(4,3) * t350;
t140 = t217 * rSges(4,1) + t216 * rSges(4,2) + rSges(4,3) * t352;
t130 = -t199 * t286 - t247 * t392;
t129 = t198 * t286 + t247 * t395;
t123 = t209 * rSges(5,1) + t208 * rSges(5,2) + rSges(5,3) * t352;
t115 = -t165 * t286 - t224 * t392;
t113 = t244 * t392 + t245 * t258 + t246 * t259;
t111 = t312 * t284;
t99 = t341 * t281;
t98 = t341 * t280;
t96 = t280 * t288 + t281 * t339;
t95 = t356 * t281;
t94 = t356 * t280;
t93 = t226 * t392 + t227 * t242 + t228 * t243;
t91 = -t165 * t395 + t149;
t90 = t198 * t280 + t199 * t281 + t373;
t85 = t376 * t392 + t345;
t84 = t172 * t286 + t229 * t395 + t378;
t77 = t140 * t280 + t141 * t281 - t375;
t76 = -t204 * t392 - t141 * t286 + (t199 * t284 - t247 * t391) * qJD(2);
t75 = t204 * t395 + t140 * t286 + (-t198 * t284 + t247 * t394) * qJD(2);
t72 = t343 * t281;
t71 = t343 * t280;
t66 = t178 + (t172 * t281 + t280 * t380) * t284;
t61 = t172 * t280 + t173 * t281 + t342;
t60 = (t140 * t281 - t141 * t280) * t284 + t312 * t369;
t59 = -t110 * t286 + t150 + (-t157 * t284 - t224 * t369) * t281;
t58 = -t164 * t370 + t360;
t57 = t355 * t392 + t338;
t56 = t131 * t286 + t205 * t395 + t114 + t378;
t55 = (-t110 * t284 - t165 * t369) * t280 + t385;
t54 = t123 * t280 + t124 * t281 + t344;
t53 = t149 + t178 + (t131 * t281 + t280 * t359) * t284;
t52 = t280 * t382 + t281 * t381 + t342;
t50 = t173 * t370 + t181 + t384 * t286 + (t284 * t379 + t369 * t376) * t281;
t49 = t179 * t395 + t123 * t286 + (t229 * t394 + (-t172 - t196) * t284) * qJD(2) + t358;
t48 = t258 * t137 + t259 * t139 + t218 * t189 + t219 * t191 + t281 * t299;
t47 = t258 * t136 + t259 * t138 + t218 * t188 + t219 * t190 + t281 * t300;
t39 = t242 * t120 + t243 * t122 + t210 * t169 + t211 * t171 + t281 * t301;
t38 = t242 * t119 + t243 * t121 + t210 * t168 + t211 * t170 + t281 * t302;
t35 = (t123 * t284 + t172 * t369) * t281 + (qJD(2) * t345 + t284 * t384) * t280 + t383;
t26 = t410 * t281 + (t109 + t96) * t280 + t344;
t23 = t132 * t370 + t150 + t181 + t363 * t286 + (t284 * t357 + t355 * t369) * t281;
t22 = t142 * t395 + t286 * t96 + (t205 * t394 + (-t196 - t382) * t284) * qJD(2) + t358 + t360;
t21 = t280 * t48 - t281 * t47;
t20 = t280 * t46 - t281 * t45;
t19 = (t131 * t369 + t284 * t96) * t281 + (qJD(2) * t338 + t284 * t363) * t280 + t383 + t385;
t18 = t280 * t39 - t281 * t38;
t17 = t280 * t37 - t281 * t36;
t11 = -(t258 * t202 + t259 * t203 + t218 * t245 + t219 * t246) * t286 + (t47 * t280 + (t48 - t386) * t281) * t284 + (t113 * t284 + (t414 + (t83 - t399) * t281) * t286) * qJD(2);
t9 = -(t242 * t176 + t243 * t177 + t210 * t227 + t211 * t228) * t286 + (t38 * t280 + (t39 - t387) * t281) * t284 + (t93 * t284 + (t415 + (t70 - t400) * t281) * t286) * qJD(2);
t1 = [0; -m(3) * t185 + m(4) * t77 + m(5) * t54 + m(6) * t26; (t26 * t52 + t71 * t98 + t72 * t99) * t429 + (t145 * t94 + t146 * t95 + t54 * t61) * t430 + (t152 * t206 + t153 * t207 + t77 * t90) * t431 + 0.2e1 * m(3) * (-t185 + t434) * t436 * (rSges(3,1) * t286 - rSges(3,2) * t284) + (-t277 * t250 - t15 - t17 - t20) * t281 + (t276 * t251 + t16 + t18 + t21 + (-t280 * t250 + t281 * t251) * t281) * t280; m(4) * t60 + m(5) * t35 + m(6) * t19; t295 + m(4) * (t111 * t77 + t129 * t153 + t130 * t152 + t206 * t76 + t207 * t75 + t60 * t90) + m(5) * (t145 * t50 + t146 * t49 + t35 * t61 + t54 * t66 + t84 * t95 + t85 * t94) + m(6) * (t19 * t52 + t22 * t99 + t23 * t98 + t26 * t53 + t56 * t72 + t57 * t71) + (((t457 * t280 - t441 * t281) * t425 + ((-t69 - t82) * t281 + t440 * t280) * t281 / 0.2e1) * t286 + (t446 * t280 - t447 * t281) * t284 / 0.2e1) * qJD(2) + ((t18 / 0.2e1 + t21 / 0.2e1) * t281 + (t17 / 0.2e1 + t20 / 0.2e1) * t280) * t284 + t11 * t425 + t9 * t425 + t444 * t424 + (-t442 * t280 + t443 * t281) * t421; (t19 * t53 + t22 * t56 + t23 * t57) * t429 + (t35 * t66 + t49 * t84 + t50 * t85) * t430 + (t111 * t60 + t129 * t75 + t130 * t76) * t431 + t296 + t444 * t395 + (t9 + t11) * t392 + (t435 * t370 + (t441 * t280 + t448) * t352 + (t440 * t281 + t414 + t415) * t350) * t284 + (t450 * t286 + t437 * t370 - t451 * t352 + (-t113 - t93) * t350 + ((-t176 * t271 + t177 * t272 - t202 * t283 + t203 * t285 - t245 * t366 - t246 * t368 + (-t227 * t272 - t228 * t271) * qJD(3)) * t286 + t442 * t281 + t443 * t280) * t284 + ((-t445 * t286 - t435) * t286 + (t456 * t286 + t437) * t284) * qJD(2)) * t286; (m(5) + m(6)) * t370; m(6) * (-t26 * t286 + t392 * t72 + t395 * t71) + m(5) * (-t286 * t54 + t392 * t95 + t395 * t94) + ((t284 * t52 + t391 * t99 + t394 * t98) * t426 + (t145 * t394 + t146 * t391 + t284 * t61) * t427) * t428; m(6) * (-t19 * t286 + t22 * t392 + t23 * t395) + m(5) * (-t286 * t35 + t392 * t49 + t395 * t50) + ((t284 * t53 + t391 * t56 + t394 * t57) * t426 + (t284 * t66 + t391 * t84 + t394 * t85) * t427) * t428; 0.4e1 * (t427 + t426) * (-0.1e1 + t436) * t284 * t369; m(6) * t55; m(6) * (t114 * t72 + t115 * t71 + t26 * t91 + t52 * t55 + t58 * t99 + t59 * t98) + t295; m(6) * (t114 * t22 + t115 * t23 + t19 * t91 + t53 * t55 + t56 * t58 + t57 * t59) + t296; m(6) * (-t286 * t55 + (t280 * t59 + t281 * t58) * t284 + (t284 * t91 + (t114 * t281 + t115 * t280) * t286) * qJD(2)); (t114 * t58 + t115 * t59 + t55 * t91) * t429 + t296;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
