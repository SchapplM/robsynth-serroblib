% Calculate joint inertia matrix for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:22
% EndTime: 2019-03-10 02:16:43
% DurationCPUTime: 9.29s
% Computational Cost: add. (36005->697), mult. (74203->939), div. (0->0), fcn. (95231->12), ass. (0->330)
t343 = cos(pkin(6));
t345 = sin(qJ(3));
t346 = sin(qJ(2));
t342 = sin(pkin(6));
t423 = cos(qJ(3));
t379 = t342 * t423;
t325 = t343 * t345 + t346 * t379;
t406 = qJ(4) + qJ(5);
t340 = sin(t406);
t373 = cos(t406);
t349 = cos(qJ(2));
t412 = t342 * t349;
t293 = t325 * t340 + t373 * t412;
t294 = t325 * t373 - t340 * t412;
t414 = t342 * t346;
t324 = -t343 * t423 + t345 * t414;
t222 = Icges(7,5) * t294 + Icges(7,6) * t324 + Icges(7,3) * t293;
t224 = Icges(7,4) * t294 + Icges(7,2) * t324 + Icges(7,6) * t293;
t226 = Icges(7,1) * t294 + Icges(7,4) * t324 + Icges(7,5) * t293;
t130 = t293 * t222 + t324 * t224 + t294 * t226;
t223 = Icges(6,5) * t294 - Icges(6,6) * t293 + Icges(6,3) * t324;
t225 = Icges(6,4) * t294 - Icges(6,2) * t293 + Icges(6,6) * t324;
t227 = Icges(6,1) * t294 - Icges(6,4) * t293 + Icges(6,5) * t324;
t131 = t324 * t223 - t293 * t225 + t294 * t227;
t446 = t130 + t131;
t347 = sin(qJ(1));
t408 = t347 * t349;
t350 = cos(qJ(1));
t409 = t346 * t350;
t327 = t343 * t409 + t408;
t411 = t342 * t350;
t306 = t327 * t423 - t345 * t411;
t407 = t349 * t350;
t410 = t346 * t347;
t326 = -t343 * t407 + t410;
t261 = t306 * t340 - t326 * t373;
t262 = t306 * t373 + t326 * t340;
t443 = rSges(7,3) + qJ(6);
t444 = rSges(7,1) + pkin(5);
t445 = -t443 * t261 - t444 * t262;
t305 = t327 * t345 + t350 * t379;
t113 = t222 * t261 + t224 * t305 + t226 * t262;
t114 = t223 * t305 - t225 * t261 + t227 * t262;
t442 = t113 + t114;
t329 = -t343 * t410 + t407;
t413 = t342 * t347;
t308 = t329 * t423 + t345 * t413;
t328 = t343 * t408 + t409;
t263 = t308 * t340 - t328 * t373;
t264 = t308 * t373 + t328 * t340;
t307 = t329 * t345 - t347 * t379;
t115 = t222 * t263 + t224 * t307 + t226 * t264;
t116 = t223 * t307 - t225 * t263 + t227 * t264;
t441 = t115 + t116;
t440 = t446 * t324;
t173 = Icges(7,5) * t262 + Icges(7,6) * t305 + Icges(7,3) * t261;
t177 = Icges(7,4) * t262 + Icges(7,2) * t305 + Icges(7,6) * t261;
t181 = Icges(7,1) * t262 + Icges(7,4) * t305 + Icges(7,5) * t261;
t81 = t173 * t261 + t177 * t305 + t181 * t262;
t174 = Icges(7,5) * t264 + Icges(7,6) * t307 + Icges(7,3) * t263;
t178 = Icges(7,4) * t264 + Icges(7,2) * t307 + Icges(7,6) * t263;
t182 = Icges(7,1) * t264 + Icges(7,4) * t307 + Icges(7,5) * t263;
t82 = t174 * t261 + t178 * t305 + t182 * t262;
t175 = Icges(6,5) * t262 - Icges(6,6) * t261 + Icges(6,3) * t305;
t179 = Icges(6,4) * t262 - Icges(6,2) * t261 + Icges(6,6) * t305;
t183 = Icges(6,1) * t262 - Icges(6,4) * t261 + Icges(6,5) * t305;
t83 = t175 * t305 - t179 * t261 + t183 * t262;
t176 = Icges(6,5) * t264 - Icges(6,6) * t263 + Icges(6,3) * t307;
t180 = Icges(6,4) * t264 - Icges(6,2) * t263 + Icges(6,6) * t307;
t184 = Icges(6,1) * t264 - Icges(6,4) * t263 + Icges(6,5) * t307;
t84 = t176 * t305 - t180 * t261 + t184 * t262;
t439 = t442 * t324 + (t82 + t84) * t307 + (t81 + t83) * t305;
t85 = t173 * t263 + t177 * t307 + t181 * t264;
t86 = t174 * t263 + t178 * t307 + t182 * t264;
t87 = t175 * t307 - t179 * t263 + t183 * t264;
t88 = t176 * t307 - t180 * t263 + t184 * t264;
t438 = t441 * t324 + (t86 + t88) * t307 + (t85 + t87) * t305;
t25 = -t113 * t412 + t326 * t81 + t328 * t82;
t26 = -t114 * t412 + t326 * t83 + t328 * t84;
t437 = t25 + t26;
t27 = -t115 * t412 + t326 * t85 + t328 * t86;
t28 = -t116 * t412 + t326 * t87 + t328 * t88;
t436 = t27 + t28;
t29 = t113 * t343 + (t347 * t82 - t350 * t81) * t342;
t30 = t114 * t343 + (t347 * t84 - t350 * t83) * t342;
t435 = t29 + t30;
t31 = t115 * t343 + (t347 * t86 - t350 * t85) * t342;
t32 = t116 * t343 + (t347 * t88 - t350 * t87) * t342;
t434 = t31 + t32;
t102 = t176 * t324 - t180 * t293 + t184 * t294;
t418 = t102 * t307;
t101 = t175 * t324 - t179 * t293 + t183 * t294;
t419 = t101 * t305;
t100 = t174 * t293 + t178 * t324 + t182 * t294;
t420 = t100 * t307;
t99 = t173 * t293 + t177 * t324 + t181 * t294;
t421 = t99 * t305;
t433 = t420 + t421 + t418 + t419 + t440;
t49 = t100 * t328 - t130 * t412 + t99 * t326;
t50 = t101 * t326 + t102 * t328 - t131 * t412;
t432 = t49 + t50;
t127 = t130 * t343;
t53 = t127 + (t100 * t347 - t99 * t350) * t342;
t128 = t131 * t343;
t54 = t128 + (-t101 * t350 + t102 * t347) * t342;
t431 = t53 + t54;
t401 = t305 * rSges(7,2) - t445;
t400 = t307 * rSges(7,2) + t443 * t263 + t264 * t444;
t397 = rSges(7,2) * t324 + t443 * t293 + t294 * t444;
t344 = sin(qJ(4));
t348 = cos(qJ(4));
t303 = -t325 * t344 - t348 * t412;
t388 = t344 * t412;
t304 = t325 * t348 - t388;
t232 = Icges(5,5) * t304 + Icges(5,6) * t303 + Icges(5,3) * t324;
t233 = Icges(5,4) * t304 + Icges(5,2) * t303 + Icges(5,6) * t324;
t234 = Icges(5,1) * t304 + Icges(5,4) * t303 + Icges(5,5) * t324;
t135 = t324 * t232 + t303 * t233 + t304 * t234;
t274 = Icges(4,5) * t325 - Icges(4,6) * t324 - Icges(4,3) * t412;
t275 = Icges(4,4) * t325 - Icges(4,2) * t324 - Icges(4,6) * t412;
t276 = Icges(4,1) * t325 - Icges(4,4) * t324 - Icges(4,5) * t412;
t160 = -t274 * t412 - t324 * t275 + t325 * t276;
t430 = -t135 - t160 - t446;
t429 = t305 / 0.2e1;
t428 = t307 / 0.2e1;
t427 = t324 / 0.2e1;
t426 = t326 / 0.2e1;
t425 = t328 / 0.2e1;
t424 = t343 / 0.2e1;
t339 = pkin(4) * t348 + pkin(3);
t422 = -pkin(3) + t339;
t278 = Icges(3,5) * t327 - Icges(3,6) * t326 - Icges(3,3) * t411;
t417 = t278 * t350;
t416 = t326 * t344;
t415 = t328 * t344;
t405 = t401 * t307;
t301 = t305 * pkin(10);
t351 = -pkin(11) - pkin(10);
t389 = pkin(4) * t416;
t169 = -t305 * t351 + t306 * t422 - t301 + t389;
t252 = t306 * pkin(3) + t301;
t244 = t328 * t252;
t404 = t328 * t169 + t244;
t403 = t400 * t324;
t253 = t308 * pkin(3) + pkin(10) * t307;
t381 = pkin(4) * t415 - t307 * t351 + t308 * t339;
t170 = -t253 + t381;
t190 = t264 * rSges(6,1) - t263 * rSges(6,2) + t307 * rSges(6,3);
t402 = -t170 - t190;
t272 = -t308 * t344 + t328 * t348;
t273 = t308 * t348 + t415;
t199 = t273 * rSges(5,1) + t272 * rSges(5,2) + t307 * rSges(5,3);
t399 = -t199 - t253;
t398 = t397 * t305;
t229 = rSges(6,1) * t294 - rSges(6,2) * t293 + rSges(6,3) * t324;
t230 = -pkin(4) * t388 + t422 * t325 + (-pkin(10) - t351) * t324;
t396 = -t229 - t230;
t235 = rSges(5,1) * t304 + rSges(5,2) * t303 + rSges(5,3) * t324;
t289 = t325 * pkin(3) + t324 * pkin(10);
t395 = -t235 - t289;
t394 = t252 * t412 + t326 * t289;
t291 = t329 * pkin(2) + pkin(9) * t328;
t288 = t343 * t291;
t393 = t343 * t253 + t288;
t290 = t327 * pkin(2) + t326 * pkin(9);
t392 = -t252 - t290;
t391 = t290 * t413 + t291 * t411;
t390 = t350 * pkin(1) + pkin(8) * t413;
t387 = t343 * t170 + t393;
t386 = -t169 + t392;
t385 = -t170 - t400;
t384 = -t253 + t402;
t383 = -t230 - t397;
t382 = -t289 + t396;
t243 = t308 * rSges(4,1) - t307 * rSges(4,2) + t328 * rSges(4,3);
t314 = Icges(3,3) * t343 + (Icges(3,5) * t346 + Icges(3,6) * t349) * t342;
t315 = Icges(3,6) * t343 + (Icges(3,4) * t346 + Icges(3,2) * t349) * t342;
t316 = Icges(3,5) * t343 + (Icges(3,1) * t346 + Icges(3,4) * t349) * t342;
t380 = t343 * t314 + t315 * t412 + t316 * t414;
t285 = t329 * rSges(3,1) - t328 * rSges(3,2) + rSges(3,3) * t413;
t377 = -t412 / 0.2e1;
t270 = -t306 * t344 + t326 * t348;
t271 = t306 * t348 + t416;
t192 = Icges(5,5) * t271 + Icges(5,6) * t270 + Icges(5,3) * t305;
t194 = Icges(5,4) * t271 + Icges(5,2) * t270 + Icges(5,6) * t305;
t196 = Icges(5,1) * t271 + Icges(5,4) * t270 + Icges(5,5) * t305;
t106 = t192 * t324 + t194 * t303 + t196 * t304;
t121 = t232 * t305 + t233 * t270 + t234 * t271;
t375 = t121 / 0.2e1 + t106 / 0.2e1;
t193 = Icges(5,5) * t273 + Icges(5,6) * t272 + Icges(5,3) * t307;
t195 = Icges(5,4) * t273 + Icges(5,2) * t272 + Icges(5,6) * t307;
t197 = Icges(5,1) * t273 + Icges(5,4) * t272 + Icges(5,5) * t307;
t107 = t193 * t324 + t195 * t303 + t197 * t304;
t122 = t232 * t307 + t233 * t272 + t234 * t273;
t374 = t122 / 0.2e1 + t107 / 0.2e1;
t372 = -t347 * pkin(1) + pkin(8) * t411;
t277 = rSges(4,1) * t325 - rSges(4,2) * t324 - rSges(4,3) * t412;
t330 = (pkin(2) * t346 - pkin(9) * t349) * t342;
t371 = t342 * (-t277 - t330);
t370 = t169 * t412 + t326 * t230 + t394;
t369 = -t253 + t385;
t368 = -t289 + t383;
t367 = t252 * t413 + t253 * t411 + t391;
t366 = t342 * (-t330 + t395);
t365 = t305 * t439 + t438 * t307 + t433 * t324;
t364 = -t262 * rSges(6,1) + t261 * rSges(6,2);
t363 = t291 + t390;
t362 = t342 * (-t330 + t382);
t361 = t169 * t413 + t170 * t411 + t367;
t360 = t342 * (-t330 + t368);
t359 = -t290 + t372;
t242 = rSges(4,1) * t306 - rSges(4,2) * t305 + rSges(4,3) * t326;
t198 = rSges(5,1) * t271 + rSges(5,2) * t270 + rSges(5,3) * t305;
t284 = rSges(3,1) * t327 - rSges(3,2) * t326 - rSges(3,3) * t411;
t358 = t363 + t381;
t357 = t421 / 0.2e1 + t420 / 0.2e1 + t419 / 0.2e1 + t418 / 0.2e1 + t442 * t429 + t441 * t428 + t440;
t356 = t433 * t377 + t438 * t425 + t426 * t439 + t432 * t427 + t436 * t428 + t437 * t429;
t355 = t435 * t429 + t434 * t428 + t431 * t427 + t433 * t424 + t438 * t413 / 0.2e1 - t439 * t411 / 0.2e1;
t236 = Icges(4,5) * t306 - Icges(4,6) * t305 + Icges(4,3) * t326;
t238 = Icges(4,4) * t306 - Icges(4,2) * t305 + Icges(4,6) * t326;
t240 = Icges(4,1) * t306 - Icges(4,4) * t305 + Icges(4,5) * t326;
t142 = -t236 * t412 - t238 * t324 + t240 * t325;
t151 = t274 * t326 - t275 * t305 + t276 * t306;
t354 = t99 / 0.2e1 + t113 / 0.2e1 + t114 / 0.2e1 + t101 / 0.2e1 + t142 / 0.2e1 + t151 / 0.2e1 + t375;
t353 = -t306 * t339 + t359 - t389;
t237 = Icges(4,5) * t308 - Icges(4,6) * t307 + Icges(4,3) * t328;
t239 = Icges(4,4) * t308 - Icges(4,2) * t307 + Icges(4,6) * t328;
t241 = Icges(4,1) * t308 - Icges(4,4) * t307 + Icges(4,5) * t328;
t143 = -t237 * t412 - t239 * t324 + t241 * t325;
t152 = t274 * t328 - t275 * t307 + t276 * t308;
t352 = t102 / 0.2e1 + t100 / 0.2e1 + t115 / 0.2e1 + t116 / 0.2e1 + t152 / 0.2e1 + t143 / 0.2e1 + t374;
t332 = rSges(2,1) * t350 - rSges(2,2) * t347;
t331 = -rSges(2,1) * t347 - rSges(2,2) * t350;
t317 = t343 * rSges(3,3) + (rSges(3,1) * t346 + rSges(3,2) * t349) * t342;
t283 = Icges(3,1) * t329 - Icges(3,4) * t328 + Icges(3,5) * t413;
t282 = Icges(3,1) * t327 - Icges(3,4) * t326 - Icges(3,5) * t411;
t281 = Icges(3,4) * t329 - Icges(3,2) * t328 + Icges(3,6) * t413;
t280 = Icges(3,4) * t327 - Icges(3,2) * t326 - Icges(3,6) * t411;
t279 = Icges(3,5) * t329 - Icges(3,6) * t328 + Icges(3,3) * t413;
t266 = t285 + t390;
t265 = -t284 + t372;
t246 = -t284 * t343 - t317 * t411;
t245 = t285 * t343 - t317 * t413;
t231 = t380 * t343;
t220 = (t284 * t347 + t285 * t350) * t342;
t218 = t314 * t413 - t315 * t328 + t316 * t329;
t217 = -t314 * t411 - t315 * t326 + t316 * t327;
t209 = t305 * t230;
t208 = t305 * t229;
t202 = t363 + t243;
t201 = -t242 + t359;
t188 = t305 * rSges(6,3) - t364;
t186 = -t243 * t412 - t277 * t328;
t185 = t242 * t412 + t277 * t326;
t172 = t343 * t279 + (t281 * t349 + t283 * t346) * t342;
t171 = t343 * t278 + (t280 * t349 + t282 * t346) * t342;
t164 = t324 * t190;
t161 = t324 * t170;
t159 = t307 * t188;
t157 = t160 * t343;
t156 = t307 * t169;
t155 = t242 * t328 - t243 * t326;
t154 = (-t242 - t290) * t343 + t350 * t371;
t153 = t343 * t243 + t347 * t371 + t288;
t150 = t363 - t399;
t149 = -t198 - t252 + t359;
t148 = (t242 * t347 + t243 * t350) * t342 + t391;
t147 = t199 * t324 - t235 * t307;
t146 = -t198 * t324 + t235 * t305;
t145 = t358 + t190;
t144 = (-rSges(6,3) + t351) * t305 + t353 + t364;
t141 = -t229 * t307 + t164;
t140 = -t188 * t324 + t208;
t139 = t237 * t328 - t239 * t307 + t241 * t308;
t138 = t236 * t328 - t238 * t307 + t240 * t308;
t137 = t237 * t326 - t239 * t305 + t241 * t306;
t136 = t236 * t326 - t238 * t305 + t240 * t306;
t134 = t135 * t343;
t133 = t135 * t324;
t132 = t198 * t307 - t199 * t305;
t129 = -t190 * t305 + t159;
t126 = t328 * t395 + t399 * t412;
t125 = t198 * t412 + t235 * t326 + t394;
t120 = t358 + t400;
t119 = (-rSges(7,2) + t351) * t305 + t353 + t445;
t118 = (-t198 + t392) * t343 + t350 * t366;
t117 = t343 * t199 + t347 * t366 + t393;
t108 = t328 * t198 + t326 * t399 + t244;
t105 = (t198 * t347 + t199 * t350) * t342 + t367;
t104 = -t307 * t397 + t403;
t103 = -t324 * t401 + t398;
t98 = t193 * t307 + t195 * t272 + t197 * t273;
t97 = t192 * t307 + t194 * t272 + t196 * t273;
t96 = t193 * t305 + t195 * t270 + t197 * t271;
t95 = t192 * t305 + t194 * t270 + t196 * t271;
t90 = t307 * t396 + t161 + t164;
t89 = t208 + t209 + (-t169 - t188) * t324;
t80 = t328 * t382 + t384 * t412;
t79 = t188 * t412 + t229 * t326 + t370;
t78 = -t305 * t400 + t405;
t77 = (-t188 + t386) * t343 + t350 * t362;
t76 = t343 * t190 + t347 * t362 + t387;
t75 = t305 * t402 + t156 + t159;
t74 = t157 + (-t142 * t350 + t143 * t347) * t342;
t73 = t142 * t326 + t143 * t328 - t160 * t412;
t72 = t307 * t383 + t161 + t403;
t71 = t209 + (-t169 - t401) * t324 + t398;
t70 = t328 * t188 + t326 * t384 + t404;
t69 = (t188 * t347 + t190 * t350) * t342 + t361;
t68 = t328 * t368 + t369 * t412;
t67 = t326 * t397 + t401 * t412 + t370;
t66 = (t386 - t401) * t343 + t350 * t360;
t65 = t343 * t400 + t347 * t360 + t387;
t64 = t152 * t343 + (-t138 * t350 + t139 * t347) * t342;
t63 = t151 * t343 + (-t136 * t350 + t137 * t347) * t342;
t62 = t138 * t326 + t139 * t328 - t152 * t412;
t61 = t136 * t326 + t137 * t328 - t151 * t412;
t60 = t305 * t385 + t156 + t405;
t59 = t326 * t369 + t328 * t401 + t404;
t58 = (t347 * t401 + t350 * t400) * t342 + t361;
t57 = t134 + (-t106 * t350 + t107 * t347) * t342;
t56 = t106 * t326 + t107 * t328 - t135 * t412;
t55 = t106 * t305 + t107 * t307 + t133;
t38 = t122 * t343 + (t347 * t98 - t350 * t97) * t342;
t37 = t121 * t343 + (t347 * t96 - t350 * t95) * t342;
t36 = -t122 * t412 + t326 * t97 + t328 * t98;
t35 = -t121 * t412 + t326 * t95 + t328 * t96;
t34 = t122 * t324 + t305 * t97 + t307 * t98;
t33 = t121 * t324 + t305 * t95 + t307 * t96;
t1 = [m(6) * (t144 ^ 2 + t145 ^ 2) + m(7) * (t119 ^ 2 + t120 ^ 2) + m(5) * (t149 ^ 2 + t150 ^ 2) + m(4) * (t201 ^ 2 + t202 ^ 2) + m(3) * (t265 ^ 2 + t266 ^ 2) + m(2) * (t331 ^ 2 + t332 ^ 2) + Icges(2,3) + t380 - t430; t134 + t128 + t231 + t127 + t157 + m(7) * (t119 * t66 + t120 * t65) + m(6) * (t144 * t77 + t145 * t76) + m(5) * (t117 * t150 + t118 * t149) + m(4) * (t153 * t202 + t154 * t201) + m(3) * (t245 * t266 + t246 * t265) + ((-t171 / 0.2e1 - t217 / 0.2e1 - t354) * t350 + (t172 / 0.2e1 + t218 / 0.2e1 + t352) * t347) * t342; (t57 + t74 + t231 + t431) * t343 + m(7) * (t58 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t69 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t105 ^ 2 + t117 ^ 2 + t118 ^ 2) + m(4) * (t148 ^ 2 + t153 ^ 2 + t154 ^ 2) + m(3) * (t220 ^ 2 + t245 ^ 2 + t246 ^ 2) + ((-t37 - t63 + (-t280 * t326 + t282 * t327 - t342 * t417) * t411 - t435) * t350 + (t38 + t64 + ((-t281 * t328 + t283 * t329 + (t279 * t347 - t417) * t342) * t347 + (t279 * t411 + t280 * t328 + t281 * t326 - t282 * t329 - t283 * t327) * t350) * t342 + t434) * t347 + ((-t171 - t217) * t350 + (t172 + t218) * t347) * t343) * t342; t430 * t412 + m(6) * (t144 * t79 + t145 * t80) + m(7) * (t119 * t67 + t120 * t68) + m(5) * (t125 * t149 + t126 * t150) + m(4) * (t185 * t201 + t186 * t202) + t352 * t328 + t354 * t326; (t49 / 0.2e1 + t50 / 0.2e1 + t56 / 0.2e1 + t73 / 0.2e1) * t343 + (t31 / 0.2e1 + t32 / 0.2e1 + t38 / 0.2e1 + t64 / 0.2e1) * t328 + (t29 / 0.2e1 + t30 / 0.2e1 + t37 / 0.2e1 + t63 / 0.2e1) * t326 + m(7) * (t59 * t58 + t65 * t68 + t66 * t67) + m(6) * (t69 * t70 + t76 * t80 + t77 * t79) + m(5) * (t105 * t108 + t117 * t126 + t118 * t125) + m(4) * (t148 * t155 + t153 * t186 + t154 * t185) + ((-t35 / 0.2e1 - t25 / 0.2e1 - t26 / 0.2e1 - t61 / 0.2e1) * t350 + (-t53 / 0.2e1 - t54 / 0.2e1 - t57 / 0.2e1 - t74 / 0.2e1) * t349 + (t36 / 0.2e1 + t27 / 0.2e1 + t28 / 0.2e1 + t62 / 0.2e1) * t347) * t342; (-t56 - t73 - t432) * t412 + (t36 + t62 + t436) * t328 + (t35 + t61 + t437) * t326 + m(7) * (t59 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t70 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t108 ^ 2 + t125 ^ 2 + t126 ^ 2) + m(4) * (t155 ^ 2 + t185 ^ 2 + t186 ^ 2); t357 + t133 + t374 * t307 + t375 * t305 + m(6) * (t144 * t89 + t145 * t90) + m(7) * (t119 * t71 + t120 * t72) + m(5) * (t146 * t149 + t147 * t150); m(7) * (t60 * t58 + t65 * t72 + t66 * t71) + m(6) * (t69 * t75 + t76 * t90 + t77 * t89) + m(5) * (t105 * t132 + t117 * t147 + t118 * t146) + (t347 * t34 / 0.2e1 - t350 * t33 / 0.2e1) * t342 + t37 * t429 + t38 * t428 + t57 * t427 + t55 * t424 + t355; t55 * t377 + m(7) * (t60 * t59 + t67 * t71 + t68 * t72) + m(6) * (t70 * t75 + t79 * t89 + t80 * t90) + m(5) * (t108 * t132 + t125 * t146 + t126 * t147) + t356 + t35 * t429 + t36 * t428 + t56 * t427 + t33 * t426 + t34 * t425; t305 * t33 + t307 * t34 + t324 * t55 + m(7) * (t60 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(5) * (t132 ^ 2 + t146 ^ 2 + t147 ^ 2) + t365; m(6) * (t140 * t144 + t141 * t145) + m(7) * (t103 * t119 + t104 * t120) + t357; t355 + m(7) * (t103 * t66 + t104 * t65 + t58 * t78) + m(6) * (t129 * t69 + t140 * t77 + t141 * t76); t356 + m(7) * (t103 * t67 + t104 * t68 + t59 * t78) + m(6) * (t129 * t70 + t140 * t79 + t141 * t80); m(7) * (t103 * t71 + t104 * t72 + t60 * t78) + m(6) * (t129 * t75 + t140 * t89 + t141 * t90) + t365; m(7) * (t103 ^ 2 + t104 ^ 2 + t78 ^ 2) + m(6) * (t129 ^ 2 + t140 ^ 2 + t141 ^ 2) + t365; m(7) * (t119 * t263 + t120 * t261); m(7) * (t261 * t65 + t263 * t66 + t293 * t58); m(7) * (t261 * t68 + t263 * t67 + t293 * t59); m(7) * (t261 * t72 + t263 * t71 + t293 * t60); m(7) * (t103 * t263 + t104 * t261 + t293 * t78); m(7) * (t261 ^ 2 + t263 ^ 2 + t293 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
