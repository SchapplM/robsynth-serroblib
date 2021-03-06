% Calculate joint inertia matrix for
% S6RRRRRP9
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:24
% EndTime: 2019-03-10 02:00:47
% DurationCPUTime: 9.46s
% Computational Cost: add. (37234->708), mult. (76024->950), div. (0->0), fcn. (97103->12), ass. (0->334)
t357 = -pkin(11) - pkin(10);
t454 = rSges(7,3) + qJ(6) - t357;
t349 = cos(pkin(6));
t351 = sin(qJ(3));
t352 = sin(qJ(2));
t348 = sin(pkin(6));
t431 = cos(qJ(3));
t384 = t348 * t431;
t325 = t349 * t351 + t352 * t384;
t347 = qJ(4) + qJ(5);
t342 = sin(t347);
t343 = cos(t347);
t355 = cos(qJ(2));
t419 = t348 * t355;
t294 = -t325 * t342 - t343 * t419;
t295 = t325 * t343 - t342 * t419;
t421 = t348 * t352;
t324 = -t349 * t431 + t351 * t421;
t223 = Icges(7,5) * t295 + Icges(7,6) * t294 + Icges(7,3) * t324;
t225 = Icges(7,4) * t295 + Icges(7,2) * t294 + Icges(7,6) * t324;
t227 = Icges(7,1) * t295 + Icges(7,4) * t294 + Icges(7,5) * t324;
t128 = t324 * t223 + t294 * t225 + t295 * t227;
t224 = Icges(6,5) * t295 + Icges(6,6) * t294 + Icges(6,3) * t324;
t226 = Icges(6,4) * t295 + Icges(6,2) * t294 + Icges(6,6) * t324;
t228 = Icges(6,1) * t295 + Icges(6,4) * t294 + Icges(6,5) * t324;
t129 = t324 * t224 + t294 * t226 + t295 * t228;
t453 = t128 + t129;
t452 = t357 + t454;
t353 = sin(qJ(1));
t415 = t353 * t355;
t356 = cos(qJ(1));
t416 = t352 * t356;
t327 = t349 * t416 + t415;
t418 = t348 * t356;
t307 = t327 * t431 - t351 * t418;
t414 = t355 * t356;
t417 = t352 * t353;
t326 = -t349 * t414 + t417;
t260 = -t307 * t342 + t326 * t343;
t261 = t307 * t343 + t326 * t342;
t306 = t327 * t351 + t356 * t384;
t113 = t223 * t306 + t225 * t260 + t227 * t261;
t114 = t224 * t306 + t226 * t260 + t228 * t261;
t451 = t113 + t114;
t329 = -t349 * t417 + t414;
t420 = t348 * t353;
t309 = t329 * t431 + t351 * t420;
t328 = t349 * t415 + t416;
t262 = -t309 * t342 + t328 * t343;
t263 = t309 * t343 + t328 * t342;
t308 = t329 * t351 - t353 * t384;
t115 = t223 * t308 + t225 * t262 + t227 * t263;
t116 = t224 * t308 + t226 * t262 + t228 * t263;
t450 = t115 + t116;
t449 = t453 * t324;
t354 = cos(qJ(4));
t341 = t354 * pkin(4) + pkin(3);
t331 = pkin(5) * t343 + t341;
t350 = sin(qJ(4));
t430 = pkin(4) * t350;
t332 = pkin(5) * t342 + t430;
t448 = t263 * rSges(7,1) + t262 * rSges(7,2) + t308 * t454 + t309 * t331 + t328 * t332;
t178 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t306;
t182 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t306;
t186 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t306;
t83 = t178 * t306 + t182 * t260 + t186 * t261;
t179 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t308;
t183 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t308;
t187 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t308;
t84 = t179 * t306 + t183 * t260 + t187 * t261;
t180 = Icges(6,5) * t261 + Icges(6,6) * t260 + Icges(6,3) * t306;
t184 = Icges(6,4) * t261 + Icges(6,2) * t260 + Icges(6,6) * t306;
t188 = Icges(6,1) * t261 + Icges(6,4) * t260 + Icges(6,5) * t306;
t85 = t180 * t306 + t184 * t260 + t188 * t261;
t181 = Icges(6,5) * t263 + Icges(6,6) * t262 + Icges(6,3) * t308;
t185 = Icges(6,4) * t263 + Icges(6,2) * t262 + Icges(6,6) * t308;
t189 = Icges(6,1) * t263 + Icges(6,4) * t262 + Icges(6,5) * t308;
t86 = t181 * t306 + t185 * t260 + t189 * t261;
t447 = t451 * t324 + (t84 + t86) * t308 + (t83 + t85) * t306;
t87 = t178 * t308 + t182 * t262 + t186 * t263;
t88 = t179 * t308 + t183 * t262 + t187 * t263;
t89 = t180 * t308 + t184 * t262 + t188 * t263;
t90 = t181 * t308 + t185 * t262 + t189 * t263;
t446 = t450 * t324 + (t88 + t90) * t308 + (t87 + t89) * t306;
t25 = -t113 * t419 + t326 * t83 + t328 * t84;
t26 = -t114 * t419 + t326 * t85 + t328 * t86;
t445 = t25 + t26;
t27 = -t115 * t419 + t326 * t87 + t328 * t88;
t28 = -t116 * t419 + t326 * t89 + t328 * t90;
t444 = t27 + t28;
t29 = t113 * t349 + (t353 * t84 - t356 * t83) * t348;
t30 = t114 * t349 + (t353 * t86 - t356 * t85) * t348;
t443 = t29 + t30;
t31 = t115 * t349 + (t353 * t88 - t356 * t87) * t348;
t32 = t116 * t349 + (t353 * t90 - t356 * t89) * t348;
t442 = t31 + t32;
t104 = t181 * t324 + t185 * t294 + t189 * t295;
t425 = t104 * t308;
t103 = t180 * t324 + t184 * t294 + t188 * t295;
t426 = t103 * t306;
t102 = t179 * t324 + t183 * t294 + t187 * t295;
t427 = t102 * t308;
t101 = t178 * t324 + t182 * t294 + t186 * t295;
t428 = t101 * t306;
t441 = t427 + t428 + t425 + t426 + t449;
t49 = t101 * t326 + t102 * t328 - t128 * t419;
t50 = t103 * t326 + t104 * t328 - t129 * t419;
t440 = t49 + t50;
t125 = t128 * t349;
t53 = t125 + (-t101 * t356 + t102 * t353) * t348;
t126 = t129 * t349;
t54 = t126 + (-t103 * t356 + t104 * t353) * t348;
t439 = t53 + t54;
t368 = -rSges(7,1) * t261 - rSges(7,2) * t260;
t377 = -t332 + t430;
t398 = t331 - t341;
t411 = t306 * t452 + t307 * t398 - t326 * t377 - t368;
t422 = t328 * t350;
t386 = -pkin(4) * t422 + t308 * t357 - t309 * t341;
t410 = t386 + t448;
t405 = rSges(7,1) * t295 + rSges(7,2) * t294 + t324 * t452 + t325 * t398 + t377 * t419;
t304 = -t325 * t350 - t354 * t419;
t394 = t350 * t419;
t305 = t325 * t354 - t394;
t233 = Icges(5,5) * t305 + Icges(5,6) * t304 + Icges(5,3) * t324;
t234 = Icges(5,4) * t305 + Icges(5,2) * t304 + Icges(5,6) * t324;
t235 = Icges(5,1) * t305 + Icges(5,4) * t304 + Icges(5,5) * t324;
t133 = t324 * t233 + t304 * t234 + t305 * t235;
t273 = Icges(4,5) * t325 - Icges(4,6) * t324 - Icges(4,3) * t419;
t274 = Icges(4,4) * t325 - Icges(4,2) * t324 - Icges(4,6) * t419;
t275 = Icges(4,1) * t325 - Icges(4,4) * t324 - Icges(4,5) * t419;
t164 = -t273 * t419 - t324 * t274 + t325 * t275;
t438 = -t133 - t164 - t453;
t437 = t306 / 0.2e1;
t436 = t308 / 0.2e1;
t435 = t324 / 0.2e1;
t434 = t326 / 0.2e1;
t433 = t328 / 0.2e1;
t432 = t349 / 0.2e1;
t429 = -pkin(3) + t341;
t277 = Icges(3,5) * t327 - Icges(3,6) * t326 - Icges(3,3) * t418;
t424 = t277 * t356;
t423 = t326 * t350;
t413 = t411 * t308;
t412 = t410 * t324;
t302 = t306 * pkin(10);
t395 = pkin(4) * t423;
t173 = -t306 * t357 + t307 * t429 - t302 + t395;
t252 = t307 * pkin(3) + t302;
t245 = t328 * t252;
t409 = t328 * t173 + t245;
t253 = t309 * pkin(3) + pkin(10) * t308;
t174 = -t253 - t386;
t195 = t263 * rSges(6,1) + t262 * rSges(6,2) + t308 * rSges(6,3);
t408 = -t174 - t195;
t407 = t405 * t306;
t271 = -t309 * t350 + t328 * t354;
t272 = t309 * t354 + t422;
t203 = t272 * rSges(5,1) + t271 * rSges(5,2) + t308 * rSges(5,3);
t406 = -t203 - t253;
t230 = rSges(6,1) * t295 + rSges(6,2) * t294 + rSges(6,3) * t324;
t231 = -pkin(4) * t394 + t429 * t325 + (-pkin(10) - t357) * t324;
t404 = -t230 - t231;
t236 = rSges(5,1) * t305 + rSges(5,2) * t304 + rSges(5,3) * t324;
t289 = t325 * pkin(3) + t324 * pkin(10);
t403 = -t236 - t289;
t402 = t252 * t419 + t326 * t289;
t291 = t329 * pkin(2) + pkin(9) * t328;
t288 = t349 * t291;
t401 = t349 * t253 + t288;
t290 = t327 * pkin(2) + t326 * pkin(9);
t400 = -t252 - t290;
t399 = t290 * t420 + t291 * t418;
t397 = t356 * pkin(1) + pkin(8) * t420;
t393 = -t174 - t410;
t392 = t349 * t174 + t401;
t391 = -t173 + t400;
t390 = -t253 + t408;
t389 = -t231 - t405;
t388 = -t289 + t404;
t244 = t309 * rSges(4,1) - t308 * rSges(4,2) + t328 * rSges(4,3);
t314 = Icges(3,3) * t349 + (Icges(3,5) * t352 + Icges(3,6) * t355) * t348;
t315 = Icges(3,6) * t349 + (Icges(3,4) * t352 + Icges(3,2) * t355) * t348;
t316 = Icges(3,5) * t349 + (Icges(3,1) * t352 + Icges(3,4) * t355) * t348;
t385 = t349 * t314 + t315 * t419 + t316 * t421;
t284 = t329 * rSges(3,1) - t328 * rSges(3,2) + rSges(3,3) * t420;
t382 = -t419 / 0.2e1;
t197 = Icges(5,5) * t272 + Icges(5,6) * t271 + Icges(5,3) * t308;
t199 = Icges(5,4) * t272 + Icges(5,2) * t271 + Icges(5,6) * t308;
t201 = Icges(5,1) * t272 + Icges(5,4) * t271 + Icges(5,5) * t308;
t107 = t197 * t324 + t199 * t304 + t201 * t305;
t120 = t233 * t308 + t234 * t271 + t235 * t272;
t380 = t107 / 0.2e1 + t120 / 0.2e1;
t269 = -t307 * t350 + t326 * t354;
t270 = t307 * t354 + t423;
t196 = Icges(5,5) * t270 + Icges(5,6) * t269 + Icges(5,3) * t306;
t198 = Icges(5,4) * t270 + Icges(5,2) * t269 + Icges(5,6) * t306;
t200 = Icges(5,1) * t270 + Icges(5,4) * t269 + Icges(5,5) * t306;
t106 = t196 * t324 + t198 * t304 + t200 * t305;
t119 = t233 * t306 + t234 * t269 + t235 * t270;
t379 = t119 / 0.2e1 + t106 / 0.2e1;
t378 = -t353 * pkin(1) + pkin(8) * t418;
t276 = rSges(4,1) * t325 - rSges(4,2) * t324 - rSges(4,3) * t419;
t330 = (pkin(2) * t352 - pkin(9) * t355) * t348;
t376 = t348 * (-t276 - t330);
t375 = -t253 + t393;
t374 = t173 * t419 + t326 * t231 + t402;
t373 = -t289 + t389;
t372 = t252 * t420 + t253 * t418 + t399;
t371 = t348 * (-t330 + t403);
t370 = t306 * t447 + t446 * t308 + t441 * t324;
t369 = -t261 * rSges(6,1) - t260 * rSges(6,2);
t367 = t291 + t397;
t366 = t348 * (-t330 + t388);
t365 = t173 * t420 + t174 * t418 + t372;
t364 = t348 * (-t330 + t373);
t363 = -t290 + t378;
t243 = rSges(4,1) * t307 - rSges(4,2) * t306 + rSges(4,3) * t326;
t202 = rSges(5,1) * t270 + rSges(5,2) * t269 + rSges(5,3) * t306;
t283 = rSges(3,1) * t327 - rSges(3,2) * t326 - rSges(3,3) * t418;
t362 = t428 / 0.2e1 + t427 / 0.2e1 + t426 / 0.2e1 + t425 / 0.2e1 + t451 * t437 + t450 * t436 + t449;
t361 = t441 * t382 + t446 * t433 + t434 * t447 + t440 * t435 + t444 * t436 + t445 * t437;
t360 = t443 * t437 + t442 * t436 + t439 * t435 + t441 * t432 + t446 * t420 / 0.2e1 - t447 * t418 / 0.2e1;
t237 = Icges(4,5) * t307 - Icges(4,6) * t306 + Icges(4,3) * t326;
t239 = Icges(4,4) * t307 - Icges(4,2) * t306 + Icges(4,6) * t326;
t241 = Icges(4,1) * t307 - Icges(4,4) * t306 + Icges(4,5) * t326;
t142 = -t237 * t419 - t239 * t324 + t241 * t325;
t152 = t273 * t326 - t274 * t306 + t275 * t307;
t359 = t101 / 0.2e1 + t114 / 0.2e1 + t113 / 0.2e1 + t103 / 0.2e1 + t152 / 0.2e1 + t142 / 0.2e1 + t379;
t238 = Icges(4,5) * t309 - Icges(4,6) * t308 + Icges(4,3) * t328;
t240 = Icges(4,4) * t309 - Icges(4,2) * t308 + Icges(4,6) * t328;
t242 = Icges(4,1) * t309 - Icges(4,4) * t308 + Icges(4,5) * t328;
t143 = -t238 * t419 - t240 * t324 + t242 * t325;
t153 = t273 * t328 - t274 * t308 + t275 * t309;
t358 = t104 / 0.2e1 + t102 / 0.2e1 + t116 / 0.2e1 + t115 / 0.2e1 + t153 / 0.2e1 + t143 / 0.2e1 + t380;
t334 = rSges(2,1) * t356 - rSges(2,2) * t353;
t333 = -rSges(2,1) * t353 - rSges(2,2) * t356;
t317 = rSges(3,3) * t349 + (rSges(3,1) * t352 + rSges(3,2) * t355) * t348;
t282 = Icges(3,1) * t329 - Icges(3,4) * t328 + Icges(3,5) * t420;
t281 = Icges(3,1) * t327 - Icges(3,4) * t326 - Icges(3,5) * t418;
t280 = Icges(3,4) * t329 - Icges(3,2) * t328 + Icges(3,6) * t420;
t279 = Icges(3,4) * t327 - Icges(3,2) * t326 - Icges(3,6) * t418;
t278 = Icges(3,5) * t329 - Icges(3,6) * t328 + Icges(3,3) * t420;
t265 = t284 + t397;
t264 = -t283 + t378;
t247 = -t283 * t349 - t317 * t418;
t246 = t284 * t349 - t317 * t420;
t232 = t385 * t349;
t221 = (t283 * t353 + t284 * t356) * t348;
t220 = t314 * t420 - t315 * t328 + t316 * t329;
t219 = -t314 * t418 - t315 * t326 + t316 * t327;
t213 = t306 * t231;
t212 = t306 * t230;
t205 = t367 + t244;
t204 = -t243 + t363;
t193 = t306 * rSges(6,3) - t369;
t191 = -t244 * t419 - t276 * t328;
t190 = t243 * t419 + t276 * t326;
t177 = t278 * t349 + (t280 * t355 + t282 * t352) * t348;
t176 = t277 * t349 + (t279 * t355 + t281 * t352) * t348;
t168 = t324 * t195;
t165 = t324 * t174;
t163 = t308 * t193;
t161 = t164 * t349;
t160 = t308 * t173;
t157 = t243 * t328 - t244 * t326;
t156 = (-t243 - t290) * t349 + t356 * t376;
t155 = t244 * t349 + t353 * t376 + t288;
t150 = t367 - t406;
t149 = -t202 - t252 + t363;
t148 = (t243 * t353 + t244 * t356) * t348 + t399;
t147 = t203 * t324 - t308 * t236;
t146 = -t202 * t324 + t236 * t306;
t145 = t367 - t386 + t195;
t144 = -t395 - t307 * t341 + (-rSges(6,3) + t357) * t306 + t363 + t369;
t141 = -t230 * t308 + t168;
t140 = -t193 * t324 + t212;
t139 = t367 + t448;
t138 = -t306 * t454 - t307 * t331 - t326 * t332 + t363 + t368;
t137 = t238 * t328 - t240 * t308 + t242 * t309;
t136 = t237 * t328 - t239 * t308 + t241 * t309;
t135 = t238 * t326 - t240 * t306 + t242 * t307;
t134 = t237 * t326 - t239 * t306 + t241 * t307;
t132 = t133 * t349;
t131 = t133 * t324;
t130 = t202 * t308 - t203 * t306;
t127 = -t195 * t306 + t163;
t124 = t328 * t403 + t406 * t419;
t123 = t202 * t419 + t236 * t326 + t402;
t118 = (-t202 + t400) * t349 + t356 * t371;
t117 = t203 * t349 + t353 * t371 + t401;
t108 = t202 * t328 + t326 * t406 + t245;
t105 = (t202 * t353 + t203 * t356) * t348 + t372;
t100 = t197 * t308 + t199 * t271 + t201 * t272;
t99 = t196 * t308 + t198 * t271 + t200 * t272;
t98 = t197 * t306 + t199 * t269 + t201 * t270;
t97 = t196 * t306 + t198 * t269 + t200 * t270;
t92 = t308 * t404 + t165 + t168;
t91 = t212 + t213 + (-t173 - t193) * t324;
t82 = -t308 * t405 + t412;
t81 = -t324 * t411 + t407;
t80 = t328 * t388 + t390 * t419;
t79 = t193 * t419 + t230 * t326 + t374;
t78 = (-t193 + t391) * t349 + t356 * t366;
t77 = t195 * t349 + t353 * t366 + t392;
t76 = t306 * t408 + t160 + t163;
t75 = t161 + (-t142 * t356 + t143 * t353) * t348;
t74 = t142 * t326 + t143 * t328 - t164 * t419;
t73 = -t306 * t410 + t413;
t72 = t193 * t328 + t326 * t390 + t409;
t71 = (t193 * t353 + t195 * t356) * t348 + t365;
t70 = t153 * t349 + (-t136 * t356 + t137 * t353) * t348;
t69 = t152 * t349 + (-t134 * t356 + t135 * t353) * t348;
t68 = t136 * t326 + t137 * t328 - t153 * t419;
t67 = t134 * t326 + t135 * t328 - t152 * t419;
t66 = t308 * t389 + t165 + t412;
t65 = t213 + (-t173 - t411) * t324 + t407;
t64 = t328 * t373 + t375 * t419;
t63 = t326 * t405 + t411 * t419 + t374;
t62 = (t391 - t411) * t349 + t356 * t364;
t61 = t349 * t410 + t353 * t364 + t392;
t60 = t306 * t393 + t160 + t413;
t59 = t326 * t375 + t328 * t411 + t409;
t58 = (t353 * t411 + t356 * t410) * t348 + t365;
t57 = t132 + (-t106 * t356 + t107 * t353) * t348;
t56 = t106 * t326 + t107 * t328 - t133 * t419;
t55 = t106 * t306 + t107 * t308 + t131;
t38 = t120 * t349 + (t100 * t353 - t356 * t99) * t348;
t37 = t119 * t349 + (t353 * t98 - t356 * t97) * t348;
t36 = t100 * t328 - t120 * t419 + t326 * t99;
t35 = -t119 * t419 + t326 * t97 + t328 * t98;
t34 = t100 * t308 + t120 * t324 + t306 * t99;
t33 = t119 * t324 + t306 * t97 + t308 * t98;
t1 = [m(7) * (t138 ^ 2 + t139 ^ 2) + m(6) * (t144 ^ 2 + t145 ^ 2) + m(5) * (t149 ^ 2 + t150 ^ 2) + m(4) * (t204 ^ 2 + t205 ^ 2) + m(3) * (t264 ^ 2 + t265 ^ 2) + m(2) * (t333 ^ 2 + t334 ^ 2) + Icges(2,3) + t385 - t438; t232 + t126 + t161 + t125 + t132 + m(3) * (t246 * t265 + t247 * t264) + m(4) * (t155 * t205 + t156 * t204) + m(5) * (t117 * t150 + t118 * t149) + m(7) * (t138 * t62 + t139 * t61) + m(6) * (t144 * t78 + t145 * t77) + ((-t176 / 0.2e1 - t219 / 0.2e1 - t359) * t356 + (t177 / 0.2e1 + t220 / 0.2e1 + t358) * t353) * t348; (t57 + t75 + t232 + t439) * t349 + m(7) * (t58 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t71 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t105 ^ 2 + t117 ^ 2 + t118 ^ 2) + m(4) * (t148 ^ 2 + t155 ^ 2 + t156 ^ 2) + m(3) * (t221 ^ 2 + t246 ^ 2 + t247 ^ 2) + ((-t37 - t69 + (-t279 * t326 + t281 * t327 - t348 * t424) * t418 - t443) * t356 + (t38 + t70 + ((-t280 * t328 + t282 * t329 + (t278 * t353 - t424) * t348) * t353 + (t278 * t418 + t279 * t328 + t280 * t326 - t281 * t329 - t282 * t327) * t356) * t348 + t442) * t353 + ((-t176 - t219) * t356 + (t177 + t220) * t353) * t349) * t348; t438 * t419 + m(7) * (t138 * t63 + t139 * t64) + m(6) * (t144 * t79 + t145 * t80) + m(5) * (t123 * t149 + t124 * t150) + m(4) * (t190 * t204 + t191 * t205) + t358 * t328 + t359 * t326; (t50 / 0.2e1 + t56 / 0.2e1 + t74 / 0.2e1 + t49 / 0.2e1) * t349 + (t70 / 0.2e1 + t38 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1) * t328 + (t69 / 0.2e1 + t37 / 0.2e1 + t29 / 0.2e1 + t30 / 0.2e1) * t326 + m(4) * (t148 * t157 + t155 * t191 + t156 * t190) + m(5) * (t105 * t108 + t117 * t124 + t118 * t123) + m(7) * (t58 * t59 + t61 * t64 + t62 * t63) + m(6) * (t71 * t72 + t77 * t80 + t78 * t79) + ((-t25 / 0.2e1 - t26 / 0.2e1 - t35 / 0.2e1 - t67 / 0.2e1) * t356 + (-t53 / 0.2e1 - t54 / 0.2e1 - t57 / 0.2e1 - t75 / 0.2e1) * t355 + (t27 / 0.2e1 + t28 / 0.2e1 + t36 / 0.2e1 + t68 / 0.2e1) * t353) * t348; (-t56 - t74 - t440) * t419 + (t36 + t68 + t444) * t328 + (t35 + t67 + t445) * t326 + m(7) * (t59 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t72 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t108 ^ 2 + t123 ^ 2 + t124 ^ 2) + m(4) * (t157 ^ 2 + t190 ^ 2 + t191 ^ 2); t362 + m(7) * (t138 * t65 + t139 * t66) + m(6) * (t144 * t91 + t145 * t92) + m(5) * (t146 * t149 + t147 * t150) + t380 * t308 + t379 * t306 + t131; m(7) * (t58 * t60 + t61 * t66 + t62 * t65) + m(6) * (t71 * t76 + t77 * t92 + t78 * t91) + m(5) * (t105 * t130 + t117 * t147 + t118 * t146) + (t353 * t34 / 0.2e1 - t356 * t33 / 0.2e1) * t348 + t360 + t37 * t437 + t38 * t436 + t57 * t435 + t55 * t432; t361 + m(7) * (t59 * t60 + t63 * t65 + t64 * t66) + m(6) * (t72 * t76 + t79 * t91 + t80 * t92) + m(5) * (t108 * t130 + t123 * t146 + t124 * t147) + t55 * t382 + t35 * t437 + t36 * t436 + t56 * t435 + t33 * t434 + t34 * t433; t306 * t33 + t308 * t34 + t324 * t55 + m(7) * (t60 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t76 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t130 ^ 2 + t146 ^ 2 + t147 ^ 2) + t370; m(7) * (t138 * t81 + t139 * t82) + m(6) * (t140 * t144 + t141 * t145) + t362; m(6) * (t127 * t71 + t140 * t78 + t141 * t77) + m(7) * (t58 * t73 + t61 * t82 + t62 * t81) + t360; m(7) * (t59 * t73 + t63 * t81 + t64 * t82) + m(6) * (t127 * t72 + t140 * t79 + t141 * t80) + t361; m(7) * (t60 * t73 + t65 * t81 + t66 * t82) + m(6) * (t127 * t76 + t140 * t91 + t141 * t92) + t370; m(7) * (t73 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(6) * (t127 ^ 2 + t140 ^ 2 + t141 ^ 2) + t370; m(7) * (t138 * t308 + t139 * t306); m(7) * (t306 * t61 + t308 * t62 + t324 * t58); m(7) * (t306 * t64 + t308 * t63 + t324 * t59); m(7) * (t306 * t66 + t308 * t65 + t324 * t60); m(7) * (t306 * t82 + t308 * t81 + t324 * t73); m(7) * (t306 ^ 2 + t308 ^ 2 + t324 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
