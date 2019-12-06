% Calculate time derivative of joint inertia matrix for
% S5RPRRP4
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:31
% EndTime: 2019-12-05 18:06:15
% DurationCPUTime: 19.82s
% Computational Cost: add. (15761->588), mult. (22702->824), div. (0->0), fcn. (21282->8), ass. (0->282)
t470 = Icges(5,4) + Icges(6,4);
t468 = Icges(5,2) + Icges(6,2);
t467 = Icges(5,6) + Icges(6,6);
t282 = qJ(3) + qJ(4);
t276 = cos(t282);
t472 = t470 * t276;
t471 = Icges(5,1) + Icges(6,1);
t469 = Icges(5,5) + Icges(6,5);
t466 = -Icges(6,3) - Icges(5,3);
t275 = sin(t282);
t283 = sin(pkin(8));
t284 = cos(pkin(8));
t458 = -t467 * t284 + (-t468 * t275 + t472) * t283;
t465 = t470 * t275;
t281 = qJD(3) + qJD(4);
t387 = t281 * t283;
t461 = (t469 * t275 + t467 * t276) * t387;
t464 = (-t471 * t275 - t472) * t387;
t288 = cos(qJ(1));
t286 = sin(qJ(1));
t383 = t284 * t286;
t240 = t275 * t383 + t276 * t288;
t382 = t284 * t288;
t243 = t275 * t286 + t276 * t382;
t167 = qJD(1) * t240 - t243 * t281;
t300 = t275 * t382 - t276 * t286;
t301 = -t275 * t288 + t276 * t383;
t168 = -qJD(1) * t301 - t281 * t300;
t352 = qJD(1) * t286;
t335 = t283 * t352;
t463 = t467 * t167 + t469 * t168 + t466 * t335;
t454 = t464 * t276 * t283 + t461 * t284;
t455 = t469 * t284 + (-t471 * t276 + t465) * t283;
t456 = (-t468 * t276 - t465) * t387;
t389 = t276 * t281;
t457 = t458 * t389;
t444 = t454 + ((t455 * t281 - t456) * t275 - t457) * t283;
t462 = t444 * t284;
t169 = qJD(1) * t300 + t281 * t301;
t170 = -qJD(1) * t243 + t240 * t281;
t351 = qJD(1) * t288;
t334 = t283 * t351;
t441 = -t467 * t169 - t469 * t170 - t466 * t334;
t440 = -t468 * t167 - t470 * t168 + t467 * t335;
t439 = t468 * t169 + t470 * t170 - t467 * t334;
t438 = t470 * t167 + t471 * t168 - t469 * t335;
t437 = -t470 * t169 - t471 * t170 + t469 * t334;
t386 = t283 * t286;
t436 = t467 * t240 - t469 * t301 + t466 * t386;
t385 = t283 * t288;
t435 = t469 * t243 - t467 * t300 - t466 * t385;
t434 = t468 * t240 - t470 * t301 - t467 * t386;
t433 = t470 * t243 - t468 * t300 + t467 * t385;
t432 = t470 * t240 - t471 * t301 - t469 * t386;
t431 = t471 * t243 - t470 * t300 + t469 * t385;
t459 = t466 * t284 + (-t467 * t275 + t469 * t276) * t283;
t453 = t441 * t284 + ((-t434 * t281 - t437) * t276 + (-t432 * t281 - t439) * t275) * t283;
t452 = t463 * t284 + ((t433 * t281 - t438) * t276 + (t431 * t281 - t440) * t275) * t283;
t451 = -t464 * t301 + (t461 * t286 - t459 * t351) * t283 + t456 * t240 - t455 * t170 + t458 * t169;
t450 = t434 * t240 - t432 * t301 - t436 * t386;
t449 = t433 * t240 - t431 * t301 - t435 * t386;
t448 = t432 * t243 - t434 * t300 + t436 * t385;
t447 = -t431 * t243 + t433 * t300 - t435 * t385;
t446 = -t436 * t284 + (-t434 * t275 + t432 * t276) * t283;
t445 = -t435 * t284 + (-t433 * t275 + t431 * t276) * t283;
t443 = t458 * t240 + t455 * t301 - t459 * t386;
t442 = -t455 * t243 - t458 * t300 + t459 * t385;
t285 = sin(qJ(3));
t349 = qJD(3) * t285;
t345 = pkin(3) * t349;
t406 = pkin(4) * t275;
t251 = -t281 * t406 - t345;
t306 = qJD(5) * t283 + t251 * t284;
t430 = -rSges(6,1) * t168 - rSges(6,2) * t167 - t306 * t288;
t407 = pkin(3) * t285;
t263 = t406 + t407;
t289 = -pkin(7) - pkin(6);
t279 = -qJ(5) + t289;
t429 = -rSges(6,1) * t301 + t240 * rSges(6,2) + t288 * t263 + t279 * t386;
t428 = t463 * t288;
t427 = t456 * t300 + (t461 * t288 + t459 * t352) * t283 - t464 * t243 + t455 * t168 - t458 * t167;
t287 = cos(qJ(3));
t348 = qJD(3) * t287;
t344 = pkin(3) * t348;
t252 = pkin(4) * t389 + t344;
t326 = -t263 + t407;
t354 = t279 - t289;
t272 = t287 * pkin(3) + pkin(2);
t343 = pkin(3) * qJD(3) * t288;
t317 = t285 * t343;
t332 = t284 * t352;
t358 = t272 * t332 + t284 * t317;
t259 = pkin(4) * t276 + t272;
t391 = t259 * t284;
t426 = -rSges(6,3) * t335 + (t252 - t344) * t286 + (-t326 * t288 + (t283 * t354 - t391) * t286) * qJD(1) + t358 - t430;
t357 = t259 - t272;
t325 = t357 * t284;
t379 = t285 * t288;
t384 = t283 * t289;
t356 = pkin(3) * t379 + t286 * t384;
t425 = -rSges(6,3) * t386 - t286 * t325 - t356 + t429;
t267 = t288 * t384;
t307 = -rSges(6,1) * t243 + rSges(6,2) * t300;
t413 = t326 * t286;
t424 = t267 - t413 + (-t279 * t283 + t325) * t288 + rSges(6,3) * t385 - t307;
t423 = (t354 - rSges(6,3)) * t284 + (rSges(6,1) * t276 - rSges(6,2) * t275 + t357) * t283;
t422 = -qJD(5) * t284 + (t251 + t345) * t283 + (-rSges(6,1) * t275 - rSges(6,2) * t276) * t387;
t404 = pkin(2) - t272;
t302 = pkin(6) * t283 + t284 * t404;
t380 = t285 * t286;
t421 = -pkin(3) * t380 + t288 * t302;
t420 = t451 * t284 + (t449 * t352 + t450 * t351 + (-t169 * t433 - t170 * t431 + t240 * t440 + t301 * t438 + t334 * t435) * t288 + (t437 * t301 + (t286 * t441 - t351 * t436 + t428) * t283 + t439 * t240 + t432 * t170 + t434 * t169) * t286) * t283;
t419 = t462 + (t452 * t288 + t453 * t286 + (t286 * t445 + t288 * t446) * qJD(1)) * t283;
t418 = t443 * t284 + (t286 * t450 - t288 * t449) * t283;
t417 = t442 * t284 + (t286 * t448 + t288 * t447) * t283;
t416 = t283 * (rSges(6,3) - t279);
t412 = t170 * rSges(6,1) + t169 * rSges(6,2) + t288 * t252 + t279 * t334 - t286 * t306;
t411 = 2 * m(4);
t410 = 2 * m(5);
t409 = 2 * m(6);
t280 = t283 ^ 2;
t408 = (t427 * t284 + (t447 * t352 - t448 * t351 + ((-t352 * t435 + t428) * t283 + t438 * t243 + t440 * t300 + t431 * t168 + t433 * t167) * t288 + ((t288 * t441 + t352 * t436) * t283 + t437 * t243 + t439 * t300 - t432 * t168 - t434 * t167) * t286) * t283) * t385;
t402 = -rSges(3,3) - qJ(2);
t318 = t286 * t345;
t353 = qJD(1) * t283;
t333 = t289 * t353;
t336 = t284 * t318 + t287 * t343 + t288 * t333;
t399 = rSges(6,3) * t334 - (-t288 * t325 + t413) * qJD(1) + t336 - t412;
t310 = -t168 * rSges(5,1) - t167 * rSges(5,2);
t112 = -rSges(5,3) * t335 - t310;
t221 = (-rSges(5,1) * t275 - rSges(5,2) * t276) * t387;
t398 = t284 * t112 + t221 * t385;
t397 = Icges(4,4) * t285;
t396 = Icges(4,4) * t287;
t350 = qJD(3) * t283;
t246 = (-Icges(4,2) * t287 - t397) * t350;
t381 = t285 * t246;
t378 = t286 * t287;
t377 = t287 * t288;
t277 = t288 * qJ(2);
t376 = -qJ(2) - t263;
t368 = t170 * rSges(5,1) + t169 * rSges(5,2);
t114 = -rSges(5,3) * t334 + t368;
t139 = qJD(1) * t421 + t336;
t375 = -t114 - t139;
t374 = t425 * t335;
t309 = -rSges(5,1) * t243 + rSges(5,2) * t300;
t161 = rSges(5,3) * t385 - t309;
t231 = -rSges(5,3) * t284 + (rSges(5,1) * t276 - rSges(5,2) * t275) * t283;
t130 = t284 * t161 + t231 * t385;
t361 = -rSges(5,1) * t301 + t240 * rSges(5,2);
t159 = -rSges(5,3) * t386 + t361;
t181 = t286 * t302 + t356;
t371 = -t159 - t181;
t182 = -t267 - t421;
t370 = -t161 - t182;
t222 = (pkin(6) + t289) * t284 - t404 * t283;
t367 = t284 * t182 + t222 * t385;
t366 = t423 * t386;
t298 = t284 * t379 - t378;
t299 = t284 * t378 - t379;
t194 = qJD(1) * t298 + qJD(3) * t299;
t253 = t284 * t380 + t377;
t256 = t284 * t377 + t380;
t195 = -qJD(1) * t256 + qJD(3) * t253;
t364 = t195 * rSges(4,1) + t194 * rSges(4,2);
t363 = t221 * t386 + t231 * t334;
t360 = -rSges(4,1) * t299 + t253 * rSges(4,2);
t355 = pkin(2) * t332 + pkin(6) * t335;
t346 = rSges(4,3) * t386;
t342 = -t139 + t399;
t338 = -t181 - t425;
t337 = -t182 - t424;
t328 = -pkin(1) - t391;
t327 = -qJ(2) - t407;
t324 = t376 * t286;
t245 = (-Icges(4,5) * t285 - Icges(4,6) * t287) * t350;
t247 = (-Icges(4,1) * t285 - t396) * t350;
t321 = t283 * t287 * t247 - t284 * t245;
t320 = t284 * t426 + t385 * t422;
t273 = pkin(1) * t352;
t319 = -qJD(2) * t286 + t273;
t49 = t284 * t424 + t385 * t423;
t316 = t334 * t423 + t386 * t422;
t313 = rSges(3,1) * t284 - rSges(3,2) * t283;
t192 = qJD(1) * t253 - qJD(3) * t256;
t193 = -qJD(1) * t299 - qJD(3) * t298;
t312 = -rSges(4,1) * t193 - rSges(4,2) * t192;
t311 = -rSges(4,1) * t256 + rSges(4,2) * t298;
t305 = -pkin(1) - t313;
t304 = -rSges(5,3) * t283 - t272 * t284 - pkin(1);
t303 = -rSges(6,3) * t283 + t328;
t138 = t286 * t333 + (t285 * t351 + t286 * t348) * pkin(3) + t355 - t358;
t297 = t284 * t138 - t280 * t317;
t296 = t222 * t334 - t280 * t318;
t295 = -pkin(2) * t284 - pkin(1) + (-rSges(4,3) - pkin(6)) * t283;
t293 = -qJ(2) * t286 + t288 * t295;
t213 = t286 * t402 + t288 * t305;
t292 = t286 * t327 + t288 * t304;
t291 = -(t451 + t453) * t386 / 0.2e1 + (-t427 - t452) * t385 / 0.2e1 - ((t442 + t445) * t286 + (t443 + t446) * t288) * t353 / 0.2e1;
t290 = t419 * t284 + (t420 * t286 + (t417 * t286 + t418 * t288) * qJD(1)) * t283 + t408;
t274 = qJD(2) * t288;
t248 = (-rSges(4,1) * t285 - rSges(4,2) * t287) * t350;
t239 = -rSges(4,3) * t284 + (rSges(4,1) * t287 - rSges(4,2) * t285) * t283;
t238 = -Icges(4,5) * t284 + (Icges(4,1) * t287 - t397) * t283;
t237 = -Icges(4,6) * t284 + (-Icges(4,2) * t285 + t396) * t283;
t236 = -Icges(4,3) * t284 + (Icges(4,5) * t287 - Icges(4,6) * t285) * t283;
t212 = rSges(3,3) * t288 + t286 * t305 + t277;
t209 = t231 * t386;
t206 = t222 * t386;
t186 = t213 * qJD(1) + t274;
t185 = (t286 * t313 + t288 * t402) * qJD(1) + t319;
t184 = rSges(4,3) * t385 - t311;
t183 = -t346 + t360;
t180 = Icges(4,1) * t256 - Icges(4,4) * t298 + Icges(4,5) * t385;
t179 = -Icges(4,1) * t299 + Icges(4,4) * t253 - Icges(4,5) * t386;
t178 = Icges(4,4) * t256 - Icges(4,2) * t298 + Icges(4,6) * t385;
t177 = -Icges(4,4) * t299 + Icges(4,2) * t253 - Icges(4,6) * t386;
t176 = Icges(4,5) * t256 - Icges(4,6) * t298 + Icges(4,3) * t385;
t175 = -Icges(4,5) * t299 + Icges(4,6) * t253 - Icges(4,3) * t386;
t162 = t181 * t335;
t143 = t159 * t335;
t141 = t293 + t311;
t140 = t286 * t295 + t277 + t360;
t134 = t184 * t284 + t239 * t385;
t133 = -t183 * t284 + t239 * t386;
t129 = -t159 * t284 + t209;
t128 = -rSges(4,3) * t334 + t364;
t127 = -rSges(4,3) * t335 - t312;
t126 = Icges(4,1) * t195 + Icges(4,4) * t194 - Icges(4,5) * t334;
t125 = Icges(4,1) * t193 + Icges(4,4) * t192 - Icges(4,5) * t335;
t124 = Icges(4,4) * t195 + Icges(4,2) * t194 - Icges(4,6) * t334;
t123 = Icges(4,4) * t193 + Icges(4,2) * t192 - Icges(4,6) * t335;
t122 = Icges(4,5) * t195 + Icges(4,6) * t194 - Icges(4,3) * t334;
t121 = Icges(4,5) * t193 + Icges(4,6) * t192 - Icges(4,3) * t335;
t120 = t267 + t292 + t309;
t119 = t286 * t304 + t277 + t356 + t361;
t118 = t236 * t385 - t237 * t298 + t238 * t256;
t117 = -t236 * t386 + t237 * t253 - t238 * t299;
t116 = t324 + (t328 - t416) * t288 + t307;
t115 = t286 * t303 + t277 + t429;
t96 = (-t159 * t288 - t161 * t286) * t283;
t95 = qJD(1) * t293 + t274 + t364;
t94 = (t346 - t277) * qJD(1) + t312 + t319 + t355;
t89 = (-t381 + (-t237 * t287 - t238 * t285) * qJD(3)) * t283 + t321;
t84 = -t128 * t284 + (t239 * t351 + t248 * t286) * t283;
t83 = t127 * t284 + (-t239 * t352 + t248 * t288) * t283;
t82 = -t176 * t284 + (-t178 * t285 + t180 * t287) * t283;
t81 = -t175 * t284 + (-t177 * t285 + t179 * t287) * t283;
t77 = t367 + t130;
t76 = t284 * t371 + t206 + t209;
t75 = t176 * t385 - t178 * t298 + t180 * t256;
t74 = t175 * t385 - t177 * t298 + t179 * t256;
t73 = -t176 * t386 + t178 * t253 - t180 * t299;
t72 = -t175 * t386 + t177 * t253 - t179 * t299;
t65 = qJD(1) * t292 + t274 + t336 + t368;
t64 = t273 + (-qJD(2) - t344) * t286 + (t327 * t288 + (rSges(5,3) - t289) * t386) * qJD(1) + t310 + t358;
t63 = -t114 * t284 + t363;
t62 = -t231 * t335 + t398;
t48 = -t284 * t425 + t366;
t47 = (t286 * t370 + t288 * t371) * t283;
t46 = t274 + (t288 * t303 + t324) * qJD(1) + t412;
t45 = t273 + (-qJD(2) - t252) * t286 + (t376 * t288 + (t391 + t416) * t286) * qJD(1) + t430;
t44 = t194 * t237 + t195 * t238 + t246 * t253 - t247 * t299 + (-t236 * t351 - t245 * t286) * t283;
t43 = t192 * t237 + t193 * t238 - t246 * t298 + t247 * t256 + (-t236 * t352 + t245 * t288) * t283;
t42 = (-t286 * t424 - t288 * t425) * t283;
t41 = t49 + t367;
t40 = t284 * t338 + t206 + t366;
t39 = t284 * t375 + t296 + t363;
t38 = (-t222 - t231) * t335 + t297 + t398;
t29 = t143 + (-t112 * t286 + (-qJD(1) * t161 - t114) * t288) * t283;
t28 = (t286 * t337 + t288 * t338) * t283;
t27 = -t121 * t284 + (-t123 * t285 + t125 * t287 + (-t178 * t287 - t180 * t285) * qJD(3)) * t283;
t26 = -t122 * t284 + (-t124 * t285 + t126 * t287 + (-t177 * t287 - t179 * t285) * qJD(3)) * t283;
t25 = t284 * t399 + t316;
t24 = -t335 * t423 + t320;
t15 = t284 * t342 + t296 + t316;
t14 = (-t222 - t423) * t335 + t297 + t320;
t9 = t143 + t162 + ((-t112 - t138) * t286 + (qJD(1) * t370 + t375) * t288) * t283;
t8 = (-t426 * t286 + (-qJD(1) * t424 + t399) * t288) * t283 + t374;
t7 = t162 + ((-t138 - t426) * t286 + (qJD(1) * t337 + t342) * t288) * t283 + t374;
t1 = [(t115 * t46 + t116 * t45) * t409 + 0.2e1 * m(3) * (t185 * t213 + t186 * t212) + (t140 * t95 + t141 * t94) * t411 + (t119 * t65 + t120 * t64) * t410 + t321 + (-t237 * t348 - t238 * t349 - t381 - t457) * t283 + (-t283 * t456 + t455 * t387) * t275 + t454; m(6) * (t286 * t46 + t288 * t45 + (t115 * t288 - t116 * t286) * qJD(1)) + m(3) * (t185 * t288 + t186 * t286 + (t212 * t288 - t213 * t286) * qJD(1)) + m(4) * (t286 * t95 + t288 * t94 + (t140 * t288 - t141 * t286) * qJD(1)) + m(5) * (t286 * t65 + t288 * t64 + (t119 * t288 - t120 * t286) * qJD(1)); 0; t291 + ((t27 / 0.2e1 + t43 / 0.2e1) * t288 + (-t26 / 0.2e1 - t44 / 0.2e1) * t286 + ((-t81 / 0.2e1 - t117 / 0.2e1) * t288 + (-t82 / 0.2e1 - t118 / 0.2e1) * t286) * qJD(1)) * t283 + m(6) * (t115 * t15 + t116 * t14 + t40 * t46 + t41 * t45) + m(5) * (t119 * t39 + t120 * t38 + t64 * t77 + t65 * t76) + m(4) * (t133 * t95 + t134 * t94 + t140 * t84 + t141 * t83) + (-t89 - t444) * t284; m(4) * (t286 * t84 + t288 * t83 + (t133 * t288 - t134 * t286) * qJD(1)) + m(5) * (t286 * t39 + t288 * t38 + (-t286 * t77 + t288 * t76) * qJD(1)) + m(6) * (t14 * t288 + t15 * t286 + (-t286 * t41 + t288 * t40) * qJD(1)); (t14 * t41 + t15 * t40 + t28 * t7) * t409 + (t38 * t77 + t39 * t76 + t47 * t9) * t410 + (t133 * t84 + t134 * t83 + (-t183 * t288 - t184 * t286) * (-t127 * t286 - t128 * t288 + (t183 * t286 - t184 * t288) * qJD(1)) * t280) * t411 + t408 + t420 * t386 + t417 * t335 + t418 * t334 + (-(-t286 * t72 + t288 * t73) * t334 - (-(t253 * t124 - t126 * t299 + t194 * t177 + t195 * t179) * t286 - t72 * t351 + (t253 * t123 - t125 * t299 + t194 * t178 + t195 * t180) * t288 - t73 * t352) * t386 - (-t286 * t74 + t288 * t75) * t335 + (-(-t124 * t298 + t256 * t126 + t192 * t177 + t193 * t179) * t286 - t74 * t351 + (-t123 * t298 + t256 * t125 + t192 * t178 + t193 * t180) * t288 - t75 * t352) * t385 + (-(-(-t122 * t286 - t175 * t351) * t286 + (-t121 * t286 - t176 * t351) * t288) * t386 + (-(t122 * t288 - t175 * t352) * t286 + (t121 * t288 - t176 * t352) * t288) * t385) * t283) * t283 + (t117 * t334 + t44 * t386 + t118 * t335 - t43 * t385 - (-t26 * t286 + t27 * t288 + (-t286 * t82 - t288 * t81) * qJD(1)) * t283 + t89 * t284 + t419) * t284; m(6) * (t115 * t25 + t116 * t24 + t49 * t45 + t48 * t46) + m(5) * (t119 * t63 + t120 * t62 + t129 * t65 + t130 * t64) - t462 + t291; m(5) * (t286 * t63 + t288 * t62 + (t129 * t288 - t130 * t286) * qJD(1)) + m(6) * (t24 * t288 + t25 * t286 + (-t286 * t49 + t288 * t48) * qJD(1)); m(6) * (t14 * t49 + t15 * t48 + t24 * t41 + t25 * t40 + t28 * t8 + t42 * t7) + m(5) * (t129 * t39 + t130 * t38 + t29 * t47 + t62 * t77 + t63 * t76 + t9 * t96) + t290; (t129 * t63 + t130 * t62 + t29 * t96) * t410 + (t24 * t49 + t25 * t48 + t42 * t8) * t409 + t290; m(6) * (-t286 * t45 + t288 * t46 + (-t115 * t286 - t116 * t288) * qJD(1)) * t283; 0; m(6) * (-t284 * t7 + (-t14 * t286 + t15 * t288 + (-t286 * t40 - t288 * t41) * qJD(1)) * t283); m(6) * (-t284 * t8 + (-t24 * t286 + t25 * t288 + (-t286 * t48 - t288 * t49) * qJD(1)) * t283); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
