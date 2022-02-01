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
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:31:52
% EndTime: 2022-01-23 09:32:17
% DurationCPUTime: 15.19s
% Computational Cost: add. (16353->571), mult. (22814->800), div. (0->0), fcn. (21426->8), ass. (0->268)
t458 = Icges(5,4) + Icges(6,4);
t456 = Icges(5,2) + Icges(6,2);
t455 = Icges(5,6) + Icges(6,6);
t292 = qJ(3) + qJ(4);
t284 = cos(t292);
t460 = t458 * t284;
t459 = Icges(5,1) + Icges(6,1);
t457 = Icges(5,5) + Icges(6,5);
t454 = -Icges(6,3) - Icges(5,3);
t283 = sin(t292);
t293 = sin(pkin(8));
t294 = cos(pkin(8));
t443 = -t455 * t294 + (-t283 * t456 + t460) * t293;
t453 = t458 * t283;
t291 = qJD(3) + qJD(4);
t383 = t291 * t293;
t447 = (t283 * t457 + t284 * t455) * t383;
t452 = (-t283 * t459 - t460) * t383;
t441 = t452 * t284 * t293 + t447 * t294;
t446 = (t456 * t284 + t453) * t383;
t385 = t284 * t291;
t448 = t443 * t385;
t449 = -t457 * t294 + (t459 * t284 - t453) * t293;
t434 = t441 + ((-t449 * t291 + t446) * t283 - t448) * t293;
t451 = t434 * t294;
t296 = sin(qJ(1));
t298 = cos(qJ(1));
t378 = t294 * t298;
t246 = t283 * t296 + t284 * t378;
t379 = t294 * t296;
t311 = t283 * t379 + t284 * t298;
t176 = qJD(1) * t311 - t246 * t291;
t244 = -t283 * t298 + t284 * t379;
t245 = -t283 * t378 + t284 * t296;
t177 = -qJD(1) * t244 + t245 * t291;
t355 = qJD(1) * t296;
t340 = t293 * t355;
t450 = t455 * t176 + t457 * t177 + t454 * t340;
t178 = qJD(1) * t245 - t244 * t291;
t179 = qJD(1) * t246 - t291 * t311;
t354 = qJD(1) * t298;
t339 = t293 * t354;
t423 = t455 * t178 + t457 * t179 - t454 * t339;
t422 = t456 * t176 + t458 * t177 - t455 * t340;
t421 = t456 * t178 + t458 * t179 + t455 * t339;
t420 = t458 * t176 + t459 * t177 - t457 * t340;
t419 = t458 * t178 + t459 * t179 + t457 * t339;
t382 = t293 * t296;
t417 = t457 * t244 - t455 * t311 - t454 * t382;
t380 = t293 * t298;
t416 = t455 * t245 + t457 * t246 - t454 * t380;
t415 = t458 * t244 - t456 * t311 + t455 * t382;
t414 = t456 * t245 + t458 * t246 + t455 * t380;
t413 = t459 * t244 - t458 * t311 + t457 * t382;
t412 = t458 * t245 + t459 * t246 + t457 * t380;
t444 = t454 * t294 + (-t455 * t283 + t457 * t284) * t293;
t299 = pkin(7) + pkin(6);
t289 = -qJ(5) - t299;
t432 = rSges(6,3) - t289;
t440 = t423 * t294 + ((t415 * t291 - t419) * t284 + (t413 * t291 + t421) * t283) * t293;
t439 = t450 * t294 + ((t414 * t291 - t420) * t284 + (t412 * t291 + t422) * t283) * t293;
t438 = t415 * t245 + t413 * t246 + t417 * t380;
t437 = -t414 * t245 - t412 * t246 - t416 * t380;
t436 = t417 * t294 + (t415 * t283 - t413 * t284) * t293;
t435 = -t416 * t294 + (-t414 * t283 + t412 * t284) * t293;
t433 = t443 * t245 + t449 * t246 + t444 * t380;
t286 = t298 * qJ(2);
t431 = rSges(6,1) * t244 - rSges(6,2) * t311 - t286;
t295 = sin(qJ(3));
t352 = qJD(3) * t295;
t349 = pkin(3) * t352;
t398 = pkin(4) * t283;
t256 = -t291 * t398 - t349;
t350 = qJD(5) * t293;
t430 = -rSges(6,1) * t179 - rSges(6,2) * t178 - (t256 * t294 + t350) * t296;
t429 = (t447 * t298 + t444 * t355) * t293 - t452 * t246 + t446 * t245 - t449 * t177 - t443 * t176;
t428 = -t446 * t311 + (t447 * t296 - t444 * t354) * t293 - t452 * t244 - t449 * t179 - t443 * t178;
t427 = t413 * t244 - t415 * t311 + t417 * t382;
t426 = t412 * t244 - t414 * t311 + t416 * t382;
t425 = -t244 * t449 + t443 * t311 - t444 * t382;
t297 = cos(qJ(3));
t351 = qJD(3) * t297;
t348 = pkin(3) * t351;
t257 = pkin(4) * t385 + t348;
t274 = qJD(2) + t348;
t397 = pkin(2) * t294 + pkin(1);
t399 = pkin(3) * t297;
t242 = t293 * t299 + t294 * t399 + t397;
t341 = -pkin(2) - t399;
t264 = pkin(4) * t284 - t341;
t336 = -t264 * t294 - pkin(1);
t325 = t242 + t336;
t309 = -t289 * t293 - t325;
t282 = qJD(2) * t298;
t400 = pkin(3) * t295;
t279 = qJ(2) + t400;
t329 = t296 * t349;
t358 = -t279 * t355 + t294 * t329;
t342 = t282 - t358;
t268 = t398 + t400;
t406 = (qJ(2) + t268) * t296;
t424 = rSges(6,3) * t339 + (-t257 + t274) * t298 + (t309 * t298 + t406) * qJD(1) - t342 - t430;
t418 = (-t268 + t279) * t298 + t309 * t296 + rSges(6,3) * t382 + t431;
t411 = (t299 - t432) * t294 + (rSges(6,1) * t284 - rSges(6,2) * t283 + t264 + t341) * t293;
t410 = -qJD(5) * t294 + (t256 + t349) * t293 + (-rSges(6,1) * t283 - rSges(6,2) * t284) * t383;
t409 = t450 * t298;
t408 = t451 + (t439 * t298 + t440 * t296 + (t435 * t296 + t436 * t298) * qJD(1)) * t293;
t285 = t296 * qJ(2);
t356 = t298 * pkin(1) + t285;
t119 = t246 * rSges(6,1) + t245 * rSges(6,2) + t264 * t378 + t296 * t268 + t432 * t380 + t356;
t405 = (t428 * t294 + (-t426 * t355 + t427 * t354 + (t414 * t178 + t412 * t179 + t420 * t244 - t422 * t311 + t416 * t339) * t298 + ((t423 * t296 + t417 * t354 + t409) * t293 + t419 * t244 - t421 * t311 + t413 * t179 + t415 * t178) * t296) * t293) * t382 + (t429 * t294 + (t437 * t355 + t438 * t354 + ((-t416 * t355 + t409) * t293 + t420 * t246 + t422 * t245 + t412 * t177 + t414 * t176) * t298 + ((t423 * t298 - t417 * t355) * t293 + t419 * t246 + t421 * t245 + t413 * t177 + t415 * t176) * t296) * t293) * t380 + (t425 * t294 + (t427 * t296 + t426 * t298) * t293) * t339 + (t433 * t294 + (-t438 * t296 + t437 * t298) * t293) * t340;
t404 = t177 * rSges(6,1) + t176 * rSges(6,2) + t256 * t378 + t296 * t257 + t268 * t354 + t289 * t340 + t298 * t350;
t403 = 2 * m(4);
t402 = 2 * m(5);
t401 = 2 * m(6);
t290 = t293 ^ 2;
t396 = t424 * t380;
t328 = t298 * t349;
t305 = t274 * t296 + t279 * t354 - t294 * t328;
t357 = qJ(2) * t354 + qJD(2) * t296;
t303 = t305 - t357;
t394 = rSges(6,3) * t340 - t325 * t355 + t303 - t404;
t393 = Icges(4,4) * t295;
t392 = Icges(4,4) * t297;
t387 = t279 * t298;
t381 = t293 * t297;
t353 = qJD(3) * t293;
t250 = (-Icges(4,2) * t297 - t393) * t353;
t377 = t295 * t250;
t376 = t295 * t296;
t375 = t295 * t298;
t374 = t296 * t297;
t373 = t297 * t298;
t364 = t177 * rSges(5,1) + t176 * rSges(5,2);
t115 = -rSges(5,3) * t340 + t364;
t272 = pkin(6) * t293 + t397;
t359 = t242 - t272;
t138 = -t359 * t355 + t303;
t371 = -t115 - t138;
t370 = t418 * t380;
t360 = t242 * t298 + t279 * t296;
t368 = t360 - t119;
t152 = t359 * t296 + t286 - t387;
t248 = pkin(3) * t381 + (pkin(6) - t299) * t294;
t367 = t294 * t152 + t248 * t382;
t330 = -t272 * t298 - t285;
t153 = t330 + t360;
t169 = t246 * rSges(5,1) + t245 * rSges(5,2) + rSges(5,3) * t380;
t366 = -t153 - t169;
t262 = t294 * t373 + t376;
t310 = t294 * t376 + t373;
t195 = qJD(1) * t310 - qJD(3) * t262;
t260 = t294 * t374 - t375;
t261 = -t294 * t375 + t374;
t196 = -qJD(1) * t260 + qJD(3) * t261;
t362 = t196 * rSges(4,1) + t195 * rSges(4,2);
t318 = -rSges(5,1) * t244 + rSges(5,2) * t311;
t167 = rSges(5,3) * t382 - t318;
t231 = -rSges(5,3) * t294 + (rSges(5,1) * t284 - rSges(5,2) * t283) * t293;
t134 = t294 * t167 + t231 * t382;
t361 = t248 * t340 + t290 * t328;
t347 = -t138 + t394;
t319 = -rSges(5,1) * t179 - rSges(5,2) * t178;
t117 = rSges(5,3) * t339 - t319;
t219 = (-rSges(5,1) * t283 - rSges(5,2) * t284) * t383;
t63 = t294 * t117 + t219 * t382 + t231 * t339;
t343 = -t153 + t368;
t188 = t262 * rSges(4,1) + t261 * rSges(4,2) + rSges(4,3) * t380;
t335 = -rSges(4,3) * t293 - t272;
t334 = -rSges(5,3) * t293 - t242;
t249 = (-Icges(4,5) * t295 - Icges(4,6) * t297) * t353;
t251 = (-Icges(4,1) * t295 - t392) * t353;
t331 = -t294 * t249 + t251 * t381;
t48 = t418 * t294 + t411 * t382;
t323 = -t219 * t380 + t231 * t340;
t322 = rSges(3,1) * t294 - rSges(3,2) * t293;
t197 = qJD(1) * t261 - qJD(3) * t260;
t198 = qJD(1) * t262 - qJD(3) * t310;
t321 = -rSges(4,1) * t198 - rSges(4,2) * t197;
t320 = -rSges(4,1) * t260 + rSges(4,2) * t310;
t315 = t335 * t296;
t314 = t334 * qJD(1);
t312 = -pkin(1) - t322;
t25 = t424 * t294 + t411 * t339 + t410 * t382;
t308 = t411 * t340 - t410 * t380;
t307 = -t432 * t293 + t336;
t139 = -t274 * t298 + (t359 * t298 - t285) * qJD(1) + t342;
t306 = t294 * t139 + t248 * t339 - t290 * t329;
t302 = rSges(3,3) * t298 + t296 * t312;
t301 = (-t428 - t440) * t382 / 0.2e1 - (t433 + t435) * t340 / 0.2e1 + (-t429 + (-t425 - t436) * qJD(1) - t439) * t380 / 0.2e1;
t300 = t294 * t408 + t405;
t252 = (-rSges(4,1) * t295 - rSges(4,2) * t297) * t353;
t241 = -rSges(4,3) * t294 + (rSges(4,1) * t297 - rSges(4,2) * t295) * t293;
t240 = -Icges(4,5) * t294 + (Icges(4,1) * t297 - t393) * t293;
t239 = -Icges(4,6) * t294 + (-Icges(4,2) * t295 + t392) * t293;
t238 = -Icges(4,3) * t294 + (Icges(4,5) * t297 - Icges(4,6) * t295) * t293;
t211 = rSges(3,3) * t296 + t298 * t322 + t356;
t210 = t286 + t302;
t191 = t282 + ((-rSges(3,3) - qJ(2)) * t296 + t312 * t298) * qJD(1);
t190 = qJD(1) * t302 + t357;
t187 = rSges(4,3) * t382 - t320;
t186 = Icges(4,1) * t262 + Icges(4,4) * t261 + Icges(4,5) * t380;
t185 = Icges(4,1) * t260 - Icges(4,4) * t310 + Icges(4,5) * t382;
t184 = Icges(4,4) * t262 + Icges(4,2) * t261 + Icges(4,6) * t380;
t183 = Icges(4,4) * t260 - Icges(4,2) * t310 + Icges(4,6) * t382;
t182 = Icges(4,5) * t262 + Icges(4,6) * t261 + Icges(4,3) * t380;
t181 = Icges(4,5) * t260 - Icges(4,6) * t310 + Icges(4,3) * t382;
t148 = t167 * t380;
t146 = t152 * t380;
t145 = t286 + t315 + t320;
t144 = -t330 + t188;
t143 = -t188 * t294 - t241 * t380;
t142 = t187 * t294 + t241 * t382;
t141 = t169 + t360;
t140 = t296 * t334 + t318 + t387;
t135 = -t169 * t294 - t231 * t380;
t133 = rSges(4,3) * t339 - t321;
t132 = -rSges(4,3) * t340 + t362;
t131 = Icges(4,1) * t198 + Icges(4,4) * t197 + Icges(4,5) * t339;
t130 = Icges(4,1) * t196 + Icges(4,4) * t195 - Icges(4,5) * t340;
t129 = Icges(4,4) * t198 + Icges(4,2) * t197 + Icges(4,6) * t339;
t128 = Icges(4,4) * t196 + Icges(4,2) * t195 - Icges(4,6) * t340;
t127 = Icges(4,5) * t198 + Icges(4,6) * t197 + Icges(4,3) * t339;
t126 = Icges(4,5) * t196 + Icges(4,6) * t195 - Icges(4,3) * t340;
t123 = t139 * t380;
t121 = t238 * t380 + t239 * t261 + t240 * t262;
t120 = t238 * t382 - t239 * t310 + t240 * t260;
t118 = t268 * t298 + t296 * t307 - t431;
t101 = qJD(1) * t315 + t357 + t362;
t100 = t282 + (t298 * t335 - t285) * qJD(1) + t321;
t97 = -t169 * t382 + t148;
t96 = t117 * t380;
t90 = (-t377 + (-t239 * t297 - t240 * t295) * qJD(3)) * t293 + t331;
t85 = t133 * t294 + (t241 * t354 + t252 * t296) * t293;
t84 = -t132 * t294 + (t241 * t355 - t252 * t298) * t293;
t83 = -t182 * t294 + (-t184 * t295 + t186 * t297) * t293;
t82 = -t181 * t294 + (-t183 * t295 + t185 * t297) * t293;
t81 = t366 * t294 + (-t231 - t248) * t380;
t80 = t134 + t367;
t79 = t296 * t314 + t305 + t364;
t78 = (t274 + t314) * t298 + t319 + t358;
t77 = t182 * t380 + t184 * t261 + t186 * t262;
t76 = t181 * t380 + t183 * t261 + t185 * t262;
t75 = t182 * t382 - t184 * t310 + t186 * t260;
t74 = t181 * t382 - t183 * t310 + t185 * t260;
t62 = -t115 * t294 + t323;
t49 = t368 * t294 - t380 * t411;
t47 = t257 * t298 + t282 + (t307 * t298 - t406) * qJD(1) + t430;
t46 = (-rSges(6,3) * t293 + t336) * t355 + t357 + t404;
t45 = t366 * t382 + t146 + t148;
t44 = t197 * t239 + t198 * t240 - t250 * t310 + t251 * t260 + (t238 * t354 + t249 * t296) * t293;
t43 = t195 * t239 + t196 * t240 + t250 * t261 + t251 * t262 + (-t238 * t355 + t249 * t298) * t293;
t42 = t368 * t382 + t370;
t41 = t343 * t294 + (-t248 - t411) * t380;
t40 = t48 + t367;
t39 = t306 + t63;
t38 = t371 * t294 + t323 + t361;
t29 = t96 + (-t115 * t296 + (-t167 * t296 - t169 * t298) * qJD(1)) * t293;
t28 = -t126 * t294 + (-t128 * t295 + t130 * t297 + (-t184 * t297 - t186 * t295) * qJD(3)) * t293;
t27 = -t127 * t294 + (-t129 * t295 + t131 * t297 + (-t183 * t297 - t185 * t295) * qJD(3)) * t293;
t26 = t343 * t382 + t146 + t370;
t24 = t394 * t294 + t308;
t11 = t306 + t25;
t10 = t294 * t347 + t308 + t361;
t9 = t123 + t96 + (t371 * t296 + (t366 * t298 + (-t152 - t167) * t296) * qJD(1)) * t293;
t8 = (t394 * t296 + (-t296 * t418 + t368 * t298) * qJD(1)) * t293 + t396;
t7 = t123 + (t347 * t296 + (t343 * t298 + (-t152 - t418) * t296) * qJD(1)) * t293 + t396;
t1 = [(t118 * t47 + t119 * t46) * t401 + (t140 * t78 + t141 * t79) * t402 + (t100 * t145 + t101 * t144) * t403 + 0.2e1 * m(3) * (t190 * t211 + t191 * t210) + t331 + (-t239 * t351 - t240 * t352 - t377 - t448) * t293 + (t446 * t293 - t383 * t449) * t283 + t441; m(6) * (t296 * t47 - t298 * t46 + (t118 * t298 + t119 * t296) * qJD(1)) + m(5) * (t296 * t78 - t298 * t79 + (t140 * t298 + t141 * t296) * qJD(1)) + m(4) * (t100 * t296 - t101 * t298 + (t144 * t296 + t145 * t298) * qJD(1)) + m(3) * (-t190 * t298 + t191 * t296 + (t210 * t298 + t211 * t296) * qJD(1)); 0; ((t43 / 0.2e1 + t28 / 0.2e1) * t298 + (t27 / 0.2e1 + t44 / 0.2e1) * t296 + ((t82 / 0.2e1 + t120 / 0.2e1) * t298 + (-t83 / 0.2e1 - t121 / 0.2e1) * t296) * qJD(1)) * t293 + m(6) * (t10 * t119 + t11 * t118 + t40 * t47 + t41 * t46) + m(5) * (t140 * t39 + t141 * t38 + t78 * t80 + t79 * t81) + m(4) * (t100 * t142 + t101 * t143 + t144 * t84 + t145 * t85) + (-t90 - t434) * t294 + t301; m(4) * (t296 * t85 - t298 * t84 + (t142 * t298 + t143 * t296) * qJD(1)) + m(5) * (t296 * t39 - t298 * t38 + (t296 * t81 + t298 * t80) * qJD(1)) + m(6) * (-t10 * t298 + t11 * t296 + (t296 * t41 + t298 * t40) * qJD(1)); (t142 * t85 + t143 * t84 + (t187 * t298 - t188 * t296) * (-t132 * t296 + t133 * t298 + (-t187 * t296 - t188 * t298) * qJD(1)) * t290) * t403 + (t10 * t41 + t11 * t40 + t26 * t7) * t401 + (t38 * t81 + t39 * t80 + t45 * t9) * t402 + ((t296 * t74 + t298 * t75) * t339 - (t296 * t76 + t298 * t77) * t340 + ((-t128 * t310 + t260 * t130 + t197 * t184 + t198 * t186) * t298 - t75 * t355 + (-t129 * t310 + t260 * t131 + t197 * t183 + t198 * t185) * t296 + t74 * t354) * t382 + ((t261 * t128 + t262 * t130 + t195 * t184 + t196 * t186) * t298 - t77 * t355 + (t261 * t129 + t262 * t131 + t195 * t183 + t196 * t185) * t296 + t76 * t354) * t380 + (((t126 * t296 + t182 * t354) * t298 + (t127 * t296 + t181 * t354) * t296) * t382 + ((t126 * t298 - t182 * t355) * t298 + (t127 * t298 - t181 * t355) * t296) * t380) * t293) * t293 + (-(t27 * t296 + t28 * t298 + (-t296 * t83 + t298 * t82) * qJD(1)) * t293 - t120 * t339 + t121 * t340 - t44 * t382 - t43 * t380 + t90 * t294 + t408) * t294 + t405; t301 + m(6) * (t118 * t25 + t119 * t24 + t49 * t46 + t48 * t47) + m(5) * (t134 * t78 + t135 * t79 + t140 * t63 + t141 * t62) - t451; m(5) * (t296 * t63 - t298 * t62 + (t134 * t298 + t135 * t296) * qJD(1)) + m(6) * (-t24 * t298 + t25 * t296 + (t296 * t49 + t298 * t48) * qJD(1)); m(6) * (t10 * t49 + t11 * t48 + t24 * t41 + t25 * t40 + t26 * t8 + t42 * t7) + m(5) * (t134 * t39 + t135 * t38 + t29 * t45 + t62 * t81 + t63 * t80 + t9 * t97) + t300; (t134 * t63 + t135 * t62 + t29 * t97) * t402 + (t24 * t49 + t25 * t48 + t42 * t8) * t401 + t300; m(6) * (t296 * t46 + t298 * t47 + (-t118 * t296 + t119 * t298) * qJD(1)) * t293; 0; m(6) * (-t294 * t7 + (t10 * t296 + t11 * t298 + (-t296 * t40 + t298 * t41) * qJD(1)) * t293); m(6) * (-t294 * t8 + (t24 * t296 + t25 * t298 + (-t296 * t48 + t298 * t49) * qJD(1)) * t293); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
