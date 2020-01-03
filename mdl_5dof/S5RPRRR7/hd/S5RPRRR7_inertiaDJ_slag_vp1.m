% Calculate time derivative of joint inertia matrix for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:13
% EndTime: 2019-12-31 19:03:34
% DurationCPUTime: 10.41s
% Computational Cost: add. (35208->738), mult. (36572->1027), div. (0->0), fcn. (35434->10), ass. (0->390)
t291 = sin(qJ(3));
t406 = qJD(1) * t291;
t490 = -t406 / 0.2e1;
t288 = qJ(1) + pkin(9);
t283 = cos(t288);
t294 = cos(qJ(3));
t402 = qJD(3) * t294;
t372 = t402 / 0.2e1;
t282 = sin(t288);
t383 = t282 * t406;
t489 = t283 * t372 - t383 / 0.2e1;
t407 = qJD(1) * t283;
t373 = t407 / 0.2e1;
t488 = t282 * t372 + t291 * t373;
t379 = t282 * t402;
t382 = t283 * t406;
t310 = t379 + t382;
t479 = 2 * m(4);
t462 = rSges(4,1) * t294;
t355 = -rSges(4,2) * t291 + t462;
t461 = rSges(4,3) * t283;
t220 = t282 * t355 - t461;
t427 = t283 * t294;
t429 = t283 * t291;
t482 = -rSges(4,2) * t429 + t282 * rSges(4,3);
t221 = rSges(4,1) * t427 + t482;
t267 = rSges(4,1) * t291 + rSges(4,2) * t294;
t318 = qJD(3) * t267;
t302 = rSges(4,2) * t383 + rSges(4,3) * t407 - t283 * t318;
t93 = (qJD(1) * t220 + t302) * t283 + (-t282 * t318 + (-t221 + t482) * qJD(1)) * t282;
t487 = t479 * t93;
t293 = cos(qJ(4));
t281 = pkin(4) * t293 + pkin(3);
t464 = pkin(3) - t281;
t486 = t291 * t464;
t296 = -pkin(8) - pkin(7);
t463 = pkin(7) + t296;
t485 = t294 * t463;
t484 = t294 * t464;
t457 = Icges(4,4) * t291;
t341 = Icges(4,1) * t294 - t457;
t219 = Icges(4,5) * t282 + t283 * t341;
t440 = t219 * t294;
t456 = Icges(4,4) * t294;
t337 = -Icges(4,2) * t291 + t456;
t217 = Icges(4,6) * t282 + t283 * t337;
t445 = t217 * t291;
t321 = -t440 + t445;
t483 = t321 * t283;
t286 = cos(qJ(1)) * pkin(1);
t481 = t282 * pkin(6) + t286;
t333 = Icges(4,5) * t294 - Icges(4,6) * t291;
t214 = -Icges(4,3) * t283 + t282 * t333;
t216 = -Icges(4,6) * t283 + t282 * t337;
t218 = -Icges(4,5) * t283 + t282 * t341;
t289 = qJ(4) + qJ(5);
t284 = sin(t289);
t426 = t284 * t294;
t285 = cos(t289);
t436 = t282 * t285;
t231 = -t283 * t426 + t436;
t425 = t285 * t294;
t437 = t282 * t284;
t232 = t283 * t425 + t437;
t166 = t232 * rSges(6,1) + t231 * rSges(6,2) + rSges(6,3) * t429;
t274 = pkin(3) * t427;
t248 = pkin(7) * t429 + t274;
t290 = sin(qJ(4));
t435 = t282 * t290;
t272 = pkin(4) * t435;
t319 = t281 * t427 - t296 * t429 + t272;
t187 = t319 - t248;
t416 = t166 + t187;
t431 = t283 * t285;
t229 = -t282 * t426 - t431;
t432 = t283 * t284;
t230 = t282 * t425 - t432;
t350 = -t230 * rSges(6,1) - t229 * rSges(6,2);
t434 = t282 * t291;
t165 = rSges(6,3) * t434 - t350;
t308 = -t291 * t463 - t484;
t430 = t283 * t290;
t398 = pkin(4) * t430;
t186 = t282 * t308 - t398;
t417 = t165 + t186;
t480 = -t282 * t417 - t283 * t416;
t478 = 2 * m(5);
t477 = 2 * m(6);
t476 = t282 ^ 2;
t475 = t283 ^ 2;
t474 = t282 / 0.2e1;
t473 = t283 / 0.2e1;
t472 = -t294 / 0.2e1;
t471 = -rSges(5,3) - pkin(7);
t470 = m(4) * t267;
t469 = pkin(3) * t294;
t468 = pkin(4) * t290;
t467 = pkin(7) * t291;
t466 = sin(qJ(1)) * pkin(1);
t465 = qJD(1) / 0.2e1;
t460 = rSges(6,3) * t291;
t459 = -rSges(6,3) + t296;
t287 = qJD(4) + qJD(5);
t369 = -t287 * t294 + qJD(1);
t403 = qJD(3) * t291;
t307 = t284 * t403 + t285 * t369;
t405 = qJD(1) * t294;
t368 = -t287 + t405;
t153 = t282 * t307 - t368 * t432;
t306 = t284 * t369 - t285 * t403;
t154 = t282 * t306 + t368 * t431;
t351 = t154 * rSges(6,1) + t153 * rSges(6,2);
t106 = rSges(6,3) * t310 + t351;
t377 = t283 * t402;
t458 = t106 * t429 + t165 * t377;
t455 = Icges(5,4) * t290;
t454 = Icges(5,4) * t293;
t453 = Icges(6,4) * t284;
t452 = Icges(6,4) * t285;
t423 = t290 * t294;
t428 = t283 * t293;
t239 = -t282 * t423 - t428;
t421 = t293 * t294;
t240 = t282 * t421 - t430;
t353 = -rSges(5,1) * t240 - rSges(5,2) * t239;
t181 = rSges(5,3) * t434 - t353;
t448 = t181 * t283;
t447 = t216 * t291;
t446 = t216 * t294;
t444 = t217 * t294;
t443 = t218 * t291;
t442 = t218 * t294;
t441 = t219 * t291;
t338 = Icges(6,1) * t285 - t453;
t235 = -Icges(6,5) * t294 + t291 * t338;
t439 = t235 * t285;
t339 = Icges(5,1) * t293 - t455;
t245 = -Icges(5,5) * t294 + t291 * t339;
t438 = t245 * t293;
t433 = t282 * t293;
t424 = t287 * t291;
t422 = t291 * t281;
t420 = t294 * t296;
t151 = t283 * t307 + t368 * t437;
t152 = t283 * t306 - t368 * t436;
t391 = t152 * rSges(6,1) + t151 * rSges(6,2) + rSges(6,3) * t377;
t105 = -rSges(6,3) * t383 + t391;
t262 = pkin(7) * t377;
t399 = qJD(4) * t294;
t367 = t399 * t468;
t400 = qJD(4) * t293;
t396 = pkin(4) * t400;
t386 = qJD(1) * t398 + t282 * t396 + t296 * t383;
t408 = qJD(1) * t282;
t419 = t105 - t262 + (t467 + t484) * t408 + (-t367 + (-t420 + t486) * qJD(3)) * t283 + t386;
t349 = rSges(6,1) * t285 - rSges(6,2) * t284;
t236 = -rSges(6,3) * t294 + t291 * t349;
t418 = t166 * t403 + t236 * t383;
t241 = -t283 * t423 + t433;
t242 = t283 * t421 + t435;
t182 = t242 * rSges(5,1) + t241 * rSges(5,2) + rSges(5,3) * t429;
t415 = -t182 - t248;
t194 = (-rSges(6,1) * t284 - rSges(6,2) * t285) * t424 + (t294 * t349 + t460) * qJD(3);
t401 = qJD(4) * t291;
t376 = t290 * t401;
t206 = -pkin(4) * t376 + qJD(3) * t308;
t414 = -t194 - t206;
t352 = rSges(5,1) * t293 - rSges(5,2) * t290;
t203 = (-rSges(5,1) * t290 - rSges(5,2) * t293) * t401 + (rSges(5,3) * t291 + t294 * t352) * qJD(3);
t358 = t467 + t469;
t256 = t358 * qJD(3);
t413 = -t203 - t256;
t135 = t294 * t165 + t236 * t434;
t247 = t358 * t282;
t412 = t282 * t247 + t283 * t248;
t228 = t485 - t486;
t411 = t228 + t236;
t246 = -rSges(5,3) * t294 + t291 * t352;
t271 = t291 * pkin(3) - t294 * pkin(7);
t410 = -t246 - t271;
t215 = Icges(4,3) * t282 + t283 * t333;
t409 = qJD(1) * t215;
t404 = qJD(3) * t282;
t364 = qJD(1) - t399;
t305 = t290 * t403 + t293 * t364;
t363 = -qJD(4) + t405;
t173 = t282 * t305 - t363 * t430;
t304 = t290 * t364 - t293 * t403;
t174 = t282 * t304 + t363 * t428;
t108 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t310;
t110 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t310;
t112 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t310;
t175 = Icges(5,5) * t240 + Icges(5,6) * t239 + Icges(5,3) * t434;
t177 = Icges(5,4) * t240 + Icges(5,2) * t239 + Icges(5,6) * t434;
t179 = Icges(5,1) * t240 + Icges(5,4) * t239 + Icges(5,5) * t434;
t328 = -t177 * t290 + t179 * t293;
t31 = (qJD(3) * t328 - t108) * t294 + (qJD(3) * t175 - t110 * t290 + t112 * t293 + (-t177 * t293 - t179 * t290) * qJD(4)) * t291;
t332 = Icges(5,5) * t293 - Icges(5,6) * t290;
t200 = (-Icges(5,5) * t290 - Icges(5,6) * t293) * t401 + (Icges(5,3) * t291 + t294 * t332) * qJD(3);
t335 = -Icges(5,2) * t290 + t454;
t201 = (-Icges(5,2) * t293 - t455) * t401 + (Icges(5,6) * t291 + t294 * t335) * qJD(3);
t202 = (-Icges(5,1) * t290 - t454) * t401 + (Icges(5,5) * t291 + t294 * t339) * qJD(3);
t243 = -Icges(5,3) * t294 + t291 * t332;
t244 = -Icges(5,6) * t294 + t291 * t335;
t62 = t173 * t244 + t174 * t245 + t200 * t434 + t201 * t239 + t202 * t240 + t243 * t310;
t395 = t31 / 0.2e1 + t62 / 0.2e1;
t171 = t283 * t305 + t363 * t435;
t172 = t283 * t304 - t363 * t433;
t309 = t377 - t383;
t107 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t309;
t109 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t309;
t111 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t309;
t176 = Icges(5,5) * t242 + Icges(5,6) * t241 + Icges(5,3) * t429;
t178 = Icges(5,4) * t242 + Icges(5,2) * t241 + Icges(5,6) * t429;
t180 = Icges(5,1) * t242 + Icges(5,4) * t241 + Icges(5,5) * t429;
t327 = -t178 * t290 + t180 * t293;
t32 = (qJD(3) * t327 - t107) * t294 + (qJD(3) * t176 - t109 * t290 + t111 * t293 + (-t178 * t293 - t180 * t290) * qJD(4)) * t291;
t61 = t171 * t244 + t172 * t245 + t200 * t429 + t201 * t241 + t202 * t242 + t243 * t309;
t394 = t32 / 0.2e1 + t61 / 0.2e1;
t127 = t239 * t244 + t240 * t245 + t243 * t434;
t94 = -t175 * t294 + t291 * t328;
t393 = t94 / 0.2e1 + t127 / 0.2e1;
t128 = t241 * t244 + t242 * t245 + t243 * t429;
t95 = -t176 * t294 + t291 * t327;
t392 = t95 / 0.2e1 + t128 / 0.2e1;
t390 = t172 * rSges(5,1) + t171 * rSges(5,2) + rSges(5,3) * t377;
t389 = -t256 + t414;
t261 = t282 * pkin(3) * t403;
t378 = t283 * t403;
t388 = t282 * (pkin(7) * t310 + qJD(1) * t274 - t261) + t283 * (-pkin(7) * t383 + t262 + (-t282 * t405 - t378) * pkin(3)) + t247 * t407;
t387 = -t271 - t411;
t385 = t283 * pkin(2) + t481;
t384 = t246 * t408;
t334 = -Icges(6,2) * t284 + t452;
t234 = -Icges(6,6) * t294 + t291 * t334;
t381 = t234 * t402;
t380 = t244 * t402;
t375 = t434 / 0.2e1;
t374 = t429 / 0.2e1;
t279 = t283 * pkin(6);
t371 = t279 - t466;
t370 = -t281 * t294 - pkin(2);
t205 = t410 * t283;
t366 = t283 * t396;
t365 = t294 * t106 + t194 * t434 + t310 * t236;
t146 = t387 * t283;
t331 = Icges(6,5) * t285 - Icges(6,6) * t284;
t233 = -Icges(6,3) * t294 + t291 * t331;
t138 = -t233 * t294 + (-t234 * t284 + t439) * t291;
t159 = Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t434;
t161 = Icges(6,4) * t230 + Icges(6,2) * t229 + Icges(6,6) * t434;
t163 = Icges(6,1) * t230 + Icges(6,4) * t229 + Icges(6,5) * t434;
t330 = -t161 * t284 + t163 * t285;
t85 = -t159 * t294 + t291 * t330;
t160 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t429;
t162 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t429;
t164 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t429;
t329 = -t162 * t284 + t164 * t285;
t86 = -t160 * t294 + t291 * t329;
t343 = t85 * t282 + t86 * t283;
t122 = t229 * t234 + t230 * t235 + t233 * t434;
t73 = t159 * t434 + t161 * t229 + t163 * t230;
t74 = t160 * t434 + t162 * t229 + t164 * t230;
t348 = t282 * t73 + t283 * t74;
t41 = -t122 * t294 + t291 * t348;
t123 = t231 * t234 + t232 * t235 + t233 * t429;
t75 = t159 * t429 + t161 * t231 + t163 * t232;
t76 = t160 * t429 + t162 * t231 + t164 * t232;
t347 = t282 * t75 + t283 * t76;
t42 = -t123 * t294 + t291 * t347;
t101 = Icges(6,4) * t154 + Icges(6,2) * t153 + Icges(6,6) * t310;
t103 = Icges(6,1) * t154 + Icges(6,4) * t153 + Icges(6,5) * t310;
t99 = Icges(6,5) * t154 + Icges(6,6) * t153 + Icges(6,3) * t310;
t19 = t101 * t231 + t103 * t232 + t151 * t161 + t152 * t163 + t159 * t309 + t429 * t99;
t100 = Icges(6,4) * t152 + Icges(6,2) * t151 + Icges(6,6) * t309;
t102 = Icges(6,1) * t152 + Icges(6,4) * t151 + Icges(6,5) * t309;
t98 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t309;
t20 = t100 * t231 + t102 * t232 + t151 * t162 + t152 * t164 + t160 * t309 + t429 * t98;
t53 = t282 * t76 - t283 * t75;
t183 = (-Icges(6,5) * t284 - Icges(6,6) * t285) * t424 + (Icges(6,3) * t291 + t294 * t331) * qJD(3);
t184 = (-Icges(6,2) * t285 - t453) * t424 + (Icges(6,6) * t291 + t294 * t334) * qJD(3);
t185 = (-Icges(6,1) * t284 - t452) * t424 + (Icges(6,5) * t291 + t294 * t338) * qJD(3);
t56 = t151 * t234 + t152 * t235 + t183 * t429 + t184 * t231 + t185 * t232 + t233 * t309;
t5 = (qJD(3) * t347 - t56) * t294 + (-qJD(1) * t53 + qJD(3) * t123 + t19 * t282 + t20 * t283) * t291;
t21 = t101 * t229 + t103 * t230 + t153 * t161 + t154 * t163 + t159 * t310 + t434 * t99;
t22 = t100 * t229 + t102 * t230 + t153 * t162 + t154 * t164 + t160 * t310 + t434 * t98;
t52 = t282 * t74 - t283 * t73;
t57 = t153 * t234 + t154 * t235 + t183 * t434 + t184 * t229 + t185 * t230 + t233 * t310;
t6 = (qJD(3) * t348 - t57) * t294 + (-qJD(1) * t52 + qJD(3) * t122 + t21 * t282 + t22 * t283) * t291;
t357 = t42 * t377 + t5 * t429 + t6 * t434 + (-t138 * t294 + t291 * t343) * t403 + t310 * t41;
t354 = rSges(5,1) * t174 + rSges(5,2) * t173;
t81 = t175 * t434 + t177 * t239 + t179 * t240;
t82 = t176 * t434 + t178 * t239 + t180 * t240;
t58 = t282 * t82 - t283 * t81;
t346 = t282 * t81 + t283 * t82;
t83 = t175 * t429 + t177 * t241 + t179 * t242;
t84 = t176 * t429 + t178 * t241 + t180 * t242;
t59 = t282 * t84 - t283 * t83;
t345 = t282 * t83 + t283 * t84;
t344 = t282 * t86 - t283 * t85;
t342 = t94 * t282 + t95 * t283;
t336 = Icges(4,2) * t294 + t457;
t326 = -t182 * t282 + t448;
t325 = -t181 * t282 - t182 * t283;
t322 = -t442 + t447;
t320 = -pkin(2) - t355;
t317 = t186 * t283 - t282 * t416;
t316 = t291 * t471 - pkin(2) - t469;
t314 = qJD(3) * t336;
t313 = qJD(3) * (-Icges(4,5) * t291 - Icges(4,6) * t294);
t311 = t291 * t459 + t370;
t12 = qJD(1) * t347 - t19 * t283 + t20 * t282;
t13 = qJD(1) * t348 - t21 * t283 + t22 * t282;
t25 = (qJD(3) * t330 - t99) * t294 + (qJD(3) * t159 + (-t161 * t287 + t103) * t285 + (-t163 * t287 - t101) * t284) * t291;
t26 = (qJD(3) * t329 - t98) * t294 + (qJD(3) * t160 + (-t162 * t287 + t102) * t285 + (-t164 * t287 - t100) * t284) * t291;
t303 = -t283 * t6 / 0.2e1 + t5 * t474 + t13 * t375 + t12 * t374 + (qJD(1) * t343 - t25 * t283 + t26 * t282) * t472 + t41 * t408 / 0.2e1 + t42 * t373 + t344 * t403 / 0.2e1 + t489 * t53 + t488 * t52;
t301 = -t294 * t183 + t233 * t403 + t402 * t439 + (t185 * t291 - t234 * t424) * t285;
t300 = -t294 * t200 + t243 * t403 + t402 * t438 + (t202 * t293 - t244 * t400) * t291;
t299 = t282 * t316 - t466;
t137 = t138 * t403;
t65 = (-t381 + (-t235 * t287 - t184) * t291) * t284 + t301;
t7 = t137 + (qJD(3) * t343 - t65) * t294 + (-qJD(1) * t344 + t25 * t282 + t26 * t283) * t291;
t298 = -t294 * t7 - t383 * t42 + t357;
t297 = t137 + (t25 + t57) * t375 + (t26 + t56) * t374 + (t123 + t86) * t489 + (t122 + t85) * t488;
t276 = pkin(6) * t407;
t254 = t355 * qJD(3);
t249 = t271 * t408;
t204 = t410 * t282;
t199 = t221 + t385;
t198 = t282 * t320 + t371 + t461;
t189 = t282 * t313 + t409;
t188 = -qJD(1) * t214 + t283 * t313;
t156 = t165 * t429;
t148 = t267 * t404 + (-t286 + (-rSges(4,3) - pkin(6)) * t282 + t320 * t283) * qJD(1);
t147 = t276 + (-t466 + (-pkin(2) - t462) * t282) * qJD(1) + t302;
t145 = t387 * t282;
t144 = -t243 * t294 + (-t244 * t290 + t438) * t291;
t143 = t144 * t403;
t142 = -t182 * t294 - t246 * t429;
t141 = t181 * t294 + t246 * t434;
t140 = t385 - t415;
t139 = t279 + t299 + t353;
t136 = -t166 * t294 - t236 * t429;
t134 = t215 * t282 - t483;
t133 = t214 * t282 - t283 * t322;
t132 = -t215 * t283 - t282 * t321;
t131 = -t214 * t283 - t282 * t322;
t130 = qJD(1) * t205 + t282 * t413;
t129 = t283 * t413 + t249 + t384;
t126 = t319 + t385 + t166;
t125 = t282 * t311 + t350 + t371 + t398;
t124 = t326 * t291;
t121 = -t366 + t261 + (-t367 + (-t422 - t485) * qJD(3)) * t282 + (t283 * t308 + t272) * qJD(1);
t119 = -t166 * t434 + t156;
t114 = rSges(5,3) * t310 + t354;
t113 = -rSges(5,3) * t383 + t390;
t104 = -t325 + t412;
t92 = -t294 * t416 - t411 * t429;
t91 = t186 * t294 + t228 * t434 + t135;
t90 = qJD(1) * t146 + t282 * t389;
t89 = t283 * t389 + t408 * t411 + t249;
t88 = t261 + t471 * t379 + (t283 * t316 - t481) * qJD(1) - t354;
t87 = -pkin(3) * t378 + qJD(1) * t299 + t262 + t276 + t390;
t72 = t291 * t317 + t156;
t71 = (-t380 + (-qJD(4) * t245 - t201) * t291) * t290 + t300;
t70 = t366 + (t367 + (t294 * t459 + t422) * qJD(3)) * t282 + (-t286 + (-pkin(6) - t468) * t282 + t311 * t283) * qJD(1) - t351;
t69 = t276 + (-t367 + (-t420 - t422) * qJD(3)) * t283 + (-t466 + (t370 - t460) * t282) * qJD(1) + t386 + t391;
t68 = t412 - t480;
t67 = (t246 * t404 + t114) * t294 + (-qJD(3) * t181 + t203 * t282 + t246 * t407) * t291;
t66 = (-qJD(3) * t246 * t283 - t113) * t294 + (qJD(3) * t182 - t203 * t283 + t384) * t291;
t64 = -t165 * t403 + t365;
t63 = -t105 * t294 + (-t194 * t291 - t236 * t402) * t283 + t418;
t46 = t326 * t402 + (qJD(1) * t325 - t113 * t282 + t114 * t283) * t291;
t45 = -t128 * t294 + t291 * t345;
t44 = -t127 * t294 + t291 * t346;
t43 = t113 * t283 + t114 * t282 + (t282 * t415 + t448) * qJD(1) + t388;
t40 = -t166 * t382 + (-t166 * t402 + (-qJD(1) * t165 - t105) * t291) * t282 + t458;
t34 = (t228 * t404 + t121) * t294 + (-qJD(3) * t417 + t206 * t282 + t228 * t407) * t291 + t365;
t33 = -t419 * t294 + (qJD(3) * t187 + t228 * t408) * t291 + (t291 * t414 - t402 * t411) * t283 + t418;
t30 = t107 * t434 + t109 * t239 + t111 * t240 + t173 * t178 + t174 * t180 + t176 * t310;
t29 = t108 * t434 + t110 * t239 + t112 * t240 + t173 * t177 + t174 * t179 + t175 * t310;
t28 = t107 * t429 + t109 * t241 + t111 * t242 + t171 * t178 + t172 * t180 + t176 * t309;
t27 = t108 * t429 + t110 * t241 + t112 * t242 + t171 * t177 + t172 * t179 + t175 * t309;
t18 = t419 * t283 + (t106 + t121) * t282 + (t417 * t283 + (-t248 - t416) * t282) * qJD(1) + t388;
t17 = t317 * t402 + (t480 * qJD(1) + t121 * t283 - t419 * t282) * t291 + t458;
t16 = qJD(1) * t346 + t282 * t30 - t283 * t29;
t15 = qJD(1) * t345 - t27 * t283 + t28 * t282;
t9 = (qJD(3) * t346 - t62) * t294 + (-qJD(1) * t58 + qJD(3) * t127 + t282 * t29 + t283 * t30) * t291;
t8 = (qJD(3) * t345 - t61) * t294 + (-t59 * qJD(1) + qJD(3) * t128 + t27 * t282 + t28 * t283) * t291;
t1 = [t300 + t301 + (t125 * t70 + t126 * t69) * t477 + (t139 * t88 + t140 * t87) * t478 + (t147 * t199 + t148 * t198) * t479 - t245 * t376 + (-t336 + t341) * t403 + (Icges(4,1) * t291 + t337 + t456) * t402 + (-t201 * t291 - t380) * t290 + (-t184 * t291 - t235 * t424 - t381) * t284; 0; 0; m(6) * (t125 * t89 + t126 * t90 + t145 * t69 + t146 * t70) + m(5) * (t129 * t139 + t130 * t140 + t204 * t87 + t205 * t88) + (t475 / 0.2e1 + t476 / 0.2e1) * t333 * qJD(3) + ((qJD(1) * t217 - t282 * t314) * t472 + t219 * t490 - t57 / 0.2e1 - t25 / 0.2e1 + (t447 / 0.2e1 - t442 / 0.2e1) * qJD(3) - t395 + m(4) * (-t148 * t267 - t198 * t254) + (t444 / 0.2e1 + t441 / 0.2e1 + t86 / 0.2e1 - t199 * t470 + t123 / 0.2e1 + t392) * qJD(1)) * t283 + ((-t216 * qJD(1) - t283 * t314) * t294 / 0.2e1 + t218 * t490 + t26 / 0.2e1 + t56 / 0.2e1 + (-t445 / 0.2e1 + t440 / 0.2e1) * qJD(3) + t394 + m(4) * (-t147 * t267 - t199 * t254) + (t446 / 0.2e1 + t443 / 0.2e1 + t198 * t470 + t122 / 0.2e1 + t85 / 0.2e1 + t393) * qJD(1)) * t282; m(4) * t93 + m(5) * t43 + m(6) * t18; (t145 * t90 + t146 * t89 + t18 * t68) * t477 + (t104 * t43 + t129 * t205 + t130 * t204) * t478 + (t475 + t476) * t267 * t254 * t479 + (t221 * t487 - t13 - t16 + (-t132 * qJD(1) + (-t322 * qJD(1) - t189) * t283) * t283) * t283 + (t12 + t15 + t220 * t487 + (t133 * qJD(1) + (t321 * qJD(1) + t188) * t282) * t282 + ((-t189 + (-t441 - t444) * qJD(3) + t217 * t402 + t219 * t403 - t409) * t282 + (t216 * t402 + t218 * t403 + t188 - (t443 + t446) * qJD(3)) * t283 + (t134 - t131 + (t215 - t322) * t282 + t483) * qJD(1)) * t283) * t282 + (-t131 * t283 + t132 * t282 + t52 + t58) * t408 + (-t133 * t283 + t134 * t282 + t53 + t59) * t407; (-t65 - t71 + (t282 * t393 + t283 * t392) * qJD(3)) * t294 + t297 + t143 + (t394 * t283 + t395 * t282 + (-t282 * t392 + t283 * t393) * qJD(1)) * t291 + m(6) * (t125 * t34 + t126 * t33 + t69 * t92 + t70 * t91) + m(5) * (t139 * t67 + t140 * t66 + t141 * t88 + t142 * t87); m(5) * t46 + m(6) * t17; t303 + ((qJD(1) * t94 + t32) * t472 + t8 / 0.2e1 + t58 * t372 + t44 * t465) * t282 + m(6) * (t145 * t33 + t146 * t34 + t17 * t68 + t18 * t72 + t89 * t91 + t90 * t92) + m(5) * (t104 * t46 + t124 * t43 + t129 * t141 + t130 * t142 + t204 * t66 + t205 * t67) + ((qJD(1) * t95 - t31) * t472 - t9 / 0.2e1 + t59 * t372 + t45 * t465) * t283 + (qJD(3) * (t282 * t95 - t94 * t283) / 0.2e1 + t15 * t473 + t16 * t474 + (t58 * t473 - t282 * t59 / 0.2e1) * qJD(1)) * t291; (t17 * t72 + t33 * t92 + t34 * t91) * t477 + (t124 * t46 + t141 * t67 + t142 * t66) * t478 + (t71 * t294 - t143 - t7 + (t282 * t44 + t283 * t45 - t294 * t342) * qJD(3)) * t294 + (t282 * t9 + t283 * t8 - t294 * (t31 * t282 + t32 * t283) + (-t144 * t294 + t291 * t342) * qJD(3) + ((-t294 * t94 + t44) * t283 + (t294 * t95 - t42 - t45) * t282) * qJD(1)) * t291 + t357; t297 + m(6) * (t125 * t64 + t126 * t63 + t135 * t70 + t136 * t69) - t65 * t294; m(6) * t40; t303 + m(6) * (t119 * t18 + t135 * t89 + t136 * t90 + t145 * t63 + t146 * t64 + t40 * t68); m(6) * (t119 * t17 + t135 * t34 + t136 * t33 + t40 * t72 + t63 * t92 + t64 * t91) + t298; (t119 * t40 + t135 * t64 + t136 * t63) * t477 + t298;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
