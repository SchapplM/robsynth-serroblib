% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR11_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:47
% DurationCPUTime: 11.00s
% Computational Cost: add. (15477->540), mult. (38578->813), div. (0->0), fcn. (41731->6), ass. (0->318)
t331 = cos(qJ(4));
t328 = sin(qJ(4));
t330 = cos(qJ(5));
t332 = cos(qJ(1));
t419 = t332 * t330;
t327 = sin(qJ(5));
t329 = sin(qJ(1));
t425 = t329 * t327;
t274 = t328 * t425 - t419;
t424 = t329 * t330;
t275 = t327 * t332 + t328 * t424;
t423 = t329 * t331;
t194 = -Icges(6,5) * t275 + Icges(6,6) * t274 + Icges(6,3) * t423;
t437 = t194 * t328;
t270 = Icges(6,4) * t275;
t196 = -Icges(6,2) * t274 - Icges(6,6) * t423 + t270;
t269 = Icges(6,4) * t274;
t200 = -Icges(6,1) * t275 + Icges(6,5) * t423 + t269;
t515 = t196 * t327 + t200 * t330;
t122 = t515 * t331 + t437;
t305 = rSges(5,1) * t331 - rSges(5,2) * t328;
t325 = t329 ^ 2;
t326 = t332 ^ 2;
t309 = t325 + t326;
t239 = t309 * t305;
t450 = rSges(6,2) * t327;
t372 = rSges(6,1) * t330 - t450;
t266 = rSges(6,3) * t328 + t331 * t372;
t308 = pkin(4) * t331 + pkin(7) * t328;
t394 = t266 + t308;
t219 = t394 * t329;
t221 = t394 * t332;
t501 = t219 * t329 + t221 * t332;
t513 = m(6) / 0.2e1;
t514 = m(5) / 0.2e1;
t412 = t239 * t514 + t501 * t513;
t208 = ((-rSges(6,3) - pkin(7)) * t328 + (-pkin(4) - t372) * t331) * t329;
t428 = t328 * t329;
t287 = rSges(5,1) * t423 - rSges(5,2) * t428;
t289 = t305 * t332;
t214 = -t329 * t287 - t289 * t332;
t420 = t331 * t332;
t427 = t328 * t332;
t238 = t331 * rSges(6,1) * t419 + rSges(6,3) * t427 - t420 * t450;
t386 = -pkin(4) * t420 - pkin(7) * t427 - t238;
t509 = t386 * t332;
t414 = (t329 * t208 + t509) * t513 + t214 * t514;
t52 = t414 - t412;
t521 = t52 * qJD(1);
t276 = -t327 * t427 - t424;
t277 = t328 * t419 - t425;
t195 = Icges(6,5) * t277 + Icges(6,6) * t276 - Icges(6,3) * t420;
t360 = t274 * t196 + t275 * t200;
t440 = Icges(6,4) * t277;
t198 = Icges(6,2) * t276 - Icges(6,6) * t420 + t440;
t271 = Icges(6,4) * t276;
t201 = Icges(6,1) * t277 - Icges(6,5) * t420 + t271;
t407 = t276 * t198 + t277 * t201;
t520 = -t407 - (-t194 * t329 - t195 * t332) * t331 - t360;
t408 = t276 * t196 - t277 * t200;
t409 = -t274 * t198 + t275 * t201;
t519 = -t331 * (t194 * t332 - t195 * t329) - t409 - t408;
t204 = t277 * rSges(6,1) + t276 * rSges(6,2) - rSges(6,3) * t420;
t176 = t328 * t204 + t266 * t420;
t202 = t275 * rSges(6,1) - t274 * rSges(6,2) - rSges(6,3) * t423;
t517 = -t202 * t328 - t266 * t423;
t120 = t176 * t329 + t332 * t517;
t518 = pkin(7) * t423 - t202;
t493 = -m(6) / 0.2e1;
t364 = Icges(6,5) * t330 - Icges(6,6) * t327;
t253 = Icges(6,3) * t328 + t331 * t364;
t438 = Icges(6,4) * t330;
t367 = -Icges(6,2) * t327 + t438;
t257 = Icges(6,6) * t328 + t331 * t367;
t439 = Icges(6,4) * t327;
t369 = Icges(6,1) * t330 - t439;
t261 = Icges(6,5) * t328 + t331 * t369;
t134 = t253 * t423 + t257 * t274 - t261 * t275;
t511 = t134 * t328;
t442 = Icges(5,4) * t328;
t368 = Icges(5,2) * t331 + t442;
t258 = Icges(5,6) * t332 + t368 * t329;
t441 = Icges(5,4) * t331;
t370 = Icges(5,1) * t328 + t441;
t262 = Icges(5,5) * t332 + t370 * t329;
t510 = (t258 * t331 + t262 * t328) * t329;
t452 = pkin(1) + qJ(3);
t454 = pkin(4) * t328;
t376 = t452 + t454;
t377 = (-pkin(6) + qJ(2)) * t329;
t502 = pkin(7) * t420 - t204;
t173 = t376 * t332 + t377 - t502;
t321 = t332 * qJ(2);
t499 = -pkin(6) * t332 - t376 * t329 + t321 + t518;
t115 = t329 * t173 + t332 * t499;
t481 = -t329 / 0.2e1;
t506 = t329 / 0.2e1;
t478 = -t332 / 0.2e1;
t477 = t332 / 0.2e1;
t252 = Icges(6,3) * t331 - t328 * t364;
t301 = Icges(5,1) * t331 - t442;
t421 = t330 * t261;
t430 = t327 * t257;
t503 = t328 * (t421 / 0.2e1 - t430 / 0.2e1 + t301 / 0.2e1 - t368 / 0.2e1 - t252 / 0.2e1);
t413 = (t208 * t332 - t329 * t386) * t513 + (-t287 * t332 + t329 * t289) * t514;
t281 = (-Icges(6,2) * t330 - t439) * t331;
t284 = (-Icges(6,1) * t327 - t438) * t331;
t497 = -(t261 / 0.2e1 + t281 / 0.2e1) * t327 + (t284 / 0.2e1 - t257 / 0.2e1) * t330;
t259 = -Icges(5,6) * t329 + t332 * t368;
t310 = Icges(5,4) * t420;
t263 = Icges(5,1) * t427 - Icges(5,5) * t329 + t310;
t299 = -Icges(5,2) * t328 + t441;
t282 = t299 * t329;
t283 = -Icges(5,2) * t427 + t310;
t285 = t301 * t329;
t286 = t301 * t332;
t496 = (-t329 * (t263 + t283) + (t262 + t282) * t332) * t331 + (t329 * (t259 - t286) + (-t258 + t285) * t332) * t328;
t495 = 0.2e1 * t309;
t494 = 0.4e1 * qJD(1);
t100 = t194 * t423 - t360;
t101 = -t195 * t423 + t409;
t363 = -t100 * t329 - t101 * t332;
t42 = t331 * t363 - t511;
t492 = -t42 / 0.2e1;
t223 = Icges(6,5) * t276 - Icges(6,6) * t277;
t403 = -Icges(6,2) * t277 + t201 + t271;
t405 = Icges(6,1) * t276 - t198 - t440;
t94 = t223 * t328 + (-t403 * t327 + t405 * t330) * t331;
t491 = -t94 / 0.2e1;
t264 = rSges(6,3) * t331 - t328 * t372;
t124 = (-t264 * t329 - t202) * t331;
t125 = (t264 * t332 + t204) * t331 + (-t266 * t332 + t238) * t328;
t489 = m(6) * (t124 * t499 + t125 * t173 - t176 * t386 + t208 * t517);
t228 = -rSges(6,1) * t274 - rSges(6,2) * t275;
t229 = rSges(6,1) * t276 - rSges(6,2) * t277;
t288 = (-rSges(6,1) * t327 - rSges(6,2) * t330) * t331;
t487 = m(6) * (t115 * t288 + t219 * t229 - t221 * t228);
t484 = m(6) * (-t173 * t386 + t208 * t499);
t483 = t309 / 0.2e1;
t482 = t328 / 0.2e1;
t480 = -t329 / 0.4e1;
t479 = t331 / 0.2e1;
t476 = t332 / 0.4e1;
t475 = m(3) * ((rSges(3,3) * t332 + t321) * t332 + (rSges(3,3) + qJ(2)) * t325);
t390 = rSges(4,3) + t452;
t249 = rSges(4,2) * t332 - t390 * t329 + t321;
t250 = (rSges(4,2) + qJ(2)) * t329 + t390 * t332;
t474 = m(4) * (-t249 * t329 + t332 * t250);
t473 = m(4) * (t249 * t332 + t329 * t250);
t373 = rSges(5,1) * t328 + rSges(5,2) * t331;
t341 = t373 + t452;
t216 = t321 + (-rSges(5,3) - pkin(6)) * t332 - t341 * t329;
t322 = t329 * rSges(5,3);
t217 = t332 * t341 - t322 + t377;
t472 = m(5) * (-t216 * t287 + t217 * t289);
t471 = m(5) * (-t216 * t329 + t332 * t217);
t470 = m(5) * (t216 * t332 + t329 * t217);
t466 = m(6) * (t332 * t173 - t329 * t499);
t465 = m(6) * t115;
t464 = m(6) * (t176 * t332 - t329 * t517);
t463 = m(6) * t120;
t458 = m(6) * (t219 * t332 - t221 * t329);
t356 = -t329 * t228 - t229 * t332;
t456 = m(6) * t356;
t159 = -t228 * t332 + t329 * t229;
t455 = m(6) * t159;
t453 = qJD(4) / 0.2e1;
t451 = m(6) * qJD(5);
t448 = t329 * t42;
t135 = -t253 * t420 + t276 * t257 + t277 * t261;
t132 = t135 * t328;
t102 = t194 * t420 + t408;
t103 = -t195 * t420 + t407;
t362 = -t329 * t102 - t103 * t332;
t43 = t331 * t362 + t132;
t447 = t332 * t43;
t446 = -t100 + t520;
t445 = t101 + t519;
t444 = -t102 - t519;
t443 = t103 + t520;
t436 = t195 * t328;
t432 = t253 * t328;
t237 = t266 * t329;
t256 = Icges(6,6) * t331 - t328 * t367;
t431 = t327 * t256;
t278 = (-Icges(6,5) * t327 - Icges(6,6) * t330) * t331;
t429 = t328 * t278;
t365 = Icges(5,5) * t328 + Icges(5,6) * t331;
t254 = Icges(5,3) * t332 + t365 * t329;
t426 = t329 * t254;
t260 = Icges(6,5) * t331 - t328 * t369;
t422 = t330 * t260;
t342 = m(6) * (t124 * t332 + t125 * t329);
t352 = m(6) * t288 * t483;
t63 = -t342 / 0.2e1 + t352;
t418 = t63 * qJD(3);
t85 = t124 * t329 - t125 * t332;
t417 = t85 * qJD(2);
t416 = t85 * qJD(5);
t406 = -Icges(6,1) * t274 - t196 - t270;
t404 = -Icges(6,2) * t275 - t200 - t269;
t401 = t258 * t420 + t262 * t427;
t399 = -t257 + t284;
t397 = t261 + t281;
t395 = pkin(7) * t331 + t264 - t454;
t393 = qJD(1) * t331;
t392 = qJD(5) * t331;
t178 = (t493 - m(5) / 0.2e1 - m(4) / 0.2e1) * t495;
t391 = t178 * qJD(1);
t12 = t511 + (t443 * t329 + t444 * t332) * t331;
t389 = -t12 / 0.2e1 + t492;
t13 = t132 + (t445 * t329 + t446 * t332) * t331;
t388 = t43 / 0.2e1 - t13 / 0.2e1;
t164 = t332 * (-Icges(5,3) * t329 + t332 * t365) + t259 * t423 + t263 * t428;
t383 = -t423 / 0.4e1;
t382 = -t420 / 0.4e1;
t109 = -t397 * t274 + t399 * t275 - t278 * t423;
t110 = t397 * t276 + t399 * t277 - t278 * t420;
t222 = -Icges(6,5) * t274 - Icges(6,6) * t275;
t93 = t222 * t328 + (-t404 * t327 + t406 * t330) * t331;
t374 = t487 / 0.2e1 + (t110 + t94) * t480 + (t109 + t93) * t476;
t366 = Icges(5,5) * t331 - Icges(5,6) * t328;
t358 = -t198 * t327 + t201 * t330;
t123 = t331 * t358 + t436;
t361 = t122 * t329 - t123 * t332;
t357 = t202 * t332 - t204 * t329;
t355 = t421 - t430;
t353 = m(6) * (t173 * t229 - t228 * t499) + t429 / 0.2e1;
t163 = t332 * t254 + t510;
t166 = t401 - t426;
t18 = t446 * t329 - t445 * t332;
t50 = t102 * t332 - t103 * t329;
t351 = t166 * t478 + t401 * t477 + (-t163 + t510) * t506 - t50 / 0.2e1 + t18 / 0.2e1;
t17 = t444 * t329 - t443 * t332;
t49 = t100 * t332 - t101 * t329;
t350 = (t164 - t166 - t426) * t481 + t326 * t254 / 0.2e1 + t163 * t478 + t164 * t506 - t17 / 0.2e1 - t49 / 0.2e1;
t81 = -t222 * t423 - t404 * t274 + t406 * t275;
t82 = -t223 * t423 - t403 * t274 + t405 * t275;
t35 = -t82 * t329 + t332 * t81;
t83 = -t222 * t420 + t404 * t276 + t406 * t277;
t84 = -t223 * t420 + t403 * t276 + t405 * t277;
t36 = -t84 * t329 + t332 * t83;
t349 = t35 * t477 + t36 * t481;
t346 = t253 * t329 + t515;
t345 = t253 * t332 - t358;
t344 = t252 - t355;
t340 = t12 * t480 + t13 * t476 + t18 * t383 - t448 / 0.4e1 - t447 / 0.4e1 + t50 * t423 / 0.4e1 + (t17 + t49) * t382;
t118 = (t253 + t422 - t431) * t331 + t344 * t328;
t165 = t331 * t355 + t432;
t233 = t257 * t329;
t235 = t261 * t329;
t88 = (-t233 * t327 + t235 * t330 - t194) * t331 + t346 * t328;
t234 = t257 * t332;
t236 = t261 * t332;
t89 = (-t234 * t327 + t236 * t330 + t195) * t331 + t345 * t328;
t336 = -t331 * t344 + t432;
t98 = -t256 * t274 + t260 * t275 + t336 * t329;
t99 = t276 * t256 + t277 * t260 + t332 * t336;
t339 = t118 * t482 + t165 * t479 + t489 / 0.2e1 + (-t122 - t134) * t428 / 0.4e1 + (t123 + t135) * t427 / 0.4e1 + (t88 + t98) * t383 + (t89 + t99) * t382;
t338 = -t331 * t346 - t437;
t337 = -t331 * t345 + t436;
t335 = t422 / 0.2e1 - t431 / 0.2e1 + t253 / 0.2e1 - t370 / 0.2e1 - t299 / 0.2e1;
t280 = t332 * t366;
t279 = t366 * t329;
t220 = t395 * t332;
t218 = t395 * t329;
t180 = t328 * t229 + t288 * t420;
t179 = -t228 * t328 - t288 * t423;
t177 = (-m(6) / 0.4e1 - m(5) / 0.4e1 - m(4) / 0.4e1) * t495 + (m(4) + m(5) + m(6)) * t483;
t149 = t455 / 0.2e1;
t148 = t456 / 0.2e1;
t147 = t159 * t331;
t145 = t458 / 0.2e1;
t138 = t357 * t331;
t136 = t509 + (-t308 * t329 - t237) * t329;
t128 = (t429 + (-t397 * t327 + t399 * t330) * t331) * t328;
t126 = (-pkin(4) * t427 + t502) * t332 + (-pkin(4) * t428 + t518) * t329;
t112 = t463 / 0.2e1;
t111 = t464 / 0.2e1;
t104 = (-t237 * t332 + t238 * t329) * t331 + t357 * t328;
t79 = m(6) * t85 * t453;
t78 = t276 * t234 + t277 * t236 + t332 * t337;
t77 = t276 * t233 + t277 * t235 + t332 * t338;
t76 = -t234 * t274 + t236 * t275 + t337 * t329;
t75 = -t233 * t274 + t235 * t275 + t338 * t329;
t69 = t466 + t471 + t474;
t68 = t145 - t413;
t67 = t145 + t413;
t66 = -t458 / 0.2e1 + t413;
t62 = t342 / 0.2e1 + t352;
t61 = t126 * t356 + t501 * t288;
t53 = t412 + t414;
t51 = t331 * t497 + t353;
t46 = t465 + t470 + t473 + t475;
t45 = t165 * t328 + t331 * t361;
t34 = t112 - t456 / 0.2e1;
t33 = t111 - t455 / 0.2e1;
t32 = t149 + t111;
t31 = t149 - t464 / 0.2e1;
t30 = t148 + t112;
t29 = t148 - t463 / 0.2e1;
t28 = -t78 * t329 + t332 * t77;
t27 = -t76 * t329 + t332 * t75;
t26 = t331 * t335 + t472 + t484 - t503;
t23 = -t104 * t138 + t124 * t517 + t125 * t176;
t22 = t110 * t328 + (-t329 * t83 - t332 * t84) * t331;
t21 = t109 * t328 + (-t329 * t81 - t332 * t82) * t331;
t14 = (-t88 * t329 - t89 * t332 + t165) * t331 + (t118 - t361) * t328;
t9 = (-t329 * t77 - t332 * t78 + t135) * t331 + (-t362 + t99) * t328;
t8 = (-t329 * t75 - t332 * t76 - t134) * t331 + (-t363 + t98) * t328;
t7 = m(6) * t61 + t349;
t6 = (t388 * t329 + t389 * t332) * t331;
t5 = t350 * t329 + t332 * t351;
t4 = m(6) * t23 + (t9 * t478 + t8 * t481 + t45 / 0.2e1) * t331 + (t447 / 0.2e1 + t448 / 0.2e1 + t14 / 0.2e1) * t328;
t3 = t339 + (t110 / 0.4e1 + t94 / 0.4e1) * t329 + (-t109 / 0.4e1 - t93 / 0.4e1) * t332 - t487 / 0.2e1 + t340;
t2 = t339 + (t43 / 0.4e1 - t13 / 0.4e1 + (t49 / 0.4e1 + t17 / 0.4e1) * t331) * t332 + (t42 / 0.4e1 + t12 / 0.4e1 + (-t50 / 0.4e1 + t18 / 0.4e1) * t331) * t329 + t374;
t1 = -t489 / 0.2e1 + (-t165 / 0.2e1 + (t99 / 0.4e1 + t89 / 0.4e1) * t332 + (t88 / 0.4e1 + t98 / 0.4e1) * t329) * t331 + (-t118 / 0.2e1 + (-t135 / 0.4e1 - t123 / 0.4e1) * t332 + (t134 / 0.4e1 + t122 / 0.4e1) * t329) * t328 + t340 + t374;
t10 = [t46 * qJD(2) + t69 * qJD(3) + t26 * qJD(4) + t51 * qJD(5), qJD(1) * t46 + qJD(3) * t177 + qJD(4) * t53 + qJD(5) * t30, qJD(1) * t69 + qJD(2) * t177 + qJD(4) * t67 + qJD(5) * t32, t26 * qJD(1) + t53 * qJD(2) + t67 * qJD(3) + t2 * qJD(5) + (m(6) * (t173 * t218 + t208 * t221 - t219 * t386 + t220 * t499) + (t88 / 0.2e1 + t98 / 0.2e1 + m(5) * (-t216 * t373 - t287 * t305) - t365 * t477 + (-t258 / 0.2e1 + t285 / 0.2e1) * t331 + (-t262 / 0.2e1 - t282 / 0.2e1) * t328 - t351) * t332 + (-t99 / 0.2e1 - t89 / 0.2e1 + m(5) * (-t217 * t373 + t289 * t305) + t365 * t481 + (t259 / 0.2e1 - t286 / 0.2e1) * t331 + (t263 / 0.2e1 + t283 / 0.2e1) * t328 - t350) * t329) * qJD(4), t51 * qJD(1) + t30 * qJD(2) + t32 * qJD(3) + t2 * qJD(4) + (m(6) * (t173 * t180 + t176 * t229 + t179 * t499 - t228 * t517) + t128) * qJD(5) + ((-t110 / 0.2e1 + t491 - t389) * t332 + (-t93 / 0.2e1 - t109 / 0.2e1 - t388) * t329) * t392; t178 * qJD(3) + t52 * qJD(4) + t29 * qJD(5) + (-t465 / 0.4e1 - t470 / 0.4e1 - t475 / 0.4e1 - t473 / 0.4e1) * t494, 0, t391, t521 + 0.2e1 * ((-t218 * t332 + t220 * t329) * t453 + t416 / 0.4e1) * m(6), t29 * qJD(1) + t79 + (t179 * t329 - t180 * t332) * t451; -t178 * qJD(2) + t66 * qJD(4) + t31 * qJD(5) + (-t466 / 0.4e1 - t471 / 0.4e1 - t474 / 0.4e1) * t494, -t391, 0, t66 * qJD(1) + 0.2e1 * ((t218 * t329 + t220 * t332) * t513 - m(5) * t373 * t483) * qJD(4) + t62 * qJD(5), t31 * qJD(1) + t62 * qJD(4) + (t179 * t332 + t180 * t329) * t451; -t52 * qJD(2) + t68 * qJD(3) + t5 * qJD(4) + t1 * qJD(5) + (-t484 / 0.4e1 - t472 / 0.4e1) * t494 - t335 * t393 + t503 * qJD(1), t416 * t493 - t521, qJD(1) * t68 + qJD(5) * t63, t5 * qJD(1) + (m(5) * ((-t329 * (rSges(5,3) * t332 + t373 * t329) + (-t332 * t373 + t322) * t332) * t214 - t239 * t373) + m(6) * (t126 * t136 + t218 * t219 + t220 * t221) + (t325 * t280 + (-t329 * t279 + t496) * t332 + t28) * t481 + (t326 * t279 + (-t332 * t280 + t496) * t329 + t27) * t477) * qJD(4) + t7 * qJD(5), t1 * qJD(1) + t418 + t7 * qJD(4) + t417 * t493 + (-t45 / 0.2e1 + (-t36 / 0.2e1 + t9 / 0.2e1) * t332 + (-t35 / 0.2e1 + t8 / 0.2e1) * t329) * t392 + (t21 * t477 + t22 * t481 + (-t14 / 0.2e1 + (t93 / 0.2e1 - t43 / 0.2e1) * t332 + (t491 + t492) * t329) * t328 + (t120 * t288 + t147 * t126 - t138 * t356 + t179 * t221 + t180 * t219 - t23) * m(6)) * qJD(5); -t353 * qJD(1) + t34 * qJD(2) + t33 * qJD(3) + t3 * qJD(4) + t6 * qJD(5) - t497 * t393, qJD(1) * t34 + t79, qJD(1) * t33 - qJD(4) * t63, t3 * qJD(1) - t418 + (t9 * t481 + t50 * t427 / 0.2e1 - t28 * t420 / 0.2e1 + t8 * t477 + t49 * t428 / 0.2e1 - t27 * t423 / 0.2e1 + (-t122 * t332 - t123 * t329) * t479 + (-t89 * t329 + t88 * t332) * t482 - t349) * qJD(4) + t4 * qJD(5) + (t417 / 0.2e1 + (t104 * t126 + t124 * t221 + t125 * t219 - t136 * t138 + t176 * t218 + t220 * t517 - t61) * qJD(4)) * m(6), t6 * qJD(1) + t4 * qJD(4) + (m(6) * (-t138 * t147 + t176 * t180 + t179 * t517) + t128 * t482 + (t22 * t478 + t21 * t481 + (-t329 * t93 - t332 * t94) * t482) * t331) * qJD(5);];
Cq = t10;
