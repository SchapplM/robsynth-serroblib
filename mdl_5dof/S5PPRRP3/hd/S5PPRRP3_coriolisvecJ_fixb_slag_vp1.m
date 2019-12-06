% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:11:08
% DurationCPUTime: 29.45s
% Computational Cost: add. (14417->673), mult. (40513->995), div. (0->0), fcn. (44596->8), ass. (0->270)
t491 = Icges(5,4) - Icges(6,5);
t477 = Icges(5,1) + Icges(6,1);
t476 = Icges(6,4) + Icges(5,5);
t475 = Icges(5,2) + Icges(6,3);
t474 = Icges(6,2) + Icges(5,3);
t490 = Icges(5,6) - Icges(6,6);
t282 = cos(pkin(8));
t281 = sin(pkin(7));
t286 = cos(qJ(3));
t384 = t281 * t286;
t348 = t282 * t384;
t283 = cos(pkin(7));
t285 = sin(qJ(3));
t383 = t283 * t285;
t260 = t348 - t383;
t280 = sin(pkin(8));
t284 = sin(qJ(4));
t388 = t280 * t284;
t407 = cos(qJ(4));
t222 = t260 * t407 + t281 * t388;
t505 = t491 * t222;
t382 = t283 * t286;
t385 = t281 * t285;
t262 = t282 * t382 + t385;
t224 = t262 * t407 + t283 * t388;
t504 = t491 * t224;
t345 = t280 * t407;
t302 = -t262 * t284 + t283 * t345;
t503 = t491 * t302;
t303 = -t260 * t284 + t281 * t345;
t502 = t491 * t303;
t259 = t282 * t385 + t382;
t493 = -t490 * t259 - t475 * t303 - t505;
t261 = t282 * t383 - t384;
t492 = -t490 * t261 - t475 * t302 - t504;
t432 = t476 * t222 + t474 * t259 + t490 * t303;
t431 = t476 * t224 + t474 * t261 + t490 * t302;
t489 = t477 * t222 + t476 * t259 + t502;
t488 = t477 * t224 + t476 * t261 + t503;
t332 = t286 * t345;
t264 = -t282 * t284 + t332;
t501 = t491 * t264;
t386 = t280 * t286;
t263 = t282 * t407 + t284 * t386;
t500 = t491 * t263;
t441 = t489 * t222 + t432 * t259 - t493 * t303;
t440 = t488 * t222 + t431 * t259 - t492 * t303;
t439 = t489 * t224 + t432 * t261 - t493 * t302;
t438 = t488 * t224 + t431 * t261 - t492 * t302;
t387 = t280 * t285;
t437 = t493 * t263 + t489 * t264 + t432 * t387;
t436 = t492 * t263 + t488 * t264 + t431 * t387;
t246 = t259 * qJD(3);
t155 = qJD(4) * t222 - t246 * t284;
t156 = qJD(4) * t303 - t246 * t407;
t357 = qJD(3) * t285;
t247 = qJD(3) * t348 - t283 * t357;
t499 = -t155 * t490 + t156 * t476 + t247 * t474;
t248 = t261 * qJD(3);
t157 = qJD(4) * t224 - t248 * t284;
t158 = qJD(4) * t302 - t248 * t407;
t249 = t262 * qJD(3);
t498 = -t157 * t490 + t158 * t476 + t249 * t474;
t497 = t155 * t475 - t156 * t491 - t247 * t490;
t496 = t157 * t475 - t158 * t491 - t249 * t490;
t495 = -t491 * t155 + t156 * t477 + t476 * t247;
t494 = -t491 * t157 + t158 * t477 + t476 * t249;
t342 = t280 * t357;
t225 = -qJD(4) * t332 + (qJD(4) * t282 + t342) * t284;
t226 = -qJD(4) * t263 - t342 * t407;
t356 = qJD(3) * t286;
t341 = t280 * t356;
t487 = t225 * t490 + t226 * t476 + t341 * t474;
t486 = -t225 * t475 - t226 * t491 - t341 * t490;
t485 = t491 * t225 + t226 * t477 + t476 * t341;
t426 = -t263 * t490 + t264 * t476 + t387 * t474;
t470 = -t263 * t475 + t387 * t490 + t501;
t484 = t264 * t477 + t476 * t387 - t500;
t359 = qJD(3) * t280;
t344 = t281 * t359;
t235 = qJD(4) * t259 + t344;
t418 = t235 / 0.2e1;
t343 = t283 * t359;
t236 = qJD(4) * t261 + t343;
t416 = t236 / 0.2e1;
t358 = qJD(3) * t282;
t271 = qJD(4) * t387 - t358;
t410 = t271 / 0.2e1;
t483 = qJD(4) / 0.2e1;
t451 = t155 * t493 + t156 * t489 + t222 * t495 + t247 * t432 + t259 * t499 - t303 * t497;
t450 = t155 * t492 + t156 * t488 + t222 * t494 + t247 * t431 + t259 * t498 - t303 * t496;
t449 = t157 * t493 + t158 * t489 + t224 * t495 + t249 * t432 + t261 * t499 - t302 * t497;
t448 = t157 * t492 + t158 * t488 + t224 * t494 + t249 * t431 + t261 * t498 - t302 * t496;
t446 = (t285 * t499 + t356 * t432) * t280 + t495 * t264 + t497 * t263 + t489 * t226 - t493 * t225;
t445 = (t285 * t498 + t356 * t431) * t280 + t494 * t264 + t496 * t263 + t488 * t226 - t492 * t225;
t482 = -t155 * t470 + t156 * t484 + t222 * t485 + t247 * t426 + t259 * t487 - t303 * t486;
t481 = -t157 * t470 + t158 * t484 + t224 * t485 + t249 * t426 + t261 * t487 - t302 * t486;
t480 = (t285 * t487 + t356 * t426) * t280 + t485 * t264 + t486 * t263 + t484 * t226 + t470 * t225;
t435 = -t222 * t484 - t259 * t426 - t303 * t470;
t434 = t224 * t484 + t261 * t426 + t302 * t470;
t433 = t263 * t470 - t264 * t484 - t387 * t426;
t468 = -t284 * t475 + t407 * t491;
t467 = t284 * t490 - t407 * t476;
t466 = t284 * t491 - t407 * t477;
t465 = t247 * t437 + t249 * t436;
t464 = t247 * t439 + t249 * t438;
t463 = t247 * t441 + t440 * t249;
t462 = t481 * t410 + t448 * t416 + t449 * t418 + (t434 * t341 + t464) * t483;
t461 = t482 * t410 + t450 * t416 + t451 * t418 + (-t435 * t341 + t463) * t483;
t460 = rSges(6,1) + pkin(4);
t459 = rSges(6,3) + qJ(5);
t312 = -pkin(4) * t407 - qJ(5) * t284;
t313 = -rSges(6,1) * t407 - rSges(6,3) * t284;
t458 = t312 + t313;
t455 = t480 * t271 + t445 * t236 + t446 * t235 + (-t433 * t341 + t465) * qJD(4);
t381 = rSges(6,2) * t259 + t222 * t460 - t303 * t459;
t367 = rSges(6,2) * t387 + t263 * t459 + t264 * t460;
t190 = -pkin(3) * t246 + pkin(6) * t247;
t269 = (-pkin(3) * t285 + pkin(6) * t286) * t280;
t256 = qJD(3) * t269;
t373 = t190 * t358 + t256 * t344;
t351 = qJD(5) * t263;
t378 = rSges(6,2) * t341 - t225 * t459 + t226 * t460 + t351;
t353 = qJD(5) * t303;
t405 = rSges(6,2) * t247 + t155 * t459 + t156 * t460 - t353;
t8 = qJD(5) * t157 - t405 * t271 + t378 * t235 + (t247 * t367 - t341 * t381) * qJD(4) + t373;
t454 = t381 * t8;
t453 = t441 * t235 + t440 * t236 - t435 * t271;
t452 = t439 * t235 + t438 * t236 + t434 * t271;
t447 = t437 * t235 + t436 * t236 - t433 * t271;
t390 = t280 * t281;
t389 = t280 * t283;
t430 = t468 * t259 - t260 * t490;
t429 = t468 * t261 - t262 * t490;
t428 = t466 * t259 + t476 * t260;
t427 = t466 * t261 + t476 * t262;
t425 = (t468 * t285 - t286 * t490) * t280;
t424 = (t466 * t285 + t476 * t286) * t280;
t423 = (-t484 * t407 + t470 * t284 + (t467 * t285 + t474 * t286) * t280) * t271 + (t467 * t261 + t474 * t262 - t284 * t492 - t407 * t488) * t236 + (t467 * t259 + t474 * t260 - t284 * t493 - t407 * t489) * t235;
t422 = (t475 * t264 - t484 + t500) * t271 + (t475 * t224 - t488 - t503) * t236 + (t475 * t222 - t489 - t502) * t235;
t421 = (-t477 * t263 - t470 - t501) * t271 + (t477 * t302 + t492 - t504) * t236 + (t477 * t303 + t493 - t505) * t235;
t420 = (-t476 * t263 - t264 * t490) * t271 + (-t224 * t490 + t476 * t302) * t236 + (-t222 * t490 + t476 * t303) * t235;
t265 = (-Icges(4,5) * t285 - Icges(4,6) * t286) * t280;
t250 = qJD(3) * t265;
t419 = -t235 / 0.2e1;
t417 = -t236 / 0.2e1;
t415 = t247 / 0.2e1;
t411 = -t271 / 0.2e1;
t352 = qJD(5) * t302;
t404 = rSges(6,2) * t249 + t157 * t459 + t158 * t460 - t352;
t191 = -pkin(3) * t248 + pkin(6) * t249;
t86 = rSges(5,1) * t158 - rSges(5,2) * t157 + rSges(5,3) * t249;
t403 = -t191 - t86;
t402 = Icges(4,4) * t260;
t401 = Icges(4,4) * t262;
t400 = Icges(4,4) * t285;
t399 = Icges(4,4) * t286;
t392 = t191 * t282;
t391 = t256 * t283;
t380 = rSges(6,2) * t261 + t224 * t460 - t302 * t459;
t107 = rSges(5,1) * t224 + rSges(5,2) * t302 + rSges(5,3) * t261;
t206 = pkin(3) * t262 + pkin(6) * t261;
t379 = -t107 - t206;
t377 = t222 * t459 + t303 * t460;
t376 = t224 * t459 + t302 * t460;
t375 = t260 * rSges(6,2) + t259 * t458;
t374 = t262 * rSges(6,2) + t261 * t458;
t372 = t282 * t190 + t256 * t390;
t163 = -Icges(4,2) * t259 + Icges(4,6) * t390 + t402;
t371 = -Icges(4,1) * t259 - t163 - t402;
t164 = -Icges(4,2) * t261 + Icges(4,6) * t389 + t401;
t370 = -Icges(4,1) * t261 - t164 - t401;
t254 = Icges(4,4) * t259;
t165 = Icges(4,1) * t260 + Icges(4,5) * t390 - t254;
t369 = Icges(4,2) * t260 - t165 + t254;
t255 = Icges(4,4) * t261;
t166 = Icges(4,1) * t262 + Icges(4,5) * t389 - t255;
t368 = Icges(4,2) * t262 - t166 + t255;
t203 = -pkin(3) * t259 + pkin(6) * t260;
t366 = t203 * t358 + t269 * t344;
t204 = pkin(3) * t260 + pkin(6) * t259;
t270 = (pkin(3) * t286 + pkin(6) * t285) * t280;
t365 = t282 * t204 + t270 * t390;
t364 = -t263 * t460 + t264 * t459;
t363 = (rSges(6,2) * t286 + t285 * t313) * t280 + t312 * t387;
t242 = -Icges(4,6) * t282 + (-Icges(4,2) * t285 + t399) * t280;
t267 = (-Icges(4,1) * t285 - t399) * t280;
t362 = t242 - t267;
t243 = -Icges(4,5) * t282 + (Icges(4,1) * t286 - t400) * t280;
t266 = (-Icges(4,2) * t286 - t400) * t280;
t361 = t243 + t266;
t360 = qJD(2) * t283;
t355 = qJD(4) * t249;
t354 = qJD(4) * t286;
t350 = qJD(5) * t284;
t349 = -t191 - t404;
t347 = -t206 - t380;
t279 = qJD(2) * t281;
t346 = t204 * t358 + t270 * t344 + t279;
t337 = qJD(4) * t415;
t336 = t355 / 0.2e1;
t331 = t341 / 0.2e1;
t105 = rSges(5,1) * t222 + rSges(5,2) * t303 + rSges(5,3) * t259;
t118 = t226 * rSges(5,1) + t225 * rSges(5,2) + rSges(5,3) * t341;
t178 = rSges(5,1) * t264 - rSges(5,2) * t263 + rSges(5,3) * t387;
t84 = rSges(5,1) * t156 - rSges(5,2) * t155 + rSges(5,3) * t247;
t34 = t235 * t118 - t271 * t84 + (-t105 * t341 + t247 * t178) * qJD(4) + t373;
t59 = -t105 * t271 + t178 * t235 + t346;
t329 = t34 * t105 + t59 * t84;
t167 = rSges(4,1) * t260 - rSges(4,2) * t259 + rSges(4,3) * t390;
t245 = -t282 * rSges(4,3) + (rSges(4,1) * t286 - rSges(4,2) * t285) * t280;
t90 = t279 + (t167 * t282 + t245 * t390) * qJD(3);
t168 = rSges(4,1) * t262 - rSges(4,2) * t261 + rSges(4,3) * t389;
t91 = -t360 + (-t168 * t282 - t245 * t389) * qJD(3);
t322 = t281 * t90 - t283 * t91;
t321 = t105 * t249 - t107 * t247;
t186 = -rSges(4,1) * t246 - rSges(4,2) * t247;
t305 = (-rSges(4,1) * t285 - rSges(4,2) * t286) * t280;
t253 = qJD(3) * t305;
t109 = (t186 * t282 + t253 * t390) * qJD(3);
t187 = -rSges(4,1) * t248 - rSges(4,2) * t249;
t110 = (-t187 * t282 - t253 * t389) * qJD(3);
t320 = t109 * t281 - t110 * t283;
t319 = t167 * t283 - t168 * t281;
t318 = t186 * t283 - t187 * t281;
t317 = qJD(4) * t331;
t316 = t190 * t343 - t191 * t344;
t205 = -pkin(3) * t261 + pkin(6) * t262;
t315 = t203 * t343 - t205 * t344;
t314 = -rSges(5,1) * t407 + rSges(5,2) * t284;
t304 = t204 * t343 - t206 * t344 + qJD(1);
t33 = -t235 * t380 + t236 * t381 + t304 + t351;
t7 = -qJD(5) * t225 + t405 * t236 - t404 * t235 + (-t247 * t380 + t249 * t381) * qJD(4) + t316;
t301 = t33 * t405 + t381 * t7;
t38 = t235 * t367 - t271 * t381 + t346 - t352;
t298 = t367 * t8 + t378 * t38;
t297 = -t33 * t380 + t367 * t38;
t293 = (-t206 * t282 - t270 * t389) * qJD(3) - t360;
t39 = -t236 * t367 + t271 * t380 + t293 - t353;
t296 = t33 * t381 - t367 * t39;
t295 = -t205 * t358 - t269 * t343;
t294 = t286 * (-t38 * t381 + t380 * t39);
t252 = qJD(3) * t267;
t251 = qJD(3) * t266;
t234 = (rSges(5,3) * t286 + t285 * t314) * t280;
t215 = -rSges(5,1) * t263 - rSges(5,2) * t264;
t202 = -rSges(4,1) * t261 - rSges(4,2) * t262;
t201 = -rSges(4,1) * t259 - rSges(4,2) * t260;
t196 = -Icges(4,5) * t261 - Icges(4,6) * t262;
t195 = -Icges(4,5) * t259 - Icges(4,6) * t260;
t185 = -Icges(4,1) * t248 - Icges(4,4) * t249;
t184 = -Icges(4,1) * t246 - Icges(4,4) * t247;
t183 = -Icges(4,4) * t248 - Icges(4,2) * t249;
t182 = -Icges(4,4) * t246 - Icges(4,2) * t247;
t181 = -Icges(4,5) * t248 - Icges(4,6) * t249;
t180 = -Icges(4,5) * t246 - Icges(4,6) * t247;
t179 = t204 * t389;
t160 = t190 * t389;
t154 = t262 * rSges(5,3) + t261 * t314;
t152 = t260 * rSges(5,3) + t259 * t314;
t137 = rSges(5,1) * t302 - rSges(5,2) * t224;
t133 = rSges(5,1) * t303 - rSges(5,2) * t222;
t89 = t318 * t359;
t88 = t319 * t359 + qJD(1);
t60 = t107 * t271 - t178 * t236 + t293;
t48 = t105 * t236 - t107 * t235 + t304;
t35 = -t178 * t355 - t236 * t118 + t271 * t86 + (-t392 + (t107 * t354 - t391) * t280) * qJD(3);
t28 = qJD(4) * t321 - t235 * t86 + t236 * t84 + t316;
t9 = qJD(5) * t155 + t404 * t271 - t378 * t236 - t367 * t355 + (-t392 + (t354 * t380 - t391) * t280) * qJD(3);
t1 = [m(4) * t89 + m(5) * t28 + m(6) * t7; m(4) * t320 + m(5) * (t281 * t34 - t283 * t35) + m(6) * (t281 * t8 - t283 * t9); (-t282 * (-t242 * t249 - t243 * t248 + t250 * t389 - t251 * t261 + t252 * t262) + (t281 * (-t163 * t249 - t165 * t248 + t180 * t389 - t182 * t261 + t184 * t262) + t283 * (-t164 * t249 - t166 * t248 + t181 * t389 - t183 * t261 + t185 * t262)) * t280) * t343 + (-t282 * (-t242 * t247 - t243 * t246 + t250 * t390 - t251 * t259 + t252 * t260) + (t281 * (-t163 * t247 - t165 * t246 + t180 * t390 - t182 * t259 + t184 * t260) + t283 * (-t164 * t247 - t166 * t246 + t181 * t390 - t183 * t259 + t185 * t260)) * t280) * t344 - (-t282 * (-t282 * t250 + (-t251 * t285 + t252 * t286 + (-t242 * t286 - t243 * t285) * qJD(3)) * t280) + (t281 * (-t282 * t180 + (-t182 * t285 + t184 * t286 + (-t163 * t286 - t165 * t285) * qJD(3)) * t280) + t283 * (-t282 * t181 + (-t183 * t285 + t185 * t286 + (-t164 * t286 - t166 * t285) * qJD(3)) * t280)) * t280) * t358 + (t282 ^ 2 * t250 + (((t281 * t371 + t283 * t370) * t286 + (t281 * t369 + t283 * t368) * t285) * t280 + (-t195 * t281 - t196 * t283 + t285 * t361 + t286 * t362) * t282) * t359) * t358 / 0.2e1 + ((t222 * t424 + t260 * t426 - t303 * t425) * t271 + t423 * t259 + (t222 * t427 + t260 * t431 - t303 * t429) * t236 + (t222 * t428 + t260 * t432 - t303 * t430) * t235 + (t260 * t441 + t262 * t440 - t386 * t435) * qJD(4)) * t419 + (-t482 * t282 + (t281 * t451 + t283 * t450) * t280) * t418 + ((t224 * t424 + t262 * t426 - t302 * t425) * t271 + t423 * t261 + (t224 * t427 + t262 * t431 - t302 * t429) * t236 + (t224 * t428 + t262 * t432 - t302 * t430) * t235 + (t260 * t439 + t262 * t438 + t386 * t434) * qJD(4)) * t417 + (-t481 * t282 + (t281 * t449 + t283 * t448) * t280) * t416 + ((t263 * t425 + t264 * t424) * t271 + (t263 * t429 + t264 * t427) * t236 + (t263 * t430 + t264 * t428) * t235 + (t260 * t437 + t262 * t436) * qJD(4) + ((-qJD(4) * t433 + t235 * t432 + t236 * t431 + t271 * t426) * t286 + t423 * t285) * t280) * t411 + (-t480 * t282 + (t281 * t446 + t283 * t445) * t280) * t410 - t455 * t282 / 0.2e1 + t390 * t461 + t389 * t462 + (t435 * t282 + (t281 * t441 + t283 * t440) * t280) * t337 + (-t434 * t282 + (t281 * t439 + t283 * t438) * t280) * t336 + (t433 * t282 + (t281 * t437 + t283 * t436) * t280) * t317 + (-(t260 * t297 + t262 * t296 + t280 * t294) * qJD(4) + t8 * t365 + t7 * t179 + (t347 * t9 + t454) * t282 + ((t9 * (-t270 - t367) + t301) * t283 + (t347 * t7 + t298) * t281) * t280 + (t259 * t350 - t295 - t374 * t271 + t363 * t236 + t349 * t282 + (-t256 - t378) * t389) * t39 + (t374 * t235 - t375 * t236 + t349 * t390 + t350 * t387 + t160 - t315) * t33 + (-t363 * t235 + t261 * t350 + t375 * t271 + t405 * t282 - t366 + t372) * t38) * m(6) + (-t59 * (-t271 * t152 + t235 * t234 + t366) - t60 * (t271 * t154 - t236 * t234 + t295) - t48 * (t152 * t236 - t154 * t235 + t315) - (t59 * (-t105 * t386 + t178 * t260) + t60 * (t107 * t386 - t178 * t262) + t48 * (t105 * t262 - t107 * t260)) * qJD(4) + t34 * t365 + t59 * t372 + t28 * t179 + t48 * t160 + (t35 * t379 + t403 * t60 + t329) * t282 + ((t35 * (-t178 - t270) + t60 * (-t118 - t256) + t28 * t105 + t48 * t84) * t283 + (t59 * t118 + t34 * t178 + t28 * t379 + t403 * t48) * t281) * t280) * m(5) + ((t109 * t167 - t110 * t168 + t186 * t90 - t187 * t91) * t282 + (t245 * t320 + t253 * t322 + t318 * t88 + t319 * t89) * t280 - ((t90 * t201 - t91 * t202) * t282 + (t88 * (t201 * t283 - t202 * t281) + t322 * t305) * t280) * qJD(3)) * m(4) - ((t281 * ((t196 * t390 + t259 * t368 + t260 * t370) * t389 + (t195 * t390 + t259 * t369 + t260 * t371) * t390 - (-t259 * t361 - t260 * t362 + t265 * t390) * t282) + t283 * ((t196 * t389 + t261 * t368 + t262 * t370) * t389 + (t195 * t389 + t261 * t369 + t262 * t371) * t390 - (-t261 * t361 - t262 * t362 + t265 * t389) * t282)) * qJD(3) ^ 2 + t447 * t354) * t280 / 0.2e1 - (t260 * t453 + t262 * t452) * qJD(4) / 0.2e1; (t222 * t421 + t259 * t420 - t303 * t422) * t419 + ((t285 * t482 - t356 * t435) * t280 + t450 * t261 + t451 * t259 + t463) * t418 + (t224 * t421 + t261 * t420 - t302 * t422) * t417 + ((t285 * t481 + t356 * t434) * t280 + t448 * t261 + t449 * t259 + t464) * t416 + t453 * t415 + t452 * t249 / 0.2e1 + t259 * t461 + t261 * t462 + (t263 * t422 + t264 * t421 + t387 * t420) * t411 + ((t285 * t480 - t356 * t433) * t280 + t445 * t261 + t446 * t259 + t465) * t410 + t455 * t387 / 0.2e1 + (t259 * t441 + t261 * t440 - t387 * t435) * t337 + (t259 * t439 + t261 * t438 + t387 * t434) * t336 + t447 * t331 + (t259 * t437 + t261 * t436 - t387 * t433) * t317 + (t296 * t249 + t297 * t247 + (-t367 * t9 - t378 * t39 + t301) * t261 + (-t33 * t404 - t380 * t7 + t298) * t259 + (qJD(3) * t294 + (-t38 * t405 + t380 * t9 + t39 * t404 - t454) * t285) * t280 - (t222 * t39 + t224 * t38 + t264 * t33) * qJD(5) - (t376 * t39 - t377 * t38) * t271 - (t33 * t377 - t364 * t39) * t236 - (-t33 * t376 + t364 * t38) * t235) * m(6) + (t28 * (t105 * t261 - t107 * t259) + (t259 * t59 - t261 * t60) * t118 + (t247 * t59 - t249 * t60 + t259 * t34 - t261 * t35) * t178 + ((-t105 * t59 + t107 * t60) * t356 + (t35 * t107 + t60 * t86 - t329) * t285) * t280 - t59 * (-t133 * t271 + t215 * t235) - t60 * (t137 * t271 - t215 * t236) + (-t133 * t236 + t137 * t235 - t259 * t86 + t261 * t84 + t321) * t48) * m(5); (-t303 * t9 - t302 * t8 + t263 * t7 + (t236 * t263 + t271 * t302 + t155) * t39 + (-t235 * t263 - t271 * t303 + t157) * t38 + (-t235 * t302 + t236 * t303 - t225) * t33) * m(6);];
tauc = t1(:);
