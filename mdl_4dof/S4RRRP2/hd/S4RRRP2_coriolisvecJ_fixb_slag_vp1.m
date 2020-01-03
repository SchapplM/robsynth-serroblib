% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:13:06
% DurationCPUTime: 10.16s
% Computational Cost: add. (7941->394), mult. (8585->484), div. (0->0), fcn. (6686->6), ass. (0->244)
t241 = sin(qJ(3));
t243 = cos(qJ(3));
t184 = Icges(5,5) * t243 - Icges(5,6) * t241;
t186 = Icges(4,5) * t243 - Icges(4,6) * t241;
t481 = t184 + t186;
t498 = Icges(4,5) + Icges(5,5);
t497 = Icges(4,6) + Icges(5,6);
t496 = Icges(4,3) + Icges(5,3);
t239 = qJ(1) + qJ(2);
t234 = cos(t239);
t233 = sin(t239);
t379 = t233 * t243;
t380 = t233 * t241;
t120 = Icges(5,4) * t379 - Icges(5,2) * t380 - Icges(5,6) * t234;
t122 = Icges(4,4) * t379 - Icges(4,2) * t380 - Icges(4,6) * t234;
t486 = t120 + t122;
t209 = Icges(5,4) * t380;
t124 = Icges(5,1) * t379 - Icges(5,5) * t234 - t209;
t210 = Icges(4,4) * t380;
t126 = Icges(4,1) * t379 - Icges(4,5) * t234 - t210;
t493 = t124 + t126;
t396 = Icges(5,4) * t241;
t192 = Icges(5,1) * t243 - t396;
t286 = t192 * t234;
t125 = Icges(5,5) * t233 + t286;
t397 = Icges(4,4) * t241;
t194 = Icges(4,1) * t243 - t397;
t287 = t194 * t234;
t127 = Icges(4,5) * t233 + t287;
t484 = t125 + t127;
t187 = Icges(5,2) * t243 + t396;
t189 = Icges(4,2) * t243 + t397;
t495 = t187 + t189;
t235 = Icges(5,4) * t243;
t191 = Icges(5,1) * t241 + t235;
t236 = Icges(4,4) * t243;
t193 = Icges(4,1) * t241 + t236;
t491 = -t191 - t193;
t494 = t481 * t234;
t465 = t496 * t234 - t498 * t379 + t497 * t380;
t464 = t496 * t233 + t494;
t490 = t192 + t194;
t303 = -Icges(5,2) * t241 + t235;
t304 = -Icges(4,2) * t241 + t236;
t489 = t303 + t304;
t488 = t486 * t241;
t487 = t484 * t379;
t284 = t303 * t234;
t121 = Icges(5,6) * t233 + t284;
t285 = t304 * t234;
t123 = Icges(4,6) * t233 + t285;
t485 = t121 + t123;
t477 = t495 * t241 + t491 * t243;
t460 = -t493 * t243 + t488;
t483 = t489 * qJD(3);
t482 = t490 * qJD(3);
t183 = Icges(5,5) * t241 + Icges(5,6) * t243;
t185 = Icges(4,5) * t241 + Icges(4,6) * t243;
t480 = t185 + t183;
t238 = qJD(1) + qJD(2);
t479 = -qJD(3) * t495 + t497 * t238;
t478 = t491 * qJD(3) + t498 * t238;
t476 = -t464 * t234 + t487;
t374 = t234 * t243;
t475 = t465 * t233 - t493 * t374;
t436 = t464 * t233 + t484 * t374;
t474 = -t460 * t233 + t465 * t234;
t447 = -t485 * t380 + t476;
t375 = t234 * t241;
t446 = -t486 * t375 - t475;
t445 = -t485 * t375 + t436;
t384 = t185 * t234;
t387 = t183 * t234;
t473 = -t477 * t233 - t384 - t387;
t385 = t185 * t233;
t388 = t183 * t233;
t472 = -t477 * t234 + t385 + t388;
t230 = t234 * pkin(6);
t159 = pkin(2) * t233 - t230;
t240 = -qJ(4) - pkin(6);
t215 = t234 * t240;
t412 = pkin(3) * t243;
t232 = pkin(2) + t412;
t279 = -rSges(5,1) * t379 + rSges(5,2) * t380 + t234 * rSges(5,3) - t233 * t232 - t215;
t457 = t159 + t279;
t471 = t238 * t457;
t470 = t485 * t241;
t444 = t493 * t241 + t486 * t243;
t443 = t484 * t241 + t485 * t243;
t382 = t233 * t238;
t469 = t479 * t234 - t382 * t489;
t468 = (t284 + t285) * t238 + t479 * t233;
t467 = t478 * t234 - t490 * t382;
t466 = (-t286 - t287) * t238 - t478 * t233;
t463 = t482 * t243 - t483 * t241 + t480 * t238 + (t491 * t241 - t243 * t495) * qJD(3);
t462 = -t480 * qJD(3) + t496 * t238;
t461 = -t484 * t243 + t470;
t459 = t481 * qJD(3) + t477 * t238;
t458 = t472 * t238;
t456 = (t445 * t233 - t446 * t234) * qJD(3);
t455 = (t447 * t233 - t474 * t234) * qJD(3);
t454 = t473 * t238;
t453 = t454 + t455;
t452 = t456 + t458;
t451 = qJD(3) * t460 + t241 * t466 - t243 * t468;
t450 = -qJD(3) * t461 + t241 * t467 + t243 * t469;
t449 = t233 * t459 + t234 * t463;
t448 = t233 * t463 - t234 * t459;
t442 = -qJD(3) * t443 + t238 * t464 - t241 * t469 + t243 * t467;
t441 = qJD(3) * t444 + t465 * t238 + t468 * t241 + t466 * t243;
t351 = rSges(4,2) * t380 + t234 * rSges(4,3);
t129 = rSges(4,1) * t379 - t351;
t115 = t238 * t129;
t154 = t238 * t159;
t376 = t234 * t238;
t198 = pkin(6) * t376;
t346 = qJD(3) * t243;
t329 = t234 * t346;
t373 = t238 * t241;
t341 = t233 * t373;
t288 = rSges(4,3) * t376 + (-t329 + t341) * rSges(4,2);
t347 = qJD(3) * t241;
t330 = t234 * t347;
t440 = -rSges(4,1) * t330 + t115 + t154 + t198 + t288;
t340 = t238 * t379;
t216 = qJD(4) * t233;
t430 = rSges(5,2) * t341 + rSges(5,3) * t376 + t216;
t439 = -rSges(5,1) * t340 - t232 * t382 + t154 + t430 - t471;
t438 = t465 + t470;
t437 = t233 * t346 + t234 * t373;
t242 = sin(qJ(1));
t405 = pkin(1) * qJD(1);
t342 = t242 * t405;
t157 = rSges(3,1) * t233 + rSges(3,2) * t234;
t389 = t157 * t238;
t133 = -t342 - t389;
t435 = (t460 + t494) * t238 + t462 * t233;
t434 = t462 * t234 + t461 * t238 - t481 * t382;
t433 = 0.2e1 * qJD(3);
t350 = rSges(5,1) * t374 + t233 * rSges(5,3);
t431 = rSges(4,1) * t374 + t233 * rSges(4,3);
t403 = t243 * rSges(5,2);
t199 = rSges(5,1) * t241 + t403;
t320 = pkin(3) * t241 + t199;
t348 = qJD(3) * t234;
t277 = -t320 * t348 + t216;
t229 = t233 * pkin(6);
t160 = t234 * pkin(2) + t229;
t381 = t233 * t240;
t296 = t234 * t232 + t350 - t381;
t332 = t233 * t347;
t217 = qJD(4) * t234;
t356 = pkin(3) * t332 + t217;
t429 = -rSges(5,1) * t332 - rSges(5,2) * t437 - t238 * t381 - t356;
t200 = rSges(4,1) * t241 + rSges(4,2) * t243;
t349 = qJD(3) * t233;
t156 = t200 * t349;
t131 = -rSges(4,2) * t375 + t431;
t278 = t131 + t160;
t423 = -t238 * t278 + t156;
t367 = -rSges(5,2) * t375 - t160 + t296;
t420 = t199 * t349 - t238 * (t160 + t367) + t356;
t360 = -Icges(4,2) * t379 + t126 - t210;
t364 = t193 * t233 + t122;
t419 = -t241 * t360 - t243 * t364;
t362 = -Icges(5,2) * t379 + t124 - t209;
t366 = t191 * t233 + t120;
t418 = -t241 * t362 - t243 * t366;
t417 = t233 / 0.2e1;
t416 = -t234 / 0.2e1;
t414 = t238 / 0.2e1;
t413 = pkin(1) * t242;
t244 = cos(qJ(1));
t237 = t244 * pkin(1);
t411 = pkin(2) - t232;
t280 = -t330 - t340;
t410 = -pkin(3) * t330 - t198 + (t233 * t411 - t215) * t238 + rSges(5,1) * t280 - rSges(5,2) * t329 + t430;
t409 = t429 + (-t234 * t411 - t229 + t350) * t238;
t408 = rSges(4,1) * t243;
t407 = rSges(5,1) * t243;
t406 = rSges(5,2) * t241;
t333 = t200 * t348;
t281 = -t333 - t342;
t62 = (-t129 - t159) * t238 + t281;
t404 = t234 * t62;
t386 = t184 * t238;
t383 = t186 * t238;
t365 = -t191 * t234 - t121;
t363 = -t193 * t234 - t123;
t361 = -t187 * t234 + t125;
t359 = -t189 * t234 + t127;
t355 = -t187 + t192;
t354 = t191 + t303;
t353 = -t189 + t194;
t352 = t193 + t304;
t246 = qJD(1) ^ 2;
t345 = t246 * t413;
t344 = t246 * t237;
t343 = t244 * t405;
t335 = rSges(4,1) * t332 + rSges(4,2) * t437;
t326 = -pkin(2) - t408;
t325 = -t349 / 0.2e1;
t322 = t348 / 0.2e1;
t202 = -t406 + t407;
t319 = -t202 - t412;
t158 = t234 * rSges(3,1) - rSges(3,2) * t233;
t311 = t238 * (-pkin(2) * t382 + t198) - t345;
t310 = -pkin(3) * t375 - t199 * t234;
t136 = rSges(3,1) * t376 - rSges(3,2) * t382;
t308 = -rSges(4,2) * t241 + t408;
t63 = t343 - t423;
t307 = -t233 * t63 - t404;
t171 = t202 * qJD(3);
t294 = (-t412 * qJD(3) - t171) * qJD(3);
t66 = (t129 * t233 + t131 * t234) * qJD(3);
t276 = -t233 * t457 + t234 * t367;
t275 = -t241 * t361 + t243 * t365;
t274 = -t241 * t359 + t243 * t363;
t273 = t233 * t326 + t230 + t351;
t272 = -rSges(5,3) * t382 - t429;
t271 = (-t241 * t354 + t243 * t355) * t238;
t270 = (-t241 * t352 + t243 * t353) * t238;
t263 = t277 - t342;
t249 = (t326 * t404 + (t62 * (-rSges(4,3) - pkin(6)) + t63 * t326) * t233) * t238;
t248 = ((t436 * t233 + ((t464 + t488) * t234 + t447 + t475 - t487) * t234) * qJD(3) + t458) * t322 + (-t477 * qJD(3) + t482 * t241 + t483 * t243) * t238 + (((t234 * t438 - t436 + t445) * t234 + (t233 * t438 + t446 - t476) * t233) * qJD(3) + t453 - t454) * t325 + (t449 + t450) * t349 / 0.2e1 - (t448 - t451 + t452) * t348 / 0.2e1 + ((t444 + t473) * t233 + (t443 + t472) * t234) * qJD(3) * t414;
t22 = t294 * t233 + (t410 + t277) * t238 + t311;
t48 = (-t159 + t457) * t238 + t263;
t49 = t343 - t420;
t247 = (-t22 * t406 + (t48 * (-t232 - t407) - t49 * t240) * t238 + t49 * (-t403 + (-rSges(5,1) - pkin(3)) * t241) * qJD(3)) * t234;
t172 = t308 * qJD(3);
t153 = t200 * t234;
t151 = t200 * t233;
t150 = t199 * t233;
t137 = t160 * t238;
t134 = t158 * t238 + t343;
t109 = -t136 * t238 - t344;
t108 = -t238 * t389 - t345;
t92 = t238 * t431 - t335;
t90 = rSges(4,1) * t280 + t288;
t46 = -t344 - t172 * t348 + (-t137 - t92 + t156) * t238;
t45 = t238 * t90 + (-t172 * t233 - t200 * t376) * qJD(3) + t311;
t40 = t276 * qJD(3);
t23 = -t344 + t294 * t234 + (t320 * t349 - t137 + t217 - t409) * t238;
t1 = [m(3) * (t109 * (-t157 - t413) + t108 * (t158 + t237) + (-t136 - t343 + t134) * t133) + t248 + (t23 * (t279 - t413) + t48 * (t272 - t343) + t22 * (t237 + t296) + t247 + (t48 - t263 - t342 + t439) * t49) * m(5) + (t46 * (t273 - t413) + t62 * (t335 - t343) + t45 * (t237 + t278) + t249 + (-t342 - t281 + t62 + t440) * t63) * m(4); t248 + (t22 * t296 + t23 * t279 + t247 + (-t277 + t439) * t49 + (-t420 + t272) * t48) * m(5) + (t273 * t46 + t278 * t45 + t249 + (t333 + t440) * t63 + (-t423 + t335) * t62) * m(4) + (t108 * t158 - t109 * t157 - t133 * t136 - t134 * t389 - (-t133 * t158 - t134 * t157) * t238) * m(3); -(((t352 + t354) * t243 + (t353 + t355) * t241) * t238 + (((-t360 - t362) * t234 + (t359 + t361) * t233) * t243 + ((t364 + t366) * t234 + (t363 + t365) * t233) * t241) * qJD(3)) * t238 / 0.2e1 + ((t443 * t238 + t451) * t234 + (t238 * t444 + t450) * t233) * t414 + ((-t349 * t387 + t386) * t233 + (t271 + (-t418 * t234 + (t388 + t275) * t233) * qJD(3)) * t234 + (-t349 * t384 + t383) * t233 + (t270 + (-t419 * t234 + (t385 + t274) * t233) * qJD(3)) * t234) * t325 + ((-t348 * t388 - t386) * t234 + (t271 + (t275 * t233 + (t387 - t418) * t234) * qJD(3)) * t233 + (-t348 * t385 - t383) * t234 + (t270 + (t274 * t233 + (t384 - t419) * t234) * qJD(3)) * t233) * t322 + (0.2e1 * t66 * ((t90 + t115) * t234 + (-t131 * t238 + t92) * t233) + t307 * t172 + ((-t238 * t63 - t46) * t234 + (t238 * t62 - t45) * t233) * t200 - (t151 * t62 - t153 * t63) * t238 - (t66 * (-t151 * t233 - t153 * t234) + t307 * t308) * qJD(3)) * m(4) + (t449 * t238 + ((t234 * t441 + t238 * t445) * t234 + (t434 * t233 + t446 * t238 + (-t435 + t442) * t234) * t233) * t433) * t417 + (t448 * t238 + ((t234 * t435 + t238 * t447) * t234 + (t442 * t233 + t474 * t238 + (-t434 + t441) * t234) * t233) * t433) * t416 + ((-t22 * t320 - t49 * t171 + t40 * t409 + (t48 * t199 - t367 * t40) * t238) * t233 + (-t23 * t320 - t48 * t171 + t40 * t410 + (-t320 * t49 - t40 * t457) * t238) * t234 - (t48 * t150 + t310 * t49) * t238 + (((t410 - t471) * t234 + (-t238 * t367 + t409) * t233) * t276 + (-t233 * t49 - t234 * t48) * t412 - (t310 * t40 + t319 * t48) * t234 - (t49 * t319 + (-pkin(3) * t380 - t150) * t40) * t233) * qJD(3)) * m(5) + (t453 + t455) * t382 / 0.2e1 + (t452 + t456) * t376 / 0.2e1; 0.2e1 * (t22 * t416 + t23 * t417) * m(5);];
tauc = t1(:);
