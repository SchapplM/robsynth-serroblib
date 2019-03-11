% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:52
% EndTime: 2019-03-09 01:30:05
% DurationCPUTime: 11.17s
% Computational Cost: add. (32453->555), mult. (39746->834), div. (0->0), fcn. (42815->8), ass. (0->328)
t341 = qJ(1) + pkin(9);
t340 = cos(t341);
t339 = sin(t341);
t478 = cos(qJ(1)) * pkin(1);
t368 = t478 + (-pkin(7) + qJ(3)) * t339;
t475 = pkin(2) + qJ(4);
t343 = sin(qJ(5));
t477 = pkin(5) * t343;
t395 = t475 + t477;
t345 = cos(qJ(6));
t342 = sin(qJ(6));
t442 = t342 * t343;
t278 = -t339 * t345 - t340 * t442;
t440 = t343 * t345;
t279 = -t339 * t342 + t340 * t440;
t346 = cos(qJ(5));
t445 = t340 * t346;
t218 = t279 * rSges(7,1) + t278 * rSges(7,2) - rSges(7,3) * t445;
t528 = pkin(8) * t445 - t218;
t161 = t340 * t395 + t368 - t528;
t276 = t339 * t442 - t340 * t345;
t277 = t339 * t440 + t340 * t342;
t236 = -rSges(7,1) * t276 - rSges(7,2) * t277;
t237 = rSges(7,1) * t278 - rSges(7,2) * t279;
t396 = -sin(qJ(1)) * pkin(1) + t340 * qJ(3);
t447 = t339 * t346;
t216 = t277 * rSges(7,1) - t276 * rSges(7,2) - rSges(7,3) * t447;
t541 = pkin(8) * t447 - t216;
t525 = -pkin(7) * t340 - t339 * t395 + t396 + t541;
t546 = m(7) * (t161 * t237 - t236 * t525);
t391 = rSges(7,1) * t345 - rSges(7,2) * t342;
t295 = rSges(7,3) * t343 + t346 * t391;
t452 = t218 * t343;
t186 = t295 * t445 + t452;
t453 = t216 * t343;
t544 = -t295 * t447 - t453;
t121 = t186 * t339 + t340 * t544;
t208 = -Icges(7,5) * t277 + Icges(7,6) * t276 + Icges(7,3) * t447;
t457 = t208 * t343;
t272 = Icges(7,4) * t277;
t210 = -Icges(7,2) * t276 - Icges(7,6) * t447 + t272;
t271 = Icges(7,4) * t276;
t214 = -Icges(7,1) * t277 + Icges(7,5) * t447 + t271;
t539 = t210 * t342 + t214 * t345;
t123 = t539 * t346 + t457;
t473 = m(7) * qJD(6);
t337 = t339 ^ 2;
t338 = t340 ^ 2;
t309 = t337 + t338;
t318 = rSges(6,1) * t346 - rSges(6,2) * t343;
t263 = t309 * t318;
t517 = m(7) / 0.2e1;
t326 = pkin(5) * t346 + pkin(8) * t343;
t414 = t295 + t326;
t239 = t414 * t339;
t241 = t414 * t340;
t527 = t239 * t339 + t241 * t340;
t538 = m(6) / 0.2e1;
t429 = t263 * t538 + t527 * t517;
t220 = ((-rSges(7,3) - pkin(8)) * t343 + (-pkin(5) - t391) * t346) * t339;
t448 = t339 * t343;
t292 = rSges(6,1) * t447 - rSges(6,2) * t448;
t293 = t318 * t340;
t373 = -t292 * t339 - t293 * t340;
t531 = t373 * t538;
t446 = t340 * t343;
t259 = rSges(7,3) * t446 + t391 * t445;
t403 = -pkin(5) * t445 - pkin(8) * t446 - t259;
t535 = t403 * t340;
t432 = (t220 * t339 + t535) * t517 + t531;
t59 = t432 - t429;
t545 = t59 * qJD(1);
t209 = Icges(7,5) * t279 + Icges(7,6) * t278 - Icges(7,3) * t445;
t379 = t210 * t276 + t214 * t277;
t460 = Icges(7,4) * t279;
t212 = Icges(7,2) * t278 - Icges(7,6) * t445 + t460;
t273 = Icges(7,4) * t278;
t215 = Icges(7,1) * t279 - Icges(7,5) * t445 + t273;
t425 = t278 * t212 + t279 * t215;
t543 = -t425 - (-t208 * t339 - t209 * t340) * t346 - t379;
t426 = t278 * t210 - t214 * t279;
t427 = -t276 * t212 + t277 * t215;
t542 = -t346 * (t208 * t340 - t209 * t339) - t427 - t426;
t518 = -m(7) / 0.2e1;
t505 = -t339 / 0.4e1;
t458 = Icges(7,4) * t345;
t387 = -Icges(7,2) * t342 + t458;
t283 = Icges(7,6) * t343 + t346 * t387;
t459 = Icges(7,4) * t342;
t389 = Icges(7,1) * t345 - t459;
t285 = Icges(7,5) * t343 + t346 * t389;
t384 = Icges(7,5) * t345 - Icges(7,6) * t342;
t281 = Icges(7,3) * t343 + t346 * t384;
t437 = t346 * t281;
t150 = t276 * t283 - t277 * t285 + t339 * t437;
t537 = t150 * t343;
t462 = Icges(6,4) * t343;
t388 = Icges(6,2) * t346 + t462;
t266 = Icges(6,6) * t340 + t339 * t388;
t461 = Icges(6,4) * t346;
t390 = Icges(6,1) * t343 + t461;
t268 = Icges(6,5) * t340 + t339 * t390;
t536 = (t266 * t346 + t268 * t343) * t339;
t101 = t339 * t161 + t340 * t525;
t506 = -t339 / 0.2e1;
t532 = t339 / 0.2e1;
t504 = -t340 / 0.2e1;
t503 = t340 / 0.2e1;
t431 = (t220 * t340 - t339 * t403) * t517 + (-t292 * t340 + t293 * t339) * t538;
t306 = (-Icges(7,2) * t345 - t459) * t346;
t307 = (-Icges(7,1) * t342 - t458) * t346;
t523 = -(t285 / 0.2e1 + t306 / 0.2e1) * t342 + (t307 / 0.2e1 - t283 / 0.2e1) * t345;
t522 = 0.2e1 * m(7);
t521 = 0.2e1 * t309;
t520 = 0.4e1 * qJD(1);
t519 = 2 * qJD(5);
t105 = t208 * t447 - t379;
t106 = -t209 * t447 + t427;
t382 = -t105 * t339 - t106 * t340;
t42 = t346 * t382 - t537;
t516 = -t42 / 0.2e1;
t229 = Icges(7,5) * t278 - Icges(7,6) * t279;
t421 = -Icges(7,2) * t279 + t215 + t273;
t423 = Icges(7,1) * t278 - t212 - t460;
t97 = t229 * t343 + (-t342 * t421 + t345 * t423) * t346;
t515 = -t97 / 0.2e1;
t258 = t295 * t339;
t294 = rSges(7,3) * t346 - t343 * t391;
t372 = t294 * t346 - t295 * t343;
t133 = -t216 * t346 - t258 * t343 - t339 * t372;
t134 = t218 * t346 + t259 * t343 + t340 * t372;
t513 = m(7) * (t133 * t525 + t134 * t161 - t186 * t403 + t220 * t544);
t308 = (-rSges(7,1) * t342 - rSges(7,2) * t345) * t346;
t511 = m(7) * (t101 * t308 - t236 * t241 + t237 * t239);
t508 = m(7) * (-t161 * t403 + t220 * t525);
t507 = t309 / 0.2e1;
t502 = t340 / 0.4e1;
t501 = t343 / 0.2e1;
t500 = t346 / 0.2e1;
t499 = m(4) * ((rSges(4,3) * t340 + t396) * t340 + (t478 + (rSges(4,3) + qJ(3)) * t339) * t339);
t408 = rSges(5,3) + t475;
t242 = rSges(5,2) * t340 - t339 * t408 + t396;
t243 = t478 + (rSges(5,2) + qJ(3)) * t339 + t408 * t340;
t498 = m(5) * (-t242 * t339 + t340 * t243);
t497 = m(5) * (t242 * t340 + t339 * t243);
t392 = rSges(6,1) * t343 + rSges(6,2) * t346;
t358 = t392 + t475;
t203 = (-rSges(6,3) - pkin(7)) * t340 - t358 * t339 + t396;
t334 = t339 * rSges(6,3);
t204 = t340 * t358 - t334 + t368;
t496 = m(6) * (-t203 * t292 + t204 * t293);
t495 = m(6) * (-t203 * t339 + t340 * t204);
t494 = m(6) * (t203 * t340 + t339 * t204);
t490 = m(7) * (t340 * t161 - t339 * t525);
t489 = m(7) * t101;
t488 = m(7) * (t186 * t340 - t339 * t544);
t487 = m(7) * t121;
t376 = -t236 * t339 - t237 * t340;
t484 = m(7) * t376;
t168 = -t236 * t340 + t237 * t339;
t483 = m(7) * t168;
t480 = m(7) * (t239 * t340 - t241 * t339);
t476 = qJD(5) / 0.2e1;
t470 = t339 * t42;
t151 = t278 * t283 + t279 * t285 - t340 * t437;
t148 = t151 * t343;
t107 = t208 * t445 + t426;
t108 = -t209 * t445 + t425;
t381 = -t107 * t339 - t108 * t340;
t43 = t346 * t381 + t148;
t469 = t340 * t43;
t466 = -t105 + t543;
t465 = t106 + t542;
t464 = -t107 - t542;
t463 = t108 + t543;
t454 = t209 * t343;
t385 = Icges(6,5) * t343 + Icges(6,6) * t346;
t264 = Icges(6,3) * t340 + t339 * t385;
t450 = t264 * t339;
t449 = t281 * t343;
t282 = Icges(7,6) * t346 - t343 * t387;
t444 = t342 * t282;
t443 = t342 * t283;
t305 = (-Icges(7,5) * t342 - Icges(7,6) * t345) * t346;
t441 = t343 * t305;
t284 = Icges(7,5) * t346 - t343 * t389;
t439 = t345 * t284;
t438 = t345 * t285;
t88 = t133 * t339 - t134 * t340;
t436 = t88 * qJD(3);
t435 = t88 * qJD(6);
t424 = -Icges(7,1) * t276 - t210 - t272;
t422 = -Icges(7,2) * t277 - t214 - t271;
t419 = t266 * t445 + t268 * t446;
t417 = -t283 + t307;
t416 = t285 + t306;
t415 = pkin(8) * t346 + t294 - t477;
t413 = qJD(1) * t343;
t412 = qJD(1) * t346;
t411 = qJD(6) * t346;
t170 = (t518 - m(6) / 0.2e1 - m(5) / 0.2e1) * t521;
t409 = t170 * qJD(1);
t12 = t537 + (t339 * t463 + t340 * t464) * t346;
t406 = -t12 / 0.2e1 + t516;
t13 = t148 + (t339 * t465 + t340 * t466) * t346;
t405 = t43 / 0.2e1 - t13 / 0.2e1;
t267 = -Icges(6,6) * t339 + t340 * t388;
t320 = Icges(6,4) * t445;
t269 = Icges(6,1) * t446 - Icges(6,5) * t339 + t320;
t163 = t340 * (-Icges(6,3) * t339 + t340 * t385) + t267 * t447 + t269 * t448;
t401 = -t447 / 0.4e1;
t399 = -t445 / 0.4e1;
t125 = -t276 * t416 + t277 * t417 - t305 * t447;
t126 = t278 * t416 + t279 * t417 - t305 * t445;
t228 = -Icges(7,5) * t276 - Icges(7,6) * t277;
t96 = t228 * t343 + (-t342 * t422 + t345 * t424) * t346;
t393 = t511 / 0.2e1 + (t126 + t97) * t505 + (t125 + t96) * t502;
t315 = Icges(6,1) * t346 - t462;
t313 = -Icges(6,2) * t343 + t461;
t386 = Icges(6,5) * t346 - Icges(6,6) * t343;
t360 = m(7) * (t133 * t340 + t134 * t339);
t371 = m(7) * t308 * t507;
t74 = -t360 / 0.2e1 + t371;
t109 = (-t258 * t346 + t453) * t340 + (t259 * t346 - t452) * t339;
t83 = (t109 / 0.4e1 - t376 / 0.4e1) * t522;
t383 = t83 * qJD(2) - t74 * qJD(4);
t377 = -t212 * t342 + t215 * t345;
t124 = t346 * t377 + t454;
t380 = t123 * t339 - t124 * t340;
t374 = t438 - t443;
t162 = t264 * t340 + t536;
t164 = t419 - t450;
t17 = t339 * t466 - t340 * t465;
t51 = t107 * t340 - t108 * t339;
t370 = t164 * t504 + t419 * t503 + (-t162 + t536) * t532 - t51 / 0.2e1 + t17 / 0.2e1;
t16 = t339 * t464 - t340 * t463;
t50 = t105 * t340 - t106 * t339;
t369 = (t163 - t164 - t450) * t506 + t264 * t338 / 0.2e1 + t162 * t504 + t163 * t532 - t16 / 0.2e1 - t50 / 0.2e1;
t79 = -t228 * t447 - t276 * t422 + t277 * t424;
t80 = -t229 * t447 - t276 * t421 + t277 * t423;
t29 = -t339 * t80 + t340 * t79;
t81 = -t228 * t445 + t278 * t422 + t279 * t424;
t82 = -t229 * t445 + t278 * t421 + t279 * t423;
t30 = -t339 * t82 + t340 * t81;
t367 = t29 * t503 + t30 * t506;
t364 = t281 * t339 + t539;
t363 = t281 * t340 - t377;
t280 = Icges(7,3) * t346 - t343 * t384;
t362 = t280 - t374;
t288 = t313 * t339;
t290 = t315 * t339;
t357 = (t268 + t288) * t346 + (-t266 + t290) * t343;
t289 = -Icges(6,2) * t446 + t320;
t291 = t315 * t340;
t356 = (-t269 - t289) * t346 + (t267 - t291) * t343;
t355 = t12 * t505 + t13 * t502 + t17 * t401 - t470 / 0.4e1 - t469 / 0.4e1 + t51 * t447 / 0.4e1 + (t16 + t50) * t399;
t351 = -t346 * t362 + t449;
t112 = -t276 * t282 + t277 * t284 + t339 * t351;
t113 = t278 * t282 + t279 * t284 + t340 * t351;
t132 = (t281 + t439 - t444) * t346 + t362 * t343;
t189 = t346 * t374 + t449;
t254 = t283 * t339;
t256 = t285 * t339;
t91 = (-t254 * t342 + t256 * t345 - t208) * t346 + t364 * t343;
t255 = t283 * t340;
t257 = t285 * t340;
t92 = (-t255 * t342 + t257 * t345 + t209) * t346 + t363 * t343;
t354 = t132 * t501 + t189 * t500 + t513 / 0.2e1 + (-t123 - t150) * t448 / 0.4e1 + (t112 + t91) * t401 + (t124 + t151) * t446 / 0.4e1 + (t113 + t92) * t399;
t353 = -t346 * t364 - t457;
t352 = -t346 * t363 + t454;
t350 = t439 / 0.2e1 - t444 / 0.2e1 + t281 / 0.2e1 - t390 / 0.2e1 - t313 / 0.2e1;
t349 = t438 / 0.2e1 - t443 / 0.2e1 + t315 / 0.2e1 - t388 / 0.2e1 - t280 / 0.2e1;
t287 = t340 * t386;
t286 = t386 * t339;
t240 = t415 * t340;
t238 = t415 * t339;
t197 = t237 * t343 + t308 * t445;
t196 = -t236 * t343 - t308 * t447;
t173 = t480 / 0.2e1;
t169 = (-m(7) / 0.4e1 - m(6) / 0.4e1 - m(5) / 0.4e1) * t521 + (m(5) + m(6) + m(7)) * t507;
t158 = t483 / 0.2e1;
t157 = t484 / 0.2e1;
t155 = t168 * t346;
t147 = (t441 + (-t342 * t416 + t345 * t417) * t346) * t343;
t140 = (-t216 * t340 + t218 * t339) * t346;
t136 = t535 + (-t326 * t339 - t258) * t339;
t128 = (-pkin(5) * t446 + t528) * t340 + (-pkin(5) * t448 + t541) * t339;
t117 = t487 / 0.2e1;
t116 = t488 / 0.2e1;
t85 = m(7) * t88 * t476;
t84 = (t109 + t376) * t517;
t78 = t255 * t278 + t257 * t279 + t340 * t352;
t77 = t254 * t278 + t256 * t279 + t340 * t353;
t76 = -t255 * t276 + t257 * t277 + t339 * t352;
t75 = -t254 * t276 + t256 * t277 + t339 * t353;
t73 = t360 / 0.2e1 + t371;
t67 = t173 - t431;
t66 = t173 + t431;
t65 = -t480 / 0.2e1 + t431;
t62 = t128 * t376 + t308 * t527;
t60 = t429 + t432;
t58 = t546 + t441 / 0.2e1 + t523 * t346;
t47 = t490 + t495 + t498;
t46 = t189 * t343 + t346 * t380;
t44 = t489 + t494 + t497 + t499;
t36 = t117 - t484 / 0.2e1;
t35 = t116 - t483 / 0.2e1;
t34 = t158 + t116;
t33 = t158 - t488 / 0.2e1;
t32 = t157 + t117;
t31 = t157 - t487 / 0.2e1;
t28 = -t343 * t349 + t346 * t350 + t496 + t508;
t27 = -t339 * t78 + t340 * t77;
t26 = -t339 * t76 + t340 * t75;
t23 = t109 * t140 + t133 * t544 + t134 * t186;
t22 = t126 * t343 + (-t339 * t81 - t340 * t82) * t346;
t21 = t125 * t343 + (-t339 * t79 - t340 * t80) * t346;
t18 = (-t91 * t339 - t92 * t340 + t189) * t346 + (t132 - t380) * t343;
t9 = (-t339 * t77 - t340 * t78 + t151) * t346 + (t113 - t381) * t343;
t8 = (-t339 * t75 - t340 * t76 - t150) * t346 + (t112 - t382) * t343;
t7 = m(7) * t62 + t367;
t6 = (t339 * t405 + t340 * t406) * t346;
t5 = t339 * t369 + t340 * t370;
t4 = m(7) * t23 + (t9 * t504 + t8 * t506 + t46 / 0.2e1) * t346 + (t469 / 0.2e1 + t470 / 0.2e1 + t18 / 0.2e1) * t343;
t3 = (-t189 / 0.2e1 + (t113 / 0.4e1 + t92 / 0.4e1) * t340 + (t112 / 0.4e1 + t91 / 0.4e1) * t339) * t346 + t355 - t513 / 0.2e1 + (-t132 / 0.2e1 + (-t151 / 0.4e1 - t124 / 0.4e1) * t340 + (t150 / 0.4e1 + t123 / 0.4e1) * t339) * t343 + t393;
t2 = t355 - t511 / 0.2e1 + t354 + (-t125 / 0.4e1 - t96 / 0.4e1) * t340 + (t126 / 0.4e1 + t97 / 0.4e1) * t339;
t1 = t354 + (t12 / 0.4e1 + t42 / 0.4e1 + (-t51 / 0.4e1 + t17 / 0.4e1) * t346) * t339 + (-t13 / 0.4e1 + t43 / 0.4e1 + (t50 / 0.4e1 + t16 / 0.4e1) * t346) * t340 + t393;
t10 = [t44 * qJD(3) + t47 * qJD(4) + t28 * qJD(5) + t58 * qJD(6), 0, qJD(1) * t44 + qJD(4) * t169 + qJD(5) * t60 + qJD(6) * t32, qJD(1) * t47 + qJD(3) * t169 + qJD(5) * t66 + qJD(6) * t34, t28 * qJD(1) + t60 * qJD(3) + t66 * qJD(4) + t1 * qJD(6) + (m(7) * (t161 * t238 + t220 * t241 - t239 * t403 + t240 * t525) + (t91 / 0.2e1 + t112 / 0.2e1 + m(6) * (-t203 * t392 - t292 * t318) - t385 * t503 + (-t266 / 0.2e1 + t290 / 0.2e1) * t346 + (-t268 / 0.2e1 - t288 / 0.2e1) * t343 - t370) * t340 + (-t92 / 0.2e1 - t113 / 0.2e1 + m(6) * (-t204 * t392 + t293 * t318) + t385 * t506 + (t267 / 0.2e1 - t291 / 0.2e1) * t346 + (t269 / 0.2e1 + t289 / 0.2e1) * t343 - t369) * t339) * qJD(5), t58 * qJD(1) + t32 * qJD(3) + t34 * qJD(4) + t1 * qJD(5) + t147 * qJD(6) + (t161 * t197 + t186 * t237 + t196 * t525 - t236 * t544) * t473 + ((t515 - t126 / 0.2e1 - t406) * t340 + (-t96 / 0.2e1 - t125 / 0.2e1 - t405) * t339) * t411; 0, 0, 0, 0 (t136 * t517 + t531) * t519 + t84 * qJD(6), t84 * qJD(5) + t155 * t473; t170 * qJD(4) + t59 * qJD(5) + t31 * qJD(6) + (-t489 / 0.4e1 - t494 / 0.4e1 - t499 / 0.4e1 - t497 / 0.4e1) * t520, 0, 0, t409, t545 + ((-t238 * t340 + t240 * t339) * t476 + t435 / 0.4e1) * t522, t31 * qJD(1) + t85 + (t196 * t339 - t197 * t340) * t473; -t170 * qJD(3) + t65 * qJD(5) + t33 * qJD(6) + (-t490 / 0.4e1 - t495 / 0.4e1 - t498 / 0.4e1) * t520, 0, -t409, 0, t65 * qJD(1) + ((t238 * t339 + t240 * t340) * t517 - m(6) * t392 * t507) * t519 + t73 * qJD(6), t33 * qJD(1) + t73 * qJD(5) + (t196 * t340 + t197 * t339) * t473; -t59 * qJD(3) + t67 * qJD(4) + t5 * qJD(5) + t3 * qJD(6) + (-t508 / 0.4e1 - t496 / 0.4e1) * t520 - t350 * t412 + t349 * t413, -qJD(6) * t83, t435 * t518 - t545, qJD(1) * t67 + qJD(6) * t74, t5 * qJD(1) + (m(6) * ((-t339 * (rSges(6,3) * t340 + t339 * t392) + (-t340 * t392 + t334) * t340) * t373 - t263 * t392) + m(7) * (t128 * t136 + t238 * t239 + t240 * t241) + (t337 * t287 + (t357 * t340 + (-t286 + t356) * t339) * t340 + t27) * t506 + (t338 * t286 + (t356 * t339 + (-t287 + t357) * t340) * t339 + t26) * t503) * qJD(5) + t7 * qJD(6), t3 * qJD(1) + t7 * qJD(5) + t436 * t518 + (-t46 / 0.2e1 + (-t30 / 0.2e1 + t9 / 0.2e1) * t340 + (-t29 / 0.2e1 + t8 / 0.2e1) * t339) * t411 - t383 + (t21 * t503 + t22 * t506 + (-t18 / 0.2e1 + (t96 / 0.2e1 - t43 / 0.2e1) * t340 + (t515 + t516) * t339) * t343 + (t121 * t308 + t128 * t155 + t140 * t376 + t196 * t241 + t197 * t239 - t23) * m(7)) * qJD(6); -t305 * t413 / 0.2e1 + t36 * qJD(3) + t35 * qJD(4) + t2 * qJD(5) + t6 * qJD(6) - qJD(1) * t546 - t523 * t412, qJD(5) * t83, qJD(1) * t36 + t85, qJD(1) * t35 - qJD(5) * t74, t2 * qJD(1) + (t9 * t506 + t51 * t446 / 0.2e1 - t27 * t445 / 0.2e1 + t8 * t503 + t50 * t448 / 0.2e1 - t26 * t447 / 0.2e1 + (-t123 * t340 - t124 * t339) * t500 + (-t92 * t339 + t91 * t340) * t501 - t367) * qJD(5) + t4 * qJD(6) + (t436 / 0.2e1 + (t109 * t128 + t133 * t241 + t134 * t239 + t136 * t140 + t186 * t238 + t240 * t544 - t62) * qJD(5)) * m(7) + t383, t6 * qJD(1) + t4 * qJD(5) + (t147 * t501 + m(7) * (t140 * t155 + t186 * t197 + t196 * t544) + (t22 * t504 + t21 * t506 + (-t339 * t96 - t340 * t97) * t501) * t346) * qJD(6);];
Cq  = t10;
