% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:40
% EndTime: 2019-03-09 01:31:54
% DurationCPUTime: 11.44s
% Computational Cost: add. (46170->561), mult. (39908->833), div. (0->0), fcn. (42977->10), ass. (0->332)
t348 = qJ(1) + pkin(9);
t346 = cos(t348);
t482 = -pkin(7) - qJ(4);
t333 = t346 * t482;
t344 = sin(t348);
t486 = cos(qJ(1)) * pkin(1);
t396 = t346 * pkin(2) + t486;
t349 = sin(pkin(10));
t485 = pkin(4) * t349;
t399 = qJ(3) + t485;
t347 = pkin(10) + qJ(5);
t343 = sin(t347);
t484 = pkin(5) * t343;
t353 = cos(qJ(6));
t446 = t346 * t353;
t351 = sin(qJ(6));
t450 = t344 * t351;
t295 = -t343 * t450 + t446;
t447 = t346 * t351;
t449 = t344 * t353;
t296 = t343 * t449 + t447;
t345 = cos(t347);
t451 = t344 * t345;
t220 = t296 * rSges(7,1) + t295 * rSges(7,2) - rSges(7,3) * t451;
t542 = -pkin(8) * t451 + t220;
t166 = -t333 + (t399 + t484) * t344 + t396 + t542;
t243 = rSges(7,1) * t295 - rSges(7,2) * t296;
t297 = t343 * t447 + t449;
t298 = t343 * t446 - t450;
t244 = rSges(7,1) * t297 + rSges(7,2) * t298;
t448 = t345 * t346;
t328 = pkin(8) * t448;
t400 = -sin(qJ(1)) * pkin(1) + t346 * qJ(3);
t361 = t346 * t485 + t400 + (-pkin(2) + t482) * t344;
t478 = rSges(7,3) * t345;
t540 = t298 * rSges(7,1) - t297 * rSges(7,2);
t552 = (-t478 + t484) * t346 - t328 + t361 + t540;
t565 = m(7) * (t166 * t243 - t244 * t552);
t480 = m(7) * qJD(6);
t211 = Icges(7,5) * t296 + Icges(7,6) * t295 - Icges(7,3) * t451;
t467 = Icges(7,4) * t296;
t214 = Icges(7,2) * t295 - Icges(7,6) * t451 + t467;
t287 = Icges(7,4) * t295;
t217 = Icges(7,1) * t296 - Icges(7,5) * t451 + t287;
t112 = t211 * t448 + t297 * t214 - t298 * t217;
t564 = t112 * t344;
t563 = t112 * t346;
t213 = -Icges(7,5) * t298 + Icges(7,6) * t297 + Icges(7,3) * t448;
t457 = t213 * t343;
t289 = Icges(7,4) * t298;
t216 = Icges(7,2) * t297 + Icges(7,6) * t448 - t289;
t288 = Icges(7,4) * t297;
t218 = Icges(7,1) * t298 - Icges(7,5) * t448 - t288;
t558 = t216 * t351 + t218 * t353;
t126 = t345 * t558 - t457;
t318 = rSges(6,1) * t345 - rSges(6,2) * t343;
t341 = t344 ^ 2;
t342 = t346 ^ 2;
t321 = t341 + t342;
t526 = -m(7) / 0.2e1;
t392 = rSges(7,1) * t353 - rSges(7,2) * t351;
t280 = rSges(7,3) * t343 + t345 * t392;
t320 = pkin(5) * t345 + pkin(8) * t343;
t419 = t280 + t320;
t232 = t419 * t344;
t234 = t419 * t346;
t541 = t232 * t344 + t234 * t346;
t557 = -m(6) / 0.2e1;
t435 = t318 * t321 * t557 + t541 * t526;
t210 = ((rSges(7,3) + pkin(8)) * t343 + (pkin(5) + t392) * t345) * t346;
t525 = m(7) / 0.2e1;
t453 = t343 * t344;
t256 = rSges(7,3) * t453 + (rSges(7,1) * t449 - rSges(7,2) * t450) * t345;
t407 = -pkin(5) * t451 - pkin(8) * t453 - t256;
t554 = t407 * t344;
t290 = t318 * t344;
t291 = t318 * t346;
t224 = t290 * t344 + t291 * t346;
t527 = m(6) / 0.2e1;
t559 = t224 * t527;
t436 = (t210 * t346 - t554) * t525 + t559;
t53 = t436 - t435;
t562 = t53 * qJD(1);
t379 = -t216 * t297 - t218 * t298;
t432 = t295 * t214 + t296 * t217;
t561 = t379 + t432 + (-t211 * t344 - t213 * t346) * t345;
t222 = rSges(7,3) * t448 - t540;
t560 = -t222 * t343 + t280 * t448;
t111 = -t213 * t451 + t295 * t216 - t218 * t296;
t385 = Icges(7,5) * t353 - Icges(7,6) * t351;
t274 = Icges(7,3) * t343 + t345 * t385;
t465 = Icges(7,4) * t353;
t388 = -Icges(7,2) * t351 + t465;
t276 = Icges(7,6) * t343 + t345 * t388;
t466 = Icges(7,4) * t351;
t390 = Icges(7,1) * t353 - t466;
t278 = Icges(7,5) * t343 + t345 * t390;
t151 = t274 * t448 + t276 * t297 - t278 * t298;
t555 = t343 * t151;
t469 = Icges(6,4) * t343;
t389 = Icges(6,2) * t345 + t469;
t268 = Icges(6,6) * t346 + t344 * t389;
t322 = Icges(6,4) * t451;
t270 = Icges(6,1) * t453 + Icges(6,5) * t346 + t322;
t376 = -t268 * t345 - t270 * t343;
t553 = t346 * t376;
t269 = -Icges(6,6) * t344 + t346 * t389;
t468 = Icges(6,4) * t345;
t391 = Icges(6,1) * t343 + t468;
t271 = -Icges(6,5) * t344 + t346 * t391;
t283 = -Icges(6,2) * t453 + t322;
t313 = -Icges(6,2) * t343 + t468;
t284 = t313 * t346;
t315 = Icges(6,1) * t345 - t469;
t285 = t315 * t344;
t286 = t315 * t346;
t551 = -t343 * (t344 * (t269 - t286) - t346 * (t268 - t285)) + t345 * (t344 * (t271 + t284) - t346 * (t270 + t283));
t513 = -t344 / 0.2e1;
t512 = t344 / 0.2e1;
t509 = t346 / 0.2e1;
t549 = t111 * t344;
t548 = t111 * t346;
t547 = t280 * t346;
t545 = (t269 * t345 + t271 * t343) * t346;
t543 = t346 * t166 - t344 * t552;
t183 = t220 * t343 + t280 * t451;
t539 = -t183 * t346 + t344 * t560;
t437 = (t210 * t344 + t346 * t407) * t525 + (-t290 * t346 + t291 * t344) * t527;
t273 = Icges(7,3) * t345 - t343 * t385;
t442 = t353 * t278;
t444 = t351 * t276;
t374 = t442 - t444;
t365 = t273 - t374;
t455 = t274 * t343;
t537 = t345 * t365 - t455;
t366 = -t274 * t346 + t558;
t536 = t345 * t366 - t457;
t380 = -t214 * t351 + t217 * t353;
t367 = t274 * t344 - t380;
t459 = t211 * t343;
t535 = t345 * t367 - t459;
t300 = (-Icges(7,2) * t353 - t466) * t345;
t301 = (-Icges(7,1) * t351 - t465) * t345;
t532 = -(t278 / 0.2e1 + t300 / 0.2e1) * t351 + t353 * (t301 / 0.2e1 - t276 / 0.2e1);
t531 = 0.2e1 * m(7);
t530 = 0.2e1 * t321;
t529 = 0.4e1 * qJD(1);
t528 = 2 * qJD(5);
t149 = -t274 * t451 + t276 * t295 + t278 * t296;
t144 = t149 * t343;
t110 = -t211 * t451 + t432;
t383 = -t110 * t344 + t548;
t42 = t345 * t383 + t144;
t524 = -t42 / 0.2e1;
t238 = Icges(7,5) * t297 + Icges(7,6) * t298;
t428 = Icges(7,2) * t298 - t218 + t288;
t430 = -Icges(7,1) * t297 + t216 - t289;
t97 = t238 * t343 + (-t351 * t428 - t353 * t430) * t345;
t523 = t97 / 0.2e1;
t279 = -t343 * t392 + t478;
t132 = (t279 * t344 + t220) * t345 + (-t280 * t344 + t256) * t343;
t133 = (t279 * t346 - t222) * t345;
t521 = m(7) * (t132 * t166 + t133 * t552 - t183 * t407 + t210 * t560);
t302 = (-rSges(7,1) * t351 - rSges(7,2) * t353) * t345;
t519 = m(7) * (-t232 * t244 - t234 * t243 - t302 * t543);
t516 = m(7) * (-t166 * t407 + t210 * t552);
t515 = t321 / 0.2e1;
t514 = t343 / 0.2e1;
t511 = t344 / 0.4e1;
t510 = t345 / 0.2e1;
t508 = t346 / 0.4e1;
t507 = m(4) * ((rSges(4,3) * t346 - pkin(2) * t344 + t400) * t346 + ((rSges(4,3) + qJ(3)) * t344 + t396) * t344);
t395 = rSges(5,1) * t349 + rSges(5,2) * cos(pkin(10));
t412 = rSges(5,3) + pkin(2) + qJ(4);
t235 = -t344 * t412 + t346 * t395 + t400;
t236 = t486 + t412 * t346 + (qJ(3) + t395) * t344;
t506 = m(5) * (-t235 * t344 + t346 * t236);
t505 = m(5) * (t235 * t346 + t344 * t236);
t394 = rSges(6,1) * t343 + rSges(6,2) * t345;
t360 = -t344 * rSges(6,3) + t346 * t394;
t205 = t360 + t361;
t206 = t486 - t333 + (rSges(6,3) + pkin(2)) * t346 + (t394 + t399) * t344;
t504 = m(6) * (t205 * t291 + t206 * t290);
t503 = m(6) * (-t205 * t344 + t346 * t206);
t502 = m(6) * (t205 * t346 + t344 * t206);
t498 = m(7) * t543;
t497 = m(7) * (t344 * t166 + t346 * t552);
t496 = m(7) * t539;
t495 = m(7) * (t183 * t344 + t346 * t560);
t491 = m(7) * (t232 * t346 - t234 * t344);
t177 = -t243 * t346 - t244 * t344;
t488 = m(7) * t177;
t179 = t243 * t344 - t244 * t346;
t487 = m(7) * t179;
t483 = qJD(5) / 0.2e1;
t475 = t344 * t42;
t113 = t213 * t448 - t379;
t382 = t113 * t346 - t564;
t43 = t345 * t382 + t555;
t474 = t346 * t43;
t473 = -t110 + t561;
t470 = t113 + t561;
t299 = (-Icges(7,5) * t351 - Icges(7,6) * t353) * t345;
t454 = t343 * t299;
t452 = t343 * t346;
t275 = Icges(7,6) * t345 - t343 * t388;
t445 = t351 * t275;
t277 = Icges(7,5) * t345 - t343 * t390;
t443 = t353 * t277;
t81 = t132 * t344 + t133 * t346;
t441 = t81 * qJD(4);
t440 = t81 * qJD(6);
t431 = -Icges(7,1) * t295 + t214 + t467;
t429 = -Icges(7,2) * t296 + t217 + t287;
t422 = t276 - t301;
t421 = t278 + t300;
t420 = pkin(8) * t345 + t279 - t484;
t417 = qJD(1) * t343;
t416 = qJD(1) * t345;
t415 = qJD(6) * t345;
t188 = (t526 + t557 - m(5) / 0.2e1) * t530;
t413 = t188 * qJD(1);
t13 = -t555 + (t346 * t473 + t564) * t345;
t409 = -t13 / 0.2e1 - t43 / 0.2e1;
t12 = t144 + (-t344 * t470 + t548) * t345;
t408 = t524 + t12 / 0.2e1;
t386 = Icges(6,5) * t343 + Icges(6,6) * t345;
t155 = t346 * (Icges(6,3) * t346 + t344 * t386) + t268 * t451 + t270 * t453;
t267 = -Icges(6,3) * t344 + t346 * t386;
t156 = -t346 * t267 - t269 * t451 - t271 * t453;
t404 = -t451 / 0.4e1;
t403 = t448 / 0.4e1;
t398 = t394 * t321;
t122 = t295 * t421 - t296 * t422 - t299 * t451;
t123 = t297 * t421 + t298 * t422 + t299 * t448;
t237 = Icges(7,5) * t295 - Icges(7,6) * t296;
t96 = t237 * t343 + (-t351 * t429 - t353 * t431) * t345;
t397 = t519 / 0.2e1 + (t123 + t97) * t511 + (t122 + t96) * t508;
t387 = Icges(6,5) * t345 - Icges(6,6) * t343;
t363 = m(7) * (-t132 * t346 + t133 * t344);
t373 = m(7) * t302 * t515;
t67 = -t363 / 0.2e1 + t373;
t377 = t220 * t346 + t222 * t344;
t109 = (-t256 * t346 + t344 * t547) * t345 + t377 * t343;
t87 = (t109 / 0.4e1 + t179 / 0.4e1) * t531;
t384 = t87 * qJD(2) - t67 * qJD(3);
t125 = t345 * t380 + t459;
t381 = -t125 * t344 - t126 * t346;
t158 = -t267 * t344 + t545;
t17 = t346 * t470 + t549;
t55 = t110 * t346 + t549;
t372 = -t155 * t346 / 0.2e1 + t156 * t513 + (t156 - t553) * t512 + (t158 - t545 + (t267 + t376) * t344 + t155) * t509 - t55 / 0.2e1 + t17 / 0.2e1;
t18 = t344 * t473 - t563;
t56 = t113 * t344 + t563;
t371 = t267 * t341 / 0.2e1 + t158 * t512 + t18 / 0.2e1 + t56 / 0.2e1 + (t156 + (t267 - t376) * t346 + t553) * t509;
t82 = -t237 * t451 + t295 * t429 - t296 * t431;
t83 = -t238 * t451 + t295 * t428 - t296 * t430;
t35 = t344 * t83 + t346 * t82;
t84 = t237 * t448 + t297 * t429 + t298 * t431;
t85 = t238 * t448 + t297 * t428 + t298 * t430;
t36 = t344 * t85 + t346 * t84;
t370 = t35 * t509 + t36 * t512;
t359 = t12 * t511 + t13 * t508 + t17 * t403 - t475 / 0.4e1 + t474 / 0.4e1 - t55 * t448 / 0.4e1 + (t18 + t56) * t404;
t104 = t275 * t295 + t277 * t296 - t344 * t537;
t105 = t275 * t297 - t277 * t298 + t346 * t537;
t128 = (t274 + t443 - t445) * t345 + t365 * t343;
t186 = t345 * t374 + t455;
t252 = t276 * t344;
t254 = t278 * t344;
t91 = (-t252 * t351 + t254 * t353 + t211) * t345 + t367 * t343;
t253 = t276 * t346;
t255 = t278 * t346;
t92 = (t253 * t351 - t255 * t353 + t213) * t345 + t366 * t343;
t358 = t128 * t514 + t186 * t510 + t521 / 0.2e1 + (t125 + t149) * t453 / 0.4e1 - (-t126 + t151) * t452 / 0.4e1 + (t104 + t91) * t404 + (t105 + t92) * t403;
t357 = t443 / 0.2e1 - t445 / 0.2e1 + t274 / 0.2e1 - t391 / 0.2e1 - t313 / 0.2e1;
t356 = t442 / 0.2e1 - t444 / 0.2e1 + t315 / 0.2e1 - t389 / 0.2e1 - t273 / 0.2e1;
t282 = t387 * t346;
t281 = t344 * t387;
t233 = t420 * t346;
t231 = t420 * t344;
t190 = -t244 * t343 + t302 * t448;
t189 = t243 * t343 + t302 * t451;
t187 = (-m(7) / 0.4e1 - m(6) / 0.4e1 - m(5) / 0.4e1) * t530 + (m(5) + m(6) + m(7)) * t515;
t175 = t487 / 0.2e1;
t173 = t488 / 0.2e1;
t163 = t177 * t345;
t160 = t491 / 0.2e1;
t145 = t377 * t345;
t136 = t554 + (-t320 * t346 - t547) * t346;
t135 = (t454 + (-t351 * t421 - t353 * t422) * t345) * t343;
t131 = (-pkin(5) * t452 + t222 + t328) * t346 + (-pkin(5) * t453 - t542) * t344;
t115 = t495 / 0.2e1;
t114 = -t496 / 0.2e1;
t88 = (t109 - t179) * t525;
t79 = m(7) * t81 * t483;
t78 = -t253 * t297 + t255 * t298 + t346 * t536;
t77 = t252 * t297 - t254 * t298 + t346 * t535;
t76 = -t253 * t295 - t255 * t296 - t344 * t536;
t75 = t252 * t295 + t254 * t296 - t344 * t535;
t66 = t363 / 0.2e1 + t373;
t65 = t160 - t437;
t64 = -t491 / 0.2e1 + t437;
t63 = t160 + t437;
t62 = -t131 * t179 + t302 * t541;
t52 = t435 + t436;
t50 = t565 + t454 / 0.2e1 + t532 * t345;
t47 = t498 + t503 + t506;
t46 = t186 * t343 + t345 * t381;
t44 = t497 + t502 + t505 + t507;
t34 = t115 - t488 / 0.2e1;
t33 = t114 - t487 / 0.2e1;
t32 = t175 + t114;
t31 = t175 + t496 / 0.2e1;
t30 = t173 + t115;
t29 = t173 - t495 / 0.2e1;
t28 = t344 * t78 + t346 * t77;
t27 = t344 * t76 + t346 * t75;
t26 = -t343 * t356 + t345 * t357 + t504 + t516;
t23 = -t109 * t145 + t132 * t183 + t133 * t560;
t22 = t123 * t343 + (-t344 * t84 + t346 * t85) * t345;
t21 = t122 * t343 + (-t344 * t82 + t346 * t83) * t345;
t14 = (-t91 * t344 + t92 * t346 + t186) * t345 + (t128 - t381) * t343;
t9 = (-t344 * t77 + t346 * t78 + t151) * t345 + (t105 - t382) * t343;
t8 = (-t344 * t75 + t346 * t76 + t149) * t345 + (t104 - t383) * t343;
t7 = m(7) * t62 + t370;
t6 = (t344 * t409 + t346 * t408) * t345;
t5 = t344 * t372 + t346 * t371;
t4 = m(7) * t23 + (t8 * t513 + t9 * t509 + t46 / 0.2e1) * t345 + (t475 / 0.2e1 - t474 / 0.2e1 + t14 / 0.2e1) * t343;
t3 = t358 + (t42 / 0.4e1 - t12 / 0.4e1 + (t56 / 0.4e1 + t18 / 0.4e1) * t345) * t344 + (-t43 / 0.4e1 - t13 / 0.4e1 + (-t17 / 0.4e1 + t55 / 0.4e1) * t345) * t346 + t397;
t2 = t358 - t519 / 0.2e1 + t359 + (-t122 / 0.4e1 - t96 / 0.4e1) * t346 + (-t123 / 0.4e1 - t97 / 0.4e1) * t344;
t1 = -t521 / 0.2e1 + t359 + (-t128 / 0.2e1 + (t151 / 0.4e1 - t126 / 0.4e1) * t346 + (-t149 / 0.4e1 - t125 / 0.4e1) * t344) * t343 + (-t186 / 0.2e1 + (-t105 / 0.4e1 - t92 / 0.4e1) * t346 + (t91 / 0.4e1 + t104 / 0.4e1) * t344) * t345 + t397;
t10 = [t44 * qJD(3) + t47 * qJD(4) + t26 * qJD(5) + t50 * qJD(6), 0, qJD(1) * t44 + qJD(4) * t187 + qJD(5) * t63 + qJD(6) * t30, qJD(1) * t47 + qJD(3) * t187 + qJD(5) * t52 + qJD(6) * t32, t26 * qJD(1) + t63 * qJD(3) + t52 * qJD(4) + t3 * qJD(6) + (m(7) * (-t166 * t233 + t210 * t232 + t231 * t552 + t234 * t407) + (m(6) * (t206 * t394 - t290 * t318) + t104 / 0.2e1 + t91 / 0.2e1 - t386 * t509 + (-t268 / 0.2e1 + t285 / 0.2e1) * t345 + (-t270 / 0.2e1 - t283 / 0.2e1) * t343 - t371) * t346 + (m(6) * (-t205 * t394 + t291 * t318) + t105 / 0.2e1 + t92 / 0.2e1 - t386 * t512 + (t269 / 0.2e1 - t286 / 0.2e1) * t345 + (t271 / 0.2e1 + t284 / 0.2e1) * t343 - t372) * t344) * qJD(5), t50 * qJD(1) + t30 * qJD(3) + t32 * qJD(4) + t3 * qJD(5) + t135 * qJD(6) + (t166 * t189 + t183 * t243 + t190 * t552 - t244 * t560) * t480 + ((t123 / 0.2e1 + t523 - t408) * t346 + (-t122 / 0.2e1 - t96 / 0.2e1 - t409) * t344) * t415; 0, 0, 0, 0 (t136 * t525 - t559) * t528 + t88 * qJD(6), t88 * qJD(5) + t163 * t480; t188 * qJD(4) + t64 * qJD(5) + t29 * qJD(6) + (-t497 / 0.4e1 - t502 / 0.4e1 - t507 / 0.4e1 - t505 / 0.4e1) * t529, 0, 0, t413, t64 * qJD(1) + ((t231 * t344 + t233 * t346) * t525 - t398 * t527) * t528 + t66 * qJD(6), t29 * qJD(1) + t66 * qJD(5) + (-t189 * t346 + t190 * t344) * t480; -t188 * qJD(3) + t53 * qJD(5) + t31 * qJD(6) + (-t498 / 0.4e1 - t503 / 0.4e1 - t506 / 0.4e1) * t529, 0, -t413, 0, t562 + ((t231 * t346 - t233 * t344) * t483 + t440 / 0.4e1) * t531, t31 * qJD(1) + t79 + (t189 * t344 + t190 * t346) * t480; t65 * qJD(3) - t53 * qJD(4) + t5 * qJD(5) + t1 * qJD(6) + (-t516 / 0.4e1 - t504 / 0.4e1) * t529 - t357 * t416 + t356 * t417, -qJD(6) * t87, qJD(1) * t65 + qJD(6) * t67, t440 * t526 - t562, t5 * qJD(1) + (m(6) * (-(-t346 * t360 + (-t346 * rSges(6,3) - t344 * t394) * t344) * t224 - t318 * t398) + m(7) * (t131 * t136 + t231 * t232 + t233 * t234) + (-t341 * t282 + (t344 * t281 + t551) * t346 + t28) * t512 + (t342 * t281 + (-t346 * t282 - t551) * t344 + t27) * t509) * qJD(5) + t7 * qJD(6), t1 * qJD(1) + t7 * qJD(5) + t441 * t526 + (-t46 / 0.2e1 + (t36 / 0.2e1 - t9 / 0.2e1) * t346 + (-t35 / 0.2e1 + t8 / 0.2e1) * t344) * t415 - t384 + (t21 * t509 + t22 * t512 + (-t14 / 0.2e1 + (t96 / 0.2e1 + t43 / 0.2e1) * t346 + (t523 + t524) * t344) * t343 + (t131 * t163 + t145 * t179 - t189 * t234 + t190 * t232 + t302 * t539 - t23) * m(7)) * qJD(6); -t299 * t417 / 0.2e1 + t34 * qJD(3) + t33 * qJD(4) + t2 * qJD(5) + t6 * qJD(6) - qJD(1) * t565 - t532 * t416, qJD(5) * t87, qJD(1) * t34 - qJD(5) * t67, qJD(1) * t33 + t79, t2 * qJD(1) + (t8 * t509 + t55 * t453 / 0.2e1 - t27 * t451 / 0.2e1 + t9 * t512 - t56 * t452 / 0.2e1 + t28 * t448 / 0.2e1 + (t125 * t346 - t126 * t344) * t510 + (t92 * t344 + t91 * t346) * t514 - t370) * qJD(5) + t4 * qJD(6) + (t441 / 0.2e1 + (t109 * t131 - t132 * t234 + t133 * t232 - t136 * t145 - t183 * t233 + t231 * t560 - t62) * qJD(5)) * m(7) + t384, t6 * qJD(1) + t4 * qJD(5) + (t135 * t514 + m(7) * (-t145 * t163 + t183 * t189 + t190 * t560) + (t21 * t513 + t22 * t509 + (-t344 * t96 + t346 * t97) * t514) * t345) * qJD(6);];
Cq  = t10;
