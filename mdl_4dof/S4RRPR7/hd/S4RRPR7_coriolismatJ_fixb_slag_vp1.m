% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:58
% EndTime: 2019-12-31 17:06:17
% DurationCPUTime: 13.71s
% Computational Cost: add. (28224->558), mult. (38279->818), div. (0->0), fcn. (41677->8), ass. (0->330)
t344 = qJ(2) + pkin(7);
t332 = sin(t344);
t349 = sin(qJ(1));
t470 = t332 * t349;
t586 = Icges(3,3) + Icges(4,3);
t333 = cos(t344);
t347 = sin(qJ(4));
t350 = cos(qJ(4));
t385 = Icges(5,5) * t350 - Icges(5,6) * t347;
t233 = -Icges(5,3) * t333 + t332 * t385;
t485 = Icges(5,4) * t350;
t388 = -Icges(5,2) * t347 + t485;
t235 = -Icges(5,6) * t333 + t332 * t388;
t486 = Icges(5,4) * t347;
t389 = Icges(5,1) * t350 - t486;
t237 = -Icges(5,5) * t333 + t332 * t389;
t352 = cos(qJ(1));
t456 = t350 * t352;
t461 = t349 * t347;
t283 = t333 * t461 + t456;
t453 = t352 * t347;
t460 = t349 * t350;
t284 = t333 * t460 - t453;
t125 = t233 * t470 - t235 * t283 + t237 * t284;
t186 = Icges(5,5) * t284 - Icges(5,6) * t283 + Icges(5,3) * t470;
t276 = Icges(5,4) * t284;
t189 = -Icges(5,2) * t283 + Icges(5,6) * t470 + t276;
t275 = Icges(5,4) * t283;
t193 = -Icges(5,1) * t284 - Icges(5,5) * t470 + t275;
t102 = t186 * t470 - t189 * t283 - t193 * t284;
t285 = -t333 * t453 + t460;
t286 = t333 * t456 + t461;
t469 = t332 * t352;
t188 = Icges(5,5) * t286 + Icges(5,6) * t285 + Icges(5,3) * t469;
t487 = Icges(5,4) * t286;
t191 = Icges(5,2) * t285 + Icges(5,6) * t469 + t487;
t277 = Icges(5,4) * t285;
t194 = Icges(5,1) * t286 + Icges(5,5) * t469 + t277;
t103 = t188 * t470 - t283 * t191 + t284 * t194;
t384 = t102 * t349 + t103 * t352;
t12 = t125 * t333 - t332 * t384;
t295 = Icges(4,5) * t333 - Icges(4,6) * t332;
t348 = sin(qJ(2));
t351 = cos(qJ(2));
t305 = Icges(3,5) * t351 - Icges(3,6) * t348;
t584 = (t295 + t305) * t352 + t586 * t349;
t459 = t349 * t351;
t463 = t348 * t349;
t467 = t333 * t349;
t577 = -Icges(3,5) * t459 - Icges(4,5) * t467 + Icges(3,6) * t463 + Icges(4,6) * t470 + t586 * t352;
t127 = t233 * t469 + t285 * t235 + t286 * t237;
t104 = t186 * t469 + t285 * t189 - t286 * t193;
t105 = t188 * t469 + t285 * t191 + t286 * t194;
t383 = t349 * t104 + t105 * t352;
t583 = -t127 * t333 + t332 * t383;
t54 = -t102 * t352 + t103 * t349;
t394 = t286 * rSges(5,1) + t285 * rSges(5,2);
t197 = rSges(5,3) * t469 + t394;
t500 = rSges(5,1) * t350;
t393 = -rSges(5,2) * t347 + t500;
t241 = -rSges(5,3) * t333 + t332 * t393;
t147 = t333 * t197 + t241 * t469;
t548 = -t284 * rSges(5,1) + t283 * rSges(5,2);
t195 = rSges(5,3) * t470 - t548;
t578 = t195 * t333 + t241 * t470;
t107 = -t147 * t349 + t352 * t578;
t477 = t186 * t333;
t573 = t189 * t347 + t193 * t350;
t112 = t573 * t332 + t477;
t302 = pkin(3) * t332 - pkin(6) * t333;
t507 = pkin(2) * t348;
t397 = t241 + t302 + t507;
t182 = t397 * t349;
t184 = t397 * t352;
t550 = t349 * t182 + t184 * t352;
t300 = rSges(4,1) * t332 + rSges(4,2) * t333;
t366 = t300 + t507;
t551 = t366 * t352;
t552 = t366 * t349;
t559 = t349 * t552 + t352 * t551;
t572 = -m(5) / 0.2e1;
t448 = t550 * t572 - m(4) * t559 / 0.2e1;
t520 = rSges(5,3) + pkin(6);
t174 = (t507 - t520 * t333 + (pkin(3) + t393) * t332) * t349;
t466 = t333 * t352;
t320 = pkin(6) * t466;
t427 = t332 * rSges(5,2) * t453 + rSges(5,3) * t466;
t175 = t320 + (-t507 + (-pkin(3) - t500) * t332) * t352 + t427;
t540 = m(5) / 0.2e1;
t541 = m(4) / 0.2e1;
t450 = (t349 * t174 - t175 * t352) * t540 + t559 * t541;
t38 = t450 - t448;
t581 = t38 * qJD(1);
t488 = Icges(4,4) * t332;
t299 = Icges(4,1) * t333 - t488;
t256 = Icges(4,5) * t349 + t299 * t352;
t489 = Icges(3,4) * t348;
t309 = Icges(3,1) * t351 - t489;
t265 = Icges(3,5) * t349 + t309 * t352;
t580 = -t256 * t467 - t265 * t459;
t579 = -t104 * t352 + t105 * t349;
t576 = t584 * t352 + t580;
t316 = Icges(4,4) * t470;
t255 = Icges(4,1) * t467 - Icges(4,5) * t352 - t316;
t324 = Icges(3,4) * t463;
t264 = Icges(3,1) * t459 - Icges(3,5) * t352 - t324;
t455 = t351 * t352;
t575 = -t255 * t466 - t264 * t455 + t577 * t349;
t562 = t256 * t466 + t265 * t455 + t584 * t349;
t253 = Icges(4,4) * t467 - Icges(4,2) * t470 - Icges(4,6) * t352;
t262 = Icges(3,4) * t459 - Icges(3,2) * t463 - Icges(3,6) * t352;
t574 = t253 * t332 + t262 * t348;
t504 = -qJ(3) - pkin(5);
t328 = t349 * t504;
t506 = pkin(2) * t351;
t331 = pkin(1) + t506;
t505 = pkin(3) * t333;
t545 = -t332 * t520 - t505;
t156 = -t328 + (t331 - t545) * t352 + t394;
t426 = -t349 * t331 - t352 * t504;
t564 = t545 * t349 + t426 + t548;
t110 = t156 * t349 + t352 * t564;
t327 = Icges(4,4) * t333;
t482 = Icges(4,2) * t332;
t254 = Icges(4,6) * t349 + (t327 - t482) * t352;
t340 = Icges(3,4) * t351;
t307 = -Icges(3,2) * t348 + t340;
t263 = Icges(3,6) * t349 + t307 * t352;
t462 = t348 * t352;
t565 = -t254 * t469 - t263 * t462 + t562;
t563 = t254 * t332 + t263 * t348 + t577;
t234 = Icges(5,3) * t332 + t333 * t385;
t306 = Icges(3,2) * t351 + t489;
t458 = t350 * t237;
t465 = t347 * t235;
t490 = Icges(4,1) * t332;
t561 = -(t309 / 0.2e1 - t306 / 0.2e1) * t348 - (t327 + t490 / 0.2e1 - t482 / 0.2e1 + t458 / 0.2e1 - t465 / 0.2e1 - t234 / 0.2e1) * t333;
t560 = -t253 * t469 - t254 * t470 - t262 * t462 - t263 * t463 - t575 - t576;
t526 = t349 / 0.2e1;
t523 = -t352 / 0.2e1;
t521 = t352 / 0.2e1;
t554 = t241 * t349;
t386 = -Icges(4,5) * t332 - Icges(4,6) * t333;
t387 = -Icges(3,5) * t348 - Icges(3,6) * t351;
t549 = (t386 + t387) * t352;
t308 = Icges(3,1) * t348 + t340;
t345 = t349 ^ 2;
t346 = t352 ^ 2;
t424 = t345 + t346;
t269 = (-Icges(5,2) * t350 - t486) * t332;
t272 = (-Icges(5,1) * t347 - t485) * t332;
t544 = -t347 * (t237 / 0.2e1 + t269 / 0.2e1) + t350 * (t272 / 0.2e1 - t235 / 0.2e1);
t543 = 0.4e1 * qJD(1);
t542 = 2 * qJD(2);
t539 = -t583 / 0.2e1;
t538 = t54 / 0.2e1;
t537 = t579 / 0.2e1;
t202 = -Icges(5,5) * t283 - Icges(5,6) * t284;
t443 = -Icges(5,2) * t284 - t193 - t275;
t445 = -Icges(5,1) * t283 - t189 - t276;
t84 = -t202 * t333 + (-t347 * t443 + t350 * t445) * t332;
t536 = t84 / 0.2e1;
t242 = rSges(5,3) * t332 + t333 * t393;
t114 = (t242 * t349 - t195) * t332;
t220 = -rSges(5,1) * t332 * t456 + t427;
t115 = (-t241 * t352 - t220) * t333 + (-t242 * t352 + t197) * t332;
t533 = m(5) * (t114 * t564 + t115 * t156 - t147 * t175 + t174 * t578);
t208 = -rSges(5,1) * t283 - rSges(5,2) * t284;
t209 = rSges(5,1) * t285 - rSges(5,2) * t286;
t278 = (-rSges(5,1) * t347 - rSges(5,2) * t350) * t332;
t531 = m(5) * (-t110 * t278 - t182 * t209 + t184 * t208);
t529 = m(5) * (t156 * t175 + t174 * t564);
t528 = t332 / 0.2e1;
t527 = -t333 / 0.2e1;
t525 = t349 / 0.4e1;
t522 = -t352 / 0.4e1;
t343 = t352 * pkin(5);
t502 = rSges(3,1) * t351;
t410 = pkin(1) + t502;
t425 = rSges(3,2) * t463 + t352 * rSges(3,3);
t229 = -t349 * t410 + t343 + t425;
t326 = rSges(3,2) * t462;
t230 = -t326 + t410 * t352 + (rSges(3,3) + pkin(5)) * t349;
t310 = rSges(3,1) * t348 + rSges(3,2) * t351;
t293 = t310 * t349;
t294 = t310 * t352;
t519 = m(3) * (t229 * t293 - t230 * t294);
t374 = rSges(4,1) * t467 - rSges(4,2) * t470 - t352 * rSges(4,3);
t210 = -t374 + t426;
t408 = -rSges(4,2) * t469 + t349 * rSges(4,3);
t501 = rSges(4,1) * t333;
t211 = -t328 + (t331 + t501) * t352 + t408;
t518 = m(4) * (t210 * t552 - t211 * t551);
t517 = m(4) * (t210 * t352 + t211 * t349);
t513 = m(5) * t107;
t512 = m(5) * t110;
t149 = t349 * t208 + t209 * t352;
t508 = m(5) * t149;
t503 = m(5) * qJD(4);
t497 = t349 * t12;
t496 = t352 * t583;
t476 = t188 * t333;
t475 = t233 * t333;
t266 = (-Icges(5,5) * t347 - Icges(5,6) * t350) * t332;
t468 = t333 * t266;
t236 = Icges(5,6) * t332 + t333 * t388;
t464 = t347 * t236;
t238 = Icges(5,5) * t332 + t333 * t389;
t457 = t350 * t238;
t66 = t114 * t349 - t115 * t352;
t452 = t66 * qJD(3);
t444 = Icges(5,1) * t285 - t191 - t487;
t442 = -Icges(5,2) * t286 + t194 + t277;
t437 = -t349 * (pkin(1) * t349 - t343 + t426) + t352 * (-t349 * pkin(5) - t328 + (-pkin(1) + t331) * t352);
t436 = -t235 + t272;
t435 = t237 + t269;
t390 = -t327 - t490;
t273 = t390 * t349;
t432 = t253 - t273;
t270 = -Icges(4,2) * t467 - t316;
t431 = t255 + t270;
t291 = t308 * t349;
t430 = t262 + t291;
t289 = -Icges(3,2) * t459 - t324;
t429 = t264 + t289;
t423 = qJD(1) * t332;
t422 = qJD(4) * t332;
t421 = t66 * t540;
t418 = t539 + t583 / 0.2e1;
t417 = t470 / 0.4e1;
t411 = t295 / 0.2e1 + t305 / 0.2e1;
t409 = rSges(4,2) * t332 - t501 - t506;
t274 = t390 * t352;
t407 = (-t254 + t274) * t349;
t296 = Icges(4,2) * t333 + t488;
t271 = t296 * t352;
t406 = (-t256 + t271) * t349;
t292 = t308 * t352;
t405 = (-t263 - t292) * t349;
t290 = t306 * t352;
t404 = (-t265 + t290) * t349;
t399 = t424 * t507;
t203 = Icges(5,5) * t285 - Icges(5,6) * t286;
t85 = -t203 * t333 + (-t347 * t442 + t350 * t444) * t332;
t98 = t266 * t470 - t283 * t435 + t284 * t436;
t99 = t266 * t469 + t285 * t435 + t286 * t436;
t398 = t531 / 0.2e1 + (t85 + t99) * t525 + (t84 + t98) * t522;
t303 = pkin(6) * t332 + t505;
t396 = -t242 - t303 - t506;
t380 = -t191 * t347 + t194 * t350;
t113 = t332 * t380 - t476;
t382 = -t112 * t349 + t113 * t352;
t379 = t195 * t352 - t197 * t349;
t378 = t458 - t465;
t375 = -m(5) * (t156 * t209 - t208 * t564) + t468 / 0.2e1;
t74 = t202 * t470 - t283 * t443 + t284 * t445;
t75 = t203 * t470 - t283 * t442 + t284 * t444;
t33 = t75 * t349 - t352 * t74;
t76 = t202 * t469 + t285 * t443 + t286 * t445;
t77 = t203 * t469 + t285 * t442 + t286 * t444;
t34 = t77 * t349 - t352 * t76;
t373 = t33 * t523 + t34 * t526;
t370 = -t233 * t349 + t573;
t369 = -t233 * t352 - t380;
t368 = t234 - t378;
t365 = t12 * t525 + t583 * t522 - t497 / 0.4e1 + t496 / 0.4e1 + (t417 - t470 / 0.4e1) * t579;
t364 = -t579 / 0.2e1 + t537 - t562 * t349 / 0.2e1 + t565 * t526 + ((t574 + t584) * t352 + t560 + t575 + t580) * t523;
t363 = t538 - t54 / 0.2e1 + (t577 * t352 + (t255 * t333 + t264 * t351 - t574) * t349) * t523 + (t563 * t352 - t562 + t565) * t521 + (t563 * t349 + t560 + t576) * t526;
t136 = t332 * t378 - t475;
t214 = t235 * t349;
t216 = t237 * t349;
t78 = -t370 * t333 + (t214 * t347 - t216 * t350 + t186) * t332;
t215 = t235 * t352;
t217 = t237 * t352;
t79 = -t369 * t333 + (t215 * t347 - t217 * t350 + t188) * t332;
t355 = t332 * t368 + t475;
t89 = -t236 * t283 + t238 * t284 + t349 * t355;
t90 = t285 * t236 + t286 * t238 + t352 * t355;
t94 = -t368 * t333 + (t233 + t457 - t464) * t332;
t358 = t136 * t528 + t533 / 0.2e1 + t94 * t527 + (t78 + t89) * t417 + (t79 + t90) * t469 / 0.4e1 + (-t112 + t125) * t467 / 0.4e1 + (t113 + t127) * t466 / 0.4e1;
t357 = t332 * t370 + t477;
t356 = t332 * t369 + t476;
t353 = t299 / 0.2e1 - t296 / 0.2e1 + t457 / 0.2e1 - t464 / 0.2e1 + t233 / 0.2e1;
t312 = -rSges(3,2) * t348 + t502;
t287 = t387 * t349;
t267 = t386 * t349;
t249 = t409 * t352;
t247 = t409 * t349;
t185 = t396 * t352;
t183 = t396 * t349;
t162 = -t333 * t209 - t278 * t469;
t161 = t208 * t333 + t278 * t470;
t144 = -t508 / 0.2e1;
t137 = (t208 * t352 - t209 * t349) * t332;
t133 = t379 * t332;
t121 = (-pkin(3) * t469 + t220 + t320) * t352 - t399 + (-t302 * t349 - t554) * t349;
t118 = -t468 + (-t347 * t435 + t350 * t436) * t332;
t106 = t513 / 0.2e1;
t95 = t379 * t333 + (-t220 * t349 - t352 * t554) * t332;
t92 = (t303 * t349 + t195) * t349 + (t303 * t352 + t197) * t352 + t437;
t81 = t512 + t517;
t70 = -t285 * t215 - t286 * t217 + t352 * t356;
t69 = -t285 * t214 - t286 * t216 + t352 * t357;
t68 = t215 * t283 - t217 * t284 + t349 * t356;
t67 = t214 * t283 - t216 * t284 + t349 * t357;
t65 = qJD(2) * t421;
t49 = t544 * t332 - t375;
t48 = t92 * t149 + t550 * t278;
t45 = -t136 * t333 + t332 * t382;
t40 = t448 + t450;
t31 = t106 + t508 / 0.2e1;
t30 = t144 + t106;
t29 = t144 - t513 / 0.2e1;
t28 = t70 * t349 - t352 * t69;
t27 = t68 * t349 - t352 * t67;
t24 = -t99 * t333 + (t349 * t76 + t352 * t77) * t332;
t23 = -t98 * t333 + (t349 * t74 + t352 * t75) * t332;
t22 = t114 * t578 - t115 * t147 + t133 * t95;
t21 = (t308 / 0.2e1 + t307 / 0.2e1) * t351 + t519 + t518 + t529 + t353 * t332 - t561;
t14 = (t382 - t94) * t333 + (t78 * t349 + t79 * t352 + t136) * t332;
t9 = (t383 - t90) * t333 + (t349 * t69 + t352 * t70 + t127) * t332;
t8 = (t384 - t89) * t333 + (t349 * t67 + t352 * t68 + t125) * t332;
t7 = m(5) * t48 + t373;
t6 = t418 * t470;
t5 = m(5) * t22 + (t496 / 0.2e1 - t497 / 0.2e1 - t14 / 0.2e1) * t333 + (t9 * t521 + t8 * t526 + t45 / 0.2e1) * t332;
t4 = t349 * t363 + t352 * t364;
t3 = (-t136 / 0.2e1 + (-t90 / 0.4e1 - t79 / 0.4e1) * t352 + (-t89 / 0.4e1 - t78 / 0.4e1) * t349) * t332 + (t94 / 0.2e1 + (-t113 / 0.4e1 - t127 / 0.4e1) * t352 + (-t125 / 0.4e1 + t112 / 0.4e1) * t349) * t333 - t533 / 0.2e1 + t365 + t398;
t2 = (-t99 / 0.4e1 - t85 / 0.4e1) * t349 - t531 / 0.2e1 + (t98 / 0.4e1 + t84 / 0.4e1) * t352 + t365 + t358;
t1 = t358 + t398;
t10 = [t21 * qJD(2) + t81 * qJD(3) + t49 * qJD(4), t21 * qJD(1) + t40 * qJD(3) + t1 * qJD(4) + ((t210 * t249 + t211 * t247) * t541 + (t156 * t183 - t174 * t184 - t175 * t182 + t185 * t564) * t540) * t542 + ((-t89 / 0.2e1 - t78 / 0.2e1 + m(3) * (-t229 * t312 - t293 * t310) + t411 * t352 + (-t264 / 0.2e1 - t289 / 0.2e1) * t351 + (t262 / 0.2e1 + t291 / 0.2e1) * t348 + (-t255 / 0.2e1 - t270 / 0.2e1) * t333 + (t253 / 0.2e1 - t273 / 0.2e1) * t332 - t364) * t352 + (t90 / 0.2e1 + t79 / 0.2e1 + m(3) * (-t230 * t312 + t294 * t310) + t411 * t349 + (t265 / 0.2e1 - t290 / 0.2e1) * t351 + (-t263 / 0.2e1 - t292 / 0.2e1) * t348 + (t256 / 0.2e1 - t271 / 0.2e1) * t333 + (-t254 / 0.2e1 + t274 / 0.2e1) * t332 - t363) * t349) * qJD(2), qJD(1) * t81 + qJD(2) * t40 + qJD(4) * t30, t49 * qJD(1) + t1 * qJD(2) + t30 * qJD(3) + (-t118 * t333 + m(5) * (-t147 * t209 + t156 * t162 + t161 * t564 - t208 * t578)) * qJD(4) + ((t85 / 0.2e1 + t99 / 0.2e1) * t352 + (t536 + t98 / 0.2e1 - t418) * t349) * t422; t4 * qJD(2) - t38 * qJD(3) + t3 * qJD(4) + (-t519 / 0.4e1 - t518 / 0.4e1 - t529 / 0.4e1) * t543 - t353 * t423 + (-(t308 + t307) * t351 / 0.2e1 + t561) * qJD(1), t4 * qJD(1) + t7 * qJD(4) + (m(4) * (-t551 * t249 - t552 * t247 + (t349 * t374 + t352 * (rSges(4,1) * t466 + t408) + t437) * (-t424 * t300 - t399)) + m(3) * ((t349 * (rSges(3,1) * t459 - t425) + t352 * (rSges(3,1) * t455 + t349 * rSges(3,3) - t326)) * (-t349 * t293 - t294 * t352) + t424 * t312 * t310) + m(5) * (t121 * t92 - t182 * t183 - t184 * t185) + (t28 + (-t349 * t267 + (t352 * t432 + t407) * t333 + (t352 * t431 + t406) * t332) * t352 + (-t349 * t287 + (t352 * t430 + t405) * t351 + (t352 * t429 + t404) * t348) * t352 + t549 * t345) * t526 + (t27 + (t406 * t332 + t407 * t333 + t404 * t348 + t405 * t351 + (t431 * t332 + t432 * t333 + t429 * t348 + t430 * t351 - t549) * t352) * t349 + (t267 + t287) * t346) * t523) * qJD(2), -t581 - t66 * t503 / 0.2e1, t3 * qJD(1) + t7 * qJD(2) + t452 * t572 + (-t45 / 0.2e1 + (t34 / 0.2e1 - t9 / 0.2e1) * t352 + (t33 / 0.2e1 - t8 / 0.2e1) * t349) * t422 + (t23 * t523 + t24 * t526 + (t14 / 0.2e1 + (t536 + t539) * t352 + (-t85 / 0.2e1 + t12 / 0.2e1) * t349) * t333 + (-t107 * t278 + t133 * t149 + t137 * t92 - t161 * t184 - t162 * t182 - t22) * m(5)) * qJD(4); t38 * qJD(2) + t29 * qJD(4) + (-t517 / 0.4e1 - t512 / 0.4e1) * t543, t581 + ((-t183 * t352 + t349 * t185) * t540 + (-t247 * t352 + t349 * t249) * t541) * t542 + qJD(4) * t421, 0, t29 * qJD(1) + t65 + (t161 * t349 - t162 * t352) * t503; t375 * qJD(1) + t2 * qJD(2) + t31 * qJD(3) + t6 * qJD(4) - t544 * t423, t2 * qJD(1) + (t8 * t523 + t466 * t537 + t28 * t469 / 0.2e1 + (t112 * t352 + t113 * t349) * t528 + (t79 * t349 - t78 * t352) * t527 + t9 * t526 + t467 * t538 + t27 * t470 / 0.2e1 - t373) * qJD(2) + t5 * qJD(4) + ((-t114 * t184 - t115 * t182 + t121 * t133 - t147 * t183 + t185 * t578 + t92 * t95 - t48) * qJD(2) + t452 / 0.2e1) * m(5), qJD(1) * t31 + t65, t6 * qJD(1) + t5 * qJD(2) + (m(5) * (t133 * t137 - t147 * t162 + t161 * t578) + t333 ^ 2 * t118 / 0.2e1 + (t24 * t521 + t23 * t526 + (t84 * t349 + t85 * t352) * t527) * t332) * qJD(4);];
Cq = t10;
