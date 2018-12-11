% Calculate kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:46
% EndTime: 2018-12-10 18:09:52
% DurationCPUTime: 5.87s
% Computational Cost: add. (33927->472), mult. (34354->667), div. (0->0), fcn. (34328->30), ass. (0->217)
t480 = cos(qJ(1));
t509 = t480 / 0.2e1;
t508 = cos(qJ(5));
t471 = cos(pkin(6));
t507 = pkin(10) * t471;
t476 = sin(qJ(1));
t506 = Icges(2,4) * t476;
t467 = sin(pkin(6));
t505 = t467 * t476;
t504 = t467 * t480;
t503 = qJD(2) * t467;
t502 = pkin(8) - qJ(4);
t501 = pkin(8) + qJ(4);
t500 = V_base(5) * pkin(9) + V_base(1);
t436 = t476 * t503 + V_base(4);
t457 = V_base(6) + qJD(1);
t498 = cos(t501);
t497 = sin(t502);
t463 = pkin(6) - qJ(2);
t451 = cos(t463) / 0.2e1;
t462 = pkin(6) + qJ(2);
t456 = cos(t462);
t430 = t451 + t456 / 0.2e1;
t475 = sin(qJ(2));
t416 = -t476 * t430 - t475 * t480;
t450 = sin(t462) / 0.2e1;
t455 = sin(t463);
t428 = t450 - t455 / 0.2e1;
t479 = cos(qJ(2));
t417 = -t476 * t428 + t479 * t480;
t460 = pkin(7) + pkin(14);
t448 = sin(t460) / 0.2e1;
t461 = pkin(7) - pkin(14);
t452 = sin(t461);
t420 = t448 + t452 / 0.2e1;
t449 = cos(t461) / 0.2e1;
t453 = cos(t460);
t422 = t449 + t453 / 0.2e1;
t464 = sin(pkin(14));
t379 = t416 * t422 - t417 * t464 + t420 * t505;
t466 = sin(pkin(7));
t470 = cos(pkin(7));
t403 = -t416 * t466 + t470 * t505;
t465 = sin(pkin(8));
t469 = cos(pkin(8));
t370 = -t379 * t465 + t403 * t469;
t364 = qJD(4) * t370 + t436;
t437 = qJD(2) * t471 + t457;
t496 = cos(t502) / 0.2e1;
t495 = sin(t501) / 0.2e1;
t421 = t448 - t452 / 0.2e1;
t423 = t449 - t453 / 0.2e1;
t468 = cos(pkin(14));
t380 = t416 * t421 + t417 * t468 + t423 * t505;
t474 = sin(qJ(4));
t488 = t495 + t497 / 0.2e1;
t489 = t496 + t498 / 0.2e1;
t337 = -t379 * t489 + t380 * t474 - t403 * t488;
t324 = qJD(5) * t337 + t364;
t414 = t430 * t480 - t476 * t475;
t415 = t428 * t480 + t476 * t479;
t377 = t414 * t422 - t415 * t464 - t420 * t504;
t402 = -t414 * t466 - t470 * t504;
t369 = -t377 * t465 + t402 * t469;
t427 = t450 + t455 / 0.2e1;
t431 = t451 - t456 / 0.2e1;
t386 = t420 * t471 + t422 * t427 - t431 * t464;
t413 = -t427 * t466 + t470 * t471;
t376 = -t386 * t465 + t413 * t469;
t371 = qJD(4) * t376 + t437;
t435 = -t480 * t503 + V_base(5);
t424 = t476 * pkin(1) - pkin(10) * t504;
t494 = -t424 * t457 + V_base(5) * t507 + t500;
t387 = t421 * t427 + t423 * t471 + t431 * t468;
t352 = -t386 * t489 + t387 * t474 - t413 * t488;
t330 = qJD(5) * t352 + t371;
t425 = pkin(1) * t480 + pkin(10) * t505;
t493 = V_base(4) * t424 - t425 * V_base(5) + V_base(3);
t363 = qJD(4) * t369 + t435;
t378 = t414 * t421 + t415 * t468 - t423 * t504;
t335 = -t377 * t489 + t378 * t474 - t402 * t488;
t323 = qJD(5) * t335 + t363;
t399 = pkin(2) * t431 + qJ(3) * t413;
t492 = qJD(3) * t403 + t435 * t399 + t494;
t384 = t415 * pkin(2) + qJ(3) * t402;
t491 = qJD(3) * t413 + t436 * t384 + t493;
t490 = t457 * t425 + V_base(2) + (-pkin(9) - t507) * V_base(4);
t385 = pkin(2) * t417 + qJ(3) * t403;
t487 = qJD(3) * t402 + t437 * t385 + t490;
t341 = pkin(3) * t378 + pkin(11) * t369;
t356 = pkin(3) * t387 + pkin(11) * t376;
t486 = t435 * t356 + (-t341 - t384) * t437 + t492;
t342 = pkin(3) * t380 + pkin(11) * t370;
t485 = t436 * t341 + (-t342 - t385) * t435 + t491;
t484 = t437 * t342 + (-t356 - t399) * t436 + t487;
t426 = t495 - t497 / 0.2e1;
t429 = t496 - t498 / 0.2e1;
t478 = cos(qJ(4));
t336 = t377 * t426 + t378 * t478 + t402 * t429;
t315 = pkin(4) * t336 + pkin(12) * t335;
t353 = t386 * t426 + t387 * t478 + t413 * t429;
t329 = pkin(4) * t353 + pkin(12) * t352;
t483 = -t315 * t371 + t363 * t329 + t486;
t338 = t379 * t426 + t380 * t478 + t403 * t429;
t316 = pkin(4) * t338 + pkin(12) * t337;
t482 = t364 * t315 - t316 * t363 + t485;
t481 = t371 * t316 - t329 * t364 + t484;
t477 = cos(qJ(6));
t473 = sin(qJ(5));
t472 = sin(qJ(6));
t458 = Icges(2,4) * t480;
t445 = rSges(2,1) * t480 - t476 * rSges(2,2);
t444 = t476 * rSges(2,1) + rSges(2,2) * t480;
t443 = Icges(2,1) * t480 - t506;
t442 = Icges(2,1) * t476 + t458;
t441 = -Icges(2,2) * t476 + t458;
t440 = Icges(2,2) * t480 + t506;
t434 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t433 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t432 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t411 = V_base(5) * rSges(2,3) - t444 * t457 + t500;
t410 = t445 * t457 + V_base(2) + (-rSges(2,3) - pkin(9)) * V_base(4);
t409 = t444 * V_base(4) - t445 * V_base(5) + V_base(3);
t408 = rSges(3,1) * t431 + rSges(3,2) * t427 + rSges(3,3) * t471;
t407 = Icges(3,1) * t431 + Icges(3,4) * t427 + Icges(3,5) * t471;
t406 = Icges(3,4) * t431 + Icges(3,2) * t427 + Icges(3,6) * t471;
t405 = Icges(3,5) * t431 + Icges(3,6) * t427 + Icges(3,3) * t471;
t395 = rSges(3,1) * t417 + rSges(3,2) * t416 + rSges(3,3) * t505;
t394 = t415 * rSges(3,1) + t414 * rSges(3,2) - rSges(3,3) * t504;
t393 = Icges(3,1) * t417 + Icges(3,4) * t416 + Icges(3,5) * t505;
t392 = Icges(3,1) * t415 + Icges(3,4) * t414 - Icges(3,5) * t504;
t391 = Icges(3,4) * t417 + Icges(3,2) * t416 + Icges(3,6) * t505;
t390 = Icges(3,4) * t415 + Icges(3,2) * t414 - Icges(3,6) * t504;
t389 = Icges(3,5) * t417 + Icges(3,6) * t416 + Icges(3,3) * t505;
t388 = Icges(3,5) * t415 + Icges(3,6) * t414 - Icges(3,3) * t504;
t362 = -t394 * t437 + t408 * t435 + t494;
t361 = t395 * t437 - t408 * t436 + t490;
t360 = rSges(4,1) * t387 + rSges(4,2) * t386 + rSges(4,3) * t413;
t359 = Icges(4,1) * t387 + Icges(4,4) * t386 + Icges(4,5) * t413;
t358 = Icges(4,4) * t387 + Icges(4,2) * t386 + Icges(4,6) * t413;
t357 = Icges(4,5) * t387 + Icges(4,6) * t386 + Icges(4,3) * t413;
t354 = t394 * t436 - t395 * t435 + t493;
t350 = rSges(4,1) * t380 + rSges(4,2) * t379 + rSges(4,3) * t403;
t349 = rSges(4,1) * t378 + rSges(4,2) * t377 + rSges(4,3) * t402;
t348 = Icges(4,1) * t380 + Icges(4,4) * t379 + Icges(4,5) * t403;
t347 = Icges(4,1) * t378 + Icges(4,4) * t377 + Icges(4,5) * t402;
t346 = Icges(4,4) * t380 + Icges(4,2) * t379 + Icges(4,6) * t403;
t345 = Icges(4,4) * t378 + Icges(4,2) * t377 + Icges(4,6) * t402;
t344 = Icges(4,5) * t380 + Icges(4,6) * t379 + Icges(4,3) * t403;
t343 = Icges(4,5) * t378 + Icges(4,6) * t377 + Icges(4,3) * t402;
t332 = t353 * t508 + t376 * t473;
t331 = t353 * t473 - t376 * t508;
t328 = t338 * t508 + t370 * t473;
t327 = t338 * t473 - t370 * t508;
t326 = t336 * t508 + t369 * t473;
t325 = t336 * t473 - t369 * t508;
t322 = rSges(5,1) * t353 - rSges(5,2) * t352 + rSges(5,3) * t376;
t321 = Icges(5,1) * t353 - Icges(5,4) * t352 + Icges(5,5) * t376;
t320 = Icges(5,4) * t353 - Icges(5,2) * t352 + Icges(5,6) * t376;
t319 = Icges(5,5) * t353 - Icges(5,6) * t352 + Icges(5,3) * t376;
t318 = t332 * t477 + t352 * t472;
t317 = -t332 * t472 + t352 * t477;
t313 = pkin(5) * t332 + pkin(13) * t331;
t312 = qJD(6) * t331 + t330;
t310 = t360 * t435 + (-t349 - t384) * t437 + t492;
t309 = t350 * t437 + (-t360 - t399) * t436 + t487;
t308 = rSges(5,1) * t338 - rSges(5,2) * t337 + rSges(5,3) * t370;
t307 = rSges(5,1) * t336 - rSges(5,2) * t335 + rSges(5,3) * t369;
t306 = Icges(5,1) * t338 - Icges(5,4) * t337 + Icges(5,5) * t370;
t305 = Icges(5,1) * t336 - Icges(5,4) * t335 + Icges(5,5) * t369;
t304 = Icges(5,4) * t338 - Icges(5,2) * t337 + Icges(5,6) * t370;
t303 = Icges(5,4) * t336 - Icges(5,2) * t335 + Icges(5,6) * t369;
t302 = Icges(5,5) * t338 - Icges(5,6) * t337 + Icges(5,3) * t370;
t301 = Icges(5,5) * t336 - Icges(5,6) * t335 + Icges(5,3) * t369;
t300 = t328 * t477 + t337 * t472;
t299 = -t328 * t472 + t337 * t477;
t298 = t326 * t477 + t335 * t472;
t297 = -t326 * t472 + t335 * t477;
t295 = t349 * t436 + (-t350 - t385) * t435 + t491;
t294 = rSges(6,1) * t332 - rSges(6,2) * t331 + rSges(6,3) * t352;
t293 = Icges(6,1) * t332 - Icges(6,4) * t331 + Icges(6,5) * t352;
t292 = Icges(6,4) * t332 - Icges(6,2) * t331 + Icges(6,6) * t352;
t291 = Icges(6,5) * t332 - Icges(6,6) * t331 + Icges(6,3) * t352;
t290 = pkin(5) * t328 + pkin(13) * t327;
t289 = pkin(5) * t326 + pkin(13) * t325;
t288 = qJD(6) * t327 + t324;
t287 = qJD(6) * t325 + t323;
t286 = rSges(6,1) * t328 - rSges(6,2) * t327 + rSges(6,3) * t337;
t285 = rSges(6,1) * t326 - rSges(6,2) * t325 + rSges(6,3) * t335;
t284 = Icges(6,1) * t328 - Icges(6,4) * t327 + Icges(6,5) * t337;
t283 = Icges(6,1) * t326 - Icges(6,4) * t325 + Icges(6,5) * t335;
t282 = Icges(6,4) * t328 - Icges(6,2) * t327 + Icges(6,6) * t337;
t281 = Icges(6,4) * t326 - Icges(6,2) * t325 + Icges(6,6) * t335;
t280 = Icges(6,5) * t328 - Icges(6,6) * t327 + Icges(6,3) * t337;
t279 = Icges(6,5) * t326 - Icges(6,6) * t325 + Icges(6,3) * t335;
t278 = rSges(7,1) * t318 + rSges(7,2) * t317 + rSges(7,3) * t331;
t277 = Icges(7,1) * t318 + Icges(7,4) * t317 + Icges(7,5) * t331;
t276 = Icges(7,4) * t318 + Icges(7,2) * t317 + Icges(7,6) * t331;
t275 = Icges(7,5) * t318 + Icges(7,6) * t317 + Icges(7,3) * t331;
t274 = rSges(7,1) * t300 + rSges(7,2) * t299 + rSges(7,3) * t327;
t273 = rSges(7,1) * t298 + rSges(7,2) * t297 + rSges(7,3) * t325;
t272 = Icges(7,1) * t300 + Icges(7,4) * t299 + Icges(7,5) * t327;
t271 = Icges(7,1) * t298 + Icges(7,4) * t297 + Icges(7,5) * t325;
t270 = Icges(7,4) * t300 + Icges(7,2) * t299 + Icges(7,6) * t327;
t269 = Icges(7,4) * t298 + Icges(7,2) * t297 + Icges(7,6) * t325;
t268 = Icges(7,5) * t300 + Icges(7,6) * t299 + Icges(7,3) * t327;
t267 = Icges(7,5) * t298 + Icges(7,6) * t297 + Icges(7,3) * t325;
t266 = -t307 * t371 + t322 * t363 + t486;
t265 = t308 * t371 - t322 * t364 + t484;
t264 = t307 * t364 - t308 * t363 + t485;
t263 = -t285 * t330 + t294 * t323 + t483;
t262 = t286 * t330 - t294 * t324 + t481;
t261 = t285 * t324 - t286 * t323 + t482;
t260 = -t273 * t312 + t278 * t287 - t289 * t330 + t313 * t323 + t483;
t259 = t274 * t312 - t278 * t288 + t290 * t330 - t313 * t324 + t481;
t258 = t273 * t288 - t274 * t287 + t289 * t324 - t290 * t323 + t482;
t1 = Icges(2,3) * t457 ^ 2 / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t480 - Icges(2,6) * t476) * t457 + (Icges(1,1) / 0.2e1 - t476 * t441 / 0.2e1 + t443 * t509) * V_base(4) + ((t441 + t442) * t480 + (-t440 + t443) * t476) * V_base(5) / 0.2e1) * V_base(4) + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (Icges(2,5) * t476 + Icges(2,6) * t480) * t457 + (Icges(1,2) / 0.2e1 + t440 * t509 + t476 * t442 / 0.2e1) * V_base(5)) * V_base(5) + m(7) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(6) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(5) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(4) * (t295 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t287 * ((t268 * t325 + t270 * t297 + t272 * t298) * t288 + (t325 * t267 + t297 * t269 + t298 * t271) * t287 + (t275 * t325 + t276 * t297 + t277 * t298) * t312) / 0.2e1 + t288 * ((t327 * t268 + t299 * t270 + t300 * t272) * t288 + (t267 * t327 + t269 * t299 + t271 * t300) * t287 + (t275 * t327 + t276 * t299 + t277 * t300) * t312) / 0.2e1 + t312 * ((t268 * t331 + t270 * t317 + t272 * t318) * t288 + (t267 * t331 + t269 * t317 + t271 * t318) * t287 + (t331 * t275 + t317 * t276 + t318 * t277) * t312) / 0.2e1 + t323 * ((t280 * t335 - t282 * t325 + t284 * t326) * t324 + (t335 * t279 - t325 * t281 + t326 * t283) * t323 + (t291 * t335 - t292 * t325 + t293 * t326) * t330) / 0.2e1 + t324 * ((t337 * t280 - t327 * t282 + t328 * t284) * t324 + (t279 * t337 - t281 * t327 + t283 * t328) * t323 + (t291 * t337 - t292 * t327 + t293 * t328) * t330) / 0.2e1 + t330 * ((t280 * t352 - t282 * t331 + t284 * t332) * t324 + (t279 * t352 - t281 * t331 + t283 * t332) * t323 + (t352 * t291 - t331 * t292 + t332 * t293) * t330) / 0.2e1 + m(3) * (t354 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + t363 * ((t302 * t369 - t304 * t335 + t306 * t336) * t364 + (t369 * t301 - t335 * t303 + t336 * t305) * t363 + (t319 * t369 - t320 * t335 + t321 * t336) * t371) / 0.2e1 + t364 * ((t370 * t302 - t337 * t304 + t338 * t306) * t364 + (t301 * t370 - t303 * t337 + t305 * t338) * t363 + (t319 * t370 - t320 * t337 + t321 * t338) * t371) / 0.2e1 + t371 * ((t302 * t376 - t304 * t352 + t306 * t353) * t364 + (t301 * t376 - t303 * t352 + t305 * t353) * t363 + (t376 * t319 - t352 * t320 + t353 * t321) * t371) / 0.2e1 + m(2) * (t409 ^ 2 + t410 ^ 2 + t411 ^ 2) / 0.2e1 + m(1) * (t432 ^ 2 + t433 ^ 2 + t434 ^ 2) / 0.2e1 + ((t357 * t402 + t358 * t377 + t359 * t378 - t405 * t504 + t414 * t406 + t415 * t407) * t437 + (t344 * t402 + t346 * t377 + t348 * t378 - t389 * t504 + t414 * t391 + t415 * t393) * t436 + (t402 * t343 + t377 * t345 + t378 * t347 - t388 * t504 + t414 * t390 + t415 * t392) * t435) * t435 / 0.2e1 + ((t357 * t403 + t358 * t379 + t359 * t380 + t405 * t505 + t406 * t416 + t407 * t417) * t437 + (t403 * t344 + t379 * t346 + t380 * t348 + t389 * t505 + t416 * t391 + t417 * t393) * t436 + (t343 * t403 + t345 * t379 + t347 * t380 + t388 * t505 + t390 * t416 + t392 * t417) * t435) * t436 / 0.2e1 + ((t413 * t357 + t386 * t358 + t387 * t359 + t471 * t405 + t427 * t406 + t431 * t407) * t437 + (t344 * t413 + t346 * t386 + t348 * t387 + t389 * t471 + t391 * t427 + t393 * t431) * t436 + (t343 * t413 + t345 * t386 + t347 * t387 + t388 * t471 + t390 * t427 + t392 * t431) * t435) * t437 / 0.2e1;
T  = t1;
