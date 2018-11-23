% Calculate kinetic energy for
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:08
% EndTime: 2018-11-23 10:29:13
% DurationCPUTime: 4.49s
% Computational Cost: add. (34471->469), mult. (35026->685), div. (0->0), fcn. (35064->30), ass. (0->220)
t484 = cos(qJ(1));
t513 = t484 / 0.2e1;
t512 = cos(qJ(5));
t473 = cos(pkin(6));
t511 = pkin(10) * t473;
t479 = sin(qJ(1));
t510 = Icges(2,4) * t479;
t470 = sin(pkin(6));
t509 = t470 * t479;
t508 = t470 * t484;
t507 = qJD(2) * t470;
t506 = pkin(8) - qJ(4);
t505 = pkin(8) + qJ(4);
t504 = V_base(5) * pkin(9) + V_base(1);
t440 = t479 * t507 + V_base(4);
t461 = V_base(6) + qJD(1);
t502 = cos(t505);
t501 = sin(t506);
t467 = pkin(6) - qJ(2);
t455 = cos(t467) / 0.2e1;
t466 = pkin(6) + qJ(2);
t460 = cos(t466);
t434 = t455 + t460 / 0.2e1;
t478 = sin(qJ(2));
t420 = -t479 * t434 - t478 * t484;
t469 = sin(pkin(7));
t472 = cos(pkin(7));
t406 = -t420 * t469 + t472 * t509;
t399 = qJD(3) * t406 + t440;
t441 = qJD(2) * t473 + t461;
t500 = cos(t506) / 0.2e1;
t499 = sin(t505) / 0.2e1;
t453 = sin(t466) / 0.2e1;
t458 = sin(t467);
t430 = t453 - t458 / 0.2e1;
t483 = cos(qJ(2));
t421 = -t479 * t430 + t483 * t484;
t464 = pkin(7) + qJ(3);
t452 = sin(t464) / 0.2e1;
t465 = pkin(7) - qJ(3);
t457 = sin(t465);
t427 = t452 + t457 / 0.2e1;
t454 = cos(t465) / 0.2e1;
t459 = cos(t464);
t432 = t454 + t459 / 0.2e1;
t477 = sin(qJ(3));
t380 = t420 * t432 - t421 * t477 + t427 * t509;
t468 = sin(pkin(8));
t471 = cos(pkin(8));
t372 = -t380 * t468 + t406 * t471;
t354 = qJD(4) * t372 + t399;
t418 = t434 * t484 - t479 * t478;
t419 = t430 * t484 + t479 * t483;
t378 = t418 * t432 - t419 * t477 - t427 * t508;
t405 = -t418 * t469 - t472 * t508;
t371 = -t378 * t468 + t405 * t471;
t429 = t453 + t458 / 0.2e1;
t435 = t455 - t460 / 0.2e1;
t387 = t427 * t473 + t429 * t432 - t435 * t477;
t417 = -t429 * t469 + t472 * t473;
t377 = -t387 * t468 + t417 * t471;
t408 = qJD(3) * t417 + t441;
t439 = -t484 * t507 + V_base(5);
t428 = t452 - t457 / 0.2e1;
t433 = t454 - t459 / 0.2e1;
t482 = cos(qJ(3));
t381 = t420 * t428 + t421 * t482 + t433 * t509;
t476 = sin(qJ(4));
t494 = t499 + t501 / 0.2e1;
t495 = t500 + t502 / 0.2e1;
t340 = -t380 * t495 + t381 * t476 - t406 * t494;
t325 = qJD(5) * t340 + t354;
t424 = t479 * pkin(1) - pkin(10) * t508;
t498 = -t424 * t461 + V_base(5) * t511 + t504;
t370 = qJD(4) * t377 + t408;
t425 = pkin(1) * t484 + pkin(10) * t509;
t497 = V_base(4) * t424 - t425 * V_base(5) + V_base(3);
t398 = qJD(3) * t405 + t439;
t388 = t428 * t429 + t433 * t473 + t435 * t482;
t357 = -t387 * t495 + t388 * t476 - t417 * t494;
t331 = qJD(5) * t357 + t370;
t353 = qJD(4) * t371 + t398;
t496 = t461 * t425 + V_base(2) + (-pkin(9) - t511) * V_base(4);
t379 = t418 * t428 + t419 * t482 - t433 * t508;
t338 = -t378 * t495 + t379 * t476 - t405 * t494;
t324 = qJD(5) * t338 + t353;
t385 = t419 * pkin(2) + pkin(11) * t405;
t402 = pkin(2) * t435 + pkin(11) * t417;
t493 = -t385 * t441 + t439 * t402 + t498;
t386 = pkin(2) * t421 + pkin(11) * t406;
t492 = t440 * t385 - t386 * t439 + t497;
t491 = t441 * t386 - t402 * t440 + t496;
t342 = pkin(3) * t379 + pkin(12) * t371;
t359 = pkin(3) * t388 + pkin(12) * t377;
t490 = -t342 * t408 + t398 * t359 + t493;
t343 = pkin(3) * t381 + pkin(12) * t372;
t489 = t399 * t342 - t343 * t398 + t492;
t488 = t408 * t343 - t359 * t399 + t491;
t426 = t499 - t501 / 0.2e1;
t431 = t500 - t502 / 0.2e1;
t481 = cos(qJ(4));
t339 = t378 * t426 + t379 * t481 + t405 * t431;
t316 = pkin(4) * t339 + pkin(13) * t338;
t358 = t387 * t426 + t388 * t481 + t417 * t431;
t330 = pkin(4) * t358 + pkin(13) * t357;
t487 = -t316 * t370 + t353 * t330 + t490;
t341 = t380 * t426 + t381 * t481 + t406 * t431;
t317 = pkin(4) * t341 + pkin(13) * t340;
t486 = t354 * t316 - t317 * t353 + t489;
t485 = t370 * t317 - t330 * t354 + t488;
t480 = cos(qJ(6));
t475 = sin(qJ(5));
t474 = sin(qJ(6));
t462 = Icges(2,4) * t484;
t449 = rSges(2,1) * t484 - t479 * rSges(2,2);
t448 = t479 * rSges(2,1) + rSges(2,2) * t484;
t447 = Icges(2,1) * t484 - t510;
t446 = Icges(2,1) * t479 + t462;
t445 = -Icges(2,2) * t479 + t462;
t444 = Icges(2,2) * t484 + t510;
t438 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t437 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t436 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t415 = V_base(5) * rSges(2,3) - t448 * t461 + t504;
t414 = t449 * t461 + V_base(2) + (-rSges(2,3) - pkin(9)) * V_base(4);
t413 = t448 * V_base(4) - t449 * V_base(5) + V_base(3);
t412 = rSges(3,1) * t435 + rSges(3,2) * t429 + rSges(3,3) * t473;
t411 = Icges(3,1) * t435 + Icges(3,4) * t429 + Icges(3,5) * t473;
t410 = Icges(3,4) * t435 + Icges(3,2) * t429 + Icges(3,6) * t473;
t409 = Icges(3,5) * t435 + Icges(3,6) * t429 + Icges(3,3) * t473;
t396 = rSges(3,1) * t421 + rSges(3,2) * t420 + rSges(3,3) * t509;
t395 = t419 * rSges(3,1) + t418 * rSges(3,2) - rSges(3,3) * t508;
t394 = Icges(3,1) * t421 + Icges(3,4) * t420 + Icges(3,5) * t509;
t393 = Icges(3,1) * t419 + Icges(3,4) * t418 - Icges(3,5) * t508;
t392 = Icges(3,4) * t421 + Icges(3,2) * t420 + Icges(3,6) * t509;
t391 = Icges(3,4) * t419 + Icges(3,2) * t418 - Icges(3,6) * t508;
t390 = Icges(3,5) * t421 + Icges(3,6) * t420 + Icges(3,3) * t509;
t389 = Icges(3,5) * t419 + Icges(3,6) * t418 - Icges(3,3) * t508;
t365 = rSges(4,1) * t388 + rSges(4,2) * t387 + rSges(4,3) * t417;
t364 = Icges(4,1) * t388 + Icges(4,4) * t387 + Icges(4,5) * t417;
t363 = Icges(4,4) * t388 + Icges(4,2) * t387 + Icges(4,6) * t417;
t362 = Icges(4,5) * t388 + Icges(4,6) * t387 + Icges(4,3) * t417;
t361 = -t395 * t441 + t412 * t439 + t498;
t360 = t396 * t441 - t412 * t440 + t496;
t355 = t395 * t440 - t396 * t439 + t497;
t351 = rSges(4,1) * t381 + rSges(4,2) * t380 + rSges(4,3) * t406;
t350 = rSges(4,1) * t379 + rSges(4,2) * t378 + rSges(4,3) * t405;
t349 = Icges(4,1) * t381 + Icges(4,4) * t380 + Icges(4,5) * t406;
t348 = Icges(4,1) * t379 + Icges(4,4) * t378 + Icges(4,5) * t405;
t347 = Icges(4,4) * t381 + Icges(4,2) * t380 + Icges(4,6) * t406;
t346 = Icges(4,4) * t379 + Icges(4,2) * t378 + Icges(4,6) * t405;
t345 = Icges(4,5) * t381 + Icges(4,6) * t380 + Icges(4,3) * t406;
t344 = Icges(4,5) * t379 + Icges(4,6) * t378 + Icges(4,3) * t405;
t333 = t358 * t512 + t377 * t475;
t332 = t358 * t475 - t377 * t512;
t329 = t341 * t512 + t372 * t475;
t328 = t341 * t475 - t372 * t512;
t327 = t339 * t512 + t371 * t475;
t326 = t339 * t475 - t371 * t512;
t323 = rSges(5,1) * t358 - rSges(5,2) * t357 + rSges(5,3) * t377;
t322 = Icges(5,1) * t358 - Icges(5,4) * t357 + Icges(5,5) * t377;
t321 = Icges(5,4) * t358 - Icges(5,2) * t357 + Icges(5,6) * t377;
t320 = Icges(5,5) * t358 - Icges(5,6) * t357 + Icges(5,3) * t377;
t319 = t333 * t480 + t357 * t474;
t318 = -t333 * t474 + t357 * t480;
t314 = pkin(5) * t333 + pkin(14) * t332;
t313 = qJD(6) * t332 + t331;
t312 = rSges(5,1) * t341 - rSges(5,2) * t340 + rSges(5,3) * t372;
t311 = rSges(5,1) * t339 - rSges(5,2) * t338 + rSges(5,3) * t371;
t309 = Icges(5,1) * t341 - Icges(5,4) * t340 + Icges(5,5) * t372;
t308 = Icges(5,1) * t339 - Icges(5,4) * t338 + Icges(5,5) * t371;
t307 = Icges(5,4) * t341 - Icges(5,2) * t340 + Icges(5,6) * t372;
t306 = Icges(5,4) * t339 - Icges(5,2) * t338 + Icges(5,6) * t371;
t305 = Icges(5,5) * t341 - Icges(5,6) * t340 + Icges(5,3) * t372;
t304 = Icges(5,5) * t339 - Icges(5,6) * t338 + Icges(5,3) * t371;
t303 = -t350 * t408 + t365 * t398 + t493;
t302 = t351 * t408 - t365 * t399 + t491;
t301 = t329 * t480 + t340 * t474;
t300 = -t329 * t474 + t340 * t480;
t299 = t327 * t480 + t338 * t474;
t298 = -t327 * t474 + t338 * t480;
t296 = t350 * t399 - t351 * t398 + t492;
t295 = rSges(6,1) * t333 - rSges(6,2) * t332 + rSges(6,3) * t357;
t294 = Icges(6,1) * t333 - Icges(6,4) * t332 + Icges(6,5) * t357;
t293 = Icges(6,4) * t333 - Icges(6,2) * t332 + Icges(6,6) * t357;
t292 = Icges(6,5) * t333 - Icges(6,6) * t332 + Icges(6,3) * t357;
t291 = pkin(5) * t329 + pkin(14) * t328;
t290 = pkin(5) * t327 + pkin(14) * t326;
t289 = qJD(6) * t328 + t325;
t288 = qJD(6) * t326 + t324;
t287 = rSges(6,1) * t329 - rSges(6,2) * t328 + rSges(6,3) * t340;
t286 = rSges(6,1) * t327 - rSges(6,2) * t326 + rSges(6,3) * t338;
t285 = Icges(6,1) * t329 - Icges(6,4) * t328 + Icges(6,5) * t340;
t284 = Icges(6,1) * t327 - Icges(6,4) * t326 + Icges(6,5) * t338;
t283 = Icges(6,4) * t329 - Icges(6,2) * t328 + Icges(6,6) * t340;
t282 = Icges(6,4) * t327 - Icges(6,2) * t326 + Icges(6,6) * t338;
t281 = Icges(6,5) * t329 - Icges(6,6) * t328 + Icges(6,3) * t340;
t280 = Icges(6,5) * t327 - Icges(6,6) * t326 + Icges(6,3) * t338;
t279 = rSges(7,1) * t319 + rSges(7,2) * t318 + rSges(7,3) * t332;
t278 = Icges(7,1) * t319 + Icges(7,4) * t318 + Icges(7,5) * t332;
t277 = Icges(7,4) * t319 + Icges(7,2) * t318 + Icges(7,6) * t332;
t276 = Icges(7,5) * t319 + Icges(7,6) * t318 + Icges(7,3) * t332;
t275 = rSges(7,1) * t301 + rSges(7,2) * t300 + rSges(7,3) * t328;
t274 = rSges(7,1) * t299 + rSges(7,2) * t298 + rSges(7,3) * t326;
t273 = Icges(7,1) * t301 + Icges(7,4) * t300 + Icges(7,5) * t328;
t272 = Icges(7,1) * t299 + Icges(7,4) * t298 + Icges(7,5) * t326;
t271 = Icges(7,4) * t301 + Icges(7,2) * t300 + Icges(7,6) * t328;
t270 = Icges(7,4) * t299 + Icges(7,2) * t298 + Icges(7,6) * t326;
t269 = Icges(7,5) * t301 + Icges(7,6) * t300 + Icges(7,3) * t328;
t268 = Icges(7,5) * t299 + Icges(7,6) * t298 + Icges(7,3) * t326;
t267 = -t311 * t370 + t323 * t353 + t490;
t266 = t312 * t370 - t323 * t354 + t488;
t265 = t311 * t354 - t312 * t353 + t489;
t264 = -t286 * t331 + t295 * t324 + t487;
t263 = t287 * t331 - t295 * t325 + t485;
t262 = t286 * t325 - t287 * t324 + t486;
t261 = -t274 * t313 + t279 * t288 - t290 * t331 + t314 * t324 + t487;
t260 = t275 * t313 - t279 * t289 + t291 * t331 - t314 * t325 + t485;
t259 = t274 * t289 - t275 * t288 + t290 * t325 - t291 * t324 + t486;
t1 = Icges(2,3) * t461 ^ 2 / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t484 - Icges(2,6) * t479) * t461 + (Icges(1,1) / 0.2e1 - t479 * t445 / 0.2e1 + t447 * t513) * V_base(4) + ((t445 + t446) * t484 + (-t444 + t447) * t479) * V_base(5) / 0.2e1) * V_base(4) + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + t440 * ((t390 * t509 + t392 * t420 + t394 * t421) * t440 + (t389 * t509 + t391 * t420 + t393 * t421) * t439 + (t409 * t509 + t410 * t420 + t411 * t421) * t441) / 0.2e1 + t439 * ((-t390 * t508 + t418 * t392 + t419 * t394) * t440 + (-t389 * t508 + t418 * t391 + t419 * t393) * t439 + (-t409 * t508 + t418 * t410 + t419 * t411) * t441) / 0.2e1 + (Icges(1,6) * V_base(6) + (Icges(2,5) * t479 + Icges(2,6) * t484) * t461 + (Icges(1,2) / 0.2e1 + t444 * t513 + t479 * t446 / 0.2e1) * V_base(5)) * V_base(5) + m(7) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(4) * (t296 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + t288 * ((t269 * t326 + t271 * t298 + t273 * t299) * t289 + (t326 * t268 + t298 * t270 + t299 * t272) * t288 + (t276 * t326 + t277 * t298 + t278 * t299) * t313) / 0.2e1 + t289 * ((t328 * t269 + t300 * t271 + t301 * t273) * t289 + (t268 * t328 + t270 * t300 + t272 * t301) * t288 + (t276 * t328 + t277 * t300 + t278 * t301) * t313) / 0.2e1 + t313 * ((t269 * t332 + t271 * t318 + t273 * t319) * t289 + (t268 * t332 + t270 * t318 + t272 * t319) * t288 + (t332 * t276 + t318 * t277 + t319 * t278) * t313) / 0.2e1 + t324 * ((t281 * t338 - t283 * t326 + t285 * t327) * t325 + (t338 * t280 - t326 * t282 + t327 * t284) * t324 + (t292 * t338 - t293 * t326 + t294 * t327) * t331) / 0.2e1 + t325 * ((t340 * t281 - t328 * t283 + t329 * t285) * t325 + (t280 * t340 - t282 * t328 + t284 * t329) * t324 + (t292 * t340 - t293 * t328 + t294 * t329) * t331) / 0.2e1 + t331 * ((t281 * t357 - t283 * t332 + t285 * t333) * t325 + (t280 * t357 - t282 * t332 + t284 * t333) * t324 + (t357 * t292 - t332 * t293 + t333 * t294) * t331) / 0.2e1 + m(3) * (t355 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + t353 * ((t305 * t371 - t307 * t338 + t309 * t339) * t354 + (t304 * t371 - t306 * t338 + t308 * t339) * t353 + (t320 * t371 - t321 * t338 + t322 * t339) * t370) / 0.2e1 + t354 * ((t305 * t372 - t307 * t340 + t341 * t309) * t354 + (t304 * t372 - t306 * t340 + t308 * t341) * t353 + (t320 * t372 - t321 * t340 + t322 * t341) * t370) / 0.2e1 + t370 * ((t305 * t377 - t307 * t357 + t309 * t358) * t354 + (t304 * t377 - t306 * t357 + t308 * t358) * t353 + (t320 * t377 - t321 * t357 + t322 * t358) * t370) / 0.2e1 + t398 * ((t345 * t405 + t347 * t378 + t349 * t379) * t399 + (t344 * t405 + t346 * t378 + t348 * t379) * t398 + (t362 * t405 + t363 * t378 + t364 * t379) * t408) / 0.2e1 + t399 * ((t345 * t406 + t347 * t380 + t349 * t381) * t399 + (t344 * t406 + t346 * t380 + t348 * t381) * t398 + (t362 * t406 + t363 * t380 + t364 * t381) * t408) / 0.2e1 + m(2) * (t413 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + t408 * ((t345 * t417 + t347 * t387 + t349 * t388) * t399 + (t344 * t417 + t346 * t387 + t348 * t388) * t398 + (t362 * t417 + t363 * t387 + t364 * t388) * t408) / 0.2e1 + m(1) * (t436 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + t441 * ((t390 * t473 + t392 * t429 + t394 * t435) * t440 + (t389 * t473 + t391 * t429 + t393 * t435) * t439 + (t409 * t473 + t410 * t429 + t411 * t435) * t441) / 0.2e1;
T  = t1;
