% Calculate kinetic energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:09
% EndTime: 2024-09-27 21:45:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (2468->250), mult. (1753->383), div. (0->0), fcn. (1636->20), ass. (0->130)
t422 = sin(pkin(5));
t419 = pkin(10) + qJ(2);
t411 = sin(t419);
t412 = cos(t419);
t425 = sin(qJ(3));
t423 = cos(pkin(5));
t426 = cos(qJ(3));
t442 = t423 * t426;
t371 = -t411 * t425 + t412 * t442;
t443 = t423 * t425;
t372 = t411 * t426 + t412 * t443;
t373 = -t411 * t442 - t412 * t425;
t374 = -t411 * t443 + t412 * t426;
t444 = t412 * t422;
t445 = t411 * t422;
t432 = (Icges(4,5) * t372 + Icges(4,6) * t371 - Icges(4,3) * t444) * t412 - (Icges(4,5) * t374 + Icges(4,6) * t373 + Icges(4,3) * t445) * t411;
t447 = t432 * t422;
t427 = pkin(7) + pkin(8);
t446 = t426 * pkin(3);
t421 = qJ(3) + qJ(4);
t387 = pkin(3) * t443 - t422 * t427;
t407 = cos(qJ(4)) * pkin(4) + pkin(3);
t420 = pkin(9) + t427;
t436 = pkin(4) * sin(qJ(4)) * t426;
t441 = -t422 * t420 + (t407 * t425 + t436) * t423 - t387;
t416 = cos(t421);
t440 = pkin(4) * t416;
t438 = qJD(3) * t422;
t393 = t411 * t438;
t368 = qJD(4) * t445 + t393;
t439 = qJD(2) * (t411 * pkin(2) - pkin(7) * t444);
t400 = qJD(3) * t423 + qJD(2);
t437 = -qJD(3) - qJD(4);
t435 = t412 * t438;
t433 = pkin(7) * t422 + t387;
t343 = t446 * t411 + t433 * t412;
t344 = -t433 * t411 + t446 * t412;
t434 = t343 * t393 + t344 * t435 + qJD(1);
t390 = qJD(4) * t423 + t400;
t414 = pkin(5) - t421;
t413 = pkin(5) + t421;
t375 = qJD(2) * (t412 * pkin(2) + pkin(7) * t445);
t376 = t422 * t425 * pkin(3) + (-pkin(7) + t427) * t423;
t431 = t400 * t344 - t376 * t393 + t375;
t430 = -t400 * t343 - t376 * t435 - t439;
t429 = qJD(1) ^ 2;
t428 = qJD(2) ^ 2;
t417 = qJ(5) + t421;
t415 = sin(t421);
t406 = cos(t417);
t405 = sin(t417);
t404 = -qJ(5) + t414;
t403 = qJ(5) + t413;
t402 = cos(t413);
t401 = sin(t414);
t399 = cos(t403);
t398 = sin(t404);
t397 = cos(t414) / 0.2e1;
t396 = sin(t413) / 0.2e1;
t395 = cos(t404) / 0.2e1;
t394 = sin(t403) / 0.2e1;
t389 = t412 * rSges(3,1) - t411 * rSges(3,2);
t388 = t411 * rSges(3,1) + t412 * rSges(3,2);
t386 = t397 - t402 / 0.2e1;
t385 = t397 + t402 / 0.2e1;
t384 = t396 - t401 / 0.2e1;
t383 = t396 + t401 / 0.2e1;
t381 = t395 - t399 / 0.2e1;
t380 = t395 + t399 / 0.2e1;
t379 = t394 - t398 / 0.2e1;
t378 = t394 + t398 / 0.2e1;
t377 = qJD(5) * t423 + t390;
t370 = t423 * rSges(4,3) + (rSges(4,1) * t425 + rSges(4,2) * t426) * t422;
t369 = t437 * t444;
t367 = Icges(4,5) * t423 + (Icges(4,1) * t425 + Icges(4,4) * t426) * t422;
t366 = Icges(4,6) * t423 + (Icges(4,4) * t425 + Icges(4,2) * t426) * t422;
t365 = Icges(4,3) * t423 + (Icges(4,5) * t425 + Icges(4,6) * t426) * t422;
t363 = (-qJD(5) + t437) * t444;
t362 = qJD(5) * t445 + t368;
t361 = -t411 * t384 + t412 * t416;
t360 = -t411 * t385 - t412 * t415;
t359 = t412 * t384 + t411 * t416;
t358 = t412 * t385 - t411 * t415;
t357 = -t411 * t379 + t412 * t406;
t356 = -t411 * t380 - t412 * t405;
t355 = t412 * t379 + t411 * t406;
t354 = t412 * t380 - t411 * t405;
t353 = t386 * rSges(5,1) + t383 * rSges(5,2) + t423 * rSges(5,3);
t352 = Icges(5,1) * t386 + Icges(5,4) * t383 + Icges(5,5) * t423;
t351 = Icges(5,4) * t386 + Icges(5,2) * t383 + Icges(5,6) * t423;
t350 = Icges(5,5) * t386 + Icges(5,6) * t383 + Icges(5,3) * t423;
t349 = t381 * rSges(6,1) + t378 * rSges(6,2) + t423 * rSges(6,3);
t348 = (t420 - t427) * t423 + (t436 + (-pkin(3) + t407) * t425) * t422;
t347 = Icges(6,1) * t381 + Icges(6,4) * t378 + Icges(6,5) * t423;
t346 = Icges(6,4) * t381 + Icges(6,2) * t378 + Icges(6,6) * t423;
t345 = Icges(6,5) * t381 + Icges(6,6) * t378 + Icges(6,3) * t423;
t342 = t374 * rSges(4,1) + t373 * rSges(4,2) + rSges(4,3) * t445;
t341 = t372 * rSges(4,1) + t371 * rSges(4,2) - rSges(4,3) * t444;
t340 = Icges(4,1) * t374 + Icges(4,4) * t373 + Icges(4,5) * t445;
t339 = Icges(4,1) * t372 + Icges(4,4) * t371 - Icges(4,5) * t444;
t338 = Icges(4,4) * t374 + Icges(4,2) * t373 + Icges(4,6) * t445;
t337 = Icges(4,4) * t372 + Icges(4,2) * t371 - Icges(4,6) * t444;
t331 = t361 * rSges(5,1) + t360 * rSges(5,2) + rSges(5,3) * t445;
t330 = t359 * rSges(5,1) + t358 * rSges(5,2) - rSges(5,3) * t444;
t329 = Icges(5,1) * t361 + Icges(5,4) * t360 + Icges(5,5) * t445;
t328 = Icges(5,1) * t359 + Icges(5,4) * t358 - Icges(5,5) * t444;
t327 = Icges(5,4) * t361 + Icges(5,2) * t360 + Icges(5,6) * t445;
t326 = Icges(5,4) * t359 + Icges(5,2) * t358 - Icges(5,6) * t444;
t325 = Icges(5,5) * t361 + Icges(5,6) * t360 + Icges(5,3) * t445;
t324 = Icges(5,5) * t359 + Icges(5,6) * t358 - Icges(5,3) * t444;
t323 = t357 * rSges(6,1) + t356 * rSges(6,2) + rSges(6,3) * t445;
t322 = t355 * rSges(6,1) + t354 * rSges(6,2) - rSges(6,3) * t444;
t321 = Icges(6,1) * t357 + Icges(6,4) * t356 + Icges(6,5) * t445;
t320 = Icges(6,1) * t355 + Icges(6,4) * t354 - Icges(6,5) * t444;
t319 = Icges(6,4) * t357 + Icges(6,2) * t356 + Icges(6,6) * t445;
t318 = Icges(6,4) * t355 + Icges(6,2) * t354 - Icges(6,6) * t444;
t317 = Icges(6,5) * t357 + Icges(6,6) * t356 + Icges(6,3) * t445;
t316 = Icges(6,5) * t355 + Icges(6,6) * t354 - Icges(6,3) * t444;
t315 = -t441 * t411 + t440 * t412;
t314 = t440 * t411 + t441 * t412;
t313 = t400 * t342 - t370 * t393 + t375;
t312 = -t400 * t341 - t370 * t435 - t439;
t311 = qJD(1) + (t341 * t411 + t342 * t412) * t438;
t310 = t390 * t331 - t368 * t353 + t431;
t309 = -t390 * t330 + t369 * t353 + t430;
t308 = t368 * t330 - t369 * t331 + t434;
t307 = t390 * t315 + t377 * t323 - t368 * t348 - t362 * t349 + t431;
t306 = -t314 * t390 - t322 * t377 + t348 * t369 + t349 * t363 + t430;
t305 = t314 * t368 - t315 * t369 + t322 * t362 - t323 * t363 + t434;
t1 = m(2) * t429 / 0.2e1 + m(3) * (t429 + (t388 ^ 2 + t389 ^ 2) * t428) / 0.2e1 + t428 * Icges(3,3) / 0.2e1 + m(4) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + ((t365 * t445 + t366 * t373 + t367 * t374) * t400 + (-(t337 * t373 + t339 * t374) * t412 + (t373 * t338 + t374 * t340 - t447) * t411) * t438) * t393 / 0.2e1 - ((-t365 * t444 + t366 * t371 + t367 * t372) * t400 + ((t338 * t371 + t340 * t372) * t411 + (-t371 * t337 - t372 * t339 + t447) * t412) * t438) * t435 / 0.2e1 + t400 * ((t423 * t365 + (t366 * t426 + t367 * t425) * t422) * t400 + (((t338 * t426 + t340 * t425) * t411 - (t337 * t426 + t339 * t425) * t412) * t422 - t432 * t423) * t438) / 0.2e1 + m(5) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t368 * ((t325 * t445 + t360 * t327 + t361 * t329) * t368 + (t324 * t445 + t326 * t360 + t328 * t361) * t369 + (t350 * t445 + t351 * t360 + t352 * t361) * t390) / 0.2e1 + t369 * ((-t325 * t444 + t327 * t358 + t329 * t359) * t368 + (-t324 * t444 + t358 * t326 + t359 * t328) * t369 + (-t350 * t444 + t351 * t358 + t352 * t359) * t390) / 0.2e1 + t390 * ((t325 * t423 + t327 * t383 + t329 * t386) * t368 + (t324 * t423 + t326 * t383 + t328 * t386) * t369 + (t423 * t350 + t383 * t351 + t386 * t352) * t390) / 0.2e1 + m(6) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + t362 * ((t317 * t445 + t356 * t319 + t357 * t321) * t362 + (t316 * t445 + t318 * t356 + t320 * t357) * t363 + (t345 * t445 + t346 * t356 + t347 * t357) * t377) / 0.2e1 + t363 * ((-t317 * t444 + t319 * t354 + t321 * t355) * t362 + (-t316 * t444 + t354 * t318 + t355 * t320) * t363 + (-t345 * t444 + t346 * t354 + t347 * t355) * t377) / 0.2e1 + t377 * ((t317 * t423 + t319 * t378 + t321 * t381) * t362 + (t316 * t423 + t318 * t378 + t320 * t381) * t363 + (t423 * t345 + t378 * t346 + t381 * t347) * t377) / 0.2e1;
T = t1;
