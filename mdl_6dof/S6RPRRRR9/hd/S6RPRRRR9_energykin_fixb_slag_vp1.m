% Calculate kinetic energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:01
% EndTime: 2019-03-09 07:23:03
% DurationCPUTime: 1.98s
% Computational Cost: add. (1342->298), mult. (1968->482), div. (0->0), fcn. (1958->10), ass. (0->145)
t435 = sin(qJ(3));
t438 = cos(qJ(3));
t433 = qJ(4) + qJ(5);
t429 = cos(t433);
t462 = pkin(5) * t429;
t478 = -pkin(10) * t438 + t435 * t462;
t437 = cos(qJ(4));
t474 = t437 * pkin(4);
t477 = -pkin(9) * t438 + t435 * t474;
t472 = Icges(4,4) * t435;
t471 = Icges(4,4) * t438;
t434 = sin(qJ(4));
t436 = sin(qJ(1));
t470 = t434 * t436;
t439 = cos(qJ(1));
t469 = t434 * t439;
t468 = t435 * t436;
t467 = t435 * t439;
t466 = t436 * t437;
t465 = t436 * t438;
t464 = t437 * t439;
t463 = t438 * t439;
t406 = qJD(1) * (pkin(1) * t439 + qJ(2) * t436);
t461 = pkin(7) * qJD(1) * t439 + t406;
t425 = qJD(3) * t436;
t459 = qJD(4) * t438;
t402 = t439 * t459 + t425;
t426 = qJD(3) * t439;
t418 = qJD(4) * t435 + qJD(1);
t458 = -qJD(4) - qJD(5);
t374 = qJD(5) * t463 + t402;
t405 = qJD(5) * t435 + t418;
t428 = sin(t433);
t457 = pkin(5) * t428;
t411 = pkin(1) * t436 - qJ(2) * t439;
t456 = -pkin(7) * t436 - t411;
t455 = pkin(3) * t435 - pkin(8) * t438;
t399 = t455 * t436;
t400 = t455 * t439;
t454 = -t399 * t425 - t400 * t426;
t453 = rSges(4,1) * t435 + rSges(4,2) * t438;
t452 = Icges(4,1) * t435 + t471;
t451 = Icges(4,2) * t438 + t472;
t450 = Icges(4,5) * t435 + Icges(4,6) * t438;
t380 = Icges(4,6) * t439 + t436 * t451;
t383 = Icges(4,5) * t439 + t436 * t452;
t449 = -t380 * t438 - t383 * t435;
t381 = Icges(4,6) * t436 - t439 * t451;
t384 = Icges(4,5) * t436 - t439 * t452;
t448 = t381 * t438 + t384 * t435;
t409 = -Icges(4,2) * t435 + t471;
t410 = Icges(4,1) * t438 - t472;
t447 = t409 * t438 + t410 * t435;
t353 = pkin(4) * t469 + t436 * t477;
t354 = pkin(4) * t470 - t439 * t477;
t403 = -t436 * t459 + t426;
t446 = -t353 * t402 + t354 * t403 + t454;
t415 = pkin(3) * t438 + pkin(8) * t435;
t427 = qJD(2) * t436;
t445 = t415 * t425 + t427 + (t400 + t456) * qJD(1);
t444 = qJD(1) * t399 + (-qJD(3) * t415 - qJD(2)) * t439 + t461;
t365 = pkin(9) * t435 + t438 * t474;
t443 = -t354 * t418 + t365 * t402 + t445;
t442 = t353 * t418 - t365 * t403 + t444;
t430 = qJ(6) + t433;
t420 = cos(t430);
t419 = sin(t430);
t414 = rSges(2,1) * t439 - rSges(2,2) * t436;
t413 = rSges(4,1) * t438 - rSges(4,2) * t435;
t412 = rSges(2,1) * t436 + rSges(2,2) * t439;
t408 = Icges(4,5) * t438 - Icges(4,6) * t435;
t398 = -t435 * t464 + t470;
t397 = t434 * t467 + t466;
t396 = t435 * t466 + t469;
t395 = -t434 * t468 + t464;
t394 = qJD(6) * t435 + t405;
t391 = t428 * t436 - t429 * t467;
t390 = t428 * t467 + t429 * t436;
t389 = t428 * t439 + t429 * t468;
t388 = -t428 * t468 + t429 * t439;
t387 = rSges(4,3) * t436 - t439 * t453;
t386 = rSges(5,3) * t435 + (rSges(5,1) * t437 - rSges(5,2) * t434) * t438;
t385 = rSges(4,3) * t439 + t436 * t453;
t382 = Icges(5,5) * t435 + (Icges(5,1) * t437 - Icges(5,4) * t434) * t438;
t379 = Icges(5,6) * t435 + (Icges(5,4) * t437 - Icges(5,2) * t434) * t438;
t378 = Icges(4,3) * t436 - t439 * t450;
t377 = Icges(4,3) * t439 + t436 * t450;
t376 = Icges(5,3) * t435 + (Icges(5,5) * t437 - Icges(5,6) * t434) * t438;
t375 = t458 * t465 + t426;
t373 = t419 * t436 - t420 * t467;
t372 = t419 * t467 + t420 * t436;
t371 = t419 * t439 + t420 * t468;
t370 = -t419 * t468 + t420 * t439;
t369 = rSges(6,3) * t435 + (rSges(6,1) * t429 - rSges(6,2) * t428) * t438;
t368 = Icges(6,5) * t435 + (Icges(6,1) * t429 - Icges(6,4) * t428) * t438;
t367 = Icges(6,6) * t435 + (Icges(6,4) * t429 - Icges(6,2) * t428) * t438;
t366 = Icges(6,3) * t435 + (Icges(6,5) * t429 - Icges(6,6) * t428) * t438;
t364 = rSges(7,3) * t435 + (rSges(7,1) * t420 - rSges(7,2) * t419) * t438;
t363 = Icges(7,5) * t435 + (Icges(7,1) * t420 - Icges(7,4) * t419) * t438;
t362 = Icges(7,6) * t435 + (Icges(7,4) * t420 - Icges(7,2) * t419) * t438;
t361 = Icges(7,3) * t435 + (Icges(7,5) * t420 - Icges(7,6) * t419) * t438;
t360 = t426 + (-qJD(6) + t458) * t465;
t359 = qJD(6) * t463 + t374;
t358 = t406 - qJD(2) * t439 + qJD(1) * (-rSges(3,2) * t439 + rSges(3,3) * t436);
t357 = t427 + (rSges(3,2) * t436 + rSges(3,3) * t439 - t411) * qJD(1);
t355 = pkin(10) * t435 + t438 * t462;
t352 = rSges(5,1) * t398 + rSges(5,2) * t397 + rSges(5,3) * t463;
t351 = rSges(5,1) * t396 + rSges(5,2) * t395 - rSges(5,3) * t465;
t350 = Icges(5,1) * t398 + Icges(5,4) * t397 + Icges(5,5) * t463;
t349 = Icges(5,1) * t396 + Icges(5,4) * t395 - Icges(5,5) * t465;
t348 = Icges(5,4) * t398 + Icges(5,2) * t397 + Icges(5,6) * t463;
t347 = Icges(5,4) * t396 + Icges(5,2) * t395 - Icges(5,6) * t465;
t346 = Icges(5,5) * t398 + Icges(5,6) * t397 + Icges(5,3) * t463;
t345 = Icges(5,5) * t396 + Icges(5,6) * t395 - Icges(5,3) * t465;
t344 = (-t385 * t436 + t387 * t439) * qJD(3);
t342 = rSges(6,1) * t391 + rSges(6,2) * t390 + rSges(6,3) * t463;
t341 = rSges(6,1) * t389 + rSges(6,2) * t388 - rSges(6,3) * t465;
t340 = Icges(6,1) * t391 + Icges(6,4) * t390 + Icges(6,5) * t463;
t339 = Icges(6,1) * t389 + Icges(6,4) * t388 - Icges(6,5) * t465;
t338 = Icges(6,4) * t391 + Icges(6,2) * t390 + Icges(6,6) * t463;
t337 = Icges(6,4) * t389 + Icges(6,2) * t388 - Icges(6,6) * t465;
t336 = Icges(6,5) * t391 + Icges(6,6) * t390 + Icges(6,3) * t463;
t335 = Icges(6,5) * t389 + Icges(6,6) * t388 - Icges(6,3) * t465;
t334 = rSges(7,1) * t373 + rSges(7,2) * t372 + rSges(7,3) * t463;
t333 = rSges(7,1) * t371 + rSges(7,2) * t370 - rSges(7,3) * t465;
t331 = Icges(7,1) * t373 + Icges(7,4) * t372 + Icges(7,5) * t463;
t330 = Icges(7,1) * t371 + Icges(7,4) * t370 - Icges(7,5) * t465;
t329 = Icges(7,4) * t373 + Icges(7,2) * t372 + Icges(7,6) * t463;
t328 = Icges(7,4) * t371 + Icges(7,2) * t370 - Icges(7,6) * t465;
t327 = Icges(7,5) * t373 + Icges(7,6) * t372 + Icges(7,3) * t463;
t326 = Icges(7,5) * t371 + Icges(7,6) * t370 - Icges(7,3) * t465;
t325 = qJD(1) * t385 + (-qJD(3) * t413 - qJD(2)) * t439 + t461;
t324 = t413 * t425 + t427 + (-t387 + t456) * qJD(1);
t323 = t436 * t457 - t439 * t478;
t322 = t436 * t478 + t439 * t457;
t321 = t351 * t418 - t386 * t403 + t444;
t320 = -t352 * t418 + t386 * t402 + t445;
t319 = -t351 * t402 + t352 * t403 + t454;
t318 = t341 * t405 - t369 * t375 + t442;
t317 = -t342 * t405 + t369 * t374 + t443;
t316 = -t341 * t374 + t342 * t375 + t446;
t315 = t322 * t405 + t333 * t394 - t355 * t375 - t360 * t364 + t442;
t314 = -t323 * t405 - t334 * t394 + t355 * t374 + t359 * t364 + t443;
t313 = -t322 * t374 + t323 * t375 - t333 * t359 + t334 * t360 + t446;
t1 = ((t408 * t439 + t436 * t447) * qJD(1) + (t439 ^ 2 * t377 + (t448 * t436 + (t378 - t449) * t439) * t436) * qJD(3)) * t426 / 0.2e1 + m(7) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(6) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(5) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(4) * (t324 ^ 2 + t325 ^ 2 + t344 ^ 2) / 0.2e1 + m(3) * (t357 ^ 2 + t358 ^ 2) / 0.2e1 + t360 * ((-t326 * t465 + t370 * t328 + t371 * t330) * t360 + (-t327 * t465 + t329 * t370 + t331 * t371) * t359 + (-t361 * t465 + t362 * t370 + t363 * t371) * t394) / 0.2e1 + t359 * ((t326 * t463 + t328 * t372 + t330 * t373) * t360 + (t327 * t463 + t372 * t329 + t373 * t331) * t359 + (t361 * t463 + t362 * t372 + t363 * t373) * t394) / 0.2e1 + t394 * ((t326 * t360 + t327 * t359 + t361 * t394) * t435 + ((-t328 * t419 + t330 * t420) * t360 + (-t329 * t419 + t331 * t420) * t359 + (-t362 * t419 + t363 * t420) * t394) * t438) / 0.2e1 + t405 * ((t335 * t375 + t336 * t374 + t366 * t405) * t435 + ((-t337 * t428 + t339 * t429) * t375 + (-t338 * t428 + t340 * t429) * t374 + (-t367 * t428 + t368 * t429) * t405) * t438) / 0.2e1 + t375 * ((-t335 * t465 + t388 * t337 + t389 * t339) * t375 + (-t336 * t465 + t338 * t388 + t340 * t389) * t374 + (-t366 * t465 + t367 * t388 + t368 * t389) * t405) / 0.2e1 + t374 * ((t335 * t463 + t337 * t390 + t339 * t391) * t375 + (t336 * t463 + t390 * t338 + t391 * t340) * t374 + (t366 * t463 + t367 * t390 + t368 * t391) * t405) / 0.2e1 + t403 * ((-t345 * t465 + t395 * t347 + t396 * t349) * t403 + (-t346 * t465 + t348 * t395 + t350 * t396) * t402 + (-t376 * t465 + t379 * t395 + t382 * t396) * t418) / 0.2e1 + t402 * ((t345 * t463 + t347 * t397 + t349 * t398) * t403 + (t346 * t463 + t397 * t348 + t398 * t350) * t402 + (t376 * t463 + t379 * t397 + t382 * t398) * t418) / 0.2e1 + t418 * ((t345 * t403 + t346 * t402 + t376 * t418) * t435 + ((-t347 * t434 + t349 * t437) * t403 + (-t348 * t434 + t350 * t437) * t402 + (-t379 * t434 + t382 * t437) * t418) * t438) / 0.2e1 + ((t436 * t408 - t439 * t447) * qJD(1) + (t436 ^ 2 * t378 + (t449 * t439 + (t377 - t448) * t436) * t439) * qJD(3)) * t425 / 0.2e1 + qJD(1) * ((-t435 * t409 + t438 * t410) * qJD(1) + ((-t380 * t435 + t383 * t438) * t439 + (-t381 * t435 + t384 * t438) * t436) * qJD(3)) / 0.2e1 + (m(2) * (t412 ^ 2 + t414 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
