% Calculate kinetic energy for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:12
% EndTime: 2019-03-09 05:00:15
% DurationCPUTime: 2.45s
% Computational Cost: add. (2041->290), mult. (1897->447), div. (0->0), fcn. (1886->12), ass. (0->141)
t484 = -Icges(6,3) - Icges(5,3);
t430 = qJ(4) + pkin(11);
t423 = sin(t430);
t425 = cos(t430);
t431 = qJ(1) + pkin(10);
t426 = cos(t431);
t424 = sin(t431);
t437 = cos(qJ(3));
t466 = t437 * t424;
t378 = -t423 * t466 - t425 * t426;
t379 = -t423 * t426 + t425 * t466;
t433 = sin(qJ(4));
t436 = cos(qJ(4));
t389 = -t426 * t436 - t433 * t466;
t468 = t426 * t433;
t390 = t436 * t466 - t468;
t434 = sin(qJ(3));
t469 = t424 * t434;
t483 = Icges(5,5) * t390 + Icges(6,5) * t379 + Icges(5,6) * t389 + Icges(6,6) * t378 - t469 * t484;
t465 = t437 * t426;
t380 = -t423 * t465 + t424 * t425;
t381 = t423 * t424 + t425 * t465;
t391 = t424 * t436 - t433 * t465;
t470 = t424 * t433;
t392 = t436 * t465 + t470;
t467 = t426 * t434;
t482 = Icges(5,5) * t392 + Icges(6,5) * t381 + Icges(5,6) * t391 + Icges(6,6) * t380 - t467 * t484;
t481 = t484 * t437 + (Icges(5,5) * t436 + Icges(6,5) * t425 - Icges(5,6) * t433 - Icges(6,6) * t423) * t434;
t435 = sin(qJ(1));
t476 = pkin(1) * t435;
t474 = t436 * pkin(4);
t472 = Icges(4,4) * t434;
t471 = Icges(4,4) * t437;
t438 = cos(qJ(1));
t422 = qJD(1) * t438 * pkin(1);
t464 = qJD(1) * (pkin(2) * t426 + pkin(7) * t424) + t422;
t463 = pkin(5) * t425;
t417 = qJD(3) * t424;
t460 = qJD(4) * t434;
t399 = t426 * t460 + t417;
t461 = qJD(3) * t426;
t459 = qJD(5) * t434;
t458 = qJD(6) * t434;
t454 = pkin(3) * t437 + pkin(8) * t434;
t397 = t454 * t424;
t398 = t454 * t426;
t457 = t397 * t417 + t398 * t461 + qJD(2);
t456 = -pkin(2) * t424 + pkin(7) * t426 - t476;
t455 = pkin(5) * t423;
t400 = t424 * t460 - t461;
t453 = rSges(4,1) * t437 - rSges(4,2) * t434;
t452 = Icges(4,1) * t437 - t472;
t451 = -Icges(4,2) * t434 + t471;
t450 = Icges(4,5) * t437 - Icges(4,6) * t434;
t368 = -Icges(4,6) * t426 + t424 * t451;
t370 = -Icges(4,5) * t426 + t424 * t452;
t449 = t368 * t434 - t370 * t437;
t369 = Icges(4,6) * t424 + t426 * t451;
t371 = Icges(4,5) * t424 + t426 * t452;
t448 = -t369 * t434 + t371 * t437;
t407 = Icges(4,2) * t437 + t472;
t408 = Icges(4,1) * t434 + t471;
t447 = -t407 * t434 + t408 * t437;
t416 = pkin(3) * t434 - pkin(8) * t437;
t446 = qJD(1) * t398 - t416 * t417 + t464;
t444 = qJ(5) * t434 + t437 * t474;
t351 = -pkin(4) * t468 + t424 * t444;
t445 = -qJD(5) * t437 + t399 * t351 + t457;
t443 = pkin(9) * t434 + t437 * t463;
t352 = pkin(4) * t470 + t426 * t444;
t418 = -qJD(4) * t437 + qJD(1);
t442 = t418 * t352 + t424 * t459 + t446;
t441 = (-t397 + t456) * qJD(1) - t416 * t461;
t377 = -qJ(5) * t437 + t434 * t474;
t440 = t400 * t377 + t426 * t459 + t441;
t427 = qJ(6) + t430;
t420 = cos(t427);
t419 = sin(t427);
t415 = rSges(2,1) * t438 - rSges(2,2) * t435;
t414 = rSges(2,1) * t435 + rSges(2,2) * t438;
t413 = rSges(4,1) * t434 + rSges(4,2) * t437;
t406 = Icges(4,5) * t434 + Icges(4,6) * t437;
t404 = qJD(1) + (-qJD(4) - qJD(6)) * t437;
t396 = -rSges(5,3) * t437 + (rSges(5,1) * t436 - rSges(5,2) * t433) * t434;
t395 = -Icges(5,5) * t437 + (Icges(5,1) * t436 - Icges(5,4) * t433) * t434;
t394 = -Icges(5,6) * t437 + (Icges(5,4) * t436 - Icges(5,2) * t433) * t434;
t387 = t422 + qJD(1) * (rSges(3,1) * t426 - rSges(3,2) * t424);
t386 = (-rSges(3,1) * t424 - rSges(3,2) * t426 - t476) * qJD(1);
t385 = -rSges(6,3) * t437 + (rSges(6,1) * t425 - rSges(6,2) * t423) * t434;
t384 = -Icges(6,5) * t437 + (Icges(6,1) * t425 - Icges(6,4) * t423) * t434;
t383 = -Icges(6,6) * t437 + (Icges(6,4) * t425 - Icges(6,2) * t423) * t434;
t374 = rSges(4,3) * t424 + t426 * t453;
t373 = -rSges(4,3) * t426 + t424 * t453;
t372 = -rSges(7,3) * t437 + (rSges(7,1) * t420 - rSges(7,2) * t419) * t434;
t367 = Icges(4,3) * t424 + t426 * t450;
t366 = -Icges(4,3) * t426 + t424 * t450;
t365 = -Icges(7,5) * t437 + (Icges(7,1) * t420 - Icges(7,4) * t419) * t434;
t364 = -Icges(7,6) * t437 + (Icges(7,4) * t420 - Icges(7,2) * t419) * t434;
t363 = -Icges(7,3) * t437 + (Icges(7,5) * t420 - Icges(7,6) * t419) * t434;
t362 = t419 * t424 + t420 * t465;
t361 = -t419 * t465 + t420 * t424;
t360 = -t419 * t426 + t420 * t466;
t359 = -t419 * t466 - t420 * t426;
t358 = t424 * t458 + t400;
t357 = t426 * t458 + t399;
t356 = -pkin(9) * t437 + t434 * t463;
t354 = rSges(5,1) * t392 + rSges(5,2) * t391 + rSges(5,3) * t467;
t353 = rSges(5,1) * t390 + rSges(5,2) * t389 + rSges(5,3) * t469;
t350 = Icges(5,1) * t392 + Icges(5,4) * t391 + Icges(5,5) * t467;
t349 = Icges(5,1) * t390 + Icges(5,4) * t389 + Icges(5,5) * t469;
t348 = Icges(5,4) * t392 + Icges(5,2) * t391 + Icges(5,6) * t467;
t347 = Icges(5,4) * t390 + Icges(5,2) * t389 + Icges(5,6) * t469;
t343 = rSges(6,1) * t381 + rSges(6,2) * t380 + rSges(6,3) * t467;
t342 = rSges(6,1) * t379 + rSges(6,2) * t378 + rSges(6,3) * t469;
t341 = Icges(6,1) * t381 + Icges(6,4) * t380 + Icges(6,5) * t467;
t340 = Icges(6,1) * t379 + Icges(6,4) * t378 + Icges(6,5) * t469;
t339 = Icges(6,4) * t381 + Icges(6,2) * t380 + Icges(6,6) * t467;
t338 = Icges(6,4) * t379 + Icges(6,2) * t378 + Icges(6,6) * t469;
t335 = qJD(1) * t374 - t413 * t417 + t464;
t334 = -t413 * t461 + (-t373 + t456) * qJD(1);
t333 = rSges(7,1) * t362 + rSges(7,2) * t361 + rSges(7,3) * t467;
t332 = rSges(7,1) * t360 + rSges(7,2) * t359 + rSges(7,3) * t469;
t331 = Icges(7,1) * t362 + Icges(7,4) * t361 + Icges(7,5) * t467;
t330 = Icges(7,1) * t360 + Icges(7,4) * t359 + Icges(7,5) * t469;
t329 = Icges(7,4) * t362 + Icges(7,2) * t361 + Icges(7,6) * t467;
t328 = Icges(7,4) * t360 + Icges(7,2) * t359 + Icges(7,6) * t469;
t327 = Icges(7,5) * t362 + Icges(7,6) * t361 + Icges(7,3) * t467;
t326 = Icges(7,5) * t360 + Icges(7,6) * t359 + Icges(7,3) * t469;
t325 = qJD(2) + (t373 * t424 + t374 * t426) * qJD(3);
t323 = t424 * t455 + t426 * t443;
t322 = t424 * t443 - t426 * t455;
t321 = t354 * t418 - t396 * t399 + t446;
t320 = -t353 * t418 + t396 * t400 + t441;
t319 = t353 * t399 - t354 * t400 + t457;
t318 = t343 * t418 + (-t377 - t385) * t399 + t442;
t317 = t385 * t400 + (-t342 - t351) * t418 + t440;
t316 = t342 * t399 + (-t343 - t352) * t400 + t445;
t315 = t323 * t418 + t333 * t404 - t357 * t372 + (-t356 - t377) * t399 + t442;
t314 = -t332 * t404 + t356 * t400 + t358 * t372 + (-t322 - t351) * t418 + t440;
t313 = t322 * t399 + t332 * t357 - t333 * t358 + (-t323 - t352) * t400 + t445;
t1 = m(3) * (qJD(2) ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(6) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(7) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(5) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(4) * (t325 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + t357 * ((t327 * t467 + t361 * t329 + t362 * t331) * t357 + (t326 * t467 + t328 * t361 + t330 * t362) * t358 + (t361 * t364 + t362 * t365 + t363 * t467) * t404) / 0.2e1 + t404 * ((-t326 * t358 - t327 * t357 - t363 * t404) * t437 + ((-t329 * t419 + t331 * t420) * t357 + (-t328 * t419 + t330 * t420) * t358 + (-t364 * t419 + t365 * t420) * t404) * t434) / 0.2e1 + t358 * ((t327 * t469 + t329 * t359 + t331 * t360) * t357 + (t326 * t469 + t359 * t328 + t360 * t330) * t358 + (t359 * t364 + t360 * t365 + t363 * t469) * t404) / 0.2e1 - ((-t426 * t406 + t424 * t447) * qJD(1) + (t426 ^ 2 * t366 + (t448 * t424 + (-t367 + t449) * t426) * t424) * qJD(3)) * t461 / 0.2e1 + qJD(1) * ((t437 * t407 + t434 * t408) * qJD(1) + ((t369 * t437 + t371 * t434) * t424 - (t368 * t437 + t370 * t434) * t426) * qJD(3)) / 0.2e1 + ((t424 * t406 + t426 * t447) * qJD(1) + (t424 ^ 2 * t367 + (t449 * t426 + (-t366 + t448) * t424) * t426) * qJD(3)) * t417 / 0.2e1 + ((t380 * t383 + t381 * t384 + t391 * t394 + t392 * t395 + t467 * t481) * t418 + (t338 * t380 + t340 * t381 + t347 * t391 + t349 * t392 + t467 * t483) * t400 + (t380 * t339 + t381 * t341 + t391 * t348 + t392 * t350 + t482 * t467) * t399) * t399 / 0.2e1 + ((t378 * t383 + t379 * t384 + t389 * t394 + t390 * t395 + t469 * t481) * t418 + (t378 * t338 + t379 * t340 + t389 * t347 + t390 * t349 + t483 * t469) * t400 + (t339 * t378 + t341 * t379 + t348 * t389 + t350 * t390 + t469 * t482) * t399) * t400 / 0.2e1 + ((-t482 * t399 - t400 * t483 - t481 * t418) * t437 + ((-t383 * t423 + t384 * t425 - t394 * t433 + t395 * t436) * t418 + (-t338 * t423 + t340 * t425 - t347 * t433 + t349 * t436) * t400 + (-t339 * t423 + t341 * t425 - t348 * t433 + t350 * t436) * t399) * t434) * t418 / 0.2e1 + (m(2) * (t414 ^ 2 + t415 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
