% Calculate kinetic energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:24
% EndTime: 2019-03-09 03:44:27
% DurationCPUTime: 2.83s
% Computational Cost: add. (1445->247), mult. (1548->394), div. (0->0), fcn. (1481->10), ass. (0->132)
t512 = Icges(4,4) + Icges(5,6);
t511 = Icges(4,1) + Icges(5,2);
t510 = -Icges(4,2) - Icges(5,3);
t437 = cos(qJ(3));
t509 = t512 * t437;
t434 = sin(qJ(3));
t508 = t512 * t434;
t507 = -Icges(5,4) + Icges(4,5);
t506 = Icges(5,5) - Icges(4,6);
t505 = t510 * t434 + t509;
t504 = -t511 * t437 + t508;
t503 = Icges(5,1) + Icges(4,3);
t431 = qJ(1) + pkin(10);
t426 = sin(t431);
t427 = cos(t431);
t502 = t505 * t426 + t506 * t427;
t501 = -t506 * t426 + t505 * t427;
t500 = t504 * t426 + t507 * t427;
t499 = t507 * t426 - t504 * t427;
t498 = t510 * t437 - t508;
t497 = t511 * t434 + t509;
t496 = t506 * t434 + t507 * t437;
t495 = t496 * t426 - t427 * t503;
t494 = t426 * t503 + t496 * t427;
t493 = t507 * t434 - t506 * t437;
t492 = t434 * t498 + t437 * t497;
t491 = -t434 * t501 + t437 * t499;
t490 = t434 * t502 + t437 * t500;
t435 = sin(qJ(1));
t486 = pkin(1) * t435;
t436 = cos(qJ(5));
t484 = pkin(5) * t436;
t478 = t426 * t437;
t477 = t427 * t437;
t432 = qJ(5) + qJ(6);
t429 = sin(t432);
t476 = t429 * t434;
t430 = cos(t432);
t475 = t430 * t434;
t433 = sin(qJ(5));
t474 = t433 * t434;
t473 = t434 * t436;
t438 = cos(qJ(1));
t425 = qJD(1) * t438 * pkin(1);
t472 = qJD(1) * (pkin(2) * t427 + pkin(7) * t426) + t425;
t422 = qJD(3) * t426;
t469 = qJD(5) * t437;
t400 = t427 * t469 + t422;
t471 = qJD(3) * t427;
t470 = qJD(4) * t434;
t468 = qJD(6) * t437;
t423 = qJD(5) * t434 + qJD(1);
t465 = -pkin(2) * t426 + pkin(7) * t427 - t486;
t417 = pkin(3) * t434 - qJ(4) * t437;
t464 = qJD(3) * (rSges(5,2) * t434 + rSges(5,3) * t437 - t417);
t401 = t426 * t469 - t471;
t459 = pkin(3) * t437 + qJ(4) * t434;
t396 = t459 * t427;
t463 = qJD(1) * t396 + t426 * t470 + t472;
t395 = t459 * t426;
t462 = -t395 + t465;
t461 = rSges(4,1) * t437 - rSges(4,2) * t434;
t460 = -rSges(5,2) * t437 + rSges(5,3) * t434;
t458 = qJD(3) * (-pkin(8) * t434 - t417);
t445 = -qJD(4) * t437 + t395 * t422 + t396 * t471 + qJD(2);
t444 = pkin(5) * t474 + pkin(9) * t437;
t402 = pkin(4) * t426 + pkin(8) * t477;
t403 = -pkin(4) * t427 + pkin(8) * t478;
t443 = t402 * t471 + t403 * t422 + t445;
t442 = qJD(1) * t402 + t426 * t458 + t463;
t416 = t427 * t470;
t441 = t416 + (-t403 + t462) * qJD(1) + t427 * t458;
t421 = rSges(2,1) * t438 - rSges(2,2) * t435;
t420 = rSges(2,1) * t435 + rSges(2,2) * t438;
t419 = rSges(4,1) * t434 + rSges(4,2) * t437;
t406 = qJD(6) * t434 + t423;
t399 = -pkin(5) * t433 * t437 + pkin(9) * t434;
t397 = rSges(6,3) * t434 + (-rSges(6,1) * t433 - rSges(6,2) * t436) * t437;
t394 = Icges(6,5) * t434 + (-Icges(6,1) * t433 - Icges(6,4) * t436) * t437;
t393 = Icges(6,6) * t434 + (-Icges(6,4) * t433 - Icges(6,2) * t436) * t437;
t392 = Icges(6,3) * t434 + (-Icges(6,5) * t433 - Icges(6,6) * t436) * t437;
t391 = t426 * t474 - t427 * t436;
t390 = t426 * t473 + t427 * t433;
t389 = t426 * t436 + t427 * t474;
t388 = -t426 * t433 + t427 * t473;
t386 = t425 + qJD(1) * (rSges(3,1) * t427 - rSges(3,2) * t426);
t385 = (-rSges(3,1) * t426 - rSges(3,2) * t427 - t486) * qJD(1);
t384 = rSges(7,3) * t434 + (-rSges(7,1) * t429 - rSges(7,2) * t430) * t437;
t381 = Icges(7,5) * t434 + (-Icges(7,1) * t429 - Icges(7,4) * t430) * t437;
t380 = Icges(7,6) * t434 + (-Icges(7,4) * t429 - Icges(7,2) * t430) * t437;
t379 = Icges(7,3) * t434 + (-Icges(7,5) * t429 - Icges(7,6) * t430) * t437;
t378 = t426 * t476 - t427 * t430;
t377 = t426 * t475 + t427 * t429;
t376 = t426 * t430 + t427 * t476;
t375 = -t426 * t429 + t427 * t475;
t374 = -rSges(5,1) * t427 + t426 * t460;
t373 = rSges(5,1) * t426 + t427 * t460;
t372 = rSges(4,3) * t426 + t427 * t461;
t371 = -rSges(4,3) * t427 + t426 * t461;
t356 = t426 * t468 + t401;
t355 = t427 * t468 + t400;
t354 = t426 * t444 - t427 * t484;
t353 = t426 * t484 + t427 * t444;
t352 = rSges(6,1) * t391 + rSges(6,2) * t390 + rSges(6,3) * t478;
t351 = rSges(6,1) * t389 + rSges(6,2) * t388 + rSges(6,3) * t477;
t350 = Icges(6,1) * t391 + Icges(6,4) * t390 + Icges(6,5) * t478;
t349 = Icges(6,1) * t389 + Icges(6,4) * t388 + Icges(6,5) * t477;
t348 = Icges(6,4) * t391 + Icges(6,2) * t390 + Icges(6,6) * t478;
t347 = Icges(6,4) * t389 + Icges(6,2) * t388 + Icges(6,6) * t477;
t346 = Icges(6,5) * t391 + Icges(6,6) * t390 + Icges(6,3) * t478;
t345 = Icges(6,5) * t389 + Icges(6,6) * t388 + Icges(6,3) * t477;
t344 = rSges(7,1) * t378 + rSges(7,2) * t377 + rSges(7,3) * t478;
t343 = rSges(7,1) * t376 + rSges(7,2) * t375 + rSges(7,3) * t477;
t342 = Icges(7,1) * t378 + Icges(7,4) * t377 + Icges(7,5) * t478;
t341 = Icges(7,1) * t376 + Icges(7,4) * t375 + Icges(7,5) * t477;
t340 = Icges(7,4) * t378 + Icges(7,2) * t377 + Icges(7,6) * t478;
t339 = Icges(7,4) * t376 + Icges(7,2) * t375 + Icges(7,6) * t477;
t338 = Icges(7,5) * t378 + Icges(7,6) * t377 + Icges(7,3) * t478;
t337 = Icges(7,5) * t376 + Icges(7,6) * t375 + Icges(7,3) * t477;
t336 = qJD(1) * t372 - t419 * t422 + t472;
t335 = -t419 * t471 + (-t371 + t465) * qJD(1);
t334 = qJD(2) + (t371 * t426 + t372 * t427) * qJD(3);
t333 = qJD(1) * t373 + t426 * t464 + t463;
t332 = t416 + t427 * t464 + (-t374 + t462) * qJD(1);
t331 = (t373 * t427 + t374 * t426) * qJD(3) + t445;
t330 = t351 * t423 - t397 * t400 + t442;
t329 = -t352 * t423 + t397 * t401 + t441;
t328 = -t351 * t401 + t352 * t400 + t443;
t327 = t343 * t406 + t353 * t423 - t355 * t384 - t399 * t400 + t442;
t326 = -t344 * t406 - t354 * t423 + t356 * t384 + t399 * t401 + t441;
t325 = -t343 * t356 + t344 * t355 - t353 * t401 + t354 * t400 + t443;
t1 = m(7) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + t406 * ((t337 * t355 + t338 * t356 + t379 * t406) * t434 + ((-t339 * t430 - t341 * t429) * t355 + (-t340 * t430 - t342 * t429) * t356 + (-t380 * t430 - t381 * t429) * t406) * t437) / 0.2e1 + t355 * ((t337 * t477 + t375 * t339 + t376 * t341) * t355 + (t338 * t477 + t340 * t375 + t342 * t376) * t356 + (t375 * t380 + t376 * t381 + t379 * t477) * t406) / 0.2e1 + t356 * ((t337 * t478 + t339 * t377 + t341 * t378) * t355 + (t338 * t478 + t377 * t340 + t378 * t342) * t356 + (t377 * t380 + t378 * t381 + t379 * t478) * t406) / 0.2e1 + t401 * ((t345 * t478 + t347 * t390 + t349 * t391) * t400 + (t346 * t478 + t390 * t348 + t391 * t350) * t401 + (t390 * t393 + t391 * t394 + t392 * t478) * t423) / 0.2e1 + t400 * ((t345 * t477 + t388 * t347 + t389 * t349) * t400 + (t346 * t477 + t348 * t388 + t350 * t389) * t401 + (t388 * t393 + t389 * t394 + t392 * t477) * t423) / 0.2e1 + t423 * ((t345 * t400 + t346 * t401 + t392 * t423) * t434 + ((-t347 * t436 - t349 * t433) * t400 + (-t348 * t436 - t350 * t433) * t401 + (-t393 * t436 - t394 * t433) * t423) * t437) / 0.2e1 + m(5) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(6) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(4) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + (((t434 * t500 - t437 * t502) * t427 + (t434 * t499 + t437 * t501) * t426) * qJD(3) + (t497 * t434 - t498 * t437) * qJD(1)) * qJD(1) / 0.2e1 + ((t494 * t426 ^ 2 + (t490 * t427 + (t491 - t495) * t426) * t427) * qJD(3) + (t426 * t493 + t427 * t492) * qJD(1)) * t422 / 0.2e1 - ((t495 * t427 ^ 2 + (t491 * t426 + (t490 - t494) * t427) * t426) * qJD(3) + (t426 * t492 - t493 * t427) * qJD(1)) * t471 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t420 ^ 2 + t421 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
