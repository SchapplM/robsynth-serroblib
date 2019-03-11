% Calculate kinetic energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:22
% EndTime: 2019-03-09 15:28:25
% DurationCPUTime: 2.87s
% Computational Cost: add. (1221->238), mult. (1626->363), div. (0->0), fcn. (1471->8), ass. (0->141)
t531 = Icges(4,4) + Icges(6,4) - Icges(5,5);
t530 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t529 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t438 = qJ(2) + qJ(3);
t436 = sin(t438);
t528 = t531 * t436;
t437 = cos(t438);
t527 = t531 * t437;
t526 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t525 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t524 = t436 * t529 - t527;
t523 = t437 * t530 - t528;
t522 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t441 = sin(qJ(1));
t444 = cos(qJ(1));
t521 = t524 * t441 + t525 * t444;
t520 = -t525 * t441 + t524 * t444;
t519 = t523 * t441 - t526 * t444;
t518 = t526 * t441 + t523 * t444;
t517 = -t437 * t529 - t528;
t516 = t436 * t530 + t527;
t515 = -t525 * t436 + t526 * t437;
t435 = qJD(2) * t441;
t420 = qJD(3) * t441 + t435;
t421 = (-qJD(2) - qJD(3)) * t444;
t514 = (t521 * t436 + t519 * t437) * t421 + (t520 * t436 + t518 * t437) * t420 + (t517 * t436 + t516 * t437) * qJD(1);
t513 = (t515 * t441 - t522 * t444) * t421 + (t522 * t441 + t515 * t444) * t420 + (t526 * t436 + t525 * t437) * qJD(1);
t509 = pkin(4) * t436;
t443 = cos(qJ(2));
t507 = pkin(2) * t443;
t440 = sin(qJ(2));
t505 = Icges(3,4) * t440;
t504 = Icges(3,4) * t443;
t497 = t437 * t441;
t496 = t437 * t444;
t439 = sin(qJ(6));
t495 = t439 * t441;
t494 = t439 * t444;
t442 = cos(qJ(6));
t493 = t441 * t442;
t492 = t442 * t444;
t356 = -pkin(8) * t444 + t441 * t507;
t357 = pkin(8) * t441 + t444 * t507;
t488 = qJD(2) * t444;
t491 = t356 * t435 + t357 * t488;
t430 = pkin(1) * t441 - pkin(7) * t444;
t490 = -t356 - t430;
t475 = pkin(3) * t437 + qJ(4) * t436;
t394 = t475 * t444;
t404 = pkin(4) * t496 - qJ(5) * t441;
t489 = -t394 - t404;
t487 = qJD(4) * t436;
t486 = qJD(6) * t437;
t485 = pkin(2) * qJD(2) * t440;
t393 = t475 * t441;
t484 = -t393 + t490;
t414 = pkin(3) * t436 - qJ(4) * t437;
t483 = -t414 - t509;
t482 = t444 * t485;
t403 = pkin(4) * t497 + qJ(5) * t444;
t481 = -t403 + t484;
t480 = pkin(5) * t436 + pkin(9) * t437;
t479 = rSges(3,1) * t443 - rSges(3,2) * t440;
t478 = rSges(4,1) * t437 - rSges(4,2) * t436;
t477 = rSges(5,1) * t437 + rSges(5,3) * t436;
t476 = rSges(6,1) * t436 - rSges(6,2) * t437;
t474 = Icges(3,1) * t443 - t505;
t470 = -Icges(3,2) * t440 + t504;
t466 = Icges(3,5) * t443 - Icges(3,6) * t440;
t385 = -Icges(3,6) * t444 + t441 * t470;
t387 = -Icges(3,5) * t444 + t441 * t474;
t462 = t385 * t440 - t387 * t443;
t386 = Icges(3,6) * t441 + t444 * t470;
t388 = Icges(3,5) * t441 + t444 * t474;
t461 = -t386 * t440 + t388 * t443;
t423 = Icges(3,2) * t443 + t505;
t424 = Icges(3,1) * t440 + t504;
t460 = -t423 * t440 + t424 * t443;
t459 = -qJD(4) * t437 + t420 * t393 + t491;
t419 = qJD(1) * (pkin(1) * t444 + pkin(7) * t441);
t458 = qJD(1) * t357 - t441 * t485 + t419;
t457 = t421 * t414 + t444 * t487 - t482;
t456 = t420 * t403 + t459;
t452 = qJD(1) * t394 + t441 * t487 + t458;
t451 = -qJD(5) * t441 + t421 * t509 + t457;
t450 = qJD(1) * t404 + qJD(5) * t444 + t452;
t431 = qJD(6) * t436 + qJD(1);
t427 = rSges(2,1) * t444 - rSges(2,2) * t441;
t426 = rSges(2,1) * t441 + rSges(2,2) * t444;
t425 = rSges(3,1) * t440 + rSges(3,2) * t443;
t422 = Icges(3,5) * t440 + Icges(3,6) * t443;
t418 = -pkin(5) * t437 + pkin(9) * t436;
t417 = -rSges(6,1) * t437 - rSges(6,2) * t436;
t416 = rSges(4,1) * t436 + rSges(4,2) * t437;
t415 = rSges(5,1) * t436 - rSges(5,3) * t437;
t400 = t436 * t492 - t495;
t399 = -t436 * t494 - t493;
t398 = t436 * t493 + t494;
t397 = -t436 * t495 + t492;
t396 = t480 * t444;
t395 = t480 * t441;
t392 = rSges(3,3) * t441 + t444 * t479;
t391 = -rSges(3,3) * t444 + t441 * t479;
t390 = t441 * t486 + t421;
t389 = t444 * t486 + t420;
t384 = Icges(3,3) * t441 + t444 * t466;
t383 = -Icges(3,3) * t444 + t441 * t466;
t381 = rSges(4,3) * t441 + t444 * t478;
t380 = rSges(5,2) * t441 + t444 * t477;
t379 = -rSges(6,3) * t441 + t444 * t476;
t378 = -rSges(4,3) * t444 + t441 * t478;
t377 = -rSges(5,2) * t444 + t441 * t477;
t376 = rSges(6,3) * t444 + t441 * t476;
t353 = rSges(7,3) * t436 + (-rSges(7,1) * t442 + rSges(7,2) * t439) * t437;
t352 = Icges(7,5) * t436 + (-Icges(7,1) * t442 + Icges(7,4) * t439) * t437;
t351 = Icges(7,6) * t436 + (-Icges(7,4) * t442 + Icges(7,2) * t439) * t437;
t350 = Icges(7,3) * t436 + (-Icges(7,5) * t442 + Icges(7,6) * t439) * t437;
t345 = qJD(1) * t392 - t425 * t435 + t419;
t344 = -t425 * t488 + (-t391 - t430) * qJD(1);
t343 = rSges(7,1) * t400 + rSges(7,2) * t399 + rSges(7,3) * t496;
t342 = rSges(7,1) * t398 + rSges(7,2) * t397 + rSges(7,3) * t497;
t341 = Icges(7,1) * t400 + Icges(7,4) * t399 + Icges(7,5) * t496;
t340 = Icges(7,1) * t398 + Icges(7,4) * t397 + Icges(7,5) * t497;
t339 = Icges(7,4) * t400 + Icges(7,2) * t399 + Icges(7,6) * t496;
t338 = Icges(7,4) * t398 + Icges(7,2) * t397 + Icges(7,6) * t497;
t337 = Icges(7,5) * t400 + Icges(7,6) * t399 + Icges(7,3) * t496;
t336 = Icges(7,5) * t398 + Icges(7,6) * t397 + Icges(7,3) * t497;
t335 = (t391 * t441 + t392 * t444) * qJD(2);
t334 = qJD(1) * t381 - t416 * t420 + t458;
t333 = -t482 + t416 * t421 + (-t378 + t490) * qJD(1);
t332 = t378 * t420 - t381 * t421 + t491;
t331 = qJD(1) * t380 + (-t414 - t415) * t420 + t452;
t330 = t415 * t421 + (-t377 + t484) * qJD(1) + t457;
t329 = qJD(1) * t379 + (-t417 + t483) * t420 + t450;
t328 = t417 * t421 + (-t376 + t481) * qJD(1) + t451;
t327 = t377 * t420 + (-t380 - t394) * t421 + t459;
t326 = t376 * t420 + (-t379 + t489) * t421 + t456;
t325 = qJD(1) * t396 + t343 * t431 - t353 * t389 + (-t418 + t483) * t420 + t450;
t324 = -t342 * t431 + t353 * t390 + t418 * t421 + (-t395 + t481) * qJD(1) + t451;
t323 = t342 * t389 - t343 * t390 + t395 * t420 + (-t396 + t489) * t421 + t456;
t1 = -((-t444 * t422 + t441 * t460) * qJD(1) + (t444 ^ 2 * t383 + (t461 * t441 + (-t384 + t462) * t444) * t441) * qJD(2)) * t488 / 0.2e1 + ((t441 * t422 + t444 * t460) * qJD(1) + (t441 ^ 2 * t384 + (t462 * t444 + (-t383 + t461) * t441) * t444) * qJD(2)) * t435 / 0.2e1 + t390 * ((t337 * t497 + t339 * t397 + t341 * t398) * t389 + (t336 * t497 + t397 * t338 + t398 * t340) * t390 + (t350 * t497 + t351 * t397 + t352 * t398) * t431) / 0.2e1 + t431 * ((t336 * t390 + t337 * t389 + t350 * t431) * t436 + ((t339 * t439 - t341 * t442) * t389 + (t338 * t439 - t340 * t442) * t390 + (t351 * t439 - t352 * t442) * t431) * t437) / 0.2e1 + m(3) * (t335 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(5) * (t327 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(6) * (t326 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(7) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + t389 * ((t337 * t496 + t399 * t339 + t400 * t341) * t389 + (t336 * t496 + t338 * t399 + t340 * t400) * t390 + (t350 * t496 + t351 * t399 + t352 * t400) * t431) / 0.2e1 + (Icges(2,3) + m(2) * (t426 ^ 2 + t427 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (t513 * t441 + t514 * t444) * t420 / 0.2e1 + (t514 * t441 - t513 * t444) * t421 / 0.2e1 + (((t386 * t443 + t388 * t440) * t441 - (t385 * t443 + t387 * t440) * t444) * qJD(2) + (t519 * t436 - t521 * t437) * t421 + (t518 * t436 - t520 * t437) * t420 + (t443 * t423 + t440 * t424 + t516 * t436 - t517 * t437) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
