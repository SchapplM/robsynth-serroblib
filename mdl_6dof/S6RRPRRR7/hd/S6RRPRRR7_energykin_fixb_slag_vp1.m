% Calculate kinetic energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:40
% EndTime: 2019-03-09 13:55:44
% DurationCPUTime: 3.87s
% Computational Cost: add. (1414->303), mult. (3055->475), div. (0->0), fcn. (3394->10), ass. (0->152)
t544 = Icges(3,4) - Icges(4,5);
t543 = Icges(3,1) + Icges(4,1);
t542 = Icges(3,2) + Icges(4,3);
t468 = sin(qJ(2));
t541 = t544 * t468;
t471 = cos(qJ(2));
t540 = t544 * t471;
t539 = Icges(4,4) + Icges(3,5);
t538 = Icges(3,6) - Icges(4,6);
t537 = t542 * t468 - t540;
t536 = t543 * t471 - t541;
t535 = Icges(4,2) + Icges(3,3);
t469 = sin(qJ(1));
t472 = cos(qJ(1));
t534 = t469 * t537 + t472 * t538;
t533 = -t469 * t538 + t472 * t537;
t532 = -t536 * t469 + t472 * t539;
t531 = t469 * t539 + t536 * t472;
t530 = -t542 * t471 - t541;
t529 = t543 * t468 + t540;
t528 = -t538 * t468 + t471 * t539;
t527 = t528 * t469 - t472 * t535;
t526 = t469 * t535 + t528 * t472;
t525 = t468 * t539 + t538 * t471;
t524 = t468 * t530 + t471 * t529;
t523 = t468 * t533 + t471 * t531;
t522 = -t468 * t534 + t471 * t532;
t467 = sin(qJ(4));
t518 = cos(qJ(4));
t501 = t468 * t518;
t438 = -t471 * t467 + t501;
t470 = cos(qJ(5));
t516 = pkin(5) * t470;
t466 = sin(qJ(5));
t510 = t466 * t469;
t509 = t466 * t472;
t507 = t471 * t472;
t494 = pkin(2) * t471 + qJ(3) * t468;
t433 = t494 * t469;
t455 = pkin(1) * t469 - pkin(7) * t472;
t506 = -t433 - t455;
t462 = qJD(2) * t469;
t505 = qJD(2) * t472;
t504 = qJD(3) * t468;
t437 = t468 * t467 + t471 * t518;
t429 = qJD(5) * t437 + qJD(1);
t434 = t494 * t472;
t441 = qJD(1) * (pkin(1) * t472 + pkin(7) * t469);
t503 = qJD(1) * t434 + t469 * t504 + t441;
t439 = pkin(3) * t469 * t471 + pkin(8) * t472;
t502 = -t439 + t506;
t450 = pkin(2) * t468 - qJ(3) * t471;
t498 = qJD(2) * (-rSges(4,1) * t468 + rSges(4,3) * t471 - t450);
t443 = qJD(4) * t472 - t505;
t442 = -qJD(4) * t469 + t462;
t425 = t438 * t469;
t395 = -qJD(5) * t425 + t443;
t497 = -qJD(3) * t471 + t433 * t462 + t434 * t505;
t427 = t467 * t507 - t472 * t501;
t394 = qJD(5) * t427 + t442;
t496 = rSges(3,1) * t471 - rSges(3,2) * t468;
t495 = rSges(4,1) * t471 + rSges(4,3) * t468;
t493 = qJD(2) * (-pkin(3) * t468 - t450);
t440 = pkin(3) * t507 - pkin(8) * t469;
t480 = t439 * t462 + t440 * t505 + t497;
t459 = t472 * t504;
t479 = t472 * t493 + t459;
t426 = t437 * t469;
t387 = pkin(4) * t426 - pkin(9) * t425;
t428 = t437 * t472;
t388 = pkin(4) * t428 + pkin(9) * t427;
t478 = t442 * t387 - t388 * t443 + t480;
t477 = qJD(1) * t440 + t469 * t493 + t503;
t400 = pkin(4) * t438 + pkin(9) * t437;
t476 = qJD(1) * t388 - t400 * t442 + t477;
t475 = t443 * t400 + (-t387 + t502) * qJD(1) + t479;
t465 = qJ(5) + qJ(6);
t464 = cos(t465);
t463 = sin(t465);
t454 = rSges(2,1) * t472 - rSges(2,2) * t469;
t453 = rSges(2,1) * t469 + rSges(2,2) * t472;
t452 = rSges(3,1) * t468 + rSges(3,2) * t471;
t424 = rSges(3,3) * t469 + t472 * t496;
t423 = rSges(4,2) * t469 + t472 * t495;
t422 = -rSges(3,3) * t472 + t469 * t496;
t421 = -rSges(4,2) * t472 + t469 * t495;
t404 = t428 * t470 - t510;
t403 = -t428 * t466 - t469 * t470;
t402 = t426 * t470 + t509;
t401 = -t426 * t466 + t470 * t472;
t399 = rSges(5,1) * t438 - rSges(5,2) * t437;
t398 = Icges(5,1) * t438 - Icges(5,4) * t437;
t397 = Icges(5,4) * t438 - Icges(5,2) * t437;
t396 = Icges(5,5) * t438 - Icges(5,6) * t437;
t393 = qJD(6) * t437 + t429;
t392 = t428 * t464 - t463 * t469;
t391 = -t428 * t463 - t464 * t469;
t390 = t426 * t464 + t463 * t472;
t389 = -t426 * t463 + t464 * t472;
t384 = rSges(5,1) * t428 - rSges(5,2) * t427 - rSges(5,3) * t469;
t383 = rSges(5,1) * t426 + rSges(5,2) * t425 + rSges(5,3) * t472;
t382 = Icges(5,1) * t428 - Icges(5,4) * t427 - Icges(5,5) * t469;
t381 = Icges(5,1) * t426 + Icges(5,4) * t425 + Icges(5,5) * t472;
t380 = Icges(5,4) * t428 - Icges(5,2) * t427 - Icges(5,6) * t469;
t379 = Icges(5,4) * t426 + Icges(5,2) * t425 + Icges(5,6) * t472;
t378 = Icges(5,5) * t428 - Icges(5,6) * t427 - Icges(5,3) * t469;
t377 = Icges(5,5) * t426 + Icges(5,6) * t425 + Icges(5,3) * t472;
t376 = qJD(1) * t424 - t452 * t462 + t441;
t375 = -t452 * t505 + (-t422 - t455) * qJD(1);
t374 = (t422 * t469 + t424 * t472) * qJD(2);
t373 = rSges(6,3) * t437 + (rSges(6,1) * t470 - rSges(6,2) * t466) * t438;
t372 = Icges(6,5) * t437 + (Icges(6,1) * t470 - Icges(6,4) * t466) * t438;
t371 = Icges(6,6) * t437 + (Icges(6,4) * t470 - Icges(6,2) * t466) * t438;
t370 = Icges(6,3) * t437 + (Icges(6,5) * t470 - Icges(6,6) * t466) * t438;
t369 = -qJD(6) * t425 + t395;
t368 = qJD(6) * t427 + t394;
t366 = rSges(7,3) * t437 + (rSges(7,1) * t464 - rSges(7,2) * t463) * t438;
t365 = Icges(7,5) * t437 + (Icges(7,1) * t464 - Icges(7,4) * t463) * t438;
t364 = Icges(7,6) * t437 + (Icges(7,4) * t464 - Icges(7,2) * t463) * t438;
t363 = Icges(7,3) * t437 + (Icges(7,5) * t464 - Icges(7,6) * t463) * t438;
t362 = pkin(10) * t437 + t438 * t516;
t361 = rSges(6,1) * t404 + rSges(6,2) * t403 + rSges(6,3) * t427;
t360 = rSges(6,1) * t402 + rSges(6,2) * t401 - rSges(6,3) * t425;
t359 = Icges(6,1) * t404 + Icges(6,4) * t403 + Icges(6,5) * t427;
t358 = Icges(6,1) * t402 + Icges(6,4) * t401 - Icges(6,5) * t425;
t357 = Icges(6,4) * t404 + Icges(6,2) * t403 + Icges(6,6) * t427;
t356 = Icges(6,4) * t402 + Icges(6,2) * t401 - Icges(6,6) * t425;
t355 = Icges(6,5) * t404 + Icges(6,6) * t403 + Icges(6,3) * t427;
t354 = Icges(6,5) * t402 + Icges(6,6) * t401 - Icges(6,3) * t425;
t353 = qJD(1) * t423 + t469 * t498 + t503;
t352 = t459 + t472 * t498 + (-t421 + t506) * qJD(1);
t351 = rSges(7,1) * t392 + rSges(7,2) * t391 + rSges(7,3) * t427;
t350 = rSges(7,1) * t390 + rSges(7,2) * t389 - rSges(7,3) * t425;
t349 = Icges(7,1) * t392 + Icges(7,4) * t391 + Icges(7,5) * t427;
t348 = Icges(7,1) * t390 + Icges(7,4) * t389 - Icges(7,5) * t425;
t347 = Icges(7,4) * t392 + Icges(7,2) * t391 + Icges(7,6) * t427;
t346 = Icges(7,4) * t390 + Icges(7,2) * t389 - Icges(7,6) * t425;
t345 = Icges(7,5) * t392 + Icges(7,6) * t391 + Icges(7,3) * t427;
t344 = Icges(7,5) * t390 + Icges(7,6) * t389 - Icges(7,3) * t425;
t343 = -pkin(5) * t510 + pkin(10) * t427 + t428 * t516;
t342 = pkin(5) * t509 - pkin(10) * t425 + t426 * t516;
t341 = (t421 * t469 + t423 * t472) * qJD(2) + t497;
t340 = qJD(1) * t384 - t399 * t442 + t477;
t339 = t399 * t443 + (-t383 + t502) * qJD(1) + t479;
t338 = t383 * t442 - t384 * t443 + t480;
t337 = t361 * t429 - t373 * t394 + t476;
t336 = -t360 * t429 + t373 * t395 + t475;
t335 = t360 * t394 - t361 * t395 + t478;
t334 = t343 * t429 + t351 * t393 - t362 * t394 - t366 * t368 + t476;
t333 = -t342 * t429 - t350 * t393 + t362 * t395 + t366 * t369 + t475;
t332 = t342 * t394 - t343 * t395 + t350 * t368 - t351 * t369 + t478;
t1 = t393 * ((t344 * t369 + t345 * t368 + t363 * t393) * t437 + ((-t347 * t463 + t349 * t464) * t368 + (-t346 * t463 + t348 * t464) * t369 + (-t364 * t463 + t365 * t464) * t393) * t438) / 0.2e1 + t369 * ((-t345 * t425 + t347 * t389 + t349 * t390) * t368 + (-t425 * t344 + t389 * t346 + t390 * t348) * t369 + (-t363 * t425 + t364 * t389 + t365 * t390) * t393) / 0.2e1 + t368 * ((t427 * t345 + t391 * t347 + t392 * t349) * t368 + (t344 * t427 + t346 * t391 + t348 * t392) * t369 + (t363 * t427 + t364 * t391 + t365 * t392) * t393) / 0.2e1 + t395 * ((-t355 * t425 + t357 * t401 + t359 * t402) * t394 + (-t425 * t354 + t401 * t356 + t402 * t358) * t395 + (-t370 * t425 + t371 * t401 + t372 * t402) * t429) / 0.2e1 + t429 * ((t354 * t395 + t355 * t394 + t370 * t429) * t437 + ((-t357 * t466 + t359 * t470) * t394 + (-t356 * t466 + t358 * t470) * t395 + (-t371 * t466 + t372 * t470) * t429) * t438) / 0.2e1 + t394 * ((t427 * t355 + t403 * t357 + t404 * t359) * t394 + (t354 * t427 + t356 * t403 + t358 * t404) * t395 + (t370 * t427 + t371 * t403 + t372 * t404) * t429) / 0.2e1 + t442 * ((-t469 * t378 - t427 * t380 + t428 * t382) * t442 + (-t377 * t469 - t379 * t427 + t381 * t428) * t443 + (-t396 * t469 - t397 * t427 + t398 * t428) * qJD(1)) / 0.2e1 + t443 * ((t378 * t472 + t380 * t425 + t382 * t426) * t442 + (t472 * t377 + t425 * t379 + t426 * t381) * t443 + (t396 * t472 + t397 * t425 + t398 * t426) * qJD(1)) / 0.2e1 + m(6) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(7) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(4) * (t341 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(5) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(3) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + (m(2) * (t453 ^ 2 + t454 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t526 * t469 ^ 2 + (t522 * t472 + (t523 - t527) * t469) * t472) * qJD(2) + (t469 * t525 + t472 * t524) * qJD(1)) * t462 / 0.2e1 - ((t527 * t472 ^ 2 + (t523 * t469 + (t522 - t526) * t472) * t469) * qJD(2) + (t469 * t524 - t472 * t525) * qJD(1)) * t505 / 0.2e1 + ((-t380 * t437 + t382 * t438) * t442 + (-t379 * t437 + t381 * t438) * t443 + ((t468 * t532 + t471 * t534) * t472 + (t468 * t531 - t471 * t533) * t469) * qJD(2) + (-t437 * t397 + t438 * t398 + t529 * t468 - t530 * t471) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
