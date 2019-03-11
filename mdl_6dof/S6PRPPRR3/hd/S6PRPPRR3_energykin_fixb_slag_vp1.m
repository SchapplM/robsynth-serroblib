% Calculate kinetic energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:25
% EndTime: 2019-03-08 19:21:28
% DurationCPUTime: 3.39s
% Computational Cost: add. (2393->320), mult. (6198->489), div. (0->0), fcn. (7697->12), ass. (0->146)
t559 = Icges(3,1) + Icges(4,1);
t557 = Icges(3,4) - Icges(4,5);
t556 = Icges(4,4) + Icges(3,5);
t558 = Icges(3,2) + Icges(4,3);
t555 = Icges(4,6) - Icges(3,6);
t553 = -Icges(5,3) - Icges(3,3) - Icges(4,2);
t499 = sin(pkin(11));
t501 = sin(pkin(6));
t506 = sin(qJ(2));
t508 = cos(qJ(2));
t536 = cos(pkin(11));
t484 = (t499 * t506 + t508 * t536) * t501;
t485 = (-t499 * t508 + t506 * t536) * t501;
t503 = cos(pkin(6));
t552 = (Icges(5,5) * t485 - Icges(5,6) * t484 - (t506 * t556 - t508 * t555) * t501 + t553 * t503) * t503;
t500 = sin(pkin(10));
t502 = cos(pkin(10));
t531 = t503 * t508;
t488 = t500 * t506 - t502 * t531;
t532 = t503 * t506;
t489 = t500 * t508 + t502 * t532;
t534 = t501 * t502;
t551 = t488 * t558 - t489 * t557 - t534 * t555;
t490 = t500 * t531 + t502 * t506;
t491 = -t500 * t532 + t502 * t508;
t535 = t500 * t501;
t548 = t490 * t558 - t491 * t557 + t535 * t555;
t547 = t557 * t488 - t489 * t559 + t556 * t534;
t546 = -t557 * t490 + t491 * t559 + t556 * t535;
t461 = -t490 * t536 + t491 * t499;
t462 = t490 * t499 + t491 * t536;
t541 = Icges(5,5) * t462 - Icges(5,6) * t461 - t490 * t555 - t556 * t491 + t535 * t553;
t459 = -t488 * t536 + t489 * t499;
t460 = t488 * t499 + t489 * t536;
t542 = Icges(5,5) * t460 - Icges(5,6) * t459 - t488 * t555 - t489 * t556 - t534 * t553;
t550 = -t542 * t534 + t541 * t535 + t552;
t545 = t555 * t503 + (-t557 * t506 - t508 * t558) * t501;
t543 = t556 * t503 + (t506 * t559 + t557 * t508) * t501;
t538 = qJD(2) ^ 2;
t537 = cos(qJ(5));
t505 = sin(qJ(5));
t533 = t501 * t505;
t464 = pkin(2) * t491 + qJ(3) * t490;
t498 = qJD(2) * t503;
t530 = qJD(3) * t488 + t464 * t498;
t463 = pkin(2) * t489 + qJ(3) * t488;
t471 = pkin(3) * t489 + qJ(4) * t534;
t529 = -t463 - t471;
t492 = (pkin(2) * t506 - qJ(3) * t508) * t501;
t528 = -pkin(3) * t501 * t506 + qJ(4) * t503 - t492;
t527 = qJD(2) * t501;
t497 = t500 * t527;
t431 = qJD(5) * t461 + t497;
t470 = qJD(5) * t484 + t498;
t526 = qJD(3) * t508;
t525 = qJD(4) * t501;
t523 = t501 * t537;
t522 = t502 * t527;
t521 = t463 * t497 + t464 * t522 + qJD(1);
t519 = (-t503 * rSges(4,2) - (rSges(4,1) * t506 - rSges(4,3) * t508) * t501 - t492) * t501;
t472 = pkin(3) * t491 - qJ(4) * t535;
t518 = t472 * t498 + t502 * t525 + t530;
t515 = (-rSges(5,1) * t485 + rSges(5,2) * t484 + rSges(5,3) * t503 + t528) * t501;
t514 = (-pkin(4) * t485 - pkin(8) * t484 + t528) * t501;
t432 = qJD(5) * t459 - t522;
t487 = qJD(3) * t490;
t513 = -t500 * t525 + t487;
t512 = -qJD(4) * t503 + t471 * t497 + t472 * t522 + t521;
t416 = pkin(4) * t462 + pkin(8) * t461;
t511 = qJD(2) * t500 * t514 + t416 * t498 + t518;
t415 = pkin(4) * t460 + pkin(8) * t459;
t510 = t415 * t497 + t416 * t522 - t501 * t526 + t512;
t509 = ((-t415 + t529) * t503 + t502 * t514) * qJD(2) + t513;
t507 = cos(qJ(6));
t504 = sin(qJ(6));
t481 = t503 * rSges(3,3) + (rSges(3,1) * t506 + rSges(3,2) * t508) * t501;
t469 = t485 * t537 - t503 * t505;
t468 = t485 * t505 + t503 * t537;
t454 = rSges(3,1) * t491 - rSges(3,2) * t490 + rSges(3,3) * t535;
t453 = rSges(4,1) * t491 + rSges(4,2) * t535 + rSges(4,3) * t490;
t452 = rSges(3,1) * t489 - rSges(3,2) * t488 - rSges(3,3) * t534;
t451 = rSges(4,1) * t489 - rSges(4,2) * t534 + rSges(4,3) * t488;
t435 = Icges(5,1) * t485 - Icges(5,4) * t484 - Icges(5,5) * t503;
t434 = Icges(5,4) * t485 - Icges(5,2) * t484 - Icges(5,6) * t503;
t430 = t462 * t537 - t500 * t533;
t429 = t462 * t505 + t500 * t523;
t428 = t460 * t537 + t502 * t533;
t427 = t460 * t505 - t502 * t523;
t426 = t469 * t507 + t484 * t504;
t425 = -t469 * t504 + t484 * t507;
t424 = qJD(6) * t468 + t470;
t423 = pkin(5) * t469 + pkin(9) * t468;
t422 = (-t452 * t503 - t481 * t534) * qJD(2);
t421 = (t454 * t503 - t481 * t535) * qJD(2);
t420 = rSges(6,1) * t469 - rSges(6,2) * t468 + rSges(6,3) * t484;
t419 = Icges(6,1) * t469 - Icges(6,4) * t468 + Icges(6,5) * t484;
t418 = Icges(6,4) * t469 - Icges(6,2) * t468 + Icges(6,6) * t484;
t417 = Icges(6,5) * t469 - Icges(6,6) * t468 + Icges(6,3) * t484;
t411 = rSges(5,1) * t462 - rSges(5,2) * t461 - rSges(5,3) * t535;
t410 = rSges(5,1) * t460 - rSges(5,2) * t459 + rSges(5,3) * t534;
t409 = Icges(5,1) * t462 - Icges(5,4) * t461 - Icges(5,5) * t535;
t408 = Icges(5,1) * t460 - Icges(5,4) * t459 + Icges(5,5) * t534;
t407 = Icges(5,4) * t462 - Icges(5,2) * t461 - Icges(5,6) * t535;
t406 = Icges(5,4) * t460 - Icges(5,2) * t459 + Icges(5,6) * t534;
t403 = t430 * t507 + t461 * t504;
t402 = -t430 * t504 + t461 * t507;
t401 = t428 * t507 + t459 * t504;
t400 = -t428 * t504 + t459 * t507;
t399 = qJD(1) + (t452 * t500 + t454 * t502) * t527;
t398 = qJD(6) * t427 + t432;
t397 = qJD(6) * t429 + t431;
t396 = pkin(5) * t430 + pkin(9) * t429;
t395 = pkin(5) * t428 + pkin(9) * t427;
t394 = rSges(7,1) * t426 + rSges(7,2) * t425 + rSges(7,3) * t468;
t393 = Icges(7,1) * t426 + Icges(7,4) * t425 + Icges(7,5) * t468;
t392 = Icges(7,4) * t426 + Icges(7,2) * t425 + Icges(7,6) * t468;
t391 = Icges(7,5) * t426 + Icges(7,6) * t425 + Icges(7,3) * t468;
t390 = rSges(6,1) * t430 - rSges(6,2) * t429 + rSges(6,3) * t461;
t389 = rSges(6,1) * t428 - rSges(6,2) * t427 + rSges(6,3) * t459;
t388 = Icges(6,1) * t430 - Icges(6,4) * t429 + Icges(6,5) * t461;
t387 = Icges(6,1) * t428 - Icges(6,4) * t427 + Icges(6,5) * t459;
t386 = Icges(6,4) * t430 - Icges(6,2) * t429 + Icges(6,6) * t461;
t385 = Icges(6,4) * t428 - Icges(6,2) * t427 + Icges(6,6) * t459;
t384 = Icges(6,5) * t430 - Icges(6,6) * t429 + Icges(6,3) * t461;
t383 = Icges(6,5) * t428 - Icges(6,6) * t427 + Icges(6,3) * t459;
t382 = t487 + ((-t451 - t463) * t503 + t502 * t519) * qJD(2);
t381 = (t453 * t503 + t500 * t519) * qJD(2) + t530;
t380 = (-t526 + (t451 * t500 + t453 * t502) * qJD(2)) * t501 + t521;
t379 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t429;
t378 = rSges(7,1) * t401 + rSges(7,2) * t400 + rSges(7,3) * t427;
t377 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t429;
t376 = Icges(7,1) * t401 + Icges(7,4) * t400 + Icges(7,5) * t427;
t375 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t429;
t374 = Icges(7,4) * t401 + Icges(7,2) * t400 + Icges(7,6) * t427;
t373 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t429;
t372 = Icges(7,5) * t401 + Icges(7,6) * t400 + Icges(7,3) * t427;
t371 = ((-t410 + t529) * t503 + t502 * t515) * qJD(2) + t513;
t370 = (t411 * t503 + t500 * t515) * qJD(2) + t518;
t369 = (-t526 + (t410 * t500 + t411 * t502) * qJD(2)) * t501 + t512;
t368 = -t389 * t470 + t420 * t432 + t509;
t367 = t390 * t470 - t420 * t431 + t511;
t366 = t431 * t389 - t432 * t390 + t510;
t365 = -t378 * t424 + t394 * t398 - t395 * t470 + t423 * t432 + t509;
t364 = t379 * t424 - t394 * t397 + t396 * t470 - t423 * t431 + t511;
t363 = t397 * t378 - t398 * t379 + t431 * t395 - t432 * t396 + t510;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + t398 * ((t373 * t427 + t375 * t400 + t377 * t401) * t397 + (t427 * t372 + t400 * t374 + t401 * t376) * t398 + (t391 * t427 + t392 * t400 + t393 * t401) * t424) / 0.2e1 + t424 * ((t373 * t468 + t375 * t425 + t377 * t426) * t397 + (t372 * t468 + t374 * t425 + t376 * t426) * t398 + (t468 * t391 + t425 * t392 + t426 * t393) * t424) / 0.2e1 + t397 * ((t429 * t373 + t402 * t375 + t403 * t377) * t397 + (t372 * t429 + t374 * t402 + t376 * t403) * t398 + (t391 * t429 + t392 * t402 + t393 * t403) * t424) / 0.2e1 + t431 * ((t461 * t384 - t429 * t386 + t430 * t388) * t431 + (t383 * t461 - t385 * t429 + t387 * t430) * t432 + (t417 * t461 - t418 * t429 + t419 * t430) * t470) / 0.2e1 + t470 * ((t384 * t484 - t386 * t468 + t388 * t469) * t431 + (t383 * t484 - t385 * t468 + t387 * t469) * t432 + (t417 * t484 - t418 * t468 + t419 * t469) * t470) / 0.2e1 + t432 * ((t384 * t459 - t386 * t427 + t388 * t428) * t431 + (t459 * t383 - t427 * t385 + t428 * t387) * t432 + (t417 * t459 - t418 * t427 + t419 * t428) * t470) / 0.2e1 + m(7) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(6) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(4) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(5) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(3) * (t399 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 - ((-t407 * t459 + t409 * t460 + t488 * t548 + t489 * t546) * t535 + (-t434 * t459 + t435 * t460 + t488 * t545 + t489 * t543) * t503 + (t406 * t459 - t408 * t460 - t488 * t551 + t547 * t489 + t550) * t534) * t538 * t534 / 0.2e1 + ((((t406 * t484 - t408 * t485 + (t506 * t547 + t508 * t551) * t501) * t502 + (-t407 * t484 + t409 * t485 + (t506 * t546 - t508 * t548) * t501) * t500) * t501 + (-t434 * t484 + t485 * t435 + (-t500 * t541 + t502 * t542 + t506 * t543 - t508 * t545) * t501 - t552) * t503) * t503 + ((t406 * t461 - t408 * t462 - t490 * t551 + t547 * t491) * t534 + (-t434 * t461 + t435 * t462 + t490 * t545 + t491 * t543) * t503 + (-t407 * t461 + t409 * t462 + t490 * t548 + t491 * t546 - t550) * t535) * t535) * t538 / 0.2e1;
T  = t1;
