% Calculate kinetic energy for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:43
% EndTime: 2019-03-09 12:16:47
% DurationCPUTime: 3.18s
% Computational Cost: add. (1214->250), mult. (2949->379), div. (0->0), fcn. (3271->8), ass. (0->135)
t562 = Icges(3,4) - Icges(4,5);
t561 = Icges(3,1) + Icges(4,1);
t560 = Icges(3,2) + Icges(4,3);
t471 = cos(qJ(2));
t559 = t562 * t471;
t468 = sin(qJ(2));
t558 = t562 * t468;
t557 = Icges(4,4) + Icges(3,5);
t556 = Icges(3,6) - Icges(4,6);
t555 = t560 * t468 - t559;
t554 = t561 * t471 - t558;
t553 = Icges(6,1) + Icges(7,1);
t552 = -Icges(6,4) + Icges(7,5);
t551 = Icges(7,4) + Icges(6,5);
t550 = Icges(4,2) + Icges(3,3);
t549 = Icges(6,2) + Icges(7,3);
t548 = Icges(7,6) - Icges(6,6);
t547 = Icges(6,3) + Icges(7,2);
t469 = sin(qJ(1));
t472 = cos(qJ(1));
t546 = t555 * t469 + t556 * t472;
t545 = -t556 * t469 + t555 * t472;
t544 = -t554 * t469 + t557 * t472;
t543 = t557 * t469 + t554 * t472;
t542 = -t560 * t471 - t558;
t541 = t561 * t468 + t559;
t540 = -t556 * t468 + t557 * t471;
t539 = rSges(7,1) + pkin(5);
t538 = rSges(7,3) + qJ(6);
t467 = sin(qJ(4));
t516 = cos(qJ(4));
t440 = t468 * t467 + t471 * t516;
t430 = t440 * t469;
t466 = sin(qJ(5));
t470 = cos(qJ(5));
t407 = t430 * t466 - t472 * t470;
t408 = t430 * t470 + t466 * t472;
t500 = t468 * t516;
t441 = -t471 * t467 + t500;
t429 = t441 * t469;
t537 = t549 * t407 + t552 * t408 - t548 * t429;
t432 = t440 * t472;
t409 = t432 * t466 + t469 * t470;
t410 = t432 * t470 - t466 * t469;
t509 = t471 * t472;
t431 = t467 * t509 - t472 * t500;
t536 = t549 * t409 + t552 * t410 + t548 * t431;
t535 = t548 * t407 + t551 * t408 - t547 * t429;
t534 = t548 * t409 + t551 * t410 + t547 * t431;
t533 = t552 * t407 + t553 * t408 - t551 * t429;
t532 = t552 * t409 + t553 * t410 + t551 * t431;
t531 = (t549 * t466 + t552 * t470) * t441 + t548 * t440;
t530 = (t548 * t466 + t551 * t470) * t441 + t547 * t440;
t529 = (t552 * t466 + t553 * t470) * t441 + t551 * t440;
t528 = t540 * t469 - t550 * t472;
t527 = t550 * t469 + t540 * t472;
t526 = t557 * t468 + t556 * t471;
t525 = t542 * t468 + t541 * t471;
t524 = t545 * t468 + t543 * t471;
t523 = -t546 * t468 + t544 * t471;
t508 = -rSges(7,2) * t429 + t538 * t407 + t408 * t539;
t507 = rSges(7,2) * t431 + t538 * t409 + t410 * t539;
t506 = rSges(7,2) * t440 + (t538 * t466 + t470 * t539) * t441;
t493 = pkin(2) * t471 + qJ(3) * t468;
t437 = t493 * t469;
t458 = pkin(1) * t469 - pkin(7) * t472;
t505 = -t437 - t458;
t465 = qJD(2) * t469;
t504 = qJD(2) * t472;
t503 = qJD(3) * t468;
t438 = t493 * t472;
t444 = qJD(1) * (pkin(1) * t472 + pkin(7) * t469);
t502 = qJD(1) * t438 + t469 * t503 + t444;
t442 = pkin(3) * t469 * t471 + pkin(8) * t472;
t501 = -t442 + t505;
t453 = pkin(2) * t468 - qJ(3) * t471;
t497 = qJD(2) * (-rSges(4,1) * t468 + rSges(4,3) * t471 - t453);
t446 = qJD(4) * t472 - t504;
t445 = -qJD(4) * t469 + t465;
t496 = -qJD(3) * t471 + t437 * t465 + t438 * t504;
t495 = rSges(3,1) * t471 - rSges(3,2) * t468;
t494 = rSges(4,1) * t471 + rSges(4,3) * t468;
t492 = qJD(2) * (-pkin(3) * t468 - t453);
t443 = pkin(3) * t509 - pkin(8) * t469;
t479 = t442 * t465 + t443 * t504 + t496;
t462 = t472 * t503;
t478 = t472 * t492 + t462;
t397 = pkin(4) * t430 - pkin(9) * t429;
t398 = pkin(4) * t432 + pkin(9) * t431;
t477 = t445 * t397 - t398 * t446 + t479;
t476 = qJD(1) * t443 + t469 * t492 + t502;
t406 = pkin(4) * t441 + pkin(9) * t440;
t475 = qJD(1) * t398 - t406 * t445 + t476;
t474 = t446 * t406 + (-t397 + t501) * qJD(1) + t478;
t457 = rSges(2,1) * t472 - rSges(2,2) * t469;
t456 = rSges(2,1) * t469 + rSges(2,2) * t472;
t455 = rSges(3,1) * t468 + rSges(3,2) * t471;
t433 = qJD(5) * t440 + qJD(1);
t428 = rSges(3,3) * t469 + t472 * t495;
t427 = rSges(4,2) * t469 + t472 * t494;
t426 = -rSges(3,3) * t472 + t469 * t495;
t425 = -rSges(4,2) * t472 + t469 * t494;
t405 = rSges(5,1) * t441 - rSges(5,2) * t440;
t404 = Icges(5,1) * t441 - Icges(5,4) * t440;
t403 = Icges(5,4) * t441 - Icges(5,2) * t440;
t402 = Icges(5,5) * t441 - Icges(5,6) * t440;
t401 = -qJD(5) * t429 + t446;
t400 = qJD(5) * t431 + t445;
t394 = rSges(5,1) * t432 - rSges(5,2) * t431 - rSges(5,3) * t469;
t393 = rSges(5,1) * t430 + rSges(5,2) * t429 + rSges(5,3) * t472;
t392 = Icges(5,1) * t432 - Icges(5,4) * t431 - Icges(5,5) * t469;
t391 = Icges(5,1) * t430 + Icges(5,4) * t429 + Icges(5,5) * t472;
t390 = Icges(5,4) * t432 - Icges(5,2) * t431 - Icges(5,6) * t469;
t389 = Icges(5,4) * t430 + Icges(5,2) * t429 + Icges(5,6) * t472;
t388 = Icges(5,5) * t432 - Icges(5,6) * t431 - Icges(5,3) * t469;
t387 = Icges(5,5) * t430 + Icges(5,6) * t429 + Icges(5,3) * t472;
t386 = qJD(1) * t428 - t455 * t465 + t444;
t385 = -t455 * t504 + (-t426 - t458) * qJD(1);
t384 = (t426 * t469 + t428 * t472) * qJD(2);
t383 = rSges(6,3) * t440 + (rSges(6,1) * t470 - rSges(6,2) * t466) * t441;
t372 = rSges(6,1) * t410 - rSges(6,2) * t409 + rSges(6,3) * t431;
t370 = rSges(6,1) * t408 - rSges(6,2) * t407 - rSges(6,3) * t429;
t356 = qJD(1) * t427 + t469 * t497 + t502;
t355 = t462 + t472 * t497 + (-t425 + t505) * qJD(1);
t354 = (t425 * t469 + t427 * t472) * qJD(2) + t496;
t353 = qJD(1) * t394 - t405 * t445 + t476;
t352 = t405 * t446 + (-t393 + t501) * qJD(1) + t478;
t351 = t393 * t445 - t394 * t446 + t479;
t350 = t372 * t433 - t383 * t400 + t475;
t349 = -t370 * t433 + t383 * t401 + t474;
t348 = t370 * t400 - t372 * t401 + t477;
t347 = qJD(6) * t407 - t400 * t506 + t433 * t507 + t475;
t346 = qJD(6) * t409 + t401 * t506 - t433 * t508 + t474;
t345 = qJD(6) * t441 * t466 + t400 * t508 - t401 * t507 + t477;
t1 = m(7) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(5) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(4) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + t445 * ((-t469 * t388 - t431 * t390 + t432 * t392) * t445 + (-t387 * t469 - t389 * t431 + t391 * t432) * t446 + (-t402 * t469 - t403 * t431 + t404 * t432) * qJD(1)) / 0.2e1 + t446 * ((t388 * t472 + t390 * t429 + t392 * t430) * t445 + (t472 * t387 + t429 * t389 + t430 * t391) * t446 + (t402 * t472 + t403 * t429 + t404 * t430) * qJD(1)) / 0.2e1 + m(6) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + ((t531 * t409 + t529 * t410 + t530 * t431) * t433 + (t537 * t409 + t533 * t410 + t535 * t431) * t401 + (t536 * t409 + t532 * t410 + t534 * t431) * t400) * t400 / 0.2e1 + ((t531 * t407 + t529 * t408 - t530 * t429) * t433 + (t537 * t407 + t533 * t408 - t535 * t429) * t401 + (t536 * t407 + t532 * t408 - t534 * t429) * t400) * t401 / 0.2e1 + (((t531 * t466 + t529 * t470) * t433 + (t537 * t466 + t533 * t470) * t401 + (t536 * t466 + t532 * t470) * t400) * t441 + (t534 * t400 + t535 * t401 + t530 * t433) * t440) * t433 / 0.2e1 + (Icges(2,3) + m(2) * (t456 ^ 2 + t457 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t527 * t469 ^ 2 + (t523 * t472 + (t524 - t528) * t469) * t472) * qJD(2) + (t526 * t469 + t525 * t472) * qJD(1)) * t465 / 0.2e1 - ((t528 * t472 ^ 2 + (t524 * t469 + (t523 - t527) * t472) * t469) * qJD(2) + (t525 * t469 - t526 * t472) * qJD(1)) * t504 / 0.2e1 + ((-t390 * t440 + t392 * t441) * t445 + (-t389 * t440 + t391 * t441) * t446 + ((t544 * t468 + t546 * t471) * t472 + (t543 * t468 - t545 * t471) * t469) * qJD(2) + (-t440 * t403 + t441 * t404 + t541 * t468 - t542 * t471) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
