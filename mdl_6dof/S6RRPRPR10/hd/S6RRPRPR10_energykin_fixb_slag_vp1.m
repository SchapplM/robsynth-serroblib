% Calculate kinetic energy for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:22
% EndTime: 2019-03-09 11:04:25
% DurationCPUTime: 3.10s
% Computational Cost: add. (2723->341), mult. (4987->505), div. (0->0), fcn. (5916->12), ass. (0->160)
t579 = Icges(5,1) + Icges(6,2);
t578 = Icges(6,1) + Icges(5,3);
t577 = -Icges(5,4) - Icges(6,6);
t576 = Icges(6,4) - Icges(5,5);
t575 = Icges(6,5) - Icges(5,6);
t574 = Icges(5,2) + Icges(6,3);
t515 = sin(qJ(2));
t516 = sin(qJ(1));
t518 = cos(qJ(2));
t519 = cos(qJ(1));
t554 = cos(pkin(6));
t538 = t519 * t554;
t492 = t515 * t516 - t518 * t538;
t493 = t515 * t538 + t516 * t518;
t511 = sin(pkin(6));
t550 = t511 * t519;
t450 = Icges(3,5) * t493 - Icges(3,6) * t492 - Icges(3,3) * t550;
t539 = t516 * t554;
t494 = t519 * t515 + t518 * t539;
t495 = -t515 * t539 + t519 * t518;
t552 = t511 * t516;
t451 = Icges(3,5) * t495 - Icges(3,6) * t494 + Icges(3,3) * t552;
t573 = t511 * (t450 * t519 - t451 * t516);
t544 = pkin(11) + qJ(4);
t536 = cos(t544);
t529 = t511 * t536;
t535 = sin(t544);
t468 = t493 * t535 + t519 * t529;
t528 = t511 * t535;
t469 = t493 * t536 - t519 * t528;
t572 = t574 * t468 + t577 * t469 + t575 * t492;
t470 = t495 * t535 - t516 * t529;
t471 = t495 * t536 + t516 * t528;
t571 = t574 * t470 + t577 * t471 + t575 * t494;
t570 = t575 * t468 - t576 * t469 + t578 * t492;
t569 = t575 * t470 - t576 * t471 + t578 * t494;
t568 = t577 * t468 + t579 * t469 - t576 * t492;
t567 = t577 * t470 + t579 * t471 - t576 * t494;
t510 = sin(pkin(11));
t512 = cos(pkin(11));
t473 = -t493 * t510 - t512 * t550;
t542 = t510 * t550;
t474 = t493 * t512 - t542;
t417 = Icges(4,5) * t474 + Icges(4,6) * t473 + Icges(4,3) * t492;
t452 = Icges(3,4) * t493 - Icges(3,2) * t492 - Icges(3,6) * t550;
t566 = -t417 + t452;
t475 = -t495 * t510 + t512 * t552;
t543 = t510 * t552;
t476 = t495 * t512 + t543;
t418 = Icges(4,5) * t476 + Icges(4,6) * t475 + Icges(4,3) * t494;
t453 = Icges(3,4) * t495 - Icges(3,2) * t494 + Icges(3,6) * t552;
t565 = t418 - t453;
t483 = t515 * t528 - t536 * t554;
t484 = t515 * t529 + t535 * t554;
t551 = t511 * t518;
t564 = t574 * t483 + t577 * t484 - t575 * t551;
t563 = t577 * t483 + t579 * t484 + t576 * t551;
t562 = t575 * t483 - t576 * t484 - t578 * t551;
t553 = t511 * t515;
t490 = -t510 * t553 + t512 * t554;
t537 = t554 * t510;
t491 = t512 * t553 + t537;
t444 = Icges(4,5) * t491 + Icges(4,6) * t490 - Icges(4,3) * t551;
t481 = Icges(3,6) * t554 + (Icges(3,4) * t515 + Icges(3,2) * t518) * t511;
t561 = t444 - t481;
t555 = pkin(3) * t512;
t462 = pkin(2) * t493 + qJ(3) * t492;
t463 = pkin(2) * t495 + qJ(3) * t494;
t546 = qJD(2) * t511;
t506 = t516 * t546;
t540 = t519 * t546;
t548 = t462 * t506 + t463 * t540;
t477 = qJD(4) * t494 + t506;
t547 = qJD(1) * (pkin(1) * t516 - pkin(8) * t550);
t545 = qJD(3) * t518;
t507 = qJD(2) * t554 + qJD(1);
t498 = qJD(1) * (pkin(1) * t519 + pkin(8) * t552);
t541 = qJD(3) * t492 + t507 * t463 + t498;
t534 = qJD(3) * t494 - t547;
t496 = (pkin(2) * t515 - qJ(3) * t518) * t511;
t531 = (-rSges(4,1) * t491 - rSges(4,2) * t490 + rSges(4,3) * t551 - t496) * t546;
t530 = (-pkin(3) * t537 - (-pkin(9) * t518 + t515 * t555) * t511 - t496) * t546;
t478 = qJD(4) * t492 - t540;
t497 = -qJD(4) * t551 + t507;
t414 = -pkin(3) * t542 + pkin(9) * t492 + t493 * t555;
t415 = pkin(3) * t543 + pkin(9) * t494 + t495 * t555;
t526 = t414 * t506 + t415 * t540 - t511 * t545 + t548;
t426 = pkin(4) * t469 + qJ(5) * t468;
t525 = qJD(5) * t483 + t477 * t426 + t526;
t524 = t507 * t415 + t516 * t530 + t541;
t427 = pkin(4) * t471 + qJ(5) * t470;
t523 = qJD(5) * t468 + t497 * t427 + t524;
t522 = (-t414 - t462) * t507 + t519 * t530 + t534;
t447 = pkin(4) * t484 + qJ(5) * t483;
t521 = qJD(5) * t470 + t478 * t447 + t522;
t517 = cos(qJ(6));
t514 = sin(qJ(6));
t503 = rSges(2,1) * t519 - rSges(2,2) * t516;
t502 = rSges(2,1) * t516 + rSges(2,2) * t519;
t485 = t554 * rSges(3,3) + (rSges(3,1) * t515 + rSges(3,2) * t518) * t511;
t482 = Icges(3,5) * t554 + (Icges(3,1) * t515 + Icges(3,4) * t518) * t511;
t480 = Icges(3,3) * t554 + (Icges(3,5) * t515 + Icges(3,6) * t518) * t511;
t472 = -pkin(5) * t551 + pkin(10) * t484;
t467 = t483 * t514 - t517 * t551;
t466 = t483 * t517 + t514 * t551;
t461 = qJD(6) * t484 + t497;
t460 = rSges(3,1) * t495 - rSges(3,2) * t494 + rSges(3,3) * t552;
t459 = rSges(3,1) * t493 - rSges(3,2) * t492 - rSges(3,3) * t550;
t455 = Icges(3,1) * t495 - Icges(3,4) * t494 + Icges(3,5) * t552;
t454 = Icges(3,1) * t493 - Icges(3,4) * t492 - Icges(3,5) * t550;
t446 = Icges(4,1) * t491 + Icges(4,4) * t490 - Icges(4,5) * t551;
t445 = Icges(4,4) * t491 + Icges(4,2) * t490 - Icges(4,6) * t551;
t443 = pkin(5) * t494 + pkin(10) * t471;
t442 = pkin(5) * t492 + pkin(10) * t469;
t441 = rSges(5,1) * t484 - rSges(5,2) * t483 - rSges(5,3) * t551;
t440 = -rSges(6,1) * t551 - rSges(6,2) * t484 + rSges(6,3) * t483;
t433 = t470 * t514 + t494 * t517;
t432 = t470 * t517 - t494 * t514;
t431 = t468 * t514 + t492 * t517;
t430 = t468 * t517 - t492 * t514;
t429 = qJD(6) * t469 + t478;
t428 = qJD(6) * t471 + t477;
t424 = rSges(4,1) * t476 + rSges(4,2) * t475 + rSges(4,3) * t494;
t423 = rSges(4,1) * t474 + rSges(4,2) * t473 + rSges(4,3) * t492;
t422 = Icges(4,1) * t476 + Icges(4,4) * t475 + Icges(4,5) * t494;
t421 = Icges(4,1) * t474 + Icges(4,4) * t473 + Icges(4,5) * t492;
t420 = Icges(4,4) * t476 + Icges(4,2) * t475 + Icges(4,6) * t494;
t419 = Icges(4,4) * t474 + Icges(4,2) * t473 + Icges(4,6) * t492;
t413 = rSges(5,1) * t471 - rSges(5,2) * t470 + rSges(5,3) * t494;
t412 = rSges(5,1) * t469 - rSges(5,2) * t468 + rSges(5,3) * t492;
t411 = rSges(6,1) * t494 - rSges(6,2) * t471 + rSges(6,3) * t470;
t410 = rSges(6,1) * t492 - rSges(6,2) * t469 + rSges(6,3) * t468;
t397 = t460 * t507 - t485 * t506 + t498;
t396 = -t459 * t507 - t485 * t540 - t547;
t395 = rSges(7,1) * t467 + rSges(7,2) * t466 + rSges(7,3) * t484;
t394 = Icges(7,1) * t467 + Icges(7,4) * t466 + Icges(7,5) * t484;
t393 = Icges(7,4) * t467 + Icges(7,2) * t466 + Icges(7,6) * t484;
t392 = Icges(7,5) * t467 + Icges(7,6) * t466 + Icges(7,3) * t484;
t387 = (t459 * t516 + t460 * t519) * t546;
t386 = rSges(7,1) * t433 + rSges(7,2) * t432 + rSges(7,3) * t471;
t385 = rSges(7,1) * t431 + rSges(7,2) * t430 + rSges(7,3) * t469;
t384 = Icges(7,1) * t433 + Icges(7,4) * t432 + Icges(7,5) * t471;
t383 = Icges(7,1) * t431 + Icges(7,4) * t430 + Icges(7,5) * t469;
t382 = Icges(7,4) * t433 + Icges(7,2) * t432 + Icges(7,6) * t471;
t381 = Icges(7,4) * t431 + Icges(7,2) * t430 + Icges(7,6) * t469;
t380 = Icges(7,5) * t433 + Icges(7,6) * t432 + Icges(7,3) * t471;
t379 = Icges(7,5) * t431 + Icges(7,6) * t430 + Icges(7,3) * t469;
t378 = t424 * t507 + t516 * t531 + t541;
t377 = (-t423 - t462) * t507 + t519 * t531 + t534;
t376 = (-t545 + (t423 * t516 + t424 * t519) * qJD(2)) * t511 + t548;
t375 = t413 * t497 - t441 * t477 + t524;
t374 = -t412 * t497 + t441 * t478 + t522;
t373 = t412 * t477 - t413 * t478 + t526;
t372 = t411 * t497 + (-t440 - t447) * t477 + t523;
t371 = t440 * t478 + (-t410 - t426) * t497 + t521;
t370 = t410 * t477 + (-t411 - t427) * t478 + t525;
t369 = t386 * t461 - t395 * t428 + t443 * t497 + (-t447 - t472) * t477 + t523;
t368 = -t385 * t461 + t395 * t429 + t472 * t478 + (-t426 - t442) * t497 + t521;
t367 = t385 * t428 - t386 * t429 + t442 * t477 + (-t427 - t443) * t478 + t525;
t1 = m(3) * (t387 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 + m(4) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(5) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(6) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(7) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + t428 * ((t380 * t471 + t382 * t432 + t384 * t433) * t428 + (t379 * t471 + t381 * t432 + t383 * t433) * t429 + (t392 * t471 + t393 * t432 + t394 * t433) * t461) / 0.2e1 + t429 * ((t380 * t469 + t382 * t430 + t384 * t431) * t428 + (t379 * t469 + t381 * t430 + t383 * t431) * t429 + (t392 * t469 + t393 * t430 + t394 * t431) * t461) / 0.2e1 + t461 * ((t380 * t484 + t382 * t466 + t384 * t467) * t428 + (t379 * t484 + t381 * t466 + t383 * t467) * t429 + (t392 * t484 + t393 * t466 + t394 * t467) * t461) / 0.2e1 + ((t564 * t470 + t563 * t471 + t562 * t494) * t497 + (t572 * t470 + t568 * t471 + t570 * t494) * t478 + (t571 * t470 + t567 * t471 + t569 * t494) * t477) * t477 / 0.2e1 + ((t564 * t468 + t563 * t469 + t562 * t492) * t497 + (t572 * t468 + t568 * t469 + t570 * t492) * t478 + (t571 * t468 + t567 * t469 + t569 * t492) * t477) * t478 / 0.2e1 + ((t564 * t483 + t563 * t484 - t562 * t551) * t497 + (t572 * t483 + t568 * t484 - t570 * t551) * t478 + (t571 * t483 + t567 * t484 - t569 * t551) * t477) * t497 / 0.2e1 + ((t554 * t451 + (t453 * t518 + t455 * t515) * t511) * t506 - (t554 * t450 + (t452 * t518 + t454 * t515) * t511) * t540 + ((t420 * t490 + t422 * t491) * t516 - (t419 * t490 + t421 * t491) * t519 + (t417 * t519 - t418 * t516) * t551) * t546 + (t554 * t480 + (t481 * t518 + t482 * t515) * t511 - t444 * t551 + t445 * t490 + t446 * t491) * t507) * t507 / 0.2e1 + (Icges(2,3) + m(2) * (t502 ^ 2 + t503 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t419 * t475 - t421 * t476 - t454 * t495 + t566 * t494) * t519 + (t420 * t475 + t422 * t476 + t495 * t455 + t565 * t494 - t573) * t516) * t546 + (t445 * t475 + t446 * t476 + t480 * t552 + t482 * t495 + t561 * t494) * t507) * t506 / 0.2e1 - (((-t419 * t473 - t421 * t474 - t493 * t454 + t566 * t492 + t573) * t519 + (t420 * t473 + t422 * t474 + t455 * t493 + t565 * t492) * t516) * t546 + (t445 * t473 + t446 * t474 - t480 * t550 + t482 * t493 + t561 * t492) * t507) * t540 / 0.2e1;
T  = t1;
