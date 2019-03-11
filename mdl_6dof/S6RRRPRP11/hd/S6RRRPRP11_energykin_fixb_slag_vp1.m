% Calculate kinetic energy for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:32
% EndTime: 2019-03-09 17:41:35
% DurationCPUTime: 2.76s
% Computational Cost: add. (2293->291), mult. (5705->428), div. (0->0), fcn. (6906->10), ass. (0->141)
t584 = Icges(4,1) + Icges(5,2);
t583 = Icges(5,1) + Icges(4,3);
t582 = Icges(6,1) + Icges(7,1);
t581 = -Icges(4,4) - Icges(5,6);
t580 = Icges(5,4) - Icges(4,5);
t579 = Icges(6,4) + Icges(7,4);
t578 = Icges(5,5) - Icges(4,6);
t577 = Icges(6,5) + Icges(7,5);
t576 = Icges(4,2) + Icges(5,3);
t575 = Icges(6,2) + Icges(7,2);
t574 = Icges(6,6) + Icges(7,6);
t573 = Icges(6,3) + Icges(7,3);
t572 = rSges(7,3) + qJ(6);
t509 = sin(qJ(2));
t510 = sin(qJ(1));
t512 = cos(qJ(2));
t513 = cos(qJ(1));
t541 = cos(pkin(6));
t524 = t513 * t541;
t488 = t509 * t510 - t512 * t524;
t489 = t509 * t524 + t510 * t512;
t506 = sin(pkin(6));
t535 = t506 * t513;
t453 = Icges(3,5) * t489 - Icges(3,6) * t488 - Icges(3,3) * t535;
t525 = t510 * t541;
t490 = t513 * t509 + t512 * t525;
t491 = -t509 * t525 + t513 * t512;
t537 = t506 * t510;
t454 = Icges(3,5) * t491 - Icges(3,6) * t490 + Icges(3,3) * t537;
t571 = t506 * (t453 * t513 - t454 * t510);
t545 = cos(qJ(3));
t528 = t506 * t545;
t544 = sin(qJ(3));
t472 = t489 * t544 + t513 * t528;
t508 = sin(qJ(5));
t511 = cos(qJ(5));
t439 = t472 * t511 - t488 * t508;
t540 = t472 * t508;
t440 = t488 * t511 + t540;
t527 = t506 * t544;
t473 = t489 * t545 - t513 * t527;
t570 = t439 * t574 + t440 * t577 + t473 * t573;
t474 = t491 * t544 - t510 * t528;
t441 = t474 * t511 - t490 * t508;
t539 = t474 * t508;
t442 = t490 * t511 + t539;
t475 = t491 * t545 + t510 * t527;
t569 = t441 * t574 + t442 * t577 + t475 * t573;
t568 = t439 * t575 + t440 * t579 + t473 * t574;
t567 = t441 * t575 + t442 * t579 + t475 * t574;
t566 = t439 * t579 + t440 * t582 + t473 * t577;
t565 = t441 * t579 + t442 * t582 + t475 * t577;
t486 = t509 * t527 - t541 * t545;
t536 = t506 * t512;
t470 = t486 * t511 + t508 * t536;
t538 = t486 * t508;
t471 = -t511 * t536 + t538;
t487 = t509 * t528 + t541 * t544;
t564 = t470 * t574 + t471 * t577 + t487 * t573;
t563 = t470 * t575 + t471 * t579 + t487 * t574;
t562 = t470 * t579 + t471 * t582 + t487 * t577;
t561 = t472 * t576 + t473 * t581 + t488 * t578;
t560 = t474 * t576 + t475 * t581 + t490 * t578;
t559 = t472 * t578 - t473 * t580 + t488 * t583;
t558 = t474 * t578 - t475 * t580 + t490 * t583;
t557 = t581 * t472 + t473 * t584 - t580 * t488;
t556 = t581 * t474 + t475 * t584 - t580 * t490;
t555 = t486 * t576 + t487 * t581 - t536 * t578;
t554 = t581 * t486 + t487 * t584 + t580 * t536;
t553 = t486 * t578 - t487 * t580 - t536 * t583;
t543 = pkin(5) * t511;
t534 = rSges(7,1) * t440 + rSges(7,2) * t439 + pkin(5) * t540 + t572 * t473 + t488 * t543;
t533 = rSges(7,1) * t442 + rSges(7,2) * t441 + pkin(5) * t539 + t572 * t475 + t490 * t543;
t532 = rSges(7,1) * t471 + rSges(7,2) * t470 + pkin(5) * t538 + t572 * t487 - t536 * t543;
t465 = pkin(2) * t489 + pkin(9) * t488;
t466 = pkin(2) * t491 + pkin(9) * t490;
t529 = qJD(2) * t506;
t501 = t510 * t529;
t526 = t513 * t529;
t531 = t465 * t501 + t466 * t526;
t476 = qJD(3) * t490 + t501;
t530 = qJD(1) * (pkin(1) * t510 - pkin(8) * t535);
t502 = qJD(2) * t541 + qJD(1);
t434 = pkin(3) * t473 + qJ(4) * t472;
t523 = qJD(4) * t486 + t476 * t434 + t531;
t477 = qJD(3) * t488 - t526;
t493 = -qJD(3) * t536 + t502;
t492 = (pkin(2) * t509 - pkin(9) * t512) * t506;
t494 = qJD(1) * (pkin(1) * t513 + pkin(8) * t537);
t521 = t502 * t466 - t492 * t501 + t494;
t435 = pkin(3) * t475 + qJ(4) * t474;
t520 = qJD(4) * t472 + t493 * t435 + t521;
t443 = pkin(4) * t488 + pkin(10) * t473;
t444 = pkin(4) * t490 + pkin(10) * t475;
t519 = t476 * t443 + (-t435 - t444) * t477 + t523;
t518 = -t465 * t502 - t492 * t526 - t530;
t464 = pkin(3) * t487 + qJ(4) * t486;
t517 = qJD(4) * t474 + t477 * t464 + t518;
t478 = -pkin(4) * t536 + pkin(10) * t487;
t516 = t493 * t444 + (-t464 - t478) * t476 + t520;
t515 = t477 * t478 + (-t434 - t443) * t493 + t517;
t497 = rSges(2,1) * t513 - rSges(2,2) * t510;
t496 = rSges(2,1) * t510 + rSges(2,2) * t513;
t482 = t541 * rSges(3,3) + (rSges(3,1) * t509 + rSges(3,2) * t512) * t506;
t481 = Icges(3,5) * t541 + (Icges(3,1) * t509 + Icges(3,4) * t512) * t506;
t480 = Icges(3,6) * t541 + (Icges(3,4) * t509 + Icges(3,2) * t512) * t506;
t479 = Icges(3,3) * t541 + (Icges(3,5) * t509 + Icges(3,6) * t512) * t506;
t467 = qJD(5) * t487 + t493;
t461 = rSges(3,1) * t491 - rSges(3,2) * t490 + rSges(3,3) * t537;
t460 = rSges(3,1) * t489 - rSges(3,2) * t488 - rSges(3,3) * t535;
t458 = Icges(3,1) * t491 - Icges(3,4) * t490 + Icges(3,5) * t537;
t457 = Icges(3,1) * t489 - Icges(3,4) * t488 - Icges(3,5) * t535;
t456 = Icges(3,4) * t491 - Icges(3,2) * t490 + Icges(3,6) * t537;
t455 = Icges(3,4) * t489 - Icges(3,2) * t488 - Icges(3,6) * t535;
t452 = rSges(4,1) * t487 - rSges(4,2) * t486 - rSges(4,3) * t536;
t451 = -rSges(5,1) * t536 - rSges(5,2) * t487 + rSges(5,3) * t486;
t437 = qJD(5) * t473 + t477;
t436 = qJD(5) * t475 + t476;
t428 = rSges(4,1) * t475 - rSges(4,2) * t474 + rSges(4,3) * t490;
t427 = rSges(4,1) * t473 - rSges(4,2) * t472 + rSges(4,3) * t488;
t426 = rSges(5,1) * t490 - rSges(5,2) * t475 + rSges(5,3) * t474;
t425 = rSges(5,1) * t488 - rSges(5,2) * t473 + rSges(5,3) * t472;
t412 = rSges(6,1) * t471 + rSges(6,2) * t470 + rSges(6,3) * t487;
t403 = t461 * t502 - t482 * t501 + t494;
t402 = -t460 * t502 - t482 * t526 - t530;
t401 = (t460 * t510 + t461 * t513) * t529;
t398 = rSges(6,1) * t442 + rSges(6,2) * t441 + rSges(6,3) * t475;
t396 = rSges(6,1) * t440 + rSges(6,2) * t439 + rSges(6,3) * t473;
t382 = t428 * t493 - t452 * t476 + t521;
t381 = -t427 * t493 + t452 * t477 + t518;
t380 = t427 * t476 - t428 * t477 + t531;
t379 = t426 * t493 + (-t451 - t464) * t476 + t520;
t378 = t451 * t477 + (-t425 - t434) * t493 + t517;
t377 = t425 * t476 + (-t426 - t435) * t477 + t523;
t376 = t398 * t467 - t412 * t436 + t516;
t375 = -t396 * t467 + t412 * t437 + t515;
t374 = t396 * t436 - t398 * t437 + t519;
t373 = qJD(6) * t473 - t436 * t532 + t467 * t533 + t516;
t372 = qJD(6) * t475 + t437 * t532 - t467 * t534 + t515;
t371 = qJD(6) * t487 + t436 * t534 - t437 * t533 + t519;
t1 = m(6) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(7) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + t502 * ((t541 * t454 + (t456 * t512 + t458 * t509) * t506) * t501 - (t541 * t453 + (t455 * t512 + t457 * t509) * t506) * t526 + (t541 * t479 + (t480 * t512 + t481 * t509) * t506) * t502) / 0.2e1 + m(5) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(4) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(3) * (t401 ^ 2 + t402 ^ 2 + t403 ^ 2) / 0.2e1 - ((-t479 * t535 - t480 * t488 + t481 * t489) * t502 + ((-t456 * t488 + t458 * t489) * t510 + (t488 * t455 - t489 * t457 + t571) * t513) * t529) * t526 / 0.2e1 + ((t479 * t537 - t480 * t490 + t481 * t491) * t502 + (-(-t455 * t490 + t457 * t491) * t513 + (-t490 * t456 + t491 * t458 - t571) * t510) * t529) * t501 / 0.2e1 + ((t563 * t441 + t562 * t442 + t564 * t475) * t467 + (t568 * t441 + t566 * t442 + t570 * t475) * t437 + (t567 * t441 + t565 * t442 + t569 * t475) * t436) * t436 / 0.2e1 + ((t563 * t439 + t562 * t440 + t564 * t473) * t467 + (t568 * t439 + t566 * t440 + t570 * t473) * t437 + (t567 * t439 + t565 * t440 + t569 * t473) * t436) * t437 / 0.2e1 + ((t563 * t470 + t562 * t471 + t564 * t487) * t467 + (t568 * t470 + t566 * t471 + t570 * t487) * t437 + (t567 * t470 + t565 * t471 + t569 * t487) * t436) * t467 / 0.2e1 + ((t555 * t474 + t554 * t475 + t553 * t490) * t493 + (t561 * t474 + t557 * t475 + t559 * t490) * t477 + (t560 * t474 + t556 * t475 + t558 * t490) * t476) * t476 / 0.2e1 + ((t555 * t472 + t554 * t473 + t553 * t488) * t493 + (t561 * t472 + t557 * t473 + t559 * t488) * t477 + (t560 * t472 + t556 * t473 + t558 * t488) * t476) * t477 / 0.2e1 + ((t555 * t486 + t554 * t487 - t553 * t536) * t493 + (t561 * t486 + t557 * t487 - t559 * t536) * t477 + (t560 * t486 + t556 * t487 - t558 * t536) * t476) * t493 / 0.2e1 + (Icges(2,3) + m(2) * (t496 ^ 2 + t497 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
