% Calculate kinetic energy for
% S6RRPRRP11
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:23
% EndTime: 2019-03-09 12:45:26
% DurationCPUTime: 3.62s
% Computational Cost: add. (1286->268), mult. (2208->412), div. (0->0), fcn. (2165->8), ass. (0->147)
t593 = Icges(3,4) + Icges(4,6);
t592 = Icges(3,1) + Icges(4,2);
t591 = -Icges(3,2) - Icges(4,3);
t487 = cos(qJ(2));
t590 = t593 * t487;
t484 = sin(qJ(2));
t589 = t593 * t484;
t588 = -Icges(4,4) + Icges(3,5);
t587 = Icges(4,5) - Icges(3,6);
t586 = t591 * t484 + t590;
t585 = -t592 * t487 + t589;
t584 = Icges(4,1) + Icges(3,3);
t583 = Icges(6,1) + Icges(7,1);
t582 = Icges(6,4) + Icges(7,4);
t581 = Icges(7,5) + Icges(6,5);
t580 = Icges(6,2) + Icges(7,2);
t579 = Icges(7,6) + Icges(6,6);
t578 = Icges(7,3) + Icges(6,3);
t485 = sin(qJ(1));
t488 = cos(qJ(1));
t577 = t586 * t485 + t587 * t488;
t576 = -t587 * t485 + t586 * t488;
t575 = t585 * t485 + t588 * t488;
t574 = t588 * t485 - t585 * t488;
t573 = t591 * t487 - t589;
t572 = t592 * t484 + t590;
t571 = t587 * t484 + t588 * t487;
t482 = qJ(4) + qJ(5);
t478 = sin(t482);
t479 = cos(t482);
t536 = t484 * t488;
t435 = -t478 * t485 + t479 * t536;
t436 = t478 * t536 + t479 * t485;
t532 = t487 * t488;
t570 = t579 * t435 + t581 * t436 + t578 * t532;
t537 = t484 * t485;
t437 = t478 * t488 + t479 * t537;
t438 = t478 * t537 - t479 * t488;
t534 = t485 * t487;
t569 = t579 * t437 + t581 * t438 + t578 * t534;
t568 = t580 * t435 + t582 * t436 + t579 * t532;
t567 = t580 * t437 + t582 * t438 + t579 * t534;
t566 = t582 * t435 + t583 * t436 + t581 * t532;
t565 = t582 * t437 + t583 * t438 + t581 * t534;
t483 = sin(qJ(4));
t547 = pkin(4) * t483;
t496 = pkin(9) * t487 + t484 * t547;
t486 = cos(qJ(4));
t545 = t486 * pkin(4);
t399 = t485 * t545 + t488 * t496;
t439 = pkin(9) * t484 - t487 * t547;
t477 = qJD(2) * t485;
t523 = qJD(4) * t487;
t450 = t488 * t523 + t477;
t474 = qJD(4) * t484 + qJD(1);
t564 = t474 * t399 - t439 * t450;
t400 = t485 * t496 - t488 * t545;
t525 = qJD(2) * t488;
t451 = t485 * t523 - t525;
t563 = -t400 * t474 + t451 * t439;
t562 = (-t581 * t478 - t579 * t479) * t487 + t578 * t484;
t561 = (-t582 * t478 - t580 * t479) * t487 + t579 * t484;
t560 = (-t583 * t478 - t582 * t479) * t487 + t581 * t484;
t559 = t584 * t485 + t571 * t488;
t558 = t571 * t485 - t584 * t488;
t557 = t588 * t484 - t587 * t487;
t556 = t573 * t484 + t572 * t487;
t555 = -t576 * t484 + t574 * t487;
t554 = t577 * t484 + t575 * t487;
t535 = t485 * t486;
t533 = t486 * t488;
t518 = pkin(5) * t478;
t494 = qJ(6) * t487 + t484 * t518;
t527 = pkin(5) * t479;
t531 = rSges(7,1) * t436 + rSges(7,2) * t435 + rSges(7,3) * t532 + t485 * t527 + t488 * t494;
t530 = rSges(7,1) * t438 + rSges(7,2) * t437 + rSges(7,3) * t534 + t485 * t494 - t488 * t527;
t529 = (-rSges(7,1) * t478 - rSges(7,2) * t479 - t518) * t487 + (qJ(6) + rSges(7,3)) * t484;
t512 = pkin(2) * t487 + qJ(3) * t484;
t447 = t512 * t485;
t469 = pkin(1) * t485 - pkin(7) * t488;
t528 = -t447 - t469;
t524 = qJD(3) * t484;
t522 = qJD(5) * t487;
t448 = t512 * t488;
t456 = qJD(1) * (pkin(1) * t488 + pkin(7) * t485);
t521 = qJD(1) * t448 + t485 * t524 + t456;
t464 = pkin(2) * t484 - qJ(3) * t487;
t517 = qJD(2) * (rSges(4,2) * t484 + rSges(4,3) * t487 - t464);
t453 = pkin(3) * t485 + pkin(8) * t532;
t516 = qJD(1) * t453 + t521;
t515 = -qJD(3) * t487 + t447 * t477 + t448 * t525;
t514 = rSges(3,1) * t487 - rSges(3,2) * t484;
t513 = -rSges(4,2) * t487 + rSges(4,3) * t484;
t511 = (-pkin(8) * t484 - t464) * qJD(2);
t454 = -pkin(3) * t488 + pkin(8) * t534;
t473 = t488 * t524;
t498 = t473 + (-t454 + t528) * qJD(1);
t497 = t453 * t525 + t454 * t477 + t515;
t495 = qJD(6) * t487 + t511;
t493 = -t399 * t451 + t450 * t400 + t497;
t492 = t485 * t511 + t516;
t491 = t488 * t511 + t498;
t468 = rSges(2,1) * t488 - rSges(2,2) * t485;
t467 = rSges(2,1) * t485 + rSges(2,2) * t488;
t466 = rSges(3,1) * t484 + rSges(3,2) * t487;
t455 = qJD(5) * t484 + t474;
t446 = t483 * t537 - t533;
t445 = t483 * t488 + t484 * t535;
t444 = t483 * t536 + t535;
t443 = -t483 * t485 + t484 * t533;
t434 = -rSges(4,1) * t488 + t485 * t513;
t433 = rSges(4,1) * t485 + t488 * t513;
t432 = rSges(3,3) * t485 + t488 * t514;
t431 = rSges(5,3) * t484 + (-rSges(5,1) * t483 - rSges(5,2) * t486) * t487;
t430 = -rSges(3,3) * t488 + t485 * t514;
t419 = Icges(5,5) * t484 + (-Icges(5,1) * t483 - Icges(5,4) * t486) * t487;
t416 = Icges(5,6) * t484 + (-Icges(5,4) * t483 - Icges(5,2) * t486) * t487;
t413 = Icges(5,3) * t484 + (-Icges(5,5) * t483 - Icges(5,6) * t486) * t487;
t412 = t485 * t522 + t451;
t411 = t488 * t522 + t450;
t410 = rSges(6,3) * t484 + (-rSges(6,1) * t478 - rSges(6,2) * t479) * t487;
t398 = rSges(5,1) * t446 + rSges(5,2) * t445 + rSges(5,3) * t534;
t397 = rSges(5,1) * t444 + rSges(5,2) * t443 + rSges(5,3) * t532;
t396 = Icges(5,1) * t446 + Icges(5,4) * t445 + Icges(5,5) * t534;
t395 = Icges(5,1) * t444 + Icges(5,4) * t443 + Icges(5,5) * t532;
t394 = Icges(5,4) * t446 + Icges(5,2) * t445 + Icges(5,6) * t534;
t393 = Icges(5,4) * t444 + Icges(5,2) * t443 + Icges(5,6) * t532;
t392 = Icges(5,5) * t446 + Icges(5,6) * t445 + Icges(5,3) * t534;
t391 = Icges(5,5) * t444 + Icges(5,6) * t443 + Icges(5,3) * t532;
t389 = qJD(1) * t432 - t466 * t477 + t456;
t388 = -t466 * t525 + (-t430 - t469) * qJD(1);
t387 = (t430 * t485 + t432 * t488) * qJD(2);
t386 = rSges(6,1) * t438 + rSges(6,2) * t437 + rSges(6,3) * t534;
t384 = rSges(6,1) * t436 + rSges(6,2) * t435 + rSges(6,3) * t532;
t367 = qJD(1) * t433 + t485 * t517 + t521;
t366 = t473 + t488 * t517 + (-t434 + t528) * qJD(1);
t365 = (t433 * t488 + t434 * t485) * qJD(2) + t515;
t364 = t397 * t474 - t431 * t450 + t492;
t363 = -t398 * t474 + t431 * t451 + t491;
t362 = -t397 * t451 + t398 * t450 + t497;
t361 = t384 * t455 - t410 * t411 + t492 + t564;
t360 = -t386 * t455 + t410 * t412 + t491 + t563;
t359 = -t384 * t412 + t386 * t411 + t493;
t358 = -t411 * t529 + t455 * t531 + t485 * t495 + t516 + t564;
t357 = t412 * t529 - t455 * t530 + t488 * t495 + t498 + t563;
t356 = qJD(6) * t484 + t411 * t530 - t412 * t531 + t493;
t1 = t450 * ((t391 * t532 + t443 * t393 + t444 * t395) * t450 + (t392 * t532 + t394 * t443 + t396 * t444) * t451 + (t413 * t532 + t416 * t443 + t419 * t444) * t474) / 0.2e1 + t451 * ((t391 * t534 + t393 * t445 + t395 * t446) * t450 + (t392 * t534 + t445 * t394 + t446 * t396) * t451 + (t413 * t534 + t416 * t445 + t419 * t446) * t474) / 0.2e1 + t474 * ((t391 * t450 + t392 * t451 + t413 * t474) * t484 + ((-t393 * t486 - t395 * t483) * t450 + (-t394 * t486 - t396 * t483) * t451 + (-t416 * t486 - t419 * t483) * t474) * t487) / 0.2e1 + m(6) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(7) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(5) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(4) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(3) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + ((t561 * t435 + t560 * t436 + t562 * t532) * t455 + (t567 * t435 + t565 * t436 + t569 * t532) * t412 + (t568 * t435 + t566 * t436 + t570 * t532) * t411) * t411 / 0.2e1 + ((t561 * t437 + t560 * t438 + t562 * t534) * t455 + (t567 * t437 + t565 * t438 + t569 * t534) * t412 + (t568 * t437 + t566 * t438 + t570 * t534) * t411) * t412 / 0.2e1 + (((-t560 * t478 - t561 * t479) * t455 + (-t565 * t478 - t567 * t479) * t412 + (-t566 * t478 - t568 * t479) * t411) * t487 + (t570 * t411 + t569 * t412 + t562 * t455) * t484) * t455 / 0.2e1 + (m(2) * (t467 ^ 2 + t468 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t575 * t484 - t577 * t487) * t488 + (t574 * t484 + t576 * t487) * t485) * qJD(2) + (t572 * t484 - t573 * t487) * qJD(1)) * qJD(1) / 0.2e1 + ((t559 * t485 ^ 2 + (t554 * t488 + (t555 - t558) * t485) * t488) * qJD(2) + (t557 * t485 + t556 * t488) * qJD(1)) * t477 / 0.2e1 - ((t558 * t488 ^ 2 + (t555 * t485 + (t554 - t559) * t488) * t485) * qJD(2) + (t556 * t485 - t557 * t488) * qJD(1)) * t525 / 0.2e1;
T  = t1;
