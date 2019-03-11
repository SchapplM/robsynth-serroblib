% Calculate kinetic energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:00
% EndTime: 2019-03-08 18:49:02
% DurationCPUTime: 2.10s
% Computational Cost: add. (4226->275), mult. (11689->418), div. (0->0), fcn. (15179->14), ass. (0->134)
t589 = Icges(5,1) + Icges(6,2);
t588 = Icges(6,1) + Icges(5,3);
t587 = -Icges(5,4) - Icges(6,6);
t586 = -Icges(6,4) + Icges(5,5);
t585 = Icges(6,5) - Icges(5,6);
t584 = Icges(5,2) + Icges(6,3);
t537 = sin(pkin(11));
t539 = cos(pkin(11));
t565 = sin(pkin(12));
t569 = cos(pkin(6));
t557 = t569 * t565;
t567 = cos(pkin(12));
t528 = t537 * t567 + t539 * t557;
t542 = sin(qJ(3));
t559 = t569 * t567;
t549 = t537 * t565 - t539 * t559;
t568 = cos(pkin(7));
t546 = t549 * t568;
t538 = sin(pkin(6));
t566 = sin(pkin(7));
t560 = t538 * t566;
t571 = cos(qJ(3));
t508 = t528 * t571 + (-t539 * t560 - t546) * t542;
t561 = t538 * t568;
t521 = -t539 * t561 + t549 * t566;
t541 = sin(qJ(4));
t570 = cos(qJ(4));
t490 = t508 * t541 - t521 * t570;
t491 = t508 * t570 + t521 * t541;
t555 = t571 * t560;
t507 = t528 * t542 + t539 * t555 + t546 * t571;
t583 = t584 * t490 + t587 * t491 + t585 * t507;
t529 = -t537 * t557 + t539 * t567;
t548 = t537 * t559 + t539 * t565;
t545 = t548 * t568;
t510 = t529 * t571 + (t537 * t560 - t545) * t542;
t522 = t537 * t561 + t548 * t566;
t492 = t510 * t541 - t522 * t570;
t493 = t510 * t570 + t522 * t541;
t509 = t529 * t542 - t537 * t555 + t545 * t571;
t582 = t584 * t492 + t587 * t493 + t585 * t509;
t581 = t585 * t490 + t586 * t491 + t588 * t507;
t580 = t585 * t492 + t586 * t493 + t588 * t509;
t579 = t587 * t490 + t589 * t491 + t586 * t507;
t578 = t587 * t492 + t589 * t493 + t586 * t509;
t556 = t568 * t567;
t558 = t569 * t566;
t520 = t542 * t558 + (t542 * t556 + t565 * t571) * t538;
t527 = -t560 * t567 + t568 * t569;
t511 = t520 * t541 - t527 * t570;
t512 = t520 * t570 + t527 * t541;
t519 = -t571 * t558 + (t542 * t565 - t556 * t571) * t538;
t577 = t584 * t511 + t587 * t512 + t585 * t519;
t576 = t585 * t511 + t586 * t512 + t588 * t519;
t575 = t587 * t511 + t589 * t512 + t586 * t519;
t517 = qJD(3) * t521;
t496 = qJD(4) * t507 + t517;
t518 = qJD(3) * t522;
t497 = qJD(4) * t509 + t518;
t526 = qJD(3) * t527;
t513 = qJD(4) * t519 + t526;
t564 = qJD(2) * t538;
t534 = qJD(2) * t569 + qJD(1);
t562 = t539 * t564;
t483 = pkin(3) * t508 + pkin(9) * t507;
t503 = pkin(3) * t520 + pkin(9) * t519;
t533 = t537 * t564;
t554 = -t483 * t526 + t503 * t517 + t533;
t484 = pkin(3) * t510 + pkin(9) * t509;
t553 = t483 * t518 - t484 * t517 + t534;
t485 = pkin(4) * t512 + qJ(5) * t511;
t552 = qJD(5) * t492 + t496 * t485 + t554;
t551 = t484 * t526 - t503 * t518 - t562;
t455 = pkin(4) * t491 + qJ(5) * t490;
t550 = qJD(5) * t511 + t497 * t455 + t553;
t456 = pkin(4) * t493 + qJ(5) * t492;
t547 = qJD(5) * t490 + t513 * t456 + t551;
t543 = cos(qJ(6));
t540 = sin(qJ(6));
t502 = rSges(4,1) * t520 - rSges(4,2) * t519 + rSges(4,3) * t527;
t501 = Icges(4,1) * t520 - Icges(4,4) * t519 + Icges(4,5) * t527;
t500 = Icges(4,4) * t520 - Icges(4,2) * t519 + Icges(4,6) * t527;
t499 = Icges(4,5) * t520 - Icges(4,6) * t519 + Icges(4,3) * t527;
t498 = pkin(5) * t519 + pkin(10) * t512;
t495 = t511 * t540 + t519 * t543;
t494 = t511 * t543 - t519 * t540;
t486 = qJD(6) * t512 + t513;
t481 = rSges(5,1) * t512 - rSges(5,2) * t511 + rSges(5,3) * t519;
t480 = rSges(6,1) * t519 - rSges(6,2) * t512 + rSges(6,3) * t511;
t472 = rSges(4,1) * t510 - rSges(4,2) * t509 + rSges(4,3) * t522;
t471 = rSges(4,1) * t508 - rSges(4,2) * t507 + rSges(4,3) * t521;
t470 = Icges(4,1) * t510 - Icges(4,4) * t509 + Icges(4,5) * t522;
t469 = Icges(4,1) * t508 - Icges(4,4) * t507 + Icges(4,5) * t521;
t468 = Icges(4,4) * t510 - Icges(4,2) * t509 + Icges(4,6) * t522;
t467 = Icges(4,4) * t508 - Icges(4,2) * t507 + Icges(4,6) * t521;
t466 = Icges(4,5) * t510 - Icges(4,6) * t509 + Icges(4,3) * t522;
t465 = Icges(4,5) * t508 - Icges(4,6) * t507 + Icges(4,3) * t521;
t464 = pkin(5) * t509 + pkin(10) * t493;
t463 = pkin(5) * t507 + pkin(10) * t491;
t462 = t492 * t540 + t509 * t543;
t461 = t492 * t543 - t509 * t540;
t460 = t490 * t540 + t507 * t543;
t459 = t490 * t543 - t507 * t540;
t458 = qJD(6) * t493 + t497;
t457 = qJD(6) * t491 + t496;
t452 = rSges(7,1) * t495 + rSges(7,2) * t494 + rSges(7,3) * t512;
t451 = Icges(7,1) * t495 + Icges(7,4) * t494 + Icges(7,5) * t512;
t450 = Icges(7,4) * t495 + Icges(7,2) * t494 + Icges(7,6) * t512;
t449 = Icges(7,5) * t495 + Icges(7,6) * t494 + Icges(7,3) * t512;
t448 = rSges(5,1) * t493 - rSges(5,2) * t492 + rSges(5,3) * t509;
t447 = rSges(5,1) * t491 - rSges(5,2) * t490 + rSges(5,3) * t507;
t446 = rSges(6,1) * t509 - rSges(6,2) * t493 + rSges(6,3) * t492;
t445 = rSges(6,1) * t507 - rSges(6,2) * t491 + rSges(6,3) * t490;
t431 = -t562 + (t472 * t527 - t502 * t522) * qJD(3);
t430 = t533 + (-t471 * t527 + t502 * t521) * qJD(3);
t429 = (t471 * t522 - t472 * t521) * qJD(3) + t534;
t428 = rSges(7,1) * t462 + rSges(7,2) * t461 + rSges(7,3) * t493;
t427 = rSges(7,1) * t460 + rSges(7,2) * t459 + rSges(7,3) * t491;
t426 = Icges(7,1) * t462 + Icges(7,4) * t461 + Icges(7,5) * t493;
t425 = Icges(7,1) * t460 + Icges(7,4) * t459 + Icges(7,5) * t491;
t424 = Icges(7,4) * t462 + Icges(7,2) * t461 + Icges(7,6) * t493;
t423 = Icges(7,4) * t460 + Icges(7,2) * t459 + Icges(7,6) * t491;
t422 = Icges(7,5) * t462 + Icges(7,6) * t461 + Icges(7,3) * t493;
t421 = Icges(7,5) * t460 + Icges(7,6) * t459 + Icges(7,3) * t491;
t420 = t448 * t513 - t481 * t497 + t551;
t419 = -t447 * t513 + t481 * t496 + t554;
t418 = t447 * t497 - t448 * t496 + t553;
t417 = t446 * t513 + (-t480 - t485) * t497 + t547;
t416 = t480 * t496 + (-t445 - t455) * t513 + t552;
t415 = t445 * t497 + (-t446 - t456) * t496 + t550;
t414 = t428 * t486 - t452 * t458 + t464 * t513 + (-t485 - t498) * t497 + t547;
t413 = -t427 * t486 + t452 * t457 + t496 * t498 + (-t455 - t463) * t513 + t552;
t412 = t427 * t458 - t428 * t457 + t463 * t497 + (-t456 - t464) * t496 + t550;
t1 = t458 * ((t493 * t422 + t461 * t424 + t462 * t426) * t458 + (t421 * t493 + t423 * t461 + t425 * t462) * t457 + (t449 * t493 + t450 * t461 + t451 * t462) * t486) / 0.2e1 + t457 * ((t422 * t491 + t424 * t459 + t426 * t460) * t458 + (t491 * t421 + t459 * t423 + t460 * t425) * t457 + (t449 * t491 + t450 * t459 + t451 * t460) * t486) / 0.2e1 + t486 * ((t422 * t512 + t424 * t494 + t426 * t495) * t458 + (t421 * t512 + t423 * t494 + t425 * t495) * t457 + (t449 * t512 + t450 * t494 + t451 * t495) * t486) / 0.2e1 + m(6) * (t415 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(7) * (t412 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(5) * (t418 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(4) * (t429 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + m(3) * (t534 ^ 2 + (t537 ^ 2 + t539 ^ 2) * qJD(2) ^ 2 * t538 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + ((t577 * t490 + t575 * t491 + t576 * t507) * t513 + (t582 * t490 + t578 * t491 + t580 * t507) * t497 + (t583 * t490 + t579 * t491 + t581 * t507) * t496) * t496 / 0.2e1 + ((t577 * t492 + t575 * t493 + t576 * t509) * t513 + (t582 * t492 + t578 * t493 + t580 * t509) * t497 + (t583 * t492 + t579 * t493 + t581 * t509) * t496) * t497 / 0.2e1 + ((t577 * t511 + t575 * t512 + t576 * t519) * t513 + (t582 * t511 + t578 * t512 + t580 * t519) * t497 + (t583 * t511 + t579 * t512 + t581 * t519) * t496) * t513 / 0.2e1 + (t522 * ((t466 * t522 - t468 * t509 + t470 * t510) * t522 + (t465 * t522 - t467 * t509 + t469 * t510) * t521 + (t499 * t522 - t509 * t500 + t501 * t510) * t527) + t521 * ((t466 * t521 - t468 * t507 + t470 * t508) * t522 + (t465 * t521 - t467 * t507 + t469 * t508) * t521 + (t499 * t521 - t500 * t507 + t501 * t508) * t527) + t527 * ((t466 * t527 - t468 * t519 + t470 * t520) * t522 + (t465 * t527 - t467 * t519 + t469 * t520) * t521 + (t499 * t527 - t500 * t519 + t501 * t520) * t527)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
