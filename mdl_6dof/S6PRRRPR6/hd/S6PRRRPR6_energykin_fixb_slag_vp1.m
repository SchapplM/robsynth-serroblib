% Calculate kinetic energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:09
% EndTime: 2019-03-08 23:31:12
% DurationCPUTime: 2.92s
% Computational Cost: add. (2906->330), mult. (7433->505), div. (0->0), fcn. (9322->12), ass. (0->149)
t557 = Icges(5,1) + Icges(6,1);
t556 = -Icges(5,4) + Icges(6,5);
t555 = Icges(6,4) + Icges(5,5);
t554 = Icges(5,2) + Icges(6,3);
t553 = Icges(6,2) + Icges(5,3);
t552 = -Icges(5,6) + Icges(6,6);
t508 = sin(pkin(11));
t510 = cos(pkin(11));
t517 = cos(qJ(2));
t511 = cos(pkin(6));
t515 = sin(qJ(2));
t532 = t511 * t515;
t496 = t508 * t517 + t510 * t532;
t514 = sin(qJ(3));
t509 = sin(pkin(6));
t533 = t510 * t509;
t538 = cos(qJ(3));
t479 = t496 * t538 - t514 * t533;
t531 = t511 * t517;
t495 = t508 * t515 - t510 * t531;
t513 = sin(qJ(4));
t537 = cos(qJ(4));
t450 = t479 * t513 - t495 * t537;
t451 = t479 * t537 + t495 * t513;
t528 = t509 * t538;
t478 = t496 * t514 + t510 * t528;
t551 = t554 * t450 + t556 * t451 + t552 * t478;
t498 = -t508 * t532 + t510 * t517;
t535 = t509 * t514;
t481 = t498 * t538 + t508 * t535;
t497 = t508 * t531 + t510 * t515;
t452 = t481 * t513 - t497 * t537;
t453 = t481 * t537 + t497 * t513;
t480 = t498 * t514 - t508 * t528;
t550 = t554 * t452 + t556 * t453 + t552 * t480;
t549 = t552 * t450 + t555 * t451 + t553 * t478;
t548 = t552 * t452 + t555 * t453 + t553 * t480;
t547 = t556 * t450 + t557 * t451 + t555 * t478;
t546 = t556 * t452 + t557 * t453 + t555 * t480;
t500 = t511 * t514 + t515 * t528;
t534 = t509 * t517;
t482 = t500 * t513 + t534 * t537;
t483 = t500 * t537 - t513 * t534;
t499 = -t511 * t538 + t515 * t535;
t545 = t554 * t482 + t556 * t483 + t552 * t499;
t544 = t552 * t482 + t555 * t483 + t553 * t499;
t543 = t556 * t482 + t557 * t483 + t555 * t499;
t542 = qJD(2) ^ 2;
t536 = t508 * t509;
t530 = qJD(2) * t509;
t505 = t508 * t530;
t484 = qJD(3) * t497 + t505;
t507 = qJD(2) * t511;
t446 = qJD(4) * t480 + t484;
t527 = t510 * t530;
t471 = pkin(2) * t496 + pkin(8) * t495;
t472 = pkin(2) * t498 + pkin(8) * t497;
t526 = t471 * t505 + t472 * t527 + qJD(1);
t485 = qJD(3) * t495 - t527;
t502 = -qJD(3) * t534 + t507;
t501 = (pkin(2) * t515 - pkin(8) * t517) * t509;
t525 = t472 * t507 - t501 * t505;
t447 = qJD(4) * t478 + t485;
t476 = qJD(4) * t499 + t502;
t443 = pkin(3) * t479 + pkin(9) * t478;
t444 = pkin(3) * t481 + pkin(9) * t480;
t524 = t484 * t443 - t444 * t485 + t526;
t523 = (-t471 * t511 - t501 * t533) * qJD(2);
t413 = pkin(4) * t451 + qJ(5) * t450;
t522 = qJD(5) * t482 + t446 * t413 + t524;
t473 = pkin(3) * t500 + pkin(9) * t499;
t521 = t502 * t444 - t473 * t484 + t525;
t414 = pkin(4) * t453 + qJ(5) * t452;
t520 = qJD(5) * t450 + t476 * t414 + t521;
t519 = -t443 * t502 + t485 * t473 + t523;
t445 = pkin(4) * t483 + qJ(5) * t482;
t518 = qJD(5) * t452 + t447 * t445 + t519;
t516 = cos(qJ(6));
t512 = sin(qJ(6));
t489 = t511 * rSges(3,3) + (rSges(3,1) * t515 + rSges(3,2) * t517) * t509;
t488 = Icges(3,5) * t511 + (Icges(3,1) * t515 + Icges(3,4) * t517) * t509;
t487 = Icges(3,6) * t511 + (Icges(3,4) * t515 + Icges(3,2) * t517) * t509;
t486 = Icges(3,3) * t511 + (Icges(3,5) * t515 + Icges(3,6) * t517) * t509;
t469 = t500 * rSges(4,1) - t499 * rSges(4,2) - rSges(4,3) * t534;
t468 = Icges(4,1) * t500 - Icges(4,4) * t499 - Icges(4,5) * t534;
t467 = Icges(4,4) * t500 - Icges(4,2) * t499 - Icges(4,6) * t534;
t466 = Icges(4,5) * t500 - Icges(4,6) * t499 - Icges(4,3) * t534;
t463 = rSges(3,1) * t498 - rSges(3,2) * t497 + rSges(3,3) * t536;
t462 = rSges(3,1) * t496 - rSges(3,2) * t495 - rSges(3,3) * t533;
t461 = Icges(3,1) * t498 - Icges(3,4) * t497 + Icges(3,5) * t536;
t460 = Icges(3,1) * t496 - Icges(3,4) * t495 - Icges(3,5) * t533;
t459 = Icges(3,4) * t498 - Icges(3,2) * t497 + Icges(3,6) * t536;
t458 = Icges(3,4) * t496 - Icges(3,2) * t495 - Icges(3,6) * t533;
t457 = Icges(3,5) * t498 - Icges(3,6) * t497 + Icges(3,3) * t536;
t456 = Icges(3,5) * t496 - Icges(3,6) * t495 - Icges(3,3) * t533;
t455 = pkin(5) * t483 - pkin(10) * t499;
t454 = -qJD(6) * t499 + t476;
t441 = t482 * t512 + t483 * t516;
t440 = t482 * t516 - t483 * t512;
t438 = (-t462 * t511 - t489 * t533) * qJD(2);
t437 = (t463 * t511 - t489 * t536) * qJD(2);
t436 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t499;
t435 = rSges(6,1) * t483 + rSges(6,2) * t499 + rSges(6,3) * t482;
t428 = rSges(4,1) * t481 - rSges(4,2) * t480 + rSges(4,3) * t497;
t427 = rSges(4,1) * t479 - rSges(4,2) * t478 + rSges(4,3) * t495;
t426 = Icges(4,1) * t481 - Icges(4,4) * t480 + Icges(4,5) * t497;
t425 = Icges(4,1) * t479 - Icges(4,4) * t478 + Icges(4,5) * t495;
t424 = Icges(4,4) * t481 - Icges(4,2) * t480 + Icges(4,6) * t497;
t423 = Icges(4,4) * t479 - Icges(4,2) * t478 + Icges(4,6) * t495;
t422 = Icges(4,5) * t481 - Icges(4,6) * t480 + Icges(4,3) * t497;
t421 = Icges(4,5) * t479 - Icges(4,6) * t478 + Icges(4,3) * t495;
t420 = pkin(5) * t453 - pkin(10) * t480;
t419 = pkin(5) * t451 - pkin(10) * t478;
t418 = -qJD(6) * t478 + t447;
t417 = -qJD(6) * t480 + t446;
t415 = qJD(1) + (t462 * t508 + t463 * t510) * t530;
t412 = t452 * t512 + t453 * t516;
t411 = t452 * t516 - t453 * t512;
t410 = t450 * t512 + t451 * t516;
t409 = t450 * t516 - t451 * t512;
t407 = rSges(5,1) * t453 - rSges(5,2) * t452 + rSges(5,3) * t480;
t406 = rSges(6,1) * t453 + rSges(6,2) * t480 + rSges(6,3) * t452;
t405 = rSges(5,1) * t451 - rSges(5,2) * t450 + rSges(5,3) * t478;
t404 = rSges(6,1) * t451 + rSges(6,2) * t478 + rSges(6,3) * t450;
t390 = rSges(7,1) * t441 + rSges(7,2) * t440 - rSges(7,3) * t499;
t389 = Icges(7,1) * t441 + Icges(7,4) * t440 - Icges(7,5) * t499;
t388 = Icges(7,4) * t441 + Icges(7,2) * t440 - Icges(7,6) * t499;
t387 = Icges(7,5) * t441 + Icges(7,6) * t440 - Icges(7,3) * t499;
t385 = -t427 * t502 + t469 * t485 + t523;
t384 = t428 * t502 - t469 * t484 + t525;
t383 = rSges(7,1) * t412 + rSges(7,2) * t411 - rSges(7,3) * t480;
t382 = rSges(7,1) * t410 + rSges(7,2) * t409 - rSges(7,3) * t478;
t381 = Icges(7,1) * t412 + Icges(7,4) * t411 - Icges(7,5) * t480;
t380 = Icges(7,1) * t410 + Icges(7,4) * t409 - Icges(7,5) * t478;
t379 = Icges(7,4) * t412 + Icges(7,2) * t411 - Icges(7,6) * t480;
t378 = Icges(7,4) * t410 + Icges(7,2) * t409 - Icges(7,6) * t478;
t377 = Icges(7,5) * t412 + Icges(7,6) * t411 - Icges(7,3) * t480;
t376 = Icges(7,5) * t410 + Icges(7,6) * t409 - Icges(7,3) * t478;
t375 = t427 * t484 - t428 * t485 + t526;
t374 = -t405 * t476 + t436 * t447 + t519;
t373 = t407 * t476 - t436 * t446 + t521;
t372 = t405 * t446 - t407 * t447 + t524;
t371 = t435 * t447 + (-t404 - t413) * t476 + t518;
t370 = t406 * t476 + (-t435 - t445) * t446 + t520;
t369 = t404 * t446 + (-t406 - t414) * t447 + t522;
t368 = -t382 * t454 + t390 * t418 + t447 * t455 + (-t413 - t419) * t476 + t518;
t367 = t383 * t454 - t390 * t417 + t420 * t476 + (-t445 - t455) * t446 + t520;
t366 = t382 * t417 - t383 * t418 + t419 * t446 + (-t414 - t420) * t447 + t522;
t1 = t454 * ((-t377 * t499 + t379 * t440 + t381 * t441) * t417 + (-t376 * t499 + t378 * t440 + t380 * t441) * t418 + (-t387 * t499 + t388 * t440 + t389 * t441) * t454) / 0.2e1 + t417 * ((-t377 * t480 + t379 * t411 + t381 * t412) * t417 + (-t376 * t480 + t378 * t411 + t380 * t412) * t418 + (-t387 * t480 + t388 * t411 + t389 * t412) * t454) / 0.2e1 + t418 * ((-t377 * t478 + t379 * t409 + t381 * t410) * t417 + (-t376 * t478 + t378 * t409 + t380 * t410) * t418 + (-t387 * t478 + t388 * t409 + t389 * t410) * t454) / 0.2e1 + t484 * ((t422 * t497 - t424 * t480 + t426 * t481) * t484 + (t421 * t497 - t423 * t480 + t425 * t481) * t485 + (t466 * t497 - t467 * t480 + t468 * t481) * t502) / 0.2e1 + t485 * ((t422 * t495 - t424 * t478 + t426 * t479) * t484 + (t421 * t495 - t423 * t478 + t425 * t479) * t485 + (t466 * t495 - t467 * t478 + t468 * t479) * t502) / 0.2e1 + t502 * ((-t422 * t534 - t499 * t424 + t500 * t426) * t484 + (-t421 * t534 - t499 * t423 + t500 * t425) * t485 + (-t466 * t534 - t499 * t467 + t500 * t468) * t502) / 0.2e1 + m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(4) * (t375 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(3) * (t415 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 - t542 * ((-t457 * t533 - t459 * t495 + t461 * t496) * t536 - (-t456 * t533 - t458 * t495 + t460 * t496) * t533 + (-t486 * t533 - t487 * t495 + t488 * t496) * t511) * t533 / 0.2e1 + ((t452 * t545 + t453 * t543 + t480 * t544) * t476 + (t452 * t551 + t547 * t453 + t549 * t480) * t447 + (t550 * t452 + t546 * t453 + t548 * t480) * t446) * t446 / 0.2e1 + ((t450 * t545 + t451 * t543 + t478 * t544) * t476 + (t551 * t450 + t547 * t451 + t549 * t478) * t447 + (t450 * t550 + t451 * t546 + t478 * t548) * t446) * t447 / 0.2e1 + ((t545 * t482 + t543 * t483 + t544 * t499) * t476 + (t482 * t551 + t547 * t483 + t549 * t499) * t447 + (t482 * t550 + t483 * t546 + t499 * t548) * t446) * t476 / 0.2e1 + (t511 * (t511 ^ 2 * t486 + (((t459 * t517 + t461 * t515) * t508 - (t458 * t517 + t460 * t515) * t510) * t509 + (-t456 * t510 + t457 * t508 + t487 * t517 + t488 * t515) * t511) * t509) + ((t457 * t536 - t459 * t497 + t461 * t498) * t536 - (t456 * t536 - t458 * t497 + t460 * t498) * t533 + (t486 * t536 - t487 * t497 + t488 * t498) * t511) * t536) * t542 / 0.2e1;
T  = t1;
