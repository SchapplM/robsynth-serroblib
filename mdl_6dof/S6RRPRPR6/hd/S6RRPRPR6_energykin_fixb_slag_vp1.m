% Calculate kinetic energy for
% S6RRPRPR6
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:08
% EndTime: 2019-03-09 10:38:11
% DurationCPUTime: 3.33s
% Computational Cost: add. (3031->339), mult. (7605->506), div. (0->0), fcn. (9548->12), ass. (0->157)
t602 = Icges(5,1) + Icges(6,2);
t601 = -Icges(6,1) - Icges(5,3);
t600 = -Icges(5,4) - Icges(6,6);
t599 = Icges(6,4) - Icges(5,5);
t598 = Icges(6,5) - Icges(5,6);
t597 = Icges(5,2) + Icges(6,3);
t537 = sin(qJ(1));
t539 = cos(qJ(1));
t536 = sin(qJ(2));
t573 = sin(pkin(11));
t575 = cos(pkin(6));
t551 = t575 * t573;
t574 = cos(pkin(11));
t552 = t575 * t574;
t579 = cos(qJ(2));
t544 = -t536 * t551 + t552 * t579;
t546 = t536 * t574 + t573 * t579;
t495 = -t537 * t546 + t539 * t544;
t511 = t536 * t552 + t551 * t579;
t519 = -t536 * t573 + t579 * t574;
t496 = t511 * t539 + t519 * t537;
t534 = sin(pkin(6));
t571 = t534 * t539;
t441 = Icges(4,5) * t496 + Icges(4,6) * t495 - Icges(4,3) * t571;
t497 = -t537 * t544 - t539 * t546;
t498 = -t511 * t537 + t519 * t539;
t572 = t534 * t537;
t442 = Icges(4,5) * t498 + Icges(4,6) * t497 + Icges(4,3) * t572;
t554 = t575 * t579;
t513 = -t537 * t536 + t539 * t554;
t561 = t536 * t575;
t514 = t537 * t579 + t539 * t561;
t482 = Icges(3,5) * t514 + Icges(3,6) * t513 - Icges(3,3) * t571;
t515 = -t539 * t536 - t537 * t554;
t516 = -t537 * t561 + t539 * t579;
t483 = Icges(3,5) * t516 + Icges(3,6) * t515 + Icges(3,3) * t572;
t596 = t534 * ((t441 + t482) * t539 + (-t442 - t483) * t537);
t578 = cos(qJ(4));
t565 = t534 * t578;
t577 = sin(qJ(4));
t468 = t496 * t577 + t539 * t565;
t564 = t534 * t577;
t469 = t496 * t578 - t539 * t564;
t595 = t597 * t468 + t600 * t469 - t598 * t495;
t470 = t498 * t577 - t537 * t565;
t471 = t498 * t578 + t537 * t564;
t594 = t597 * t470 + t600 * t471 - t598 * t497;
t593 = t598 * t468 - t599 * t469 + t601 * t495;
t592 = t598 * t470 - t599 * t471 + t601 * t497;
t591 = t600 * t468 + t602 * t469 + t599 * t495;
t590 = t600 * t470 + t602 * t471 + t599 * t497;
t510 = t546 * t534;
t500 = t510 * t577 - t575 * t578;
t501 = t510 * t578 + t575 * t577;
t509 = t519 * t534;
t589 = t597 * t500 + t600 * t501 - t598 * t509;
t588 = t598 * t500 - t599 * t501 + t601 * t509;
t587 = t600 * t500 + t602 * t501 + t599 * t509;
t586 = Icges(4,5) * t510 + Icges(4,6) * t509 + (Icges(3,5) * t536 + Icges(3,6) * t579) * t534 + (Icges(4,3) + Icges(3,3)) * t575;
t576 = pkin(2) * t579;
t562 = pkin(2) * t561 - qJ(3) * t534;
t494 = -t537 * t562 + t539 * t576;
t517 = qJD(1) * (pkin(1) * t539 + pkin(8) * t572);
t529 = qJD(2) * t575 + qJD(1);
t570 = t529 * t494 + t517;
t568 = qJD(2) * t534;
t528 = t537 * t568;
t473 = -qJD(4) * t497 + t528;
t569 = qJD(1) * (pkin(1) * t537 - pkin(8) * t571);
t567 = qJD(3) * t539;
t493 = t537 * t576 + t539 * t562;
t563 = t539 * t568;
t566 = qJD(3) * t575 + t493 * t528 + t494 * t563;
t502 = -qJD(4) * t509 + t529;
t520 = t534 * t536 * pkin(2) + qJ(3) * t575;
t559 = qJD(2) * (-t510 * rSges(4,1) - t509 * rSges(4,2) - rSges(4,3) * t575 - t520);
t558 = qJD(2) * (-pkin(3) * t510 + pkin(9) * t509 - t520);
t557 = qJD(3) * t572 - t569;
t474 = -qJD(4) * t495 - t563;
t460 = pkin(3) * t496 - pkin(9) * t495;
t461 = pkin(3) * t498 - pkin(9) * t497;
t550 = t460 * t528 + t461 * t563 + t566;
t428 = pkin(4) * t469 + qJ(5) * t468;
t547 = qJD(5) * t500 + t473 * t428 + t550;
t545 = t529 * t461 + (t537 * t558 - t567) * t534 + t570;
t543 = (-t460 - t493) * t529 + t558 * t571 + t557;
t429 = pkin(4) * t471 + qJ(5) * t470;
t542 = qJD(5) * t468 + t502 * t429 + t545;
t462 = pkin(4) * t501 + qJ(5) * t500;
t541 = qJD(5) * t470 + t474 * t462 + t543;
t538 = cos(qJ(6));
t535 = sin(qJ(6));
t523 = rSges(2,1) * t539 - rSges(2,2) * t537;
t522 = rSges(2,1) * t537 + rSges(2,2) * t539;
t508 = t575 * rSges(3,3) + (rSges(3,1) * t536 + rSges(3,2) * t579) * t534;
t507 = Icges(3,5) * t575 + (Icges(3,1) * t536 + Icges(3,4) * t579) * t534;
t506 = Icges(3,6) * t575 + (Icges(3,4) * t536 + Icges(3,2) * t579) * t534;
t489 = rSges(3,1) * t516 + rSges(3,2) * t515 + rSges(3,3) * t572;
t488 = rSges(3,1) * t514 + rSges(3,2) * t513 - rSges(3,3) * t571;
t487 = Icges(3,1) * t516 + Icges(3,4) * t515 + Icges(3,5) * t572;
t486 = Icges(3,1) * t514 + Icges(3,4) * t513 - Icges(3,5) * t571;
t485 = Icges(3,4) * t516 + Icges(3,2) * t515 + Icges(3,6) * t572;
t484 = Icges(3,4) * t514 + Icges(3,2) * t513 - Icges(3,6) * t571;
t480 = Icges(4,1) * t510 + Icges(4,4) * t509 + Icges(4,5) * t575;
t479 = Icges(4,4) * t510 + Icges(4,2) * t509 + Icges(4,6) * t575;
t472 = -pkin(5) * t509 + pkin(10) * t501;
t465 = t500 * t535 - t509 * t538;
t464 = t500 * t538 + t509 * t535;
t463 = qJD(6) * t501 + t502;
t459 = rSges(5,1) * t501 - rSges(5,2) * t500 - rSges(5,3) * t509;
t458 = -rSges(6,1) * t509 - rSges(6,2) * t501 + rSges(6,3) * t500;
t449 = rSges(4,1) * t498 + rSges(4,2) * t497 + rSges(4,3) * t572;
t448 = rSges(4,1) * t496 + rSges(4,2) * t495 - rSges(4,3) * t571;
t446 = Icges(4,1) * t498 + Icges(4,4) * t497 + Icges(4,5) * t572;
t445 = Icges(4,1) * t496 + Icges(4,4) * t495 - Icges(4,5) * t571;
t444 = Icges(4,4) * t498 + Icges(4,2) * t497 + Icges(4,6) * t572;
t443 = Icges(4,4) * t496 + Icges(4,2) * t495 - Icges(4,6) * t571;
t440 = t489 * t529 - t508 * t528 + t517;
t439 = -t488 * t529 - t508 * t563 - t569;
t438 = -pkin(5) * t497 + pkin(10) * t471;
t437 = -pkin(5) * t495 + pkin(10) * t469;
t436 = t470 * t535 - t497 * t538;
t435 = t470 * t538 + t497 * t535;
t434 = t468 * t535 - t495 * t538;
t433 = t468 * t538 + t495 * t535;
t432 = (t488 * t537 + t489 * t539) * t568;
t431 = qJD(6) * t469 + t474;
t430 = qJD(6) * t471 + t473;
t425 = rSges(7,1) * t465 + rSges(7,2) * t464 + rSges(7,3) * t501;
t424 = Icges(7,1) * t465 + Icges(7,4) * t464 + Icges(7,5) * t501;
t423 = Icges(7,4) * t465 + Icges(7,2) * t464 + Icges(7,6) * t501;
t422 = Icges(7,5) * t465 + Icges(7,6) * t464 + Icges(7,3) * t501;
t421 = rSges(5,1) * t471 - rSges(5,2) * t470 - rSges(5,3) * t497;
t420 = rSges(5,1) * t469 - rSges(5,2) * t468 - rSges(5,3) * t495;
t419 = -rSges(6,1) * t497 - rSges(6,2) * t471 + rSges(6,3) * t470;
t418 = -rSges(6,1) * t495 - rSges(6,2) * t469 + rSges(6,3) * t468;
t404 = t449 * t529 + (t537 * t559 - t567) * t534 + t570;
t403 = (-t448 - t493) * t529 + t559 * t571 + t557;
t402 = rSges(7,1) * t436 + rSges(7,2) * t435 + rSges(7,3) * t471;
t401 = rSges(7,1) * t434 + rSges(7,2) * t433 + rSges(7,3) * t469;
t400 = Icges(7,1) * t436 + Icges(7,4) * t435 + Icges(7,5) * t471;
t399 = Icges(7,1) * t434 + Icges(7,4) * t433 + Icges(7,5) * t469;
t398 = Icges(7,4) * t436 + Icges(7,2) * t435 + Icges(7,6) * t471;
t397 = Icges(7,4) * t434 + Icges(7,2) * t433 + Icges(7,6) * t469;
t396 = Icges(7,5) * t436 + Icges(7,6) * t435 + Icges(7,3) * t471;
t395 = Icges(7,5) * t434 + Icges(7,6) * t433 + Icges(7,3) * t469;
t394 = (t448 * t537 + t449 * t539) * t568 + t566;
t393 = t421 * t502 - t459 * t473 + t545;
t392 = -t420 * t502 + t459 * t474 + t543;
t391 = t420 * t473 - t421 * t474 + t550;
t390 = t419 * t502 + (-t458 - t462) * t473 + t542;
t389 = t458 * t474 + (-t418 - t428) * t502 + t541;
t388 = t418 * t473 + (-t419 - t429) * t474 + t547;
t387 = t402 * t463 - t425 * t430 + t438 * t502 + (-t462 - t472) * t473 + t542;
t386 = -t401 * t463 + t425 * t431 + t472 * t474 + (-t428 - t437) * t502 + t541;
t385 = t401 * t430 - t402 * t431 + t437 * t473 + (-t429 - t438) * t474 + t547;
t1 = m(5) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + m(6) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(7) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + t430 * ((t396 * t471 + t398 * t435 + t400 * t436) * t430 + (t395 * t471 + t397 * t435 + t399 * t436) * t431 + (t422 * t471 + t423 * t435 + t424 * t436) * t463) / 0.2e1 + t431 * ((t396 * t469 + t398 * t433 + t400 * t434) * t430 + (t395 * t469 + t397 * t433 + t399 * t434) * t431 + (t422 * t469 + t423 * t433 + t424 * t434) * t463) / 0.2e1 + t463 * ((t396 * t501 + t398 * t464 + t400 * t465) * t430 + (t395 * t501 + t397 * t464 + t399 * t465) * t431 + (t422 * t501 + t423 * t464 + t424 * t465) * t463) / 0.2e1 + m(3) * (t432 ^ 2 + t439 ^ 2 + t440 ^ 2) / 0.2e1 + m(4) * (t394 ^ 2 + t403 ^ 2 + t404 ^ 2) / 0.2e1 + ((t470 * t589 + t471 * t587 - t497 * t588) * t502 + (t470 * t595 + t591 * t471 - t593 * t497) * t474 + (t470 * t594 + t471 * t590 - t497 * t592) * t473) * t473 / 0.2e1 + ((t468 * t589 + t469 * t587 - t495 * t588) * t502 + (t468 * t595 + t591 * t469 - t593 * t495) * t474 + (t468 * t594 + t469 * t590 - t495 * t592) * t473) * t474 / 0.2e1 + ((t500 * t589 + t501 * t587 - t509 * t588) * t502 + (t500 * t595 + t591 * t501 - t593 * t509) * t474 + (t500 * t594 + t501 * t590 - t509 * t592) * t473) * t502 / 0.2e1 + ((t575 * t483 + (t485 * t579 + t487 * t536) * t534) * t528 - (t575 * t482 + (t484 * t579 + t486 * t536) * t534) * t563 + ((t442 * t575 + t509 * t444 + t510 * t446) * t537 - (t441 * t575 + t509 * t443 + t510 * t445) * t539) * t568 + ((t506 * t579 + t507 * t536) * t534 + t509 * t479 + t510 * t480 + t586 * t575) * t529) * t529 / 0.2e1 + (Icges(2,3) + m(2) * (t522 ^ 2 + t523 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t443 * t497 - t445 * t498 - t484 * t515 - t486 * t516) * t539 + (t497 * t444 + t498 * t446 + t515 * t485 + t516 * t487 - t596) * t537) * t568 + (t479 * t497 + t480 * t498 + t506 * t515 + t507 * t516 + t572 * t586) * t529) * t528 / 0.2e1 - (((-t495 * t443 - t496 * t445 - t513 * t484 - t514 * t486 + t596) * t539 + (t444 * t495 + t446 * t496 + t485 * t513 + t487 * t514) * t537) * t568 + (t479 * t495 + t480 * t496 + t506 * t513 + t507 * t514 - t571 * t586) * t529) * t563 / 0.2e1;
T  = t1;
