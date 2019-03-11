% Calculate kinetic energy for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:18
% EndTime: 2019-03-09 09:28:21
% DurationCPUTime: 2.85s
% Computational Cost: add. (1600->287), mult. (4009->437), div. (0->0), fcn. (4628->10), ass. (0->136)
t550 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t549 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t548 = Icges(3,5) + Icges(5,5) - Icges(4,4);
t547 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t546 = Icges(3,6) - Icges(4,5) - Icges(5,4);
t545 = Icges(3,3) + Icges(5,1) + Icges(4,1);
t491 = sin(pkin(6));
t492 = cos(pkin(6));
t498 = cos(qJ(2));
t499 = cos(qJ(1));
t525 = t498 * t499;
t495 = sin(qJ(2));
t496 = sin(qJ(1));
t528 = t495 * t496;
t471 = -t492 * t525 + t528;
t526 = t496 * t498;
t527 = t495 * t499;
t472 = t492 * t527 + t526;
t473 = t492 * t526 + t527;
t474 = -t492 * t528 + t525;
t529 = t491 * t499;
t531 = t491 * t496;
t536 = (t546 * t471 - t548 * t472 + t545 * t529) * t499 + (-t546 * t473 + t548 * t474 + t545 * t531) * t496;
t544 = t536 * t491;
t543 = -t549 * t473 + t550 * t474 + t548 * t531;
t542 = t549 * t471 - t550 * t472 + t548 * t529;
t541 = t547 * t473 - t549 * t474 - t546 * t531;
t540 = t547 * t471 - t549 * t472 + t546 * t529;
t539 = t546 * t492 + (t549 * t495 + t547 * t498) * t491;
t538 = t548 * t492 + (t550 * t495 + t549 * t498) * t491;
t537 = t545 * t492 + (t548 * t495 + t546 * t498) * t491;
t533 = cos(qJ(5));
t532 = t491 * t495;
t530 = t491 * t498;
t429 = pkin(2) * t472 + qJ(3) * t471;
t430 = pkin(2) * t474 + qJ(3) * t473;
t520 = qJD(2) * t491;
t488 = t496 * t520;
t516 = t499 * t520;
t524 = t429 * t488 + t430 * t516;
t447 = -pkin(3) * t529 + qJ(4) * t472;
t523 = -t429 - t447;
t475 = (pkin(2) * t495 - qJ(3) * t498) * t491;
t522 = -pkin(3) * t492 - qJ(4) * t532 - t475;
t444 = -qJD(5) * t473 + t488;
t521 = qJD(1) * (pkin(1) * t496 - pkin(8) * t529);
t519 = qJD(3) * t498;
t489 = qJD(2) * t492 + qJD(1);
t477 = qJD(1) * (pkin(1) * t499 + pkin(8) * t531);
t518 = qJD(3) * t471 + t489 * t430 + t477;
t517 = t491 * t533;
t476 = qJD(5) * t530 + t489;
t515 = qJD(3) * t473 - t521;
t512 = (-rSges(4,1) * t492 - (-rSges(4,2) * t495 - rSges(4,3) * t498) * t491 - t475) * t520;
t445 = -qJD(5) * t471 - t516;
t511 = qJD(4) * t474 + t515;
t446 = pkin(3) * t531 + qJ(4) * t474;
t510 = qJD(4) * t472 + t489 * t446 + t518;
t509 = qJD(4) * t532 + t446 * t516 + t447 * t488 + t524;
t505 = (-rSges(5,1) * t492 - (-rSges(5,2) * t498 + rSges(5,3) * t495) * t491 + t522) * t520;
t504 = (-pkin(4) * t492 - pkin(9) * t530 + t522) * t520;
t448 = pkin(4) * t531 - pkin(9) * t473;
t449 = -pkin(4) * t529 - pkin(9) * t471;
t503 = t448 * t516 + t449 * t488 - t491 * t519 + t509;
t502 = t489 * t448 + t496 * t504 + t510;
t501 = (-t449 + t523) * t489 + t499 * t504 + t511;
t497 = cos(qJ(6));
t494 = sin(qJ(5));
t493 = sin(qJ(6));
t482 = rSges(2,1) * t499 - rSges(2,2) * t496;
t481 = rSges(2,1) * t496 + rSges(2,2) * t499;
t470 = t492 * t533 + t494 * t532;
t469 = t492 * t494 - t495 * t517;
t459 = rSges(3,3) * t492 + (rSges(3,1) * t495 + rSges(3,2) * t498) * t491;
t443 = t472 * t494 - t499 * t517;
t442 = t472 * t533 + t494 * t529;
t441 = t474 * t494 + t496 * t517;
t440 = -t474 * t533 + t494 * t531;
t439 = t470 * t497 + t493 * t530;
t438 = -t470 * t493 + t497 * t530;
t434 = qJD(6) * t469 + t476;
t428 = pkin(5) * t470 + pkin(10) * t469;
t427 = rSges(3,1) * t474 - rSges(3,2) * t473 + rSges(3,3) * t531;
t426 = rSges(3,1) * t472 - rSges(3,2) * t471 - rSges(3,3) * t529;
t425 = -rSges(4,1) * t529 - rSges(4,2) * t472 + rSges(4,3) * t471;
t424 = -rSges(5,1) * t529 + rSges(5,2) * t471 + rSges(5,3) * t472;
t423 = rSges(4,1) * t531 - rSges(4,2) * t474 + rSges(4,3) * t473;
t422 = rSges(5,1) * t531 + rSges(5,2) * t473 + rSges(5,3) * t474;
t400 = rSges(6,1) * t470 - rSges(6,2) * t469 + rSges(6,3) * t530;
t399 = Icges(6,1) * t470 - Icges(6,4) * t469 + Icges(6,5) * t530;
t398 = Icges(6,4) * t470 - Icges(6,2) * t469 + Icges(6,6) * t530;
t397 = Icges(6,5) * t470 - Icges(6,6) * t469 + Icges(6,3) * t530;
t396 = t443 * t497 - t471 * t493;
t395 = -t443 * t493 - t471 * t497;
t394 = t441 * t497 - t473 * t493;
t393 = -t441 * t493 - t473 * t497;
t392 = -qJD(6) * t442 + t445;
t391 = qJD(6) * t440 + t444;
t390 = pkin(5) * t443 - pkin(10) * t442;
t389 = pkin(5) * t441 + pkin(10) * t440;
t388 = rSges(6,1) * t443 + rSges(6,2) * t442 - rSges(6,3) * t471;
t387 = rSges(6,1) * t441 - rSges(6,2) * t440 - rSges(6,3) * t473;
t386 = Icges(6,1) * t443 + Icges(6,4) * t442 - Icges(6,5) * t471;
t385 = Icges(6,1) * t441 - Icges(6,4) * t440 - Icges(6,5) * t473;
t384 = Icges(6,4) * t443 + Icges(6,2) * t442 - Icges(6,6) * t471;
t383 = Icges(6,4) * t441 - Icges(6,2) * t440 - Icges(6,6) * t473;
t382 = Icges(6,5) * t443 + Icges(6,6) * t442 - Icges(6,3) * t471;
t381 = Icges(6,5) * t441 - Icges(6,6) * t440 - Icges(6,3) * t473;
t380 = rSges(7,1) * t439 + rSges(7,2) * t438 + rSges(7,3) * t469;
t379 = Icges(7,1) * t439 + Icges(7,4) * t438 + Icges(7,5) * t469;
t378 = Icges(7,4) * t439 + Icges(7,2) * t438 + Icges(7,6) * t469;
t377 = Icges(7,5) * t439 + Icges(7,6) * t438 + Icges(7,3) * t469;
t376 = t427 * t489 - t459 * t488 + t477;
t375 = -t426 * t489 - t459 * t516 - t521;
t374 = (t426 * t496 + t427 * t499) * t520;
t373 = rSges(7,1) * t396 + rSges(7,2) * t395 - rSges(7,3) * t442;
t372 = rSges(7,1) * t394 + rSges(7,2) * t393 + rSges(7,3) * t440;
t371 = Icges(7,1) * t396 + Icges(7,4) * t395 - Icges(7,5) * t442;
t370 = Icges(7,1) * t394 + Icges(7,4) * t393 + Icges(7,5) * t440;
t369 = Icges(7,4) * t396 + Icges(7,2) * t395 - Icges(7,6) * t442;
t368 = Icges(7,4) * t394 + Icges(7,2) * t393 + Icges(7,6) * t440;
t367 = Icges(7,5) * t396 + Icges(7,6) * t395 - Icges(7,3) * t442;
t366 = Icges(7,5) * t394 + Icges(7,6) * t393 + Icges(7,3) * t440;
t365 = t423 * t489 + t496 * t512 + t518;
t364 = (-t425 - t429) * t489 + t499 * t512 + t515;
t363 = (-t519 + (t423 * t499 + t425 * t496) * qJD(2)) * t491 + t524;
t362 = t422 * t489 + t496 * t505 + t510;
t361 = (-t424 + t523) * t489 + t499 * t505 + t511;
t360 = (-t519 + (t422 * t499 + t424 * t496) * qJD(2)) * t491 + t509;
t359 = t387 * t476 - t400 * t444 + t502;
t358 = -t388 * t476 + t400 * t445 + t501;
t357 = -t387 * t445 + t388 * t444 + t503;
t356 = t372 * t434 - t380 * t391 + t389 * t476 - t428 * t444 + t502;
t355 = -t373 * t434 + t380 * t392 - t390 * t476 + t428 * t445 + t501;
t354 = -t372 * t392 + t373 * t391 - t389 * t445 + t390 * t444 + t503;
t1 = m(7) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(6) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(4) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(5) * (t360 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + m(3) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + t391 * ((t440 * t366 + t393 * t368 + t394 * t370) * t391 + (t367 * t440 + t369 * t393 + t371 * t394) * t392 + (t377 * t440 + t378 * t393 + t379 * t394) * t434) / 0.2e1 + t434 * ((t366 * t469 + t368 * t438 + t370 * t439) * t391 + (t367 * t469 + t369 * t438 + t371 * t439) * t392 + (t469 * t377 + t438 * t378 + t439 * t379) * t434) / 0.2e1 + t476 * ((t381 * t530 - t383 * t469 + t385 * t470) * t444 + (t382 * t530 - t384 * t469 + t386 * t470) * t445 + (t397 * t530 - t469 * t398 + t470 * t399) * t476) / 0.2e1 + t392 * ((-t366 * t442 + t368 * t395 + t370 * t396) * t391 + (-t442 * t367 + t395 * t369 + t396 * t371) * t392 + (-t377 * t442 + t378 * t395 + t379 * t396) * t434) / 0.2e1 + t444 * ((-t473 * t381 - t440 * t383 + t441 * t385) * t444 + (-t382 * t473 - t384 * t440 + t386 * t441) * t445 + (-t397 * t473 - t398 * t440 + t399 * t441) * t476) / 0.2e1 + t445 * ((-t381 * t471 + t383 * t442 + t385 * t443) * t444 + (-t471 * t382 + t442 * t384 + t443 * t386) * t445 + (-t397 * t471 + t398 * t442 + t399 * t443) * t476) / 0.2e1 + (Icges(2,3) + m(2) * (t481 ^ 2 + t482 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t536 * t492 + ((t495 * t542 + t498 * t540) * t499 + (t495 * t543 - t541 * t498) * t496) * t491) * t520 + (t537 * t492 + (t495 * t538 + t498 * t539) * t491) * t489) * t489 / 0.2e1 + (((-t473 * t540 + t474 * t542) * t499 + (t541 * t473 + t474 * t543 + t544) * t496) * t520 + (-t473 * t539 + t474 * t538 + t531 * t537) * t489) * t488 / 0.2e1 - (((-t471 * t540 + t472 * t542 - t544) * t499 + (t541 * t471 + t472 * t543) * t496) * t520 + (-t471 * t539 + t472 * t538 - t529 * t537) * t489) * t516 / 0.2e1;
T  = t1;
