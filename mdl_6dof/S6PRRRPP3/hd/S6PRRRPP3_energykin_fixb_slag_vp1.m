% Calculate kinetic energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:25
% EndTime: 2019-03-08 22:54:27
% DurationCPUTime: 2.32s
% Computational Cost: add. (2507->272), mult. (6404->412), div. (0->0), fcn. (7909->10), ass. (0->128)
t543 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t542 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t541 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t540 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t539 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t538 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t537 = rSges(7,1) + pkin(5);
t536 = rSges(7,3) + qJ(6);
t491 = sin(pkin(10));
t493 = cos(pkin(10));
t497 = cos(qJ(2));
t494 = cos(pkin(6));
t496 = sin(qJ(2));
t515 = t494 * t496;
t478 = t491 * t497 + t493 * t515;
t492 = sin(pkin(6));
t495 = sin(qJ(3));
t517 = t492 * t495;
t522 = cos(qJ(3));
t462 = t478 * t522 - t493 * t517;
t514 = t494 * t497;
t477 = t491 * t496 - t493 * t514;
t520 = sin(qJ(4));
t521 = cos(qJ(4));
t436 = t462 * t520 - t477 * t521;
t437 = t462 * t521 + t477 * t520;
t508 = t492 * t522;
t461 = t478 * t495 + t493 * t508;
t535 = t541 * t436 + t543 * t437 + t540 * t461;
t480 = -t491 * t515 + t493 * t497;
t464 = t480 * t522 + t491 * t517;
t479 = t491 * t514 + t493 * t496;
t438 = t464 * t520 - t479 * t521;
t439 = t464 * t521 + t479 * t520;
t463 = t480 * t495 - t491 * t508;
t534 = t541 * t438 + t543 * t439 + t540 * t463;
t533 = t538 * t436 + t541 * t437 + t539 * t461;
t532 = t538 * t438 + t541 * t439 + t539 * t463;
t531 = t539 * t436 + t540 * t437 + t542 * t461;
t530 = t539 * t438 + t540 * t439 + t542 * t463;
t482 = t494 * t495 + t496 * t508;
t516 = t492 * t497;
t465 = t482 * t520 + t516 * t521;
t466 = t482 * t521 - t516 * t520;
t481 = -t494 * t522 + t496 * t517;
t529 = t541 * t465 + t543 * t466 + t540 * t481;
t528 = t538 * t465 + t541 * t466 + t539 * t481;
t527 = t539 * t465 + t540 * t466 + t542 * t481;
t526 = qJD(2) ^ 2;
t519 = t491 * t492;
t518 = t492 * t493;
t513 = rSges(7,2) * t436 + t536 * t437 + t537 * t461;
t512 = rSges(7,2) * t438 + t536 * t439 + t537 * t463;
t511 = rSges(7,2) * t465 + t536 * t466 + t537 * t481;
t510 = qJD(2) * t492;
t488 = t491 * t510;
t467 = qJD(3) * t479 + t488;
t490 = qJD(2) * t494;
t507 = t493 * t510;
t456 = pkin(2) * t478 + pkin(8) * t477;
t457 = pkin(2) * t480 + pkin(8) * t479;
t506 = t456 * t488 + t457 * t507 + qJD(1);
t468 = qJD(3) * t477 - t507;
t484 = -qJD(3) * t516 + t490;
t483 = (pkin(2) * t496 - pkin(8) * t497) * t492;
t505 = t457 * t490 - t483 * t488;
t429 = pkin(3) * t462 + pkin(9) * t461;
t430 = pkin(3) * t464 + pkin(9) * t463;
t504 = t467 * t429 - t430 * t468 + t506;
t503 = (-t456 * t494 - t483 * t518) * qJD(2);
t399 = pkin(4) * t437 + qJ(5) * t436;
t432 = qJD(4) * t463 + t467;
t502 = qJD(5) * t465 + t432 * t399 + t504;
t458 = pkin(3) * t482 + pkin(9) * t481;
t501 = t484 * t430 - t458 * t467 + t505;
t400 = pkin(4) * t439 + qJ(5) * t438;
t459 = qJD(4) * t481 + t484;
t500 = qJD(5) * t436 + t459 * t400 + t501;
t499 = -t429 * t484 + t468 * t458 + t503;
t431 = pkin(4) * t466 + qJ(5) * t465;
t433 = qJD(4) * t461 + t468;
t498 = qJD(5) * t438 + t433 * t431 + t499;
t472 = t494 * rSges(3,3) + (rSges(3,1) * t496 + rSges(3,2) * t497) * t492;
t471 = Icges(3,5) * t494 + (Icges(3,1) * t496 + Icges(3,4) * t497) * t492;
t470 = Icges(3,6) * t494 + (Icges(3,4) * t496 + Icges(3,2) * t497) * t492;
t469 = Icges(3,3) * t494 + (Icges(3,5) * t496 + Icges(3,6) * t497) * t492;
t454 = t482 * rSges(4,1) - t481 * rSges(4,2) - rSges(4,3) * t516;
t453 = Icges(4,1) * t482 - Icges(4,4) * t481 - Icges(4,5) * t516;
t452 = Icges(4,4) * t482 - Icges(4,2) * t481 - Icges(4,6) * t516;
t451 = Icges(4,5) * t482 - Icges(4,6) * t481 - Icges(4,3) * t516;
t448 = rSges(3,1) * t480 - rSges(3,2) * t479 + rSges(3,3) * t519;
t447 = rSges(3,1) * t478 - rSges(3,2) * t477 - rSges(3,3) * t518;
t446 = Icges(3,1) * t480 - Icges(3,4) * t479 + Icges(3,5) * t519;
t445 = Icges(3,1) * t478 - Icges(3,4) * t477 - Icges(3,5) * t518;
t444 = Icges(3,4) * t480 - Icges(3,2) * t479 + Icges(3,6) * t519;
t443 = Icges(3,4) * t478 - Icges(3,2) * t477 - Icges(3,6) * t518;
t442 = Icges(3,5) * t480 - Icges(3,6) * t479 + Icges(3,3) * t519;
t441 = Icges(3,5) * t478 - Icges(3,6) * t477 - Icges(3,3) * t518;
t426 = (-t447 * t494 - t472 * t518) * qJD(2);
t425 = (t448 * t494 - t472 * t519) * qJD(2);
t424 = rSges(5,1) * t466 - rSges(5,2) * t465 + rSges(5,3) * t481;
t423 = rSges(6,1) * t481 - rSges(6,2) * t466 + rSges(6,3) * t465;
t412 = rSges(4,1) * t464 - rSges(4,2) * t463 + rSges(4,3) * t479;
t411 = rSges(4,1) * t462 - rSges(4,2) * t461 + rSges(4,3) * t477;
t410 = Icges(4,1) * t464 - Icges(4,4) * t463 + Icges(4,5) * t479;
t409 = Icges(4,1) * t462 - Icges(4,4) * t461 + Icges(4,5) * t477;
t408 = Icges(4,4) * t464 - Icges(4,2) * t463 + Icges(4,6) * t479;
t407 = Icges(4,4) * t462 - Icges(4,2) * t461 + Icges(4,6) * t477;
t406 = Icges(4,5) * t464 - Icges(4,6) * t463 + Icges(4,3) * t479;
t405 = Icges(4,5) * t462 - Icges(4,6) * t461 + Icges(4,3) * t477;
t401 = qJD(1) + (t447 * t491 + t448 * t493) * t510;
t397 = rSges(5,1) * t439 - rSges(5,2) * t438 + rSges(5,3) * t463;
t396 = rSges(5,1) * t437 - rSges(5,2) * t436 + rSges(5,3) * t461;
t395 = rSges(6,1) * t463 - rSges(6,2) * t439 + rSges(6,3) * t438;
t393 = rSges(6,1) * t461 - rSges(6,2) * t437 + rSges(6,3) * t436;
t371 = -t411 * t484 + t454 * t468 + t503;
t370 = t412 * t484 - t454 * t467 + t505;
t369 = t411 * t467 - t412 * t468 + t506;
t368 = -t396 * t459 + t424 * t433 + t499;
t367 = t397 * t459 - t424 * t432 + t501;
t366 = t396 * t432 - t397 * t433 + t504;
t365 = t423 * t433 + (-t393 - t399) * t459 + t498;
t364 = t395 * t459 + (-t423 - t431) * t432 + t500;
t363 = t393 * t432 + (-t395 - t400) * t433 + t502;
t362 = qJD(6) * t439 + t511 * t433 + (-t399 - t513) * t459 + t498;
t361 = qJD(6) * t437 + t512 * t459 + (-t431 - t511) * t432 + t500;
t360 = qJD(6) * t466 + t513 * t432 + (-t400 - t512) * t433 + t502;
t1 = -t526 * ((-t442 * t518 - t444 * t477 + t446 * t478) * t519 - (-t441 * t518 - t443 * t477 + t445 * t478) * t518 + (-t469 * t518 - t470 * t477 + t471 * t478) * t494) * t518 / 0.2e1 + t467 * ((t406 * t479 - t408 * t463 + t410 * t464) * t467 + (t405 * t479 - t407 * t463 + t409 * t464) * t468 + (t451 * t479 - t452 * t463 + t453 * t464) * t484) / 0.2e1 + t468 * ((t406 * t477 - t408 * t461 + t410 * t462) * t467 + (t405 * t477 - t407 * t461 + t409 * t462) * t468 + (t451 * t477 - t452 * t461 + t453 * t462) * t484) / 0.2e1 + t484 * ((-t406 * t516 - t481 * t408 + t482 * t410) * t467 + (-t405 * t516 - t481 * t407 + t482 * t409) * t468 + (-t451 * t516 - t481 * t452 + t482 * t453) * t484) / 0.2e1 + m(6) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(7) * (t360 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + m(5) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(4) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(3) * (t401 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + (((t442 * t519 - t444 * t479 + t446 * t480) * t519 - (t441 * t519 - t443 * t479 + t445 * t480) * t518 + (t469 * t519 - t470 * t479 + t471 * t480) * t494) * t519 + t494 * (t494 ^ 2 * t469 + (((t444 * t497 + t446 * t496) * t491 - (t443 * t497 + t445 * t496) * t493) * t492 + (-t441 * t493 + t442 * t491 + t470 * t497 + t471 * t496) * t494) * t492)) * t526 / 0.2e1 + ((t438 * t528 + t439 * t529 + t463 * t527) * t459 + (t533 * t438 + t439 * t535 + t531 * t463) * t433 + (t532 * t438 + t534 * t439 + t530 * t463) * t432) * t432 / 0.2e1 + ((t436 * t528 + t437 * t529 + t461 * t527) * t459 + (t533 * t436 + t535 * t437 + t531 * t461) * t433 + (t436 * t532 + t437 * t534 + t461 * t530) * t432) * t433 / 0.2e1 + ((t528 * t465 + t529 * t466 + t527 * t481) * t459 + (t533 * t465 + t466 * t535 + t531 * t481) * t433 + (t465 * t532 + t466 * t534 + t481 * t530) * t432) * t459 / 0.2e1;
T  = t1;
