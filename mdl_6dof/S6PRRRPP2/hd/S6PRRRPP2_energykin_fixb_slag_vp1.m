% Calculate kinetic energy for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:44
% EndTime: 2019-03-08 22:48:46
% DurationCPUTime: 2.34s
% Computational Cost: add. (2502->272), mult. (6391->412), div. (0->0), fcn. (7891->10), ass. (0->128)
t539 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t538 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t537 = Icges(6,4) + Icges(5,5) - Icges(7,5);
t536 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t535 = Icges(6,2) + Icges(5,3) + Icges(7,3);
t534 = -Icges(5,6) + Icges(6,6) - Icges(7,6);
t533 = rSges(7,1) + pkin(5);
t532 = -rSges(7,3) - qJ(6);
t487 = sin(pkin(10));
t489 = cos(pkin(10));
t494 = cos(qJ(2));
t490 = cos(pkin(6));
t493 = sin(qJ(2));
t512 = t490 * t493;
t475 = t487 * t494 + t489 * t512;
t492 = sin(qJ(3));
t488 = sin(pkin(6));
t513 = t489 * t488;
t518 = cos(qJ(3));
t459 = t475 * t518 - t492 * t513;
t511 = t490 * t494;
t474 = t487 * t493 - t489 * t511;
t491 = sin(qJ(4));
t517 = cos(qJ(4));
t433 = t459 * t491 - t474 * t517;
t434 = t459 * t517 + t474 * t491;
t505 = t488 * t518;
t458 = t475 * t492 + t489 * t505;
t531 = t534 * t433 + t537 * t434 + t535 * t458;
t477 = -t487 * t512 + t489 * t494;
t515 = t488 * t492;
t461 = t477 * t518 + t487 * t515;
t476 = t487 * t511 + t489 * t493;
t435 = t461 * t491 - t476 * t517;
t436 = t461 * t517 + t476 * t491;
t460 = t477 * t492 - t487 * t505;
t530 = t534 * t435 + t537 * t436 + t535 * t460;
t529 = t536 * t433 + t538 * t434 + t534 * t458;
t528 = t536 * t435 + t538 * t436 + t534 * t460;
t527 = t538 * t433 + t539 * t434 + t537 * t458;
t526 = t538 * t435 + t539 * t436 + t537 * t460;
t479 = t490 * t492 + t493 * t505;
t514 = t488 * t494;
t462 = t479 * t491 + t514 * t517;
t463 = t479 * t517 - t491 * t514;
t478 = -t490 * t518 + t493 * t515;
t525 = t534 * t462 + t537 * t463 + t535 * t478;
t524 = t536 * t462 + t538 * t463 + t534 * t478;
t523 = t538 * t462 + t539 * t463 + t537 * t478;
t522 = qJD(2) ^ 2;
t516 = t487 * t488;
t510 = rSges(7,2) * t433 + t533 * t434 + t532 * t458;
t509 = rSges(7,2) * t435 + t533 * t436 + t532 * t460;
t508 = rSges(7,2) * t462 + t533 * t463 + t532 * t478;
t507 = qJD(2) * t488;
t484 = t487 * t507;
t464 = qJD(3) * t476 + t484;
t486 = qJD(2) * t490;
t504 = t489 * t507;
t453 = pkin(2) * t475 + pkin(8) * t474;
t454 = pkin(2) * t477 + pkin(8) * t476;
t503 = t453 * t484 + t454 * t504 + qJD(1);
t465 = qJD(3) * t474 - t504;
t481 = -qJD(3) * t514 + t486;
t480 = (pkin(2) * t493 - pkin(8) * t494) * t488;
t502 = t454 * t486 - t480 * t484;
t426 = pkin(3) * t459 + pkin(9) * t458;
t427 = pkin(3) * t461 + pkin(9) * t460;
t501 = t464 * t426 - t427 * t465 + t503;
t500 = (-t453 * t490 - t480 * t513) * qJD(2);
t396 = pkin(4) * t434 + qJ(5) * t433;
t429 = qJD(4) * t460 + t464;
t499 = qJD(5) * t462 + t429 * t396 + t501;
t455 = pkin(3) * t479 + pkin(9) * t478;
t498 = t481 * t427 - t455 * t464 + t502;
t397 = pkin(4) * t436 + qJ(5) * t435;
t456 = qJD(4) * t478 + t481;
t497 = qJD(5) * t433 + t456 * t397 + t498;
t496 = -t426 * t481 + t465 * t455 + t500;
t428 = pkin(4) * t463 + qJ(5) * t462;
t430 = qJD(4) * t458 + t465;
t495 = qJD(5) * t435 + t430 * t428 + t496;
t469 = t490 * rSges(3,3) + (rSges(3,1) * t493 + rSges(3,2) * t494) * t488;
t468 = Icges(3,5) * t490 + (Icges(3,1) * t493 + Icges(3,4) * t494) * t488;
t467 = Icges(3,6) * t490 + (Icges(3,4) * t493 + Icges(3,2) * t494) * t488;
t466 = Icges(3,3) * t490 + (Icges(3,5) * t493 + Icges(3,6) * t494) * t488;
t451 = t479 * rSges(4,1) - t478 * rSges(4,2) - rSges(4,3) * t514;
t450 = Icges(4,1) * t479 - Icges(4,4) * t478 - Icges(4,5) * t514;
t449 = Icges(4,4) * t479 - Icges(4,2) * t478 - Icges(4,6) * t514;
t448 = Icges(4,5) * t479 - Icges(4,6) * t478 - Icges(4,3) * t514;
t445 = rSges(3,1) * t477 - rSges(3,2) * t476 + rSges(3,3) * t516;
t444 = rSges(3,1) * t475 - rSges(3,2) * t474 - rSges(3,3) * t513;
t443 = Icges(3,1) * t477 - Icges(3,4) * t476 + Icges(3,5) * t516;
t442 = Icges(3,1) * t475 - Icges(3,4) * t474 - Icges(3,5) * t513;
t441 = Icges(3,4) * t477 - Icges(3,2) * t476 + Icges(3,6) * t516;
t440 = Icges(3,4) * t475 - Icges(3,2) * t474 - Icges(3,6) * t513;
t439 = Icges(3,5) * t477 - Icges(3,6) * t476 + Icges(3,3) * t516;
t438 = Icges(3,5) * t475 - Icges(3,6) * t474 - Icges(3,3) * t513;
t423 = (-t444 * t490 - t469 * t513) * qJD(2);
t422 = (t445 * t490 - t469 * t516) * qJD(2);
t421 = rSges(5,1) * t463 - rSges(5,2) * t462 + rSges(5,3) * t478;
t420 = rSges(6,1) * t463 + rSges(6,2) * t478 + rSges(6,3) * t462;
t409 = rSges(4,1) * t461 - rSges(4,2) * t460 + rSges(4,3) * t476;
t408 = rSges(4,1) * t459 - rSges(4,2) * t458 + rSges(4,3) * t474;
t407 = Icges(4,1) * t461 - Icges(4,4) * t460 + Icges(4,5) * t476;
t406 = Icges(4,1) * t459 - Icges(4,4) * t458 + Icges(4,5) * t474;
t405 = Icges(4,4) * t461 - Icges(4,2) * t460 + Icges(4,6) * t476;
t404 = Icges(4,4) * t459 - Icges(4,2) * t458 + Icges(4,6) * t474;
t403 = Icges(4,5) * t461 - Icges(4,6) * t460 + Icges(4,3) * t476;
t402 = Icges(4,5) * t459 - Icges(4,6) * t458 + Icges(4,3) * t474;
t398 = qJD(1) + (t444 * t487 + t445 * t489) * t507;
t394 = rSges(5,1) * t436 - rSges(5,2) * t435 + rSges(5,3) * t460;
t393 = rSges(6,1) * t436 + rSges(6,2) * t460 + rSges(6,3) * t435;
t391 = rSges(5,1) * t434 - rSges(5,2) * t433 + rSges(5,3) * t458;
t390 = rSges(6,1) * t434 + rSges(6,2) * t458 + rSges(6,3) * t433;
t368 = -t408 * t481 + t451 * t465 + t500;
t367 = t409 * t481 - t451 * t464 + t502;
t366 = t408 * t464 - t409 * t465 + t503;
t365 = -t391 * t456 + t421 * t430 + t496;
t364 = t394 * t456 - t421 * t429 + t498;
t363 = t391 * t429 - t394 * t430 + t501;
t362 = t420 * t430 + (-t390 - t396) * t456 + t495;
t361 = t393 * t456 + (-t420 - t428) * t429 + t497;
t360 = t390 * t429 + (-t393 - t397) * t430 + t499;
t359 = -qJD(6) * t460 + t508 * t430 + (-t396 - t510) * t456 + t495;
t358 = -qJD(6) * t458 + t509 * t456 + (-t428 - t508) * t429 + t497;
t357 = -qJD(6) * t478 + t510 * t429 + (-t397 - t509) * t430 + t499;
t1 = -t522 * ((-t439 * t513 - t441 * t474 + t443 * t475) * t516 - (-t438 * t513 - t440 * t474 + t442 * t475) * t513 + (-t466 * t513 - t467 * t474 + t468 * t475) * t490) * t513 / 0.2e1 + t481 * ((-t403 * t514 - t478 * t405 + t479 * t407) * t464 + (-t402 * t514 - t478 * t404 + t479 * t406) * t465 + (-t448 * t514 - t478 * t449 + t479 * t450) * t481) / 0.2e1 + t464 * ((t476 * t403 - t460 * t405 + t461 * t407) * t464 + (t402 * t476 - t404 * t460 + t406 * t461) * t465 + (t448 * t476 - t449 * t460 + t450 * t461) * t481) / 0.2e1 + t465 * ((t403 * t474 - t405 * t458 + t407 * t459) * t464 + (t474 * t402 - t458 * t404 + t459 * t406) * t465 + (t448 * t474 - t449 * t458 + t450 * t459) * t481) / 0.2e1 + m(7) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(5) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(6) * (t360 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + m(4) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(3) * (t398 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + (((t439 * t516 - t441 * t476 + t443 * t477) * t516 - (t438 * t516 - t440 * t476 + t442 * t477) * t513 + (t466 * t516 - t467 * t476 + t468 * t477) * t490) * t516 + t490 * (t490 ^ 2 * t466 + (((t441 * t494 + t443 * t493) * t487 - (t440 * t494 + t442 * t493) * t489) * t488 + (-t438 * t489 + t439 * t487 + t467 * t494 + t468 * t493) * t490) * t488)) * t522 / 0.2e1 + ((t524 * t435 + t523 * t436 + t525 * t460) * t456 + (t529 * t435 + t527 * t436 + t531 * t460) * t430 + (t528 * t435 + t526 * t436 + t530 * t460) * t429) * t429 / 0.2e1 + ((t524 * t433 + t523 * t434 + t525 * t458) * t456 + (t529 * t433 + t527 * t434 + t531 * t458) * t430 + (t528 * t433 + t526 * t434 + t530 * t458) * t429) * t430 / 0.2e1 + ((t524 * t462 + t523 * t463 + t525 * t478) * t456 + (t529 * t462 + t527 * t463 + t531 * t478) * t430 + (t528 * t462 + t526 * t463 + t530 * t478) * t429) * t456 / 0.2e1;
T  = t1;
