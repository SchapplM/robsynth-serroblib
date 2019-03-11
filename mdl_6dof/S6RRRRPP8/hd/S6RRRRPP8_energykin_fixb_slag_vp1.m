% Calculate kinetic energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:17
% EndTime: 2019-03-09 21:30:19
% DurationCPUTime: 2.46s
% Computational Cost: add. (2571->280), mult. (6443->421), div. (0->0), fcn. (7925->10), ass. (0->131)
t542 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t541 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t540 = Icges(6,4) + Icges(5,5) - Icges(7,5);
t539 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t538 = Icges(6,2) + Icges(5,3) + Icges(7,3);
t537 = -Icges(5,6) + Icges(6,6) - Icges(7,6);
t536 = rSges(7,1) + pkin(5);
t535 = -rSges(7,3) - qJ(6);
t491 = sin(qJ(2));
t492 = sin(qJ(1));
t493 = cos(qJ(2));
t494 = cos(qJ(1));
t518 = cos(pkin(6));
t505 = t494 * t518;
t472 = t491 * t492 - t493 * t505;
t473 = t491 * t505 + t492 * t493;
t488 = sin(pkin(6));
t515 = t488 * t494;
t438 = Icges(3,5) * t473 - Icges(3,6) * t472 - Icges(3,3) * t515;
t506 = t492 * t518;
t474 = t494 * t491 + t493 * t506;
t475 = -t491 * t506 + t494 * t493;
t517 = t488 * t492;
t439 = Icges(3,5) * t475 - Icges(3,6) * t474 + Icges(3,3) * t517;
t534 = t488 * (t438 * t494 - t439 * t492);
t490 = sin(qJ(3));
t520 = cos(qJ(3));
t457 = t473 * t520 - t490 * t515;
t489 = sin(qJ(4));
t519 = cos(qJ(4));
t429 = t457 * t489 - t472 * t519;
t430 = t457 * t519 + t472 * t489;
t508 = t488 * t520;
t456 = t473 * t490 + t494 * t508;
t533 = t537 * t429 + t540 * t430 + t538 * t456;
t459 = t475 * t520 + t490 * t517;
t431 = t459 * t489 - t474 * t519;
t432 = t459 * t519 + t474 * t489;
t458 = t475 * t490 - t492 * t508;
t532 = t537 * t431 + t540 * t432 + t538 * t458;
t531 = t539 * t429 + t541 * t430 + t537 * t456;
t530 = t539 * t431 + t541 * t432 + t537 * t458;
t529 = t541 * t429 + t542 * t430 + t540 * t456;
t528 = t541 * t431 + t542 * t432 + t540 * t458;
t471 = t490 * t518 + t491 * t508;
t516 = t488 * t493;
t454 = t471 * t489 + t516 * t519;
t455 = t471 * t519 - t489 * t516;
t470 = t488 * t490 * t491 - t518 * t520;
t527 = t537 * t454 + t540 * t455 + t538 * t470;
t526 = t539 * t454 + t541 * t455 + t537 * t470;
t525 = t541 * t454 + t542 * t455 + t540 * t470;
t514 = rSges(7,2) * t429 + t536 * t430 + t535 * t456;
t513 = rSges(7,2) * t431 + t536 * t432 + t535 * t458;
t512 = rSges(7,2) * t454 + t536 * t455 + t535 * t470;
t450 = pkin(2) * t473 + pkin(9) * t472;
t451 = pkin(2) * t475 + pkin(9) * t474;
t509 = qJD(2) * t488;
t484 = t492 * t509;
t507 = t494 * t509;
t511 = t450 * t484 + t451 * t507;
t460 = qJD(3) * t474 + t484;
t510 = qJD(1) * (pkin(1) * t492 - pkin(8) * t515);
t485 = qJD(2) * t518 + qJD(1);
t461 = qJD(3) * t472 - t507;
t423 = pkin(3) * t457 + pkin(10) * t456;
t424 = pkin(3) * t459 + pkin(10) * t458;
t503 = t460 * t423 - t424 * t461 + t511;
t477 = -qJD(3) * t516 + t485;
t476 = (pkin(2) * t491 - pkin(9) * t493) * t488;
t478 = qJD(1) * (pkin(1) * t494 + pkin(8) * t517);
t502 = t485 * t451 - t476 * t484 + t478;
t392 = pkin(4) * t430 + qJ(5) * t429;
t425 = qJD(4) * t458 + t460;
t501 = qJD(5) * t454 + t425 * t392 + t503;
t500 = -t450 * t485 - t476 * t507 - t510;
t449 = pkin(3) * t471 + pkin(10) * t470;
t499 = t477 * t424 - t449 * t460 + t502;
t393 = pkin(4) * t432 + qJ(5) * t431;
t452 = qJD(4) * t470 + t477;
t498 = qJD(5) * t429 + t452 * t393 + t499;
t497 = -t423 * t477 + t461 * t449 + t500;
t422 = pkin(4) * t455 + qJ(5) * t454;
t426 = qJD(4) * t456 + t461;
t496 = qJD(5) * t431 + t426 * t422 + t497;
t481 = rSges(2,1) * t494 - rSges(2,2) * t492;
t480 = rSges(2,1) * t492 + rSges(2,2) * t494;
t465 = t518 * rSges(3,3) + (rSges(3,1) * t491 + rSges(3,2) * t493) * t488;
t464 = Icges(3,5) * t518 + (Icges(3,1) * t491 + Icges(3,4) * t493) * t488;
t463 = Icges(3,6) * t518 + (Icges(3,4) * t491 + Icges(3,2) * t493) * t488;
t462 = Icges(3,3) * t518 + (Icges(3,5) * t491 + Icges(3,6) * t493) * t488;
t446 = rSges(3,1) * t475 - rSges(3,2) * t474 + rSges(3,3) * t517;
t445 = rSges(3,1) * t473 - rSges(3,2) * t472 - rSges(3,3) * t515;
t443 = Icges(3,1) * t475 - Icges(3,4) * t474 + Icges(3,5) * t517;
t442 = Icges(3,1) * t473 - Icges(3,4) * t472 - Icges(3,5) * t515;
t441 = Icges(3,4) * t475 - Icges(3,2) * t474 + Icges(3,6) * t517;
t440 = Icges(3,4) * t473 - Icges(3,2) * t472 - Icges(3,6) * t515;
t437 = rSges(4,1) * t471 - rSges(4,2) * t470 - rSges(4,3) * t516;
t436 = Icges(4,1) * t471 - Icges(4,4) * t470 - Icges(4,5) * t516;
t435 = Icges(4,4) * t471 - Icges(4,2) * t470 - Icges(4,6) * t516;
t434 = Icges(4,5) * t471 - Icges(4,6) * t470 - Icges(4,3) * t516;
t419 = rSges(4,1) * t459 - rSges(4,2) * t458 + rSges(4,3) * t474;
t418 = rSges(4,1) * t457 - rSges(4,2) * t456 + rSges(4,3) * t472;
t417 = Icges(4,1) * t459 - Icges(4,4) * t458 + Icges(4,5) * t474;
t416 = Icges(4,1) * t457 - Icges(4,4) * t456 + Icges(4,5) * t472;
t415 = Icges(4,4) * t459 - Icges(4,2) * t458 + Icges(4,6) * t474;
t414 = Icges(4,4) * t457 - Icges(4,2) * t456 + Icges(4,6) * t472;
t413 = Icges(4,5) * t459 - Icges(4,6) * t458 + Icges(4,3) * t474;
t412 = Icges(4,5) * t457 - Icges(4,6) * t456 + Icges(4,3) * t472;
t411 = rSges(5,1) * t455 - rSges(5,2) * t454 + rSges(5,3) * t470;
t410 = rSges(6,1) * t455 + rSges(6,2) * t470 + rSges(6,3) * t454;
t396 = t446 * t485 - t465 * t484 + t478;
t395 = -t445 * t485 - t465 * t507 - t510;
t394 = (t445 * t492 + t446 * t494) * t509;
t390 = rSges(5,1) * t432 - rSges(5,2) * t431 + rSges(5,3) * t458;
t389 = rSges(6,1) * t432 + rSges(6,2) * t458 + rSges(6,3) * t431;
t387 = rSges(5,1) * t430 - rSges(5,2) * t429 + rSges(5,3) * t456;
t386 = rSges(6,1) * t430 + rSges(6,2) * t456 + rSges(6,3) * t429;
t364 = t419 * t477 - t437 * t460 + t502;
t363 = -t418 * t477 + t437 * t461 + t500;
t362 = t418 * t460 - t419 * t461 + t511;
t361 = t390 * t452 - t411 * t425 + t499;
t360 = -t387 * t452 + t411 * t426 + t497;
t359 = t387 * t425 - t390 * t426 + t503;
t358 = t389 * t452 + (-t410 - t422) * t425 + t498;
t357 = t410 * t426 + (-t386 - t392) * t452 + t496;
t356 = t386 * t425 + (-t389 - t393) * t426 + t501;
t355 = -qJD(6) * t456 + t513 * t452 + (-t422 - t512) * t425 + t498;
t354 = -qJD(6) * t458 + t512 * t426 + (-t392 - t514) * t452 + t496;
t353 = -qJD(6) * t470 + t514 * t425 + (-t393 - t513) * t426 + t501;
t1 = t461 * ((t413 * t472 - t415 * t456 + t417 * t457) * t460 + (t412 * t472 - t414 * t456 + t416 * t457) * t461 + (t434 * t472 - t435 * t456 + t436 * t457) * t477) / 0.2e1 + t477 * ((-t413 * t516 - t415 * t470 + t417 * t471) * t460 + (-t412 * t516 - t414 * t470 + t416 * t471) * t461 + (-t434 * t516 - t435 * t470 + t436 * t471) * t477) / 0.2e1 + t460 * ((t413 * t474 - t415 * t458 + t417 * t459) * t460 + (t412 * t474 - t414 * t458 + t416 * t459) * t461 + (t434 * t474 - t435 * t458 + t436 * t459) * t477) / 0.2e1 + t485 * ((t518 * t439 + (t441 * t493 + t443 * t491) * t488) * t484 - (t518 * t438 + (t440 * t493 + t442 * t491) * t488) * t507 + (t518 * t462 + (t463 * t493 + t464 * t491) * t488) * t485) / 0.2e1 + m(6) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(7) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(5) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(4) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(3) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 - ((-t462 * t515 - t463 * t472 + t464 * t473) * t485 + ((-t441 * t472 + t443 * t473) * t492 + (t472 * t440 - t473 * t442 + t534) * t494) * t509) * t507 / 0.2e1 + ((t462 * t517 - t463 * t474 + t464 * t475) * t485 + (-(-t440 * t474 + t442 * t475) * t494 + (-t474 * t441 + t475 * t443 - t534) * t492) * t509) * t484 / 0.2e1 + (m(2) * (t480 ^ 2 + t481 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t431 * t526 + t432 * t525 + t458 * t527) * t452 + (t431 * t531 + t432 * t529 + t458 * t533) * t426 + (t530 * t431 + t528 * t432 + t532 * t458) * t425) * t425 / 0.2e1 + ((t429 * t526 + t430 * t525 + t456 * t527) * t452 + (t531 * t429 + t529 * t430 + t533 * t456) * t426 + (t429 * t530 + t430 * t528 + t456 * t532) * t425) * t426 / 0.2e1 + ((t526 * t454 + t525 * t455 + t527 * t470) * t452 + (t454 * t531 + t455 * t529 + t470 * t533) * t426 + (t454 * t530 + t455 * t528 + t470 * t532) * t425) * t452 / 0.2e1;
T  = t1;
