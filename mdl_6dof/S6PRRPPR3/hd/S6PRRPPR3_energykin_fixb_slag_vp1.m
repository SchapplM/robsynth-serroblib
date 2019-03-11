% Calculate kinetic energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:08
% EndTime: 2019-03-08 21:08:10
% DurationCPUTime: 2.32s
% Computational Cost: add. (1932->282), mult. (4947->422), div. (0->0), fcn. (5925->10), ass. (0->131)
t532 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t531 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t530 = -Icges(4,4) - Icges(6,4) + Icges(5,5);
t529 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t528 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t527 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t483 = sin(pkin(10));
t485 = cos(pkin(10));
t489 = cos(qJ(2));
t487 = sin(qJ(2));
t512 = cos(pkin(6));
t499 = t487 * t512;
t469 = t483 * t489 + t485 * t499;
t484 = sin(pkin(6));
t514 = cos(qJ(3));
t503 = t484 * t514;
t513 = sin(qJ(3));
t450 = t469 * t513 + t485 * t503;
t502 = t484 * t513;
t451 = t469 * t514 - t485 * t502;
t498 = t489 * t512;
t468 = t483 * t487 - t485 * t498;
t526 = t528 * t450 - t529 * t451 - t527 * t468;
t471 = -t483 * t499 + t485 * t489;
t452 = t471 * t513 - t483 * t503;
t453 = t471 * t514 + t483 * t502;
t470 = t483 * t498 + t485 * t487;
t525 = t528 * t452 - t529 * t453 - t527 * t470;
t524 = t531 * t450 + t530 * t451 - t528 * t468;
t523 = t531 * t452 + t530 * t453 - t528 * t470;
t522 = t530 * t450 + t532 * t451 + t529 * t468;
t521 = t530 * t452 + t532 * t453 + t529 * t470;
t472 = t487 * t502 - t512 * t514;
t473 = t487 * t503 + t512 * t513;
t509 = t484 * t489;
t520 = t528 * t472 - t529 * t473 + t527 * t509;
t519 = t531 * t472 + t530 * t473 + t528 * t509;
t518 = t530 * t472 + t532 * t473 - t529 * t509;
t511 = t483 * t484;
t510 = t484 * t485;
t407 = pkin(3) * t451 + qJ(4) * t450;
t418 = pkin(4) * t451 - qJ(5) * t468;
t508 = -t407 - t418;
t408 = pkin(3) * t453 + qJ(4) * t452;
t419 = pkin(4) * t453 - qJ(5) * t470;
t507 = -t408 - t419;
t445 = pkin(3) * t473 + qJ(4) * t472;
t458 = t473 * pkin(4) + qJ(5) * t509;
t506 = -t445 - t458;
t505 = qJD(2) * t484;
t479 = t483 * t505;
t456 = qJD(3) * t470 + t479;
t501 = t485 * t505;
t443 = pkin(2) * t469 + pkin(8) * t468;
t444 = pkin(2) * t471 + pkin(8) * t470;
t500 = t443 * t479 + t444 * t501 + qJD(1);
t482 = qJD(2) * t512;
t457 = qJD(3) * t468 - t501;
t475 = -qJD(3) * t509 + t482;
t497 = qJD(4) * t472 + t456 * t407 + t500;
t474 = (pkin(2) * t487 - pkin(8) * t489) * t484;
t496 = t444 * t482 - t474 * t479;
t495 = qJD(5) * t509 + t456 * t418 + t497;
t494 = qJD(4) * t450 + t475 * t408 + t496;
t493 = (-t443 * t512 - t474 * t510) * qJD(2);
t492 = -qJD(5) * t468 + t475 * t419 + t494;
t491 = qJD(4) * t452 + t457 * t445 + t493;
t490 = -qJD(5) * t470 + t457 * t458 + t491;
t488 = cos(qJ(6));
t486 = sin(qJ(6));
t462 = t512 * rSges(3,3) + (rSges(3,1) * t487 + rSges(3,2) * t489) * t484;
t461 = Icges(3,5) * t512 + (Icges(3,1) * t487 + Icges(3,4) * t489) * t484;
t460 = Icges(3,6) * t512 + (Icges(3,4) * t487 + Icges(3,2) * t489) * t484;
t459 = Icges(3,3) * t512 + (Icges(3,5) * t487 + Icges(3,6) * t489) * t484;
t455 = t472 * t488 + t486 * t509;
t454 = -t472 * t486 + t488 * t509;
t449 = qJD(6) * t473 + t475;
t446 = pkin(5) * t472 + pkin(9) * t473;
t441 = t473 * rSges(4,1) - t472 * rSges(4,2) - rSges(4,3) * t509;
t440 = t473 * rSges(5,1) - rSges(5,2) * t509 + t472 * rSges(5,3);
t439 = t472 * rSges(6,1) - t473 * rSges(6,2) + rSges(6,3) * t509;
t427 = rSges(3,1) * t471 - rSges(3,2) * t470 + rSges(3,3) * t511;
t426 = rSges(3,1) * t469 - rSges(3,2) * t468 - rSges(3,3) * t510;
t425 = Icges(3,1) * t471 - Icges(3,4) * t470 + Icges(3,5) * t511;
t424 = Icges(3,1) * t469 - Icges(3,4) * t468 - Icges(3,5) * t510;
t423 = Icges(3,4) * t471 - Icges(3,2) * t470 + Icges(3,6) * t511;
t422 = Icges(3,4) * t469 - Icges(3,2) * t468 - Icges(3,6) * t510;
t421 = Icges(3,5) * t471 - Icges(3,6) * t470 + Icges(3,3) * t511;
t420 = Icges(3,5) * t469 - Icges(3,6) * t468 - Icges(3,3) * t510;
t417 = t452 * t488 - t470 * t486;
t416 = -t452 * t486 - t470 * t488;
t415 = t450 * t488 - t468 * t486;
t414 = -t450 * t486 - t468 * t488;
t412 = qJD(6) * t451 + t457;
t411 = qJD(6) * t453 + t456;
t410 = pkin(5) * t452 + pkin(9) * t453;
t409 = pkin(5) * t450 + pkin(9) * t451;
t402 = (-t426 * t512 - t462 * t510) * qJD(2);
t401 = (t427 * t512 - t462 * t511) * qJD(2);
t400 = rSges(7,1) * t455 + rSges(7,2) * t454 + rSges(7,3) * t473;
t399 = Icges(7,1) * t455 + Icges(7,4) * t454 + Icges(7,5) * t473;
t398 = Icges(7,4) * t455 + Icges(7,2) * t454 + Icges(7,6) * t473;
t397 = Icges(7,5) * t455 + Icges(7,6) * t454 + Icges(7,3) * t473;
t396 = rSges(4,1) * t453 - rSges(4,2) * t452 + rSges(4,3) * t470;
t395 = rSges(5,1) * t453 + rSges(5,2) * t470 + rSges(5,3) * t452;
t394 = rSges(6,1) * t452 - rSges(6,2) * t453 - rSges(6,3) * t470;
t393 = rSges(4,1) * t451 - rSges(4,2) * t450 + rSges(4,3) * t468;
t392 = rSges(5,1) * t451 + rSges(5,2) * t468 + rSges(5,3) * t450;
t391 = rSges(6,1) * t450 - rSges(6,2) * t451 - rSges(6,3) * t468;
t371 = qJD(1) + (t426 * t483 + t427 * t485) * t505;
t370 = rSges(7,1) * t417 + rSges(7,2) * t416 + rSges(7,3) * t453;
t369 = rSges(7,1) * t415 + rSges(7,2) * t414 + rSges(7,3) * t451;
t368 = Icges(7,1) * t417 + Icges(7,4) * t416 + Icges(7,5) * t453;
t367 = Icges(7,1) * t415 + Icges(7,4) * t414 + Icges(7,5) * t451;
t366 = Icges(7,4) * t417 + Icges(7,2) * t416 + Icges(7,6) * t453;
t365 = Icges(7,4) * t415 + Icges(7,2) * t414 + Icges(7,6) * t451;
t364 = Icges(7,5) * t417 + Icges(7,6) * t416 + Icges(7,3) * t453;
t363 = Icges(7,5) * t415 + Icges(7,6) * t414 + Icges(7,3) * t451;
t362 = -t475 * t393 + t457 * t441 + t493;
t361 = t396 * t475 - t441 * t456 + t496;
t360 = t393 * t456 - t396 * t457 + t500;
t359 = t457 * t440 + (-t392 - t407) * t475 + t491;
t358 = t395 * t475 + (-t440 - t445) * t456 + t494;
t357 = t392 * t456 + (-t395 - t408) * t457 + t497;
t356 = t457 * t439 + (-t391 + t508) * t475 + t490;
t355 = t394 * t475 + (-t439 + t506) * t456 + t492;
t354 = t391 * t456 + (-t394 + t507) * t457 + t495;
t353 = -t449 * t369 + t412 * t400 + t457 * t446 + (-t409 + t508) * t475 + t490;
t352 = t370 * t449 - t400 * t411 + t410 * t475 + (-t446 + t506) * t456 + t492;
t351 = t369 * t411 - t370 * t412 + t409 * t456 + (-t410 + t507) * t457 + t495;
t1 = t411 * ((t453 * t364 + t416 * t366 + t417 * t368) * t411 + (t363 * t453 + t365 * t416 + t367 * t417) * t412 + (t397 * t453 + t398 * t416 + t399 * t417) * t449) / 0.2e1 + t412 * ((t364 * t451 + t366 * t414 + t368 * t415) * t411 + (t451 * t363 + t414 * t365 + t415 * t367) * t412 + (t397 * t451 + t398 * t414 + t399 * t415) * t449) / 0.2e1 + t449 * ((t364 * t473 + t366 * t454 + t368 * t455) * t411 + (t363 * t473 + t365 * t454 + t367 * t455) * t412 + (t397 * t473 + t398 * t454 + t399 * t455) * t449) / 0.2e1 + m(3) * (t371 ^ 2 + t401 ^ 2 + t402 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(4) * (t360 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + m(5) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(6) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(7) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + ((t512 * t421 + (t423 * t489 + t425 * t487) * t484) * t479 - (t512 * t420 + (t422 * t489 + t424 * t487) * t484) * t501 + (t512 * t459 + (t460 * t489 + t461 * t487) * t484) * t482) * t482 / 0.2e1 + (-t485 * ((-t421 * t510 - t423 * t468 + t425 * t469) * t511 - (-t420 * t510 - t422 * t468 + t424 * t469) * t510 + (-t459 * t510 - t460 * t468 + t461 * t469) * t512) / 0.2e1 + t483 * ((t421 * t511 - t423 * t470 + t425 * t471) * t511 - (t420 * t511 - t422 * t470 + t424 * t471) * t510 + (t459 * t511 - t460 * t470 + t461 * t471) * t512) / 0.2e1) * t484 * qJD(2) ^ 2 + ((t452 * t519 + t453 * t518 - t470 * t520) * t475 + (t524 * t452 + t522 * t453 - t470 * t526) * t457 + (t452 * t523 + t453 * t521 - t470 * t525) * t456) * t456 / 0.2e1 + ((t450 * t519 + t451 * t518 - t468 * t520) * t475 + (t524 * t450 + t522 * t451 - t468 * t526) * t457 + (t450 * t523 + t451 * t521 - t468 * t525) * t456) * t457 / 0.2e1 + ((t472 * t519 + t473 * t518 + t509 * t520) * t475 + (t524 * t472 + t522 * t473 + t509 * t526) * t457 + (t472 * t523 + t473 * t521 + t509 * t525) * t456) * t475 / 0.2e1;
T  = t1;
