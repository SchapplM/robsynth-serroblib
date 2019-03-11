% Calculate kinetic energy for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:43
% EndTime: 2019-03-09 16:47:46
% DurationCPUTime: 3.23s
% Computational Cost: add. (1990->303), mult. (2604->459), div. (0->0), fcn. (2627->10), ass. (0->151)
t550 = Icges(6,1) + Icges(7,1);
t549 = -Icges(6,4) + Icges(7,5);
t548 = Icges(7,4) + Icges(6,5);
t547 = Icges(6,2) + Icges(7,3);
t546 = -Icges(7,6) + Icges(6,6);
t545 = -Icges(5,3) - Icges(4,3);
t544 = -Icges(6,3) - Icges(7,2);
t543 = rSges(7,1) + pkin(5);
t542 = rSges(7,3) + qJ(6);
t474 = qJ(3) + pkin(10);
t470 = qJ(5) + t474;
t465 = sin(t470);
t466 = cos(t470);
t481 = cos(qJ(1));
t478 = sin(qJ(1));
t480 = cos(qJ(2));
t513 = t480 * t478;
t413 = t465 * t513 + t466 * t481;
t414 = -t465 * t481 + t466 * t513;
t477 = sin(qJ(2));
t515 = t477 * t478;
t541 = t547 * t413 + t549 * t414 - t546 * t515;
t512 = t480 * t481;
t415 = t465 * t512 - t478 * t466;
t416 = t465 * t478 + t466 * t512;
t514 = t477 * t481;
t540 = t547 * t415 + t549 * t416 - t546 * t514;
t539 = -t546 * t413 + t548 * t414 - t544 * t515;
t538 = -t546 * t415 + t548 * t416 - t544 * t514;
t537 = t549 * t413 + t550 * t414 + t548 * t515;
t536 = t549 * t415 + t550 * t416 + t548 * t514;
t468 = sin(t474);
t469 = cos(t474);
t420 = -t468 * t513 - t469 * t481;
t421 = -t468 * t481 + t469 * t513;
t476 = sin(qJ(3));
t479 = cos(qJ(3));
t439 = -t476 * t513 - t479 * t481;
t516 = t476 * t481;
t440 = t479 * t513 - t516;
t535 = Icges(4,5) * t440 + Icges(5,5) * t421 + Icges(4,6) * t439 + Icges(5,6) * t420 - t545 * t515;
t422 = -t468 * t512 + t469 * t478;
t423 = t468 * t478 + t469 * t512;
t441 = -t476 * t512 + t478 * t479;
t517 = t476 * t478;
t442 = t479 * t512 + t517;
t534 = Icges(4,5) * t442 + Icges(5,5) * t423 + Icges(4,6) * t441 + Icges(5,6) * t422 - t545 * t514;
t533 = t546 * t480 + (t547 * t465 + t549 * t466) * t477;
t532 = t544 * t480 + (-t546 * t465 + t548 * t466) * t477;
t531 = -t548 * t480 + (t549 * t465 + t550 * t466) * t477;
t530 = t545 * t480 + (Icges(4,5) * t479 + Icges(5,5) * t469 - Icges(4,6) * t476 - Icges(5,6) * t468) * t477;
t521 = t479 * pkin(3);
t519 = Icges(3,4) * t477;
t518 = Icges(3,4) * t480;
t511 = rSges(7,2) * t515 + t542 * t413 + t543 * t414;
t510 = rSges(7,2) * t514 + t542 * t415 + t543 * t416;
t509 = -rSges(7,2) * t480 + (t542 * t465 + t543 * t466) * t477;
t500 = pkin(2) * t480 + pkin(8) * t477;
t443 = t500 * t478;
t444 = t500 * t481;
t471 = qJD(2) * t478;
t505 = qJD(2) * t481;
t508 = t443 * t471 + t444 * t505;
t507 = pkin(4) * t469;
t504 = qJD(3) * t477;
t445 = t481 * t504 + t471;
t503 = qJD(4) * t477;
t502 = qJD(5) * t477;
t501 = pkin(4) * t468;
t446 = t478 * t504 - t505;
t499 = rSges(3,1) * t480 - rSges(3,2) * t477;
t498 = Icges(3,1) * t480 - t519;
t497 = -Icges(3,2) * t477 + t518;
t496 = Icges(3,5) * t480 - Icges(3,6) * t477;
t428 = -Icges(3,6) * t481 + t478 * t497;
t431 = -Icges(3,5) * t481 + t478 * t498;
t495 = t428 * t477 - t431 * t480;
t429 = Icges(3,6) * t478 + t481 * t497;
t432 = Icges(3,5) * t478 + t481 * t498;
t494 = -t429 * t477 + t432 * t480;
t452 = Icges(3,2) * t480 + t519;
t453 = Icges(3,1) * t477 + t518;
t493 = -t452 * t477 + t453 * t480;
t450 = qJD(1) * (pkin(1) * t481 + pkin(7) * t478);
t458 = pkin(2) * t477 - pkin(8) * t480;
t492 = qJD(1) * t444 - t458 * t471 + t450;
t490 = qJ(4) * t477 + t480 * t521;
t394 = -pkin(3) * t516 + t478 * t490;
t491 = -qJD(4) * t480 + t445 * t394 + t508;
t395 = pkin(3) * t517 + t481 * t490;
t464 = -qJD(3) * t480 + qJD(1);
t489 = t464 * t395 + t478 * t503 + t492;
t459 = pkin(1) * t478 - pkin(7) * t481;
t488 = (-t443 - t459) * qJD(1) - t458 * t505;
t487 = pkin(9) * t477 + t480 * t507;
t408 = -qJ(4) * t480 + t477 * t521;
t486 = t446 * t408 + t481 * t503 + t488;
t354 = t478 * t487 - t481 * t501;
t355 = t478 * t501 + t481 * t487;
t485 = t445 * t354 + (-t355 - t395) * t446 + t491;
t398 = -pkin(9) * t480 + t477 * t507;
t484 = t464 * t355 + (-t398 - t408) * t445 + t489;
t483 = t446 * t398 + (-t354 - t394) * t464 + t486;
t456 = rSges(2,1) * t481 - rSges(2,2) * t478;
t455 = rSges(2,1) * t478 + rSges(2,2) * t481;
t454 = rSges(3,1) * t477 + rSges(3,2) * t480;
t451 = Icges(3,5) * t477 + Icges(3,6) * t480;
t448 = qJD(1) + (-qJD(3) - qJD(5)) * t480;
t435 = rSges(3,3) * t478 + t481 * t499;
t434 = -rSges(3,3) * t481 + t478 * t499;
t433 = -rSges(4,3) * t480 + (rSges(4,1) * t479 - rSges(4,2) * t476) * t477;
t430 = -Icges(4,5) * t480 + (Icges(4,1) * t479 - Icges(4,4) * t476) * t477;
t427 = -Icges(4,6) * t480 + (Icges(4,4) * t479 - Icges(4,2) * t476) * t477;
t426 = Icges(3,3) * t478 + t481 * t496;
t425 = -Icges(3,3) * t481 + t478 * t496;
t419 = t478 * t502 + t446;
t418 = t481 * t502 + t445;
t412 = -rSges(5,3) * t480 + (rSges(5,1) * t469 - rSges(5,2) * t468) * t477;
t411 = -Icges(5,5) * t480 + (Icges(5,1) * t469 - Icges(5,4) * t468) * t477;
t410 = -Icges(5,6) * t480 + (Icges(5,4) * t469 - Icges(5,2) * t468) * t477;
t407 = -rSges(6,3) * t480 + (rSges(6,1) * t466 - rSges(6,2) * t465) * t477;
t397 = rSges(4,1) * t442 + rSges(4,2) * t441 + rSges(4,3) * t514;
t396 = rSges(4,1) * t440 + rSges(4,2) * t439 + rSges(4,3) * t515;
t393 = Icges(4,1) * t442 + Icges(4,4) * t441 + Icges(4,5) * t514;
t392 = Icges(4,1) * t440 + Icges(4,4) * t439 + Icges(4,5) * t515;
t391 = Icges(4,4) * t442 + Icges(4,2) * t441 + Icges(4,6) * t514;
t390 = Icges(4,4) * t440 + Icges(4,2) * t439 + Icges(4,6) * t515;
t387 = qJD(1) * t435 - t454 * t471 + t450;
t386 = -t454 * t505 + (-t434 - t459) * qJD(1);
t385 = (t434 * t478 + t435 * t481) * qJD(2);
t381 = rSges(5,1) * t423 + rSges(5,2) * t422 + rSges(5,3) * t514;
t380 = rSges(5,1) * t421 + rSges(5,2) * t420 + rSges(5,3) * t515;
t379 = Icges(5,1) * t423 + Icges(5,4) * t422 + Icges(5,5) * t514;
t378 = Icges(5,1) * t421 + Icges(5,4) * t420 + Icges(5,5) * t515;
t377 = Icges(5,4) * t423 + Icges(5,2) * t422 + Icges(5,6) * t514;
t376 = Icges(5,4) * t421 + Icges(5,2) * t420 + Icges(5,6) * t515;
t371 = rSges(6,1) * t416 - rSges(6,2) * t415 + rSges(6,3) * t514;
t369 = rSges(6,1) * t414 - rSges(6,2) * t413 + rSges(6,3) * t515;
t351 = t397 * t464 - t433 * t445 + t492;
t350 = -t396 * t464 + t433 * t446 + t488;
t349 = t396 * t445 - t397 * t446 + t508;
t348 = t381 * t464 + (-t408 - t412) * t445 + t489;
t347 = t412 * t446 + (-t380 - t394) * t464 + t486;
t346 = t380 * t445 + (-t381 - t395) * t446 + t491;
t345 = t371 * t448 - t407 * t418 + t484;
t344 = -t369 * t448 + t407 * t419 + t483;
t343 = t369 * t418 - t371 * t419 + t485;
t342 = qJD(6) * t413 - t418 * t509 + t448 * t510 + t484;
t341 = qJD(6) * t415 + t419 * t509 - t448 * t511 + t483;
t340 = qJD(6) * t465 * t477 + t418 * t511 - t419 * t510 + t485;
t1 = m(6) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(5) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(4) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(3) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + qJD(1) * ((t452 * t480 + t453 * t477) * qJD(1) + ((t429 * t480 + t432 * t477) * t478 - (t428 * t480 + t431 * t477) * t481) * qJD(2)) / 0.2e1 - ((-t481 * t451 + t478 * t493) * qJD(1) + (t481 ^ 2 * t425 + (t494 * t478 + (-t426 + t495) * t481) * t478) * qJD(2)) * t505 / 0.2e1 + ((t478 * t451 + t481 * t493) * qJD(1) + (t478 ^ 2 * t426 + (t495 * t481 + (-t425 + t494) * t478) * t481) * qJD(2)) * t471 / 0.2e1 + m(7) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + ((t533 * t415 + t531 * t416 + t532 * t514) * t448 + (t541 * t415 + t537 * t416 + t539 * t514) * t419 + (t540 * t415 + t536 * t416 + t538 * t514) * t418) * t418 / 0.2e1 + ((t533 * t413 + t531 * t414 + t532 * t515) * t448 + (t541 * t413 + t537 * t414 + t539 * t515) * t419 + (t540 * t413 + t536 * t414 + t538 * t515) * t418) * t419 / 0.2e1 + ((t410 * t422 + t411 * t423 + t427 * t441 + t430 * t442 + t530 * t514) * t464 + (t376 * t422 + t378 * t423 + t390 * t441 + t392 * t442 + t535 * t514) * t446 + (t422 * t377 + t423 * t379 + t441 * t391 + t442 * t393 + t534 * t514) * t445) * t445 / 0.2e1 + ((t410 * t420 + t411 * t421 + t427 * t439 + t430 * t440 + t530 * t515) * t464 + (t420 * t376 + t421 * t378 + t439 * t390 + t440 * t392 + t535 * t515) * t446 + (t377 * t420 + t379 * t421 + t391 * t439 + t393 * t440 + t534 * t515) * t445) * t446 / 0.2e1 + ((-t538 * t418 - t539 * t419 - t532 * t448) * t480 + ((t533 * t465 + t531 * t466) * t448 + (t541 * t465 + t537 * t466) * t419 + (t540 * t465 + t536 * t466) * t418) * t477) * t448 / 0.2e1 + ((-t534 * t445 - t535 * t446 - t530 * t464) * t480 + ((-t410 * t468 + t411 * t469 - t427 * t476 + t430 * t479) * t464 + (-t376 * t468 + t378 * t469 - t390 * t476 + t392 * t479) * t446 + (-t377 * t468 + t379 * t469 - t391 * t476 + t393 * t479) * t445) * t477) * t464 / 0.2e1 + (m(2) * (t455 ^ 2 + t456 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
