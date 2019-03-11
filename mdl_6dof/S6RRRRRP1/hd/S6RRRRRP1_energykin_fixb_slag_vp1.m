% Calculate kinetic energy for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:32
% EndTime: 2019-03-10 00:55:34
% DurationCPUTime: 2.47s
% Computational Cost: add. (1975->268), mult. (1970->416), div. (0->0), fcn. (1871->10), ass. (0->157)
t549 = Icges(6,1) + Icges(7,1);
t548 = Icges(6,4) + Icges(7,4);
t547 = -Icges(7,5) - Icges(6,5);
t546 = Icges(6,2) + Icges(7,2);
t545 = -Icges(7,6) - Icges(6,6);
t544 = -Icges(7,3) - Icges(6,3);
t458 = qJ(2) + qJ(3);
t455 = qJ(4) + t458;
t447 = cos(t455);
t463 = cos(qJ(5));
t465 = cos(qJ(1));
t510 = t463 * t465;
t460 = sin(qJ(5));
t462 = sin(qJ(1));
t513 = t460 * t462;
t418 = -t447 * t513 - t510;
t511 = t462 * t463;
t512 = t460 * t465;
t419 = t447 * t511 - t512;
t446 = sin(t455);
t515 = t446 * t462;
t543 = -t545 * t418 - t547 * t419 - t544 * t515;
t420 = -t447 * t512 + t511;
t421 = t447 * t510 + t513;
t514 = t446 * t465;
t542 = -t545 * t420 - t547 * t421 - t544 * t514;
t541 = t546 * t418 + t548 * t419 - t545 * t515;
t540 = t546 * t420 + t548 * t421 - t545 * t514;
t539 = t548 * t418 + t549 * t419 - t547 * t515;
t538 = t548 * t420 + t549 * t421 - t547 * t514;
t454 = cos(t458);
t503 = pkin(3) * t454;
t372 = pkin(9) * t462 + t465 * t503;
t452 = qJD(2) * t462;
t436 = qJD(3) * t462 + t452;
t453 = sin(t458);
t528 = pkin(3) * t453;
t537 = qJD(1) * t372 - t436 * t528;
t536 = t544 * t447 + (t545 * t460 - t547 * t463) * t446;
t535 = t545 * t447 + (-t546 * t460 + t548 * t463) * t446;
t534 = t547 * t447 + (-t548 * t460 + t549 * t463) * t446;
t494 = pkin(4) * t447 + pkin(10) * t446;
t409 = t494 * t465;
t427 = pkin(4) * t446 - pkin(10) * t447;
t428 = qJD(4) * t462 + t436;
t533 = qJD(1) * t409 - t427 * t428;
t464 = cos(qJ(2));
t526 = t464 * pkin(2);
t525 = pkin(5) * t463;
t461 = sin(qJ(2));
t522 = Icges(3,4) * t461;
t521 = Icges(3,4) * t464;
t520 = Icges(4,4) * t453;
t519 = Icges(4,4) * t454;
t518 = Icges(5,4) * t446;
t517 = Icges(5,4) * t447;
t475 = qJ(6) * t446 + t447 * t525;
t509 = rSges(7,1) * t419 + rSges(7,2) * t418 + rSges(7,3) * t515 - pkin(5) * t512 + t462 * t475;
t508 = rSges(7,1) * t421 + rSges(7,2) * t420 + rSges(7,3) * t514 + pkin(5) * t513 + t465 * t475;
t507 = (-qJ(6) - rSges(7,3)) * t447 + (rSges(7,1) * t463 - rSges(7,2) * t460 + t525) * t446;
t397 = -pkin(8) * t465 + t462 * t526;
t398 = pkin(8) * t462 + t465 * t526;
t501 = qJD(2) * t465;
t506 = t397 * t452 + t398 * t501;
t435 = qJD(1) * (pkin(1) * t465 + pkin(7) * t462);
t505 = qJD(1) * t398 + t435;
t445 = pkin(1) * t462 - pkin(7) * t465;
t504 = -t397 - t445;
t500 = qJD(5) * t446;
t499 = -qJD(2) - qJD(3);
t497 = pkin(2) * qJD(2) * t461;
t371 = -pkin(9) * t465 + t462 * t503;
t496 = -t371 + t504;
t495 = t465 * t497;
t493 = rSges(3,1) * t464 - rSges(3,2) * t461;
t492 = rSges(4,1) * t454 - rSges(4,2) * t453;
t491 = rSges(5,1) * t447 - rSges(5,2) * t446;
t429 = (-qJD(4) + t499) * t465;
t490 = Icges(3,1) * t464 - t522;
t489 = Icges(4,1) * t454 - t520;
t488 = Icges(5,1) * t447 - t518;
t487 = -Icges(3,2) * t461 + t521;
t486 = -Icges(4,2) * t453 + t519;
t485 = -Icges(5,2) * t446 + t517;
t484 = Icges(3,5) * t464 - Icges(3,6) * t461;
t483 = Icges(4,5) * t454 - Icges(4,6) * t453;
t482 = Icges(5,5) * t447 - Icges(5,6) * t446;
t412 = -Icges(3,6) * t465 + t462 * t487;
t414 = -Icges(3,5) * t465 + t462 * t490;
t481 = t412 * t461 - t414 * t464;
t413 = Icges(3,6) * t462 + t465 * t487;
t415 = Icges(3,5) * t462 + t465 * t490;
t480 = -t413 * t461 + t415 * t464;
t440 = Icges(3,2) * t464 + t522;
t441 = Icges(3,1) * t461 + t521;
t479 = -t440 * t461 + t441 * t464;
t437 = t499 * t465;
t478 = t436 * t371 - t372 * t437 + t506;
t477 = qJD(6) * t446 - t497;
t476 = -t462 * t497 + t505;
t474 = (Icges(5,5) * t446 + Icges(5,6) * t447) * qJD(1) + (-Icges(5,3) * t465 + t462 * t482) * t429 + (Icges(5,3) * t462 + t465 * t482) * t428;
t473 = (Icges(4,5) * t453 + Icges(4,6) * t454) * qJD(1) + (-Icges(4,3) * t465 + t462 * t483) * t437 + (Icges(4,3) * t462 + t465 * t483) * t436;
t408 = t494 * t462;
t422 = t437 * t528;
t472 = t429 * t427 + t422 + (-t408 + t496) * qJD(1);
t471 = t428 * t408 - t409 * t429 + t478;
t470 = t476 + t537;
t389 = -Icges(5,6) * t465 + t462 * t485;
t390 = Icges(5,6) * t462 + t465 * t485;
t391 = -Icges(5,5) * t465 + t462 * t488;
t392 = Icges(5,5) * t462 + t465 * t488;
t424 = Icges(5,2) * t447 + t518;
t425 = Icges(5,1) * t446 + t517;
t469 = (-t390 * t446 + t392 * t447) * t428 + (-t389 * t446 + t391 * t447) * t429 + (-t424 * t446 + t425 * t447) * qJD(1);
t401 = -Icges(4,6) * t465 + t462 * t486;
t402 = Icges(4,6) * t462 + t465 * t486;
t403 = -Icges(4,5) * t465 + t462 * t489;
t404 = Icges(4,5) * t462 + t465 * t489;
t431 = Icges(4,2) * t454 + t520;
t432 = Icges(4,1) * t453 + t519;
t468 = (-t402 * t453 + t404 * t454) * t436 + (-t401 * t453 + t403 * t454) * t437 + (-t431 * t453 + t432 * t454) * qJD(1);
t444 = rSges(2,1) * t465 - rSges(2,2) * t462;
t443 = rSges(2,1) * t462 + rSges(2,2) * t465;
t442 = rSges(3,1) * t461 + rSges(3,2) * t464;
t439 = Icges(3,5) * t461 + Icges(3,6) * t464;
t438 = -qJD(5) * t447 + qJD(1);
t433 = rSges(4,1) * t453 + rSges(4,2) * t454;
t426 = rSges(5,1) * t446 + rSges(5,2) * t447;
t417 = rSges(3,3) * t462 + t465 * t493;
t416 = -rSges(3,3) * t465 + t462 * t493;
t411 = Icges(3,3) * t462 + t465 * t484;
t410 = -Icges(3,3) * t465 + t462 * t484;
t406 = rSges(4,3) * t462 + t465 * t492;
t405 = -rSges(4,3) * t465 + t462 * t492;
t396 = rSges(5,3) * t462 + t465 * t491;
t395 = -rSges(5,3) * t465 + t462 * t491;
t394 = t462 * t500 + t429;
t393 = t465 * t500 + t428;
t385 = -rSges(6,3) * t447 + (rSges(6,1) * t463 - rSges(6,2) * t460) * t446;
t369 = qJD(1) * t417 - t442 * t452 + t435;
t368 = -t442 * t501 + (-t416 - t445) * qJD(1);
t366 = (t416 * t462 + t417 * t465) * qJD(2);
t365 = rSges(6,1) * t421 + rSges(6,2) * t420 + rSges(6,3) * t514;
t363 = rSges(6,1) * t419 + rSges(6,2) * t418 + rSges(6,3) * t515;
t347 = qJD(1) * t406 - t433 * t436 + t476;
t346 = -t495 + t433 * t437 + (-t405 + t504) * qJD(1);
t345 = t405 * t436 - t406 * t437 + t506;
t344 = qJD(1) * t396 - t426 * t428 + t470;
t343 = -t495 + t426 * t429 + t422 + (-t395 + t496) * qJD(1);
t342 = t395 * t428 - t396 * t429 + t478;
t341 = t365 * t438 - t385 * t393 + t470 + t533;
t340 = -t363 * t438 + t385 * t394 + t472 - t495;
t339 = t363 * t393 - t365 * t394 + t471;
t338 = -t393 * t507 + t438 * t508 + t462 * t477 + t505 + t533 + t537;
t337 = t394 * t507 - t438 * t509 + t465 * t477 + t472;
t336 = -qJD(6) * t447 + t393 * t509 - t394 * t508 + t471;
t1 = ((t462 * t439 + t465 * t479) * qJD(1) + (t462 ^ 2 * t411 + (t481 * t465 + (-t410 + t480) * t462) * t465) * qJD(2)) * t452 / 0.2e1 - ((-t465 * t439 + t462 * t479) * qJD(1) + (t465 ^ 2 * t410 + (t480 * t462 + (-t411 + t481) * t465) * t462) * qJD(2)) * t501 / 0.2e1 + m(3) * (t366 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(4) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(5) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(6) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(7) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + t436 * (t473 * t462 + t468 * t465) / 0.2e1 + t437 * (t468 * t462 - t473 * t465) / 0.2e1 + t428 * (t474 * t462 + t469 * t465) / 0.2e1 + t429 * (t469 * t462 - t474 * t465) / 0.2e1 + ((t420 * t535 + t421 * t534 + t514 * t536) * t438 + (t541 * t420 + t539 * t421 + t514 * t543) * t394 + (t540 * t420 + t538 * t421 + t542 * t514) * t393) * t393 / 0.2e1 + ((t418 * t535 + t419 * t534 + t515 * t536) * t438 + (t541 * t418 + t539 * t419 + t543 * t515) * t394 + (t418 * t540 + t419 * t538 + t515 * t542) * t393) * t394 / 0.2e1 + ((-t542 * t393 - t394 * t543 - t536 * t438) * t447 + ((-t460 * t535 + t463 * t534) * t438 + (-t460 * t541 + t463 * t539) * t394 + (-t460 * t540 + t463 * t538) * t393) * t446) * t438 / 0.2e1 + (Icges(2,3) + m(2) * (t443 ^ 2 + t444 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t413 * t464 + t415 * t461) * t462 - (t412 * t464 + t414 * t461) * t465) * qJD(2) + (t402 * t454 + t404 * t453) * t436 + (t401 * t454 + t403 * t453) * t437 + (t390 * t447 + t392 * t446) * t428 + (t389 * t447 + t391 * t446) * t429 + (t447 * t424 + t446 * t425 + t454 * t431 + t453 * t432 + t464 * t440 + t461 * t441) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
