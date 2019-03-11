% Calculate kinetic energy for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:44:00
% EndTime: 2019-03-09 20:44:03
% DurationCPUTime: 2.47s
% Computational Cost: add. (1982->277), mult. (2205->417), div. (0->0), fcn. (2172->10), ass. (0->151)
t545 = Icges(6,1) + Icges(7,1);
t544 = -Icges(6,4) + Icges(7,5);
t543 = Icges(7,4) + Icges(6,5);
t542 = Icges(6,2) + Icges(7,3);
t541 = -Icges(7,6) + Icges(6,6);
t540 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t539 = rSges(7,1) + pkin(5);
t538 = rSges(7,3) + qJ(6);
t464 = qJ(4) + pkin(10);
t459 = cos(t464);
t465 = qJ(2) + qJ(3);
t463 = cos(t465);
t472 = cos(qJ(1));
t458 = sin(t464);
t469 = sin(qJ(1));
t515 = t458 * t469;
t416 = t459 * t472 + t463 * t515;
t509 = t469 * t459;
t417 = -t458 * t472 + t463 * t509;
t462 = sin(t465);
t514 = t462 * t469;
t537 = t542 * t416 + t544 * t417 - t541 * t514;
t512 = t463 * t472;
t418 = t458 * t512 - t509;
t419 = t459 * t512 + t515;
t513 = t462 * t472;
t536 = t542 * t418 + t544 * t419 - t541 * t513;
t535 = t544 * t416 + t545 * t417 + t543 * t514;
t534 = t544 * t418 + t545 * t419 + t543 * t513;
t533 = t541 * t463 + (t542 * t458 + t544 * t459) * t462;
t532 = -t543 * t463 + (t544 * t458 + t545 * t459) * t462;
t470 = cos(qJ(4));
t507 = t470 * t472;
t467 = sin(qJ(4));
t511 = t467 * t469;
t433 = -t463 * t511 - t507;
t508 = t469 * t470;
t510 = t467 * t472;
t434 = t463 * t508 - t510;
t531 = Icges(5,5) * t434 + Icges(5,6) * t433 - t541 * t416 + t543 * t417 - t540 * t514;
t435 = -t463 * t510 + t508;
t436 = t463 * t507 + t511;
t530 = Icges(5,5) * t436 + Icges(5,6) * t435 - t541 * t418 + t543 * t419 - t540 * t513;
t529 = t540 * t463 + (Icges(5,5) * t470 - Icges(5,6) * t467 - t541 * t458 + t543 * t459) * t462;
t471 = cos(qJ(2));
t523 = pkin(2) * t471;
t522 = pkin(4) * t470;
t468 = sin(qJ(2));
t519 = Icges(3,4) * t468;
t518 = Icges(3,4) * t471;
t517 = Icges(4,4) * t462;
t516 = Icges(4,4) * t463;
t506 = rSges(7,2) * t514 + t538 * t416 + t539 * t417;
t505 = rSges(7,2) * t513 + t538 * t418 + t539 * t419;
t504 = -rSges(7,2) * t463 + (t538 * t458 + t539 * t459) * t462;
t404 = -pkin(8) * t472 + t469 * t523;
t405 = pkin(8) * t469 + t472 * t523;
t461 = qJD(2) * t469;
t501 = qJD(2) * t472;
t503 = t404 * t461 + t405 * t501;
t453 = pkin(1) * t469 - pkin(7) * t472;
t502 = -t404 - t453;
t443 = qJD(3) * t469 + t461;
t500 = qJD(4) * t462;
t499 = qJD(5) * t462;
t498 = pkin(2) * qJD(2) * t468;
t444 = (-qJD(2) - qJD(3)) * t472;
t497 = t472 * t498;
t496 = pkin(3) * t463 + pkin(9) * t462;
t495 = rSges(3,1) * t471 - rSges(3,2) * t468;
t494 = rSges(4,1) * t463 - rSges(4,2) * t462;
t493 = Icges(3,1) * t471 - t519;
t492 = Icges(4,1) * t463 - t517;
t491 = -Icges(3,2) * t468 + t518;
t490 = -Icges(4,2) * t462 + t516;
t489 = Icges(3,5) * t471 - Icges(3,6) * t468;
t488 = Icges(4,5) * t463 - Icges(4,6) * t462;
t423 = -Icges(3,6) * t472 + t469 * t491;
t425 = -Icges(3,5) * t472 + t469 * t493;
t487 = t423 * t468 - t425 * t471;
t424 = Icges(3,6) * t469 + t472 * t491;
t426 = Icges(3,5) * t469 + t472 * t493;
t486 = -t424 * t468 + t426 * t471;
t446 = Icges(3,2) * t471 + t519;
t447 = Icges(3,1) * t468 + t518;
t485 = -t446 * t468 + t447 * t471;
t431 = t496 * t469;
t432 = t496 * t472;
t484 = t443 * t431 - t432 * t444 + t503;
t442 = qJD(1) * (pkin(1) * t472 + pkin(7) * t469);
t483 = qJD(1) * t405 - t469 * t498 + t442;
t482 = qJ(5) * t462 + t463 * t522;
t481 = (Icges(4,5) * t462 + Icges(4,6) * t463) * qJD(1) + (-Icges(4,3) * t472 + t469 * t488) * t444 + (Icges(4,3) * t469 + t472 * t488) * t443;
t372 = -pkin(4) * t510 + t469 * t482;
t427 = t472 * t500 + t443;
t480 = -qJD(5) * t463 + t427 * t372 + t484;
t441 = pkin(3) * t462 - pkin(9) * t463;
t479 = qJD(1) * t432 - t441 * t443 + t483;
t478 = t444 * t441 + (-t431 + t502) * qJD(1) - t497;
t373 = pkin(4) * t511 + t472 * t482;
t454 = -qJD(4) * t463 + qJD(1);
t477 = t454 * t373 + t469 * t499 + t479;
t387 = -qJ(5) * t463 + t462 * t522;
t428 = t469 * t500 + t444;
t476 = t428 * t387 + t472 * t499 + t478;
t409 = -Icges(4,6) * t472 + t469 * t490;
t410 = Icges(4,6) * t469 + t472 * t490;
t411 = -Icges(4,5) * t472 + t469 * t492;
t412 = Icges(4,5) * t469 + t472 * t492;
t438 = Icges(4,2) * t463 + t517;
t439 = Icges(4,1) * t462 + t516;
t475 = (-t410 * t462 + t412 * t463) * t443 + (-t409 * t462 + t411 * t463) * t444 + (-t438 * t462 + t439 * t463) * qJD(1);
t450 = rSges(2,1) * t472 - rSges(2,2) * t469;
t449 = rSges(2,1) * t469 + rSges(2,2) * t472;
t448 = rSges(3,1) * t468 + rSges(3,2) * t471;
t445 = Icges(3,5) * t468 + Icges(3,6) * t471;
t440 = rSges(4,1) * t462 + rSges(4,2) * t463;
t430 = rSges(3,3) * t469 + t472 * t495;
t429 = -rSges(3,3) * t472 + t469 * t495;
t422 = Icges(3,3) * t469 + t472 * t489;
t421 = -Icges(3,3) * t472 + t469 * t489;
t415 = rSges(4,3) * t469 + t472 * t494;
t414 = -rSges(4,3) * t472 + t469 * t494;
t403 = -rSges(5,3) * t463 + (rSges(5,1) * t470 - rSges(5,2) * t467) * t462;
t402 = -Icges(5,5) * t463 + (Icges(5,1) * t470 - Icges(5,4) * t467) * t462;
t401 = -Icges(5,6) * t463 + (Icges(5,4) * t470 - Icges(5,2) * t467) * t462;
t395 = -rSges(6,3) * t463 + (rSges(6,1) * t459 - rSges(6,2) * t458) * t462;
t386 = qJD(1) * t430 - t448 * t461 + t442;
t385 = -t448 * t501 + (-t429 - t453) * qJD(1);
t384 = rSges(5,1) * t436 + rSges(5,2) * t435 + rSges(5,3) * t513;
t383 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t514;
t380 = Icges(5,1) * t436 + Icges(5,4) * t435 + Icges(5,5) * t513;
t379 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t514;
t378 = Icges(5,4) * t436 + Icges(5,2) * t435 + Icges(5,6) * t513;
t377 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t514;
t374 = (t429 * t469 + t430 * t472) * qJD(2);
t370 = rSges(6,1) * t419 - rSges(6,2) * t418 + rSges(6,3) * t513;
t368 = rSges(6,1) * t417 - rSges(6,2) * t416 + rSges(6,3) * t514;
t352 = qJD(1) * t415 - t440 * t443 + t483;
t351 = -t497 + t440 * t444 + (-t414 + t502) * qJD(1);
t350 = t414 * t443 - t415 * t444 + t503;
t349 = t384 * t454 - t403 * t427 + t479;
t348 = -t383 * t454 + t403 * t428 + t478;
t347 = t383 * t427 - t384 * t428 + t484;
t346 = t370 * t454 + (-t387 - t395) * t427 + t477;
t345 = t395 * t428 + (-t368 - t372) * t454 + t476;
t344 = t368 * t427 + (-t370 - t373) * t428 + t480;
t343 = qJD(6) * t416 + t505 * t454 + (-t387 - t504) * t427 + t477;
t342 = qJD(6) * t418 + t504 * t428 + (-t372 - t506) * t454 + t476;
t341 = qJD(6) * t458 * t462 + t506 * t427 + (-t373 - t505) * t428 + t480;
t1 = ((t469 * t445 + t472 * t485) * qJD(1) + (t469 ^ 2 * t422 + (t487 * t472 + (-t421 + t486) * t469) * t472) * qJD(2)) * t461 / 0.2e1 - ((-t472 * t445 + t469 * t485) * qJD(1) + (t472 ^ 2 * t421 + (t486 * t469 + (-t422 + t487) * t472) * t469) * qJD(2)) * t501 / 0.2e1 + m(3) * (t374 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(4) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(5) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(6) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(7) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + t443 * (t481 * t469 + t475 * t472) / 0.2e1 + t444 * (t475 * t469 - t481 * t472) / 0.2e1 + (Icges(2,3) + m(2) * (t449 ^ 2 + t450 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t424 * t471 + t426 * t468) * t469 - (t423 * t471 + t425 * t468) * t472) * qJD(2) + (t410 * t463 + t412 * t462) * t443 + (t409 * t463 + t411 * t462) * t444 + (t438 * t463 + t462 * t439 + t471 * t446 + t468 * t447) * qJD(1)) * qJD(1) / 0.2e1 + ((t401 * t435 + t402 * t436 + t533 * t418 + t532 * t419 + t529 * t513) * t454 + (t377 * t435 + t379 * t436 + t537 * t418 + t535 * t419 + t531 * t513) * t428 + (t435 * t378 + t436 * t380 + t536 * t418 + t534 * t419 + t530 * t513) * t427) * t427 / 0.2e1 + ((t401 * t433 + t402 * t434 + t533 * t416 + t532 * t417 + t529 * t514) * t454 + (t433 * t377 + t434 * t379 + t537 * t416 + t535 * t417 + t531 * t514) * t428 + (t378 * t433 + t380 * t434 + t536 * t416 + t534 * t417 + t530 * t514) * t427) * t428 / 0.2e1 + ((-t530 * t427 - t531 * t428 - t529 * t454) * t463 + ((-t401 * t467 + t402 * t470 + t533 * t458 + t532 * t459) * t454 + (-t377 * t467 + t379 * t470 + t537 * t458 + t535 * t459) * t428 + (-t378 * t467 + t380 * t470 + t536 * t458 + t534 * t459) * t427) * t462) * t454 / 0.2e1;
T  = t1;
