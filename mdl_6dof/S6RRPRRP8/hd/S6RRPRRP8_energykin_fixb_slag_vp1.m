% Calculate kinetic energy for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:21:56
% EndTime: 2019-03-09 12:21:59
% DurationCPUTime: 2.98s
% Computational Cost: add. (1957->309), mult. (2549->478), div. (0->0), fcn. (2572->10), ass. (0->154)
t544 = Icges(6,1) + Icges(7,1);
t543 = -Icges(6,4) + Icges(7,5);
t542 = Icges(7,4) + Icges(6,5);
t541 = Icges(6,2) + Icges(7,3);
t540 = -Icges(7,6) + Icges(6,6);
t539 = -Icges(6,3) - Icges(7,2);
t538 = rSges(7,1) + pkin(5);
t537 = rSges(7,3) + qJ(6);
t470 = pkin(10) + qJ(4);
t466 = qJ(5) + t470;
t461 = sin(t466);
t462 = cos(t466);
t477 = cos(qJ(1));
t475 = sin(qJ(1));
t476 = cos(qJ(2));
t514 = t475 * t476;
t409 = t461 * t514 + t462 * t477;
t410 = -t461 * t477 + t462 * t514;
t474 = sin(qJ(2));
t516 = t474 * t475;
t536 = t541 * t409 + t543 * t410 - t540 * t516;
t513 = t476 * t477;
t411 = t461 * t513 - t475 * t462;
t412 = t461 * t475 + t462 * t513;
t515 = t474 * t477;
t535 = t541 * t411 + t543 * t412 - t540 * t515;
t534 = -t540 * t409 + t542 * t410 - t539 * t516;
t533 = -t540 * t411 + t542 * t412 - t539 * t515;
t532 = t543 * t409 + t544 * t410 + t542 * t516;
t531 = t543 * t411 + t544 * t412 + t542 * t515;
t530 = t540 * t476 + (t541 * t461 + t543 * t462) * t474;
t529 = t539 * t476 + (-t540 * t461 + t542 * t462) * t474;
t528 = -t542 * t476 + (t543 * t461 + t544 * t462) * t474;
t472 = cos(pkin(10));
t521 = t472 * pkin(3);
t520 = Icges(3,4) * t474;
t519 = Icges(3,4) * t476;
t471 = sin(pkin(10));
t518 = t471 * t475;
t517 = t471 * t477;
t511 = rSges(7,2) * t516 + t537 * t409 + t538 * t410;
t510 = rSges(7,2) * t515 + t537 * t411 + t538 * t412;
t509 = -rSges(7,2) * t476 + (t537 * t461 + t538 * t462) * t474;
t493 = pkin(2) * t476 + qJ(3) * t474;
t439 = t493 * t475;
t455 = pkin(1) * t475 - pkin(7) * t477;
t508 = -t439 - t455;
t465 = cos(t470);
t507 = pkin(4) * t465;
t467 = qJD(2) * t475;
t503 = qJD(4) * t474;
t441 = t477 * t503 + t467;
t505 = qJD(2) * t477;
t504 = qJD(3) * t474;
t502 = qJD(5) * t474;
t440 = t493 * t477;
t446 = qJD(1) * (pkin(1) * t477 + pkin(7) * t475);
t501 = qJD(1) * t440 + t475 * t504 + t446;
t464 = sin(t470);
t498 = pkin(4) * t464;
t450 = pkin(2) * t474 - qJ(3) * t476;
t497 = qJD(2) * (pkin(8) * t476 - t474 * t521 - t450);
t496 = qJD(2) * (rSges(4,3) * t476 - (rSges(4,1) * t472 - rSges(4,2) * t471) * t474 - t450);
t442 = t475 * t503 - t505;
t495 = -qJD(3) * t476 + t439 * t467 + t440 * t505;
t494 = rSges(3,1) * t476 - rSges(3,2) * t474;
t492 = Icges(3,1) * t476 - t520;
t491 = -Icges(3,2) * t474 + t519;
t490 = Icges(3,5) * t476 - Icges(3,6) * t474;
t426 = -Icges(3,6) * t477 + t475 * t491;
t428 = -Icges(3,5) * t477 + t475 * t492;
t489 = t426 * t474 - t428 * t476;
t427 = Icges(3,6) * t475 + t477 * t491;
t429 = Icges(3,5) * t475 + t477 * t492;
t488 = -t427 * t474 + t429 * t476;
t448 = Icges(3,2) * t476 + t520;
t449 = Icges(3,1) * t474 + t519;
t487 = -t448 * t474 + t449 * t476;
t485 = pkin(8) * t474 + t476 * t521;
t393 = -pkin(3) * t517 + t475 * t485;
t394 = pkin(3) * t518 + t477 * t485;
t486 = t393 * t467 + t394 * t505 + t495;
t484 = pkin(9) * t474 + t476 * t507;
t483 = qJD(1) * t394 + t475 * t497 + t501;
t350 = t475 * t484 - t477 * t498;
t351 = t475 * t498 + t477 * t484;
t482 = t441 * t350 - t351 * t442 + t486;
t459 = t477 * t504;
t481 = t459 + (-t393 + t508) * qJD(1) + t477 * t497;
t395 = -pkin(9) * t476 + t474 * t507;
t460 = -qJD(4) * t476 + qJD(1);
t480 = t460 * t351 - t395 * t441 + t483;
t479 = -t350 * t460 + t442 * t395 + t481;
t453 = rSges(2,1) * t477 - rSges(2,2) * t475;
t452 = rSges(2,1) * t475 + rSges(2,2) * t477;
t451 = rSges(3,1) * t474 + rSges(3,2) * t476;
t447 = Icges(3,5) * t474 + Icges(3,6) * t476;
t445 = qJD(1) + (-qJD(4) - qJD(5)) * t476;
t438 = t472 * t513 + t518;
t437 = -t471 * t513 + t472 * t475;
t436 = t472 * t514 - t517;
t435 = -t471 * t514 - t472 * t477;
t433 = rSges(3,3) * t475 + t477 * t494;
t432 = -rSges(3,3) * t477 + t475 * t494;
t425 = Icges(3,3) * t475 + t477 * t490;
t424 = -Icges(3,3) * t477 + t475 * t490;
t423 = t464 * t475 + t465 * t513;
t422 = -t464 * t513 + t465 * t475;
t421 = -t464 * t477 + t465 * t514;
t420 = -t464 * t514 - t465 * t477;
t418 = t475 * t502 + t442;
t417 = t477 * t502 + t441;
t416 = -Icges(4,5) * t476 + (Icges(4,1) * t472 - Icges(4,4) * t471) * t474;
t415 = -Icges(4,6) * t476 + (Icges(4,4) * t472 - Icges(4,2) * t471) * t474;
t414 = -Icges(4,3) * t476 + (Icges(4,5) * t472 - Icges(4,6) * t471) * t474;
t408 = -rSges(5,3) * t476 + (rSges(5,1) * t465 - rSges(5,2) * t464) * t474;
t407 = -Icges(5,5) * t476 + (Icges(5,1) * t465 - Icges(5,4) * t464) * t474;
t406 = -Icges(5,6) * t476 + (Icges(5,4) * t465 - Icges(5,2) * t464) * t474;
t405 = -Icges(5,3) * t476 + (Icges(5,5) * t465 - Icges(5,6) * t464) * t474;
t403 = -rSges(6,3) * t476 + (rSges(6,1) * t462 - rSges(6,2) * t461) * t474;
t392 = rSges(4,1) * t438 + rSges(4,2) * t437 + rSges(4,3) * t515;
t391 = rSges(4,1) * t436 + rSges(4,2) * t435 + rSges(4,3) * t516;
t390 = Icges(4,1) * t438 + Icges(4,4) * t437 + Icges(4,5) * t515;
t389 = Icges(4,1) * t436 + Icges(4,4) * t435 + Icges(4,5) * t516;
t388 = Icges(4,4) * t438 + Icges(4,2) * t437 + Icges(4,6) * t515;
t387 = Icges(4,4) * t436 + Icges(4,2) * t435 + Icges(4,6) * t516;
t386 = Icges(4,5) * t438 + Icges(4,6) * t437 + Icges(4,3) * t515;
t385 = Icges(4,5) * t436 + Icges(4,6) * t435 + Icges(4,3) * t516;
t383 = qJD(1) * t433 - t451 * t467 + t446;
t382 = -t451 * t505 + (-t432 - t455) * qJD(1);
t379 = (t432 * t475 + t433 * t477) * qJD(2);
t376 = rSges(5,1) * t423 + rSges(5,2) * t422 + rSges(5,3) * t515;
t375 = rSges(5,1) * t421 + rSges(5,2) * t420 + rSges(5,3) * t516;
t374 = Icges(5,1) * t423 + Icges(5,4) * t422 + Icges(5,5) * t515;
t373 = Icges(5,1) * t421 + Icges(5,4) * t420 + Icges(5,5) * t516;
t372 = Icges(5,4) * t423 + Icges(5,2) * t422 + Icges(5,6) * t515;
t371 = Icges(5,4) * t421 + Icges(5,2) * t420 + Icges(5,6) * t516;
t370 = Icges(5,5) * t423 + Icges(5,6) * t422 + Icges(5,3) * t515;
t369 = Icges(5,5) * t421 + Icges(5,6) * t420 + Icges(5,3) * t516;
t367 = rSges(6,1) * t412 - rSges(6,2) * t411 + rSges(6,3) * t515;
t365 = rSges(6,1) * t410 - rSges(6,2) * t409 + rSges(6,3) * t516;
t347 = qJD(1) * t392 + t475 * t496 + t501;
t346 = t459 + t477 * t496 + (-t391 + t508) * qJD(1);
t345 = (t391 * t475 + t392 * t477) * qJD(2) + t495;
t344 = t376 * t460 - t408 * t441 + t483;
t343 = -t375 * t460 + t408 * t442 + t481;
t342 = t375 * t441 - t376 * t442 + t486;
t341 = t367 * t445 - t403 * t417 + t480;
t340 = -t365 * t445 + t403 * t418 + t479;
t339 = t365 * t417 - t367 * t418 + t482;
t338 = qJD(6) * t409 - t417 * t509 + t445 * t510 + t480;
t337 = qJD(6) * t411 + t418 * t509 - t445 * t511 + t479;
t336 = qJD(6) * t461 * t474 + t417 * t511 - t418 * t510 + t482;
t1 = m(3) * (t379 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(4) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(5) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(6) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(7) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + t441 * ((t370 * t515 + t422 * t372 + t423 * t374) * t441 + (t369 * t515 + t371 * t422 + t373 * t423) * t442 + (t405 * t515 + t406 * t422 + t407 * t423) * t460) / 0.2e1 + t442 * ((t370 * t516 + t372 * t420 + t374 * t421) * t441 + (t369 * t516 + t420 * t371 + t421 * t373) * t442 + (t405 * t516 + t406 * t420 + t407 * t421) * t460) / 0.2e1 + t460 * ((-t369 * t442 - t370 * t441 - t405 * t460) * t476 + ((-t372 * t464 + t374 * t465) * t441 + (-t371 * t464 + t373 * t465) * t442 + (-t406 * t464 + t407 * t465) * t460) * t474) / 0.2e1 + ((t411 * t530 + t412 * t528 + t515 * t529) * t445 + (t411 * t536 + t532 * t412 + t534 * t515) * t418 + (t535 * t411 + t531 * t412 + t533 * t515) * t417) * t417 / 0.2e1 + ((t409 * t530 + t410 * t528 + t516 * t529) * t445 + (t536 * t409 + t532 * t410 + t534 * t516) * t418 + (t409 * t535 + t410 * t531 + t516 * t533) * t417) * t418 / 0.2e1 + ((-t417 * t533 - t418 * t534 - t445 * t529) * t476 + ((t461 * t530 + t462 * t528) * t445 + (t461 * t536 + t532 * t462) * t418 + (t461 * t535 + t462 * t531) * t417) * t474) * t445 / 0.2e1 + (m(2) * (t452 ^ 2 + t453 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t427 * t476 + t429 * t474) * t475 - (t426 * t476 + t428 * t474) * t477 + (t385 * t477 - t386 * t475) * t476 + ((-t388 * t471 + t390 * t472) * t475 - (-t387 * t471 + t389 * t472) * t477) * t474) * qJD(2) + ((t448 - t414) * t476 + (-t415 * t471 + t416 * t472 + t449) * t474) * qJD(1)) * qJD(1) / 0.2e1 + (((-t385 * t515 - t387 * t437 - t389 * t438 + t489 * t477) * t477 + ((-t424 + t488) * t477 + t386 * t515 + t388 * t437 + t390 * t438 + t425 * t475) * t475) * qJD(2) + (t414 * t515 + t415 * t437 + t416 * t438 + t475 * t447 + t477 * t487) * qJD(1)) * t467 / 0.2e1 - (((-t385 * t516 - t387 * t435 - t389 * t436 + t424 * t477) * t477 + ((-t425 + t489) * t477 + t386 * t516 + t388 * t435 + t390 * t436 + t488 * t475) * t475) * qJD(2) + (t414 * t516 + t415 * t435 + t416 * t436 - t477 * t447 + t475 * t487) * qJD(1)) * t505 / 0.2e1;
T  = t1;
