% Calculate kinetic energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:54
% EndTime: 2019-03-09 15:24:57
% DurationCPUTime: 3.08s
% Computational Cost: add. (1633->263), mult. (1659->393), div. (0->0), fcn. (1504->10), ass. (0->154)
t548 = Icges(5,4) + Icges(6,6);
t547 = Icges(5,1) + Icges(6,2);
t546 = -Icges(5,2) - Icges(6,3);
t453 = qJ(2) + qJ(3);
t445 = pkin(10) + t453;
t443 = cos(t445);
t545 = t548 * t443;
t442 = sin(t445);
t544 = t548 * t442;
t543 = Icges(6,4) - Icges(5,5);
t542 = Icges(6,5) - Icges(5,6);
t541 = t546 * t442 + t545;
t540 = t547 * t443 - t544;
t456 = sin(qJ(1));
t459 = cos(qJ(1));
t539 = t541 * t456 + t542 * t459;
t538 = -t542 * t456 + t541 * t459;
t537 = t540 * t456 + t543 * t459;
t536 = -t543 * t456 + t540 * t459;
t535 = t546 * t443 - t544;
t534 = t547 * t442 + t545;
t533 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t449 = sin(t453);
t450 = cos(t453);
t532 = Icges(4,5) * t450 - Icges(4,6) * t449 + t542 * t442 - t543 * t443;
t518 = Icges(4,4) * t450;
t484 = -Icges(4,2) * t449 + t518;
t390 = -Icges(4,6) * t459 + t456 * t484;
t391 = Icges(4,6) * t456 + t459 * t484;
t519 = Icges(4,4) * t449;
t487 = Icges(4,1) * t450 - t519;
t392 = -Icges(4,5) * t459 + t456 * t487;
t393 = Icges(4,5) * t456 + t459 * t487;
t425 = Icges(4,2) * t450 + t519;
t426 = Icges(4,1) * t449 + t518;
t448 = qJD(2) * t456;
t433 = qJD(3) * t456 + t448;
t434 = (-qJD(2) - qJD(3)) * t459;
t531 = (-t390 * t449 + t392 * t450 - t442 * t539 + t537 * t443) * t434 + (-t391 * t449 + t393 * t450 - t442 * t538 + t443 * t536) * t433 + (-t425 * t449 + t426 * t450 + t442 * t535 + t443 * t534) * qJD(1);
t530 = (t456 * t532 - t459 * t533) * t434 + (t456 * t533 + t459 * t532) * t433 + (Icges(4,5) * t449 + Icges(4,6) * t450 - t543 * t442 - t542 * t443) * qJD(1);
t526 = pkin(3) * t449;
t525 = pkin(9) * t442;
t458 = cos(qJ(2));
t523 = t458 * pkin(2);
t455 = sin(qJ(2));
t521 = Icges(3,4) * t455;
t520 = Icges(3,4) * t458;
t513 = t443 * t456;
t512 = t443 * t459;
t454 = sin(qJ(6));
t511 = t454 * t456;
t510 = t454 * t459;
t457 = cos(qJ(6));
t509 = t456 * t457;
t508 = t457 * t459;
t504 = pkin(3) * t450;
t359 = qJ(4) * t456 + t459 * t504;
t489 = pkin(4) * t443 + qJ(5) * t442;
t397 = t489 * t459;
t507 = -t359 - t397;
t385 = -pkin(8) * t459 + t456 * t523;
t386 = pkin(8) * t456 + t459 * t523;
t502 = qJD(2) * t459;
t506 = t385 * t448 + t386 * t502;
t441 = pkin(1) * t456 - pkin(7) * t459;
t505 = -t385 - t441;
t501 = qJD(5) * t442;
t500 = qJD(6) * t443;
t499 = pkin(2) * qJD(2) * t455;
t358 = -qJ(4) * t459 + t456 * t504;
t498 = t433 * t358 + t506;
t497 = -t358 + t505;
t418 = pkin(4) * t442 - qJ(5) * t443;
t496 = -t418 - t526;
t495 = t459 * t499;
t396 = t489 * t456;
t494 = -t396 + t497;
t493 = rSges(3,1) * t458 - rSges(3,2) * t455;
t492 = rSges(4,1) * t450 - rSges(4,2) * t449;
t491 = rSges(5,1) * t443 - rSges(5,2) * t442;
t490 = -rSges(6,2) * t443 + rSges(6,3) * t442;
t488 = Icges(3,1) * t458 - t521;
t485 = -Icges(3,2) * t455 + t520;
t481 = Icges(3,5) * t458 - Icges(3,6) * t455;
t402 = -Icges(3,6) * t459 + t456 * t485;
t404 = -Icges(3,5) * t459 + t456 * t488;
t476 = t402 * t455 - t404 * t458;
t403 = Icges(3,6) * t456 + t459 * t485;
t405 = Icges(3,5) * t456 + t459 * t488;
t475 = -t403 * t455 + t405 * t458;
t436 = Icges(3,2) * t458 + t521;
t437 = Icges(3,1) * t455 + t520;
t474 = -t436 * t455 + t437 * t458;
t429 = qJD(1) * (pkin(1) * t459 + pkin(7) * t456);
t473 = qJD(1) * t386 - t456 * t499 + t429;
t472 = qJD(4) * t456 + t434 * t526 - t495;
t471 = -qJD(5) * t443 + t433 * t396 + t498;
t467 = t434 * t418 + t459 * t501 + t472;
t466 = qJD(1) * t359 - qJD(4) * t459 + t473;
t465 = qJD(1) * t397 + t456 * t501 + t466;
t440 = rSges(2,1) * t459 - rSges(2,2) * t456;
t439 = rSges(2,1) * t456 + rSges(2,2) * t459;
t438 = rSges(3,1) * t455 + rSges(3,2) * t458;
t435 = Icges(3,5) * t455 + Icges(3,6) * t458;
t432 = qJD(6) * t442 + qJD(1);
t427 = rSges(4,1) * t449 + rSges(4,2) * t450;
t423 = -pkin(5) * t459 + pkin(9) * t513;
t422 = pkin(5) * t456 + pkin(9) * t512;
t420 = rSges(5,1) * t442 + rSges(5,2) * t443;
t419 = -rSges(6,2) * t442 - rSges(6,3) * t443;
t411 = t442 * t511 - t508;
t410 = t442 * t509 + t510;
t409 = t442 * t510 + t509;
t408 = t442 * t508 - t511;
t407 = rSges(3,3) * t456 + t459 * t493;
t406 = -rSges(3,3) * t459 + t456 * t493;
t401 = Icges(3,3) * t456 + t459 * t481;
t400 = -Icges(3,3) * t459 + t456 * t481;
t399 = t456 * t500 + t434;
t398 = t459 * t500 + t433;
t395 = rSges(4,3) * t456 + t459 * t492;
t394 = -rSges(4,3) * t459 + t456 * t492;
t383 = -rSges(6,1) * t459 + t456 * t490;
t382 = rSges(6,1) * t456 + t459 * t490;
t381 = rSges(5,3) * t456 + t459 * t491;
t380 = -rSges(5,3) * t459 + t456 * t491;
t364 = rSges(7,3) * t442 + (-rSges(7,1) * t454 - rSges(7,2) * t457) * t443;
t363 = Icges(7,5) * t442 + (-Icges(7,1) * t454 - Icges(7,4) * t457) * t443;
t362 = Icges(7,6) * t442 + (-Icges(7,4) * t454 - Icges(7,2) * t457) * t443;
t361 = Icges(7,3) * t442 + (-Icges(7,5) * t454 - Icges(7,6) * t457) * t443;
t356 = qJD(1) * t407 - t438 * t448 + t429;
t355 = -t438 * t502 + (-t406 - t441) * qJD(1);
t354 = (t406 * t456 + t407 * t459) * qJD(2);
t352 = rSges(7,1) * t411 + rSges(7,2) * t410 + rSges(7,3) * t513;
t351 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t512;
t350 = Icges(7,1) * t411 + Icges(7,4) * t410 + Icges(7,5) * t513;
t349 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t512;
t348 = Icges(7,4) * t411 + Icges(7,2) * t410 + Icges(7,6) * t513;
t347 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t512;
t346 = Icges(7,5) * t411 + Icges(7,6) * t410 + Icges(7,3) * t513;
t345 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t512;
t344 = qJD(1) * t395 - t427 * t433 + t473;
t343 = -t495 + t427 * t434 + (-t394 + t505) * qJD(1);
t342 = t394 * t433 - t395 * t434 + t506;
t341 = qJD(1) * t381 + (-t420 - t526) * t433 + t466;
t340 = t420 * t434 + (-t380 + t497) * qJD(1) + t472;
t339 = qJD(1) * t382 + (-t419 + t496) * t433 + t465;
t338 = t419 * t434 + (-t383 + t494) * qJD(1) + t467;
t337 = t380 * t433 + (-t359 - t381) * t434 + t498;
t336 = t383 * t433 + (-t382 + t507) * t434 + t471;
t335 = qJD(1) * t422 + t351 * t432 - t364 * t398 + (t496 - t525) * t433 + t465;
t334 = t434 * t525 - t352 * t432 + t364 * t399 + (-t423 + t494) * qJD(1) + t467;
t333 = -t351 * t399 + t352 * t398 + t423 * t433 + (-t422 + t507) * t434 + t471;
t1 = m(4) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(3) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(5) * (t337 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + t398 * ((t345 * t512 + t408 * t347 + t409 * t349) * t398 + (t346 * t512 + t348 * t408 + t350 * t409) * t399 + (t361 * t512 + t362 * t408 + t363 * t409) * t432) / 0.2e1 + t399 * ((t345 * t513 + t347 * t410 + t349 * t411) * t398 + (t346 * t513 + t410 * t348 + t411 * t350) * t399 + (t361 * t513 + t362 * t410 + t363 * t411) * t432) / 0.2e1 + t432 * ((t345 * t398 + t346 * t399 + t361 * t432) * t442 + ((-t347 * t457 - t349 * t454) * t398 + (-t348 * t457 - t350 * t454) * t399 + (-t362 * t457 - t363 * t454) * t432) * t443) / 0.2e1 - ((-t459 * t435 + t456 * t474) * qJD(1) + (t459 ^ 2 * t400 + (t475 * t456 + (-t401 + t476) * t459) * t456) * qJD(2)) * t502 / 0.2e1 + ((t456 * t435 + t459 * t474) * qJD(1) + (t456 ^ 2 * t401 + (t476 * t459 + (-t400 + t475) * t456) * t459) * qJD(2)) * t448 / 0.2e1 + (m(2) * (t439 ^ 2 + t440 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (t530 * t456 + t531 * t459) * t433 / 0.2e1 + (t531 * t456 - t530 * t459) * t434 / 0.2e1 + (((t403 * t458 + t405 * t455) * t456 - (t402 * t458 + t404 * t455) * t459) * qJD(2) + (t390 * t450 + t392 * t449 + t537 * t442 + t443 * t539) * t434 + (t391 * t450 + t393 * t449 + t442 * t536 + t443 * t538) * t433 + (t450 * t425 + t449 * t426 + t458 * t436 + t455 * t437 + t534 * t442 - t535 * t443) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
