% Calculate kinetic energy for
% S6RRRRRP2
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:39
% EndTime: 2019-03-10 01:00:41
% DurationCPUTime: 2.33s
% Computational Cost: add. (1933->259), mult. (1961->409), div. (0->0), fcn. (1872->10), ass. (0->154)
t542 = Icges(6,1) + Icges(7,1);
t541 = -Icges(6,4) + Icges(7,5);
t540 = Icges(7,4) + Icges(6,5);
t539 = Icges(6,2) + Icges(7,3);
t538 = -Icges(7,6) + Icges(6,6);
t537 = -Icges(6,3) - Icges(7,2);
t536 = rSges(7,1) + pkin(5);
t535 = rSges(7,3) + qJ(6);
t457 = qJ(2) + qJ(3);
t454 = qJ(4) + t457;
t447 = cos(t454);
t461 = cos(qJ(5));
t463 = cos(qJ(1));
t506 = t461 * t463;
t458 = sin(qJ(5));
t460 = sin(qJ(1));
t509 = t458 * t460;
t417 = t447 * t509 + t506;
t507 = t460 * t461;
t508 = t458 * t463;
t418 = t447 * t507 - t508;
t446 = sin(t454);
t511 = t446 * t460;
t534 = t539 * t417 + t541 * t418 - t538 * t511;
t419 = t447 * t508 - t507;
t420 = t447 * t506 + t509;
t510 = t446 * t463;
t533 = t539 * t419 + t541 * t420 - t538 * t510;
t532 = -t538 * t417 + t540 * t418 - t537 * t511;
t531 = -t538 * t419 + t540 * t420 - t537 * t510;
t530 = t541 * t417 + t542 * t418 + t540 * t511;
t529 = t541 * t419 + t542 * t420 + t540 * t510;
t528 = t538 * t447 + (t539 * t458 + t541 * t461) * t446;
t527 = t537 * t447 + (-t538 * t458 + t540 * t461) * t446;
t526 = -t540 * t447 + (t541 * t458 + t542 * t461) * t446;
t452 = sin(t457);
t521 = pkin(3) * t452;
t462 = cos(qJ(2));
t519 = t462 * pkin(2);
t459 = sin(qJ(2));
t517 = Icges(3,4) * t459;
t516 = Icges(3,4) * t462;
t515 = Icges(4,4) * t452;
t453 = cos(t457);
t514 = Icges(4,4) * t453;
t513 = Icges(5,4) * t446;
t512 = Icges(5,4) * t447;
t505 = rSges(7,2) * t511 + t535 * t417 + t418 * t536;
t504 = rSges(7,2) * t510 + t535 * t419 + t420 * t536;
t395 = -pkin(8) * t463 + t460 * t519;
t396 = pkin(8) * t460 + t463 * t519;
t451 = qJD(2) * t460;
t498 = qJD(2) * t463;
t503 = t395 * t451 + t396 * t498;
t502 = -rSges(7,2) * t447 + (t535 * t458 + t461 * t536) * t446;
t444 = pkin(1) * t460 - pkin(7) * t463;
t501 = -t395 - t444;
t500 = pkin(3) * t453;
t435 = qJD(3) * t460 + t451;
t497 = qJD(5) * t446;
t496 = -qJD(2) - qJD(3);
t495 = pkin(2) * qJD(2) * t459;
t370 = -pkin(9) * t463 + t460 * t500;
t494 = -t370 + t501;
t427 = qJD(4) * t460 + t435;
t493 = t463 * t495;
t492 = pkin(4) * t447 + pkin(10) * t446;
t491 = rSges(3,1) * t462 - rSges(3,2) * t459;
t490 = rSges(4,1) * t453 - rSges(4,2) * t452;
t489 = rSges(5,1) * t447 - rSges(5,2) * t446;
t428 = (-qJD(4) + t496) * t463;
t488 = Icges(3,1) * t462 - t517;
t487 = Icges(4,1) * t453 - t515;
t486 = Icges(5,1) * t447 - t513;
t485 = -Icges(3,2) * t459 + t516;
t484 = -Icges(4,2) * t452 + t514;
t483 = -Icges(5,2) * t446 + t512;
t482 = Icges(3,5) * t462 - Icges(3,6) * t459;
t481 = Icges(4,5) * t453 - Icges(4,6) * t452;
t480 = Icges(5,5) * t447 - Icges(5,6) * t446;
t411 = -Icges(3,6) * t463 + t460 * t485;
t413 = -Icges(3,5) * t463 + t460 * t488;
t479 = t411 * t459 - t413 * t462;
t412 = Icges(3,6) * t460 + t463 * t485;
t414 = Icges(3,5) * t460 + t463 * t488;
t478 = -t412 * t459 + t414 * t462;
t439 = Icges(3,2) * t462 + t517;
t440 = Icges(3,1) * t459 + t516;
t477 = -t439 * t459 + t440 * t462;
t436 = t496 * t463;
t476 = t436 * t521 - t493;
t371 = pkin(9) * t460 + t463 * t500;
t475 = t435 * t370 - t371 * t436 + t503;
t434 = qJD(1) * (pkin(1) * t463 + pkin(7) * t460);
t474 = qJD(1) * t396 - t460 * t495 + t434;
t473 = (Icges(5,5) * t446 + Icges(5,6) * t447) * qJD(1) + (-Icges(5,3) * t463 + t460 * t480) * t428 + (Icges(5,3) * t460 + t463 * t480) * t427;
t472 = (Icges(4,5) * t452 + Icges(4,6) * t453) * qJD(1) + (-Icges(4,3) * t463 + t460 * t481) * t436 + (Icges(4,3) * t460 + t463 * t481) * t435;
t407 = t492 * t460;
t408 = t492 * t463;
t471 = t427 * t407 - t408 * t428 + t475;
t470 = qJD(1) * t371 - t435 * t521 + t474;
t426 = pkin(4) * t446 - pkin(10) * t447;
t469 = t428 * t426 + (-t407 + t494) * qJD(1) + t476;
t468 = qJD(1) * t408 - t426 * t427 + t470;
t387 = -Icges(5,6) * t463 + t460 * t483;
t388 = Icges(5,6) * t460 + t463 * t483;
t389 = -Icges(5,5) * t463 + t460 * t486;
t390 = Icges(5,5) * t460 + t463 * t486;
t423 = Icges(5,2) * t447 + t513;
t424 = Icges(5,1) * t446 + t512;
t467 = (-t388 * t446 + t390 * t447) * t427 + (-t387 * t446 + t389 * t447) * t428 + (-t423 * t446 + t424 * t447) * qJD(1);
t399 = -Icges(4,6) * t463 + t460 * t484;
t400 = Icges(4,6) * t460 + t463 * t484;
t401 = -Icges(4,5) * t463 + t460 * t487;
t402 = Icges(4,5) * t460 + t463 * t487;
t430 = Icges(4,2) * t453 + t515;
t431 = Icges(4,1) * t452 + t514;
t466 = (-t400 * t452 + t402 * t453) * t435 + (-t399 * t452 + t401 * t453) * t436 + (-t430 * t452 + t431 * t453) * qJD(1);
t443 = rSges(2,1) * t463 - rSges(2,2) * t460;
t442 = rSges(2,1) * t460 + rSges(2,2) * t463;
t441 = rSges(3,1) * t459 + rSges(3,2) * t462;
t438 = Icges(3,5) * t459 + Icges(3,6) * t462;
t437 = -qJD(5) * t447 + qJD(1);
t432 = rSges(4,1) * t452 + rSges(4,2) * t453;
t425 = rSges(5,1) * t446 + rSges(5,2) * t447;
t416 = rSges(3,3) * t460 + t463 * t491;
t415 = -rSges(3,3) * t463 + t460 * t491;
t410 = Icges(3,3) * t460 + t463 * t482;
t409 = -Icges(3,3) * t463 + t460 * t482;
t404 = rSges(4,3) * t460 + t463 * t490;
t403 = -rSges(4,3) * t463 + t460 * t490;
t394 = rSges(5,3) * t460 + t463 * t489;
t393 = -rSges(5,3) * t463 + t460 * t489;
t392 = t460 * t497 + t428;
t391 = t463 * t497 + t427;
t383 = -rSges(6,3) * t447 + (rSges(6,1) * t461 - rSges(6,2) * t458) * t446;
t366 = qJD(1) * t416 - t441 * t451 + t434;
t365 = -t441 * t498 + (-t415 - t444) * qJD(1);
t363 = (t415 * t460 + t416 * t463) * qJD(2);
t362 = rSges(6,1) * t420 - rSges(6,2) * t419 + rSges(6,3) * t510;
t360 = rSges(6,1) * t418 - rSges(6,2) * t417 + rSges(6,3) * t511;
t346 = qJD(1) * t404 - t432 * t435 + t474;
t345 = -t493 + t432 * t436 + (-t403 + t501) * qJD(1);
t344 = t403 * t435 - t404 * t436 + t503;
t343 = qJD(1) * t394 - t425 * t427 + t470;
t342 = t425 * t428 + (-t393 + t494) * qJD(1) + t476;
t341 = t393 * t427 - t394 * t428 + t475;
t340 = t362 * t437 - t383 * t391 + t468;
t339 = -t360 * t437 + t383 * t392 + t469;
t338 = t360 * t391 - t362 * t392 + t471;
t337 = qJD(6) * t417 - t391 * t502 + t437 * t504 + t468;
t336 = qJD(6) * t419 + t392 * t502 - t437 * t505 + t469;
t335 = qJD(6) * t446 * t458 + t391 * t505 - t392 * t504 + t471;
t1 = t427 * (t473 * t460 + t467 * t463) / 0.2e1 + t428 * (t467 * t460 - t473 * t463) / 0.2e1 + t435 * (t472 * t460 + t466 * t463) / 0.2e1 + t436 * (t466 * t460 - t472 * t463) / 0.2e1 + m(7) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(5) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(4) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t363 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 - ((-t463 * t438 + t460 * t477) * qJD(1) + (t463 ^ 2 * t409 + (t478 * t460 + (-t410 + t479) * t463) * t460) * qJD(2)) * t498 / 0.2e1 + ((t460 * t438 + t463 * t477) * qJD(1) + (t460 ^ 2 * t410 + (t479 * t463 + (-t409 + t478) * t460) * t463) * qJD(2)) * t451 / 0.2e1 + ((t528 * t419 + t526 * t420 + t527 * t510) * t437 + (t534 * t419 + t530 * t420 + t532 * t510) * t392 + (t533 * t419 + t529 * t420 + t531 * t510) * t391) * t391 / 0.2e1 + ((t528 * t417 + t526 * t418 + t527 * t511) * t437 + (t534 * t417 + t530 * t418 + t532 * t511) * t392 + (t533 * t417 + t529 * t418 + t531 * t511) * t391) * t392 / 0.2e1 + ((-t531 * t391 - t532 * t392 - t527 * t437) * t447 + ((t528 * t458 + t526 * t461) * t437 + (t534 * t458 + t530 * t461) * t392 + (t533 * t458 + t529 * t461) * t391) * t446) * t437 / 0.2e1 + (m(2) * (t442 ^ 2 + t443 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t388 * t447 + t390 * t446) * t427 + (t387 * t447 + t389 * t446) * t428 + (t400 * t453 + t402 * t452) * t435 + (t399 * t453 + t401 * t452) * t436 + ((t412 * t462 + t414 * t459) * t460 - (t411 * t462 + t413 * t459) * t463) * qJD(2) + (t447 * t423 + t446 * t424 + t453 * t430 + t452 * t431 + t462 * t439 + t459 * t440) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
