% Calculate kinetic energy for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:25
% EndTime: 2019-03-09 13:12:27
% DurationCPUTime: 2.16s
% Computational Cost: add. (1953->283), mult. (1697->447), div. (0->0), fcn. (1542->12), ass. (0->166)
t539 = Icges(3,3) + Icges(4,3);
t452 = qJ(2) + pkin(11);
t442 = sin(t452);
t443 = cos(t452);
t455 = sin(qJ(2));
t458 = cos(qJ(2));
t538 = Icges(3,5) * t458 + Icges(4,5) * t443 - Icges(3,6) * t455 - Icges(4,6) * t442;
t456 = sin(qJ(1));
t459 = cos(qJ(1));
t537 = t538 * t456 - t539 * t459;
t536 = t539 * t456 + t538 * t459;
t535 = Icges(3,5) * t455 + Icges(4,5) * t442 + Icges(3,6) * t458 + Icges(4,6) * t443;
t521 = Icges(4,4) * t442;
t420 = Icges(4,2) * t443 + t521;
t520 = Icges(4,4) * t443;
t421 = Icges(4,1) * t442 + t520;
t523 = Icges(3,4) * t455;
t429 = Icges(3,2) * t458 + t523;
t522 = Icges(3,4) * t458;
t430 = Icges(3,1) * t455 + t522;
t534 = -t420 * t442 + t421 * t443 - t429 * t455 + t430 * t458;
t483 = -Icges(4,2) * t442 + t520;
t387 = Icges(4,6) * t456 + t459 * t483;
t487 = Icges(4,1) * t443 - t521;
t389 = Icges(4,5) * t456 + t459 * t487;
t484 = -Icges(3,2) * t455 + t522;
t402 = Icges(3,6) * t456 + t459 * t484;
t488 = Icges(3,1) * t458 - t523;
t404 = Icges(3,5) * t456 + t459 * t488;
t533 = -t387 * t442 + t389 * t443 - t402 * t455 + t404 * t458;
t386 = -Icges(4,6) * t459 + t456 * t483;
t388 = -Icges(4,5) * t459 + t456 * t487;
t401 = -Icges(3,6) * t459 + t456 * t484;
t403 = -Icges(3,5) * t459 + t456 * t488;
t532 = t386 * t442 - t388 * t443 + t401 * t455 - t403 * t458;
t528 = pkin(2) * t455;
t444 = qJ(4) + t452;
t438 = sin(t444);
t527 = pkin(4) * t438;
t525 = t458 * pkin(2);
t519 = Icges(5,4) * t438;
t439 = cos(t444);
t518 = Icges(5,4) * t439;
t440 = qJ(5) + t444;
t435 = sin(t440);
t517 = Icges(6,4) * t435;
t436 = cos(t440);
t516 = Icges(6,4) * t436;
t515 = t435 * t456;
t514 = t435 * t459;
t454 = sin(qJ(6));
t513 = t454 * t456;
t512 = t454 * t459;
t457 = cos(qJ(6));
t511 = t456 * t457;
t510 = t457 * t459;
t382 = -qJ(3) * t459 + t456 * t525;
t383 = qJ(3) * t456 + t459 * t525;
t448 = qJD(2) * t456;
t503 = qJD(2) * t459;
t509 = t382 * t448 + t383 * t503;
t434 = pkin(1) * t456 - pkin(7) * t459;
t508 = -t382 - t434;
t507 = pkin(4) * t439;
t506 = pkin(3) * t443;
t426 = qJD(4) * t456 + t448;
t502 = qJD(6) * t435;
t501 = -qJD(2) - qJD(4);
t355 = -pkin(8) * t459 + t456 * t506;
t500 = -t355 + t508;
t417 = qJD(5) * t456 + t426;
t347 = -pkin(9) * t459 + t456 * t507;
t497 = -t347 + t500;
t356 = pkin(8) * t456 + t459 * t506;
t496 = t355 * t448 + t356 * t503 + t509;
t495 = pkin(5) * t436 + pkin(10) * t435;
t425 = qJD(1) * (pkin(1) * t459 + pkin(7) * t456);
t494 = qJD(1) * t383 - qJD(3) * t459 + t425;
t493 = rSges(3,1) * t458 - rSges(3,2) * t455;
t492 = rSges(4,1) * t443 - rSges(4,2) * t442;
t491 = rSges(5,1) * t439 - rSges(5,2) * t438;
t490 = rSges(6,1) * t436 - rSges(6,2) * t435;
t489 = qJD(2) * (-rSges(4,1) * t442 - rSges(4,2) * t443 - t528);
t418 = (-qJD(5) + t501) * t459;
t486 = Icges(5,1) * t439 - t519;
t485 = Icges(6,1) * t436 - t517;
t482 = -Icges(5,2) * t438 + t518;
t481 = -Icges(6,2) * t435 + t516;
t478 = Icges(5,5) * t439 - Icges(5,6) * t438;
t477 = Icges(6,5) * t436 - Icges(6,6) * t435;
t470 = qJD(2) * (-pkin(3) * t442 - t528);
t348 = pkin(9) * t456 + t459 * t507;
t427 = t501 * t459;
t469 = t426 * t347 - t348 * t427 + t496;
t468 = qJD(1) * (Icges(6,5) * t435 + Icges(6,6) * t436) + (-Icges(6,3) * t459 + t456 * t477) * t418 + (Icges(6,3) * t456 + t459 * t477) * t417;
t467 = qJD(1) * (Icges(5,5) * t438 + Icges(5,6) * t439) + (-Icges(5,3) * t459 + t456 * t478) * t427 + (Icges(5,3) * t456 + t459 * t478) * t426;
t447 = qJD(3) * t456;
t466 = t459 * t470 + t447;
t465 = t427 * t527 + t466;
t464 = qJD(1) * t356 + t456 * t470 + t494;
t463 = qJD(1) * t348 - t426 * t527 + t464;
t363 = -Icges(6,6) * t459 + t456 * t481;
t364 = Icges(6,6) * t456 + t459 * t481;
t365 = -Icges(6,5) * t459 + t456 * t485;
t366 = Icges(6,5) * t456 + t459 * t485;
t406 = Icges(6,2) * t436 + t517;
t407 = Icges(6,1) * t435 + t516;
t462 = (-t364 * t435 + t366 * t436) * t417 + (-t363 * t435 + t365 * t436) * t418 + (-t406 * t435 + t407 * t436) * qJD(1);
t375 = -Icges(5,6) * t459 + t456 * t482;
t376 = Icges(5,6) * t456 + t459 * t482;
t377 = -Icges(5,5) * t459 + t456 * t486;
t378 = Icges(5,5) * t456 + t459 * t486;
t414 = Icges(5,2) * t439 + t519;
t415 = Icges(5,1) * t438 + t518;
t461 = (-t376 * t438 + t378 * t439) * t426 + (-t375 * t438 + t377 * t439) * t427 + (-t414 * t438 + t415 * t439) * qJD(1);
t433 = rSges(2,1) * t459 - rSges(2,2) * t456;
t432 = rSges(2,1) * t456 + rSges(2,2) * t459;
t431 = rSges(3,1) * t455 + rSges(3,2) * t458;
t424 = -qJD(6) * t436 + qJD(1);
t416 = rSges(5,1) * t438 + rSges(5,2) * t439;
t411 = pkin(5) * t435 - pkin(10) * t436;
t410 = rSges(6,1) * t435 + rSges(6,2) * t436;
t409 = rSges(3,3) * t456 + t459 * t493;
t408 = -rSges(3,3) * t459 + t456 * t493;
t398 = t436 * t510 + t513;
t397 = -t436 * t512 + t511;
t396 = t436 * t511 - t512;
t395 = -t436 * t513 - t510;
t393 = rSges(4,3) * t456 + t459 * t492;
t392 = -rSges(4,3) * t459 + t456 * t492;
t391 = t495 * t459;
t390 = t495 * t456;
t381 = rSges(5,3) * t456 + t459 * t491;
t380 = -rSges(5,3) * t459 + t456 * t491;
t372 = t456 * t502 + t418;
t371 = t459 * t502 + t417;
t370 = rSges(6,3) * t456 + t459 * t490;
t369 = -rSges(6,3) * t459 + t456 * t490;
t360 = -rSges(7,3) * t436 + (rSges(7,1) * t457 - rSges(7,2) * t454) * t435;
t359 = -Icges(7,5) * t436 + (Icges(7,1) * t457 - Icges(7,4) * t454) * t435;
t358 = -Icges(7,6) * t436 + (Icges(7,4) * t457 - Icges(7,2) * t454) * t435;
t357 = -Icges(7,3) * t436 + (Icges(7,5) * t457 - Icges(7,6) * t454) * t435;
t351 = qJD(1) * t409 - t431 * t448 + t425;
t350 = -t431 * t503 + (-t408 - t434) * qJD(1);
t349 = (t408 * t456 + t409 * t459) * qJD(2);
t346 = rSges(7,1) * t398 + rSges(7,2) * t397 + rSges(7,3) * t514;
t345 = rSges(7,1) * t396 + rSges(7,2) * t395 + rSges(7,3) * t515;
t344 = Icges(7,1) * t398 + Icges(7,4) * t397 + Icges(7,5) * t514;
t343 = Icges(7,1) * t396 + Icges(7,4) * t395 + Icges(7,5) * t515;
t342 = Icges(7,4) * t398 + Icges(7,2) * t397 + Icges(7,6) * t514;
t341 = Icges(7,4) * t396 + Icges(7,2) * t395 + Icges(7,6) * t515;
t340 = Icges(7,5) * t398 + Icges(7,6) * t397 + Icges(7,3) * t514;
t339 = Icges(7,5) * t396 + Icges(7,6) * t395 + Icges(7,3) * t515;
t336 = qJD(1) * t393 + t456 * t489 + t494;
t335 = t447 + t459 * t489 + (-t392 + t508) * qJD(1);
t334 = (t392 * t456 + t393 * t459) * qJD(2) + t509;
t333 = qJD(1) * t381 - t416 * t426 + t464;
t332 = t416 * t427 + (-t380 + t500) * qJD(1) + t466;
t331 = t380 * t426 - t381 * t427 + t496;
t330 = qJD(1) * t370 - t410 * t417 + t463;
t329 = t410 * t418 + (-t369 + t497) * qJD(1) + t465;
t328 = t369 * t417 - t370 * t418 + t469;
t327 = qJD(1) * t391 + t346 * t424 - t360 * t371 - t411 * t417 + t463;
t326 = -t345 * t424 + t360 * t372 + t411 * t418 + (-t390 + t497) * qJD(1) + t465;
t325 = t345 * t371 - t346 * t372 + t390 * t417 - t391 * t418 + t469;
t1 = m(7) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(6) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(5) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + t424 * ((-t339 * t372 - t340 * t371 - t357 * t424) * t436 + ((-t342 * t454 + t344 * t457) * t371 + (-t341 * t454 + t343 * t457) * t372 + (-t358 * t454 + t359 * t457) * t424) * t435) / 0.2e1 + t371 * ((t340 * t514 + t397 * t342 + t398 * t344) * t371 + (t339 * t514 + t341 * t397 + t343 * t398) * t372 + (t357 * t514 + t358 * t397 + t359 * t398) * t424) / 0.2e1 + t372 * ((t340 * t515 + t342 * t395 + t344 * t396) * t371 + (t339 * t515 + t395 * t341 + t396 * t343) * t372 + (t357 * t515 + t358 * t395 + t359 * t396) * t424) / 0.2e1 + t417 * (t468 * t456 + t462 * t459) / 0.2e1 + t418 * (t462 * t456 - t468 * t459) / 0.2e1 + t426 * (t467 * t456 + t461 * t459) / 0.2e1 + t427 * (t461 * t456 - t467 * t459) / 0.2e1 + m(4) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(3) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + (m(2) * (t432 ^ 2 + t433 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t536 * t456 ^ 2 + (t532 * t459 + (t533 - t537) * t456) * t459) * qJD(2) + (t535 * t456 + t534 * t459) * qJD(1)) * t448 / 0.2e1 - ((t537 * t459 ^ 2 + (t533 * t456 + (t532 - t536) * t459) * t456) * qJD(2) + (t534 * t456 - t535 * t459) * qJD(1)) * t503 / 0.2e1 + ((t364 * t436 + t366 * t435) * t417 + (t363 * t436 + t365 * t435) * t418 + (t376 * t439 + t378 * t438) * t426 + (t375 * t439 + t377 * t438) * t427 + ((-t386 * t443 - t388 * t442 - t401 * t458 - t403 * t455) * t459 + (t387 * t443 + t389 * t442 + t402 * t458 + t404 * t455) * t456) * qJD(2) + (t436 * t406 + t435 * t407 + t439 * t414 + t438 * t415 + t443 * t420 + t442 * t421 + t458 * t429 + t455 * t430) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
