% Calculate kinetic energy for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:55
% EndTime: 2019-03-09 20:48:57
% DurationCPUTime: 2.21s
% Computational Cost: add. (1611->239), mult. (2172->362), div. (0->0), fcn. (2153->8), ass. (0->136)
t509 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t508 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t507 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t506 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t505 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t504 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t503 = rSges(7,1) + pkin(5);
t502 = rSges(7,3) + qJ(6);
t432 = qJ(2) + qJ(3);
t431 = cos(t432);
t435 = sin(qJ(1));
t436 = cos(qJ(4));
t473 = t435 * t436;
t433 = sin(qJ(4));
t438 = cos(qJ(1));
t474 = t433 * t438;
t406 = t431 * t474 - t473;
t472 = t436 * t438;
t475 = t433 * t435;
t407 = t431 * t472 + t475;
t358 = pkin(4) * t407 + qJ(5) * t406;
t404 = t431 * t475 + t472;
t425 = -qJD(4) * t431 + qJD(1);
t501 = qJD(5) * t404 + t425 * t358;
t416 = (-qJD(2) - qJD(3)) * t438;
t430 = sin(t432);
t464 = qJD(4) * t430;
t396 = t435 * t464 + t416;
t401 = (pkin(4) * t436 + qJ(5) * t433) * t430;
t500 = qJD(5) * t406 + t396 * t401;
t461 = pkin(3) * t431 + pkin(9) * t430;
t403 = t461 * t438;
t413 = pkin(3) * t430 - pkin(9) * t431;
t429 = qJD(2) * t435;
t415 = qJD(3) * t435 + t429;
t499 = qJD(1) * t403 - t413 * t415;
t405 = t431 * t473 - t474;
t477 = t430 * t435;
t498 = t505 * t404 + t507 * t405 - t504 * t477;
t476 = t430 * t438;
t497 = t505 * t406 + t507 * t407 - t504 * t476;
t496 = t506 * t404 + t508 * t405 - t505 * t477;
t495 = t506 * t406 + t508 * t407 - t505 * t476;
t494 = t508 * t404 + t509 * t405 - t507 * t477;
t493 = t508 * t406 + t509 * t407 - t507 * t476;
t492 = t504 * t431 + (t505 * t433 + t507 * t436) * t430;
t491 = t505 * t431 + (t506 * t433 + t508 * t436) * t430;
t490 = t507 * t431 + (t508 * t433 + t509 * t436) * t430;
t437 = cos(qJ(2));
t484 = pkin(2) * t437;
t434 = sin(qJ(2));
t482 = Icges(3,4) * t434;
t481 = Icges(3,4) * t437;
t480 = Icges(4,4) * t430;
t479 = Icges(4,4) * t431;
t471 = rSges(7,2) * t404 + t405 * t503 - t502 * t477;
t470 = rSges(7,2) * t406 + t407 * t503 - t502 * t476;
t377 = -pkin(8) * t438 + t435 * t484;
t378 = pkin(8) * t435 + t438 * t484;
t465 = qJD(2) * t438;
t469 = t377 * t429 + t378 * t465;
t414 = qJD(1) * (pkin(1) * t438 + pkin(7) * t435);
t468 = qJD(1) * t378 + t414;
t467 = t502 * t431 + (rSges(7,2) * t433 + t436 * t503) * t430;
t424 = pkin(1) * t435 - pkin(7) * t438;
t466 = -t377 - t424;
t463 = pkin(2) * qJD(2) * t434;
t462 = t438 * t463;
t460 = rSges(3,1) * t437 - rSges(3,2) * t434;
t459 = rSges(4,1) * t431 - rSges(4,2) * t430;
t458 = Icges(3,1) * t437 - t482;
t457 = Icges(4,1) * t431 - t480;
t456 = -Icges(3,2) * t434 + t481;
t455 = -Icges(4,2) * t430 + t479;
t454 = Icges(3,5) * t437 - Icges(3,6) * t434;
t453 = Icges(4,5) * t431 - Icges(4,6) * t430;
t391 = -Icges(3,6) * t438 + t435 * t456;
t393 = -Icges(3,5) * t438 + t435 * t458;
t452 = t391 * t434 - t393 * t437;
t392 = Icges(3,6) * t435 + t438 * t456;
t394 = Icges(3,5) * t435 + t438 * t458;
t451 = -t392 * t434 + t394 * t437;
t418 = Icges(3,2) * t437 + t482;
t419 = Icges(3,1) * t434 + t481;
t450 = -t418 * t434 + t419 * t437;
t402 = t461 * t435;
t449 = t415 * t402 - t403 * t416 + t469;
t448 = -qJD(6) * t430 - t463;
t447 = t416 * t413 + (-t402 + t466) * qJD(1);
t446 = -t435 * t463 + t468;
t357 = pkin(4) * t405 + qJ(5) * t404;
t395 = t438 * t464 + t415;
t445 = qJD(5) * t430 * t433 + t395 * t357 + t449;
t444 = (Icges(4,5) * t430 + Icges(4,6) * t431) * qJD(1) + (-Icges(4,3) * t438 + t435 * t453) * t416 + (Icges(4,3) * t435 + t438 * t453) * t415;
t443 = t446 + t499;
t442 = t447 - t462;
t382 = -Icges(4,6) * t438 + t435 * t455;
t383 = Icges(4,6) * t435 + t438 * t455;
t384 = -Icges(4,5) * t438 + t435 * t457;
t385 = Icges(4,5) * t435 + t438 * t457;
t410 = Icges(4,2) * t431 + t480;
t411 = Icges(4,1) * t430 + t479;
t441 = (-t383 * t430 + t385 * t431) * t415 + (-t382 * t430 + t384 * t431) * t416 + (-t410 * t430 + t411 * t431) * qJD(1);
t422 = rSges(2,1) * t438 - rSges(2,2) * t435;
t421 = rSges(2,1) * t435 + rSges(2,2) * t438;
t420 = rSges(3,1) * t434 + rSges(3,2) * t437;
t417 = Icges(3,5) * t434 + Icges(3,6) * t437;
t412 = rSges(4,1) * t430 + rSges(4,2) * t431;
t400 = rSges(3,3) * t435 + t438 * t460;
t399 = -rSges(3,3) * t438 + t435 * t460;
t390 = Icges(3,3) * t435 + t438 * t454;
t389 = -Icges(3,3) * t438 + t435 * t454;
t387 = rSges(4,3) * t435 + t438 * t459;
t386 = -rSges(4,3) * t438 + t435 * t459;
t376 = -rSges(5,3) * t431 + (rSges(5,1) * t436 - rSges(5,2) * t433) * t430;
t375 = -rSges(6,2) * t431 + (rSges(6,1) * t436 + rSges(6,3) * t433) * t430;
t355 = qJD(1) * t400 - t420 * t429 + t414;
t354 = -t420 * t465 + (-t399 - t424) * qJD(1);
t353 = rSges(5,1) * t407 - rSges(5,2) * t406 + rSges(5,3) * t476;
t352 = rSges(6,1) * t407 + rSges(6,2) * t476 + rSges(6,3) * t406;
t350 = rSges(5,1) * t405 - rSges(5,2) * t404 + rSges(5,3) * t477;
t349 = rSges(6,1) * t405 + rSges(6,2) * t477 + rSges(6,3) * t404;
t328 = (t399 * t435 + t400 * t438) * qJD(2);
t326 = qJD(1) * t387 - t412 * t415 + t446;
t325 = -t462 + t412 * t416 + (-t386 + t466) * qJD(1);
t324 = t386 * t415 - t387 * t416 + t469;
t323 = t353 * t425 - t376 * t395 + t443;
t322 = -t350 * t425 + t376 * t396 + t442;
t321 = t350 * t395 - t353 * t396 + t449;
t320 = t352 * t425 + (-t375 - t401) * t395 + t443 + t501;
t319 = t375 * t396 + (-t349 - t357) * t425 + t442 + t500;
t318 = t349 * t395 + (-t352 - t358) * t396 + t445;
t317 = t448 * t435 + t470 * t425 + (-t401 - t467) * t395 + t468 + t499 + t501;
t316 = t448 * t438 + t467 * t396 + (-t357 - t471) * t425 + t447 + t500;
t315 = qJD(6) * t431 + t471 * t395 + (-t358 - t470) * t396 + t445;
t1 = ((t417 * t435 + t438 * t450) * qJD(1) + (t390 * t435 ^ 2 + (t452 * t438 + (-t389 + t451) * t435) * t438) * qJD(2)) * t429 / 0.2e1 - ((-t417 * t438 + t435 * t450) * qJD(1) + (t389 * t438 ^ 2 + (t451 * t435 + (-t390 + t452) * t438) * t435) * qJD(2)) * t465 / 0.2e1 + t415 * (t444 * t435 + t441 * t438) / 0.2e1 + t416 * (t441 * t435 - t444 * t438) / 0.2e1 + m(6) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(7) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(5) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(4) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(3) * (t328 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t421 ^ 2 + t422 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t383 * t431 + t385 * t430) * t415 + (t382 * t431 + t384 * t430) * t416 + ((t392 * t437 + t394 * t434) * t435 - (t391 * t437 + t393 * t434) * t438) * qJD(2) + (t410 * t431 + t411 * t430 + t418 * t437 + t419 * t434) * qJD(1)) * qJD(1) / 0.2e1 + ((t491 * t406 + t490 * t407 - t492 * t476) * t425 + (t496 * t406 + t494 * t407 - t498 * t476) * t396 + (t495 * t406 + t493 * t407 - t497 * t476) * t395) * t395 / 0.2e1 + ((t491 * t404 + t490 * t405 - t492 * t477) * t425 + (t496 * t404 + t494 * t405 - t498 * t477) * t396 + (t495 * t404 + t493 * t405 - t497 * t477) * t395) * t396 / 0.2e1 + ((t497 * t395 + t498 * t396 + t492 * t425) * t431 + ((t491 * t433 + t490 * t436) * t425 + (t496 * t433 + t494 * t436) * t396 + (t495 * t433 + t493 * t436) * t395) * t430) * t425 / 0.2e1;
T  = t1;
