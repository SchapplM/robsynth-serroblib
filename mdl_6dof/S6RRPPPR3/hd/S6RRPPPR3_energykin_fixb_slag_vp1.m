% Calculate kinetic energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:50
% EndTime: 2019-03-09 08:13:53
% DurationCPUTime: 3.17s
% Computational Cost: add. (898->259), mult. (1785->390), div. (0->0), fcn. (1686->8), ass. (0->138)
t549 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t548 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t547 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t456 = sin(qJ(2));
t546 = t549 * t456;
t458 = cos(qJ(2));
t545 = t549 * t458;
t544 = Icges(4,4) + Icges(3,5) + Icges(5,6);
t543 = -Icges(5,5) - Icges(3,6) + Icges(4,6);
t542 = t547 * t456 - t545;
t541 = -t548 * t458 + t546;
t540 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t457 = sin(qJ(1));
t459 = cos(qJ(1));
t539 = -t457 * t542 + t459 * t543;
t538 = t457 * t543 + t459 * t542;
t537 = t541 * t457 + t459 * t544;
t536 = t457 * t544 - t541 * t459;
t535 = -t547 * t458 - t546;
t534 = t548 * t456 + t545;
t533 = t543 * t456 + t458 * t544;
t532 = t540 * t457 + t533 * t459;
t531 = t533 * t457 - t540 * t459;
t530 = -t456 * t544 + t543 * t458;
t529 = t535 * t456 + t534 * t458;
t528 = t538 * t456 + t536 * t458;
t527 = t539 * t456 + t537 * t458;
t454 = cos(pkin(9));
t522 = pkin(5) * t454;
t453 = sin(pkin(9));
t515 = t453 * t457;
t514 = t453 * t459;
t513 = t456 * t457;
t512 = t456 * t459;
t511 = t457 * t458;
t510 = t458 * t459;
t488 = pkin(2) * t458 + qJ(3) * t456;
t416 = t488 * t457;
t441 = pkin(1) * t457 - pkin(7) * t459;
t508 = -t416 - t441;
t507 = qJD(2) * t457;
t506 = qJD(2) * t459;
t505 = qJD(3) * t456;
t504 = qJD(5) * t458;
t503 = qJD(6) * t458;
t418 = t488 * t459;
t424 = qJD(1) * (pkin(1) * t459 + pkin(7) * t457);
t502 = qJD(1) * t418 + t457 * t505 + t424;
t422 = pkin(3) * t511 + qJ(4) * t459;
t501 = -t422 + t508;
t434 = pkin(2) * t456 - qJ(3) * t458;
t498 = -pkin(3) * t456 - t434;
t497 = qJD(2) * (-rSges(4,1) * t456 + rSges(4,3) * t458 - t434);
t445 = t459 * t505;
t496 = -qJD(4) * t457 + t445;
t487 = pkin(4) * t456 + qJ(5) * t458;
t415 = t487 * t457;
t495 = -t415 + t501;
t494 = pkin(4) * t458 - qJ(5) * t456 + t498;
t493 = -qJD(3) * t458 + t416 * t507 + t418 * t506;
t492 = t459 * t504 + t496;
t491 = rSges(3,1) * t458 - rSges(3,2) * t456;
t490 = rSges(4,1) * t458 + rSges(4,3) * t456;
t489 = rSges(5,1) * t456 - rSges(5,2) * t458;
t423 = pkin(3) * t510 - qJ(4) * t457;
t486 = qJD(1) * t423 + qJD(4) * t459 + t502;
t467 = qJD(2) * (rSges(5,1) * t458 + rSges(5,2) * t456 + t498);
t466 = t422 * t507 + t423 * t506 + t493;
t465 = qJD(2) * (-pkin(8) * t456 + t458 * t522 + t494);
t464 = qJD(2) * (-rSges(6,3) * t456 - (-rSges(6,1) * t454 + rSges(6,2) * t453) * t458 + t494);
t417 = t487 * t459;
t463 = qJD(1) * t417 + t457 * t504 + t486;
t462 = pkin(8) * t458 + t456 * t522;
t461 = qJD(5) * t456 + t415 * t507 + t417 * t506 + t466;
t452 = pkin(9) + qJ(6);
t449 = cos(t452);
t448 = sin(t452);
t446 = qJD(6) * t456 + qJD(1);
t440 = rSges(2,1) * t459 - rSges(2,2) * t457;
t437 = rSges(2,1) * t457 + rSges(2,2) * t459;
t436 = rSges(3,1) * t456 + rSges(3,2) * t458;
t421 = t457 * t503 - t506;
t420 = t459 * t503 + t507;
t414 = t454 * t512 - t515;
t413 = -t453 * t512 - t454 * t457;
t412 = t454 * t513 + t514;
t411 = -t453 * t513 + t454 * t459;
t406 = rSges(3,3) * t457 + t459 * t491;
t405 = rSges(4,2) * t457 + t459 * t490;
t404 = -rSges(5,3) * t457 + t459 * t489;
t403 = -rSges(3,3) * t459 + t457 * t491;
t402 = -rSges(4,2) * t459 + t457 * t490;
t401 = rSges(5,3) * t459 + t457 * t489;
t378 = -t448 * t457 + t449 * t512;
t377 = -t448 * t512 - t449 * t457;
t376 = t448 * t459 + t449 * t513;
t375 = -t448 * t513 + t449 * t459;
t373 = Icges(6,5) * t456 + (-Icges(6,1) * t454 + Icges(6,4) * t453) * t458;
t372 = Icges(6,6) * t456 + (-Icges(6,4) * t454 + Icges(6,2) * t453) * t458;
t371 = Icges(6,3) * t456 + (-Icges(6,5) * t454 + Icges(6,6) * t453) * t458;
t370 = rSges(7,3) * t456 + (-rSges(7,1) * t449 + rSges(7,2) * t448) * t458;
t369 = Icges(7,5) * t456 + (-Icges(7,1) * t449 + Icges(7,4) * t448) * t458;
t368 = Icges(7,6) * t456 + (-Icges(7,4) * t449 + Icges(7,2) * t448) * t458;
t367 = Icges(7,3) * t456 + (-Icges(7,5) * t449 + Icges(7,6) * t448) * t458;
t365 = -pkin(5) * t515 + t459 * t462;
t364 = pkin(5) * t514 + t457 * t462;
t363 = rSges(6,1) * t414 + rSges(6,2) * t413 + rSges(6,3) * t510;
t362 = rSges(6,1) * t412 + rSges(6,2) * t411 + rSges(6,3) * t511;
t361 = Icges(6,1) * t414 + Icges(6,4) * t413 + Icges(6,5) * t510;
t360 = Icges(6,1) * t412 + Icges(6,4) * t411 + Icges(6,5) * t511;
t359 = Icges(6,4) * t414 + Icges(6,2) * t413 + Icges(6,6) * t510;
t358 = Icges(6,4) * t412 + Icges(6,2) * t411 + Icges(6,6) * t511;
t357 = Icges(6,5) * t414 + Icges(6,6) * t413 + Icges(6,3) * t510;
t356 = Icges(6,5) * t412 + Icges(6,6) * t411 + Icges(6,3) * t511;
t355 = qJD(1) * t406 - t436 * t507 + t424;
t354 = -t436 * t506 + (-t403 - t441) * qJD(1);
t353 = (t403 * t457 + t406 * t459) * qJD(2);
t352 = rSges(7,1) * t378 + rSges(7,2) * t377 + rSges(7,3) * t510;
t351 = rSges(7,1) * t376 + rSges(7,2) * t375 + rSges(7,3) * t511;
t350 = Icges(7,1) * t378 + Icges(7,4) * t377 + Icges(7,5) * t510;
t349 = Icges(7,1) * t376 + Icges(7,4) * t375 + Icges(7,5) * t511;
t348 = Icges(7,4) * t378 + Icges(7,2) * t377 + Icges(7,6) * t510;
t347 = Icges(7,4) * t376 + Icges(7,2) * t375 + Icges(7,6) * t511;
t346 = Icges(7,5) * t378 + Icges(7,6) * t377 + Icges(7,3) * t510;
t345 = Icges(7,5) * t376 + Icges(7,6) * t375 + Icges(7,3) * t511;
t344 = qJD(1) * t405 + t457 * t497 + t502;
t343 = t445 + t459 * t497 + (-t402 + t508) * qJD(1);
t342 = (t402 * t457 + t405 * t459) * qJD(2) + t493;
t341 = qJD(1) * t404 + t457 * t467 + t486;
t340 = t459 * t467 + (-t401 + t501) * qJD(1) + t496;
t339 = (t401 * t457 + t404 * t459) * qJD(2) + t466;
t338 = qJD(1) * t363 + t457 * t464 + t463;
t337 = t459 * t464 + (-t362 + t495) * qJD(1) + t492;
t336 = (t362 * t457 + t363 * t459) * qJD(2) + t461;
t335 = qJD(1) * t365 + t352 * t446 - t370 * t420 + t457 * t465 + t463;
t334 = -t351 * t446 + t370 * t421 + t459 * t465 + (-t364 + t495) * qJD(1) + t492;
t333 = t351 * t420 - t352 * t421 + (t364 * t457 + t365 * t459) * qJD(2) + t461;
t1 = m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(4) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(3) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + t420 * ((t346 * t510 + t377 * t348 + t378 * t350) * t420 + (t345 * t510 + t347 * t377 + t349 * t378) * t421 + (t367 * t510 + t368 * t377 + t369 * t378) * t446) / 0.2e1 + t446 * ((t345 * t421 + t346 * t420 + t367 * t446) * t456 + ((t348 * t448 - t350 * t449) * t420 + (t347 * t448 - t349 * t449) * t421 + (t368 * t448 - t369 * t449) * t446) * t458) / 0.2e1 + t421 * ((t346 * t511 + t348 * t375 + t350 * t376) * t420 + (t345 * t511 + t375 * t347 + t376 * t349) * t421 + (t367 * t511 + t368 * t375 + t369 * t376) * t446) / 0.2e1 + (m(2) * (t437 ^ 2 + t440 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((-t358 * t453 + t360 * t454 - t539) * t458 + (-t356 + t537) * t456) * t459 + ((t359 * t453 - t361 * t454 - t538) * t458 + (t357 + t536) * t456) * t457) * qJD(2) + ((t372 * t453 - t373 * t454 - t535) * t458 + (t371 + t534) * t456) * qJD(1)) * qJD(1) / 0.2e1 + (((-t356 * t510 - t358 * t413 - t360 * t414 + t527 * t459) * t459 + (t357 * t510 + t413 * t359 + t414 * t361 + (t528 - t531) * t459 + t532 * t457) * t457) * qJD(2) + (t371 * t510 + t372 * t413 + t373 * t414 - t457 * t530 + t459 * t529) * qJD(1)) * t507 / 0.2e1 - (((t357 * t511 + t359 * t411 + t361 * t412 + t528 * t457) * t457 + (-t356 * t511 - t411 * t358 - t412 * t360 + (t527 - t532) * t457 + t531 * t459) * t459) * qJD(2) + (t371 * t511 + t372 * t411 + t373 * t412 + t457 * t529 + t459 * t530) * qJD(1)) * t506 / 0.2e1;
T  = t1;
