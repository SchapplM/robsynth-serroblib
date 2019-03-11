% Calculate kinetic energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:32
% EndTime: 2019-03-09 08:33:35
% DurationCPUTime: 2.63s
% Computational Cost: add. (779->212), mult. (1820->318), div. (0->0), fcn. (1721->6), ass. (0->124)
t565 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t564 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t563 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t456 = cos(qJ(2));
t562 = t565 * t456;
t453 = sin(qJ(2));
t561 = t565 * t453;
t560 = Icges(4,4) + Icges(3,5) + Icges(5,6);
t559 = Icges(5,5) + Icges(3,6) - Icges(4,6);
t558 = t453 * t563 - t562;
t557 = -t456 * t564 + t561;
t556 = Icges(6,1) + Icges(7,1);
t555 = Icges(6,4) + Icges(7,4);
t554 = Icges(7,5) + Icges(6,5);
t553 = Icges(6,2) + Icges(7,2);
t552 = Icges(7,6) + Icges(6,6);
t551 = Icges(7,3) + Icges(6,3);
t550 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t454 = sin(qJ(1));
t457 = cos(qJ(1));
t549 = t454 * t558 + t457 * t559;
t548 = -t454 * t559 + t457 * t558;
t547 = t454 * t557 + t457 * t560;
t546 = t454 * t560 - t457 * t557;
t545 = -t456 * t563 - t561;
t544 = t453 * t564 + t562;
t543 = -t453 * t559 + t456 * t560;
t455 = cos(qJ(5));
t508 = t455 * t457;
t452 = sin(qJ(5));
t512 = t452 * t454;
t415 = -t453 * t512 + t508;
t510 = t454 * t455;
t511 = t452 * t457;
t416 = t453 * t510 + t511;
t509 = t454 * t456;
t542 = t415 * t552 + t416 * t554 + t509 * t551;
t417 = -t453 * t511 - t510;
t418 = t453 * t508 - t512;
t507 = t456 * t457;
t541 = t417 * t552 + t418 * t554 + t507 * t551;
t540 = t415 * t553 + t416 * t555 + t509 * t552;
t539 = t417 * t553 + t418 * t555 + t507 * t552;
t538 = t415 * t555 + t416 * t556 + t509 * t554;
t537 = t417 * t555 + t418 * t556 + t507 * t554;
t536 = (t452 * t552 - t455 * t554) * t456 + t551 * t453;
t535 = (t452 * t553 - t455 * t555) * t456 + t552 * t453;
t534 = (t452 * t555 - t455 * t556) * t456 + t554 * t453;
t533 = t454 * t543 - t457 * t550;
t532 = t454 * t550 + t457 * t543;
t531 = -t453 * t560 - t456 * t559;
t530 = t453 * t545 + t456 * t544;
t529 = t453 * t548 + t456 * t546;
t528 = -t453 * t549 + t456 * t547;
t520 = pkin(5) * t455;
t462 = qJ(6) * t456 + t453 * t520;
t506 = rSges(7,1) * t416 + rSges(7,2) * t415 + rSges(7,3) * t509 + pkin(5) * t511 + t454 * t462;
t505 = rSges(7,1) * t418 + rSges(7,2) * t417 + rSges(7,3) * t507 - pkin(5) * t512 + t457 * t462;
t504 = (-rSges(7,1) * t455 + rSges(7,2) * t452 - t520) * t456 + (qJ(6) + rSges(7,3)) * t453;
t486 = pkin(2) * t456 + qJ(3) * t453;
t419 = t486 * t454;
t444 = pkin(1) * t454 - pkin(7) * t457;
t503 = -t419 - t444;
t502 = qJD(2) * t454;
t501 = qJD(2) * t457;
t500 = qJD(3) * t453;
t499 = qJD(5) * t456;
t420 = t486 * t457;
t428 = qJD(1) * (pkin(1) * t457 + pkin(7) * t454);
t498 = qJD(1) * t420 + t454 * t500 + t428;
t426 = pkin(3) * t509 + qJ(4) * t457;
t497 = -t426 + t503;
t438 = pkin(2) * t453 - qJ(3) * t456;
t494 = -pkin(3) * t453 - t438;
t493 = qJD(2) * (-rSges(4,1) * t453 + rSges(4,3) * t456 - t438);
t447 = t457 * t500;
t492 = -qJD(4) * t454 + t447;
t491 = pkin(4) * t453 + pkin(8) * t456;
t490 = -qJD(3) * t456 + t419 * t502 + t420 * t501;
t489 = rSges(3,1) * t456 - rSges(3,2) * t453;
t488 = rSges(4,1) * t456 + rSges(4,3) * t453;
t487 = rSges(5,1) * t453 - rSges(5,2) * t456;
t427 = pkin(3) * t507 - qJ(4) * t454;
t485 = qJD(1) * t427 + qJD(4) * t457 + t498;
t466 = qJD(2) * (rSges(5,1) * t456 + rSges(5,2) * t453 + t494);
t465 = (pkin(4) * t456 - pkin(8) * t453 + t494) * qJD(2);
t422 = t491 * t457;
t464 = qJD(1) * t422 + t485;
t463 = t426 * t502 + t427 * t501 + t490;
t421 = t491 * t454;
t461 = t421 * t502 + t422 * t501 + t463;
t460 = qJD(6) * t456 + t465;
t459 = (-t421 + t497) * qJD(1) + t492;
t448 = qJD(5) * t453 + qJD(1);
t443 = rSges(2,1) * t457 - rSges(2,2) * t454;
t441 = rSges(2,1) * t454 + rSges(2,2) * t457;
t440 = rSges(3,1) * t453 + rSges(3,2) * t456;
t425 = t454 * t499 - t501;
t424 = t457 * t499 + t502;
t408 = rSges(3,3) * t454 + t457 * t489;
t407 = rSges(4,2) * t454 + t457 * t488;
t406 = -rSges(5,3) * t454 + t457 * t487;
t405 = rSges(6,3) * t453 + (-rSges(6,1) * t455 + rSges(6,2) * t452) * t456;
t403 = -rSges(3,3) * t457 + t454 * t489;
t402 = -rSges(4,2) * t457 + t454 * t488;
t401 = rSges(5,3) * t457 + t454 * t487;
t373 = rSges(6,1) * t418 + rSges(6,2) * t417 + rSges(6,3) * t507;
t371 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t509;
t355 = qJD(1) * t408 - t440 * t502 + t428;
t354 = -t440 * t501 + (-t403 - t444) * qJD(1);
t353 = (t403 * t454 + t408 * t457) * qJD(2);
t352 = qJD(1) * t407 + t454 * t493 + t498;
t351 = t447 + t457 * t493 + (-t402 + t503) * qJD(1);
t350 = (t402 * t454 + t407 * t457) * qJD(2) + t490;
t349 = qJD(1) * t406 + t454 * t466 + t485;
t348 = t457 * t466 + (-t401 + t497) * qJD(1) + t492;
t347 = (t401 * t454 + t406 * t457) * qJD(2) + t463;
t346 = t373 * t448 - t405 * t424 + t454 * t465 + t464;
t345 = -t371 * t448 + t405 * t425 + t457 * t465 + t459;
t344 = t371 * t424 - t373 * t425 + t461;
t343 = -t424 * t504 + t448 * t505 + t454 * t460 + t464;
t342 = t425 * t504 - t448 * t506 + t457 * t460 + t459;
t341 = qJD(6) * t453 + t424 * t506 - t425 * t505 + t461;
t1 = m(7) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(5) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(6) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(4) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + ((t417 * t535 + t418 * t534 + t507 * t536) * t448 + (t417 * t540 + t418 * t538 + t507 * t542) * t425 + (t539 * t417 + t537 * t418 + t541 * t507) * t424) * t424 / 0.2e1 + ((t415 * t535 + t416 * t534 + t509 * t536) * t448 + (t540 * t415 + t538 * t416 + t542 * t509) * t425 + (t415 * t539 + t416 * t537 + t509 * t541) * t424) * t425 / 0.2e1 + (((t452 * t535 - t455 * t534) * t448 + (t452 * t540 - t455 * t538) * t425 + (t452 * t539 - t455 * t537) * t424) * t456 + (t424 * t541 + t425 * t542 + t448 * t536) * t453) * t448 / 0.2e1 + (Icges(2,3) + m(2) * (t441 ^ 2 + t443 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t453 * t547 + t456 * t549) * t457 + (t453 * t546 - t456 * t548) * t454) * qJD(2) + (t544 * t453 - t545 * t456) * qJD(1)) * qJD(1) / 0.2e1 + ((t532 * t454 ^ 2 + (t528 * t457 + (t529 - t533) * t454) * t457) * qJD(2) + (-t454 * t531 + t457 * t530) * qJD(1)) * t502 / 0.2e1 - ((t533 * t457 ^ 2 + (t529 * t454 + (t528 - t532) * t457) * t454) * qJD(2) + (t454 * t530 + t457 * t531) * qJD(1)) * t501 / 0.2e1;
T  = t1;
