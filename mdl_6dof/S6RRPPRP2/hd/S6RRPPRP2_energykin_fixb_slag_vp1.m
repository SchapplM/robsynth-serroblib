% Calculate kinetic energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:37
% EndTime: 2019-03-09 08:29:40
% DurationCPUTime: 2.93s
% Computational Cost: add. (1360->237), mult. (1825->349), div. (0->0), fcn. (1726->8), ass. (0->138)
t579 = Icges(4,4) + Icges(5,6);
t578 = Icges(4,1) + Icges(5,2);
t577 = -Icges(4,2) - Icges(5,3);
t460 = qJ(2) + pkin(9);
t458 = cos(t460);
t576 = t579 * t458;
t457 = sin(t460);
t575 = t579 * t457;
t574 = -Icges(5,4) + Icges(4,5);
t573 = Icges(5,5) - Icges(4,6);
t572 = t577 * t457 + t576;
t571 = -t578 * t458 + t575;
t570 = Icges(6,1) + Icges(7,1);
t569 = Icges(6,4) + Icges(7,4);
t568 = Icges(7,5) + Icges(6,5);
t567 = Icges(6,2) + Icges(7,2);
t566 = Icges(7,6) + Icges(6,6);
t565 = Icges(7,3) + Icges(6,3);
t465 = sin(qJ(1));
t468 = cos(qJ(1));
t564 = t572 * t465 + t573 * t468;
t563 = -t573 * t465 + t572 * t468;
t562 = t571 * t465 + t574 * t468;
t561 = t574 * t465 - t571 * t468;
t560 = t577 * t458 - t575;
t559 = t578 * t457 + t576;
t558 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t464 = sin(qJ(2));
t467 = cos(qJ(2));
t557 = Icges(3,5) * t467 - Icges(3,6) * t464 + t573 * t457 + t574 * t458;
t466 = cos(qJ(5));
t517 = t466 * t468;
t463 = sin(qJ(5));
t520 = t463 * t465;
t426 = t457 * t517 - t520;
t518 = t465 * t466;
t519 = t463 * t468;
t427 = t457 * t519 + t518;
t521 = t458 * t468;
t556 = t566 * t426 + t568 * t427 + t565 * t521;
t428 = t457 * t518 + t519;
t429 = t457 * t520 - t517;
t522 = t458 * t465;
t555 = t566 * t428 + t568 * t429 + t565 * t522;
t554 = t567 * t426 + t569 * t427 + t566 * t521;
t553 = t567 * t428 + t569 * t429 + t566 * t522;
t552 = t569 * t426 + t570 * t427 + t568 * t521;
t551 = t569 * t428 + t570 * t429 + t568 * t522;
t550 = (-t568 * t463 - t566 * t466) * t458 + t565 * t457;
t549 = (-t569 * t463 - t567 * t466) * t458 + t566 * t457;
t548 = (-t570 * t463 - t569 * t466) * t458 + t568 * t457;
t547 = t557 * t465 - t558 * t468;
t546 = t558 * t465 + t557 * t468;
t545 = Icges(3,5) * t464 + Icges(3,6) * t467 + t574 * t457 - t573 * t458;
t528 = Icges(3,4) * t464;
t446 = Icges(3,2) * t467 + t528;
t527 = Icges(3,4) * t467;
t447 = Icges(3,1) * t464 + t527;
t544 = -t446 * t464 + t447 * t467 + t560 * t457 + t559 * t458;
t494 = -Icges(3,2) * t464 + t527;
t417 = Icges(3,6) * t465 + t468 * t494;
t496 = Icges(3,1) * t467 - t528;
t419 = Icges(3,5) * t465 + t468 * t496;
t543 = -t417 * t464 + t419 * t467 - t563 * t457 + t561 * t458;
t416 = -Icges(3,6) * t468 + t465 * t494;
t418 = -Icges(3,5) * t468 + t465 * t496;
t542 = t416 * t464 - t418 * t467 + t564 * t457 + t562 * t458;
t535 = pkin(2) * t464;
t534 = pkin(5) * t463;
t532 = pkin(2) * t467;
t531 = pkin(5) * t466;
t475 = qJ(6) * t458 + t457 * t534;
t516 = rSges(7,1) * t427 + rSges(7,2) * t426 + rSges(7,3) * t521 + t465 * t531 + t468 * t475;
t515 = rSges(7,1) * t429 + rSges(7,2) * t428 + rSges(7,3) * t522 + t465 * t475 - t468 * t531;
t393 = -qJ(3) * t468 + t465 * t532;
t394 = qJ(3) * t465 + t468 * t532;
t509 = qJD(2) * t468;
t510 = qJD(2) * t465;
t514 = t393 * t510 + t394 * t509;
t513 = (-rSges(7,1) * t463 - rSges(7,2) * t466 - t534) * t458 + (rSges(7,3) + qJ(6)) * t457;
t453 = pkin(1) * t465 - pkin(7) * t468;
t512 = -t393 - t453;
t459 = qJD(3) * t465;
t508 = qJD(4) * t457;
t511 = t468 * t508 + t459;
t507 = qJD(5) * t458;
t498 = pkin(3) * t458 + qJ(4) * t457;
t420 = t498 * t465;
t506 = -t420 + t512;
t503 = -pkin(3) * t457 + qJ(4) * t458 - t535;
t444 = qJD(1) * (pkin(1) * t468 + pkin(7) * t465);
t502 = qJD(1) * t394 - qJD(3) * t468 + t444;
t501 = rSges(3,1) * t467 - rSges(3,2) * t464;
t500 = rSges(4,1) * t458 - rSges(4,2) * t457;
t499 = -rSges(5,2) * t458 + rSges(5,3) * t457;
t497 = qJD(2) * (-rSges(4,1) * t457 - rSges(4,2) * t458 - t535);
t478 = qJD(2) * (rSges(5,2) * t457 + rSges(5,3) * t458 + t503);
t421 = t498 * t468;
t477 = qJD(1) * t421 + t465 * t508 + t502;
t476 = -qJD(4) * t458 + t420 * t510 + t421 * t509 + t514;
t474 = (-pkin(8) * t457 + t503) * qJD(2);
t442 = pkin(4) * t465 + pkin(8) * t521;
t473 = qJD(1) * t442 + t477;
t443 = -pkin(4) * t468 + pkin(8) * t522;
t472 = (-t443 + t506) * qJD(1) + t511;
t471 = t442 * t509 + t443 * t510 + t476;
t470 = qJD(6) * t458 + t474;
t454 = qJD(5) * t457 + qJD(1);
t452 = rSges(2,1) * t468 - rSges(2,2) * t465;
t451 = rSges(2,1) * t465 + rSges(2,2) * t468;
t450 = rSges(3,1) * t464 + rSges(3,2) * t467;
t432 = t465 * t507 - t509;
t431 = t468 * t507 + t510;
t423 = rSges(3,3) * t465 + t468 * t501;
t422 = -rSges(3,3) * t468 + t465 * t501;
t412 = -rSges(5,1) * t468 + t465 * t499;
t411 = rSges(5,1) * t465 + t468 * t499;
t410 = rSges(4,3) * t465 + t468 * t500;
t409 = -rSges(4,3) * t468 + t465 * t500;
t391 = rSges(6,3) * t457 + (-rSges(6,1) * t463 - rSges(6,2) * t466) * t458;
t378 = qJD(1) * t423 - t450 * t510 + t444;
t377 = -t450 * t509 + (-t422 - t453) * qJD(1);
t376 = (t422 * t465 + t423 * t468) * qJD(2);
t375 = rSges(6,1) * t429 + rSges(6,2) * t428 + rSges(6,3) * t522;
t373 = rSges(6,1) * t427 + rSges(6,2) * t426 + rSges(6,3) * t521;
t359 = qJD(1) * t410 + t465 * t497 + t502;
t358 = t459 + t468 * t497 + (-t409 + t512) * qJD(1);
t357 = (t409 * t465 + t410 * t468) * qJD(2) + t514;
t356 = qJD(1) * t411 + t465 * t478 + t477;
t355 = t468 * t478 + (-t412 + t506) * qJD(1) + t511;
t354 = (t411 * t468 + t412 * t465) * qJD(2) + t476;
t353 = t373 * t454 - t391 * t431 + t465 * t474 + t473;
t352 = -t375 * t454 + t391 * t432 + t468 * t474 + t472;
t351 = -t373 * t432 + t375 * t431 + t471;
t350 = -t431 * t513 + t454 * t516 + t465 * t470 + t473;
t349 = t432 * t513 - t454 * t515 + t468 * t470 + t472;
t348 = qJD(6) * t457 + t431 * t515 - t432 * t516 + t471;
t1 = m(3) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(4) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(5) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(6) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(7) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + ((t549 * t426 + t548 * t427 + t550 * t521) * t454 + (t553 * t426 + t551 * t427 + t555 * t521) * t432 + (t554 * t426 + t552 * t427 + t556 * t521) * t431) * t431 / 0.2e1 + ((t549 * t428 + t548 * t429 + t550 * t522) * t454 + (t553 * t428 + t551 * t429 + t555 * t522) * t432 + (t554 * t428 + t552 * t429 + t556 * t522) * t431) * t432 / 0.2e1 + (((-t548 * t463 - t549 * t466) * t454 + (-t551 * t463 - t553 * t466) * t432 + (-t552 * t463 - t554 * t466) * t431) * t458 + (t556 * t431 + t555 * t432 + t550 * t454) * t457) * t454 / 0.2e1 + (m(2) * (t451 ^ 2 + t452 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t416 * t467 - t418 * t464 + t562 * t457 - t564 * t458) * t468 + (t417 * t467 + t419 * t464 + t561 * t457 + t563 * t458) * t465) * qJD(2) + (t467 * t446 + t464 * t447 + t559 * t457 - t560 * t458) * qJD(1)) * qJD(1) / 0.2e1 + ((t546 * t465 ^ 2 + (t542 * t468 + (t543 - t547) * t465) * t468) * qJD(2) + (t545 * t465 + t544 * t468) * qJD(1)) * t510 / 0.2e1 - ((t547 * t468 ^ 2 + (t543 * t465 + (t542 - t546) * t468) * t465) * qJD(2) + (t544 * t465 - t545 * t468) * qJD(1)) * t509 / 0.2e1;
T  = t1;
