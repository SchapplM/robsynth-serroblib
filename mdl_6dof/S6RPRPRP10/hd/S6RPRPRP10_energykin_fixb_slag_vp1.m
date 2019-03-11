% Calculate kinetic energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:55
% EndTime: 2019-03-09 03:30:58
% DurationCPUTime: 2.32s
% Computational Cost: add. (649->201), mult. (1531->298), div. (0->0), fcn. (1475->6), ass. (0->111)
t520 = -Icges(4,4) - Icges(5,6);
t519 = Icges(4,1) + Icges(5,2);
t518 = Icges(4,2) + Icges(5,3);
t429 = cos(qJ(3));
t517 = t520 * t429;
t426 = sin(qJ(3));
t516 = t520 * t426;
t515 = -Icges(5,4) + Icges(4,5);
t514 = Icges(5,5) - Icges(4,6);
t513 = -t429 * t518 + t516;
t512 = t426 * t519 - t517;
t511 = Icges(5,1) + Icges(4,3);
t510 = Icges(6,1) + Icges(7,1);
t509 = Icges(6,4) - Icges(7,5);
t508 = Icges(7,4) + Icges(6,5);
t507 = Icges(6,2) + Icges(7,3);
t506 = Icges(7,6) - Icges(6,6);
t505 = Icges(6,3) + Icges(7,2);
t427 = sin(qJ(1));
t430 = cos(qJ(1));
t504 = t427 * t513 + t430 * t514;
t503 = -t427 * t514 + t430 * t513;
t502 = t512 * t427 + t430 * t515;
t501 = t427 * t515 - t512 * t430;
t500 = t426 * t518 + t517;
t499 = t429 * t519 + t516;
t498 = t426 * t515 - t514 * t429;
t497 = rSges(7,1) + pkin(5);
t496 = rSges(7,3) + qJ(6);
t425 = sin(qJ(5));
t428 = cos(qJ(5));
t465 = t429 * t430;
t392 = t425 * t427 - t428 * t465;
t393 = t425 * t465 + t427 * t428;
t467 = t426 * t430;
t495 = t507 * t392 - t509 * t393 - t506 * t467;
t466 = t427 * t429;
t394 = t425 * t430 + t428 * t466;
t395 = -t425 * t466 + t428 * t430;
t468 = t426 * t427;
t494 = t507 * t394 - t509 * t395 + t506 * t468;
t493 = t506 * t392 + t508 * t393 - t505 * t467;
t492 = t506 * t394 + t508 * t395 + t505 * t468;
t491 = -t509 * t392 + t510 * t393 - t508 * t467;
t490 = -t509 * t394 + t510 * t395 + t508 * t468;
t489 = t506 * t429 + (-t509 * t425 - t507 * t428) * t426;
t488 = t505 * t429 + (t508 * t425 - t506 * t428) * t426;
t487 = t508 * t429 + (t510 * t425 + t509 * t428) * t426;
t486 = t511 * t427 - t498 * t430;
t485 = t498 * t427 + t511 * t430;
t484 = t514 * t426 + t429 * t515;
t483 = t499 * t426 - t500 * t429;
t482 = -t502 * t426 + t504 * t429;
t481 = t501 * t426 + t503 * t429;
t474 = pkin(8) * t429;
t464 = -rSges(7,2) * t467 + t496 * t392 + t497 * t393;
t463 = rSges(7,2) * t468 + t496 * t394 + t497 * t395;
t447 = pkin(3) * t426 - qJ(4) * t429;
t398 = t447 * t430;
t457 = qJD(3) * t430;
t462 = qJD(4) * t426 - t398 * t457;
t461 = rSges(7,2) * t429 + (t497 * t425 - t496 * t428) * t426;
t414 = pkin(3) * t429 + qJ(4) * t426;
t424 = qJD(2) * t427;
t458 = qJD(3) * t427;
t460 = t414 * t458 + t424;
t405 = qJD(1) * (pkin(1) * t430 + qJ(2) * t427);
t459 = qJD(1) * t430 * pkin(7) + t405;
t456 = qJD(4) * t429;
t455 = qJD(5) * t426;
t412 = pkin(1) * t427 - qJ(2) * t430;
t452 = -pkin(7) * t427 - t412;
t397 = t447 * t427;
t451 = qJD(1) * t397 + t430 * t456 + t459;
t450 = t398 + t452;
t449 = rSges(4,1) * t426 + rSges(4,2) * t429;
t448 = rSges(5,2) * t426 + rSges(5,3) * t429;
t403 = pkin(4) * t427 - pkin(8) * t467;
t404 = pkin(4) * t430 + pkin(8) * t468;
t434 = t403 * t457 + (-t397 - t404) * t458 + t462;
t433 = qJD(1) * t404 + (-qJD(2) + (-t414 - t474) * qJD(3)) * t430 + t451;
t432 = t458 * t474 + (-t403 + t450) * qJD(1) - t427 * t456 + t460;
t421 = qJD(5) * t429 + qJD(1);
t417 = rSges(2,1) * t430 - rSges(2,2) * t427;
t416 = rSges(4,1) * t429 - rSges(4,2) * t426;
t415 = -rSges(5,2) * t429 + rSges(5,3) * t426;
t413 = rSges(2,1) * t427 + rSges(2,2) * t430;
t402 = t427 * t455 + t457;
t401 = -t430 * t455 + t458;
t389 = rSges(5,1) * t430 - t427 * t448;
t388 = rSges(5,1) * t427 + t430 * t448;
t387 = rSges(4,3) * t427 - t430 * t449;
t386 = rSges(4,3) * t430 + t427 * t449;
t385 = rSges(6,3) * t429 + (rSges(6,1) * t425 + rSges(6,2) * t428) * t426;
t364 = t405 - qJD(2) * t430 + qJD(1) * (-rSges(3,2) * t430 + rSges(3,3) * t427);
t363 = t424 + (rSges(3,2) * t427 + rSges(3,3) * t430 - t412) * qJD(1);
t360 = rSges(6,1) * t395 - rSges(6,2) * t394 + rSges(6,3) * t468;
t358 = rSges(6,1) * t393 - rSges(6,2) * t392 - rSges(6,3) * t467;
t344 = (-t386 * t427 + t387 * t430) * qJD(3);
t343 = qJD(1) * t386 + (-qJD(3) * t416 - qJD(2)) * t430 + t459;
t342 = t416 * t458 + t424 + (-t387 + t452) * qJD(1);
t341 = (t388 * t430 + (-t389 - t397) * t427) * qJD(3) + t462;
t340 = qJD(1) * t389 + (-qJD(2) + (-t414 - t415) * qJD(3)) * t430 + t451;
t339 = (qJD(3) * t415 - t456) * t427 + (-t388 + t450) * qJD(1) + t460;
t338 = t360 * t421 - t385 * t402 + t433;
t337 = -t358 * t421 + t385 * t401 + t432;
t336 = t358 * t402 - t360 * t401 + t434;
t335 = qJD(6) * t392 - t402 * t461 + t421 * t463 + t433;
t334 = qJD(6) * t394 + t401 * t461 - t421 * t464 + t432;
t333 = -qJD(6) * t426 * t428 - t401 * t463 + t402 * t464 + t434;
t1 = m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(3) * (t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(4) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + ((t489 * t392 + t487 * t393 - t488 * t467) * t421 + (t494 * t392 + t490 * t393 - t492 * t467) * t402 + (t495 * t392 + t491 * t393 - t493 * t467) * t401) * t401 / 0.2e1 + ((t489 * t394 + t487 * t395 + t488 * t468) * t421 + (t494 * t394 + t490 * t395 + t492 * t468) * t402 + (t495 * t394 + t491 * t395 + t493 * t468) * t401) * t402 / 0.2e1 + ((t493 * t401 + t492 * t402 + t488 * t421) * t429 + ((t487 * t425 - t489 * t428) * t421 + (t490 * t425 - t494 * t428) * t402 + (t491 * t425 - t495 * t428) * t401) * t426) * t421 / 0.2e1 + (((t504 * t426 + t502 * t429) * t430 + (-t503 * t426 + t501 * t429) * t427) * qJD(3) + (t500 * t426 + t499 * t429) * qJD(1)) * qJD(1) / 0.2e1 + ((t486 * t427 ^ 2 + (t482 * t430 + (-t481 + t485) * t427) * t430) * qJD(3) + (t484 * t427 - t483 * t430) * qJD(1)) * t458 / 0.2e1 + ((t485 * t430 ^ 2 + (t481 * t427 + (-t482 + t486) * t430) * t427) * qJD(3) + (t483 * t427 + t484 * t430) * qJD(1)) * t457 / 0.2e1 + (m(2) * (t413 ^ 2 + t417 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
