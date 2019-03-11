% Calculate kinetic energy for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:00
% EndTime: 2019-03-09 04:35:02
% DurationCPUTime: 1.72s
% Computational Cost: add. (1482->195), mult. (1811->296), div. (0->0), fcn. (1826->8), ass. (0->105)
t448 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t447 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t446 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t445 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t444 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t443 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t442 = rSges(7,1) + pkin(5);
t441 = rSges(7,3) + qJ(6);
t390 = qJ(1) + pkin(9);
t388 = sin(t390);
t389 = cos(t390);
t394 = cos(qJ(4));
t391 = sin(qJ(4));
t395 = cos(qJ(3));
t422 = t391 * t395;
t352 = t388 * t422 + t389 * t394;
t420 = t394 * t395;
t353 = t388 * t420 - t389 * t391;
t392 = sin(qJ(3));
t424 = t388 * t392;
t440 = t447 * t352 + t448 * t353 - t446 * t424;
t354 = -t388 * t394 + t389 * t422;
t355 = t388 * t391 + t389 * t420;
t423 = t389 * t392;
t439 = t447 * t354 + t448 * t355 - t446 * t423;
t438 = t445 * t352 + t447 * t353 - t444 * t424;
t437 = t445 * t354 + t447 * t355 - t444 * t423;
t436 = -t444 * t352 - t446 * t353 - t443 * t424;
t435 = -t444 * t354 - t446 * t355 - t443 * t423;
t434 = t443 * t395 + (-t444 * t391 - t446 * t394) * t392;
t433 = t444 * t395 + (t445 * t391 + t447 * t394) * t392;
t432 = t446 * t395 + (t447 * t391 + t448 * t394) * t392;
t393 = sin(qJ(1));
t427 = pkin(1) * t393;
t426 = Icges(4,4) * t392;
t425 = Icges(4,4) * t395;
t421 = t392 * t394;
t419 = rSges(7,2) * t352 + t441 * t353 + t442 * t424;
t418 = rSges(7,2) * t354 + t441 * t355 + t442 * t423;
t417 = (rSges(7,2) * t391 + rSges(7,3) * t394) * t392 + qJ(6) * t421 - t442 * t395;
t396 = cos(qJ(1));
t387 = qJD(1) * t396 * pkin(1);
t416 = qJD(1) * (pkin(2) * t389 + pkin(7) * t388) + t387;
t415 = qJD(3) * t388;
t414 = qJD(3) * t389;
t413 = qJD(4) * t392;
t410 = pkin(3) * t395 + pkin(8) * t392;
t368 = t410 * t388;
t369 = t410 * t389;
t412 = t368 * t415 + t369 * t414 + qJD(2);
t411 = -pkin(2) * t388 + pkin(7) * t389 - t427;
t409 = rSges(4,1) * t395 - rSges(4,2) * t392;
t332 = pkin(4) * t353 + qJ(5) * t352;
t370 = t389 * t413 + t415;
t408 = qJD(5) * t392 * t391 + t370 * t332 + t412;
t407 = Icges(4,1) * t395 - t426;
t406 = -Icges(4,2) * t392 + t425;
t405 = Icges(4,5) * t395 - Icges(4,6) * t392;
t339 = -Icges(4,6) * t389 + t388 * t406;
t341 = -Icges(4,5) * t389 + t388 * t407;
t404 = t339 * t392 - t341 * t395;
t340 = Icges(4,6) * t388 + t389 * t406;
t342 = Icges(4,5) * t388 + t389 * t407;
t403 = -t340 * t392 + t342 * t395;
t377 = Icges(4,2) * t395 + t426;
t378 = Icges(4,1) * t392 + t425;
t402 = -t377 * t392 + t378 * t395;
t382 = pkin(3) * t392 - pkin(8) * t395;
t401 = qJD(1) * t369 - t382 * t415 + t416;
t333 = pkin(4) * t355 + qJ(5) * t354;
t386 = -qJD(4) * t395 + qJD(1);
t400 = qJD(5) * t352 + t386 * t333 + t401;
t399 = (-t368 + t411) * qJD(1) - t382 * t414;
t371 = t388 * t413 - t414;
t373 = (pkin(4) * t394 + qJ(5) * t391) * t392;
t398 = qJD(5) * t354 + t371 * t373 + t399;
t381 = rSges(2,1) * t396 - rSges(2,2) * t393;
t380 = rSges(2,1) * t393 + rSges(2,2) * t396;
t379 = rSges(4,1) * t392 + rSges(4,2) * t395;
t376 = Icges(4,5) * t392 + Icges(4,6) * t395;
t367 = -rSges(6,1) * t395 + (-rSges(6,2) * t394 + rSges(6,3) * t391) * t392;
t365 = -rSges(5,3) * t395 + (rSges(5,1) * t394 - rSges(5,2) * t391) * t392;
t348 = t387 + qJD(1) * (rSges(3,1) * t389 - rSges(3,2) * t388);
t347 = (-rSges(3,1) * t388 - rSges(3,2) * t389 - t427) * qJD(1);
t344 = rSges(4,3) * t388 + t389 * t409;
t343 = -rSges(4,3) * t389 + t388 * t409;
t338 = Icges(4,3) * t388 + t389 * t405;
t337 = -Icges(4,3) * t389 + t388 * t405;
t330 = rSges(5,1) * t355 - rSges(5,2) * t354 + rSges(5,3) * t423;
t329 = rSges(5,1) * t353 - rSges(5,2) * t352 + rSges(5,3) * t424;
t328 = rSges(6,1) * t423 - rSges(6,2) * t355 + rSges(6,3) * t354;
t326 = rSges(6,1) * t424 - rSges(6,2) * t353 + rSges(6,3) * t352;
t305 = qJD(1) * t344 - t379 * t415 + t416;
t304 = -t379 * t414 + (-t343 + t411) * qJD(1);
t303 = qJD(2) + (t343 * t388 + t344 * t389) * qJD(3);
t302 = t330 * t386 - t365 * t370 + t401;
t301 = -t329 * t386 + t365 * t371 + t399;
t300 = t329 * t370 - t330 * t371 + t412;
t299 = t328 * t386 + (-t367 - t373) * t370 + t400;
t298 = t367 * t371 + (-t326 - t332) * t386 + t398;
t297 = t326 * t370 + (-t328 - t333) * t371 + t408;
t296 = qJD(6) * t353 + t418 * t386 + (-t373 - t417) * t370 + t400;
t295 = qJD(6) * t355 + t417 * t371 + (-t332 - t419) * t386 + t398;
t294 = qJD(6) * t421 + t419 * t370 + (-t333 - t418) * t371 + t408;
t1 = ((t388 * t376 + t389 * t402) * qJD(1) + (t388 ^ 2 * t338 + (t404 * t389 + (-t337 + t403) * t388) * t389) * qJD(3)) * t415 / 0.2e1 + m(6) * (t297 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + m(7) * (t294 ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(5) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 - ((-t389 * t376 + t388 * t402) * qJD(1) + (t389 ^ 2 * t337 + (t403 * t388 + (-t338 + t404) * t389) * t388) * qJD(3)) * t414 / 0.2e1 + qJD(1) * ((t395 * t377 + t392 * t378) * qJD(1) + ((t340 * t395 + t342 * t392) * t388 - (t339 * t395 + t341 * t392) * t389) * qJD(3)) / 0.2e1 + m(4) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + ((t433 * t354 + t432 * t355 + t434 * t423) * t386 + (t438 * t354 + t440 * t355 + t436 * t423) * t371 + (t437 * t354 + t439 * t355 + t435 * t423) * t370) * t370 / 0.2e1 + ((t433 * t352 + t432 * t353 + t434 * t424) * t386 + (t438 * t352 + t440 * t353 + t436 * t424) * t371 + (t437 * t352 + t439 * t353 + t435 * t424) * t370) * t371 / 0.2e1 + ((-t435 * t370 - t436 * t371 - t434 * t386) * t395 + ((t433 * t391 + t432 * t394) * t386 + (t438 * t391 + t440 * t394) * t371 + (t437 * t391 + t439 * t394) * t370) * t392) * t386 / 0.2e1 + (m(2) * (t380 ^ 2 + t381 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
