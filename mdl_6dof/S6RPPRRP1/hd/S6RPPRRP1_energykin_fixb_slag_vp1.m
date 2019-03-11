% Calculate kinetic energy for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:35
% EndTime: 2019-03-09 01:57:36
% DurationCPUTime: 1.38s
% Computational Cost: add. (1546->191), mult. (1312->293), div. (0->0), fcn. (1273->10), ass. (0->106)
t459 = Icges(6,1) + Icges(7,1);
t458 = Icges(6,4) + Icges(7,4);
t457 = -Icges(7,5) - Icges(6,5);
t456 = Icges(6,2) + Icges(7,2);
t455 = -Icges(7,6) - Icges(6,6);
t454 = -Icges(7,3) - Icges(6,3);
t392 = pkin(10) + qJ(4);
t390 = cos(t392);
t393 = qJ(1) + pkin(9);
t391 = cos(t393);
t400 = cos(qJ(5));
t429 = t391 * t400;
t389 = sin(t393);
t398 = sin(qJ(5));
t432 = t389 * t398;
t368 = -t390 * t432 - t429;
t430 = t391 * t398;
t431 = t389 * t400;
t369 = t390 * t431 - t430;
t388 = sin(t392);
t434 = t388 * t389;
t453 = -t368 * t455 - t369 * t457 - t434 * t454;
t370 = -t390 * t430 + t431;
t371 = t390 * t429 + t432;
t433 = t388 * t391;
t452 = -t370 * t455 - t371 * t457 - t433 * t454;
t451 = t368 * t456 + t369 * t458 - t434 * t455;
t450 = t370 * t456 + t371 * t458 - t433 * t455;
t449 = t458 * t368 + t369 * t459 - t457 * t434;
t448 = t458 * t370 + t371 * t459 - t457 * t433;
t447 = t454 * t390 + (t398 * t455 - t400 * t457) * t388;
t446 = t455 * t390 + (-t398 * t456 + t400 * t458) * t388;
t445 = t457 * t390 + (-t458 * t398 + t400 * t459) * t388;
t399 = sin(qJ(1));
t440 = t399 * pkin(1);
t395 = cos(pkin(10));
t439 = pkin(3) * t395;
t438 = pkin(5) * t400;
t436 = Icges(5,4) * t388;
t435 = Icges(5,4) * t390;
t405 = qJ(6) * t388 + t390 * t438;
t427 = rSges(7,1) * t369 + rSges(7,2) * t368 + rSges(7,3) * t434 - pkin(5) * t430 + t389 * t405;
t426 = rSges(7,1) * t371 + rSges(7,2) * t370 + rSges(7,3) * t433 + pkin(5) * t432 + t391 * t405;
t425 = (-qJ(6) - rSges(7,3)) * t390 + (rSges(7,1) * t400 - rSges(7,2) * t398 + t438) * t388;
t401 = cos(qJ(1));
t387 = qJD(1) * t401 * pkin(1);
t424 = qJD(1) * (pkin(2) * t391 + qJ(3) * t389) + t387;
t423 = qJD(4) * t389;
t422 = qJD(4) * t391;
t421 = qJD(5) * t388;
t417 = pkin(4) * t390 + pkin(8) * t388;
t364 = t417 * t389;
t365 = t417 * t391;
t420 = t364 * t423 + t365 * t422 + qJD(2);
t419 = -pkin(2) * t389 + qJ(3) * t391 - t440;
t418 = pkin(7) * t391 - t389 * t439 + t419;
t394 = sin(pkin(10));
t416 = rSges(4,1) * t395 - rSges(4,2) * t394;
t415 = rSges(5,1) * t390 - rSges(5,2) * t388;
t414 = Icges(5,1) * t390 - t436;
t413 = -Icges(5,2) * t388 + t435;
t412 = Icges(5,5) * t390 - Icges(5,6) * t388;
t347 = -Icges(5,6) * t391 + t389 * t413;
t349 = -Icges(5,5) * t391 + t389 * t414;
t411 = t347 * t388 - t349 * t390;
t348 = Icges(5,6) * t389 + t391 * t413;
t350 = Icges(5,5) * t389 + t391 * t414;
t410 = -t348 * t388 + t350 * t390;
t376 = Icges(5,2) * t390 + t436;
t377 = Icges(5,1) * t388 + t435;
t409 = -t376 * t388 + t377 * t390;
t380 = pkin(4) * t388 - pkin(8) * t390;
t408 = -qJD(4) * t380 + qJD(6) * t388;
t407 = -qJD(3) * t391 + qJD(1) * (pkin(7) * t389 + t391 * t439) + t424;
t406 = qJD(1) * t365 + t407;
t384 = qJD(3) * t389;
t404 = t384 + (-t364 + t418) * qJD(1);
t402 = qJD(2) ^ 2;
t383 = -qJD(5) * t390 + qJD(1);
t382 = rSges(2,1) * t401 - rSges(2,2) * t399;
t381 = rSges(2,1) * t399 + rSges(2,2) * t401;
t378 = rSges(5,1) * t388 + rSges(5,2) * t390;
t375 = Icges(5,5) * t388 + Icges(5,6) * t390;
t373 = t389 * t421 - t422;
t372 = t391 * t421 + t423;
t367 = t387 + qJD(1) * (rSges(3,1) * t391 - rSges(3,2) * t389);
t366 = (-rSges(3,1) * t389 - rSges(3,2) * t391 - t440) * qJD(1);
t362 = -rSges(6,3) * t390 + (rSges(6,1) * t400 - rSges(6,2) * t398) * t388;
t352 = rSges(5,3) * t389 + t391 * t415;
t351 = -rSges(5,3) * t391 + t389 * t415;
t346 = Icges(5,3) * t389 + t391 * t412;
t345 = -Icges(5,3) * t391 + t389 * t412;
t341 = qJD(1) * t389 * rSges(4,3) + (qJD(1) * t416 - qJD(3)) * t391 + t424;
t340 = t384 + (t391 * rSges(4,3) - t389 * t416 + t419) * qJD(1);
t339 = rSges(6,1) * t371 + rSges(6,2) * t370 + rSges(6,3) * t433;
t337 = rSges(6,1) * t369 + rSges(6,2) * t368 + rSges(6,3) * t434;
t321 = qJD(2) + (t351 * t389 + t352 * t391) * qJD(4);
t320 = qJD(1) * t352 - t378 * t423 + t407;
t319 = -t378 * t422 + t384 + (-t351 + t418) * qJD(1);
t318 = t337 * t372 - t339 * t373 + t420;
t317 = t339 * t383 - t362 * t372 - t380 * t423 + t406;
t316 = -t337 * t383 + t362 * t373 - t380 * t422 + t404;
t315 = -t372 * t425 + t383 * t426 + t389 * t408 + t406;
t314 = t373 * t425 - t383 * t427 + t391 * t408 + t404;
t313 = -qJD(6) * t390 + t372 * t427 - t373 * t426 + t420;
t1 = ((t389 * t375 + t391 * t409) * qJD(1) + (t389 ^ 2 * t346 + (t411 * t391 + (-t345 + t410) * t389) * t391) * qJD(4)) * t423 / 0.2e1 - ((-t391 * t375 + t389 * t409) * qJD(1) + (t391 ^ 2 * t345 + (t410 * t389 + (-t346 + t411) * t391) * t389) * qJD(4)) * t422 / 0.2e1 + qJD(1) * ((t390 * t376 + t388 * t377) * qJD(1) + ((t348 * t390 + t350 * t388) * t389 - (t347 * t390 + t388 * t349) * t391) * qJD(4)) / 0.2e1 + m(6) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(7) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(5) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(3) * (t366 ^ 2 + t367 ^ 2 + t402) / 0.2e1 + m(4) * (t340 ^ 2 + t341 ^ 2 + t402) / 0.2e1 + ((t446 * t370 + t445 * t371 + t447 * t433) * t383 + (t451 * t370 + t449 * t371 + t453 * t433) * t373 + (t450 * t370 + t448 * t371 + t452 * t433) * t372) * t372 / 0.2e1 + ((t446 * t368 + t445 * t369 + t447 * t434) * t383 + (t451 * t368 + t449 * t369 + t453 * t434) * t373 + (t450 * t368 + t448 * t369 + t452 * t434) * t372) * t373 / 0.2e1 + ((-t452 * t372 - t453 * t373 - t447 * t383) * t390 + ((-t446 * t398 + t445 * t400) * t383 + (-t451 * t398 + t449 * t400) * t373 + (-t450 * t398 + t448 * t400) * t372) * t388) * t383 / 0.2e1 + (m(2) * (t381 ^ 2 + t382 ^ 2) + Icges(2,3) + Icges(4,2) * t395 ^ 2 + (Icges(4,1) * t394 + 0.2e1 * Icges(4,4) * t395) * t394 + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
