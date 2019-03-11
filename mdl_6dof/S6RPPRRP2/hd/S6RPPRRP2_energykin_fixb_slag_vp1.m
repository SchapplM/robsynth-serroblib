% Calculate kinetic energy for
% S6RPPRRP2
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:04
% EndTime: 2019-03-09 02:00:05
% DurationCPUTime: 1.37s
% Computational Cost: add. (1514->185), mult. (1303->285), div. (0->0), fcn. (1274->10), ass. (0->105)
t455 = Icges(6,1) + Icges(7,1);
t454 = -Icges(6,4) + Icges(7,5);
t453 = Icges(7,4) + Icges(6,5);
t452 = Icges(6,2) + Icges(7,3);
t451 = -Icges(7,6) + Icges(6,6);
t450 = -Icges(6,3) - Icges(7,2);
t449 = rSges(7,1) + pkin(5);
t448 = rSges(7,3) + qJ(6);
t391 = pkin(10) + qJ(4);
t389 = cos(t391);
t392 = qJ(1) + pkin(9);
t390 = cos(t392);
t398 = cos(qJ(5));
t425 = t390 * t398;
t388 = sin(t392);
t396 = sin(qJ(5));
t428 = t388 * t396;
t366 = t389 * t428 + t425;
t426 = t390 * t396;
t427 = t388 * t398;
t367 = t389 * t427 - t426;
t387 = sin(t391);
t430 = t387 * t388;
t447 = t452 * t366 + t454 * t367 - t451 * t430;
t368 = t389 * t426 - t427;
t369 = t389 * t425 + t428;
t429 = t387 * t390;
t446 = t452 * t368 + t454 * t369 - t451 * t429;
t445 = -t451 * t366 + t453 * t367 - t450 * t430;
t444 = -t451 * t368 + t453 * t369 - t450 * t429;
t443 = t454 * t366 + t455 * t367 + t453 * t430;
t442 = t454 * t368 + t455 * t369 + t453 * t429;
t441 = t451 * t389 + (t452 * t396 + t454 * t398) * t387;
t440 = t450 * t389 + (-t451 * t396 + t453 * t398) * t387;
t439 = -t453 * t389 + (t454 * t396 + t455 * t398) * t387;
t397 = sin(qJ(1));
t434 = t397 * pkin(1);
t394 = cos(pkin(10));
t433 = pkin(3) * t394;
t432 = Icges(5,4) * t387;
t431 = Icges(5,4) * t389;
t423 = rSges(7,2) * t430 + t448 * t366 + t367 * t449;
t422 = rSges(7,2) * t429 + t448 * t368 + t369 * t449;
t421 = -rSges(7,2) * t389 + (t448 * t396 + t398 * t449) * t387;
t399 = cos(qJ(1));
t386 = qJD(1) * t399 * pkin(1);
t420 = qJD(1) * (pkin(2) * t390 + qJ(3) * t388) + t386;
t419 = qJD(4) * t388;
t418 = qJD(4) * t390;
t417 = qJD(5) * t387;
t413 = pkin(4) * t389 + pkin(8) * t387;
t362 = t413 * t388;
t363 = t413 * t390;
t416 = t362 * t419 + t363 * t418 + qJD(2);
t415 = -pkin(2) * t388 + qJ(3) * t390 - t434;
t414 = pkin(7) * t390 - t388 * t433 + t415;
t393 = sin(pkin(10));
t412 = rSges(4,1) * t394 - rSges(4,2) * t393;
t411 = rSges(5,1) * t389 - rSges(5,2) * t387;
t410 = Icges(5,1) * t389 - t432;
t409 = -Icges(5,2) * t387 + t431;
t408 = Icges(5,5) * t389 - Icges(5,6) * t387;
t345 = -Icges(5,6) * t390 + t388 * t409;
t347 = -Icges(5,5) * t390 + t388 * t410;
t407 = t345 * t387 - t347 * t389;
t346 = Icges(5,6) * t388 + t390 * t409;
t348 = Icges(5,5) * t388 + t390 * t410;
t406 = -t346 * t387 + t348 * t389;
t375 = Icges(5,2) * t389 + t432;
t376 = Icges(5,1) * t387 + t431;
t405 = -t375 * t387 + t376 * t389;
t404 = -qJD(3) * t390 + qJD(1) * (pkin(7) * t388 + t390 * t433) + t420;
t379 = pkin(4) * t387 - pkin(8) * t389;
t403 = qJD(1) * t363 - t379 * t419 + t404;
t384 = qJD(3) * t388;
t402 = t384 + (-t362 + t414) * qJD(1) - t379 * t418;
t400 = qJD(2) ^ 2;
t382 = -qJD(5) * t389 + qJD(1);
t381 = rSges(2,1) * t399 - rSges(2,2) * t397;
t380 = rSges(2,1) * t397 + rSges(2,2) * t399;
t377 = rSges(5,1) * t387 + rSges(5,2) * t389;
t374 = Icges(5,5) * t387 + Icges(5,6) * t389;
t371 = t388 * t417 - t418;
t370 = t390 * t417 + t419;
t365 = t386 + qJD(1) * (rSges(3,1) * t390 - rSges(3,2) * t388);
t364 = (-rSges(3,1) * t388 - rSges(3,2) * t390 - t434) * qJD(1);
t360 = -rSges(6,3) * t389 + (rSges(6,1) * t398 - rSges(6,2) * t396) * t387;
t350 = rSges(5,3) * t388 + t390 * t411;
t349 = -rSges(5,3) * t390 + t388 * t411;
t344 = Icges(5,3) * t388 + t390 * t408;
t343 = -Icges(5,3) * t390 + t388 * t408;
t338 = qJD(1) * t388 * rSges(4,3) + (qJD(1) * t412 - qJD(3)) * t390 + t420;
t337 = t384 + (t390 * rSges(4,3) - t388 * t412 + t415) * qJD(1);
t336 = rSges(6,1) * t369 - rSges(6,2) * t368 + rSges(6,3) * t429;
t334 = rSges(6,1) * t367 - rSges(6,2) * t366 + rSges(6,3) * t430;
t320 = qJD(2) + (t349 * t388 + t350 * t390) * qJD(4);
t319 = qJD(1) * t350 - t377 * t419 + t404;
t318 = -t377 * t418 + t384 + (-t349 + t414) * qJD(1);
t317 = t334 * t370 - t336 * t371 + t416;
t316 = t336 * t382 - t360 * t370 + t403;
t315 = -t334 * t382 + t360 * t371 + t402;
t314 = qJD(6) * t366 - t370 * t421 + t382 * t422 + t403;
t313 = qJD(6) * t368 + t371 * t421 - t382 * t423 + t402;
t312 = qJD(6) * t387 * t396 + t370 * t423 - t371 * t422 + t416;
t1 = ((t388 * t374 + t390 * t405) * qJD(1) + (t388 ^ 2 * t344 + (t407 * t390 + (-t343 + t406) * t388) * t390) * qJD(4)) * t419 / 0.2e1 + m(6) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(5) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(3) * (t364 ^ 2 + t365 ^ 2 + t400) / 0.2e1 + m(4) * (t337 ^ 2 + t338 ^ 2 + t400) / 0.2e1 + m(7) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + qJD(1) * ((t389 * t375 + t387 * t376) * qJD(1) + ((t346 * t389 + t348 * t387) * t388 - (t345 * t389 + t347 * t387) * t390) * qJD(4)) / 0.2e1 - ((-t390 * t374 + t388 * t405) * qJD(1) + (t390 ^ 2 * t343 + (t406 * t388 + (-t344 + t407) * t390) * t388) * qJD(4)) * t418 / 0.2e1 + ((t441 * t368 + t439 * t369 + t440 * t429) * t382 + (t447 * t368 + t443 * t369 + t445 * t429) * t371 + (t446 * t368 + t442 * t369 + t444 * t429) * t370) * t370 / 0.2e1 + ((t441 * t366 + t439 * t367 + t440 * t430) * t382 + (t447 * t366 + t443 * t367 + t445 * t430) * t371 + (t446 * t366 + t442 * t367 + t444 * t430) * t370) * t371 / 0.2e1 + ((-t444 * t370 - t445 * t371 - t440 * t382) * t389 + ((t441 * t396 + t439 * t398) * t382 + (t447 * t396 + t443 * t398) * t371 + (t446 * t396 + t442 * t398) * t370) * t387) * t382 / 0.2e1 + (m(2) * (t380 ^ 2 + t381 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t394 ^ 2 + (Icges(4,1) * t393 + 0.2e1 * Icges(4,4) * t394) * t393) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
