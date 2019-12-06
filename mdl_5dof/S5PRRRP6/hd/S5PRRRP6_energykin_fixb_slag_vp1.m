% Calculate kinetic energy for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:44
% EndTime: 2019-12-05 16:50:46
% DurationCPUTime: 1.98s
% Computational Cost: add. (1039->186), mult. (1730->304), div. (0->0), fcn. (1776->8), ass. (0->103)
t431 = Icges(5,1) + Icges(6,1);
t430 = -Icges(5,4) + Icges(6,5);
t429 = Icges(6,4) + Icges(5,5);
t428 = Icges(5,2) + Icges(6,3);
t427 = -Icges(6,6) + Icges(5,6);
t426 = -Icges(5,3) - Icges(6,2);
t368 = cos(pkin(8));
t410 = t368 ^ 2;
t367 = sin(pkin(8));
t411 = t367 ^ 2;
t412 = t410 + t411;
t425 = qJD(2) * t412;
t424 = rSges(6,1) + pkin(4);
t423 = rSges(6,3) + qJ(5);
t422 = t367 * t368;
t366 = qJ(3) + qJ(4);
t364 = sin(t366);
t365 = cos(t366);
t372 = cos(qJ(2));
t401 = t367 * t372;
t340 = t364 * t401 + t365 * t368;
t341 = -t364 * t368 + t365 * t401;
t370 = sin(qJ(2));
t402 = t367 * t370;
t421 = t428 * t340 + t430 * t341 - t427 * t402;
t398 = t368 * t372;
t342 = t364 * t398 - t367 * t365;
t343 = t364 * t367 + t365 * t398;
t399 = t368 * t370;
t420 = t428 * t342 + t430 * t343 - t427 * t399;
t419 = -t427 * t340 + t429 * t341 - t426 * t402;
t418 = -t427 * t342 + t429 * t343 - t426 * t399;
t417 = t430 * t340 + t431 * t341 + t429 * t402;
t416 = t430 * t342 + t431 * t343 + t429 * t399;
t415 = t427 * t372 + (t428 * t364 + t430 * t365) * t370;
t414 = t426 * t372 + (-t427 * t364 + t429 * t365) * t370;
t413 = -t429 * t372 + (t430 * t364 + t431 * t365) * t370;
t409 = qJD(2) ^ 2;
t371 = cos(qJ(3));
t405 = pkin(3) * t371;
t369 = sin(qJ(3));
t403 = t367 * t369;
t400 = t368 * t369;
t397 = t369 * t372;
t396 = t371 * t372;
t395 = rSges(6,2) * t402 + t423 * t340 + t424 * t341;
t394 = rSges(6,2) * t399 + t423 * t342 + t424 * t343;
t393 = -rSges(6,2) * t372 + (t423 * t364 + t424 * t365) * t370;
t363 = qJD(2) * t367;
t391 = qJD(3) * t370;
t353 = t368 * t391 + t363;
t392 = qJD(2) * t368;
t390 = qJD(3) * t372;
t389 = qJD(4) * t370;
t358 = pkin(2) * t370 - pkin(6) * t372;
t388 = t358 * t363;
t387 = t358 * t392;
t386 = qJD(1) + (pkin(2) * t372 + pkin(6) * t370) * t425;
t354 = t367 * t391 - t392;
t383 = Icges(3,5) * t372 - Icges(3,6) * t370;
t376 = pkin(7) * t370 + t372 * t405;
t318 = -pkin(3) * t400 + t367 * t376;
t321 = -pkin(7) * t372 + t370 * t405;
t380 = t318 * t390 + t354 * t321 - t387;
t319 = pkin(3) * t403 + t368 * t376;
t377 = t353 * t318 - t319 * t354 + t386;
t375 = -t319 * t390 - t321 * t353 - t388;
t357 = rSges(3,1) * t370 + rSges(3,2) * t372;
t355 = (-qJD(3) - qJD(4)) * t372;
t352 = t368 * t396 + t403;
t351 = t367 * t371 - t368 * t397;
t350 = t367 * t396 - t400;
t349 = -t367 * t397 - t368 * t371;
t347 = -rSges(4,3) * t372 + (rSges(4,1) * t371 - rSges(4,2) * t369) * t370;
t346 = -Icges(4,5) * t372 + (Icges(4,1) * t371 - Icges(4,4) * t369) * t370;
t345 = -Icges(4,6) * t372 + (Icges(4,4) * t371 - Icges(4,2) * t369) * t370;
t344 = -Icges(4,3) * t372 + (Icges(4,5) * t371 - Icges(4,6) * t369) * t370;
t333 = Icges(3,3) * t367 + t368 * t383;
t332 = -Icges(3,3) * t368 + t367 * t383;
t331 = t367 * t389 + t354;
t330 = t368 * t389 + t353;
t329 = -rSges(5,3) * t372 + (rSges(5,1) * t365 - rSges(5,2) * t364) * t370;
t317 = rSges(4,1) * t352 + rSges(4,2) * t351 + rSges(4,3) * t399;
t316 = rSges(4,1) * t350 + rSges(4,2) * t349 + rSges(4,3) * t402;
t313 = Icges(4,1) * t352 + Icges(4,4) * t351 + Icges(4,5) * t399;
t312 = Icges(4,1) * t350 + Icges(4,4) * t349 + Icges(4,5) * t402;
t311 = Icges(4,4) * t352 + Icges(4,2) * t351 + Icges(4,6) * t399;
t310 = Icges(4,4) * t350 + Icges(4,2) * t349 + Icges(4,6) * t402;
t309 = Icges(4,5) * t352 + Icges(4,6) * t351 + Icges(4,3) * t399;
t308 = Icges(4,5) * t350 + Icges(4,6) * t349 + Icges(4,3) * t402;
t306 = rSges(5,1) * t343 - rSges(5,2) * t342 + rSges(5,3) * t399;
t304 = rSges(5,1) * t341 - rSges(5,2) * t340 + rSges(5,3) * t402;
t290 = qJD(1) + (rSges(3,1) * t372 - rSges(3,2) * t370) * t425;
t288 = t316 * t390 + t347 * t354 - t387;
t287 = -t317 * t390 - t347 * t353 - t388;
t286 = t316 * t353 - t317 * t354 + t386;
t285 = -t304 * t355 + t329 * t331 + t380;
t284 = t306 * t355 - t329 * t330 + t375;
t283 = t304 * t330 - t306 * t331 + t377;
t282 = qJD(5) * t342 + t331 * t393 - t355 * t395 + t380;
t281 = qJD(5) * t340 - t330 * t393 + t355 * t394 + t375;
t280 = qJD(5) * t364 * t370 + t330 * t395 - t331 * t394 + t377;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t412 * t409 * t357 ^ 2 + t290 ^ 2) / 0.2e1 + t409 * t367 * (-t332 * t422 + t411 * t333) / 0.2e1 - t409 * t368 * (t410 * t332 - t333 * t422) / 0.2e1 + m(4) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + t353 * ((t309 * t399 + t351 * t311 + t352 * t313) * t353 + (t308 * t399 + t310 * t351 + t312 * t352) * t354 - (t344 * t399 + t345 * t351 + t346 * t352) * t390) / 0.2e1 + t354 * ((t309 * t402 + t311 * t349 + t313 * t350) * t353 + (t308 * t402 + t349 * t310 + t350 * t312) * t354 - (t344 * t402 + t345 * t349 + t346 * t350) * t390) / 0.2e1 - ((-t308 * t354 - t309 * t353 + t344 * t390) * t372 + ((-t311 * t369 + t313 * t371) * t353 + (-t310 * t369 + t312 * t371) * t354 - (-t345 * t369 + t346 * t371) * t390) * t370) * t390 / 0.2e1 + m(5) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(6) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + ((t415 * t342 + t413 * t343 + t414 * t399) * t355 + (t421 * t342 + t417 * t343 + t419 * t399) * t331 + (t420 * t342 + t416 * t343 + t418 * t399) * t330) * t330 / 0.2e1 + ((t415 * t340 + t413 * t341 + t414 * t402) * t355 + (t421 * t340 + t417 * t341 + t419 * t402) * t331 + (t420 * t340 + t416 * t341 + t418 * t402) * t330) * t331 / 0.2e1 + ((-t418 * t330 - t419 * t331 - t414 * t355) * t372 + ((t415 * t364 + t413 * t365) * t355 + (t421 * t364 + t417 * t365) * t331 + (t420 * t364 + t416 * t365) * t330) * t370) * t355 / 0.2e1;
T = t1;
