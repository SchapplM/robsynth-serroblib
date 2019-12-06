% Calculate kinetic energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:28
% EndTime: 2019-12-05 18:05:30
% DurationCPUTime: 1.55s
% Computational Cost: add. (1026->191), mult. (1629->307), div. (0->0), fcn. (1667->8), ass. (0->101)
t424 = Icges(5,1) + Icges(6,1);
t423 = Icges(5,4) + Icges(6,4);
t422 = -Icges(6,5) - Icges(5,5);
t421 = Icges(5,2) + Icges(6,2);
t420 = -Icges(6,6) - Icges(5,6);
t419 = -Icges(6,3) - Icges(5,3);
t362 = sin(pkin(8));
t363 = cos(pkin(8));
t366 = cos(qJ(3));
t367 = cos(qJ(1));
t390 = t367 * t366;
t364 = sin(qJ(3));
t365 = sin(qJ(1));
t396 = t365 * t364;
t341 = t363 * t396 + t390;
t391 = t367 * t364;
t395 = t365 * t366;
t342 = -t363 * t395 + t391;
t343 = -t363 * t391 + t395;
t344 = t363 * t390 + t396;
t392 = t367 * t362;
t397 = t365 * t362;
t373 = (Icges(4,5) * t342 + Icges(4,6) * t341 - Icges(4,3) * t397) * t365 - (Icges(4,5) * t344 + Icges(4,6) * t343 + Icges(4,3) * t392) * t367;
t418 = t362 * t373;
t361 = qJ(3) + qJ(4);
t358 = cos(t361);
t393 = t367 * t358;
t357 = sin(t361);
t399 = t365 * t357;
t334 = t363 * t399 + t393;
t394 = t367 * t357;
t398 = t365 * t358;
t335 = -t363 * t398 + t394;
t417 = -t420 * t334 - t422 * t335 + t419 * t397;
t336 = -t363 * t394 + t398;
t337 = t363 * t393 + t399;
t416 = -t420 * t336 - t422 * t337 - t419 * t392;
t415 = t421 * t334 + t423 * t335 + t420 * t397;
t414 = t421 * t336 + t423 * t337 - t420 * t392;
t413 = t423 * t334 + t424 * t335 + t422 * t397;
t412 = t423 * t336 + t424 * t337 - t422 * t392;
t411 = t419 * t363 + (t420 * t357 - t422 * t358) * t362;
t410 = t420 * t363 + (-t421 * t357 + t423 * t358) * t362;
t409 = t422 * t363 + (-t423 * t357 + t424 * t358) * t362;
t386 = pkin(4) * t358;
t408 = qJ(5) * t362 + t363 * t386;
t401 = t366 * pkin(3);
t407 = pkin(7) * t362 + t363 * t401;
t377 = pkin(4) * t357;
t389 = t335 * rSges(6,1) + t334 * rSges(6,2) - rSges(6,3) * t397 - t408 * t365 + t377 * t367;
t388 = t337 * rSges(6,1) + t336 * rSges(6,2) + rSges(6,3) * t392 + t377 * t365 + t408 * t367;
t387 = (-qJ(5) - rSges(6,3)) * t363 + (rSges(6,1) * t358 - rSges(6,2) * t357 + t386) * t362;
t385 = qJD(1) * (-t365 * pkin(1) + t367 * qJ(2)) + qJD(2) * t365;
t383 = qJD(3) * t362;
t382 = qJD(5) * t362;
t381 = qJD(3) + qJD(4);
t376 = pkin(2) * t363 + pkin(6) * t362;
t380 = -qJD(1) * t376 * t365 + t385;
t379 = t365 * t383;
t378 = t367 * t383;
t375 = -rSges(3,1) * t363 + rSges(3,2) * t362;
t312 = pkin(3) * t391 - t407 * t365;
t321 = -pkin(7) * t363 + t362 * t401;
t353 = -qJD(3) * t363 + qJD(1);
t374 = t353 * t312 + t321 * t379 + t380;
t351 = t367 * pkin(1) + t365 * qJ(2);
t356 = qJD(2) * t367;
t372 = t356 + (-t376 * t367 - t351) * qJD(1);
t313 = pkin(3) * t396 + t407 * t367;
t371 = (-t312 * t367 - t313 * t365) * t383;
t370 = -t353 * t313 + t321 * t378 + t372;
t352 = t367 * rSges(2,1) - t365 * rSges(2,2);
t350 = -t365 * rSges(2,1) - t367 * rSges(2,2);
t346 = -t363 * t381 + qJD(1);
t340 = t381 * t392;
t339 = t381 * t397;
t333 = -t363 * rSges(4,3) + (rSges(4,1) * t366 - rSges(4,2) * t364) * t362;
t332 = -Icges(4,5) * t363 + (Icges(4,1) * t366 - Icges(4,4) * t364) * t362;
t331 = -Icges(4,6) * t363 + (Icges(4,4) * t366 - Icges(4,2) * t364) * t362;
t330 = -Icges(4,3) * t363 + (Icges(4,5) * t366 - Icges(4,6) * t364) * t362;
t329 = -t363 * rSges(5,3) + (rSges(5,1) * t358 - rSges(5,2) * t357) * t362;
t318 = t356 + (-t365 * rSges(3,3) + t367 * t375 - t351) * qJD(1);
t317 = qJD(1) * (t367 * rSges(3,3) + t365 * t375) + t385;
t315 = t344 * rSges(4,1) + t343 * rSges(4,2) + rSges(4,3) * t392;
t314 = t342 * rSges(4,1) + t341 * rSges(4,2) - rSges(4,3) * t397;
t311 = Icges(4,1) * t344 + Icges(4,4) * t343 + Icges(4,5) * t392;
t310 = Icges(4,1) * t342 + Icges(4,4) * t341 - Icges(4,5) * t397;
t309 = Icges(4,4) * t344 + Icges(4,2) * t343 + Icges(4,6) * t392;
t308 = Icges(4,4) * t342 + Icges(4,2) * t341 - Icges(4,6) * t397;
t305 = t337 * rSges(5,1) + t336 * rSges(5,2) + rSges(5,3) * t392;
t303 = t335 * rSges(5,1) + t334 * rSges(5,2) - rSges(5,3) * t397;
t286 = (-t314 * t367 - t315 * t365) * t383;
t285 = -t353 * t315 + t333 * t378 + t372;
t284 = t353 * t314 + t333 * t379 + t380;
t283 = -t346 * t305 + t340 * t329 + t370;
t282 = t346 * t303 + t339 * t329 + t374;
t281 = -t340 * t303 - t339 * t305 + t371;
t280 = t340 * t387 - t346 * t388 - t365 * t382 + t370;
t279 = t339 * t387 + t346 * t389 + t367 * t382 + t374;
t278 = -qJD(5) * t363 - t339 * t388 - t340 * t389 + t371;
t1 = m(3) * (t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(4) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + t353 * ((-t363 * t330 + (-t331 * t364 + t332 * t366) * t362) * t353 + ((-(-t308 * t364 + t310 * t366) * t365 + (-t309 * t364 + t311 * t366) * t367) * t362 + t373 * t363) * t383) / 0.2e1 - ((-t330 * t397 + t341 * t331 + t342 * t332) * t353 + ((t341 * t309 + t342 * t311) * t367 + (-t341 * t308 - t342 * t310 + t418) * t365) * t383) * t379 / 0.2e1 + ((t330 * t392 + t343 * t331 + t344 * t332) * t353 + (-(t343 * t308 + t344 * t310) * t365 + (t343 * t309 + t344 * t311 - t418) * t367) * t383) * t378 / 0.2e1 + m(5) * (t281 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(6) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 - ((t410 * t334 + t409 * t335 - t411 * t397) * t346 + (t414 * t334 + t412 * t335 - t416 * t397) * t340 - (t415 * t334 + t413 * t335 - t417 * t397) * t339) * t339 / 0.2e1 + ((t410 * t336 + t409 * t337 + t411 * t392) * t346 + (t414 * t336 + t412 * t337 + t416 * t392) * t340 - (t415 * t336 + t413 * t337 + t417 * t392) * t339) * t340 / 0.2e1 + ((t417 * t339 - t416 * t340 - t411 * t346) * t363 + ((-t410 * t357 + t409 * t358) * t346 + (-t414 * t357 + t412 * t358) * t340 - (-t415 * t357 + t413 * t358) * t339) * t362) * t346 / 0.2e1 + (m(2) * (t350 ^ 2 + t352 ^ 2) + Icges(2,3) + Icges(3,2) * t363 ^ 2 + (Icges(3,1) * t362 + 0.2e1 * Icges(3,4) * t363) * t362) * qJD(1) ^ 2 / 0.2e1;
T = t1;
