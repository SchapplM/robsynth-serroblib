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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:49:13
% EndTime: 2020-01-03 11:49:14
% DurationCPUTime: 1.73s
% Computational Cost: add. (1026->196), mult. (1629->314), div. (0->0), fcn. (1667->8), ass. (0->103)
t421 = Icges(5,1) + Icges(6,1);
t420 = Icges(5,4) + Icges(6,4);
t419 = -Icges(6,5) - Icges(5,5);
t418 = Icges(5,2) + Icges(6,2);
t417 = -Icges(6,6) - Icges(5,6);
t416 = -Icges(6,3) - Icges(5,3);
t356 = sin(pkin(8));
t357 = cos(pkin(8));
t360 = cos(qJ(3));
t361 = cos(qJ(1));
t387 = t361 * t360;
t358 = sin(qJ(3));
t359 = sin(qJ(1));
t393 = t359 * t358;
t337 = -t357 * t393 - t387;
t388 = t361 * t358;
t392 = t359 * t360;
t338 = t357 * t392 - t388;
t339 = t357 * t388 - t392;
t340 = -t357 * t387 - t393;
t389 = t361 * t356;
t394 = t359 * t356;
t365 = (Icges(4,5) * t338 + Icges(4,6) * t337 + Icges(4,3) * t394) * t359 - (Icges(4,5) * t340 + Icges(4,6) * t339 - Icges(4,3) * t389) * t361;
t415 = t356 * t365;
t355 = qJ(3) + qJ(4);
t352 = cos(t355);
t390 = t361 * t352;
t351 = sin(t355);
t396 = t359 * t351;
t330 = -t357 * t396 - t390;
t391 = t361 * t351;
t395 = t359 * t352;
t331 = t357 * t395 - t391;
t414 = -t417 * t330 - t419 * t331 - t416 * t394;
t332 = t357 * t391 - t395;
t333 = -t357 * t390 - t396;
t413 = -t417 * t332 - t419 * t333 + t416 * t389;
t412 = t418 * t330 + t420 * t331 - t417 * t394;
t411 = t418 * t332 + t420 * t333 + t417 * t389;
t410 = t420 * t330 + t421 * t331 - t419 * t394;
t409 = t420 * t332 + t421 * t333 + t419 * t389;
t408 = t416 * t357 + (t417 * t351 - t419 * t352) * t356;
t407 = t417 * t357 + (-t418 * t351 + t420 * t352) * t356;
t406 = t419 * t357 + (-t420 * t351 + t421 * t352) * t356;
t381 = pkin(4) * t352;
t405 = qJ(5) * t356 + t357 * t381;
t398 = t360 * pkin(3);
t404 = pkin(7) * t356 + t357 * t398;
t371 = pkin(4) * t351;
t386 = t331 * rSges(6,1) + t330 * rSges(6,2) + rSges(6,3) * t394 + t405 * t359 - t371 * t361;
t385 = t333 * rSges(6,1) + t332 * rSges(6,2) - rSges(6,3) * t389 - t371 * t359 - t405 * t361;
t310 = -pkin(3) * t388 + t404 * t359;
t311 = -pkin(3) * t393 - t404 * t361;
t376 = qJD(3) * t361;
t372 = t356 * t376;
t377 = qJD(3) * t359;
t373 = t356 * t377;
t384 = t310 * t372 + t311 * t373;
t383 = (-qJ(5) - rSges(6,3)) * t357 + (rSges(6,1) * t352 - rSges(6,2) * t351 + t381) * t356;
t344 = qJD(1) * (t359 * pkin(1) - t361 * qJ(2));
t369 = pkin(2) * t357 + pkin(6) * t356;
t382 = qJD(1) * t369 * t359 + t344;
t379 = qJD(2) * t361;
t378 = qJD(3) * t356;
t375 = qJD(3) + qJD(4);
t349 = -qJD(3) * t357 + qJD(1);
t374 = t349 * t310 + t382;
t347 = -t361 * pkin(1) - t359 * qJ(2);
t370 = (t369 * t361 - t347) * qJD(1);
t368 = rSges(3,1) * t357 - rSges(3,2) * t356;
t317 = -pkin(7) * t357 + t356 * t398;
t367 = -t317 * t378 - qJD(2);
t366 = -(-t357 * rSges(4,3) + (rSges(4,1) * t360 - rSges(4,2) * t358) * t356) * t378 - qJD(2);
t364 = -t349 * t311 + t370;
t348 = -t361 * rSges(2,1) + t359 * rSges(2,2);
t346 = t359 * rSges(2,1) + t361 * rSges(2,2);
t342 = -t357 * t375 + qJD(1);
t336 = t375 * t389;
t335 = t375 * t394;
t328 = -Icges(4,5) * t357 + (Icges(4,1) * t360 - Icges(4,4) * t358) * t356;
t327 = -Icges(4,6) * t357 + (Icges(4,4) * t360 - Icges(4,2) * t358) * t356;
t326 = -Icges(4,3) * t357 + (Icges(4,5) * t360 - Icges(4,6) * t358) * t356;
t325 = -t357 * rSges(5,3) + (rSges(5,1) * t352 - rSges(5,2) * t351) * t356;
t316 = -t379 + (t359 * rSges(3,3) + t361 * t368 - t347) * qJD(1);
t315 = -qJD(1) * t361 * rSges(3,3) + t344 + (qJD(1) * t368 - qJD(2)) * t359;
t313 = t340 * rSges(4,1) + t339 * rSges(4,2) - rSges(4,3) * t389;
t312 = t338 * rSges(4,1) + t337 * rSges(4,2) + rSges(4,3) * t394;
t309 = Icges(4,1) * t340 + Icges(4,4) * t339 - Icges(4,5) * t389;
t308 = Icges(4,1) * t338 + Icges(4,4) * t337 + Icges(4,5) * t394;
t307 = Icges(4,4) * t340 + Icges(4,2) * t339 - Icges(4,6) * t389;
t306 = Icges(4,4) * t338 + Icges(4,2) * t337 + Icges(4,6) * t394;
t301 = t333 * rSges(5,1) + t332 * rSges(5,2) - rSges(5,3) * t389;
t299 = t331 * rSges(5,1) + t330 * rSges(5,2) + rSges(5,3) * t394;
t282 = (t312 * t361 + t313 * t359) * t378;
t281 = -t349 * t313 + t361 * t366 + t370;
t280 = t349 * t312 + t359 * t366 + t382;
t279 = -t342 * t301 - t336 * t325 + t361 * t367 + t364;
t278 = t342 * t299 - t335 * t325 + t359 * t367 + t374;
t277 = t336 * t299 + t335 * t301 + t384;
t276 = -t379 + (qJD(5) * t359 - t317 * t376) * t356 - t385 * t342 - t383 * t336 + t364;
t275 = -qJD(2) * t359 + (-qJD(5) * t361 - t317 * t377) * t356 + t386 * t342 - t383 * t335 + t374;
t274 = -qJD(5) * t357 + t335 * t385 + t336 * t386 + t384;
t1 = m(3) * (t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(4) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + t349 * ((-t357 * t326 + (-t327 * t358 + t328 * t360) * t356) * t349 + (((-t306 * t358 + t308 * t360) * t359 - (-t307 * t358 + t309 * t360) * t361) * t356 - t365 * t357) * t378) / 0.2e1 + ((t326 * t394 + t337 * t327 + t338 * t328) * t349 + (-(t337 * t307 + t338 * t309) * t361 + (t337 * t306 + t338 * t308 + t415) * t359) * t378) * t373 / 0.2e1 - ((-t326 * t389 + t339 * t327 + t340 * t328) * t349 + ((t339 * t306 + t340 * t308) * t359 + (-t339 * t307 - t340 * t309 - t415) * t361) * t378) * t372 / 0.2e1 + m(5) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(6) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t407 * t330 + t406 * t331 + t408 * t394) * t342 - (t411 * t330 + t409 * t331 + t413 * t394) * t336 + (t412 * t330 + t410 * t331 + t414 * t394) * t335) * t335 / 0.2e1 - ((t407 * t332 + t406 * t333 - t408 * t389) * t342 - (t411 * t332 + t409 * t333 - t413 * t389) * t336 + (t412 * t332 + t410 * t333 - t414 * t389) * t335) * t336 / 0.2e1 + ((-t414 * t335 + t413 * t336 - t408 * t342) * t357 + ((-t407 * t351 + t406 * t352) * t342 - (-t411 * t351 + t409 * t352) * t336 + (-t412 * t351 + t410 * t352) * t335) * t356) * t342 / 0.2e1 + (m(2) * (t346 ^ 2 + t348 ^ 2) + Icges(2,3) + Icges(3,2) * t357 ^ 2 + (Icges(3,1) * t356 + 0.2e1 * Icges(3,4) * t357) * t356) * qJD(1) ^ 2 / 0.2e1;
T = t1;
