% Calculate kinetic energy for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP13_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP13_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:30
% EndTime: 2019-12-31 18:58:31
% DurationCPUTime: 1.28s
% Computational Cost: add. (520->166), mult. (1232->262), div. (0->0), fcn. (1224->6), ass. (0->92)
t383 = Icges(5,1) + Icges(6,1);
t382 = Icges(5,4) - Icges(6,5);
t381 = Icges(6,4) + Icges(5,5);
t380 = Icges(5,2) + Icges(6,3);
t379 = Icges(6,6) - Icges(5,6);
t378 = Icges(5,3) + Icges(6,2);
t377 = rSges(6,1) + pkin(4);
t376 = rSges(6,3) + qJ(5);
t330 = sin(qJ(3));
t332 = cos(qJ(4));
t334 = cos(qJ(1));
t355 = t334 * t332;
t329 = sin(qJ(4));
t331 = sin(qJ(1));
t360 = t331 * t329;
t306 = t330 * t360 - t355;
t356 = t334 * t329;
t359 = t331 * t332;
t307 = t330 * t359 + t356;
t333 = cos(qJ(3));
t358 = t331 * t333;
t375 = t380 * t306 - t382 * t307 - t379 * t358;
t308 = t330 * t356 + t359;
t309 = -t330 * t355 + t360;
t357 = t333 * t334;
t374 = -t380 * t308 - t382 * t309 + t379 * t357;
t373 = t379 * t306 + t381 * t307 - t378 * t358;
t372 = -t379 * t308 + t381 * t309 + t378 * t357;
t371 = -t382 * t306 + t383 * t307 - t381 * t358;
t370 = t382 * t308 + t383 * t309 + t381 * t357;
t369 = (t380 * t329 - t382 * t332) * t333 + t379 * t330;
t368 = (t379 * t329 + t381 * t332) * t333 + t378 * t330;
t367 = (-t382 * t329 + t383 * t332) * t333 + t381 * t330;
t362 = Icges(4,4) * t330;
t361 = Icges(4,4) * t333;
t354 = -rSges(6,2) * t358 + t376 * t306 + t377 * t307;
t353 = rSges(6,2) * t357 - t376 * t308 + t377 * t309;
t352 = t330 * rSges(6,2) + (t376 * t329 + t377 * t332) * t333;
t316 = qJD(1) * (t334 * pkin(1) + t331 * qJ(2));
t351 = qJD(1) * t334 * pkin(6) + t316;
t350 = qJD(3) * t331;
t349 = qJD(3) * t334;
t348 = qJD(4) * t333;
t320 = t331 * pkin(1) - t334 * qJ(2);
t347 = -pkin(6) * t331 - t320;
t346 = pkin(3) * t330 - pkin(7) * t333;
t311 = t346 * t331;
t312 = t346 * t334;
t345 = -t311 * t350 - t312 * t349;
t344 = rSges(4,1) * t330 + rSges(4,2) * t333;
t343 = Icges(4,1) * t330 + t361;
t342 = Icges(4,2) * t333 + t362;
t341 = Icges(4,5) * t330 + Icges(4,6) * t333;
t294 = Icges(4,6) * t334 + t342 * t331;
t298 = Icges(4,5) * t334 + t343 * t331;
t340 = -t294 * t333 - t298 * t330;
t295 = Icges(4,6) * t331 - t342 * t334;
t299 = Icges(4,5) * t331 - t343 * t334;
t339 = t295 * t333 + t299 * t330;
t318 = -Icges(4,2) * t330 + t361;
t319 = Icges(4,1) * t333 - t362;
t338 = t318 * t333 + t319 * t330;
t324 = t333 * pkin(3) + t330 * pkin(7);
t328 = qJD(2) * t331;
t337 = t324 * t350 + t328 + (t312 + t347) * qJD(1);
t336 = qJD(1) * t311 + (-qJD(3) * t324 - qJD(2)) * t334 + t351;
t325 = qJD(4) * t330 + qJD(1);
t323 = t334 * rSges(2,1) - t331 * rSges(2,2);
t322 = t333 * rSges(4,1) - t330 * rSges(4,2);
t321 = t331 * rSges(2,1) + t334 * rSges(2,2);
t317 = Icges(4,5) * t333 - Icges(4,6) * t330;
t315 = -t331 * t348 + t349;
t314 = t334 * t348 + t350;
t303 = t331 * rSges(4,3) - t344 * t334;
t302 = t330 * rSges(5,3) + (rSges(5,1) * t332 - rSges(5,2) * t329) * t333;
t300 = t334 * rSges(4,3) + t344 * t331;
t291 = Icges(4,3) * t331 - t341 * t334;
t290 = Icges(4,3) * t334 + t341 * t331;
t287 = t316 - qJD(2) * t334 + qJD(1) * (-t334 * rSges(3,2) + t331 * rSges(3,3));
t286 = t328 + (t331 * rSges(3,2) + t334 * rSges(3,3) - t320) * qJD(1);
t283 = t309 * rSges(5,1) + t308 * rSges(5,2) + rSges(5,3) * t357;
t281 = t307 * rSges(5,1) - t306 * rSges(5,2) - rSges(5,3) * t358;
t267 = (-t300 * t331 + t303 * t334) * qJD(3);
t266 = qJD(1) * t300 + (-qJD(3) * t322 - qJD(2)) * t334 + t351;
t265 = t322 * t350 + t328 + (-t303 + t347) * qJD(1);
t264 = t325 * t281 - t315 * t302 + t336;
t263 = -t325 * t283 + t314 * t302 + t337;
t262 = -t314 * t281 + t315 * t283 + t345;
t261 = -qJD(5) * t308 - t352 * t315 + t354 * t325 + t336;
t260 = qJD(5) * t306 + t352 * t314 - t353 * t325 + t337;
t259 = qJD(5) * t333 * t329 - t354 * t314 + t353 * t315 + t345;
t1 = m(3) * (t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(4) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + ((t334 * t317 + t338 * t331) * qJD(1) + (t334 ^ 2 * t290 + (t339 * t331 + (t291 - t340) * t334) * t331) * qJD(3)) * t349 / 0.2e1 + ((t331 * t317 - t338 * t334) * qJD(1) + (t331 ^ 2 * t291 + (t340 * t334 + (t290 - t339) * t331) * t334) * qJD(3)) * t350 / 0.2e1 + qJD(1) * ((-t330 * t318 + t333 * t319) * qJD(1) + ((-t330 * t294 + t333 * t298) * t334 + (-t330 * t295 + t333 * t299) * t331) * qJD(3)) / 0.2e1 + m(5) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + m(6) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + ((-t308 * t369 + t309 * t367 + t357 * t368) * t325 + (-t308 * t375 + t371 * t309 + t373 * t357) * t315 + (-t374 * t308 + t370 * t309 + t372 * t357) * t314) * t314 / 0.2e1 + ((t306 * t369 + t307 * t367 - t358 * t368) * t325 + (t375 * t306 + t371 * t307 - t373 * t358) * t315 + (t306 * t374 + t307 * t370 - t358 * t372) * t314) * t315 / 0.2e1 + (((t329 * t369 + t332 * t367) * t325 + (t329 * t375 + t371 * t332) * t315 + (t329 * t374 + t332 * t370) * t314) * t333 + (t314 * t372 + t315 * t373 + t325 * t368) * t330) * t325 / 0.2e1 + (m(2) * (t321 ^ 2 + t323 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
