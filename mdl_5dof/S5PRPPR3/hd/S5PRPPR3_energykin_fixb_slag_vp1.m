% Calculate kinetic energy for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:16
% EndTime: 2019-12-05 15:26:17
% DurationCPUTime: 1.07s
% Computational Cost: add. (635->119), mult. (920->199), div. (0->0), fcn. (868->8), ass. (0->72)
t333 = cos(pkin(7));
t386 = t333 ^ 2;
t332 = sin(pkin(7));
t387 = t332 ^ 2;
t367 = t386 + t387;
t388 = qJD(2) * t367;
t395 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t331 = qJ(2) + pkin(8);
t329 = sin(t331);
t330 = cos(t331);
t336 = sin(qJ(2));
t338 = cos(qJ(2));
t394 = Icges(3,5) * t338 - Icges(3,6) * t336 + (-Icges(5,4) + Icges(4,5)) * t330 + (Icges(5,5) - Icges(4,6)) * t329;
t393 = t332 * t333;
t392 = t394 * t332 - t395 * t333;
t391 = t395 * t332 + t394 * t333;
t385 = qJD(2) ^ 2;
t384 = pkin(2) * t336;
t383 = t338 * pkin(2);
t381 = t330 * t332;
t380 = t330 * t333;
t335 = sin(qJ(5));
t379 = t332 * t335;
t337 = cos(qJ(5));
t378 = t332 * t337;
t377 = t333 * t335;
t376 = t333 * t337;
t328 = qJD(3) * t332;
t371 = qJD(4) * t329;
t375 = t333 * t371 + t328;
t374 = qJD(2) * t332;
t373 = qJD(2) * t333;
t372 = qJD(3) * t333;
t370 = qJD(5) * t329;
t369 = qJD(5) * t330;
t368 = qJD(1) + (-qJ(3) * t333 + t383 * t332) * t374 + (qJ(3) * t332 + t383 * t333) * t373;
t364 = -t329 * pkin(3) + t330 * qJ(4) - t384;
t363 = t332 * t371 - t372;
t362 = qJD(2) * (-t329 * rSges(4,1) - t330 * rSges(4,2) - t384);
t361 = t368 + (pkin(3) * t330 + qJ(4) * t329) * t388;
t345 = qJD(2) * (t329 * rSges(5,2) + t330 * rSges(5,3) + t364);
t340 = qJD(2) * (-pkin(6) * t329 + t364);
t326 = t336 * rSges(3,1) + t338 * rSges(3,2);
t319 = t332 * t369 - t373;
t318 = t333 * t369 + t374;
t317 = t329 * t379 - t376;
t316 = t329 * t378 + t377;
t315 = t329 * t377 + t378;
t314 = t329 * t376 - t379;
t293 = t329 * rSges(6,3) + (-rSges(6,1) * t335 - rSges(6,2) * t337) * t330;
t292 = Icges(6,5) * t329 + (-Icges(6,1) * t335 - Icges(6,4) * t337) * t330;
t291 = Icges(6,6) * t329 + (-Icges(6,4) * t335 - Icges(6,2) * t337) * t330;
t290 = Icges(6,3) * t329 + (-Icges(6,5) * t335 - Icges(6,6) * t337) * t330;
t287 = t333 * t362 + t328;
t286 = t332 * t362 - t372;
t285 = t317 * rSges(6,1) + t316 * rSges(6,2) + rSges(6,3) * t381;
t284 = t315 * rSges(6,1) + t314 * rSges(6,2) + rSges(6,3) * t380;
t283 = Icges(6,1) * t317 + Icges(6,4) * t316 + Icges(6,5) * t381;
t282 = Icges(6,1) * t315 + Icges(6,4) * t314 + Icges(6,5) * t380;
t281 = Icges(6,4) * t317 + Icges(6,2) * t316 + Icges(6,6) * t381;
t280 = Icges(6,4) * t315 + Icges(6,2) * t314 + Icges(6,6) * t380;
t279 = Icges(6,5) * t317 + Icges(6,6) * t316 + Icges(6,3) * t381;
t278 = Icges(6,5) * t315 + Icges(6,6) * t314 + Icges(6,3) * t380;
t277 = qJD(1) + (rSges(3,1) * t338 - rSges(3,2) * t336) * t388;
t276 = t333 * t345 + t375;
t275 = t332 * t345 + t363;
t274 = t368 + (rSges(4,1) * t330 - rSges(4,2) * t329) * t388;
t273 = -t285 * t370 + t319 * t293 + t333 * t340 + t375;
t272 = t284 * t370 - t318 * t293 + t332 * t340 + t363;
t271 = -qJD(4) * t330 + t361 + (-rSges(5,2) * t330 + rSges(5,3) * t329) * t388;
t270 = -t319 * t284 + t318 * t285 + (pkin(6) * t388 - qJD(4)) * t330 + t361;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t367 * t385 * t326 ^ 2 + t277 ^ 2) / 0.2e1 + m(4) * (t274 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(5) * (t271 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(6) * (t270 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + t318 * ((t278 * t380 + t314 * t280 + t315 * t282) * t318 + (t279 * t380 + t314 * t281 + t315 * t283) * t319 + (t290 * t380 + t314 * t291 + t315 * t292) * t370) / 0.2e1 + t319 * ((t278 * t381 + t316 * t280 + t317 * t282) * t318 + (t279 * t381 + t316 * t281 + t317 * t283) * t319 + (t290 * t381 + t316 * t291 + t317 * t292) * t370) / 0.2e1 + ((t278 * t318 + t279 * t319 + t290 * t370) * t329 + ((-t280 * t337 - t282 * t335) * t318 + (-t281 * t337 - t283 * t335) * t319 + (-t291 * t337 - t292 * t335) * t370) * t330) * t370 / 0.2e1 + (t391 * t387 - t392 * t393) * t332 * t385 / 0.2e1 - (t392 * t386 - t391 * t393) * t333 * t385 / 0.2e1;
T = t1;
