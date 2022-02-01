% Calculate kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:14
% EndTime: 2022-01-23 09:16:15
% DurationCPUTime: 0.86s
% Computational Cost: add. (994->213), mult. (1152->323), div. (0->0), fcn. (1161->10), ass. (0->99)
t353 = cos(pkin(8));
t383 = t353 ^ 2;
t350 = sin(pkin(9));
t381 = pkin(3) * t350;
t354 = qJ(3) + pkin(6);
t352 = cos(pkin(9));
t339 = t352 * pkin(3) + pkin(2);
t355 = sin(qJ(1));
t380 = t353 * t355;
t356 = cos(qJ(1));
t379 = t353 * t356;
t351 = sin(pkin(8));
t378 = t355 * t351;
t377 = t356 * t351;
t329 = pkin(2) * t353 + qJ(3) * t351 + pkin(1);
t346 = t355 * pkin(1);
t330 = -qJ(2) * t356 + t346;
t376 = -t329 * t355 - t330 + t346;
t322 = t339 * t353 + t351 * t354 + pkin(1);
t375 = t322 - t329;
t336 = qJ(2) + t381;
t349 = pkin(9) + qJ(4);
t340 = sin(t349);
t374 = pkin(4) * t340 - t336 + t381;
t343 = qJD(2) * t355;
t371 = qJD(3) * t351;
t373 = t356 * t371 + t343;
t344 = t355 * qJ(2);
t347 = t356 * pkin(1);
t372 = t347 + t344;
t370 = qJD(3) * t353;
t369 = qJD(4) * t351;
t368 = qJD(4) + qJD(5);
t328 = qJD(1) * t372;
t367 = qJD(1) * (t329 * t356 - t347) + t328 + t355 * t371;
t366 = t355 * t369;
t364 = t351 * t368;
t363 = qJD(1) * (t336 * t355 + t356 * t375 - t344) + t367;
t362 = rSges(3,1) * t353 - rSges(3,2) * t351;
t341 = cos(t349);
t317 = -t340 * t380 - t341 * t356;
t318 = -t340 * t356 + t341 * t380;
t319 = -t340 * t379 + t341 * t355;
t320 = t340 * t355 + t341 * t379;
t361 = (Icges(5,5) * t318 + Icges(5,6) * t317 + Icges(5,3) * t378) * t355 + (Icges(5,5) * t320 + Icges(5,6) * t319 + Icges(5,3) * t377) * t356;
t325 = pkin(4) * t341 + t339;
t348 = -pkin(7) - t354;
t360 = t325 * t353 - t348 * t351 - t322;
t359 = t361 * t351;
t358 = (-(qJ(2) - t336) * t356 - t375 * t355 + t376) * qJD(1) + t373;
t342 = qJ(5) + t349;
t338 = cos(t342);
t337 = sin(t342);
t335 = -t353 * qJD(4) + qJD(1);
t332 = rSges(2,1) * t356 - rSges(2,2) * t355;
t331 = rSges(2,1) * t355 + rSges(2,2) * t356;
t326 = -t353 * t368 + qJD(1);
t324 = t356 * t364;
t323 = t355 * t364;
t315 = t337 * t355 + t338 * t379;
t314 = -t337 * t379 + t338 * t355;
t313 = -t337 * t356 + t338 * t380;
t312 = -t337 * t380 - t338 * t356;
t311 = -rSges(5,3) * t353 + (rSges(5,1) * t341 - rSges(5,2) * t340) * t351;
t310 = -Icges(5,5) * t353 + (Icges(5,1) * t341 - Icges(5,4) * t340) * t351;
t309 = -Icges(5,6) * t353 + (Icges(5,4) * t341 - Icges(5,2) * t340) * t351;
t308 = -Icges(5,3) * t353 + (Icges(5,5) * t341 - Icges(5,6) * t340) * t351;
t307 = -rSges(6,3) * t353 + (rSges(6,1) * t338 - rSges(6,2) * t337) * t351;
t306 = -Icges(6,5) * t353 + (Icges(6,1) * t338 - Icges(6,4) * t337) * t351;
t305 = -Icges(6,6) * t353 + (Icges(6,4) * t338 - Icges(6,2) * t337) * t351;
t304 = -Icges(6,3) * t353 + (Icges(6,5) * t338 - Icges(6,6) * t337) * t351;
t303 = qJD(1) * t355 * rSges(3,3) + t328 + (qJD(1) * t362 - qJD(2)) * t356;
t302 = t343 + (t356 * rSges(3,3) - t355 * t362 - t330) * qJD(1);
t301 = (t348 + t354) * t353 + (t325 - t339) * t351;
t299 = rSges(5,1) * t320 + rSges(5,2) * t319 + rSges(5,3) * t377;
t298 = rSges(5,1) * t318 + rSges(5,2) * t317 + rSges(5,3) * t378;
t297 = Icges(5,1) * t320 + Icges(5,4) * t319 + Icges(5,5) * t377;
t296 = Icges(5,1) * t318 + Icges(5,4) * t317 + Icges(5,5) * t378;
t295 = Icges(5,4) * t320 + Icges(5,2) * t319 + Icges(5,6) * t377;
t294 = Icges(5,4) * t318 + Icges(5,2) * t317 + Icges(5,6) * t378;
t290 = rSges(6,1) * t315 + rSges(6,2) * t314 + rSges(6,3) * t377;
t289 = rSges(6,1) * t313 + rSges(6,2) * t312 + rSges(6,3) * t378;
t288 = Icges(6,1) * t315 + Icges(6,4) * t314 + Icges(6,5) * t377;
t287 = Icges(6,1) * t313 + Icges(6,4) * t312 + Icges(6,5) * t378;
t286 = Icges(6,4) * t315 + Icges(6,2) * t314 + Icges(6,6) * t377;
t285 = Icges(6,4) * t313 + Icges(6,2) * t312 + Icges(6,6) * t378;
t284 = Icges(6,5) * t315 + Icges(6,6) * t314 + Icges(6,3) * t377;
t283 = Icges(6,5) * t313 + Icges(6,6) * t312 + Icges(6,3) * t378;
t282 = t355 * t374 + t356 * t360 + t372;
t281 = t346 + (-qJ(2) - t374) * t356 + t360 * t355;
t280 = -qJD(2) * t356 + qJD(1) * ((t350 * t355 + t352 * t379) * rSges(4,1) + (-t350 * t379 + t352 * t355) * rSges(4,2) + rSges(4,3) * t377) + t367;
t279 = (-(-t350 * t356 + t352 * t380) * rSges(4,1) - (-t350 * t380 - t352 * t356) * rSges(4,2) - rSges(4,3) * t378 + t376) * qJD(1) + t373;
t278 = -t370 + (t298 * t356 - t299 * t355) * t369;
t277 = t299 * t335 + (-t311 * t369 - qJD(2)) * t356 + t363;
t276 = -t298 * t335 + t311 * t366 + t358;
t275 = -t370 + t289 * t324 - t290 * t323 + (t281 * t356 - t282 * t355) * t369;
t274 = t282 * t335 + t290 * t326 - t307 * t324 + (-t301 * t369 - qJD(2)) * t356 + t363;
t273 = -t281 * t335 - t289 * t326 + t301 * t366 + t307 * t323 + t358;
t1 = m(3) * (t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t383 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + m(5) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t335 * ((-t353 * t308 + (-t309 * t340 + t310 * t341) * t351) * t335 + (((-t295 * t340 + t297 * t341) * t356 + (-t294 * t340 + t296 * t341) * t355) * t351 - t361 * t353) * t369) / 0.2e1 + m(6) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t324 * ((t284 * t377 + t286 * t314 + t288 * t315) * t324 + (t283 * t377 + t285 * t314 + t287 * t315) * t323 + (t304 * t377 + t305 * t314 + t306 * t315) * t326) / 0.2e1 + t323 * ((t284 * t378 + t286 * t312 + t288 * t313) * t324 + (t283 * t378 + t285 * t312 + t287 * t313) * t323 + (t304 * t378 + t305 * t312 + t306 * t313) * t326) / 0.2e1 + t326 * ((-t283 * t323 - t284 * t324 - t304 * t326) * t353 + ((-t286 * t337 + t288 * t338) * t324 + (-t285 * t337 + t287 * t338) * t323 + (-t305 * t337 + t306 * t338) * t326) * t351) / 0.2e1 + (t356 * ((t308 * t377 + t309 * t319 + t310 * t320) * t335 + ((t294 * t319 + t296 * t320) * t355 + (t319 * t295 + t320 * t297 + t359) * t356) * t369) + t355 * ((t308 * t378 + t309 * t317 + t310 * t318) * t335 + ((t295 * t317 + t297 * t318) * t356 + (t317 * t294 + t318 * t296 + t359) * t355) * t369)) * t369 / 0.2e1 + (m(2) * (t331 ^ 2 + t332 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t383 + ((Icges(3,1) + Icges(4,1) * t352 ^ 2 + (-0.2e1 * Icges(4,4) * t352 + Icges(4,2) * t350) * t350) * t351 + 0.2e1 * (-Icges(4,5) * t352 + Icges(4,6) * t350 + Icges(3,4)) * t353) * t351) * qJD(1) ^ 2 / 0.2e1;
T = t1;
