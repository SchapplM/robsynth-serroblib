% Calculate kinetic energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:17
% EndTime: 2019-03-09 01:35:17
% DurationCPUTime: 0.70s
% Computational Cost: add. (644->169), mult. (1378->275), div. (0->0), fcn. (1604->8), ass. (0->87)
t377 = cos(qJ(1));
t376 = sin(qJ(1));
t373 = sin(pkin(9));
t374 = cos(pkin(9));
t324 = -t373 * t376 - t374 * t377;
t375 = pkin(7) * t324;
t342 = sin(qJ(5));
t372 = Icges(6,4) * t342;
t344 = cos(qJ(5));
t371 = Icges(6,4) * t344;
t370 = t324 * t344;
t325 = t373 * t377 - t374 * t376;
t369 = t325 * t344;
t341 = sin(qJ(6));
t368 = t341 * t342;
t343 = cos(qJ(6));
t367 = t342 * t343;
t340 = qJD(2) * t376;
t366 = qJD(4) * t325 + t340;
t365 = qJD(5) * t324;
t364 = qJD(5) * t325;
t363 = qJD(5) * (-rSges(6,1) * t344 + rSges(6,2) * t342);
t362 = qJD(5) * (-pkin(5) * t344 - pkin(8) * t342);
t361 = qJD(6) * t344;
t360 = pkin(7) * qJD(1) * t325 + t366;
t330 = pkin(1) * t376 - qJ(2) * t377;
t359 = -pkin(2) * t376 - t330;
t358 = pkin(5) * t342 - pkin(8) * t344;
t357 = -qJD(2) * t377 + qJD(1) * (pkin(1) * t377 + qJ(2) * t376);
t356 = rSges(6,1) * t342 + rSges(6,2) * t344;
t355 = Icges(6,1) * t342 + t371;
t354 = Icges(6,2) * t344 + t372;
t353 = Icges(6,5) * t342 + Icges(6,6) * t344;
t299 = -Icges(6,6) * t324 + t325 * t354;
t301 = -Icges(6,5) * t324 + t325 * t355;
t352 = t299 * t344 + t301 * t342;
t300 = -Icges(6,6) * t325 - t324 * t354;
t302 = -Icges(6,5) * t325 - t324 * t355;
t351 = -t300 * t344 - t302 * t342;
t328 = Icges(6,2) * t342 - t371;
t329 = -Icges(6,1) * t344 + t372;
t350 = t328 * t344 + t329 * t342;
t349 = pkin(3) * t325 + qJ(4) * t324 + t359;
t348 = pkin(2) * qJD(1) * t377 + t357;
t347 = qJD(1) * (-pkin(3) * t324 + qJ(4) * t325) - qJD(4) * t324 + t348;
t345 = qJD(3) ^ 2;
t335 = -qJD(6) * t342 + qJD(1);
t333 = rSges(2,1) * t377 - rSges(2,2) * t376;
t331 = rSges(2,1) * t376 + rSges(2,2) * t377;
t327 = -Icges(6,5) * t344 + Icges(6,6) * t342;
t320 = -t342 * rSges(7,3) + (-rSges(7,1) * t343 + rSges(7,2) * t341) * t344;
t319 = -Icges(7,5) * t342 + (-Icges(7,1) * t343 + Icges(7,4) * t341) * t344;
t318 = -Icges(7,6) * t342 + (-Icges(7,4) * t343 + Icges(7,2) * t341) * t344;
t317 = -Icges(7,3) * t342 + (-Icges(7,5) * t343 + Icges(7,6) * t341) * t344;
t316 = qJD(1) * (rSges(3,1) * t377 + rSges(3,3) * t376) + t357;
t315 = t340 + (-rSges(3,1) * t376 + rSges(3,3) * t377 - t330) * qJD(1);
t312 = -t325 * t361 - t365;
t311 = t324 * t361 - t364;
t310 = t358 * t324;
t309 = t358 * t325;
t308 = -t324 * t367 - t325 * t341;
t307 = t324 * t368 - t325 * t343;
t306 = -t324 * t341 + t325 * t367;
t305 = -t324 * t343 - t325 * t368;
t304 = -t325 * rSges(6,3) - t324 * t356;
t303 = -t324 * rSges(6,3) + t325 * t356;
t298 = -Icges(6,3) * t325 - t324 * t353;
t297 = -Icges(6,3) * t324 + t325 * t353;
t296 = qJD(1) * (-rSges(4,1) * t324 - rSges(4,2) * t325) + t348;
t295 = t340 + (rSges(4,1) * t325 - rSges(4,2) * t324 + t359) * qJD(1);
t294 = rSges(7,1) * t308 + rSges(7,2) * t307 + rSges(7,3) * t370;
t293 = rSges(7,1) * t306 + rSges(7,2) * t305 - rSges(7,3) * t369;
t292 = Icges(7,1) * t308 + Icges(7,4) * t307 + Icges(7,5) * t370;
t291 = Icges(7,1) * t306 + Icges(7,4) * t305 - Icges(7,5) * t369;
t290 = Icges(7,4) * t308 + Icges(7,2) * t307 + Icges(7,6) * t370;
t289 = Icges(7,4) * t306 + Icges(7,2) * t305 - Icges(7,6) * t369;
t288 = Icges(7,5) * t308 + Icges(7,6) * t307 + Icges(7,3) * t370;
t287 = Icges(7,5) * t306 + Icges(7,6) * t305 - Icges(7,3) * t369;
t286 = qJD(1) * (rSges(5,2) * t324 + rSges(5,3) * t325) + t347;
t285 = (-rSges(5,2) * t325 + rSges(5,3) * t324 + t349) * qJD(1) + t366;
t284 = -qJD(3) + (t303 * t325 - t304 * t324) * qJD(5);
t283 = t324 * t363 + (t303 - t375) * qJD(1) + t347;
t282 = -t325 * t363 + (-t304 + t349) * qJD(1) + t360;
t281 = t324 * t362 + t335 * t293 - t312 * t320 + (t309 - t375) * qJD(1) + t347;
t280 = -t325 * t362 - t335 * t294 + t311 * t320 + (t310 + t349) * qJD(1) + t360;
t279 = -t311 * t293 + t312 * t294 - qJD(3) + (t309 * t325 + t310 * t324) * qJD(5);
t1 = t312 * ((-t287 * t369 + t305 * t289 + t306 * t291) * t312 + (-t288 * t369 + t290 * t305 + t292 * t306) * t311 + (t305 * t318 + t306 * t319 - t317 * t369) * t335) / 0.2e1 + t311 * ((t287 * t370 + t289 * t307 + t291 * t308) * t312 + (t288 * t370 + t307 * t290 + t308 * t292) * t311 + (t307 * t318 + t308 * t319 + t317 * t370) * t335) / 0.2e1 + t335 * ((-t287 * t312 - t288 * t311 - t317 * t335) * t342 + ((t289 * t341 - t291 * t343) * t312 + (t290 * t341 - t292 * t343) * t311 + (t318 * t341 - t319 * t343) * t335) * t344) / 0.2e1 - ((-t324 * t327 + t325 * t350) * qJD(1) + (t324 ^ 2 * t297 + (t351 * t325 + (t298 - t352) * t324) * t325) * qJD(5)) * t365 / 0.2e1 - ((-t324 * t350 - t325 * t327) * qJD(1) + (t325 ^ 2 * t298 + (t352 * t324 + (t297 - t351) * t325) * t324) * qJD(5)) * t364 / 0.2e1 + qJD(1) * ((t342 * t328 - t344 * t329) * qJD(1) + (-(t299 * t342 - t301 * t344) * t324 - (t342 * t300 - t344 * t302) * t325) * qJD(5)) / 0.2e1 + m(7) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(3) * (t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(4) * (t295 ^ 2 + t296 ^ 2 + t345) / 0.2e1 + m(5) * (t285 ^ 2 + t286 ^ 2 + t345) / 0.2e1 + m(6) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + (Icges(5,1) + Icges(2,3) + Icges(3,2) + Icges(4,3) + m(2) * (t331 ^ 2 + t333 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
