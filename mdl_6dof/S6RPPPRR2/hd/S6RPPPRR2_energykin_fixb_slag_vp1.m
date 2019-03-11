% Calculate kinetic energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:29
% EndTime: 2019-03-09 01:31:30
% DurationCPUTime: 0.76s
% Computational Cost: add. (918->175), mult. (801->276), div. (0->0), fcn. (730->10), ass. (0->93)
t348 = qJ(1) + pkin(9);
t344 = sin(t348);
t346 = cos(t348);
t388 = qJ(4) * qJD(1) * t346 + qJD(4) * t344;
t349 = sin(pkin(10));
t386 = pkin(4) * t349;
t353 = sin(qJ(1));
t385 = t353 * pkin(1);
t347 = pkin(10) + qJ(5);
t343 = sin(t347);
t384 = Icges(6,4) * t343;
t345 = cos(t347);
t383 = Icges(6,4) * t345;
t382 = t344 * t345;
t352 = sin(qJ(6));
t381 = t344 * t352;
t354 = cos(qJ(6));
t380 = t344 * t354;
t379 = t345 * t346;
t378 = t346 * t352;
t377 = t346 * t354;
t355 = cos(qJ(1));
t342 = qJD(1) * t355 * pkin(1);
t375 = qJD(1) * (pkin(2) * t346 + qJ(3) * t344) + t342;
t341 = qJD(3) * t344;
t374 = qJD(4) * t346 + t341;
t373 = qJD(5) * t344;
t372 = qJD(5) * t346;
t371 = qJD(6) * t345;
t370 = -pkin(2) * t344 + qJ(3) * t346 - t385;
t369 = pkin(5) * t343 - pkin(8) * t345;
t368 = -qJD(3) * t346 + t375;
t350 = cos(pkin(10));
t367 = rSges(5,1) * t349 + rSges(5,2) * t350;
t366 = rSges(6,1) * t343 + rSges(6,2) * t345;
t365 = qJD(1) * (pkin(7) * t346 + t344 * t386) + t375 + t388;
t364 = Icges(6,1) * t343 + t383;
t363 = Icges(6,2) * t345 + t384;
t362 = Icges(6,5) * t343 + Icges(6,6) * t345;
t306 = Icges(6,6) * t346 + t344 * t363;
t308 = Icges(6,5) * t346 + t344 * t364;
t361 = -t306 * t345 - t308 * t343;
t307 = Icges(6,6) * t344 - t346 * t363;
t309 = Icges(6,5) * t344 - t346 * t364;
t360 = t307 * t345 + t309 * t343;
t330 = -Icges(6,2) * t343 + t383;
t331 = Icges(6,1) * t345 - t384;
t359 = t330 * t345 + t331 * t343;
t358 = t346 * t386 + t370 + (-pkin(7) - qJ(4)) * t344;
t356 = qJD(2) ^ 2;
t337 = qJD(6) * t343 + qJD(1);
t336 = rSges(2,1) * t355 - rSges(2,2) * t353;
t335 = rSges(2,1) * t353 + rSges(2,2) * t355;
t334 = pkin(5) * t345 + pkin(8) * t343;
t333 = rSges(6,1) * t345 - rSges(6,2) * t343;
t329 = Icges(6,5) * t345 - Icges(6,6) * t343;
t327 = -t344 * t371 + t372;
t326 = t346 * t371 + t373;
t325 = -t343 * t377 + t381;
t324 = t343 * t378 + t380;
t323 = t343 * t380 + t378;
t322 = -t343 * t381 + t377;
t321 = t342 + qJD(1) * (rSges(3,1) * t346 - rSges(3,2) * t344);
t320 = (-rSges(3,1) * t344 - rSges(3,2) * t346 - t385) * qJD(1);
t319 = t369 * t346;
t318 = t369 * t344;
t316 = t343 * rSges(7,3) + (rSges(7,1) * t354 - rSges(7,2) * t352) * t345;
t315 = Icges(7,5) * t343 + (Icges(7,1) * t354 - Icges(7,4) * t352) * t345;
t314 = Icges(7,6) * t343 + (Icges(7,4) * t354 - Icges(7,2) * t352) * t345;
t313 = Icges(7,3) * t343 + (Icges(7,5) * t354 - Icges(7,6) * t352) * t345;
t311 = t344 * rSges(6,3) - t346 * t366;
t310 = t346 * rSges(6,3) + t344 * t366;
t305 = Icges(6,3) * t344 - t346 * t362;
t304 = Icges(6,3) * t346 + t344 * t362;
t303 = qJD(1) * (-rSges(4,2) * t346 + rSges(4,3) * t344) + t368;
t302 = t341 + (rSges(4,2) * t344 + rSges(4,3) * t346 + t370) * qJD(1);
t301 = rSges(7,1) * t325 + rSges(7,2) * t324 + rSges(7,3) * t379;
t300 = rSges(7,1) * t323 + rSges(7,2) * t322 - rSges(7,3) * t382;
t299 = Icges(7,1) * t325 + Icges(7,4) * t324 + Icges(7,5) * t379;
t298 = Icges(7,1) * t323 + Icges(7,4) * t322 - Icges(7,5) * t382;
t297 = Icges(7,4) * t325 + Icges(7,2) * t324 + Icges(7,6) * t379;
t296 = Icges(7,4) * t323 + Icges(7,2) * t322 - Icges(7,6) * t382;
t295 = Icges(7,5) * t325 + Icges(7,6) * t324 + Icges(7,3) * t379;
t294 = Icges(7,5) * t323 + Icges(7,6) * t322 - Icges(7,3) * t382;
t293 = qJD(1) * (t346 * rSges(5,3) + t344 * t367) + t368 + t388;
t292 = (t367 * t346 + (-rSges(5,3) - qJ(4)) * t344 + t370) * qJD(1) + t374;
t291 = qJD(2) + (-t310 * t344 + t311 * t346) * qJD(5);
t290 = qJD(1) * t310 + (-qJD(5) * t333 - qJD(3)) * t346 + t365;
t289 = t333 * t373 + (-t311 + t358) * qJD(1) + t374;
t288 = -t326 * t300 + t327 * t301 + qJD(2) + (-t318 * t344 - t319 * t346) * qJD(5);
t287 = qJD(1) * t318 + t337 * t300 - t327 * t316 + (-qJD(5) * t334 - qJD(3)) * t346 + t365;
t286 = t334 * t373 - t337 * t301 + t326 * t316 + (t319 + t358) * qJD(1) + t374;
t1 = m(6) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(7) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(3) * (t320 ^ 2 + t321 ^ 2 + t356) / 0.2e1 + m(4) * (t302 ^ 2 + t303 ^ 2 + t356) / 0.2e1 + m(5) * (t292 ^ 2 + t293 ^ 2 + t356) / 0.2e1 + t326 * ((t294 * t379 + t296 * t324 + t298 * t325) * t327 + (t295 * t379 + t324 * t297 + t325 * t299) * t326 + (t313 * t379 + t314 * t324 + t315 * t325) * t337) / 0.2e1 + t337 * ((t294 * t327 + t295 * t326 + t313 * t337) * t343 + ((-t296 * t352 + t298 * t354) * t327 + (-t297 * t352 + t299 * t354) * t326 + (-t314 * t352 + t315 * t354) * t337) * t345) / 0.2e1 + ((t346 * t329 + t344 * t359) * qJD(1) + (t346 ^ 2 * t304 + (t360 * t344 + (t305 - t361) * t346) * t344) * qJD(5)) * t372 / 0.2e1 + ((t344 * t329 - t346 * t359) * qJD(1) + (t344 ^ 2 * t305 + (t361 * t346 + (t304 - t360) * t344) * t346) * qJD(5)) * t373 / 0.2e1 + qJD(1) * ((-t343 * t330 + t345 * t331) * qJD(1) + ((-t306 * t343 + t308 * t345) * t346 + (-t343 * t307 + t345 * t309) * t344) * qJD(5)) / 0.2e1 + t327 * ((-t294 * t382 + t322 * t296 + t323 * t298) * t327 + (-t295 * t382 + t297 * t322 + t299 * t323) * t326 + (-t313 * t382 + t314 * t322 + t315 * t323) * t337) / 0.2e1 + (m(2) * (t335 ^ 2 + t336 ^ 2) + Icges(3,3) + Icges(4,1) + Icges(5,1) * t350 ^ 2 + (-0.2e1 * Icges(5,4) * t350 + Icges(5,2) * t349) * t349 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
