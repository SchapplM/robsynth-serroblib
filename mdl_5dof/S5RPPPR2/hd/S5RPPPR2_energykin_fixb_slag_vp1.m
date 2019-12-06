% Calculate kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:56
% EndTime: 2019-12-05 17:30:57
% DurationCPUTime: 0.48s
% Computational Cost: add. (634->144), mult. (1594->235), div. (0->0), fcn. (1929->10), ass. (0->75)
t346 = cos(pkin(7));
t373 = t346 ^ 2;
t342 = sin(pkin(8));
t343 = sin(pkin(7));
t371 = t342 * t343;
t345 = cos(pkin(8));
t370 = t343 * t345;
t348 = sin(qJ(1));
t369 = t343 * t348;
t350 = cos(qJ(1));
t368 = t343 * t350;
t367 = t348 * t342;
t366 = t348 * t345;
t365 = t350 * t342;
t364 = t350 * t345;
t333 = t350 * pkin(1) + t348 * qJ(2);
t354 = pkin(2) * t346 + qJ(3) * t343;
t363 = -t354 * t350 - t333;
t362 = qJD(1) * (-t348 * pkin(1) + t350 * qJ(2)) + qJD(2) * t348;
t361 = qJD(3) * t343;
t326 = -t346 * t366 + t365;
t341 = sin(pkin(9));
t344 = cos(pkin(9));
t315 = t326 * t341 + t344 * t369;
t360 = qJD(5) * t315;
t328 = t346 * t364 + t367;
t317 = t328 * t341 - t344 * t368;
t359 = qJD(5) * t317;
t327 = t346 * t365 - t366;
t358 = -t328 * pkin(3) - t327 * qJ(4) + t363;
t330 = -qJD(3) * t346 + qJD(4) * t371;
t357 = -qJD(1) * t354 * t348 + t350 * t361 + t362;
t340 = qJD(2) * t350;
t356 = -t348 * t361 + t340;
t355 = -rSges(3,1) * t346 + rSges(3,2) * t343;
t325 = t346 * t367 + t364;
t353 = -qJD(4) * t325 + t356;
t352 = qJD(1) * (t326 * pkin(3) - t325 * qJ(4)) + qJD(4) * t327 + t357;
t349 = cos(qJ(5));
t347 = sin(qJ(5));
t334 = t350 * rSges(2,1) - t348 * rSges(2,2);
t332 = -t348 * rSges(2,1) - t350 * rSges(2,2);
t324 = -t346 * t341 + t344 * t370;
t323 = t341 * t370 + t346 * t344;
t319 = qJD(5) * t323 + qJD(1);
t318 = t328 * t344 + t341 * t368;
t316 = t326 * t344 - t341 * t369;
t314 = t324 * t349 + t347 * t371;
t313 = -t324 * t347 + t349 * t371;
t310 = t340 + (-t348 * rSges(3,3) + t350 * t355 - t333) * qJD(1);
t309 = qJD(1) * (t350 * rSges(3,3) + t348 * t355) + t362;
t308 = t318 * t349 + t327 * t347;
t307 = -t318 * t347 + t327 * t349;
t306 = t316 * t349 - t325 * t347;
t305 = -t316 * t347 - t325 * t349;
t304 = t314 * rSges(6,1) + t313 * rSges(6,2) + t323 * rSges(6,3);
t303 = Icges(6,1) * t314 + Icges(6,4) * t313 + Icges(6,5) * t323;
t302 = Icges(6,4) * t314 + Icges(6,2) * t313 + Icges(6,6) * t323;
t301 = Icges(6,5) * t314 + Icges(6,6) * t313 + Icges(6,3) * t323;
t300 = (-t328 * rSges(4,1) + t327 * rSges(4,2) - rSges(4,3) * t368 + t363) * qJD(1) + t356;
t299 = qJD(1) * (t326 * rSges(4,1) + t325 * rSges(4,2) - rSges(4,3) * t369) + t357;
t298 = t308 * rSges(6,1) + t307 * rSges(6,2) + t317 * rSges(6,3);
t297 = t306 * rSges(6,1) + t305 * rSges(6,2) + t315 * rSges(6,3);
t296 = Icges(6,1) * t308 + Icges(6,4) * t307 + Icges(6,5) * t317;
t295 = Icges(6,1) * t306 + Icges(6,4) * t305 + Icges(6,5) * t315;
t294 = Icges(6,4) * t308 + Icges(6,2) * t307 + Icges(6,6) * t317;
t293 = Icges(6,4) * t306 + Icges(6,2) * t305 + Icges(6,6) * t315;
t292 = Icges(6,5) * t308 + Icges(6,6) * t307 + Icges(6,3) * t317;
t291 = Icges(6,5) * t306 + Icges(6,6) * t305 + Icges(6,3) * t315;
t290 = (-t318 * rSges(5,1) + t317 * rSges(5,2) - t327 * rSges(5,3) + t358) * qJD(1) + t353;
t289 = qJD(1) * (t316 * rSges(5,1) - t315 * rSges(5,2) - t325 * rSges(5,3)) + t352;
t288 = (-t297 * t317 + t298 * t315) * qJD(5) + t330;
t287 = t304 * t359 - t319 * t298 + (-t318 * pkin(4) - t317 * pkin(6) + t358) * qJD(1) + t353;
t286 = qJD(1) * (t316 * pkin(4) + t315 * pkin(6)) + t319 * t297 - t304 * t360 + t352;
t1 = m(3) * (t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t373 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + m(5) * (t289 ^ 2 + t290 ^ 2 + t330 ^ 2) / 0.2e1 + m(6) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + t319 * ((t323 * t301 + t313 * t302 + t314 * t303) * t319 + ((t323 * t291 + t313 * t293 + t314 * t295) * t315 + (t323 * t292 + t313 * t294 + t314 * t296) * t317) * qJD(5)) / 0.2e1 + ((t315 * t301 + t305 * t302 + t306 * t303) * t319 + ((t315 * t291 + t305 * t293 + t306 * t295) * t315 + (t315 * t292 + t305 * t294 + t306 * t296) * t317) * qJD(5)) * t360 / 0.2e1 + ((t317 * t301 + t307 * t302 + t308 * t303) * t319 + ((t317 * t291 + t307 * t293 + t308 * t295) * t315 + (t317 * t292 + t307 * t294 + t308 * t296) * t317) * qJD(5)) * t359 / 0.2e1 + (m(2) * (t332 ^ 2 + t334 ^ 2) + Icges(2,3) + (Icges(5,1) * t324 + 0.2e1 * Icges(5,5) * t371) * t324 + (-0.2e1 * Icges(5,4) * t324 + Icges(5,2) * t323 - 0.2e1 * Icges(5,6) * t371) * t323 + (Icges(3,2) + Icges(4,3)) * t373 + ((Icges(4,1) * t345 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t345 + (Icges(4,2) + Icges(5,3)) * t342) * t342) * t343 + 0.2e1 * (-Icges(4,5) * t345 + Icges(4,6) * t342 + Icges(3,4)) * t346) * t343) * qJD(1) ^ 2 / 0.2e1;
T = t1;
