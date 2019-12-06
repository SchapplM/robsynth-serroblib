% Calculate kinetic energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:24
% EndTime: 2019-12-05 14:59:25
% DurationCPUTime: 0.48s
% Computational Cost: add. (845->157), mult. (2180->258), div. (0->0), fcn. (2703->10), ass. (0->74)
t359 = cos(qJ(4));
t340 = sin(pkin(9));
t341 = sin(pkin(8));
t358 = t341 * t340;
t347 = sin(qJ(4));
t357 = t341 * t347;
t342 = sin(pkin(7));
t344 = cos(pkin(8));
t356 = t342 * t344;
t345 = cos(pkin(7));
t355 = t345 * t340;
t343 = cos(pkin(9));
t354 = t345 * t343;
t353 = qJD(3) * t341;
t332 = qJD(2) * t342 + t345 * t353;
t351 = t341 * t359;
t333 = -qJD(2) * t345 + t342 * t353;
t338 = -qJD(3) * t344 + qJD(1);
t350 = qJD(1) ^ 2;
t348 = cos(qJ(5));
t346 = sin(qJ(5));
t331 = t343 * t351 - t344 * t347;
t330 = t343 * t357 + t344 * t359;
t329 = t342 * t340 + t344 * t354;
t328 = -t342 * t343 + t344 * t355;
t327 = t343 * t356 - t355;
t326 = t340 * t356 + t354;
t325 = qJD(4) * t358 + qJD(5) * t330;
t324 = t331 * t348 + t346 * t358;
t323 = -t331 * t346 + t348 * t358;
t322 = t329 * t359 + t345 * t357;
t321 = t329 * t347 - t345 * t351;
t320 = t327 * t359 + t342 * t357;
t319 = t327 * t347 - t342 * t351;
t318 = t331 * pkin(4) + t330 * pkin(6);
t317 = t331 * rSges(5,1) - t330 * rSges(5,2) + rSges(5,3) * t358;
t316 = Icges(5,1) * t331 - Icges(5,4) * t330 + Icges(5,5) * t358;
t315 = Icges(5,4) * t331 - Icges(5,2) * t330 + Icges(5,6) * t358;
t314 = Icges(5,5) * t331 - Icges(5,6) * t330 + Icges(5,3) * t358;
t313 = qJD(4) * t328 + qJD(5) * t321;
t312 = qJD(4) * t326 + qJD(5) * t319;
t311 = t322 * t348 + t328 * t346;
t310 = -t322 * t346 + t328 * t348;
t309 = t320 * t348 + t326 * t346;
t308 = -t320 * t346 + t326 * t348;
t307 = t322 * pkin(4) + t321 * pkin(6);
t306 = t320 * pkin(4) + t319 * pkin(6);
t305 = t324 * rSges(6,1) + t323 * rSges(6,2) + t330 * rSges(6,3);
t304 = Icges(6,1) * t324 + Icges(6,4) * t323 + Icges(6,5) * t330;
t303 = Icges(6,4) * t324 + Icges(6,2) * t323 + Icges(6,6) * t330;
t302 = Icges(6,5) * t324 + Icges(6,6) * t323 + Icges(6,3) * t330;
t301 = t322 * rSges(5,1) - t321 * rSges(5,2) + t328 * rSges(5,3);
t300 = t320 * rSges(5,1) - t319 * rSges(5,2) + t326 * rSges(5,3);
t299 = Icges(5,1) * t322 - Icges(5,4) * t321 + Icges(5,5) * t328;
t298 = Icges(5,1) * t320 - Icges(5,4) * t319 + Icges(5,5) * t326;
t297 = Icges(5,4) * t322 - Icges(5,2) * t321 + Icges(5,6) * t328;
t296 = Icges(5,4) * t320 - Icges(5,2) * t319 + Icges(5,6) * t326;
t295 = Icges(5,5) * t322 - Icges(5,6) * t321 + Icges(5,3) * t328;
t294 = Icges(5,5) * t320 - Icges(5,6) * t319 + Icges(5,3) * t326;
t293 = t311 * rSges(6,1) + t310 * rSges(6,2) + t321 * rSges(6,3);
t292 = t309 * rSges(6,1) + t308 * rSges(6,2) + t319 * rSges(6,3);
t291 = Icges(6,1) * t311 + Icges(6,4) * t310 + Icges(6,5) * t321;
t290 = Icges(6,1) * t309 + Icges(6,4) * t308 + Icges(6,5) * t319;
t289 = Icges(6,4) * t311 + Icges(6,2) * t310 + Icges(6,6) * t321;
t288 = Icges(6,4) * t309 + Icges(6,2) * t308 + Icges(6,6) * t319;
t287 = Icges(6,5) * t311 + Icges(6,6) * t310 + Icges(6,3) * t321;
t286 = Icges(6,5) * t309 + Icges(6,6) * t308 + Icges(6,3) * t319;
t285 = (t301 * t358 - t317 * t328) * qJD(4) + t333;
t284 = (-t300 * t358 + t317 * t326) * qJD(4) + t332;
t283 = (t300 * t328 - t301 * t326) * qJD(4) + t338;
t282 = t325 * t293 - t313 * t305 + (t307 * t358 - t318 * t328) * qJD(4) + t333;
t281 = -t325 * t292 + t312 * t305 + (-t306 * t358 + t318 * t326) * qJD(4) + t332;
t280 = t313 * t292 - t312 * t293 + (t306 * t328 - t307 * t326) * qJD(4) + t338;
t1 = m(2) * t350 / 0.2e1 + m(3) * (t350 + (t342 ^ 2 + t345 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t332 ^ 2 + t333 ^ 2 + t338 ^ 2) / 0.2e1 + m(5) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(6) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + t313 * ((t321 * t287 + t310 * t289 + t311 * t291) * t313 + (t321 * t286 + t310 * t288 + t311 * t290) * t312 + (t321 * t302 + t310 * t303 + t311 * t304) * t325) / 0.2e1 + t312 * ((t319 * t287 + t308 * t289 + t309 * t291) * t313 + (t319 * t286 + t308 * t288 + t309 * t290) * t312 + (t319 * t302 + t308 * t303 + t309 * t304) * t325) / 0.2e1 + t325 * ((t330 * t287 + t323 * t289 + t324 * t291) * t313 + (t330 * t286 + t323 * t288 + t324 * t290) * t312 + (t330 * t302 + t323 * t303 + t324 * t304) * t325) / 0.2e1 + (t328 * ((t328 * t295 - t321 * t297 + t322 * t299) * t328 + (t328 * t294 - t321 * t296 + t322 * t298) * t326 + (t328 * t314 - t321 * t315 + t322 * t316) * t358) + t326 * ((t326 * t295 - t319 * t297 + t320 * t299) * t328 + (t326 * t294 - t319 * t296 + t320 * t298) * t326 + (t326 * t314 - t319 * t315 + t320 * t316) * t358) + ((t295 * t358 - t330 * t297 + t331 * t299) * t328 + (t294 * t358 - t330 * t296 + t331 * t298) * t326 + (t314 * t358 - t330 * t315 + t331 * t316) * t358) * t358) * qJD(4) ^ 2 / 0.2e1;
T = t1;
