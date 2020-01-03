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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:19
% EndTime: 2020-01-03 11:22:19
% DurationCPUTime: 0.55s
% Computational Cost: add. (634->145), mult. (1594->235), div. (0->0), fcn. (1929->10), ass. (0->75)
t342 = cos(pkin(7));
t369 = t342 ^ 2;
t338 = sin(pkin(8));
t339 = sin(pkin(7));
t367 = t338 * t339;
t341 = cos(pkin(8));
t366 = t339 * t341;
t344 = sin(qJ(1));
t365 = t339 * t344;
t346 = cos(qJ(1));
t364 = t339 * t346;
t363 = t344 * t338;
t362 = t344 * t341;
t361 = t346 * t338;
t360 = t346 * t341;
t331 = -t346 * pkin(1) - t344 * qJ(2);
t350 = pkin(2) * t342 + qJ(3) * t339;
t359 = t350 * t346 - t331;
t358 = qJD(2) * t346;
t357 = qJD(3) * t339;
t324 = t342 * t362 - t361;
t337 = sin(pkin(9));
t340 = cos(pkin(9));
t313 = t324 * t337 - t340 * t365;
t356 = qJD(5) * t313;
t326 = -t342 * t360 - t363;
t315 = t326 * t337 + t340 * t364;
t355 = qJD(5) * t315;
t325 = t342 * t361 - t362;
t354 = -t326 * pkin(3) + t325 * qJ(4) + t359;
t353 = t344 * t357 - t358;
t328 = -qJD(3) * t342 + qJD(4) * t367;
t323 = t342 * t363 + t360;
t352 = qJD(4) * t323 + t353;
t351 = rSges(3,1) * t342 - rSges(3,2) * t339;
t329 = qJD(1) * (t344 * pkin(1) - t346 * qJ(2));
t349 = -t346 * t357 + t329 + (qJD(1) * t350 - qJD(2)) * t344;
t348 = qJD(1) * (t324 * pkin(3) + t323 * qJ(4)) - qJD(4) * t325 + t349;
t345 = cos(qJ(5));
t343 = sin(qJ(5));
t332 = -t346 * rSges(2,1) + t344 * rSges(2,2);
t330 = t344 * rSges(2,1) + t346 * rSges(2,2);
t322 = -t342 * t337 + t340 * t366;
t321 = t337 * t366 + t342 * t340;
t317 = qJD(5) * t321 + qJD(1);
t316 = t326 * t340 - t337 * t364;
t314 = t324 * t340 + t337 * t365;
t312 = t322 * t345 + t343 * t367;
t311 = -t322 * t343 + t345 * t367;
t308 = -t358 + (t344 * rSges(3,3) + t351 * t346 - t331) * qJD(1);
t307 = -qJD(1) * t346 * rSges(3,3) + t329 + (qJD(1) * t351 - qJD(2)) * t344;
t306 = t316 * t345 - t325 * t343;
t305 = -t316 * t343 - t325 * t345;
t304 = t314 * t345 + t323 * t343;
t303 = -t314 * t343 + t323 * t345;
t302 = t312 * rSges(6,1) + t311 * rSges(6,2) + t321 * rSges(6,3);
t301 = Icges(6,1) * t312 + Icges(6,4) * t311 + Icges(6,5) * t321;
t300 = Icges(6,4) * t312 + Icges(6,2) * t311 + Icges(6,6) * t321;
t299 = Icges(6,5) * t312 + Icges(6,6) * t311 + Icges(6,3) * t321;
t298 = (-t326 * rSges(4,1) - t325 * rSges(4,2) + rSges(4,3) * t364 + t359) * qJD(1) + t353;
t297 = qJD(1) * (t324 * rSges(4,1) - t323 * rSges(4,2) + rSges(4,3) * t365) + t349;
t296 = t306 * rSges(6,1) + t305 * rSges(6,2) + t315 * rSges(6,3);
t295 = t304 * rSges(6,1) + t303 * rSges(6,2) + t313 * rSges(6,3);
t294 = Icges(6,1) * t306 + Icges(6,4) * t305 + Icges(6,5) * t315;
t293 = Icges(6,1) * t304 + Icges(6,4) * t303 + Icges(6,5) * t313;
t292 = Icges(6,4) * t306 + Icges(6,2) * t305 + Icges(6,6) * t315;
t291 = Icges(6,4) * t304 + Icges(6,2) * t303 + Icges(6,6) * t313;
t290 = Icges(6,5) * t306 + Icges(6,6) * t305 + Icges(6,3) * t315;
t289 = Icges(6,5) * t304 + Icges(6,6) * t303 + Icges(6,3) * t313;
t288 = (-t316 * rSges(5,1) + t315 * rSges(5,2) + t325 * rSges(5,3) + t354) * qJD(1) + t352;
t287 = qJD(1) * (t314 * rSges(5,1) - t313 * rSges(5,2) + t323 * rSges(5,3)) + t348;
t286 = (-t295 * t315 + t296 * t313) * qJD(5) + t328;
t285 = t302 * t355 - t317 * t296 + (-t316 * pkin(4) - t315 * pkin(6) + t354) * qJD(1) + t352;
t284 = qJD(1) * (t314 * pkin(4) + t313 * pkin(6)) + t317 * t295 - t302 * t356 + t348;
t1 = m(3) * (t307 ^ 2 + t308 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t369 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + m(5) * (t287 ^ 2 + t288 ^ 2 + t328 ^ 2) / 0.2e1 + m(6) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + t317 * ((t321 * t299 + t311 * t300 + t312 * t301) * t317 + ((t321 * t289 + t311 * t291 + t312 * t293) * t313 + (t321 * t290 + t311 * t292 + t312 * t294) * t315) * qJD(5)) / 0.2e1 + ((t313 * t299 + t303 * t300 + t304 * t301) * t317 + ((t313 * t289 + t303 * t291 + t304 * t293) * t313 + (t313 * t290 + t303 * t292 + t304 * t294) * t315) * qJD(5)) * t356 / 0.2e1 + ((t315 * t299 + t305 * t300 + t306 * t301) * t317 + ((t315 * t289 + t305 * t291 + t306 * t293) * t313 + (t315 * t290 + t305 * t292 + t306 * t294) * t315) * qJD(5)) * t355 / 0.2e1 + (m(2) * (t330 ^ 2 + t332 ^ 2) + Icges(2,3) + (Icges(5,1) * t322 + 0.2e1 * Icges(5,5) * t367) * t322 + (-0.2e1 * Icges(5,4) * t322 + Icges(5,2) * t321 - 0.2e1 * Icges(5,6) * t367) * t321 + (Icges(3,2) + Icges(4,3)) * t369 + ((Icges(4,1) * t341 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t341 + (Icges(4,2) + Icges(5,3)) * t338) * t338) * t339 + 0.2e1 * (-Icges(4,5) * t341 + Icges(4,6) * t338 + Icges(3,4)) * t342) * t339) * qJD(1) ^ 2 / 0.2e1;
T = t1;
