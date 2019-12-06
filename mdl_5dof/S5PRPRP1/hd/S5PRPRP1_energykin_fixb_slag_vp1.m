% Calculate kinetic energy for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:17
% EndTime: 2019-12-05 15:28:18
% DurationCPUTime: 0.93s
% Computational Cost: add. (645->99), mult. (516->154), div. (0->0), fcn. (423->6), ass. (0->63)
t362 = Icges(5,4) - Icges(6,5);
t361 = Icges(5,1) + Icges(6,1);
t360 = Icges(5,2) + Icges(6,3);
t295 = pkin(8) + qJ(4);
t293 = cos(t295);
t359 = t362 * t293;
t291 = sin(t295);
t358 = t362 * t291;
t357 = Icges(6,4) + Icges(5,5);
t356 = Icges(5,6) - Icges(6,6);
t355 = t360 * t291 - t359;
t354 = t361 * t293 - t358;
t353 = rSges(6,1) + pkin(4);
t352 = rSges(6,3) + qJ(5);
t351 = Icges(6,2) + Icges(5,3);
t296 = pkin(7) + qJ(2);
t292 = sin(t296);
t294 = cos(t296);
t350 = t355 * t292 + t356 * t294;
t349 = -t356 * t292 + t355 * t294;
t348 = -t354 * t292 + t357 * t294;
t347 = t357 * t292 + t354 * t294;
t346 = -t360 * t293 - t358;
t345 = t361 * t291 + t359;
t344 = -t356 * t291 + t357 * t293;
t343 = t352 * t291 + t353 * t293;
t342 = t344 * t292 - t351 * t294;
t341 = t351 * t292 + t344 * t294;
t340 = t357 * t291 + t356 * t293;
t339 = t346 * t291 + t345 * t293;
t338 = t349 * t291 + t347 * t293;
t337 = -t350 * t291 + t348 * t293;
t298 = cos(pkin(8));
t332 = pkin(3) * t298;
t286 = pkin(2) * t292 - qJ(3) * t294;
t326 = pkin(6) * t294 - t332 * t292 - t286;
t325 = -rSges(6,2) * t294 + t343 * t292;
t324 = rSges(6,2) * t292 + t343 * t294;
t323 = qJD(4) * t292;
t322 = qJD(4) * t294;
t276 = qJD(2) * (pkin(2) * t294 + qJ(3) * t292);
t319 = -qJD(3) * t294 + qJD(2) * (pkin(6) * t292 + t332 * t294) + t276;
t297 = sin(pkin(8));
t318 = rSges(4,1) * t298 - rSges(4,2) * t297;
t317 = rSges(5,1) * t293 - rSges(5,2) * t291;
t302 = t352 * qJD(4) * t293 + (-t353 * qJD(4) + qJD(5)) * t291;
t301 = qJD(1) ^ 2;
t300 = qJD(2) ^ 2;
t289 = qJD(3) * t292;
t288 = rSges(3,1) * t294 - rSges(3,2) * t292;
t287 = rSges(3,1) * t292 + rSges(3,2) * t294;
t285 = rSges(5,1) * t291 + rSges(5,2) * t293;
t273 = rSges(5,3) * t292 + t317 * t294;
t271 = -rSges(5,3) * t294 + t317 * t292;
t255 = qJD(2) * t292 * rSges(4,3) + t276 + (qJD(2) * t318 - qJD(3)) * t294;
t254 = t289 + (t294 * rSges(4,3) - t318 * t292 - t286) * qJD(2);
t253 = qJD(1) + (t271 * t292 + t273 * t294) * qJD(4);
t252 = qJD(2) * t273 - t285 * t323 + t319;
t251 = -t285 * t322 + t289 + (-t271 + t326) * qJD(2);
t250 = -qJD(5) * t293 + qJD(1) + (t325 * t292 + t324 * t294) * qJD(4);
t249 = t324 * qJD(2) + t302 * t292 + t319;
t248 = t289 + t302 * t294 + (-t325 + t326) * qJD(2);
t1 = m(2) * t301 / 0.2e1 + m(3) * (t301 + (t287 ^ 2 + t288 ^ 2) * t300) / 0.2e1 + m(4) * (t254 ^ 2 + t255 ^ 2 + t301) / 0.2e1 + m(5) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + m(6) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t298 ^ 2 + (Icges(4,1) * t297 + 0.2e1 * Icges(4,4) * t298) * t297) * t300 / 0.2e1 + (((t348 * t291 + t350 * t293) * t294 + (t347 * t291 - t349 * t293) * t292) * qJD(4) + (t345 * t291 - t346 * t293) * qJD(2)) * qJD(2) / 0.2e1 + ((t341 * t292 ^ 2 + (t337 * t294 + (t338 - t342) * t292) * t294) * qJD(4) + (t292 * t340 + t294 * t339) * qJD(2)) * t323 / 0.2e1 - ((t342 * t294 ^ 2 + (t338 * t292 + (t337 - t341) * t294) * t292) * qJD(4) + (t292 * t339 - t294 * t340) * qJD(2)) * t322 / 0.2e1;
T = t1;
