% Calculate kinetic energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:34
% EndTime: 2019-12-05 17:35:35
% DurationCPUTime: 0.93s
% Computational Cost: add. (790->129), mult. (1002->202), div. (0->0), fcn. (1007->8), ass. (0->69)
t360 = Icges(5,1) + Icges(6,1);
t359 = Icges(5,4) + Icges(6,4);
t358 = -Icges(6,5) - Icges(5,5);
t357 = Icges(5,2) + Icges(6,2);
t356 = -Icges(6,6) - Icges(5,6);
t355 = -Icges(6,3) - Icges(5,3);
t308 = sin(pkin(8));
t307 = qJ(1) + pkin(7);
t305 = sin(t307);
t306 = cos(t307);
t313 = cos(qJ(4));
t309 = cos(pkin(8));
t311 = sin(qJ(4));
t334 = t309 * t311;
t292 = t305 * t334 + t306 * t313;
t333 = t309 * t313;
t335 = t306 * t311;
t293 = -t305 * t333 + t335;
t294 = t305 * t313 - t306 * t334;
t337 = t305 * t311;
t295 = t306 * t333 + t337;
t336 = t306 * t308;
t338 = t305 * t308;
t346 = (t356 * t294 + t358 * t295 + t355 * t336) * t306 + (-t356 * t292 - t358 * t293 + t355 * t338) * t305;
t354 = t346 * t308;
t353 = t357 * t292 + t359 * t293 + t356 * t338;
t352 = t357 * t294 + t359 * t295 - t356 * t336;
t351 = -t359 * t292 - t360 * t293 - t358 * t338;
t350 = t359 * t294 + t360 * t295 - t358 * t336;
t349 = -t355 * t309 + (-t356 * t311 + t358 * t313) * t308;
t348 = t356 * t309 + (-t357 * t311 + t359 * t313) * t308;
t347 = t358 * t309 + (-t359 * t311 + t360 * t313) * t308;
t340 = t313 * pkin(4);
t345 = qJ(5) * t308 + t309 * t340;
t312 = sin(qJ(1));
t342 = t312 * pkin(1);
t314 = cos(qJ(1));
t341 = t314 * pkin(1);
t332 = t293 * rSges(6,1) + t292 * rSges(6,2) - rSges(6,3) * t338 + pkin(4) * t335 - t345 * t305;
t331 = -t295 * rSges(6,1) - t294 * rSges(6,2) - rSges(6,3) * t336 - pkin(4) * t337 - t345 * t306;
t330 = qJD(1) * (-t305 * pkin(2) + t306 * qJ(3)) + qJD(3) * t305;
t329 = qJD(4) * t308;
t328 = t305 * t329;
t327 = t306 * t329;
t326 = -t306 * pkin(2) - t305 * qJ(3) - t341;
t325 = qJD(4) * ((-qJ(5) - rSges(6,3)) * t309 + (rSges(6,1) * t313 - rSges(6,2) * t311 + t340) * t308);
t322 = pkin(3) * t309 + pkin(6) * t308;
t321 = -rSges(4,1) * t309 + rSges(4,2) * t308;
t318 = t330 + (-t305 * t322 - t342) * qJD(1);
t303 = qJD(3) * t306;
t317 = t303 + (-t322 * t306 + t326) * qJD(1);
t315 = qJD(2) ^ 2;
t301 = -qJD(4) * t309 + qJD(1);
t300 = t314 * rSges(2,1) - t312 * rSges(2,2);
t299 = -t312 * rSges(2,1) - t314 * rSges(2,2);
t291 = -t309 * rSges(5,3) + (rSges(5,1) * t313 - rSges(5,2) * t311) * t308;
t282 = (-t306 * rSges(3,1) + t305 * rSges(3,2) - t341) * qJD(1);
t281 = (-t305 * rSges(3,1) - t306 * rSges(3,2) - t342) * qJD(1);
t279 = t295 * rSges(5,1) + t294 * rSges(5,2) + rSges(5,3) * t336;
t277 = t293 * rSges(5,1) + t292 * rSges(5,2) - rSges(5,3) * t338;
t261 = t303 + (-t305 * rSges(4,3) + t306 * t321 + t326) * qJD(1);
t260 = (t306 * rSges(4,3) + t305 * t321 - t342) * qJD(1) + t330;
t259 = qJD(2) + (-t277 * t306 - t279 * t305) * t329;
t258 = -t301 * t279 + t291 * t327 + t317;
t257 = t301 * t277 + t291 * t328 + t318;
t256 = t331 * t301 + (-qJD(5) * t305 + t306 * t325) * t308 + t317;
t255 = t332 * t301 + (qJD(5) * t306 + t305 * t325) * t308 + t318;
t254 = -qJD(5) * t309 + qJD(2) + (t305 * t331 - t306 * t332) * t329;
t1 = m(3) * (t281 ^ 2 + t282 ^ 2 + t315) / 0.2e1 + m(4) * (t260 ^ 2 + t261 ^ 2 + t315) / 0.2e1 + m(5) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + m(6) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + ((t346 * t309 + ((-t352 * t311 + t350 * t313) * t306 + (t353 * t311 + t351 * t313) * t305) * t308) * t329 + (t349 * t309 + (-t348 * t311 + t347 * t313) * t308) * t301) * t301 / 0.2e1 - (((t352 * t292 + t350 * t293) * t306 + (-t353 * t292 + t351 * t293 + t354) * t305) * t329 + (t348 * t292 + t347 * t293 + t349 * t338) * t301) * t328 / 0.2e1 + (((t352 * t294 + t350 * t295 - t354) * t306 + (-t353 * t294 + t351 * t295) * t305) * t329 + (t348 * t294 + t347 * t295 - t349 * t336) * t301) * t327 / 0.2e1 + (m(2) * (t299 ^ 2 + t300 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t309 ^ 2 + (Icges(4,1) * t308 + 0.2e1 * Icges(4,4) * t309) * t308) * qJD(1) ^ 2 / 0.2e1;
T = t1;
