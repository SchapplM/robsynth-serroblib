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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:12:10
% EndTime: 2022-01-23 09:12:11
% DurationCPUTime: 0.95s
% Computational Cost: add. (790->130), mult. (1002->201), div. (0->0), fcn. (1007->8), ass. (0->67)
t358 = -Icges(6,5) - Icges(5,5);
t356 = -Icges(6,6) - Icges(5,6);
t361 = Icges(5,3) + Icges(6,3);
t360 = Icges(5,1) + Icges(6,1);
t359 = Icges(5,4) + Icges(6,4);
t357 = Icges(5,2) + Icges(6,2);
t307 = qJ(1) + pkin(7);
t305 = sin(t307);
t306 = cos(t307);
t313 = cos(qJ(4));
t309 = cos(pkin(8));
t311 = sin(qJ(4));
t337 = t309 * t311;
t292 = -t305 * t337 - t306 * t313;
t336 = t309 * t313;
t338 = t306 * t311;
t293 = t305 * t336 - t338;
t294 = t305 * t313 - t306 * t337;
t340 = t305 * t311;
t295 = t306 * t336 + t340;
t308 = sin(pkin(8));
t339 = t306 * t308;
t341 = t305 * t308;
t355 = (-t356 * t294 - t358 * t295 + t361 * t339) * t306 + (-t356 * t292 - t358 * t293 + t361 * t341) * t305;
t354 = t357 * t292 + t359 * t293 - t356 * t341;
t353 = t357 * t294 + t359 * t295 - t356 * t339;
t352 = t359 * t292 + t360 * t293 - t358 * t341;
t351 = t359 * t294 + t360 * t295 - t358 * t339;
t350 = t361 * t309 + (-t356 * t311 + t358 * t313) * t308;
t349 = t356 * t309 + (-t357 * t311 + t359 * t313) * t308;
t348 = t358 * t309 + (-t359 * t311 + t360 * t313) * t308;
t347 = t355 * t308;
t312 = sin(qJ(1));
t344 = t312 * pkin(1);
t343 = t313 * pkin(4);
t317 = qJ(5) * t308 + t343 * t309;
t335 = t293 * rSges(6,1) + t292 * rSges(6,2) + rSges(6,3) * t341 - pkin(4) * t338 + t317 * t305;
t334 = t295 * rSges(6,1) + t294 * rSges(6,2) + rSges(6,3) * t339 + pkin(4) * t340 + t317 * t306;
t333 = (-qJ(5) - rSges(6,3)) * t309 + (rSges(6,1) * t313 - rSges(6,2) * t311 + t343) * t308;
t314 = cos(qJ(1));
t304 = qJD(1) * t314 * pkin(1);
t332 = qJD(1) * (t306 * pkin(2) + t305 * qJ(3)) + t304;
t331 = qJD(4) * t308;
t324 = pkin(3) * t309 + pkin(6) * t308;
t330 = qJD(1) * t324 * t306 + t332;
t329 = (-t309 * rSges(5,3) + (rSges(5,1) * t313 - rSges(5,2) * t311) * t308) * t331;
t327 = -t305 * pkin(2) + t306 * qJ(3) - t344;
t323 = rSges(4,1) * t309 - rSges(4,2) * t308;
t302 = qJD(3) * t305;
t318 = t302 + (-t324 * t305 + t327) * qJD(1);
t315 = qJD(2) ^ 2;
t301 = -qJD(4) * t309 + qJD(1);
t300 = t314 * rSges(2,1) - t312 * rSges(2,2);
t299 = t312 * rSges(2,1) + t314 * rSges(2,2);
t282 = t304 + qJD(1) * (t306 * rSges(3,1) - t305 * rSges(3,2));
t281 = (-t305 * rSges(3,1) - t306 * rSges(3,2) - t344) * qJD(1);
t279 = t295 * rSges(5,1) + t294 * rSges(5,2) + rSges(5,3) * t339;
t277 = t293 * rSges(5,1) + t292 * rSges(5,2) + rSges(5,3) * t341;
t261 = qJD(1) * t305 * rSges(4,3) + (qJD(1) * t323 - qJD(3)) * t306 + t332;
t260 = t302 + (t306 * rSges(4,3) - t323 * t305 + t327) * qJD(1);
t259 = qJD(2) + (t277 * t306 - t279 * t305) * t331;
t258 = t301 * t279 + (-qJD(3) - t329) * t306 + t330;
t257 = -t301 * t277 + t305 * t329 + t318;
t256 = qJD(5) * t341 + t334 * t301 + (-t333 * t331 - qJD(3)) * t306 + t330;
t255 = -t335 * t301 + (t333 * t305 * qJD(4) + qJD(5) * t306) * t308 + t318;
t254 = -qJD(5) * t309 + qJD(2) + (-t334 * t305 + t335 * t306) * t331;
t1 = m(3) * (t281 ^ 2 + t282 ^ 2 + t315) / 0.2e1 + m(4) * (t260 ^ 2 + t261 ^ 2 + t315) / 0.2e1 + m(5) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + m(6) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + ((-t355 * t309 + ((-t353 * t311 + t351 * t313) * t306 + (-t354 * t311 + t352 * t313) * t305) * t308) * t331 + (t350 * t309 + (-t349 * t311 + t348 * t313) * t308) * t301) * t301 / 0.2e1 + (m(2) * (t299 ^ 2 + t300 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t309 ^ 2 + (Icges(4,1) * t308 + 0.2e1 * Icges(4,4) * t309) * t308) * qJD(1) ^ 2 / 0.2e1 + ((((t353 * t292 + t351 * t293) * t306 + (t354 * t292 + t352 * t293 + t347) * t305) * t331 + (t349 * t292 + t348 * t293 - t350 * t341) * t301) * t305 + (((t353 * t294 + t351 * t295 + t347) * t306 + (t354 * t294 + t352 * t295) * t305) * t331 + (t349 * t294 + t348 * t295 - t350 * t339) * t301) * t306) * t331 / 0.2e1;
T = t1;
