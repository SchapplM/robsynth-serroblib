% Calculate kinetic energy for
% S5PRPRP2
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:29
% EndTime: 2019-12-05 15:30:30
% DurationCPUTime: 0.92s
% Computational Cost: add. (778->122), mult. (976->193), div. (0->0), fcn. (995->6), ass. (0->62)
t351 = -Icges(6,5) - Icges(5,5);
t349 = -Icges(6,6) - Icges(5,6);
t354 = Icges(5,3) + Icges(6,3);
t353 = Icges(5,1) + Icges(6,1);
t352 = Icges(5,4) + Icges(6,4);
t350 = Icges(5,2) + Icges(6,2);
t305 = pkin(7) + qJ(2);
t303 = sin(t305);
t304 = cos(t305);
t310 = cos(qJ(4));
t307 = cos(pkin(8));
t309 = sin(qJ(4));
t331 = t307 * t309;
t291 = -t303 * t331 - t304 * t310;
t330 = t307 * t310;
t332 = t304 * t309;
t292 = t303 * t330 - t332;
t293 = t303 * t310 - t304 * t331;
t334 = t303 * t309;
t294 = t304 * t330 + t334;
t306 = sin(pkin(8));
t333 = t304 * t306;
t335 = t303 * t306;
t348 = (-t349 * t293 - t351 * t294 + t354 * t333) * t304 + (-t349 * t291 - t351 * t292 + t354 * t335) * t303;
t347 = t350 * t291 + t352 * t292 - t349 * t335;
t346 = t350 * t293 + t352 * t294 - t349 * t333;
t345 = t352 * t291 + t353 * t292 - t351 * t335;
t344 = t352 * t293 + t353 * t294 - t351 * t333;
t343 = t354 * t307 + (-t349 * t309 + t351 * t310) * t306;
t342 = t349 * t307 + (-t350 * t309 + t352 * t310) * t306;
t341 = t351 * t307 + (-t352 * t309 + t353 * t310) * t306;
t340 = t348 * t306;
t337 = t310 * pkin(4);
t313 = qJ(5) * t306 + t307 * t337;
t329 = t292 * rSges(6,1) + t291 * rSges(6,2) + rSges(6,3) * t335 - pkin(4) * t332 + t303 * t313;
t328 = t294 * rSges(6,1) + t293 * rSges(6,2) + rSges(6,3) * t333 + pkin(4) * t334 + t304 * t313;
t327 = (-qJ(5) - rSges(6,3)) * t307 + (rSges(6,1) * t310 - rSges(6,2) * t309 + t337) * t306;
t296 = qJD(2) * (t304 * pkin(2) + t303 * qJ(3));
t320 = pkin(3) * t307 + pkin(6) * t306;
t326 = qJD(2) * t320 * t304 + t296;
t325 = qJD(4) * t306;
t324 = (-t307 * rSges(5,3) + (rSges(5,1) * t310 - rSges(5,2) * t309) * t306) * t325;
t319 = rSges(4,1) * t307 - rSges(4,2) * t306;
t297 = t303 * pkin(2) - t304 * qJ(3);
t301 = qJD(3) * t303;
t316 = t301 + (-t320 * t303 - t297) * qJD(2);
t312 = qJD(1) ^ 2;
t311 = qJD(2) ^ 2;
t300 = -qJD(4) * t307 + qJD(2);
t299 = t304 * rSges(3,1) - t303 * rSges(3,2);
t298 = t303 * rSges(3,1) + t304 * rSges(3,2);
t280 = qJD(2) * t303 * rSges(4,3) + t296 + (qJD(2) * t319 - qJD(3)) * t304;
t279 = t301 + (t304 * rSges(4,3) - t303 * t319 - t297) * qJD(2);
t278 = t294 * rSges(5,1) + t293 * rSges(5,2) + rSges(5,3) * t333;
t276 = t292 * rSges(5,1) + t291 * rSges(5,2) + rSges(5,3) * t335;
t260 = qJD(1) + (t276 * t304 - t278 * t303) * t325;
t259 = t300 * t278 + (-qJD(3) - t324) * t304 + t326;
t258 = -t300 * t276 + t303 * t324 + t316;
t257 = qJD(5) * t335 + t328 * t300 + (-t325 * t327 - qJD(3)) * t304 + t326;
t256 = -t329 * t300 + (qJD(4) * t303 * t327 + qJD(5) * t304) * t306 + t316;
t255 = -qJD(5) * t307 + qJD(1) + (-t303 * t328 + t304 * t329) * t325;
t1 = m(2) * t312 / 0.2e1 + m(3) * (t312 + (t298 ^ 2 + t299 ^ 2) * t311) / 0.2e1 + m(4) * (t279 ^ 2 + t280 ^ 2 + t312) / 0.2e1 + m(5) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(6) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + ((-t348 * t307 + ((-t346 * t309 + t344 * t310) * t304 + (-t347 * t309 + t345 * t310) * t303) * t306) * t325 + (t343 * t307 + (-t342 * t309 + t341 * t310) * t306) * t300) * t300 / 0.2e1 + (Icges(3,3) + Icges(4,2) * t307 ^ 2 + (Icges(4,1) * t306 + 0.2e1 * Icges(4,4) * t307) * t306) * t311 / 0.2e1 + ((((t346 * t291 + t344 * t292) * t304 + (t347 * t291 + t345 * t292 + t340) * t303) * t325 + (t342 * t291 + t341 * t292 - t343 * t335) * t300) * t303 + (((t346 * t293 + t344 * t294 + t340) * t304 + (t347 * t293 + t345 * t294) * t303) * t325 + (t342 * t293 + t341 * t294 - t343 * t333) * t300) * t304) * t325 / 0.2e1;
T = t1;
