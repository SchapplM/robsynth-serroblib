% Calculate kinetic energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:14
% EndTime: 2019-12-31 18:56:16
% DurationCPUTime: 1.30s
% Computational Cost: add. (538->172), mult. (1241->269), div. (0->0), fcn. (1223->6), ass. (0->94)
t389 = Icges(5,1) + Icges(6,1);
t388 = Icges(5,4) + Icges(6,4);
t387 = Icges(6,5) + Icges(5,5);
t386 = Icges(5,2) + Icges(6,2);
t385 = Icges(6,6) + Icges(5,6);
t384 = Icges(6,3) + Icges(5,3);
t333 = sin(qJ(3));
t335 = cos(qJ(4));
t337 = cos(qJ(1));
t361 = t335 * t337;
t332 = sin(qJ(4));
t334 = sin(qJ(1));
t365 = t332 * t334;
t309 = -t333 * t365 + t361;
t363 = t334 * t335;
t364 = t332 * t337;
t310 = t333 * t363 + t364;
t336 = cos(qJ(3));
t362 = t334 * t336;
t383 = t385 * t309 + t387 * t310 - t384 * t362;
t311 = t333 * t364 + t363;
t312 = -t333 * t361 + t365;
t360 = t336 * t337;
t382 = t385 * t311 + t387 * t312 + t384 * t360;
t381 = t386 * t309 + t388 * t310 - t385 * t362;
t380 = t386 * t311 + t388 * t312 + t385 * t360;
t379 = t388 * t309 + t389 * t310 - t387 * t362;
t378 = t388 * t311 + t389 * t312 + t387 * t360;
t377 = (-t385 * t332 + t387 * t335) * t336 + t384 * t333;
t376 = (-t386 * t332 + t388 * t335) * t336 + t385 * t333;
t375 = (-t388 * t332 + t389 * t335) * t336 + t387 * t333;
t369 = pkin(4) * t335;
t374 = -qJ(5) * t336 + t333 * t369;
t367 = Icges(4,4) * t333;
t366 = Icges(4,4) * t336;
t359 = rSges(6,1) * t310 + rSges(6,2) * t309 - rSges(6,3) * t362 + pkin(4) * t364 + t374 * t334;
t358 = rSges(6,1) * t312 + rSges(6,2) * t311 + rSges(6,3) * t360 + pkin(4) * t365 - t374 * t337;
t357 = (rSges(6,1) * t335 - rSges(6,2) * t332 + t369) * t336 + (qJ(5) + rSges(6,3)) * t333;
t318 = qJD(1) * (pkin(1) * t337 + qJ(2) * t334);
t356 = qJD(1) * t337 * pkin(6) + t318;
t355 = qJD(3) * t334;
t354 = qJD(3) * t337;
t353 = qJD(4) * t336;
t352 = qJD(5) * t336;
t348 = pkin(3) * t333 - pkin(7) * t336;
t313 = t348 * t334;
t351 = qJD(1) * t313 + t356;
t322 = pkin(1) * t334 - qJ(2) * t337;
t350 = -pkin(6) * t334 - t322;
t326 = pkin(3) * t336 + pkin(7) * t333;
t349 = -qJD(3) * t326 - qJD(2);
t314 = t348 * t337;
t347 = -t313 * t355 - t314 * t354;
t346 = rSges(4,1) * t333 + rSges(4,2) * t336;
t345 = Icges(4,1) * t333 + t366;
t344 = Icges(4,2) * t336 + t367;
t343 = Icges(4,5) * t333 + Icges(4,6) * t336;
t297 = Icges(4,6) * t337 + t334 * t344;
t301 = Icges(4,5) * t337 + t334 * t345;
t342 = -t297 * t336 - t301 * t333;
t298 = Icges(4,6) * t334 - t337 * t344;
t302 = Icges(4,5) * t334 - t337 * t345;
t341 = t298 * t336 + t302 * t333;
t320 = -Icges(4,2) * t333 + t366;
t321 = Icges(4,1) * t336 - t367;
t340 = t320 * t336 + t321 * t333;
t330 = qJD(2) * t334;
t339 = t326 * t355 + t330 + (t314 + t350) * qJD(1);
t327 = qJD(4) * t333 + qJD(1);
t325 = rSges(2,1) * t337 - rSges(2,2) * t334;
t324 = rSges(4,1) * t336 - rSges(4,2) * t333;
t323 = rSges(2,1) * t334 + rSges(2,2) * t337;
t319 = Icges(4,5) * t336 - Icges(4,6) * t333;
t317 = -t334 * t353 + t354;
t316 = t337 * t353 + t355;
t306 = rSges(4,3) * t334 - t337 * t346;
t305 = rSges(5,3) * t333 + (rSges(5,1) * t335 - rSges(5,2) * t332) * t336;
t303 = rSges(4,3) * t337 + t334 * t346;
t294 = Icges(4,3) * t334 - t337 * t343;
t293 = Icges(4,3) * t337 + t334 * t343;
t289 = t318 - qJD(2) * t337 + qJD(1) * (-rSges(3,2) * t337 + rSges(3,3) * t334);
t288 = t330 + (rSges(3,2) * t334 + rSges(3,3) * t337 - t322) * qJD(1);
t287 = rSges(5,1) * t312 + rSges(5,2) * t311 + rSges(5,3) * t360;
t285 = rSges(5,1) * t310 + rSges(5,2) * t309 - rSges(5,3) * t362;
t269 = (-t303 * t334 + t306 * t337) * qJD(3);
t268 = qJD(1) * t303 + (-qJD(3) * t324 - qJD(2)) * t337 + t356;
t267 = t324 * t355 + t330 + (-t306 + t350) * qJD(1);
t266 = t285 * t327 - t305 * t317 + t337 * t349 + t351;
t265 = -t287 * t327 + t305 * t316 + t339;
t264 = -t285 * t316 + t287 * t317 + t347;
t263 = t359 * t327 - t357 * t317 + (t349 + t352) * t337 + t351;
t262 = t316 * t357 - t327 * t358 - t334 * t352 + t339;
t261 = qJD(5) * t333 - t316 * t359 + t317 * t358 + t347;
t1 = m(3) * (t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(4) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + ((t337 * t319 + t334 * t340) * qJD(1) + (t337 ^ 2 * t293 + (t341 * t334 + (t294 - t342) * t337) * t334) * qJD(3)) * t354 / 0.2e1 + ((t334 * t319 - t337 * t340) * qJD(1) + (t334 ^ 2 * t294 + (t342 * t337 + (t293 - t341) * t334) * t337) * qJD(3)) * t355 / 0.2e1 + qJD(1) * ((-t333 * t320 + t336 * t321) * qJD(1) + ((-t297 * t333 + t301 * t336) * t337 + (-t298 * t333 + t302 * t336) * t334) * qJD(3)) / 0.2e1 + m(5) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(6) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + ((t376 * t311 + t375 * t312 + t377 * t360) * t327 + (t381 * t311 + t379 * t312 + t383 * t360) * t317 + (t380 * t311 + t378 * t312 + t382 * t360) * t316) * t316 / 0.2e1 + ((t376 * t309 + t375 * t310 - t377 * t362) * t327 + (t381 * t309 + t379 * t310 - t383 * t362) * t317 + (t380 * t309 + t378 * t310 - t382 * t362) * t316) * t317 / 0.2e1 + (((-t376 * t332 + t375 * t335) * t327 + (-t381 * t332 + t379 * t335) * t317 + (-t380 * t332 + t378 * t335) * t316) * t336 + (t382 * t316 + t383 * t317 + t377 * t327) * t333) * t327 / 0.2e1 + (m(2) * (t323 ^ 2 + t325 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
