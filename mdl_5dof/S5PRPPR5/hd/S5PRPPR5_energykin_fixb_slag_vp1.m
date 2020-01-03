% Calculate kinetic energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:01
% EndTime: 2019-12-31 17:38:02
% DurationCPUTime: 1.09s
% Computational Cost: add. (558->147), mult. (1439->243), div. (0->0), fcn. (1579->8), ass. (0->75)
t389 = Icges(4,2) + Icges(3,3);
t344 = sin(qJ(2));
t346 = cos(qJ(2));
t388 = (Icges(4,4) + Icges(3,5)) * t346 + (-Icges(3,6) + Icges(4,6)) * t344;
t342 = cos(pkin(7));
t380 = t342 ^ 2;
t341 = sin(pkin(7));
t381 = t341 ^ 2;
t383 = t380 + t381;
t382 = qJD(2) * t383;
t387 = t388 * t341 - t389 * t342;
t386 = t389 * t341 + t388 * t342;
t340 = sin(pkin(8));
t377 = cos(pkin(8));
t330 = -t346 * t340 + t344 * t377;
t379 = qJD(2) ^ 2;
t378 = pkin(3) * t346;
t372 = qJD(3) * t344;
t337 = t341 * t372;
t375 = qJD(4) * t342 + t337;
t374 = qJD(2) * t341;
t373 = qJD(2) * t342;
t329 = t344 * t340 + t346 * t377;
t371 = qJD(5) * t329;
t332 = pkin(2) * t344 - qJ(3) * t346;
t368 = -pkin(3) * t344 - t332;
t366 = qJD(2) * (-rSges(4,1) * t344 + rSges(4,3) * t346 - t332);
t338 = t342 * t372;
t365 = -qJD(4) * t341 + t338;
t354 = -qJD(3) * t346 + qJD(1) + (pkin(2) * t346 + qJ(3) * t344) * t382;
t353 = qJD(2) * (-rSges(5,1) * t330 + rSges(5,2) * t329 + t368);
t352 = qJD(2) * (-pkin(4) * t330 - pkin(6) * t329 + t368);
t348 = (qJ(4) * t342 + t341 * t378) * t374 + (-qJ(4) * t341 + t342 * t378) * t373 + t354;
t345 = cos(qJ(5));
t343 = sin(qJ(5));
t334 = rSges(3,1) * t344 + rSges(3,2) * t346;
t326 = t329 * t342;
t325 = t330 * t342;
t324 = t329 * t341;
t323 = t330 * t341;
t308 = -t323 * qJD(5) - t373;
t307 = -t325 * qJD(5) + t374;
t306 = t326 * t345 - t341 * t343;
t305 = -t326 * t343 - t341 * t345;
t304 = t324 * t345 + t342 * t343;
t303 = -t324 * t343 + t342 * t345;
t300 = t342 * t366 + t338;
t299 = t341 * t366 + t337;
t298 = Icges(5,1) * t326 + Icges(5,4) * t325 - Icges(5,5) * t341;
t297 = Icges(5,1) * t324 + Icges(5,4) * t323 + Icges(5,5) * t342;
t296 = Icges(5,4) * t326 + Icges(5,2) * t325 - Icges(5,6) * t341;
t295 = Icges(5,4) * t324 + Icges(5,2) * t323 + Icges(5,6) * t342;
t294 = Icges(5,5) * t326 + Icges(5,6) * t325 - Icges(5,3) * t341;
t293 = Icges(5,5) * t324 + Icges(5,6) * t323 + Icges(5,3) * t342;
t292 = rSges(6,3) * t329 + (rSges(6,1) * t345 - rSges(6,2) * t343) * t330;
t291 = Icges(6,5) * t329 + (Icges(6,1) * t345 - Icges(6,4) * t343) * t330;
t290 = Icges(6,6) * t329 + (Icges(6,4) * t345 - Icges(6,2) * t343) * t330;
t289 = Icges(6,3) * t329 + (Icges(6,5) * t345 - Icges(6,6) * t343) * t330;
t288 = qJD(1) + (rSges(3,1) * t346 - rSges(3,2) * t344) * t382;
t287 = t342 * t353 + t365;
t286 = t341 * t353 + t375;
t285 = rSges(6,1) * t306 + rSges(6,2) * t305 - rSges(6,3) * t325;
t284 = rSges(6,1) * t304 + rSges(6,2) * t303 - rSges(6,3) * t323;
t283 = Icges(6,1) * t306 + Icges(6,4) * t305 - Icges(6,5) * t325;
t282 = Icges(6,1) * t304 + Icges(6,4) * t303 - Icges(6,5) * t323;
t281 = Icges(6,4) * t306 + Icges(6,2) * t305 - Icges(6,6) * t325;
t280 = Icges(6,4) * t304 + Icges(6,2) * t303 - Icges(6,6) * t323;
t279 = Icges(6,5) * t306 + Icges(6,6) * t305 - Icges(6,3) * t325;
t278 = Icges(6,5) * t304 + Icges(6,6) * t303 - Icges(6,3) * t323;
t277 = t354 + (rSges(4,1) * t346 + rSges(4,3) * t344) * t382;
t276 = (t341 * (rSges(5,1) * t324 + rSges(5,2) * t323) + t342 * (rSges(5,1) * t326 + rSges(5,2) * t325)) * qJD(2) + t348;
t275 = -t284 * t371 + t292 * t308 + t342 * t352 + t365;
t274 = t285 * t371 - t292 * t307 + t341 * t352 + t375;
t273 = t307 * t284 - t308 * t285 + (t341 * (pkin(4) * t324 - pkin(6) * t323) + t342 * (pkin(4) * t326 - pkin(6) * t325)) * qJD(2) + t348;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t383 * t379 * t334 ^ 2 + t288 ^ 2) / 0.2e1 + m(4) * (t277 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + m(5) * (t276 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(6) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t307 * ((-t279 * t325 + t281 * t305 + t283 * t306) * t307 + (-t278 * t325 + t280 * t305 + t282 * t306) * t308 + (-t289 * t325 + t290 * t305 + t291 * t306) * t371) / 0.2e1 + t308 * ((-t279 * t323 + t281 * t303 + t283 * t304) * t307 + (-t278 * t323 + t280 * t303 + t282 * t304) * t308 + (-t289 * t323 + t290 * t303 + t291 * t304) * t371) / 0.2e1 + ((t278 * t308 + t279 * t307 + t289 * t371) * t329 + ((-t281 * t343 + t283 * t345) * t307 + (-t280 * t343 + t282 * t345) * t308 + (-t290 * t343 + t291 * t345) * t371) * t330) * t371 / 0.2e1 + ((-t294 * t341 + t296 * t325 + t298 * t326) * t341 + t386 * t381 + (-t295 * t325 - t297 * t326 + (t293 - t387) * t341) * t342) * t341 * t379 / 0.2e1 - (-(t293 * t342 + t295 * t323 + t297 * t324) * t342 + t387 * t380 + (t296 * t323 + t298 * t324 + (t294 - t386) * t342) * t341) * t342 * t379 / 0.2e1;
T = t1;
