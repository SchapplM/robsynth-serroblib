% Calculate kinetic energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:35
% EndTime: 2019-12-05 15:12:36
% DurationCPUTime: 0.77s
% Computational Cost: add. (781->128), mult. (772->224), div. (0->0), fcn. (746->8), ass. (0->74)
t318 = sin(pkin(8));
t319 = cos(pkin(8));
t358 = t318 * t319;
t355 = t319 ^ 2;
t356 = t318 ^ 2;
t357 = t355 + t356;
t317 = pkin(9) + qJ(3);
t314 = qJ(4) + t317;
t306 = sin(t314);
t353 = t306 * t318;
t352 = t306 * t319;
t321 = sin(qJ(5));
t351 = t318 * t321;
t322 = cos(qJ(5));
t350 = t318 * t322;
t349 = t319 * t321;
t348 = t319 * t322;
t313 = cos(t317);
t347 = pkin(3) * t313;
t310 = qJD(3) * t318;
t304 = qJD(4) * t318 + t310;
t345 = qJD(2) * t319;
t344 = qJD(3) * t319;
t343 = qJD(5) * t306;
t307 = cos(t314);
t342 = qJD(5) * t307;
t312 = sin(t317);
t341 = pkin(3) * qJD(3) * t312;
t340 = qJD(1) + (-pkin(6) * t319 + t318 * t347) * t310 + (pkin(6) * t318 + t319 * t347) * t344;
t305 = (-qJD(3) - qJD(4)) * t319;
t338 = rSges(5,1) * t307 - rSges(5,2) * t306;
t336 = Icges(5,1) * t307 - Icges(5,4) * t306;
t334 = Icges(5,4) * t307 - Icges(5,2) * t306;
t333 = Icges(4,5) * t313 - Icges(4,6) * t312;
t332 = Icges(5,5) * t307 - Icges(5,6) * t306;
t331 = (-Icges(5,3) * t319 + t318 * t332) * t305 + (Icges(5,3) * t318 + t332 * t319) * t304;
t311 = qJD(2) * t318;
t328 = -t319 * t341 + t311;
t326 = -t318 * t341 - t345;
t325 = (-(Icges(5,6) * t318 + t319 * t334) * t306 + (Icges(5,5) * t318 + t319 * t336) * t307) * t304 + (-(-Icges(5,6) * t319 + t318 * t334) * t306 + (-Icges(5,5) * t319 + t318 * t336) * t307) * t305;
t324 = qJD(1) ^ 2;
t302 = rSges(4,1) * t312 + rSges(4,2) * t313;
t301 = pkin(4) * t306 - pkin(7) * t307;
t300 = rSges(5,1) * t306 + rSges(5,2) * t307;
t299 = t307 * t348 + t351;
t298 = -t307 * t349 + t350;
t297 = t307 * t350 - t349;
t296 = -t307 * t351 - t348;
t295 = t318 * t343 + t305;
t294 = t319 * t343 + t304;
t293 = -t302 * t310 - t345;
t292 = -t302 * t344 + t311;
t287 = Icges(4,3) * t318 + t319 * t333;
t286 = -Icges(4,3) * t319 + t318 * t333;
t279 = -t307 * rSges(6,3) + (rSges(6,1) * t322 - rSges(6,2) * t321) * t306;
t278 = -Icges(6,5) * t307 + (Icges(6,1) * t322 - Icges(6,4) * t321) * t306;
t277 = -Icges(6,6) * t307 + (Icges(6,4) * t322 - Icges(6,2) * t321) * t306;
t276 = -Icges(6,3) * t307 + (Icges(6,5) * t322 - Icges(6,6) * t321) * t306;
t275 = -t300 * t304 + t326;
t274 = t300 * t305 + t328;
t271 = rSges(6,1) * t299 + rSges(6,2) * t298 + rSges(6,3) * t352;
t270 = rSges(6,1) * t297 + rSges(6,2) * t296 + rSges(6,3) * t353;
t269 = Icges(6,1) * t299 + Icges(6,4) * t298 + Icges(6,5) * t352;
t268 = Icges(6,1) * t297 + Icges(6,4) * t296 + Icges(6,5) * t353;
t267 = Icges(6,4) * t299 + Icges(6,2) * t298 + Icges(6,6) * t352;
t266 = Icges(6,4) * t297 + Icges(6,2) * t296 + Icges(6,6) * t353;
t265 = Icges(6,5) * t299 + Icges(6,6) * t298 + Icges(6,3) * t352;
t264 = Icges(6,5) * t297 + Icges(6,6) * t296 + Icges(6,3) * t353;
t263 = qJD(1) + t357 * qJD(3) * (rSges(4,1) * t313 - rSges(4,2) * t312);
t262 = -t271 * t342 - t279 * t294 - t301 * t304 + t326;
t261 = t270 * t342 + t279 * t295 + t301 * t305 + t328;
t260 = t304 * (-t319 * rSges(5,3) + t318 * t338) - t305 * (t318 * rSges(5,3) + t338 * t319) + t340;
t259 = t294 * t270 - t295 * t271 + t340 + (t304 * t318 - t305 * t319) * (pkin(4) * t307 + pkin(7) * t306);
t1 = m(2) * t324 / 0.2e1 + m(3) * (qJD(2) ^ 2 * t357 + t324) / 0.2e1 + m(4) * (t263 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(5) * (t260 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t304 * (t318 * t331 + t325 * t319) / 0.2e1 + t305 * (t325 * t318 - t319 * t331) / 0.2e1 + m(6) * (t259 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t294 * ((t265 * t352 + t298 * t267 + t299 * t269) * t294 + (t264 * t352 + t266 * t298 + t268 * t299) * t295 - (t276 * t352 + t277 * t298 + t278 * t299) * t342) / 0.2e1 + t295 * ((t265 * t353 + t267 * t296 + t269 * t297) * t294 + (t264 * t353 + t296 * t266 + t297 * t268) * t295 - (t276 * t353 + t277 * t296 + t278 * t297) * t342) / 0.2e1 - ((-t264 * t295 - t265 * t294 + t276 * t342) * t307 + ((-t267 * t321 + t269 * t322) * t294 + (-t266 * t321 + t268 * t322) * t295 - (-t277 * t321 + t278 * t322) * t342) * t306) * t342 / 0.2e1 + (t318 * (-t286 * t358 + t356 * t287) / 0.2e1 - t319 * (t355 * t286 - t287 * t358) / 0.2e1) * qJD(3) ^ 2;
T = t1;
