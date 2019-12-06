% Calculate kinetic energy for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:37
% EndTime: 2019-12-05 15:06:38
% DurationCPUTime: 1.06s
% Computational Cost: add. (832->122), mult. (1094->192), div. (0->0), fcn. (1113->6), ass. (0->71)
t384 = Icges(5,1) + Icges(6,1);
t383 = Icges(5,4) + Icges(6,4);
t382 = Icges(5,5) + Icges(6,5);
t381 = Icges(5,2) + Icges(6,2);
t380 = -Icges(6,6) - Icges(5,6);
t379 = -Icges(6,3) - Icges(5,3);
t329 = cos(pkin(7));
t365 = t329 ^ 2;
t328 = sin(pkin(7));
t366 = t328 ^ 2;
t367 = t365 + t366;
t378 = qJD(3) * t367;
t377 = t328 * t329;
t327 = pkin(8) + qJ(3);
t326 = cos(t327);
t332 = cos(qJ(4));
t354 = t329 * t332;
t331 = sin(qJ(4));
t357 = t328 * t331;
t315 = -t326 * t357 - t354;
t355 = t329 * t331;
t356 = t328 * t332;
t316 = t326 * t356 - t355;
t325 = sin(t327);
t359 = t325 * t328;
t376 = -t380 * t315 + t382 * t316 - t379 * t359;
t317 = -t326 * t355 + t356;
t318 = t326 * t354 + t357;
t358 = t325 * t329;
t375 = -t380 * t317 + t382 * t318 - t379 * t358;
t374 = t381 * t315 + t383 * t316 - t380 * t359;
t373 = t381 * t317 + t383 * t318 - t380 * t358;
t372 = t383 * t315 + t384 * t316 + t382 * t359;
t371 = t383 * t317 + t384 * t318 + t382 * t358;
t370 = t379 * t326 + (t380 * t331 + t382 * t332) * t325;
t369 = t380 * t326 + (-t381 * t331 + t383 * t332) * t325;
t368 = t382 * t326 + (t383 * t331 - t384 * t332) * t325;
t361 = pkin(4) * t332;
t335 = qJ(5) * t325 + t361 * t326;
t353 = rSges(6,1) * t316 + rSges(6,2) * t315 + rSges(6,3) * t359 - pkin(4) * t355 + t335 * t328;
t352 = -rSges(6,1) * t318 - rSges(6,2) * t317 - rSges(6,3) * t358 - pkin(4) * t357 - t335 * t329;
t351 = (-qJ(5) - rSges(6,3)) * t326 + (rSges(6,1) * t332 - rSges(6,2) * t331 + t361) * t325;
t350 = qJD(2) * t329;
t349 = qJD(3) * t328;
t348 = qJD(3) * t329;
t347 = qJD(4) * t325;
t346 = qJD(4) * t326;
t345 = qJD(1) + (pkin(3) * t326 + pkin(6) * t325) * t378;
t341 = Icges(4,5) * t326 - Icges(4,6) * t325;
t322 = pkin(3) * t325 - pkin(6) * t326;
t338 = -qJD(3) * t322 + qJD(5) * t325;
t334 = qJD(1) ^ 2;
t324 = qJD(2) * t328;
t321 = rSges(4,1) * t325 + rSges(4,2) * t326;
t320 = t328 * t347 - t348;
t319 = t329 * t347 + t349;
t312 = -t321 * t349 - t350;
t311 = -t321 * t348 + t324;
t306 = Icges(4,3) * t328 + t341 * t329;
t305 = -Icges(4,3) * t329 + t341 * t328;
t304 = -rSges(5,3) * t326 + (rSges(5,1) * t332 - rSges(5,2) * t331) * t325;
t295 = rSges(5,1) * t318 + rSges(5,2) * t317 + rSges(5,3) * t358;
t293 = rSges(5,1) * t316 + rSges(5,2) * t315 + rSges(5,3) * t359;
t277 = qJD(1) + (rSges(4,1) * t326 - rSges(4,2) * t325) * t378;
t276 = -t295 * t346 - t304 * t319 - t322 * t349 - t350;
t275 = t293 * t346 + t304 * t320 - t322 * t348 + t324;
t274 = t293 * t319 - t295 * t320 + t345;
t273 = -t351 * t319 + t338 * t328 + t352 * t346 - t350;
t272 = t351 * t320 + t338 * t329 + t353 * t346 + t324;
t271 = -qJD(5) * t326 + t353 * t319 + t352 * t320 + t345;
t1 = m(2) * t334 / 0.2e1 + m(3) * (t367 * qJD(2) ^ 2 + t334) / 0.2e1 + m(4) * (t277 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(5) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(6) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + (t328 * (-t305 * t377 + t366 * t306) / 0.2e1 - t329 * (t365 * t305 - t306 * t377) / 0.2e1) * qJD(3) ^ 2 + ((-t369 * t317 + t368 * t318 - t370 * t358) * t346 + (t374 * t317 + t372 * t318 + t376 * t358) * t320 + (t373 * t317 + t371 * t318 + t375 * t358) * t319) * t319 / 0.2e1 + ((-t369 * t315 + t368 * t316 - t370 * t359) * t346 + (t374 * t315 + t372 * t316 + t376 * t359) * t320 + (t373 * t315 + t371 * t316 + t375 * t359) * t319) * t320 / 0.2e1 - ((-t375 * t319 - t376 * t320 + t370 * t346) * t326 + ((t369 * t331 + t368 * t332) * t346 + (-t374 * t331 + t372 * t332) * t320 + (-t373 * t331 + t371 * t332) * t319) * t325) * t346 / 0.2e1;
T = t1;
