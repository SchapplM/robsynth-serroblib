% Calculate kinetic energy for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:30
% EndTime: 2019-12-05 15:08:31
% DurationCPUTime: 1.04s
% Computational Cost: add. (802->116), mult. (1085->184), div. (0->0), fcn. (1114->6), ass. (0->72)
t382 = Icges(5,1) + Icges(6,1);
t381 = Icges(5,4) - Icges(6,5);
t380 = Icges(6,4) + Icges(5,5);
t379 = Icges(5,2) + Icges(6,3);
t378 = Icges(6,6) - Icges(5,6);
t377 = -Icges(5,3) - Icges(6,2);
t328 = cos(pkin(7));
t361 = t328 ^ 2;
t327 = sin(pkin(7));
t362 = t327 ^ 2;
t363 = t361 + t362;
t376 = qJD(3) * t363;
t375 = rSges(6,1) + pkin(4);
t374 = rSges(6,3) + qJ(5);
t373 = t327 * t328;
t326 = pkin(8) + qJ(3);
t325 = cos(t326);
t330 = cos(qJ(4));
t352 = t328 * t330;
t329 = sin(qJ(4));
t355 = t327 * t329;
t314 = t325 * t355 + t352;
t353 = t328 * t329;
t354 = t327 * t330;
t315 = t325 * t354 - t353;
t324 = sin(t326);
t357 = t324 * t327;
t372 = t379 * t314 - t381 * t315 + t378 * t357;
t316 = t325 * t353 - t354;
t317 = t325 * t352 + t355;
t356 = t324 * t328;
t371 = t379 * t316 - t381 * t317 + t378 * t356;
t370 = t378 * t314 + t380 * t315 - t377 * t357;
t369 = t378 * t316 + t380 * t317 - t377 * t356;
t368 = -t381 * t314 + t382 * t315 + t380 * t357;
t367 = -t381 * t316 + t382 * t317 + t380 * t356;
t366 = t378 * t325 + (-t379 * t329 + t381 * t330) * t324;
t365 = t377 * t325 + (t378 * t329 + t380 * t330) * t324;
t364 = t380 * t325 + (t381 * t329 - t382 * t330) * t324;
t351 = rSges(6,2) * t357 + t374 * t314 + t375 * t315;
t350 = -rSges(6,2) * t356 - t374 * t316 - t375 * t317;
t349 = -t325 * rSges(6,2) + (t374 * t329 + t375 * t330) * t324;
t348 = qJD(2) * t328;
t347 = qJD(3) * t327;
t346 = qJD(3) * t328;
t345 = qJD(4) * t324;
t344 = qJD(4) * t325;
t343 = qJD(1) + (pkin(3) * t325 + pkin(6) * t324) * t376;
t321 = t324 * pkin(3) - t325 * pkin(6);
t323 = qJD(2) * t327;
t341 = -t321 * t346 + t323;
t338 = Icges(4,5) * t325 - Icges(4,6) * t324;
t334 = -t321 * t347 - t348;
t332 = qJD(1) ^ 2;
t320 = t324 * rSges(4,1) + t325 * rSges(4,2);
t319 = t327 * t345 - t346;
t318 = t328 * t345 + t347;
t310 = -t320 * t347 - t348;
t309 = -t320 * t346 + t323;
t304 = Icges(4,3) * t327 + t338 * t328;
t303 = -Icges(4,3) * t328 + t338 * t327;
t302 = -t325 * rSges(5,3) + (rSges(5,1) * t330 - rSges(5,2) * t329) * t324;
t292 = t317 * rSges(5,1) - t316 * rSges(5,2) + rSges(5,3) * t356;
t290 = t315 * rSges(5,1) - t314 * rSges(5,2) + rSges(5,3) * t357;
t276 = qJD(1) + (rSges(4,1) * t325 - rSges(4,2) * t324) * t376;
t275 = -t292 * t344 - t318 * t302 + t334;
t274 = t290 * t344 + t319 * t302 + t341;
t273 = t318 * t290 - t319 * t292 + t343;
t272 = qJD(5) * t314 - t349 * t318 + t350 * t344 + t334;
t271 = qJD(5) * t316 + t349 * t319 + t351 * t344 + t341;
t270 = qJD(5) * t324 * t329 + t351 * t318 + t350 * t319 + t343;
t1 = m(2) * t332 / 0.2e1 + m(3) * (t363 * qJD(2) ^ 2 + t332) / 0.2e1 + m(4) * (t276 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(5) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(6) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + (t327 * (-t303 * t373 + t362 * t304) / 0.2e1 - t328 * (t361 * t303 - t304 * t373) / 0.2e1) * qJD(3) ^ 2 + ((t366 * t316 + t364 * t317 - t365 * t356) * t344 + (t372 * t316 + t368 * t317 + t370 * t356) * t319 + (t371 * t316 + t367 * t317 + t369 * t356) * t318) * t318 / 0.2e1 + ((t366 * t314 + t364 * t315 - t365 * t357) * t344 + (t372 * t314 + t368 * t315 + t370 * t357) * t319 + (t371 * t314 + t367 * t315 + t369 * t357) * t318) * t319 / 0.2e1 - ((-t369 * t318 - t370 * t319 + t365 * t344) * t325 + ((t366 * t329 + t364 * t330) * t344 + (t372 * t329 + t368 * t330) * t319 + (t371 * t329 + t367 * t330) * t318) * t324) * t344 / 0.2e1;
T = t1;
