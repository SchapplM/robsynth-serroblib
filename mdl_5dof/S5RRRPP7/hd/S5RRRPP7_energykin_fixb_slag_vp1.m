% Calculate kinetic energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:55
% EndTime: 2019-12-31 21:03:57
% DurationCPUTime: 1.58s
% Computational Cost: add. (742->184), mult. (1772->283), div. (0->0), fcn. (1805->6), ass. (0->100)
t393 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t392 = -Icges(4,4) + Icges(6,4) + Icges(5,5);
t391 = Icges(6,5) - Icges(5,4) - Icges(4,5);
t390 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t389 = -Icges(5,6) + Icges(6,6) + Icges(4,6);
t388 = Icges(6,3) + Icges(4,3) + Icges(5,2);
t387 = rSges(6,1) + pkin(4);
t386 = rSges(6,3) + qJ(5);
t337 = sin(qJ(3));
t339 = sin(qJ(1));
t340 = cos(qJ(3));
t341 = cos(qJ(2));
t342 = cos(qJ(1));
t365 = t341 * t342;
t317 = t337 * t365 - t339 * t340;
t318 = t337 * t339 + t340 * t365;
t286 = pkin(3) * t318 + qJ(4) * t317;
t366 = t339 * t341;
t315 = t337 * t366 + t340 * t342;
t335 = -qJD(3) * t341 + qJD(1);
t385 = qJD(4) * t315 + t335 * t286;
t338 = sin(qJ(2));
t319 = (pkin(3) * t340 + qJ(4) * t337) * t338;
t357 = qJD(3) * t338;
t358 = qJD(2) * t342;
t323 = t339 * t357 - t358;
t384 = qJD(4) * t317 + t323 * t319;
t316 = -t337 * t342 + t340 * t366;
t368 = t338 * t339;
t383 = t389 * t315 + t391 * t316 - t388 * t368;
t367 = t338 * t342;
t382 = t389 * t317 + t391 * t318 - t388 * t367;
t381 = t390 * t315 + t392 * t316 - t389 * t368;
t380 = t390 * t317 + t392 * t318 - t389 * t367;
t379 = t392 * t315 + t393 * t316 - t391 * t368;
t378 = t392 * t317 + t393 * t318 - t391 * t367;
t377 = t388 * t341 + (t389 * t337 + t391 * t340) * t338;
t376 = t389 * t341 + (t390 * t337 + t392 * t340) * t338;
t375 = t391 * t341 + (t392 * t337 + t393 * t340) * t338;
t370 = Icges(3,4) * t338;
t369 = Icges(3,4) * t341;
t364 = rSges(6,2) * t315 + t387 * t316 - t386 * t368;
t363 = rSges(6,2) * t317 + t387 * t318 - t386 * t367;
t362 = t386 * t341 + (rSges(6,2) * t337 + t387 * t340) * t338;
t354 = pkin(2) * t341 + pkin(7) * t338;
t320 = t354 * t339;
t321 = t354 * t342;
t359 = qJD(2) * t339;
t361 = t320 * t359 + t321 * t358;
t325 = qJD(1) * (pkin(1) * t342 + pkin(6) * t339);
t360 = qJD(1) * t321 + t325;
t333 = pkin(1) * t339 - pkin(6) * t342;
t356 = (-t320 - t333) * qJD(1);
t285 = pkin(3) * t316 + qJ(4) * t315;
t322 = t342 * t357 + t359;
t355 = qJD(4) * t338 * t337 + t322 * t285 + t361;
t353 = rSges(3,1) * t341 - rSges(3,2) * t338;
t352 = Icges(3,1) * t341 - t370;
t351 = -Icges(3,2) * t338 + t369;
t350 = Icges(3,5) * t341 - Icges(3,6) * t338;
t298 = -Icges(3,6) * t342 + t351 * t339;
t303 = -Icges(3,5) * t342 + t352 * t339;
t349 = t298 * t338 - t303 * t341;
t299 = Icges(3,6) * t339 + t351 * t342;
t304 = Icges(3,5) * t339 + t352 * t342;
t348 = -t299 * t338 + t304 * t341;
t327 = Icges(3,2) * t341 + t370;
t328 = Icges(3,1) * t338 + t369;
t347 = -t327 * t338 + t328 * t341;
t332 = pkin(2) * t338 - pkin(7) * t341;
t346 = -qJD(2) * t332 - qJD(5) * t338;
t345 = -t332 * t359 + t360;
t344 = -t332 * t358 + t356;
t331 = rSges(2,1) * t342 - rSges(2,2) * t339;
t330 = rSges(2,1) * t339 + rSges(2,2) * t342;
t329 = rSges(3,1) * t338 + rSges(3,2) * t341;
t326 = Icges(3,5) * t338 + Icges(3,6) * t341;
t309 = rSges(3,3) * t339 + t353 * t342;
t308 = -rSges(3,3) * t342 + t353 * t339;
t307 = -rSges(4,3) * t341 + (rSges(4,1) * t340 - rSges(4,2) * t337) * t338;
t306 = -rSges(5,2) * t341 + (rSges(5,1) * t340 + rSges(5,3) * t337) * t338;
t294 = Icges(3,3) * t339 + t350 * t342;
t293 = -Icges(3,3) * t342 + t350 * t339;
t284 = rSges(4,1) * t318 - rSges(4,2) * t317 + rSges(4,3) * t367;
t283 = rSges(5,1) * t318 + rSges(5,2) * t367 + rSges(5,3) * t317;
t281 = rSges(4,1) * t316 - rSges(4,2) * t315 + rSges(4,3) * t368;
t280 = rSges(5,1) * t316 + rSges(5,2) * t368 + rSges(5,3) * t315;
t259 = qJD(1) * t309 - t329 * t359 + t325;
t258 = -t329 * t358 + (-t308 - t333) * qJD(1);
t256 = (t308 * t339 + t309 * t342) * qJD(2);
t255 = t284 * t335 - t307 * t322 + t345;
t254 = -t281 * t335 + t307 * t323 + t344;
t253 = t281 * t322 - t284 * t323 + t361;
t252 = t283 * t335 + (-t306 - t319) * t322 + t345 + t385;
t251 = t306 * t323 + (-t280 - t285) * t335 + t344 + t384;
t250 = t280 * t322 + (-t283 - t286) * t323 + t355;
t249 = t346 * t339 + t363 * t335 + (-t319 - t362) * t322 + t360 + t385;
t248 = t346 * t342 + t362 * t323 + t356 + (-t285 - t364) * t335 + t384;
t247 = qJD(5) * t341 + t364 * t322 + (-t286 - t363) * t323 + t355;
t1 = m(3) * (t256 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + ((t339 * t326 + t347 * t342) * qJD(1) + (t339 ^ 2 * t294 + (t349 * t342 + (-t293 + t348) * t339) * t342) * qJD(2)) * t359 / 0.2e1 - ((-t342 * t326 + t347 * t339) * qJD(1) + (t342 ^ 2 * t293 + (t348 * t339 + (-t294 + t349) * t342) * t339) * qJD(2)) * t358 / 0.2e1 + qJD(1) * ((t341 * t327 + t338 * t328) * qJD(1) + ((t299 * t341 + t304 * t338) * t339 - (t298 * t341 + t303 * t338) * t342) * qJD(2)) / 0.2e1 + m(4) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + (m(2) * (t330 ^ 2 + t331 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t317 * t376 + t318 * t375 - t367 * t377) * t335 + (t317 * t381 + t318 * t379 - t367 * t383) * t323 + (t380 * t317 + t378 * t318 - t382 * t367) * t322) * t322 / 0.2e1 + ((t315 * t376 + t316 * t375 - t368 * t377) * t335 + (t381 * t315 + t379 * t316 - t383 * t368) * t323 + (t315 * t380 + t316 * t378 - t368 * t382) * t322) * t323 / 0.2e1 + ((t322 * t382 + t323 * t383 + t335 * t377) * t341 + ((t337 * t376 + t340 * t375) * t335 + (t337 * t381 + t340 * t379) * t323 + (t337 * t380 + t340 * t378) * t322) * t338) * t335 / 0.2e1;
T = t1;
