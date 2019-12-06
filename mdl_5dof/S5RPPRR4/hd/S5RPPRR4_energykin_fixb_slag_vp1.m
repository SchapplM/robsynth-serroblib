% Calculate kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:02
% EndTime: 2019-12-05 17:44:04
% DurationCPUTime: 1.29s
% Computational Cost: add. (950->192), mult. (1152->320), div. (0->0), fcn. (1161->10), ass. (0->98)
t348 = sin(pkin(8));
t350 = cos(pkin(8));
t346 = pkin(9) + qJ(4);
t340 = cos(t346);
t353 = cos(qJ(1));
t376 = t353 * t340;
t339 = sin(t346);
t352 = sin(qJ(1));
t384 = t352 * t339;
t319 = t350 * t384 + t376;
t377 = t353 * t339;
t383 = t352 * t340;
t320 = -t350 * t383 + t377;
t321 = -t350 * t377 + t383;
t322 = t350 * t376 + t384;
t374 = t353 * t348;
t381 = t352 * t348;
t356 = (Icges(5,5) * t320 + Icges(5,6) * t319 - Icges(5,3) * t381) * t352 - (Icges(5,5) * t322 + Icges(5,6) * t321 + Icges(5,3) * t374) * t353;
t393 = t348 * t356;
t370 = pkin(4) * t340;
t392 = pkin(7) * t348 + t350 * t370;
t349 = cos(pkin(9));
t391 = pkin(3) * t349 * t350 + pkin(6) * t348;
t390 = t350 ^ 2;
t341 = qJ(5) + t346;
t336 = sin(t341);
t386 = t352 * t336;
t337 = cos(t341);
t385 = t352 * t337;
t347 = sin(pkin(9));
t382 = t352 * t347;
t380 = t352 * t349;
t379 = t353 * t336;
t378 = t353 * t337;
t375 = t353 * t347;
t373 = t353 * t349;
t332 = t353 * pkin(1) + t352 * qJ(2);
t358 = pkin(2) * t350 + qJ(3) * t348;
t371 = -t358 * t353 - t332;
t369 = qJD(1) * (-t352 * pkin(1) + t353 * qJ(2)) + qJD(2) * t352;
t367 = qJD(3) * t350;
t366 = qJD(3) * t352;
t365 = qJD(4) * t348;
t364 = qJD(4) * t353;
t363 = qJD(4) + qJD(5);
t362 = t352 * t365;
t361 = pkin(4) * t339;
t360 = -qJD(1) * t358 * t352 + qJD(3) * t374 + t369;
t359 = -rSges(3,1) * t350 + rSges(3,2) * t348;
t357 = qJD(1) * (pkin(3) * t375 - t391 * t352) + t360;
t343 = qJD(2) * t353;
t355 = t343 + (-pkin(3) * t382 - t391 * t353 + t371) * qJD(1);
t335 = -qJD(4) * t350 + qJD(1);
t333 = t353 * rSges(2,1) - t352 * rSges(2,2);
t331 = -t352 * rSges(2,1) - t353 * rSges(2,2);
t328 = -t350 * t363 + qJD(1);
t325 = t363 * t374;
t324 = t363 * t381;
t318 = t350 * t378 + t386;
t317 = -t350 * t379 + t385;
t316 = -t350 * t385 + t379;
t315 = t350 * t386 + t378;
t314 = -t350 * rSges(5,3) + (rSges(5,1) * t340 - rSges(5,2) * t339) * t348;
t313 = -Icges(5,5) * t350 + (Icges(5,1) * t340 - Icges(5,4) * t339) * t348;
t312 = -Icges(5,6) * t350 + (Icges(5,4) * t340 - Icges(5,2) * t339) * t348;
t311 = -Icges(5,3) * t350 + (Icges(5,5) * t340 - Icges(5,6) * t339) * t348;
t310 = -t350 * rSges(6,3) + (rSges(6,1) * t337 - rSges(6,2) * t336) * t348;
t309 = -Icges(6,5) * t350 + (Icges(6,1) * t337 - Icges(6,4) * t336) * t348;
t308 = -Icges(6,6) * t350 + (Icges(6,4) * t337 - Icges(6,2) * t336) * t348;
t307 = -Icges(6,3) * t350 + (Icges(6,5) * t337 - Icges(6,6) * t336) * t348;
t306 = t343 + (-t352 * rSges(3,3) + t353 * t359 - t332) * qJD(1);
t305 = qJD(1) * (t353 * rSges(3,3) + t352 * t359) + t369;
t304 = -pkin(7) * t350 + t348 * t370;
t301 = t322 * rSges(5,1) + t321 * rSges(5,2) + rSges(5,3) * t374;
t300 = t320 * rSges(5,1) + t319 * rSges(5,2) - rSges(5,3) * t381;
t299 = Icges(5,1) * t322 + Icges(5,4) * t321 + Icges(5,5) * t374;
t298 = Icges(5,1) * t320 + Icges(5,4) * t319 - Icges(5,5) * t381;
t297 = Icges(5,4) * t322 + Icges(5,2) * t321 + Icges(5,6) * t374;
t296 = Icges(5,4) * t320 + Icges(5,2) * t319 - Icges(5,6) * t381;
t293 = t318 * rSges(6,1) + t317 * rSges(6,2) + rSges(6,3) * t374;
t292 = t316 * rSges(6,1) + t315 * rSges(6,2) - rSges(6,3) * t381;
t291 = Icges(6,1) * t318 + Icges(6,4) * t317 + Icges(6,5) * t374;
t290 = Icges(6,1) * t316 + Icges(6,4) * t315 - Icges(6,5) * t381;
t289 = Icges(6,4) * t318 + Icges(6,2) * t317 + Icges(6,6) * t374;
t288 = Icges(6,4) * t316 + Icges(6,2) * t315 - Icges(6,6) * t381;
t287 = Icges(6,5) * t318 + Icges(6,6) * t317 + Icges(6,3) * t374;
t286 = Icges(6,5) * t316 + Icges(6,6) * t315 - Icges(6,3) * t381;
t285 = t361 * t352 + t392 * t353;
t284 = -t392 * t352 + t361 * t353;
t283 = -t348 * t366 + t343 + (-(t350 * t373 + t382) * rSges(4,1) - (-t350 * t375 + t380) * rSges(4,2) - rSges(4,3) * t374 + t371) * qJD(1);
t282 = qJD(1) * ((-t350 * t380 + t375) * rSges(4,1) + (t350 * t382 + t373) * rSges(4,2) - rSges(4,3) * t381) + t360;
t281 = -t367 + (-t300 * t353 - t301 * t352) * t365;
t280 = -t335 * t301 + (t314 * t364 - t366) * t348 + t355;
t279 = t335 * t300 + t314 * t362 + t357;
t278 = -t367 - t325 * t292 - t324 * t293 + (-t284 * t353 - t285 * t352) * t365;
t277 = -t335 * t285 - t328 * t293 + t325 * t310 + (t304 * t364 - t366) * t348 + t355;
t276 = t335 * t284 + t328 * t292 + t304 * t362 + t324 * t310 + t357;
t1 = m(3) * (t305 ^ 2 + t306 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t390 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(5) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t335 * ((-t350 * t311 + (-t312 * t339 + t313 * t340) * t348) * t335 + ((-(-t296 * t339 + t298 * t340) * t352 + (-t297 * t339 + t299 * t340) * t353) * t348 + t356 * t350) * t365) / 0.2e1 - ((-t311 * t381 + t319 * t312 + t320 * t313) * t335 + ((t319 * t297 + t320 * t299) * t353 + (-t319 * t296 - t320 * t298 + t393) * t352) * t365) * t362 / 0.2e1 + t348 * ((t311 * t374 + t321 * t312 + t322 * t313) * t335 + (-(t321 * t296 + t322 * t298) * t352 + (t321 * t297 + t322 * t299 - t393) * t353) * t365) * t364 / 0.2e1 + m(6) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t328 * ((t286 * t324 - t287 * t325 - t307 * t328) * t350 + ((-t308 * t336 + t309 * t337) * t328 - (-t288 * t336 + t290 * t337) * t324 + (-t289 * t336 + t291 * t337) * t325) * t348) / 0.2e1 - t324 * ((-t307 * t381 + t315 * t308 + t316 * t309) * t328 - (-t286 * t381 + t315 * t288 + t316 * t290) * t324 + (-t287 * t381 + t315 * t289 + t316 * t291) * t325) / 0.2e1 + t325 * ((t307 * t374 + t317 * t308 + t318 * t309) * t328 - (t286 * t374 + t317 * t288 + t318 * t290) * t324 + (t287 * t374 + t317 * t289 + t318 * t291) * t325) / 0.2e1 + (m(2) * (t331 ^ 2 + t333 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t390 + ((Icges(3,1) + Icges(4,1) * t349 ^ 2 + (-0.2e1 * Icges(4,4) * t349 + Icges(4,2) * t347) * t347) * t348 + 0.2e1 * (-Icges(4,5) * t349 + Icges(4,6) * t347 + Icges(3,4)) * t350) * t348) * qJD(1) ^ 2 / 0.2e1;
T = t1;
