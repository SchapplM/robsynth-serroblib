% Calculate kinetic energy for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:08
% EndTime: 2019-12-31 18:42:09
% DurationCPUTime: 1.26s
% Computational Cost: add. (1020->170), mult. (1232->268), div. (0->0), fcn. (1213->8), ass. (0->95)
t384 = Icges(5,1) + Icges(6,1);
t383 = Icges(5,4) + Icges(6,4);
t382 = -Icges(6,5) - Icges(5,5);
t381 = Icges(5,2) + Icges(6,2);
t380 = -Icges(6,6) - Icges(5,6);
t379 = -Icges(6,3) - Icges(5,3);
t325 = qJ(1) + pkin(8);
t324 = cos(t325);
t327 = sin(qJ(4));
t330 = cos(qJ(4));
t323 = sin(t325);
t331 = cos(qJ(3));
t356 = t331 * t323;
t295 = -t324 * t330 - t327 * t356;
t358 = t324 * t327;
t296 = t330 * t356 - t358;
t328 = sin(qJ(3));
t359 = t323 * t328;
t378 = -t380 * t295 - t382 * t296 - t379 * t359;
t355 = t331 * t324;
t297 = t323 * t330 - t327 * t355;
t360 = t323 * t327;
t298 = t330 * t355 + t360;
t357 = t324 * t328;
t377 = -t380 * t297 - t382 * t298 - t379 * t357;
t376 = t381 * t295 + t383 * t296 - t380 * t359;
t375 = t381 * t297 + t383 * t298 - t380 * t357;
t374 = t383 * t295 + t384 * t296 - t382 * t359;
t373 = t383 * t297 + t384 * t298 - t382 * t357;
t372 = t379 * t331 + (t380 * t327 - t382 * t330) * t328;
t371 = t380 * t331 + (-t381 * t327 + t383 * t330) * t328;
t370 = t382 * t331 + (-t383 * t327 + t384 * t330) * t328;
t329 = sin(qJ(1));
t365 = pkin(1) * t329;
t364 = pkin(4) * t330;
t362 = Icges(4,4) * t328;
t361 = Icges(4,4) * t331;
t334 = qJ(5) * t328 + t364 * t331;
t354 = rSges(6,1) * t296 + rSges(6,2) * t295 + rSges(6,3) * t359 - pkin(4) * t358 + t334 * t323;
t353 = rSges(6,1) * t298 + rSges(6,2) * t297 + rSges(6,3) * t357 + pkin(4) * t360 + t334 * t324;
t352 = (-qJ(5) - rSges(6,3)) * t331 + (rSges(6,1) * t330 - rSges(6,2) * t327 + t364) * t328;
t332 = cos(qJ(1));
t322 = qJD(1) * t332 * pkin(1);
t351 = qJD(1) * (pkin(2) * t324 + pkin(6) * t323) + t322;
t350 = qJD(3) * t323;
t349 = qJD(3) * t324;
t348 = qJD(4) * t328;
t344 = pkin(3) * t331 + pkin(7) * t328;
t308 = t344 * t324;
t347 = qJD(1) * t308 + t351;
t307 = t344 * t323;
t346 = t307 * t350 + t308 * t349 + qJD(2);
t345 = -pkin(2) * t323 + pkin(6) * t324 - t365;
t343 = rSges(4,1) * t331 - rSges(4,2) * t328;
t342 = Icges(4,1) * t331 - t362;
t341 = -Icges(4,2) * t328 + t361;
t340 = Icges(4,5) * t331 - Icges(4,6) * t328;
t283 = -Icges(4,6) * t324 + t341 * t323;
t285 = -Icges(4,5) * t324 + t342 * t323;
t339 = t283 * t328 - t285 * t331;
t284 = Icges(4,6) * t323 + t341 * t324;
t286 = Icges(4,5) * t323 + t342 * t324;
t338 = -t284 * t328 + t286 * t331;
t314 = Icges(4,2) * t331 + t362;
t315 = Icges(4,1) * t328 + t361;
t337 = -t314 * t328 + t315 * t331;
t319 = pkin(3) * t328 - pkin(7) * t331;
t336 = -qJD(3) * t319 + qJD(5) * t328;
t335 = (-t307 + t345) * qJD(1);
t320 = -qJD(4) * t331 + qJD(1);
t318 = rSges(2,1) * t332 - rSges(2,2) * t329;
t317 = rSges(2,1) * t329 + rSges(2,2) * t332;
t316 = rSges(4,1) * t328 + rSges(4,2) * t331;
t313 = Icges(4,5) * t328 + Icges(4,6) * t331;
t310 = t323 * t348 - t349;
t309 = t324 * t348 + t350;
t306 = -rSges(5,3) * t331 + (rSges(5,1) * t330 - rSges(5,2) * t327) * t328;
t293 = t322 + qJD(1) * (rSges(3,1) * t324 - rSges(3,2) * t323);
t292 = (-rSges(3,1) * t323 - rSges(3,2) * t324 - t365) * qJD(1);
t288 = rSges(4,3) * t323 + t343 * t324;
t287 = -rSges(4,3) * t324 + t343 * t323;
t282 = Icges(4,3) * t323 + t340 * t324;
t281 = -Icges(4,3) * t324 + t340 * t323;
t280 = rSges(5,1) * t298 + rSges(5,2) * t297 + rSges(5,3) * t357;
t278 = rSges(5,1) * t296 + rSges(5,2) * t295 + rSges(5,3) * t359;
t262 = qJD(1) * t288 - t316 * t350 + t351;
t261 = -t316 * t349 + (-t287 + t345) * qJD(1);
t260 = qJD(2) + (t287 * t323 + t288 * t324) * qJD(3);
t259 = t280 * t320 - t306 * t309 - t319 * t350 + t347;
t258 = -t278 * t320 + t306 * t310 - t319 * t349 + t335;
t257 = t278 * t309 - t280 * t310 + t346;
t256 = -t352 * t309 + t353 * t320 + t336 * t323 + t347;
t255 = t352 * t310 - t354 * t320 + t336 * t324 + t335;
t254 = -qJD(5) * t331 + t354 * t309 - t353 * t310 + t346;
t1 = m(3) * (qJD(2) ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(4) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + ((t323 * t313 + t337 * t324) * qJD(1) + (t323 ^ 2 * t282 + (t339 * t324 + (-t281 + t338) * t323) * t324) * qJD(3)) * t350 / 0.2e1 - ((-t324 * t313 + t337 * t323) * qJD(1) + (t324 ^ 2 * t281 + (t338 * t323 + (-t282 + t339) * t324) * t323) * qJD(3)) * t349 / 0.2e1 + qJD(1) * ((t331 * t314 + t328 * t315) * qJD(1) + ((t284 * t331 + t328 * t286) * t323 - (t283 * t331 + t328 * t285) * t324) * qJD(3)) / 0.2e1 + m(5) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + m(6) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + ((t297 * t371 + t298 * t370 + t357 * t372) * t320 + (t376 * t297 + t374 * t298 + t357 * t378) * t310 + (t297 * t375 + t298 * t373 + t357 * t377) * t309) * t309 / 0.2e1 + ((t295 * t371 + t296 * t370 + t359 * t372) * t320 + (t376 * t295 + t374 * t296 + t359 * t378) * t310 + (t295 * t375 + t296 * t373 + t359 * t377) * t309) * t310 / 0.2e1 + ((-t377 * t309 - t310 * t378 - t372 * t320) * t331 + ((-t327 * t371 + t330 * t370) * t320 + (-t327 * t376 + t330 * t374) * t310 + (-t327 * t375 + t330 * t373) * t309) * t328) * t320 / 0.2e1 + (m(2) * (t317 ^ 2 + t318 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
