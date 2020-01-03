% Calculate kinetic energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:40
% EndTime: 2019-12-31 18:23:42
% DurationCPUTime: 1.71s
% Computational Cost: add. (806->169), mult. (963->270), div. (0->0), fcn. (888->8), ass. (0->96)
t392 = Icges(4,4) + Icges(5,6);
t391 = Icges(4,1) + Icges(5,2);
t390 = -Icges(4,2) - Icges(5,3);
t326 = cos(qJ(3));
t389 = t392 * t326;
t323 = sin(qJ(3));
t388 = t392 * t323;
t387 = -Icges(5,4) + Icges(4,5);
t386 = Icges(5,5) - Icges(4,6);
t385 = t390 * t323 + t389;
t384 = -t391 * t326 + t388;
t383 = Icges(5,1) + Icges(4,3);
t321 = qJ(1) + pkin(8);
t319 = sin(t321);
t320 = cos(t321);
t382 = t385 * t319 + t386 * t320;
t381 = -t386 * t319 + t385 * t320;
t380 = t384 * t319 + t387 * t320;
t379 = t387 * t319 - t384 * t320;
t378 = t390 * t326 - t388;
t377 = t391 * t323 + t389;
t376 = t386 * t323 + t387 * t326;
t375 = t376 * t319 - t383 * t320;
t374 = t383 * t319 + t376 * t320;
t373 = t387 * t323 - t386 * t326;
t372 = t378 * t323 + t377 * t326;
t371 = -t381 * t323 + t379 * t326;
t370 = t382 * t323 + t380 * t326;
t324 = sin(qJ(1));
t366 = t324 * pkin(1);
t360 = t319 * t326;
t359 = t320 * t326;
t322 = sin(qJ(5));
t358 = t322 * t323;
t325 = cos(qJ(5));
t357 = t323 * t325;
t327 = cos(qJ(1));
t318 = qJD(1) * t327 * pkin(1);
t356 = qJD(1) * (t320 * pkin(2) + t319 * pkin(6)) + t318;
t355 = qJD(3) * t319;
t354 = qJD(3) * t320;
t353 = qJD(4) * t323;
t352 = qJD(5) * t326;
t349 = -t319 * pkin(2) + t320 * pkin(6) - t366;
t312 = t323 * pkin(3) - t326 * qJ(4);
t348 = qJD(3) * (t323 * rSges(5,2) + t326 * rSges(5,3) - t312);
t343 = pkin(3) * t326 + qJ(4) * t323;
t296 = t343 * t320;
t347 = qJD(1) * t296 + t319 * t353 + t356;
t295 = t343 * t319;
t346 = -t295 + t349;
t345 = rSges(4,1) * t326 - rSges(4,2) * t323;
t344 = -rSges(5,2) * t326 + rSges(5,3) * t323;
t342 = qJD(3) * (-pkin(7) * t323 - t312);
t329 = -qJD(4) * t326 + t295 * t355 + t296 * t354 + qJD(2);
t317 = qJD(5) * t323 + qJD(1);
t316 = t327 * rSges(2,1) - t324 * rSges(2,2);
t315 = t324 * rSges(2,1) + t327 * rSges(2,2);
t314 = t323 * rSges(4,1) + t326 * rSges(4,2);
t311 = t320 * t353;
t301 = -t320 * pkin(4) + pkin(7) * t360;
t300 = t319 * pkin(4) + pkin(7) * t359;
t299 = t319 * t352 - t354;
t298 = t320 * t352 + t355;
t297 = t323 * rSges(6,3) + (-rSges(6,1) * t322 - rSges(6,2) * t325) * t326;
t294 = Icges(6,5) * t323 + (-Icges(6,1) * t322 - Icges(6,4) * t325) * t326;
t293 = Icges(6,6) * t323 + (-Icges(6,4) * t322 - Icges(6,2) * t325) * t326;
t292 = Icges(6,3) * t323 + (-Icges(6,5) * t322 - Icges(6,6) * t325) * t326;
t291 = t319 * t358 - t320 * t325;
t290 = t319 * t357 + t320 * t322;
t289 = t319 * t325 + t320 * t358;
t288 = -t319 * t322 + t320 * t357;
t286 = t318 + qJD(1) * (t320 * rSges(3,1) - t319 * rSges(3,2));
t285 = (-t319 * rSges(3,1) - t320 * rSges(3,2) - t366) * qJD(1);
t284 = -t320 * rSges(5,1) + t319 * t344;
t283 = t319 * rSges(5,1) + t320 * t344;
t282 = t319 * rSges(4,3) + t320 * t345;
t281 = -t320 * rSges(4,3) + t319 * t345;
t266 = t291 * rSges(6,1) + t290 * rSges(6,2) + rSges(6,3) * t360;
t265 = t289 * rSges(6,1) + t288 * rSges(6,2) + rSges(6,3) * t359;
t264 = Icges(6,1) * t291 + Icges(6,4) * t290 + Icges(6,5) * t360;
t263 = Icges(6,1) * t289 + Icges(6,4) * t288 + Icges(6,5) * t359;
t262 = Icges(6,4) * t291 + Icges(6,2) * t290 + Icges(6,6) * t360;
t261 = Icges(6,4) * t289 + Icges(6,2) * t288 + Icges(6,6) * t359;
t260 = Icges(6,5) * t291 + Icges(6,6) * t290 + Icges(6,3) * t360;
t259 = Icges(6,5) * t289 + Icges(6,6) * t288 + Icges(6,3) * t359;
t258 = qJD(1) * t282 - t314 * t355 + t356;
t257 = -t314 * t354 + (-t281 + t349) * qJD(1);
t256 = qJD(2) + (t281 * t319 + t282 * t320) * qJD(3);
t255 = qJD(1) * t283 + t319 * t348 + t347;
t254 = t311 + t320 * t348 + (-t284 + t346) * qJD(1);
t253 = (t283 * t320 + t284 * t319) * qJD(3) + t329;
t252 = qJD(1) * t300 + t317 * t265 - t298 * t297 + t319 * t342 + t347;
t251 = -t317 * t266 + t299 * t297 + t311 + t320 * t342 + (-t301 + t346) * qJD(1);
t250 = -t299 * t265 + t298 * t266 + (t300 * t320 + t301 * t319) * qJD(3) + t329;
t1 = m(3) * (qJD(2) ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + m(4) * (t256 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + m(5) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + m(6) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t298 * ((t259 * t359 + t288 * t261 + t289 * t263) * t298 + (t260 * t359 + t288 * t262 + t289 * t264) * t299 + (t288 * t293 + t289 * t294 + t292 * t359) * t317) / 0.2e1 + t299 * ((t259 * t360 + t290 * t261 + t291 * t263) * t298 + (t260 * t360 + t290 * t262 + t291 * t264) * t299 + (t290 * t293 + t291 * t294 + t292 * t360) * t317) / 0.2e1 + t317 * ((t259 * t298 + t260 * t299 + t292 * t317) * t323 + ((-t261 * t325 - t263 * t322) * t298 + (-t262 * t325 - t264 * t322) * t299 + (-t293 * t325 - t294 * t322) * t317) * t326) / 0.2e1 + (((t380 * t323 - t382 * t326) * t320 + (t379 * t323 + t381 * t326) * t319) * qJD(3) + (t377 * t323 - t378 * t326) * qJD(1)) * qJD(1) / 0.2e1 + ((t374 * t319 ^ 2 + (t370 * t320 + (t371 - t375) * t319) * t320) * qJD(3) + (t373 * t319 + t372 * t320) * qJD(1)) * t355 / 0.2e1 - ((t375 * t320 ^ 2 + (t371 * t319 + (t370 - t374) * t320) * t319) * qJD(3) + (t372 * t319 - t373 * t320) * qJD(1)) * t354 / 0.2e1 + (m(2) * (t315 ^ 2 + t316 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
