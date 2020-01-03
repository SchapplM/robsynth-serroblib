% Calculate kinetic energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR16_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR16_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:42
% EndTime: 2019-12-31 18:38:44
% DurationCPUTime: 1.65s
% Computational Cost: add. (412->175), mult. (972->274), div. (0->0), fcn. (898->6), ass. (0->93)
t390 = -Icges(4,4) - Icges(5,6);
t389 = Icges(4,1) + Icges(5,2);
t388 = Icges(4,2) + Icges(5,3);
t326 = cos(qJ(3));
t387 = t390 * t326;
t323 = sin(qJ(3));
t386 = t390 * t323;
t385 = -Icges(5,4) + Icges(4,5);
t384 = Icges(5,5) - Icges(4,6);
t383 = -t388 * t326 + t386;
t382 = t389 * t323 - t387;
t381 = Icges(5,1) + Icges(4,3);
t324 = sin(qJ(1));
t327 = cos(qJ(1));
t380 = t383 * t324 + t384 * t327;
t379 = -t384 * t324 + t383 * t327;
t378 = t382 * t324 + t385 * t327;
t377 = t385 * t324 - t382 * t327;
t376 = t388 * t323 + t387;
t375 = t389 * t326 + t386;
t374 = t385 * t323 - t384 * t326;
t373 = t374 * t324 + t327 * t381;
t372 = t324 * t381 - t374 * t327;
t371 = t384 * t323 + t385 * t326;
t370 = t323 * t375 - t326 * t376;
t369 = t323 * t377 + t326 * t379;
t368 = -t323 * t378 + t326 * t380;
t359 = t323 * t324;
t358 = t323 * t327;
t357 = t324 * t326;
t356 = t326 * t327;
t341 = pkin(3) * t323 - qJ(4) * t326;
t298 = t341 * t327;
t351 = qJD(3) * t327;
t355 = qJD(4) * t323 - t298 * t351;
t313 = pkin(3) * t326 + qJ(4) * t323;
t321 = qJD(2) * t324;
t352 = qJD(3) * t324;
t354 = t313 * t352 + t321;
t304 = qJD(1) * (pkin(1) * t327 + qJ(2) * t324);
t353 = qJD(1) * t327 * pkin(6) + t304;
t350 = qJD(4) * t326;
t349 = qJD(5) * t323;
t311 = pkin(1) * t324 - qJ(2) * t327;
t346 = -pkin(6) * t324 - t311;
t297 = t341 * t324;
t345 = qJD(1) * t297 + t327 * t350 + t353;
t344 = t298 + t346;
t343 = rSges(4,1) * t323 + rSges(4,2) * t326;
t342 = rSges(5,2) * t323 + rSges(5,3) * t326;
t325 = cos(qJ(5));
t322 = sin(qJ(5));
t318 = qJD(5) * t326 + qJD(1);
t316 = rSges(2,1) * t327 - rSges(2,2) * t324;
t315 = rSges(4,1) * t326 - rSges(4,2) * t323;
t314 = -rSges(5,2) * t326 + rSges(5,3) * t323;
t312 = rSges(2,1) * t324 + rSges(2,2) * t327;
t303 = pkin(4) * t327 + pkin(7) * t359;
t302 = pkin(4) * t324 - pkin(7) * t358;
t301 = t324 * t349 + t351;
t300 = -t327 * t349 + t352;
t296 = -t322 * t357 + t325 * t327;
t295 = -t322 * t327 - t325 * t357;
t294 = t322 * t356 + t324 * t325;
t293 = -t322 * t324 + t325 * t356;
t291 = rSges(5,1) * t327 - t324 * t342;
t290 = rSges(5,1) * t324 + t327 * t342;
t289 = rSges(4,3) * t324 - t327 * t343;
t288 = rSges(4,3) * t327 + t324 * t343;
t287 = rSges(6,3) * t326 + (rSges(6,1) * t322 + rSges(6,2) * t325) * t323;
t277 = Icges(6,5) * t326 + (Icges(6,1) * t322 + Icges(6,4) * t325) * t323;
t274 = Icges(6,6) * t326 + (Icges(6,4) * t322 + Icges(6,2) * t325) * t323;
t271 = Icges(6,3) * t326 + (Icges(6,5) * t322 + Icges(6,6) * t325) * t323;
t270 = t304 - qJD(2) * t327 + qJD(1) * (-rSges(3,2) * t327 + rSges(3,3) * t324);
t269 = t321 + (rSges(3,2) * t324 + rSges(3,3) * t327 - t311) * qJD(1);
t268 = rSges(6,1) * t296 + rSges(6,2) * t295 + rSges(6,3) * t359;
t267 = rSges(6,1) * t294 + rSges(6,2) * t293 - rSges(6,3) * t358;
t266 = Icges(6,1) * t296 + Icges(6,4) * t295 + Icges(6,5) * t359;
t265 = Icges(6,1) * t294 + Icges(6,4) * t293 - Icges(6,5) * t358;
t264 = Icges(6,4) * t296 + Icges(6,2) * t295 + Icges(6,6) * t359;
t263 = Icges(6,4) * t294 + Icges(6,2) * t293 - Icges(6,6) * t358;
t262 = Icges(6,5) * t296 + Icges(6,6) * t295 + Icges(6,3) * t359;
t261 = Icges(6,5) * t294 + Icges(6,6) * t293 - Icges(6,3) * t358;
t260 = (-t288 * t324 + t289 * t327) * qJD(3);
t259 = qJD(1) * t288 + (-qJD(3) * t315 - qJD(2)) * t327 + t353;
t258 = t315 * t352 + t321 + (-t289 + t346) * qJD(1);
t257 = (t290 * t327 + (-t291 - t297) * t324) * qJD(3) + t355;
t256 = qJD(1) * t291 + (-qJD(2) + (-t313 - t314) * qJD(3)) * t327 + t345;
t255 = (qJD(3) * t314 - t350) * t324 + (-t290 + t344) * qJD(1) + t354;
t254 = qJD(1) * t303 + t268 * t318 - t287 * t301 + (-qJD(2) + (-pkin(7) * t326 - t313) * qJD(3)) * t327 + t345;
t253 = -t267 * t318 + t287 * t300 + (pkin(7) * qJD(3) - qJD(4)) * t357 + (-t302 + t344) * qJD(1) + t354;
t252 = t267 * t301 - t268 * t300 + (t302 * t327 + (-t297 - t303) * t324) * qJD(3) + t355;
t1 = m(3) * (t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(4) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(5) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + m(6) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + t301 * ((t262 * t359 + t295 * t264 + t296 * t266) * t301 + (t261 * t359 + t263 * t295 + t265 * t296) * t300 + (t271 * t359 + t274 * t295 + t277 * t296) * t318) / 0.2e1 + t300 * ((-t262 * t358 + t264 * t293 + t266 * t294) * t301 + (-t261 * t358 + t293 * t263 + t294 * t265) * t300 + (-t271 * t358 + t274 * t293 + t277 * t294) * t318) / 0.2e1 + t318 * ((t261 * t300 + t262 * t301 + t271 * t318) * t326 + ((t264 * t325 + t266 * t322) * t301 + (t263 * t325 + t265 * t322) * t300 + (t274 * t325 + t277 * t322) * t318) * t323) / 0.2e1 + (((t323 * t380 + t326 * t378) * t327 + (-t323 * t379 + t326 * t377) * t324) * qJD(3) + (t376 * t323 + t375 * t326) * qJD(1)) * qJD(1) / 0.2e1 + ((t372 * t324 ^ 2 + (t368 * t327 + (-t369 + t373) * t324) * t327) * qJD(3) + (t371 * t324 - t370 * t327) * qJD(1)) * t352 / 0.2e1 + ((t373 * t327 ^ 2 + (t369 * t324 + (-t368 + t372) * t327) * t324) * qJD(3) + (t370 * t324 + t371 * t327) * qJD(1)) * t351 / 0.2e1 + (m(2) * (t312 ^ 2 + t316 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
