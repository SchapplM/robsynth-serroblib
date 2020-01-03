% Calculate kinetic energy for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:10:58
% EndTime: 2019-12-31 21:10:59
% DurationCPUTime: 1.50s
% Computational Cost: add. (822->157), mult. (958->250), div. (0->0), fcn. (906->8), ass. (0->96)
t389 = Icges(4,4) - Icges(5,5);
t388 = Icges(4,1) + Icges(5,1);
t387 = Icges(4,2) + Icges(5,3);
t324 = cos(qJ(3));
t386 = t389 * t324;
t321 = sin(qJ(3));
t385 = t389 * t321;
t384 = Icges(5,4) + Icges(4,5);
t383 = Icges(4,6) - Icges(5,6);
t382 = t387 * t321 - t386;
t381 = t388 * t324 - t385;
t380 = Icges(5,2) + Icges(4,3);
t319 = qJ(1) + qJ(2);
t316 = sin(t319);
t317 = cos(t319);
t379 = t382 * t316 + t383 * t317;
t378 = -t383 * t316 + t382 * t317;
t377 = -t381 * t316 + t384 * t317;
t376 = t384 * t316 + t381 * t317;
t375 = -t387 * t324 - t385;
t374 = t388 * t321 + t386;
t373 = -t383 * t321 + t384 * t324;
t372 = t373 * t316 - t380 * t317;
t371 = t380 * t316 + t373 * t317;
t370 = t384 * t321 + t383 * t324;
t369 = t375 * t321 + t374 * t324;
t368 = t378 * t321 + t376 * t324;
t367 = -t379 * t321 + t377 * t324;
t362 = pkin(4) * t324;
t361 = pkin(1) * qJD(1);
t341 = pkin(3) * t324 + qJ(4) * t321;
t292 = t341 * t316;
t299 = t316 * pkin(2) - t317 * pkin(7);
t356 = -t292 - t299;
t325 = cos(qJ(1));
t315 = t325 * t361;
t318 = qJD(1) + qJD(2);
t355 = t318 * (t317 * pkin(2) + t316 * pkin(7)) + t315;
t354 = qJD(3) * t316;
t353 = qJD(3) * t317;
t352 = qJD(4) * t321;
t351 = qJD(3) - qJD(5);
t322 = sin(qJ(1));
t350 = t322 * t361;
t308 = t321 * pkin(3) - t324 * qJ(4);
t347 = qJD(3) * (-t321 * rSges(5,1) + t324 * rSges(5,3) - t308);
t293 = t341 * t317;
t346 = t318 * t293 + t316 * t352 + t355;
t345 = t317 * t352 - t350;
t344 = -qJD(4) * t324 + t292 * t354 + t293 * t353;
t343 = rSges(4,1) * t324 - rSges(4,2) * t321;
t342 = rSges(5,1) * t324 + rSges(5,3) * t321;
t340 = qJD(3) * (-pkin(4) * t321 - t308);
t320 = sin(qJ(5));
t323 = cos(qJ(5));
t301 = -t324 * t320 + t321 * t323;
t327 = t321 * t320 + t324 * t323;
t312 = t325 * rSges(2,1) - t322 * rSges(2,2);
t311 = t322 * rSges(2,1) + t325 * rSges(2,2);
t310 = t321 * rSges(4,1) + t324 * rSges(4,2);
t298 = t351 * t317;
t297 = t351 * t316;
t296 = -t316 * pkin(8) + t317 * t362;
t295 = t317 * pkin(8) + t316 * t362;
t291 = t327 * t317;
t290 = t301 * t317;
t289 = t327 * t316;
t288 = t301 * t316;
t286 = t315 + t318 * (t317 * rSges(3,1) - t316 * rSges(3,2));
t285 = -t350 - t318 * (t316 * rSges(3,1) + t317 * rSges(3,2));
t284 = t316 * rSges(4,3) + t343 * t317;
t283 = t316 * rSges(5,2) + t342 * t317;
t282 = -t317 * rSges(4,3) + t343 * t316;
t281 = -t317 * rSges(5,2) + t342 * t316;
t266 = t301 * rSges(6,1) - rSges(6,2) * t327;
t265 = Icges(6,1) * t301 - Icges(6,4) * t327;
t264 = Icges(6,4) * t301 - Icges(6,2) * t327;
t263 = Icges(6,5) * t301 - Icges(6,6) * t327;
t262 = t291 * rSges(6,1) + t290 * rSges(6,2) - t316 * rSges(6,3);
t261 = t289 * rSges(6,1) + t288 * rSges(6,2) + t317 * rSges(6,3);
t260 = Icges(6,1) * t291 + Icges(6,4) * t290 - Icges(6,5) * t316;
t259 = Icges(6,1) * t289 + Icges(6,4) * t288 + Icges(6,5) * t317;
t258 = Icges(6,4) * t291 + Icges(6,2) * t290 - Icges(6,6) * t316;
t257 = Icges(6,4) * t289 + Icges(6,2) * t288 + Icges(6,6) * t317;
t256 = Icges(6,5) * t291 + Icges(6,6) * t290 - Icges(6,3) * t316;
t255 = Icges(6,5) * t289 + Icges(6,6) * t288 + Icges(6,3) * t317;
t254 = (t282 * t316 + t284 * t317) * qJD(3);
t253 = t318 * t284 - t310 * t354 + t355;
t252 = -t350 - t310 * t353 + (-t282 - t299) * t318;
t251 = t318 * t283 + t316 * t347 + t346;
t250 = t317 * t347 + (-t281 + t356) * t318 + t345;
t249 = (t281 * t316 + t283 * t317) * qJD(3) + t344;
t248 = -t297 * t266 + (t262 + t296) * t318 + t316 * t340 + t346;
t247 = -t298 * t266 + t317 * t340 + (-t261 - t295 + t356) * t318 + t345;
t246 = t297 * t261 + t298 * t262 + (t295 * t316 + t296 * t317) * qJD(3) + t344;
t1 = m(3) * (t285 ^ 2 + t286 ^ 2) / 0.2e1 + t318 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(5) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(6) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t297 * ((-t316 * t256 + t290 * t258 + t291 * t260) * t297 - (-t316 * t255 + t290 * t257 + t291 * t259) * t298 + (-t316 * t263 + t290 * t264 + t291 * t265) * t318) / 0.2e1 - t298 * ((t317 * t256 + t288 * t258 + t289 * t260) * t297 - (t317 * t255 + t288 * t257 + t289 * t259) * t298 + (t317 * t263 + t288 * t264 + t289 * t265) * t318) / 0.2e1 + (m(2) * (t311 ^ 2 + t312 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t370 * t316 + t369 * t317) * t318 + (t371 * t316 ^ 2 + (t367 * t317 + (t368 - t372) * t316) * t317) * qJD(3)) * t354 / 0.2e1 - ((t369 * t316 - t370 * t317) * t318 + (t372 * t317 ^ 2 + (t368 * t316 + (t367 - t371) * t317) * t316) * qJD(3)) * t353 / 0.2e1 + ((-t258 * t327 + t301 * t260) * t297 - (-t257 * t327 + t301 * t259) * t298 + ((t377 * t321 + t379 * t324) * t317 + (t376 * t321 - t378 * t324) * t316) * qJD(3) + (-t327 * t264 + t301 * t265 + t374 * t321 - t375 * t324) * t318) * t318 / 0.2e1;
T = t1;
