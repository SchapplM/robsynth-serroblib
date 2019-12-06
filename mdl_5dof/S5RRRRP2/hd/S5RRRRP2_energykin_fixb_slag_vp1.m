% Calculate kinetic energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:31
% EndTime: 2019-12-05 18:47:33
% DurationCPUTime: 1.34s
% Computational Cost: add. (964->138), mult. (778->218), div. (0->0), fcn. (648->8), ass. (0->93)
t406 = Icges(5,4) + Icges(6,4);
t405 = Icges(5,1) + Icges(6,1);
t404 = Icges(5,2) + Icges(6,2);
t322 = qJ(3) + qJ(4);
t317 = cos(t322);
t403 = t406 * t317;
t315 = sin(t322);
t402 = t406 * t315;
t401 = Icges(5,5) + Icges(6,5);
t400 = Icges(5,6) + Icges(6,6);
t399 = -t404 * t315 + t403;
t398 = t405 * t317 - t402;
t397 = rSges(6,1) + pkin(4);
t396 = Icges(5,3) + Icges(6,3);
t323 = qJ(1) + qJ(2);
t316 = sin(t323);
t318 = cos(t323);
t395 = -t399 * t316 + t400 * t318;
t394 = t400 * t316 + t399 * t318;
t393 = -t398 * t316 + t401 * t318;
t392 = t401 * t316 + t398 * t318;
t391 = t404 * t317 + t402;
t390 = t405 * t315 + t403;
t389 = -t400 * t315 + t401 * t317;
t388 = rSges(6,3) + qJ(5);
t387 = -rSges(6,2) * t315 + t397 * t317;
t361 = qJD(3) + qJD(4);
t295 = t361 * t316;
t296 = t361 * t318;
t321 = qJD(1) + qJD(2);
t386 = (t315 * t394 - t317 * t392) * t295 + (t315 * t395 - t317 * t393) * t296 + (t315 * t391 - t317 * t390) * t321;
t385 = (t401 * t315 + t400 * t317) * t321 + (-t389 * t316 + t318 * t396) * t296 + (t316 * t396 + t389 * t318) * t295;
t326 = cos(qJ(3));
t377 = t326 * pkin(3);
t375 = pkin(1) * qJD(1);
t324 = sin(qJ(3));
t374 = Icges(4,4) * t324;
t373 = Icges(4,4) * t326;
t368 = -t387 * t316 + t388 * t318;
t367 = t388 * t316 + t387 * t318;
t267 = pkin(8) * t316 + t377 * t318;
t305 = pkin(2) * t318 + pkin(7) * t316;
t366 = -t267 - t305;
t363 = qJD(3) * t316;
t362 = qJD(3) * t318;
t360 = pkin(3) * qJD(3) * t324;
t325 = sin(qJ(1));
t359 = t325 * t375;
t327 = cos(qJ(1));
t358 = t327 * t375;
t357 = rSges(6,2) * t317 + t397 * t315;
t356 = t321 * (-pkin(2) * t316 + pkin(7) * t318) - t359;
t355 = t316 * t360 - t358;
t266 = pkin(8) * t318 - t377 * t316;
t354 = -t266 * t363 + t267 * t362;
t353 = rSges(4,1) * t326 - rSges(4,2) * t324;
t352 = rSges(5,1) * t317 - rSges(5,2) * t315;
t350 = Icges(4,1) * t326 - t374;
t347 = -Icges(4,2) * t324 + t373;
t344 = Icges(4,5) * t326 - Icges(4,6) * t324;
t286 = Icges(4,6) * t318 - t347 * t316;
t288 = Icges(4,5) * t318 - t350 * t316;
t337 = -t286 * t324 + t288 * t326;
t287 = Icges(4,6) * t316 + t347 * t318;
t289 = Icges(4,5) * t316 + t350 * t318;
t336 = t287 * t324 - t289 * t326;
t309 = Icges(4,2) * t326 + t374;
t310 = Icges(4,1) * t324 + t373;
t333 = t309 * t324 - t310 * t326;
t330 = t321 * t266 - t318 * t360 + t356;
t313 = rSges(2,1) * t327 - rSges(2,2) * t325;
t312 = -rSges(2,1) * t325 - rSges(2,2) * t327;
t311 = rSges(4,1) * t324 + rSges(4,2) * t326;
t308 = Icges(4,5) * t324 + Icges(4,6) * t326;
t304 = rSges(5,1) * t315 + rSges(5,2) * t317;
t293 = -t358 - t321 * (rSges(3,1) * t318 - rSges(3,2) * t316);
t292 = -t359 + t321 * (-rSges(3,1) * t316 - rSges(3,2) * t318);
t291 = rSges(4,3) * t316 + t353 * t318;
t290 = rSges(4,3) * t318 - t353 * t316;
t285 = Icges(4,3) * t316 + t344 * t318;
t284 = Icges(4,3) * t318 - t344 * t316;
t283 = rSges(5,3) * t316 + t352 * t318;
t281 = rSges(5,3) * t318 - t352 * t316;
t261 = (-t290 * t316 + t291 * t318) * qJD(3);
t260 = -t358 + t311 * t363 + (-t291 - t305) * t321;
t259 = t290 * t321 - t311 * t362 + t356;
t258 = t295 * t304 + (-t283 + t366) * t321 + t355;
t257 = t281 * t321 - t296 * t304 + t330;
t256 = -t281 * t295 + t283 * t296 + t354;
t255 = qJD(5) * t318 + t357 * t295 + (t366 - t367) * t321 + t355;
t254 = qJD(5) * t316 - t357 * t296 + t368 * t321 + t330;
t253 = -t368 * t295 + t367 * t296 + t354;
t1 = m(3) * (t292 ^ 2 + t293 ^ 2) / 0.2e1 + t321 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + ((t318 * t308 + t333 * t316) * t321 + (t318 ^ 2 * t284 + (t336 * t316 + (t285 - t337) * t318) * t316) * qJD(3)) * t362 / 0.2e1 + ((t316 * t308 - t333 * t318) * t321 + (t316 ^ 2 * t285 + (t337 * t318 + (t284 - t336) * t316) * t318) * qJD(3)) * t363 / 0.2e1 + m(5) * (t256 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + m(6) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + (t385 * t316 - t386 * t318) * t295 / 0.2e1 + (t386 * t316 + t385 * t318) * t296 / 0.2e1 + (m(2) * (t312 ^ 2 + t313 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t286 * t326 + t288 * t324) * t318 + (t287 * t326 + t289 * t324) * t316) * qJD(3) + (t315 * t393 + t317 * t395) * t296 + (t315 * t392 + t317 * t394) * t295 + (t326 * t309 + t324 * t310 + t315 * t390 + t391 * t317) * t321) * t321 / 0.2e1;
T = t1;
