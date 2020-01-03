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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:03
% EndTime: 2020-01-03 12:11:04
% DurationCPUTime: 1.48s
% Computational Cost: add. (964->138), mult. (778->218), div. (0->0), fcn. (648->8), ass. (0->93)
t404 = Icges(5,4) + Icges(6,4);
t403 = Icges(5,1) + Icges(6,1);
t402 = Icges(5,2) + Icges(6,2);
t322 = qJ(3) + qJ(4);
t317 = cos(t322);
t401 = t404 * t317;
t315 = sin(t322);
t400 = t404 * t315;
t399 = Icges(5,5) + Icges(6,5);
t398 = -Icges(5,6) - Icges(6,6);
t397 = -t402 * t315 + t401;
t396 = t403 * t317 - t400;
t395 = rSges(6,1) + pkin(4);
t394 = Icges(5,3) + Icges(6,3);
t323 = qJ(1) + qJ(2);
t316 = sin(t323);
t318 = cos(t323);
t393 = t397 * t316 + t398 * t318;
t392 = t398 * t316 - t397 * t318;
t391 = t396 * t316 - t399 * t318;
t390 = -t399 * t316 - t396 * t318;
t389 = t402 * t317 + t400;
t388 = t403 * t315 + t401;
t387 = t398 * t315 + t399 * t317;
t386 = -rSges(6,3) - qJ(5);
t385 = -rSges(6,2) * t315 + t395 * t317;
t358 = -qJD(3) - qJD(4);
t293 = t358 * t316;
t294 = t358 * t318;
t321 = qJD(1) + qJD(2);
t384 = t293 * (t392 * t315 - t390 * t317) + t294 * (t393 * t315 - t391 * t317) + t321 * (t389 * t315 - t388 * t317);
t383 = (-t399 * t315 + t398 * t317) * t321 + (-t387 * t316 + t394 * t318) * t294 + (t394 * t316 + t387 * t318) * t293;
t326 = cos(qJ(3));
t375 = t326 * pkin(3);
t373 = pkin(1) * qJD(1);
t324 = sin(qJ(3));
t372 = Icges(4,4) * t324;
t371 = Icges(4,4) * t326;
t366 = t385 * t316 + t386 * t318;
t365 = t386 * t316 - t385 * t318;
t265 = -pkin(8) * t316 - t375 * t318;
t303 = -pkin(2) * t318 - pkin(7) * t316;
t364 = -t265 - t303;
t325 = sin(qJ(1));
t313 = t325 * t373;
t363 = t321 * (pkin(2) * t316 - pkin(7) * t318) + t313;
t360 = qJD(3) * t316;
t359 = qJD(3) * t318;
t357 = pkin(3) * qJD(3) * t324;
t356 = rSges(6,2) * t317 + t395 * t315;
t264 = -pkin(8) * t318 + t375 * t316;
t355 = t321 * t264 + t318 * t357 + t363;
t354 = t264 * t360 - t265 * t359;
t353 = rSges(4,1) * t326 - rSges(4,2) * t324;
t352 = rSges(5,1) * t317 - rSges(5,2) * t315;
t350 = Icges(4,1) * t326 - t372;
t347 = -Icges(4,2) * t324 + t371;
t344 = Icges(4,5) * t326 - Icges(4,6) * t324;
t284 = -Icges(4,6) * t318 + t347 * t316;
t286 = -Icges(4,5) * t318 + t350 * t316;
t337 = -t284 * t324 + t286 * t326;
t285 = -Icges(4,6) * t316 - t347 * t318;
t287 = -Icges(4,5) * t316 - t350 * t318;
t336 = t285 * t324 - t287 * t326;
t307 = Icges(4,2) * t326 + t372;
t308 = Icges(4,1) * t324 + t371;
t333 = t307 * t324 - t308 * t326;
t327 = cos(qJ(1));
t314 = t327 * t373;
t332 = -t316 * t357 + t314;
t311 = -rSges(2,1) * t327 + rSges(2,2) * t325;
t310 = rSges(2,1) * t325 + rSges(2,2) * t327;
t309 = rSges(4,1) * t324 + rSges(4,2) * t326;
t306 = Icges(4,5) * t324 + Icges(4,6) * t326;
t302 = rSges(5,1) * t315 + rSges(5,2) * t317;
t291 = t314 - t321 * (-rSges(3,1) * t318 + rSges(3,2) * t316);
t290 = t313 + t321 * (rSges(3,1) * t316 + rSges(3,2) * t318);
t289 = -rSges(4,3) * t316 - t353 * t318;
t288 = -rSges(4,3) * t318 + t353 * t316;
t283 = -Icges(4,3) * t316 - t344 * t318;
t282 = -Icges(4,3) * t318 + t344 * t316;
t281 = -rSges(5,3) * t316 - t352 * t318;
t279 = -rSges(5,3) * t318 + t352 * t316;
t259 = (t288 * t316 - t289 * t318) * qJD(3);
t258 = -t309 * t360 + t314 + (-t289 - t303) * t321;
t257 = t288 * t321 + t309 * t359 + t363;
t256 = t293 * t302 + (-t281 + t364) * t321 + t332;
t255 = t279 * t321 - t294 * t302 + t355;
t254 = -t279 * t293 + t281 * t294 + t354;
t253 = -qJD(5) * t318 + t356 * t293 + (t364 - t365) * t321 + t332;
t252 = -qJD(5) * t316 - t356 * t294 + t366 * t321 + t355;
t251 = -t366 * t293 + t365 * t294 + t354;
t1 = m(3) * (t290 ^ 2 + t291 ^ 2) / 0.2e1 + t321 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 - ((-t318 * t306 - t333 * t316) * t321 + (t318 ^ 2 * t282 + (t336 * t316 + (t283 - t337) * t318) * t316) * qJD(3)) * t359 / 0.2e1 - ((-t316 * t306 + t333 * t318) * t321 + (t316 ^ 2 * t283 + (t337 * t318 + (t282 - t336) * t316) * t318) * qJD(3)) * t360 / 0.2e1 + m(5) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + m(6) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + (t383 * t316 + t384 * t318) * t293 / 0.2e1 + (-t384 * t316 + t383 * t318) * t294 / 0.2e1 + (m(2) * (t310 ^ 2 + t311 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((-(t284 * t326 + t286 * t324) * t318 - (t285 * t326 + t287 * t324) * t316) * qJD(3) + (t391 * t315 + t393 * t317) * t294 + (t390 * t315 + t392 * t317) * t293 + (t307 * t326 + t308 * t324 + t388 * t315 + t389 * t317) * t321) * t321 / 0.2e1;
T = t1;
