% Calculate kinetic energy for
% S5RPRRP3
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:08
% EndTime: 2019-12-05 18:03:09
% DurationCPUTime: 1.39s
% Computational Cost: add. (935->141), mult. (779->215), div. (0->0), fcn. (648->8), ass. (0->91)
t404 = Icges(5,4) + Icges(6,4);
t403 = Icges(5,1) + Icges(6,1);
t402 = Icges(5,2) + Icges(6,2);
t323 = qJ(3) + qJ(4);
t319 = cos(t323);
t401 = t404 * t319;
t318 = sin(t323);
t400 = t404 * t318;
t399 = Icges(5,5) + Icges(6,5);
t398 = Icges(5,6) + Icges(6,6);
t397 = -t402 * t318 + t401;
t396 = t403 * t319 - t400;
t395 = rSges(6,1) + pkin(4);
t394 = Icges(5,3) + Icges(6,3);
t322 = qJ(1) + pkin(8);
t316 = sin(t322);
t317 = cos(t322);
t393 = -t397 * t316 + t398 * t317;
t392 = t398 * t316 + t397 * t317;
t391 = -t396 * t316 + t399 * t317;
t390 = t399 * t316 + t396 * t317;
t389 = t402 * t319 + t400;
t388 = t403 * t318 + t401;
t387 = -t398 * t318 + t399 * t319;
t386 = rSges(6,3) + qJ(5);
t385 = -rSges(6,2) * t318 + t395 * t319;
t359 = qJD(3) + qJD(4);
t296 = t359 * t316;
t297 = t359 * t317;
t384 = (t389 * t318 - t388 * t319) * qJD(1) + (t392 * t318 - t390 * t319) * t296 + (t393 * t318 - t391 * t319) * t297;
t383 = (-t387 * t316 + t394 * t317) * t297 + (t394 * t316 + t387 * t317) * t296 + (t399 * t318 + t398 * t319) * qJD(1);
t325 = sin(qJ(1));
t376 = pkin(1) * t325;
t327 = cos(qJ(1));
t375 = pkin(1) * t327;
t326 = cos(qJ(3));
t373 = t326 * pkin(3);
t324 = sin(qJ(3));
t371 = Icges(4,4) * t324;
t370 = Icges(4,4) * t326;
t365 = t385 * t316 - t386 * t317;
t364 = t386 * t316 + t385 * t317;
t361 = qJD(3) * t316;
t360 = qJD(3) * t317;
t358 = pkin(3) * qJD(3) * t324;
t357 = -pkin(2) * t317 - pkin(6) * t316 - t375;
t356 = rSges(6,2) * t319 + t395 * t318;
t268 = pkin(7) * t316 + t373 * t317;
t355 = -t268 + t357;
t354 = rSges(4,1) * t326 - rSges(4,2) * t324;
t353 = rSges(5,1) * t319 - rSges(5,2) * t318;
t351 = Icges(4,1) * t326 - t371;
t348 = -Icges(4,2) * t324 + t370;
t345 = Icges(4,5) * t326 - Icges(4,6) * t324;
t287 = Icges(4,6) * t317 - t348 * t316;
t289 = Icges(4,5) * t317 - t351 * t316;
t338 = -t287 * t324 + t289 * t326;
t288 = Icges(4,6) * t316 + t348 * t317;
t290 = Icges(4,5) * t316 + t351 * t317;
t337 = t288 * t324 - t290 * t326;
t310 = Icges(4,2) * t326 + t371;
t311 = Icges(4,1) * t324 + t370;
t334 = t310 * t324 - t311 * t326;
t267 = pkin(7) * t317 - t373 * t316;
t333 = -t267 * t361 + t268 * t360 + qJD(2);
t295 = qJD(1) * (-pkin(2) * t316 + pkin(6) * t317);
t332 = qJD(1) * t267 - t317 * t358 + t295;
t314 = rSges(2,1) * t327 - rSges(2,2) * t325;
t313 = -rSges(2,1) * t325 - rSges(2,2) * t327;
t312 = rSges(4,1) * t324 + rSges(4,2) * t326;
t309 = Icges(4,5) * t324 + Icges(4,6) * t326;
t308 = t316 * t358;
t306 = rSges(5,1) * t318 + rSges(5,2) * t319;
t294 = (-rSges(3,1) * t317 + rSges(3,2) * t316 - t375) * qJD(1);
t293 = (-rSges(3,1) * t316 - rSges(3,2) * t317 - t376) * qJD(1);
t292 = rSges(4,3) * t316 + t354 * t317;
t291 = rSges(4,3) * t317 - t354 * t316;
t286 = Icges(4,3) * t316 + t345 * t317;
t285 = Icges(4,3) * t317 - t345 * t316;
t284 = rSges(5,3) * t316 + t353 * t317;
t282 = rSges(5,3) * t317 - t353 * t316;
t262 = t312 * t361 + (-t292 + t357) * qJD(1);
t261 = -t312 * t360 + t295 + (t291 - t376) * qJD(1);
t260 = qJD(2) + (-t291 * t316 + t292 * t317) * qJD(3);
t259 = t296 * t306 + t308 + (-t284 + t355) * qJD(1);
t258 = -t297 * t306 + (t282 - t376) * qJD(1) + t332;
t257 = -t282 * t296 + t284 * t297 + t333;
t256 = qJD(5) * t317 + t308 + t356 * t296 + (t355 - t364) * qJD(1);
t255 = qJD(5) * t316 - t356 * t297 + (-t365 - t376) * qJD(1) + t332;
t254 = t365 * t296 + t364 * t297 + t333;
t1 = m(3) * (qJD(2) ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + m(4) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + ((t317 * t309 + t334 * t316) * qJD(1) + (t317 ^ 2 * t285 + (t337 * t316 + (t286 - t338) * t317) * t316) * qJD(3)) * t360 / 0.2e1 + ((t316 * t309 - t334 * t317) * qJD(1) + (t316 ^ 2 * t286 + (t338 * t317 + (t285 - t337) * t316) * t317) * qJD(3)) * t361 / 0.2e1 + m(5) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + m(6) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + (t383 * t316 - t384 * t317) * t296 / 0.2e1 + (t384 * t316 + t383 * t317) * t297 / 0.2e1 + (m(2) * (t313 ^ 2 + t314 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + (((t287 * t326 + t289 * t324) * t317 + (t288 * t326 + t290 * t324) * t316) * qJD(3) + (t391 * t318 + t393 * t319) * t297 + (t390 * t318 + t392 * t319) * t296 + (t326 * t310 + t324 * t311 + t388 * t318 + t389 * t319) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
