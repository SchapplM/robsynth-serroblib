% Calculate kinetic energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:12
% EndTime: 2019-12-05 17:59:14
% DurationCPUTime: 1.43s
% Computational Cost: add. (541->141), mult. (780->218), div. (0->0), fcn. (650->6), ass. (0->89)
t405 = Icges(5,4) + Icges(6,4);
t404 = Icges(5,1) + Icges(6,1);
t403 = Icges(5,2) + Icges(6,2);
t324 = qJ(3) + qJ(4);
t322 = cos(t324);
t402 = t405 * t322;
t321 = sin(t324);
t401 = t405 * t321;
t400 = Icges(5,5) + Icges(6,5);
t399 = Icges(5,6) + Icges(6,6);
t398 = t403 * t322 + t401;
t397 = t404 * t321 + t402;
t396 = rSges(6,1) + pkin(4);
t395 = Icges(5,3) + Icges(6,3);
t326 = sin(qJ(1));
t328 = cos(qJ(1));
t394 = t398 * t326 + t399 * t328;
t393 = t399 * t326 - t398 * t328;
t392 = t397 * t326 + t400 * t328;
t391 = t400 * t326 - t397 * t328;
t390 = -t403 * t321 + t402;
t389 = t404 * t322 - t401;
t388 = t400 * t321 + t399 * t322;
t387 = rSges(6,3) + qJ(5);
t386 = rSges(6,2) * t322 + t396 * t321;
t361 = qJD(3) + qJD(4);
t309 = t361 * t326;
t310 = t361 * t328;
t385 = (t389 * t321 + t390 * t322) * qJD(1) + (t391 * t321 + t393 * t322) * t309 + (t392 * t321 + t394 * t322) * t310;
t384 = (t388 * t326 + t395 * t328) * t310 + (t395 * t326 - t388 * t328) * t309 + (-t399 * t321 + t400 * t322) * qJD(1);
t325 = sin(qJ(3));
t377 = pkin(3) * t325;
t374 = Icges(4,4) * t325;
t327 = cos(qJ(3));
t373 = Icges(4,4) * t327;
t368 = t387 * t326 - t386 * t328;
t367 = t386 * t326 + t387 * t328;
t307 = qJD(1) * (pkin(1) * t328 + qJ(2) * t326);
t366 = qJD(1) * t328 * pkin(6) + t307;
t320 = qJD(2) * t326;
t360 = pkin(3) * qJD(3) * t327;
t365 = t326 * t360 + t320;
t363 = qJD(3) * t326;
t362 = qJD(3) * t328;
t358 = -rSges(6,2) * t321 + t396 * t322;
t314 = pkin(1) * t326 - qJ(2) * t328;
t357 = -pkin(6) * t326 - t314;
t297 = pkin(7) * t326 - t328 * t377;
t356 = -t297 + t357;
t298 = pkin(7) * t328 + t326 * t377;
t355 = t297 * t362 - t298 * t363;
t354 = rSges(4,1) * t325 + rSges(4,2) * t327;
t353 = rSges(5,1) * t321 + rSges(5,2) * t322;
t351 = Icges(4,1) * t325 + t373;
t348 = Icges(4,2) * t327 + t374;
t345 = Icges(4,5) * t325 + Icges(4,6) * t327;
t290 = Icges(4,6) * t328 + t348 * t326;
t292 = Icges(4,5) * t328 + t351 * t326;
t338 = -t290 * t327 - t292 * t325;
t291 = Icges(4,6) * t326 - t348 * t328;
t293 = Icges(4,5) * t326 - t351 * t328;
t337 = t291 * t327 + t293 * t325;
t312 = -Icges(4,2) * t325 + t373;
t313 = Icges(4,1) * t327 - t374;
t334 = t312 * t327 + t313 * t325;
t331 = qJD(1) * t298 + (-qJD(2) - t360) * t328 + t366;
t317 = rSges(2,1) * t328 - rSges(2,2) * t326;
t316 = rSges(4,1) * t327 - rSges(4,2) * t325;
t315 = rSges(2,1) * t326 + rSges(2,2) * t328;
t311 = Icges(4,5) * t327 - Icges(4,6) * t325;
t306 = rSges(5,1) * t322 - rSges(5,2) * t321;
t296 = rSges(4,3) * t326 - t354 * t328;
t295 = rSges(4,3) * t328 + t354 * t326;
t289 = Icges(4,3) * t326 - t345 * t328;
t288 = Icges(4,3) * t328 + t345 * t326;
t286 = rSges(5,3) * t326 - t353 * t328;
t284 = rSges(5,3) * t328 + t353 * t326;
t270 = t307 - qJD(2) * t328 + qJD(1) * (-rSges(3,2) * t328 + rSges(3,3) * t326);
t269 = t320 + (rSges(3,2) * t326 + rSges(3,3) * t328 - t314) * qJD(1);
t266 = (-t295 * t326 + t296 * t328) * qJD(3);
t265 = qJD(1) * t295 + (-qJD(3) * t316 - qJD(2)) * t328 + t366;
t264 = t316 * t363 + t320 + (-t296 + t357) * qJD(1);
t263 = qJD(1) * t284 - t306 * t310 + t331;
t262 = t306 * t309 + (-t286 + t356) * qJD(1) + t365;
t261 = -t284 * t309 + t286 * t310 + t355;
t260 = t367 * qJD(1) + qJD(5) * t326 - t358 * t310 + t331;
t259 = qJD(5) * t328 + t358 * t309 + (t356 - t368) * qJD(1) + t365;
t258 = -t367 * t309 + t368 * t310 + t355;
t1 = m(3) * (t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(4) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + ((t328 * t311 + t334 * t326) * qJD(1) + (t328 ^ 2 * t288 + (t337 * t326 + (t289 - t338) * t328) * t326) * qJD(3)) * t362 / 0.2e1 + ((t326 * t311 - t334 * t328) * qJD(1) + (t326 ^ 2 * t289 + (t338 * t328 + (t288 - t337) * t326) * t328) * qJD(3)) * t363 / 0.2e1 + m(5) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(6) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + (t384 * t326 - t385 * t328) * t309 / 0.2e1 + (t385 * t326 + t384 * t328) * t310 / 0.2e1 + (m(2) * (t315 ^ 2 + t317 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1 + (((-t290 * t325 + t292 * t327) * t328 + (-t291 * t325 + t293 * t327) * t326) * qJD(3) + (-t394 * t321 + t392 * t322) * t310 + (-t393 * t321 + t391 * t322) * t309 + (-t325 * t312 + t327 * t313 - t390 * t321 + t389 * t322) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
