% Calculate kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:50
% EndTime: 2019-12-05 18:23:52
% DurationCPUTime: 2.24s
% Computational Cost: add. (775->201), mult. (1217->329), div. (0->0), fcn. (1122->8), ass. (0->117)
t415 = Icges(3,4) + Icges(4,4);
t414 = Icges(3,1) + Icges(4,1);
t413 = Icges(3,2) + Icges(4,2);
t336 = cos(qJ(2));
t412 = t415 * t336;
t333 = sin(qJ(2));
t411 = t415 * t333;
t410 = Icges(3,5) + Icges(4,5);
t409 = Icges(3,6) + Icges(4,6);
t408 = -t413 * t333 + t412;
t407 = t414 * t336 - t411;
t406 = Icges(3,3) + Icges(4,3);
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t405 = t408 * t334 - t409 * t337;
t404 = t409 * t334 + t408 * t337;
t403 = -t407 * t334 + t410 * t337;
t402 = t410 * t334 + t407 * t337;
t401 = t413 * t336 + t411;
t400 = t414 * t333 + t412;
t399 = -t409 * t333 + t410 * t336;
t398 = t399 * t334 - t406 * t337;
t397 = t406 * t334 + t399 * t337;
t396 = t410 * t333 + t409 * t336;
t395 = -t401 * t333 + t400 * t336;
t394 = -t404 * t333 + t402 * t336;
t393 = t405 * t333 + t403 * t336;
t389 = pkin(1) * t333;
t388 = pkin(1) * t336;
t330 = qJ(2) + qJ(4);
t328 = sin(t330);
t381 = Icges(5,4) * t328;
t329 = cos(t330);
t380 = Icges(5,4) * t329;
t379 = t328 * t334;
t378 = t328 * t337;
t332 = sin(qJ(5));
t377 = t332 * t334;
t376 = t332 * t337;
t335 = cos(qJ(5));
t375 = t334 * t335;
t374 = t335 * t337;
t368 = pkin(2) * t336;
t266 = -pkin(3) * t337 + t334 * t368;
t309 = -qJ(3) * t337 + t334 * t388;
t372 = -t266 - t309;
t310 = qJ(3) * t334 + t337 * t388;
t327 = qJD(2) * t334;
t370 = qJD(2) * t337;
t371 = t309 * t327 + t310 * t370;
t312 = qJD(4) * t334 + t327;
t369 = qJD(5) * t328;
t365 = qJD(1) * t310 - qJD(3) * t337;
t313 = (-qJD(2) - qJD(4)) * t337;
t267 = pkin(3) * t334 + t337 * t368;
t364 = t266 * t327 + t267 * t370 + t371;
t363 = rSges(3,1) * t336 - rSges(3,2) * t333;
t362 = rSges(4,1) * t336 - rSges(4,2) * t333;
t361 = rSges(5,1) * t329 - rSges(5,2) * t328;
t360 = qJD(2) * (-pkin(2) * t333 - t389);
t359 = qJD(2) * (-rSges(4,1) * t333 - rSges(4,2) * t336 - t389);
t356 = Icges(5,1) * t329 - t381;
t353 = -Icges(5,2) * t328 + t380;
t350 = Icges(5,5) * t329 - Icges(5,6) * t328;
t326 = qJD(3) * t334;
t343 = t337 * t360 + t326;
t342 = (Icges(5,5) * t328 + Icges(5,6) * t329) * qJD(1) + (-Icges(5,3) * t337 + t334 * t350) * t313 + (Icges(5,3) * t334 + t337 * t350) * t312;
t341 = qJD(1) * t267 + t334 * t360 + t365;
t274 = -Icges(5,6) * t337 + t334 * t353;
t275 = Icges(5,6) * t334 + t337 * t353;
t276 = -Icges(5,5) * t337 + t356 * t334;
t277 = Icges(5,5) * t334 + t356 * t337;
t306 = Icges(5,2) * t329 + t381;
t307 = Icges(5,1) * t328 + t380;
t340 = (-t275 * t328 + t277 * t329) * t312 + (-t274 * t328 + t276 * t329) * t313 + (-t306 * t328 + t307 * t329) * qJD(1);
t324 = -qJD(5) * t329 + qJD(1);
t323 = rSges(2,1) * t337 - rSges(2,2) * t334;
t322 = rSges(2,1) * t334 + rSges(2,2) * t337;
t321 = rSges(3,1) * t333 + rSges(3,2) * t336;
t308 = rSges(5,1) * t328 + rSges(5,2) * t329;
t301 = t329 * t374 + t377;
t300 = -t329 * t376 + t375;
t299 = t329 * t375 - t376;
t298 = -t329 * t377 - t374;
t297 = rSges(3,3) * t334 + t363 * t337;
t296 = rSges(4,3) * t334 + t362 * t337;
t295 = -rSges(3,3) * t337 + t363 * t334;
t294 = -rSges(4,3) * t337 + t362 * t334;
t293 = t334 * t369 + t313;
t292 = t337 * t369 + t312;
t279 = rSges(5,3) * t334 + t361 * t337;
t278 = -rSges(5,3) * t337 + t361 * t334;
t271 = -rSges(6,3) * t329 + (rSges(6,1) * t335 - rSges(6,2) * t332) * t328;
t270 = -Icges(6,5) * t329 + (Icges(6,1) * t335 - Icges(6,4) * t332) * t328;
t269 = -Icges(6,6) * t329 + (Icges(6,4) * t335 - Icges(6,2) * t332) * t328;
t268 = -Icges(6,3) * t329 + (Icges(6,5) * t335 - Icges(6,6) * t332) * t328;
t262 = -qJD(1) * t295 - t321 * t370;
t261 = qJD(1) * t297 - t321 * t327;
t260 = rSges(6,1) * t301 + rSges(6,2) * t300 + rSges(6,3) * t378;
t259 = rSges(6,1) * t299 + rSges(6,2) * t298 + rSges(6,3) * t379;
t258 = Icges(6,1) * t301 + Icges(6,4) * t300 + Icges(6,5) * t378;
t257 = Icges(6,1) * t299 + Icges(6,4) * t298 + Icges(6,5) * t379;
t256 = Icges(6,4) * t301 + Icges(6,2) * t300 + Icges(6,6) * t378;
t255 = Icges(6,4) * t299 + Icges(6,2) * t298 + Icges(6,6) * t379;
t254 = Icges(6,5) * t301 + Icges(6,6) * t300 + Icges(6,3) * t378;
t253 = Icges(6,5) * t299 + Icges(6,6) * t298 + Icges(6,3) * t379;
t252 = (t295 * t334 + t297 * t337) * qJD(2);
t251 = t326 + t337 * t359 + (-t294 - t309) * qJD(1);
t250 = qJD(1) * t296 + t334 * t359 + t365;
t249 = (t294 * t334 + t296 * t337) * qJD(2) + t371;
t248 = t308 * t313 + (-t278 + t372) * qJD(1) + t343;
t247 = qJD(1) * t279 - t308 * t312 + t341;
t246 = t278 * t312 - t279 * t313 + t364;
t245 = -pkin(4) * t313 * t329 - t259 * t324 + t271 * t293 + (-pkin(4) * t379 + t372) * qJD(1) + t343;
t244 = t260 * t324 - t271 * t292 + (qJD(1) * t378 + t329 * t312) * pkin(4) + t341;
t243 = t259 * t292 - t260 * t293 + (t334 * t312 - t337 * t313) * t328 * pkin(4) + t364;
t1 = m(3) * (t252 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + m(4) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(5) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t312 * (t334 * t342 + t337 * t340) / 0.2e1 + t313 * (t334 * t340 - t337 * t342) / 0.2e1 + m(6) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + t292 * ((t254 * t378 + t256 * t300 + t258 * t301) * t292 + (t253 * t378 + t255 * t300 + t257 * t301) * t293 + (t268 * t378 + t269 * t300 + t270 * t301) * t324) / 0.2e1 + t293 * ((t254 * t379 + t256 * t298 + t258 * t299) * t292 + (t253 * t379 + t255 * t298 + t257 * t299) * t293 + (t268 * t379 + t269 * t298 + t270 * t299) * t324) / 0.2e1 + t324 * ((-t253 * t293 - t254 * t292 - t268 * t324) * t329 + ((-t256 * t332 + t258 * t335) * t292 + (-t255 * t332 + t257 * t335) * t293 + (-t269 * t332 + t270 * t335) * t324) * t328) / 0.2e1 + (m(2) * (t322 ^ 2 + t323 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t397 * t334 ^ 2 + (t393 * t337 + (t394 - t398) * t334) * t337) * qJD(2) + (t334 * t396 + t337 * t395) * qJD(1)) * t327 / 0.2e1 - ((t398 * t337 ^ 2 + (t394 * t334 + (t393 - t397) * t337) * t334) * qJD(2) + (t334 * t395 - t396 * t337) * qJD(1)) * t370 / 0.2e1 + ((t275 * t329 + t277 * t328) * t312 + (t274 * t329 + t276 * t328) * t313 + ((t403 * t333 - t405 * t336) * t337 + (t402 * t333 + t404 * t336) * t334) * qJD(2) + (t329 * t306 + t328 * t307 + t400 * t333 + t401 * t336) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
