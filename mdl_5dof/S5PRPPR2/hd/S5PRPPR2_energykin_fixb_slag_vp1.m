% Calculate kinetic energy for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:56
% EndTime: 2019-12-05 15:23:57
% DurationCPUTime: 1.25s
% Computational Cost: add. (865->152), mult. (1040->253), div. (0->0), fcn. (1024->10), ass. (0->87)
t399 = Icges(3,3) + Icges(4,3);
t336 = qJ(2) + pkin(8);
t332 = sin(t336);
t334 = cos(t336);
t343 = sin(qJ(2));
t344 = cos(qJ(2));
t398 = Icges(3,5) * t344 + Icges(4,5) * t334 - Icges(3,6) * t343 - Icges(4,6) * t332;
t340 = cos(pkin(7));
t390 = t340 ^ 2;
t338 = sin(pkin(7));
t391 = t338 ^ 2;
t393 = t390 + t391;
t392 = qJD(2) * t393;
t397 = t398 * t338 - t399 * t340;
t396 = t399 * t338 + t398 * t340;
t389 = qJD(2) ^ 2;
t388 = pkin(2) * t343;
t387 = pkin(2) * t344;
t339 = cos(pkin(9));
t386 = pkin(4) * t339;
t384 = t332 * t338;
t383 = t332 * t340;
t382 = t334 * t338;
t381 = t334 * t340;
t337 = sin(pkin(9));
t380 = t337 * t338;
t379 = t337 * t340;
t378 = t338 * t339;
t377 = t339 * t340;
t330 = qJD(3) * t338;
t371 = qJD(4) * t332;
t375 = t340 * t371 + t330;
t374 = qJD(2) * t338;
t373 = qJD(2) * t340;
t372 = qJD(3) * t340;
t370 = qJD(5) * t332;
t369 = qJD(5) * t334;
t368 = qJD(1) + (-qJ(3) * t340 + t338 * t387) * t374 + (qJ(3) * t338 + t340 * t387) * t373;
t365 = -pkin(3) * t332 + qJ(4) * t334 - t388;
t364 = t338 * t371 - t372;
t363 = qJD(2) * (-rSges(4,1) * t332 - rSges(4,2) * t334 - t388);
t352 = qJD(2) * (pkin(6) * t334 - t332 * t386 + t365);
t351 = qJD(2) * (rSges(5,3) * t334 - (rSges(5,1) * t339 - rSges(5,2) * t337) * t332 + t365);
t347 = -qJD(4) * t334 + t368 + (pkin(3) * t334 + qJ(4) * t332) * t392;
t335 = pkin(9) + qJ(5);
t333 = cos(t335);
t331 = sin(t335);
t327 = rSges(3,1) * t343 + rSges(3,2) * t344;
t321 = t338 * t370 - t373;
t320 = t340 * t370 + t374;
t319 = t334 * t377 + t380;
t318 = -t334 * t379 + t378;
t317 = t334 * t378 - t379;
t316 = -t334 * t380 - t377;
t309 = t331 * t338 + t333 * t381;
t308 = -t331 * t381 + t333 * t338;
t307 = -t331 * t340 + t333 * t382;
t306 = -t331 * t382 - t333 * t340;
t296 = -rSges(6,3) * t334 + (rSges(6,1) * t333 - rSges(6,2) * t331) * t332;
t293 = -Icges(6,5) * t334 + (Icges(6,1) * t333 - Icges(6,4) * t331) * t332;
t292 = -Icges(6,6) * t334 + (Icges(6,4) * t333 - Icges(6,2) * t331) * t332;
t291 = -Icges(6,3) * t334 + (Icges(6,5) * t333 - Icges(6,6) * t331) * t332;
t290 = t340 * t363 + t330;
t289 = t338 * t363 - t372;
t287 = Icges(5,1) * t319 + Icges(5,4) * t318 + Icges(5,5) * t383;
t286 = Icges(5,1) * t317 + Icges(5,4) * t316 + Icges(5,5) * t384;
t285 = Icges(5,4) * t319 + Icges(5,2) * t318 + Icges(5,6) * t383;
t284 = Icges(5,4) * t317 + Icges(5,2) * t316 + Icges(5,6) * t384;
t283 = Icges(5,5) * t319 + Icges(5,6) * t318 + Icges(5,3) * t383;
t282 = Icges(5,5) * t317 + Icges(5,6) * t316 + Icges(5,3) * t384;
t281 = qJD(1) + (rSges(3,1) * t344 - rSges(3,2) * t343) * t392;
t280 = rSges(6,1) * t309 + rSges(6,2) * t308 + rSges(6,3) * t383;
t279 = rSges(6,1) * t307 + rSges(6,2) * t306 + rSges(6,3) * t384;
t278 = Icges(6,1) * t309 + Icges(6,4) * t308 + Icges(6,5) * t383;
t277 = Icges(6,1) * t307 + Icges(6,4) * t306 + Icges(6,5) * t384;
t276 = Icges(6,4) * t309 + Icges(6,2) * t308 + Icges(6,6) * t383;
t275 = Icges(6,4) * t307 + Icges(6,2) * t306 + Icges(6,6) * t384;
t274 = Icges(6,5) * t309 + Icges(6,6) * t308 + Icges(6,3) * t383;
t273 = Icges(6,5) * t307 + Icges(6,6) * t306 + Icges(6,3) * t384;
t272 = t340 * t351 + t375;
t271 = t338 * t351 + t364;
t270 = t368 + (rSges(4,1) * t334 - rSges(4,2) * t332) * t392;
t269 = t279 * t369 + t296 * t321 + t340 * t352 + t375;
t268 = -t280 * t369 - t296 * t320 + t338 * t352 + t364;
t267 = (t338 * (rSges(5,1) * t317 + rSges(5,2) * t316 + rSges(5,3) * t384) + t340 * (rSges(5,1) * t319 + rSges(5,2) * t318 + rSges(5,3) * t383)) * qJD(2) + t347;
t266 = t320 * t279 - t321 * t280 + t347 + (pkin(6) * t332 + t334 * t386) * t392;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t393 * t389 * t327 ^ 2 + t281 ^ 2) / 0.2e1 + m(4) * (t270 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + m(5) * (t267 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + m(6) * (t266 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + t320 * ((t274 * t383 + t276 * t308 + t278 * t309) * t320 + (t273 * t383 + t275 * t308 + t277 * t309) * t321 - (t291 * t383 + t292 * t308 + t293 * t309) * t369) / 0.2e1 + t321 * ((t274 * t384 + t276 * t306 + t278 * t307) * t320 + (t273 * t384 + t275 * t306 + t277 * t307) * t321 - (t291 * t384 + t292 * t306 + t293 * t307) * t369) / 0.2e1 - ((-t273 * t321 - t274 * t320 + t291 * t369) * t334 + ((-t276 * t331 + t278 * t333) * t320 + (-t275 * t331 + t277 * t333) * t321 - (-t292 * t331 + t293 * t333) * t369) * t332) * t369 / 0.2e1 + ((t283 * t383 + t285 * t318 + t287 * t319) * t338 + t396 * t391 + (-t282 * t383 - t284 * t318 - t286 * t319 - t397 * t338) * t340) * t338 * t389 / 0.2e1 - (-(t282 * t384 + t284 * t316 + t286 * t317) * t340 + t397 * t390 + (t283 * t384 + t285 * t316 + t287 * t317 - t396 * t340) * t338) * t340 * t389 / 0.2e1;
T = t1;
