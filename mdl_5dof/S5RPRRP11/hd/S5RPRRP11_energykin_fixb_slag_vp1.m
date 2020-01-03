% Calculate kinetic energy for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:25
% EndTime: 2019-12-31 18:53:26
% DurationCPUTime: 1.32s
% Computational Cost: add. (966->171), mult. (1268->272), div. (0->0), fcn. (1260->8), ass. (0->97)
t401 = Icges(5,1) + Icges(6,1);
t400 = -Icges(5,4) + Icges(6,5);
t399 = Icges(6,4) + Icges(5,5);
t398 = Icges(5,2) + Icges(6,3);
t397 = -Icges(6,6) + Icges(5,6);
t396 = -Icges(5,3) - Icges(6,2);
t395 = rSges(6,1) + pkin(4);
t394 = rSges(6,3) + qJ(5);
t342 = pkin(8) + qJ(3);
t340 = cos(t342);
t348 = cos(qJ(4));
t349 = cos(qJ(1));
t372 = t348 * t349;
t346 = sin(qJ(4));
t347 = sin(qJ(1));
t375 = t346 * t347;
t321 = t340 * t375 + t372;
t373 = t347 * t348;
t374 = t346 * t349;
t322 = t340 * t373 - t374;
t339 = sin(t342);
t377 = t339 * t347;
t393 = t398 * t321 + t400 * t322 - t397 * t377;
t323 = t340 * t374 - t373;
t324 = t340 * t372 + t375;
t376 = t339 * t349;
t392 = t398 * t323 + t400 * t324 - t397 * t376;
t391 = -t397 * t321 + t399 * t322 - t396 * t377;
t390 = -t397 * t323 + t399 * t324 - t396 * t376;
t389 = t400 * t321 + t401 * t322 + t399 * t377;
t388 = t400 * t323 + t401 * t324 + t399 * t376;
t387 = t397 * t340 + (t398 * t346 + t400 * t348) * t339;
t386 = t396 * t340 + (-t397 * t346 + t399 * t348) * t339;
t385 = -t399 * t340 + (t400 * t346 + t401 * t348) * t339;
t344 = cos(pkin(8));
t380 = pkin(2) * t344;
t379 = Icges(4,4) * t339;
t378 = Icges(4,4) * t340;
t370 = rSges(6,2) * t377 + t394 * t321 + t395 * t322;
t369 = rSges(6,2) * t376 + t394 * t323 + t395 * t324;
t368 = -rSges(6,2) * t340 + (t394 * t346 + t395 * t348) * t339;
t333 = pkin(1) * t347 - qJ(2) * t349;
t367 = pkin(6) * t349 - t347 * t380 - t333;
t362 = pkin(3) * t340 + pkin(7) * t339;
t319 = t362 * t347;
t320 = t362 * t349;
t364 = qJD(3) * t349;
t365 = qJD(3) * t347;
t366 = t319 * t365 + t320 * t364;
t363 = qJD(4) * t339;
t332 = qJD(1) * (pkin(1) * t349 + qJ(2) * t347);
t361 = -qJD(2) * t349 + qJD(1) * (pkin(6) * t347 + t349 * t380) + t332;
t343 = sin(pkin(8));
t360 = rSges(3,1) * t344 - rSges(3,2) * t343;
t359 = rSges(4,1) * t340 - rSges(4,2) * t339;
t358 = Icges(4,1) * t340 - t379;
t357 = -Icges(4,2) * t339 + t378;
t356 = Icges(4,5) * t340 - Icges(4,6) * t339;
t309 = -Icges(4,6) * t349 + t347 * t357;
t311 = -Icges(4,5) * t349 + t347 * t358;
t355 = t309 * t339 - t311 * t340;
t310 = Icges(4,6) * t347 + t349 * t357;
t312 = Icges(4,5) * t347 + t349 * t358;
t354 = -t310 * t339 + t312 * t340;
t328 = Icges(4,2) * t340 + t379;
t329 = Icges(4,1) * t339 + t378;
t353 = -t328 * t339 + t329 * t340;
t331 = pkin(3) * t339 - pkin(7) * t340;
t352 = qJD(1) * t320 - t331 * t365 + t361;
t341 = qJD(2) * t347;
t351 = t341 + (-t319 + t367) * qJD(1) - t331 * t364;
t336 = -qJD(4) * t340 + qJD(1);
t335 = rSges(2,1) * t349 - rSges(2,2) * t347;
t334 = rSges(2,1) * t347 + rSges(2,2) * t349;
t330 = rSges(4,1) * t339 + rSges(4,2) * t340;
t327 = Icges(4,5) * t339 + Icges(4,6) * t340;
t326 = t347 * t363 - t364;
t325 = t349 * t363 + t365;
t314 = rSges(4,3) * t347 + t349 * t359;
t313 = -rSges(4,3) * t349 + t347 * t359;
t308 = Icges(4,3) * t347 + t349 * t356;
t307 = -Icges(4,3) * t349 + t347 * t356;
t305 = -rSges(5,3) * t340 + (rSges(5,1) * t348 - rSges(5,2) * t346) * t339;
t296 = qJD(1) * t347 * rSges(3,3) + t332 + (qJD(1) * t360 - qJD(2)) * t349;
t295 = t341 + (t349 * rSges(3,3) - t347 * t360 - t333) * qJD(1);
t292 = rSges(5,1) * t324 - rSges(5,2) * t323 + rSges(5,3) * t376;
t290 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t377;
t276 = (t313 * t347 + t314 * t349) * qJD(3);
t275 = qJD(1) * t314 - t330 * t365 + t361;
t274 = -t330 * t364 + t341 + (-t313 + t367) * qJD(1);
t273 = t290 * t325 - t292 * t326 + t366;
t272 = t292 * t336 - t305 * t325 + t352;
t271 = -t290 * t336 + t305 * t326 + t351;
t270 = qJD(5) * t321 - t325 * t368 + t336 * t369 + t352;
t269 = qJD(5) * t323 + t326 * t368 - t336 * t370 + t351;
t268 = qJD(5) * t339 * t346 + t325 * t370 - t326 * t369 + t366;
t1 = m(3) * (t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(4) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t347 * t327 + t349 * t353) * qJD(1) + (t347 ^ 2 * t308 + (t355 * t349 + (-t307 + t354) * t347) * t349) * qJD(3)) * t365 / 0.2e1 - ((-t349 * t327 + t347 * t353) * qJD(1) + (t349 ^ 2 * t307 + (t354 * t347 + (-t308 + t355) * t349) * t347) * qJD(3)) * t364 / 0.2e1 + qJD(1) * ((t340 * t328 + t339 * t329) * qJD(1) + ((t310 * t340 + t312 * t339) * t347 - (t309 * t340 + t311 * t339) * t349) * qJD(3)) / 0.2e1 + m(5) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + m(6) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + ((t387 * t323 + t385 * t324 + t386 * t376) * t336 + (t393 * t323 + t389 * t324 + t391 * t376) * t326 + (t392 * t323 + t388 * t324 + t390 * t376) * t325) * t325 / 0.2e1 + ((t387 * t321 + t385 * t322 + t386 * t377) * t336 + (t393 * t321 + t389 * t322 + t391 * t377) * t326 + (t392 * t321 + t388 * t322 + t390 * t377) * t325) * t326 / 0.2e1 + ((-t390 * t325 - t391 * t326 - t386 * t336) * t340 + ((t387 * t346 + t385 * t348) * t336 + (t393 * t346 + t389 * t348) * t326 + (t392 * t346 + t388 * t348) * t325) * t339) * t336 / 0.2e1 + (m(2) * (t334 ^ 2 + t335 ^ 2) + Icges(2,3) + Icges(3,2) * t344 ^ 2 + (Icges(3,1) * t343 + 0.2e1 * Icges(3,4) * t344) * t343) * qJD(1) ^ 2 / 0.2e1;
T = t1;
