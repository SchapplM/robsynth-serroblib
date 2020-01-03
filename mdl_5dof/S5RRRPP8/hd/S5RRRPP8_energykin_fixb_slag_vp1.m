% Calculate kinetic energy for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:28
% EndTime: 2019-12-31 21:07:30
% DurationCPUTime: 1.64s
% Computational Cost: add. (744->182), mult. (1777->282), div. (0->0), fcn. (1812->6), ass. (0->99)
t398 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t397 = -Icges(4,4) - Icges(5,6) + Icges(6,6);
t396 = -Icges(4,5) - Icges(6,5) + Icges(5,4);
t395 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t394 = Icges(4,6) - Icges(5,5) - Icges(6,4);
t393 = -Icges(4,3) - Icges(6,1) - Icges(5,1);
t392 = rSges(6,1) + pkin(4);
t391 = rSges(6,3) + qJ(5);
t343 = sin(qJ(3));
t346 = cos(qJ(3));
t348 = cos(qJ(1));
t371 = t346 * t348;
t345 = sin(qJ(1));
t347 = cos(qJ(2));
t372 = t345 * t347;
t320 = t343 * t372 + t371;
t370 = t348 * t343;
t321 = t346 * t372 - t370;
t344 = sin(qJ(2));
t375 = t344 * t345;
t390 = t397 * t320 + t398 * t321 - t396 * t375;
t322 = -t345 * t346 + t347 * t370;
t323 = t343 * t345 + t347 * t371;
t373 = t344 * t348;
t389 = t397 * t322 + t398 * t323 - t396 * t373;
t388 = t395 * t320 + t397 * t321 - t394 * t375;
t387 = t395 * t322 + t397 * t323 - t394 * t373;
t386 = -t394 * t320 - t396 * t321 - t393 * t375;
t385 = -t394 * t322 - t396 * t323 - t393 * t373;
t384 = t393 * t347 + (-t394 * t343 - t396 * t346) * t344;
t383 = t394 * t347 + (t395 * t343 + t397 * t346) * t344;
t382 = t396 * t347 + (t397 * t343 + t398 * t346) * t344;
t377 = Icges(3,4) * t344;
t376 = Icges(3,4) * t347;
t374 = t344 * t346;
t369 = rSges(6,2) * t320 + t391 * t321 + t392 * t375;
t368 = rSges(6,2) * t322 + t391 * t323 + t392 * t373;
t367 = (rSges(6,2) * t343 + rSges(6,3) * t346) * t344 + qJ(5) * t374 - t392 * t347;
t361 = pkin(2) * t347 + pkin(7) * t344;
t325 = t361 * t345;
t326 = t361 * t348;
t364 = qJD(2) * t348;
t365 = qJD(2) * t345;
t366 = t325 * t365 + t326 * t364;
t363 = qJD(3) * t344;
t290 = pkin(3) * t321 + qJ(4) * t320;
t327 = t348 * t363 + t365;
t362 = qJD(4) * t344 * t343 + t327 * t290 + t366;
t360 = rSges(3,1) * t347 - rSges(3,2) * t344;
t359 = Icges(3,1) * t347 - t377;
t358 = -Icges(3,2) * t344 + t376;
t357 = Icges(3,5) * t347 - Icges(3,6) * t344;
t299 = -Icges(3,6) * t348 + t345 * t358;
t302 = -Icges(3,5) * t348 + t345 * t359;
t356 = t299 * t344 - t302 * t347;
t300 = Icges(3,6) * t345 + t348 * t358;
t303 = Icges(3,5) * t345 + t348 * t359;
t355 = -t300 * t344 + t303 * t347;
t332 = Icges(3,2) * t347 + t377;
t333 = Icges(3,1) * t344 + t376;
t354 = -t332 * t344 + t333 * t347;
t330 = qJD(1) * (pkin(1) * t348 + pkin(6) * t345);
t337 = pkin(2) * t344 - pkin(7) * t347;
t353 = qJD(1) * t326 - t337 * t365 + t330;
t291 = pkin(3) * t323 + qJ(4) * t322;
t340 = -qJD(3) * t347 + qJD(1);
t352 = qJD(4) * t320 + t340 * t291 + t353;
t338 = pkin(1) * t345 - pkin(6) * t348;
t351 = (-t325 - t338) * qJD(1) - t337 * t364;
t324 = (pkin(3) * t346 + qJ(4) * t343) * t344;
t328 = t345 * t363 - t364;
t350 = qJD(4) * t322 + t328 * t324 + t351;
t336 = rSges(2,1) * t348 - rSges(2,2) * t345;
t335 = rSges(2,1) * t345 + rSges(2,2) * t348;
t334 = rSges(3,1) * t344 + rSges(3,2) * t347;
t331 = Icges(3,5) * t344 + Icges(3,6) * t347;
t314 = -rSges(5,1) * t347 + (-rSges(5,2) * t346 + rSges(5,3) * t343) * t344;
t312 = rSges(3,3) * t345 + t348 * t360;
t311 = -rSges(3,3) * t348 + t345 * t360;
t310 = -rSges(4,3) * t347 + (rSges(4,1) * t346 - rSges(4,2) * t343) * t344;
t297 = Icges(3,3) * t345 + t348 * t357;
t296 = -Icges(3,3) * t348 + t345 * t357;
t289 = rSges(4,1) * t323 - rSges(4,2) * t322 + rSges(4,3) * t373;
t288 = rSges(4,1) * t321 - rSges(4,2) * t320 + rSges(4,3) * t375;
t287 = rSges(5,1) * t373 - rSges(5,2) * t323 + rSges(5,3) * t322;
t285 = rSges(5,1) * t375 - rSges(5,2) * t321 + rSges(5,3) * t320;
t264 = qJD(1) * t312 - t334 * t365 + t330;
t263 = -t334 * t364 + (-t311 - t338) * qJD(1);
t261 = (t311 * t345 + t312 * t348) * qJD(2);
t260 = t289 * t340 - t310 * t327 + t353;
t259 = -t288 * t340 + t310 * t328 + t351;
t258 = t288 * t327 - t289 * t328 + t366;
t257 = t287 * t340 + (-t314 - t324) * t327 + t352;
t256 = t314 * t328 + (-t285 - t290) * t340 + t350;
t255 = t285 * t327 + (-t287 - t291) * t328 + t362;
t254 = qJD(5) * t321 + t368 * t340 + (-t324 - t367) * t327 + t352;
t253 = qJD(5) * t323 + t367 * t328 + (-t290 - t369) * t340 + t350;
t252 = qJD(5) * t374 + t369 * t327 + (-t291 - t368) * t328 + t362;
t1 = m(3) * (t261 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + ((t345 * t331 + t348 * t354) * qJD(1) + (t345 ^ 2 * t297 + (t356 * t348 + (-t296 + t355) * t345) * t348) * qJD(2)) * t365 / 0.2e1 - ((-t348 * t331 + t345 * t354) * qJD(1) + (t348 ^ 2 * t296 + (t355 * t345 + (-t297 + t356) * t348) * t345) * qJD(2)) * t364 / 0.2e1 + qJD(1) * ((t347 * t332 + t344 * t333) * qJD(1) + ((t300 * t347 + t303 * t344) * t345 - (t299 * t347 + t302 * t344) * t348) * qJD(2)) / 0.2e1 + m(4) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(5) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + m(6) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + (m(2) * (t335 ^ 2 + t336 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t383 * t322 + t382 * t323 + t384 * t373) * t340 + (t388 * t322 + t390 * t323 + t386 * t373) * t328 + (t387 * t322 + t389 * t323 + t385 * t373) * t327) * t327 / 0.2e1 + ((t383 * t320 + t382 * t321 + t384 * t375) * t340 + (t388 * t320 + t390 * t321 + t386 * t375) * t328 + (t387 * t320 + t389 * t321 + t385 * t375) * t327) * t328 / 0.2e1 + ((-t385 * t327 - t386 * t328 - t384 * t340) * t347 + ((t383 * t343 + t382 * t346) * t340 + (t388 * t343 + t390 * t346) * t328 + (t387 * t343 + t389 * t346) * t327) * t344) * t340 / 0.2e1;
T = t1;
