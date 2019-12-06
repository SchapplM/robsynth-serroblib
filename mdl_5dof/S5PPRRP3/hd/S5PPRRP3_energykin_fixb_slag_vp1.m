% Calculate kinetic energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:32
% EndTime: 2019-12-05 15:10:33
% DurationCPUTime: 1.19s
% Computational Cost: add. (851->171), mult. (2182->270), div. (0->0), fcn. (2537->8), ass. (0->87)
t418 = Icges(5,1) + Icges(6,1);
t417 = -Icges(5,4) + Icges(6,5);
t416 = Icges(6,4) + Icges(5,5);
t415 = Icges(5,2) + Icges(6,3);
t414 = Icges(6,2) + Icges(5,3);
t413 = -Icges(5,6) + Icges(6,6);
t412 = rSges(6,1) + pkin(4);
t411 = rSges(6,3) + qJ(5);
t370 = cos(pkin(8));
t371 = cos(pkin(7));
t373 = sin(qJ(3));
t391 = t371 * t373;
t369 = sin(pkin(7));
t374 = cos(qJ(3));
t393 = t369 * t374;
t357 = t370 * t393 - t391;
t372 = sin(qJ(4));
t368 = sin(pkin(8));
t397 = cos(qJ(4));
t382 = t368 * t397;
t345 = t357 * t372 - t369 * t382;
t395 = t369 * t368;
t346 = t357 * t397 + t372 * t395;
t390 = t371 * t374;
t394 = t369 * t373;
t356 = t370 * t394 + t390;
t410 = t415 * t345 + t417 * t346 + t413 * t356;
t359 = t370 * t390 + t394;
t347 = t359 * t372 - t371 * t382;
t392 = t371 * t368;
t348 = t359 * t397 + t372 * t392;
t358 = t370 * t391 - t393;
t409 = t415 * t347 + t417 * t348 + t413 * t358;
t408 = t413 * t345 + t416 * t346 + t414 * t356;
t407 = t413 * t347 + t416 * t348 + t414 * t358;
t406 = t417 * t345 + t418 * t346 + t416 * t356;
t405 = t417 * t347 + t418 * t348 + t416 * t358;
t360 = t368 * t374 * t372 + t370 * t397;
t361 = -t370 * t372 + t374 * t382;
t396 = t368 * t373;
t404 = t415 * t360 + t417 * t361 + t413 * t396;
t403 = t413 * t360 + t416 * t361 + t414 * t396;
t402 = t417 * t360 + t418 * t361 + t416 * t396;
t389 = rSges(6,2) * t356 + t411 * t345 + t412 * t346;
t388 = rSges(6,2) * t358 + t411 * t347 + t412 * t348;
t387 = rSges(6,2) * t396 + t411 * t360 + t412 * t361;
t386 = qJD(2) * t371;
t385 = qJD(3) * t368;
t384 = qJD(3) * t370;
t342 = pkin(3) * t357 + pkin(6) * t356;
t362 = (pkin(3) * t374 + pkin(6) * t373) * t368;
t367 = qJD(2) * t369;
t381 = t369 * t385;
t383 = t342 * t384 + t362 * t381 + t367;
t380 = t371 * t385;
t343 = pkin(3) * t359 + pkin(6) * t358;
t378 = t342 * t380 - t343 * t381 + qJD(1);
t377 = (-t343 * t370 - t362 * t392) * qJD(3) - t386;
t376 = qJD(1) ^ 2;
t363 = qJD(4) * t396 - t384;
t355 = -rSges(4,3) * t370 + (rSges(4,1) * t374 - rSges(4,2) * t373) * t368;
t354 = -Icges(4,5) * t370 + (Icges(4,1) * t374 - Icges(4,4) * t373) * t368;
t353 = -Icges(4,6) * t370 + (Icges(4,4) * t374 - Icges(4,2) * t373) * t368;
t352 = -Icges(4,3) * t370 + (Icges(4,5) * t374 - Icges(4,6) * t373) * t368;
t350 = qJD(4) * t358 + t380;
t349 = qJD(4) * t356 + t381;
t340 = rSges(5,1) * t361 - rSges(5,2) * t360 + rSges(5,3) * t396;
t331 = rSges(4,1) * t359 - rSges(4,2) * t358 + rSges(4,3) * t392;
t330 = rSges(4,1) * t357 - rSges(4,2) * t356 + rSges(4,3) * t395;
t329 = Icges(4,1) * t359 - Icges(4,4) * t358 + Icges(4,5) * t392;
t328 = Icges(4,1) * t357 - Icges(4,4) * t356 + Icges(4,5) * t395;
t327 = Icges(4,4) * t359 - Icges(4,2) * t358 + Icges(4,6) * t392;
t326 = Icges(4,4) * t357 - Icges(4,2) * t356 + Icges(4,6) * t395;
t325 = Icges(4,5) * t359 - Icges(4,6) * t358 + Icges(4,3) * t392;
t324 = Icges(4,5) * t357 - Icges(4,6) * t356 + Icges(4,3) * t395;
t321 = rSges(5,1) * t348 - rSges(5,2) * t347 + rSges(5,3) * t358;
t319 = rSges(5,1) * t346 - rSges(5,2) * t345 + rSges(5,3) * t356;
t305 = -t386 + (-t331 * t370 - t355 * t392) * qJD(3);
t304 = t367 + (t330 * t370 + t355 * t395) * qJD(3);
t303 = qJD(1) + (t330 * t371 - t331 * t369) * t385;
t302 = t321 * t363 - t340 * t350 + t377;
t301 = -t319 * t363 + t340 * t349 + t383;
t300 = t319 * t350 - t321 * t349 + t378;
t299 = qJD(5) * t345 - t387 * t350 + t388 * t363 + t377;
t298 = qJD(5) * t347 + t387 * t349 - t389 * t363 + t383;
t297 = qJD(5) * t360 - t388 * t349 + t389 * t350 + t378;
t1 = m(2) * t376 / 0.2e1 + m(3) * (t376 + (t369 ^ 2 + t371 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(5) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + m(6) * (t297 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + ((t404 * t345 + t402 * t346 + t403 * t356) * t363 + (t409 * t345 + t405 * t346 + t407 * t356) * t350 + (t410 * t345 + t406 * t346 + t408 * t356) * t349) * t349 / 0.2e1 + ((t404 * t347 + t402 * t348 + t403 * t358) * t363 + (t409 * t347 + t405 * t348 + t407 * t358) * t350 + (t410 * t347 + t406 * t348 + t408 * t358) * t349) * t350 / 0.2e1 + ((t404 * t360 + t402 * t361 + t403 * t396) * t363 + (t409 * t360 + t405 * t361 + t407 * t396) * t350 + (t410 * t360 + t406 * t361 + t408 * t396) * t349) * t363 / 0.2e1 + (-t370 * (t370 ^ 2 * t352 + (((-t327 * t373 + t329 * t374) * t371 + (-t326 * t373 + t328 * t374) * t369) * t368 + (-t324 * t369 - t325 * t371 + t353 * t373 - t354 * t374) * t370) * t368) / 0.2e1 + (t371 * ((t325 * t392 - t327 * t358 + t329 * t359) * t392 + (t324 * t392 - t326 * t358 + t328 * t359) * t395 - (t352 * t392 - t353 * t358 + t354 * t359) * t370) + t369 * ((t325 * t395 - t327 * t356 + t329 * t357) * t392 + (t324 * t395 - t326 * t356 + t328 * t357) * t395 - (t352 * t395 - t353 * t356 + t354 * t357) * t370)) * t368 / 0.2e1) * qJD(3) ^ 2;
T = t1;
