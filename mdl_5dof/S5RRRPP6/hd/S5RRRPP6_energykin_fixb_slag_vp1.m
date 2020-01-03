% Calculate kinetic energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:57
% EndTime: 2019-12-31 20:59:59
% DurationCPUTime: 1.92s
% Computational Cost: add. (1089->223), mult. (1805->336), div. (0->0), fcn. (1824->8), ass. (0->115)
t426 = Icges(5,1) + Icges(6,1);
t425 = -Icges(5,4) + Icges(6,5);
t424 = Icges(6,4) + Icges(5,5);
t423 = Icges(5,2) + Icges(6,3);
t422 = -Icges(6,6) + Icges(5,6);
t421 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t420 = rSges(6,1) + pkin(4);
t419 = rSges(6,3) + qJ(5);
t365 = qJ(3) + pkin(8);
t363 = sin(t365);
t364 = cos(t365);
t372 = cos(qJ(1));
t369 = sin(qJ(1));
t371 = cos(qJ(2));
t397 = t369 * t371;
t321 = t363 * t397 + t364 * t372;
t322 = -t363 * t372 + t364 * t397;
t368 = sin(qJ(2));
t399 = t368 * t369;
t418 = t423 * t321 + t425 * t322 - t422 * t399;
t396 = t371 * t372;
t323 = t363 * t396 - t369 * t364;
t324 = t363 * t369 + t364 * t396;
t398 = t368 * t372;
t417 = t423 * t323 + t425 * t324 - t422 * t398;
t416 = t425 * t321 + t426 * t322 + t424 * t399;
t415 = t425 * t323 + t426 * t324 + t424 * t398;
t414 = t422 * t371 + (t423 * t363 + t425 * t364) * t368;
t413 = -t424 * t371 + (t425 * t363 + t426 * t364) * t368;
t367 = sin(qJ(3));
t370 = cos(qJ(3));
t341 = -t367 * t397 - t370 * t372;
t400 = t367 * t372;
t342 = t370 * t397 - t400;
t412 = Icges(4,5) * t342 + Icges(4,6) * t341 - t422 * t321 + t424 * t322 - t421 * t399;
t343 = -t367 * t396 + t369 * t370;
t401 = t367 * t369;
t344 = t370 * t396 + t401;
t411 = Icges(4,5) * t344 + Icges(4,6) * t343 - t422 * t323 + t424 * t324 - t421 * t398;
t410 = t421 * t371 + (Icges(4,5) * t370 - Icges(4,6) * t367 - t422 * t363 + t424 * t364) * t368;
t405 = pkin(3) * t370;
t403 = Icges(3,4) * t368;
t402 = Icges(3,4) * t371;
t395 = rSges(6,2) * t399 + t419 * t321 + t420 * t322;
t394 = rSges(6,2) * t398 + t419 * t323 + t420 * t324;
t393 = -rSges(6,2) * t371 + (t419 * t363 + t420 * t364) * t368;
t387 = pkin(2) * t371 + pkin(7) * t368;
t345 = t387 * t369;
t346 = t387 * t372;
t390 = qJD(2) * t372;
t391 = qJD(2) * t369;
t392 = t345 * t391 + t346 * t390;
t389 = qJD(3) * t368;
t388 = qJD(4) * t368;
t386 = rSges(3,1) * t371 - rSges(3,2) * t368;
t385 = Icges(3,1) * t371 - t403;
t384 = -Icges(3,2) * t368 + t402;
t383 = Icges(3,5) * t371 - Icges(3,6) * t368;
t329 = -Icges(3,6) * t372 + t369 * t384;
t332 = -Icges(3,5) * t372 + t369 * t385;
t382 = t329 * t368 - t332 * t371;
t330 = Icges(3,6) * t369 + t372 * t384;
t333 = Icges(3,5) * t369 + t372 * t385;
t381 = -t330 * t368 + t333 * t371;
t351 = Icges(3,2) * t371 + t403;
t352 = Icges(3,1) * t368 + t402;
t380 = -t351 * t368 + t352 * t371;
t349 = qJD(1) * (pkin(1) * t372 + pkin(6) * t369);
t356 = pkin(2) * t368 - pkin(7) * t371;
t379 = qJD(1) * t346 - t356 * t391 + t349;
t377 = qJ(4) * t368 + t371 * t405;
t307 = -pkin(3) * t400 + t369 * t377;
t347 = t372 * t389 + t391;
t378 = -qJD(4) * t371 + t347 * t307 + t392;
t308 = pkin(3) * t401 + t372 * t377;
t361 = -qJD(3) * t371 + qJD(1);
t376 = t361 * t308 + t369 * t388 + t379;
t357 = pkin(1) * t369 - pkin(6) * t372;
t375 = (-t345 - t357) * qJD(1) - t356 * t390;
t312 = -qJ(4) * t371 + t368 * t405;
t348 = t369 * t389 - t390;
t374 = t348 * t312 + t372 * t388 + t375;
t355 = rSges(2,1) * t372 - rSges(2,2) * t369;
t354 = rSges(2,1) * t369 + rSges(2,2) * t372;
t353 = rSges(3,1) * t368 + rSges(3,2) * t371;
t350 = Icges(3,5) * t368 + Icges(3,6) * t371;
t337 = rSges(3,3) * t369 + t372 * t386;
t336 = -rSges(3,3) * t372 + t369 * t386;
t335 = -rSges(4,3) * t371 + (rSges(4,1) * t370 - rSges(4,2) * t367) * t368;
t331 = -Icges(4,5) * t371 + (Icges(4,1) * t370 - Icges(4,4) * t367) * t368;
t328 = -Icges(4,6) * t371 + (Icges(4,4) * t370 - Icges(4,2) * t367) * t368;
t327 = Icges(3,3) * t369 + t372 * t383;
t326 = -Icges(3,3) * t372 + t369 * t383;
t320 = -rSges(5,3) * t371 + (rSges(5,1) * t364 - rSges(5,2) * t363) * t368;
t310 = rSges(4,1) * t344 + rSges(4,2) * t343 + rSges(4,3) * t398;
t309 = rSges(4,1) * t342 + rSges(4,2) * t341 + rSges(4,3) * t399;
t306 = Icges(4,1) * t344 + Icges(4,4) * t343 + Icges(4,5) * t398;
t305 = Icges(4,1) * t342 + Icges(4,4) * t341 + Icges(4,5) * t399;
t304 = Icges(4,4) * t344 + Icges(4,2) * t343 + Icges(4,6) * t398;
t303 = Icges(4,4) * t342 + Icges(4,2) * t341 + Icges(4,6) * t399;
t298 = qJD(1) * t337 - t353 * t391 + t349;
t297 = -t353 * t390 + (-t336 - t357) * qJD(1);
t296 = (t336 * t369 + t337 * t372) * qJD(2);
t294 = rSges(5,1) * t324 - rSges(5,2) * t323 + rSges(5,3) * t398;
t292 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t399;
t277 = t310 * t361 - t335 * t347 + t379;
t276 = -t309 * t361 + t335 * t348 + t375;
t275 = t309 * t347 - t310 * t348 + t392;
t274 = t294 * t361 + (-t312 - t320) * t347 + t376;
t273 = t320 * t348 + (-t292 - t307) * t361 + t374;
t272 = t292 * t347 + (-t294 - t308) * t348 + t378;
t271 = qJD(5) * t321 + t394 * t361 + (-t312 - t393) * t347 + t376;
t270 = qJD(5) * t323 + t393 * t348 + (-t307 - t395) * t361 + t374;
t269 = qJD(5) * t363 * t368 + t395 * t347 + (-t308 - t394) * t348 + t378;
t1 = m(3) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + ((t369 * t350 + t372 * t380) * qJD(1) + (t369 ^ 2 * t327 + (t382 * t372 + (-t326 + t381) * t369) * t372) * qJD(2)) * t391 / 0.2e1 - ((-t372 * t350 + t369 * t380) * qJD(1) + (t372 ^ 2 * t326 + (t381 * t369 + (-t327 + t382) * t372) * t369) * qJD(2)) * t390 / 0.2e1 + qJD(1) * ((t371 * t351 + t368 * t352) * qJD(1) + ((t330 * t371 + t333 * t368) * t369 - (t329 * t371 + t332 * t368) * t372) * qJD(2)) / 0.2e1 + m(4) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(5) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + (m(2) * (t354 ^ 2 + t355 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t414 * t323 + t413 * t324 + t328 * t343 + t331 * t344 + t410 * t398) * t361 + (t303 * t343 + t305 * t344 + t418 * t323 + t416 * t324 + t412 * t398) * t348 + (t343 * t304 + t344 * t306 + t417 * t323 + t415 * t324 + t411 * t398) * t347) * t347 / 0.2e1 + ((t414 * t321 + t413 * t322 + t328 * t341 + t331 * t342 + t410 * t399) * t361 + (t341 * t303 + t342 * t305 + t418 * t321 + t416 * t322 + t412 * t399) * t348 + (t304 * t341 + t306 * t342 + t417 * t321 + t415 * t322 + t411 * t399) * t347) * t348 / 0.2e1 + ((-t411 * t347 - t412 * t348 - t410 * t361) * t371 + ((-t328 * t367 + t331 * t370 + t414 * t363 + t413 * t364) * t361 + (-t303 * t367 + t305 * t370 + t418 * t363 + t416 * t364) * t348 + (-t304 * t367 + t306 * t370 + t417 * t363 + t415 * t364) * t347) * t368) * t361 / 0.2e1;
T = t1;
