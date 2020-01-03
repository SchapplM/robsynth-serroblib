% Calculate kinetic energy for
% S5RRRRP7
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:06
% EndTime: 2019-12-31 21:56:07
% DurationCPUTime: 1.79s
% Computational Cost: add. (1118->205), mult. (1537->327), div. (0->0), fcn. (1500->8), ass. (0->120)
t431 = Icges(5,1) + Icges(6,1);
t430 = -Icges(5,4) + Icges(6,5);
t429 = Icges(6,4) + Icges(5,5);
t428 = Icges(5,2) + Icges(6,3);
t427 = -Icges(6,6) + Icges(5,6);
t426 = -Icges(5,3) - Icges(6,2);
t425 = rSges(6,1) + pkin(4);
t424 = rSges(6,3) + qJ(5);
t362 = qJ(2) + qJ(3);
t361 = cos(t362);
t366 = cos(qJ(4));
t368 = cos(qJ(1));
t398 = t366 * t368;
t363 = sin(qJ(4));
t365 = sin(qJ(1));
t401 = t363 * t365;
t336 = t361 * t401 + t398;
t399 = t365 * t366;
t400 = t363 * t368;
t337 = t361 * t399 - t400;
t360 = sin(t362);
t403 = t360 * t365;
t423 = t428 * t336 + t430 * t337 - t427 * t403;
t338 = t361 * t400 - t399;
t339 = t361 * t398 + t401;
t402 = t360 * t368;
t422 = t428 * t338 + t430 * t339 - t427 * t402;
t421 = -t427 * t336 + t429 * t337 - t426 * t403;
t420 = -t427 * t338 + t429 * t339 - t426 * t402;
t419 = t430 * t336 + t431 * t337 + t429 * t403;
t418 = t430 * t338 + t431 * t339 + t429 * t402;
t417 = t427 * t361 + (t428 * t363 + t430 * t366) * t360;
t416 = t426 * t361 + (-t427 * t363 + t429 * t366) * t360;
t415 = -t429 * t361 + (t430 * t363 + t431 * t366) * t360;
t367 = cos(qJ(2));
t409 = pkin(2) * t367;
t364 = sin(qJ(2));
t407 = Icges(3,4) * t364;
t406 = Icges(3,4) * t367;
t405 = Icges(4,4) * t360;
t404 = Icges(4,4) * t361;
t397 = rSges(6,2) * t403 + t424 * t336 + t425 * t337;
t396 = rSges(6,2) * t402 + t424 * t338 + t425 * t339;
t311 = -pkin(7) * t368 + t365 * t409;
t312 = pkin(7) * t365 + t368 * t409;
t359 = qJD(2) * t365;
t392 = qJD(2) * t368;
t395 = t311 * t359 + t312 * t392;
t394 = -rSges(6,2) * t361 + (t424 * t363 + t425 * t366) * t360;
t354 = pkin(1) * t365 - pkin(6) * t368;
t393 = -t311 - t354;
t346 = qJD(3) * t365 + t359;
t391 = qJD(4) * t360;
t390 = pkin(2) * qJD(2) * t364;
t347 = (-qJD(2) - qJD(3)) * t368;
t389 = t368 * t390;
t388 = pkin(3) * t361 + pkin(8) * t360;
t387 = rSges(3,1) * t367 - rSges(3,2) * t364;
t386 = rSges(4,1) * t361 - rSges(4,2) * t360;
t385 = Icges(3,1) * t367 - t407;
t384 = Icges(4,1) * t361 - t405;
t383 = -Icges(3,2) * t364 + t406;
t382 = -Icges(4,2) * t360 + t404;
t381 = Icges(3,5) * t367 - Icges(3,6) * t364;
t380 = Icges(4,5) * t361 - Icges(4,6) * t360;
t325 = -Icges(3,6) * t368 + t365 * t383;
t327 = -Icges(3,5) * t368 + t365 * t385;
t379 = t325 * t364 - t327 * t367;
t326 = Icges(3,6) * t365 + t368 * t383;
t328 = Icges(3,5) * t365 + t368 * t385;
t378 = -t326 * t364 + t328 * t367;
t349 = Icges(3,2) * t367 + t407;
t350 = Icges(3,1) * t364 + t406;
t377 = -t349 * t364 + t350 * t367;
t334 = t388 * t365;
t335 = t388 * t368;
t376 = t346 * t334 - t335 * t347 + t395;
t345 = qJD(1) * (pkin(1) * t368 + pkin(6) * t365);
t375 = qJD(1) * t312 - t365 * t390 + t345;
t374 = (Icges(4,5) * t360 + Icges(4,6) * t361) * qJD(1) + (-Icges(4,3) * t368 + t365 * t380) * t347 + (Icges(4,3) * t365 + t368 * t380) * t346;
t344 = pkin(3) * t360 - pkin(8) * t361;
t373 = qJD(1) * t335 - t344 * t346 + t375;
t372 = t347 * t344 + (-t334 + t393) * qJD(1) - t389;
t316 = -Icges(4,6) * t368 + t365 * t382;
t317 = Icges(4,6) * t365 + t368 * t382;
t318 = -Icges(4,5) * t368 + t365 * t384;
t319 = Icges(4,5) * t365 + t368 * t384;
t341 = Icges(4,2) * t361 + t405;
t342 = Icges(4,1) * t360 + t404;
t371 = (-t317 * t360 + t319 * t361) * t346 + (-t316 * t360 + t318 * t361) * t347 + (-t341 * t360 + t342 * t361) * qJD(1);
t355 = -qJD(4) * t361 + qJD(1);
t353 = rSges(2,1) * t368 - rSges(2,2) * t365;
t352 = rSges(2,1) * t365 + rSges(2,2) * t368;
t351 = rSges(3,1) * t364 + rSges(3,2) * t367;
t348 = Icges(3,5) * t364 + Icges(3,6) * t367;
t343 = rSges(4,1) * t360 + rSges(4,2) * t361;
t332 = rSges(3,3) * t365 + t368 * t387;
t331 = -rSges(3,3) * t368 + t365 * t387;
t330 = t365 * t391 + t347;
t329 = t368 * t391 + t346;
t324 = Icges(3,3) * t365 + t368 * t381;
t323 = -Icges(3,3) * t368 + t365 * t381;
t321 = rSges(4,3) * t365 + t368 * t386;
t320 = -rSges(4,3) * t368 + t365 * t386;
t310 = -rSges(5,3) * t361 + (rSges(5,1) * t366 - rSges(5,2) * t363) * t360;
t296 = qJD(1) * t332 - t351 * t359 + t345;
t295 = -t351 * t392 + (-t331 - t354) * qJD(1);
t294 = rSges(5,1) * t339 - rSges(5,2) * t338 + rSges(5,3) * t402;
t292 = rSges(5,1) * t337 - rSges(5,2) * t336 + rSges(5,3) * t403;
t278 = (t331 * t365 + t332 * t368) * qJD(2);
t277 = qJD(1) * t321 - t343 * t346 + t375;
t276 = -t389 + t343 * t347 + (-t320 + t393) * qJD(1);
t275 = t320 * t346 - t321 * t347 + t395;
t274 = t294 * t355 - t310 * t329 + t373;
t273 = -t292 * t355 + t310 * t330 + t372;
t272 = t292 * t329 - t294 * t330 + t376;
t271 = qJD(5) * t336 - t329 * t394 + t355 * t396 + t373;
t270 = qJD(5) * t338 + t330 * t394 - t355 * t397 + t372;
t269 = qJD(5) * t360 * t363 + t329 * t397 - t330 * t396 + t376;
t1 = m(3) * (t278 ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + ((t365 * t348 + t368 * t377) * qJD(1) + (t365 ^ 2 * t324 + (t379 * t368 + (-t323 + t378) * t365) * t368) * qJD(2)) * t359 / 0.2e1 - ((-t368 * t348 + t365 * t377) * qJD(1) + (t368 ^ 2 * t323 + (t378 * t365 + (-t324 + t379) * t368) * t365) * qJD(2)) * t392 / 0.2e1 + m(4) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t346 * (t374 * t365 + t371 * t368) / 0.2e1 + t347 * (t371 * t365 - t374 * t368) / 0.2e1 + m(5) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + ((t417 * t338 + t415 * t339 + t416 * t402) * t355 + (t423 * t338 + t419 * t339 + t421 * t402) * t330 + (t422 * t338 + t418 * t339 + t420 * t402) * t329) * t329 / 0.2e1 + ((t417 * t336 + t415 * t337 + t416 * t403) * t355 + (t423 * t336 + t419 * t337 + t421 * t403) * t330 + (t422 * t336 + t418 * t337 + t420 * t403) * t329) * t330 / 0.2e1 + ((-t420 * t329 - t421 * t330 - t416 * t355) * t361 + ((t417 * t363 + t415 * t366) * t355 + (t423 * t363 + t419 * t366) * t330 + (t422 * t363 + t418 * t366) * t329) * t360) * t355 / 0.2e1 + (m(2) * (t352 ^ 2 + t353 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t326 * t367 + t328 * t364) * t365 - (t325 * t367 + t327 * t364) * t368) * qJD(2) + (t317 * t361 + t319 * t360) * t346 + (t316 * t361 + t318 * t360) * t347 + (t361 * t341 + t360 * t342 + t367 * t349 + t364 * t350) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
