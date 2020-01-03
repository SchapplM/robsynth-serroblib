% Calculate kinetic energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:54
% EndTime: 2019-12-31 19:34:56
% DurationCPUTime: 2.26s
% Computational Cost: add. (893->205), mult. (1235->317), div. (0->0), fcn. (1132->8), ass. (0->119)
t447 = Icges(4,4) + Icges(5,6);
t446 = Icges(4,1) + Icges(5,2);
t445 = -Icges(4,2) - Icges(5,3);
t358 = qJ(2) + pkin(8);
t356 = cos(t358);
t444 = t447 * t356;
t355 = sin(t358);
t443 = t447 * t355;
t442 = -Icges(5,4) + Icges(4,5);
t441 = Icges(5,5) - Icges(4,6);
t440 = t445 * t355 + t444;
t439 = -t446 * t356 + t443;
t362 = sin(qJ(1));
t365 = cos(qJ(1));
t438 = t440 * t362 + t441 * t365;
t437 = -t441 * t362 + t440 * t365;
t436 = t439 * t362 + t442 * t365;
t435 = t442 * t362 - t439 * t365;
t434 = t445 * t356 - t443;
t433 = t446 * t355 + t444;
t432 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t361 = sin(qJ(2));
t364 = cos(qJ(2));
t431 = Icges(3,5) * t364 - Icges(3,6) * t361 + t441 * t355 + t442 * t356;
t430 = t431 * t362 - t432 * t365;
t429 = t432 * t362 + t431 * t365;
t428 = Icges(3,5) * t361 + Icges(3,6) * t364 + t442 * t355 - t441 * t356;
t417 = Icges(3,4) * t361;
t345 = Icges(3,2) * t364 + t417;
t416 = Icges(3,4) * t364;
t346 = Icges(3,1) * t361 + t416;
t427 = -t345 * t361 + t346 * t364 + t434 * t355 + t433 * t356;
t386 = -Icges(3,2) * t361 + t416;
t319 = Icges(3,6) * t362 + t386 * t365;
t388 = Icges(3,1) * t364 - t417;
t321 = Icges(3,5) * t362 + t388 * t365;
t426 = -t319 * t361 + t321 * t364 - t437 * t355 + t435 * t356;
t318 = -Icges(3,6) * t365 + t386 * t362;
t320 = -Icges(3,5) * t365 + t388 * t362;
t425 = t318 * t361 - t320 * t364 + t438 * t355 + t436 * t356;
t421 = pkin(2) * t361;
t419 = t364 * pkin(2);
t411 = t356 * t362;
t410 = t356 * t365;
t360 = sin(qJ(5));
t409 = t362 * t360;
t363 = cos(qJ(5));
t408 = t362 * t363;
t407 = t365 * t360;
t406 = t365 * t363;
t295 = -qJ(3) * t365 + t419 * t362;
t296 = qJ(3) * t362 + t419 * t365;
t401 = qJD(2) * t365;
t402 = qJD(2) * t362;
t405 = t295 * t402 + t296 * t401;
t352 = t362 * pkin(1) - t365 * pkin(6);
t404 = -t295 - t352;
t357 = qJD(3) * t362;
t400 = qJD(4) * t355;
t403 = t365 * t400 + t357;
t399 = qJD(5) * t356;
t390 = pkin(3) * t356 + qJ(4) * t355;
t322 = t390 * t362;
t398 = -t322 + t404;
t395 = -t355 * pkin(3) + t356 * qJ(4) - t421;
t343 = qJD(1) * (t365 * pkin(1) + t362 * pkin(6));
t394 = qJD(1) * t296 - qJD(3) * t365 + t343;
t393 = rSges(3,1) * t364 - rSges(3,2) * t361;
t392 = rSges(4,1) * t356 - rSges(4,2) * t355;
t391 = -rSges(5,2) * t356 + rSges(5,3) * t355;
t389 = qJD(2) * (-t355 * rSges(4,1) - t356 * rSges(4,2) - t421);
t370 = qJD(2) * (t355 * rSges(5,2) + t356 * rSges(5,3) + t395);
t323 = t390 * t365;
t369 = qJD(1) * t323 + t362 * t400 + t394;
t368 = -qJD(4) * t356 + t322 * t402 + t323 * t401 + t405;
t367 = qJD(2) * (-pkin(7) * t355 + t395);
t353 = qJD(5) * t355 + qJD(1);
t351 = t365 * rSges(2,1) - t362 * rSges(2,2);
t350 = t362 * rSges(2,1) + t365 * rSges(2,2);
t349 = t361 * rSges(3,1) + t364 * rSges(3,2);
t342 = -t365 * pkin(4) + pkin(7) * t411;
t341 = t362 * pkin(4) + pkin(7) * t410;
t331 = t362 * t399 - t401;
t330 = t365 * t399 + t402;
t329 = t355 * t409 - t406;
t328 = t355 * t408 + t407;
t327 = t355 * t407 + t408;
t326 = t355 * t406 - t409;
t325 = t362 * rSges(3,3) + t393 * t365;
t324 = -t365 * rSges(3,3) + t393 * t362;
t314 = -t365 * rSges(5,1) + t391 * t362;
t313 = t362 * rSges(5,1) + t391 * t365;
t312 = t362 * rSges(4,3) + t392 * t365;
t311 = -t365 * rSges(4,3) + t392 * t362;
t294 = t355 * rSges(6,3) + (-rSges(6,1) * t360 - rSges(6,2) * t363) * t356;
t292 = Icges(6,5) * t355 + (-Icges(6,1) * t360 - Icges(6,4) * t363) * t356;
t291 = Icges(6,6) * t355 + (-Icges(6,4) * t360 - Icges(6,2) * t363) * t356;
t290 = Icges(6,3) * t355 + (-Icges(6,5) * t360 - Icges(6,6) * t363) * t356;
t287 = qJD(1) * t325 - t349 * t402 + t343;
t286 = -t349 * t401 + (-t324 - t352) * qJD(1);
t285 = (t324 * t362 + t325 * t365) * qJD(2);
t284 = t329 * rSges(6,1) + t328 * rSges(6,2) + rSges(6,3) * t411;
t283 = t327 * rSges(6,1) + t326 * rSges(6,2) + rSges(6,3) * t410;
t282 = Icges(6,1) * t329 + Icges(6,4) * t328 + Icges(6,5) * t411;
t281 = Icges(6,1) * t327 + Icges(6,4) * t326 + Icges(6,5) * t410;
t280 = Icges(6,4) * t329 + Icges(6,2) * t328 + Icges(6,6) * t411;
t279 = Icges(6,4) * t327 + Icges(6,2) * t326 + Icges(6,6) * t410;
t278 = Icges(6,5) * t329 + Icges(6,6) * t328 + Icges(6,3) * t411;
t277 = Icges(6,5) * t327 + Icges(6,6) * t326 + Icges(6,3) * t410;
t276 = qJD(1) * t312 + t362 * t389 + t394;
t275 = t357 + t365 * t389 + (-t311 + t404) * qJD(1);
t274 = (t311 * t362 + t312 * t365) * qJD(2) + t405;
t273 = qJD(1) * t313 + t362 * t370 + t369;
t272 = t365 * t370 + (-t314 + t398) * qJD(1) + t403;
t271 = (t313 * t365 + t314 * t362) * qJD(2) + t368;
t270 = qJD(1) * t341 + t353 * t283 - t330 * t294 + t362 * t367 + t369;
t269 = -t353 * t284 + t331 * t294 + t365 * t367 + (-t342 + t398) * qJD(1) + t403;
t268 = -t331 * t283 + t330 * t284 + (t341 * t365 + t342 * t362) * qJD(2) + t368;
t1 = m(3) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(4) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(5) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + m(6) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + t330 * ((t277 * t410 + t326 * t279 + t327 * t281) * t330 + (t278 * t410 + t326 * t280 + t327 * t282) * t331 + (t290 * t410 + t326 * t291 + t327 * t292) * t353) / 0.2e1 + t331 * ((t277 * t411 + t328 * t279 + t329 * t281) * t330 + (t278 * t411 + t328 * t280 + t329 * t282) * t331 + (t290 * t411 + t328 * t291 + t329 * t292) * t353) / 0.2e1 + t353 * ((t277 * t330 + t278 * t331 + t290 * t353) * t355 + ((-t279 * t363 - t281 * t360) * t330 + (-t280 * t363 - t282 * t360) * t331 + (-t291 * t363 - t292 * t360) * t353) * t356) / 0.2e1 + (m(2) * (t350 ^ 2 + t351 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t364 * t318 - t361 * t320 + t436 * t355 - t438 * t356) * t365 + (t364 * t319 + t361 * t321 + t435 * t355 + t437 * t356) * t362) * qJD(2) + (t364 * t345 + t361 * t346 + t433 * t355 - t434 * t356) * qJD(1)) * qJD(1) / 0.2e1 + ((t429 * t362 ^ 2 + (t425 * t365 + (t426 - t430) * t362) * t365) * qJD(2) + (t428 * t362 + t427 * t365) * qJD(1)) * t402 / 0.2e1 - ((t430 * t365 ^ 2 + (t426 * t362 + (t425 - t429) * t365) * t362) * qJD(2) + (t427 * t362 - t428 * t365) * qJD(1)) * t401 / 0.2e1;
T = t1;
