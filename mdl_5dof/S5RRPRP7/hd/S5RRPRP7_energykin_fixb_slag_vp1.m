% Calculate kinetic energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:53
% EndTime: 2019-12-31 19:59:55
% DurationCPUTime: 1.71s
% Computational Cost: add. (1076->201), mult. (1495->308), div. (0->0), fcn. (1458->8), ass. (0->116)
t447 = Icges(5,1) + Icges(6,1);
t446 = -Icges(5,4) + Icges(6,5);
t445 = Icges(6,4) + Icges(5,5);
t444 = Icges(5,2) + Icges(6,3);
t443 = -Icges(6,6) + Icges(5,6);
t442 = Icges(3,3) + Icges(4,3);
t441 = -Icges(5,3) - Icges(6,2);
t363 = qJ(2) + pkin(8);
t360 = sin(t363);
t361 = cos(t363);
t366 = sin(qJ(2));
t369 = cos(qJ(2));
t440 = Icges(3,5) * t369 + Icges(4,5) * t361 - Icges(3,6) * t366 - Icges(4,6) * t360;
t439 = rSges(6,1) + pkin(4);
t438 = rSges(6,3) + qJ(5);
t368 = cos(qJ(4));
t370 = cos(qJ(1));
t403 = t368 * t370;
t365 = sin(qJ(4));
t367 = sin(qJ(1));
t406 = t365 * t367;
t338 = t361 * t406 + t403;
t404 = t367 * t368;
t405 = t365 * t370;
t339 = t361 * t404 - t405;
t408 = t360 * t367;
t437 = t444 * t338 + t446 * t339 - t443 * t408;
t340 = t361 * t405 - t404;
t341 = t361 * t403 + t406;
t407 = t360 * t370;
t436 = t444 * t340 + t446 * t341 - t443 * t407;
t435 = -t443 * t338 + t445 * t339 - t441 * t408;
t434 = -t443 * t340 + t445 * t341 - t441 * t407;
t433 = t446 * t338 + t447 * t339 + t445 * t408;
t432 = t446 * t340 + t447 * t341 + t445 * t407;
t431 = t443 * t361 + (t444 * t365 + t446 * t368) * t360;
t430 = t441 * t361 + (-t443 * t365 + t445 * t368) * t360;
t429 = -t445 * t361 + (t446 * t365 + t447 * t368) * t360;
t428 = t440 * t367 - t442 * t370;
t427 = t442 * t367 + t440 * t370;
t426 = Icges(3,5) * t366 + Icges(4,5) * t360 + Icges(3,6) * t369 + Icges(4,6) * t361;
t410 = Icges(4,4) * t360;
t345 = Icges(4,2) * t361 + t410;
t409 = Icges(4,4) * t361;
t346 = Icges(4,1) * t360 + t409;
t412 = Icges(3,4) * t366;
t351 = Icges(3,2) * t369 + t412;
t411 = Icges(3,4) * t369;
t352 = Icges(3,1) * t366 + t411;
t425 = -t345 * t360 + t346 * t361 - t351 * t366 + t352 * t369;
t382 = -Icges(4,2) * t360 + t409;
t319 = Icges(4,6) * t367 + t370 * t382;
t384 = Icges(4,1) * t361 - t410;
t321 = Icges(4,5) * t367 + t370 * t384;
t383 = -Icges(3,2) * t366 + t411;
t330 = Icges(3,6) * t367 + t370 * t383;
t385 = Icges(3,1) * t369 - t412;
t332 = Icges(3,5) * t367 + t370 * t385;
t424 = -t319 * t360 + t321 * t361 - t330 * t366 + t332 * t369;
t318 = -Icges(4,6) * t370 + t367 * t382;
t320 = -Icges(4,5) * t370 + t367 * t384;
t329 = -Icges(3,6) * t370 + t367 * t383;
t331 = -Icges(3,5) * t370 + t367 * t385;
t423 = t318 * t360 - t320 * t361 + t329 * t366 - t331 * t369;
t416 = pkin(2) * t366;
t414 = pkin(2) * t369;
t402 = rSges(6,2) * t408 + t438 * t338 + t439 * t339;
t401 = rSges(6,2) * t407 + t438 * t340 + t439 * t341;
t314 = -qJ(3) * t370 + t367 * t414;
t315 = qJ(3) * t367 + t370 * t414;
t396 = qJD(2) * t370;
t397 = qJD(2) * t367;
t400 = t314 * t397 + t315 * t396;
t399 = -rSges(6,2) * t361 + (t438 * t365 + t439 * t368) * t360;
t356 = pkin(1) * t367 - pkin(6) * t370;
t398 = -t314 - t356;
t395 = qJD(4) * t360;
t391 = pkin(3) * t361 + pkin(7) * t360;
t336 = t391 * t367;
t337 = t391 * t370;
t392 = t336 * t397 + t337 * t396 + t400;
t349 = qJD(1) * (pkin(1) * t370 + pkin(6) * t367);
t390 = qJD(1) * t315 - qJD(3) * t370 + t349;
t389 = rSges(3,1) * t369 - rSges(3,2) * t366;
t388 = rSges(4,1) * t361 - rSges(4,2) * t360;
t387 = qJD(2) * (-rSges(4,1) * t360 - rSges(4,2) * t361 - t416);
t386 = qJD(2) * (-pkin(3) * t360 + pkin(7) * t361 - t416);
t373 = qJD(1) * t337 + t367 * t386 + t390;
t362 = qJD(3) * t367;
t372 = t362 + (-t336 + t398) * qJD(1) + t370 * t386;
t357 = -qJD(4) * t361 + qJD(1);
t355 = rSges(2,1) * t370 - rSges(2,2) * t367;
t354 = rSges(2,1) * t367 + rSges(2,2) * t370;
t353 = rSges(3,1) * t366 + rSges(3,2) * t369;
t343 = t367 * t395 - t396;
t342 = t370 * t395 + t397;
t335 = rSges(3,3) * t367 + t370 * t389;
t334 = -rSges(3,3) * t370 + t367 * t389;
t323 = rSges(4,3) * t367 + t370 * t388;
t322 = -rSges(4,3) * t370 + t367 * t388;
t313 = -rSges(5,3) * t361 + (rSges(5,1) * t368 - rSges(5,2) * t365) * t360;
t300 = qJD(1) * t335 - t353 * t397 + t349;
t299 = -t353 * t396 + (-t334 - t356) * qJD(1);
t298 = (t334 * t367 + t335 * t370) * qJD(2);
t297 = rSges(5,1) * t341 - rSges(5,2) * t340 + rSges(5,3) * t407;
t295 = rSges(5,1) * t339 - rSges(5,2) * t338 + rSges(5,3) * t408;
t281 = qJD(1) * t323 + t367 * t387 + t390;
t280 = t362 + t370 * t387 + (-t322 + t398) * qJD(1);
t279 = (t322 * t367 + t323 * t370) * qJD(2) + t400;
t278 = t297 * t357 - t313 * t342 + t373;
t277 = -t295 * t357 + t313 * t343 + t372;
t276 = t295 * t342 - t297 * t343 + t392;
t275 = qJD(5) * t338 - t342 * t399 + t357 * t401 + t373;
t274 = qJD(5) * t340 + t343 * t399 - t357 * t402 + t372;
t273 = qJD(5) * t360 * t365 + t342 * t402 - t343 * t401 + t392;
t1 = m(3) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + m(4) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(5) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + m(6) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + ((t431 * t340 + t429 * t341 + t430 * t407) * t357 + (t437 * t340 + t433 * t341 + t435 * t407) * t343 + (t436 * t340 + t432 * t341 + t434 * t407) * t342) * t342 / 0.2e1 + ((t431 * t338 + t429 * t339 + t430 * t408) * t357 + (t437 * t338 + t433 * t339 + t435 * t408) * t343 + (t436 * t338 + t432 * t339 + t434 * t408) * t342) * t343 / 0.2e1 + ((-t434 * t342 - t435 * t343 - t430 * t357) * t361 + ((t431 * t365 + t429 * t368) * t357 + (t437 * t365 + t433 * t368) * t343 + (t436 * t365 + t432 * t368) * t342) * t360) * t357 / 0.2e1 + (m(2) * (t354 ^ 2 + t355 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t318 * t361 - t360 * t320 - t329 * t369 - t331 * t366) * t370 + (t319 * t361 + t321 * t360 + t330 * t369 + t332 * t366) * t367) * qJD(2) + (t361 * t345 + t360 * t346 + t369 * t351 + t366 * t352) * qJD(1)) * qJD(1) / 0.2e1 + ((t427 * t367 ^ 2 + (t423 * t370 + (t424 - t428) * t367) * t370) * qJD(2) + (t426 * t367 + t425 * t370) * qJD(1)) * t397 / 0.2e1 - ((t428 * t370 ^ 2 + (t424 * t367 + (t423 - t427) * t370) * t367) * qJD(2) + (t425 * t367 - t426 * t370) * qJD(1)) * t396 / 0.2e1;
T = t1;
