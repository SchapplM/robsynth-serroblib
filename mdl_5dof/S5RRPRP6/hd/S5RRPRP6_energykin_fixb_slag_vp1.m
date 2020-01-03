% Calculate kinetic energy for
% S5RRPRP6
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:44
% EndTime: 2019-12-31 19:56:46
% DurationCPUTime: 1.73s
% Computational Cost: add. (1106->207), mult. (1504->315), div. (0->0), fcn. (1457->8), ass. (0->117)
t451 = Icges(5,1) + Icges(6,1);
t450 = Icges(5,4) + Icges(6,4);
t449 = -Icges(6,5) - Icges(5,5);
t448 = Icges(5,2) + Icges(6,2);
t447 = -Icges(6,6) - Icges(5,6);
t446 = Icges(3,3) + Icges(4,3);
t445 = -Icges(6,3) - Icges(5,3);
t364 = qJ(2) + pkin(8);
t361 = sin(t364);
t362 = cos(t364);
t368 = sin(qJ(2));
t371 = cos(qJ(2));
t444 = Icges(3,5) * t371 + Icges(4,5) * t362 - Icges(3,6) * t368 - Icges(4,6) * t361;
t370 = cos(qJ(4));
t372 = cos(qJ(1));
t407 = t370 * t372;
t367 = sin(qJ(4));
t369 = sin(qJ(1));
t410 = t367 * t369;
t339 = -t362 * t410 - t407;
t408 = t369 * t370;
t409 = t367 * t372;
t340 = t362 * t408 - t409;
t412 = t361 * t369;
t443 = -t447 * t339 - t449 * t340 - t445 * t412;
t341 = -t362 * t409 + t408;
t342 = t362 * t407 + t410;
t411 = t361 * t372;
t442 = -t447 * t341 - t449 * t342 - t445 * t411;
t441 = t448 * t339 + t450 * t340 - t447 * t412;
t440 = t448 * t341 + t450 * t342 - t447 * t411;
t439 = t450 * t339 + t451 * t340 - t449 * t412;
t438 = t450 * t341 + t451 * t342 - t449 * t411;
t437 = t445 * t362 + (t447 * t367 - t449 * t370) * t361;
t436 = t447 * t362 + (-t448 * t367 + t450 * t370) * t361;
t435 = t449 * t362 + (-t450 * t367 + t451 * t370) * t361;
t434 = t444 * t369 - t446 * t372;
t433 = t446 * t369 + t444 * t372;
t432 = Icges(3,5) * t368 + Icges(4,5) * t361 + Icges(3,6) * t371 + Icges(4,6) * t362;
t414 = Icges(4,4) * t361;
t346 = Icges(4,2) * t362 + t414;
t413 = Icges(4,4) * t362;
t347 = Icges(4,1) * t361 + t413;
t416 = Icges(3,4) * t368;
t352 = Icges(3,2) * t371 + t416;
t415 = Icges(3,4) * t371;
t353 = Icges(3,1) * t368 + t415;
t431 = -t346 * t361 + t347 * t362 - t352 * t368 + t353 * t371;
t386 = -Icges(4,2) * t361 + t413;
t321 = Icges(4,6) * t369 + t372 * t386;
t388 = Icges(4,1) * t362 - t414;
t323 = Icges(4,5) * t369 + t372 * t388;
t387 = -Icges(3,2) * t368 + t415;
t332 = Icges(3,6) * t369 + t372 * t387;
t389 = Icges(3,1) * t371 - t416;
t334 = Icges(3,5) * t369 + t372 * t389;
t430 = -t321 * t361 + t323 * t362 - t332 * t368 + t334 * t371;
t320 = -Icges(4,6) * t372 + t369 * t386;
t322 = -Icges(4,5) * t372 + t369 * t388;
t331 = -Icges(3,6) * t372 + t369 * t387;
t333 = -Icges(3,5) * t372 + t369 * t389;
t429 = t320 * t361 - t322 * t362 + t331 * t368 - t333 * t371;
t422 = pkin(2) * t368;
t420 = pkin(2) * t371;
t419 = pkin(4) * t370;
t375 = qJ(5) * t361 + t362 * t419;
t406 = rSges(6,1) * t340 + rSges(6,2) * t339 + rSges(6,3) * t412 - pkin(4) * t409 + t369 * t375;
t405 = rSges(6,1) * t342 + rSges(6,2) * t341 + rSges(6,3) * t411 + pkin(4) * t410 + t372 * t375;
t404 = (-qJ(5) - rSges(6,3)) * t362 + (rSges(6,1) * t370 - rSges(6,2) * t367 + t419) * t361;
t316 = -qJ(3) * t372 + t369 * t420;
t317 = qJ(3) * t369 + t372 * t420;
t400 = qJD(2) * t372;
t401 = qJD(2) * t369;
t403 = t316 * t401 + t317 * t400;
t357 = t369 * pkin(1) - t372 * pkin(6);
t402 = -t316 - t357;
t399 = qJD(4) * t361;
t395 = pkin(3) * t362 + pkin(7) * t361;
t337 = t395 * t369;
t338 = t395 * t372;
t396 = t337 * t401 + t338 * t400 + t403;
t350 = qJD(1) * (t372 * pkin(1) + t369 * pkin(6));
t394 = qJD(1) * t317 - qJD(3) * t372 + t350;
t393 = rSges(3,1) * t371 - rSges(3,2) * t368;
t392 = rSges(4,1) * t362 - rSges(4,2) * t361;
t391 = qJD(2) * (-rSges(4,1) * t361 - rSges(4,2) * t362 - t422);
t390 = (-t361 * pkin(3) + t362 * pkin(7) - t422) * qJD(2);
t377 = qJD(1) * t338 + t394;
t363 = qJD(3) * t369;
t376 = t363 + (-t337 + t402) * qJD(1);
t374 = qJD(5) * t361 + t390;
t358 = -qJD(4) * t362 + qJD(1);
t356 = rSges(2,1) * t372 - rSges(2,2) * t369;
t355 = rSges(2,1) * t369 + rSges(2,2) * t372;
t354 = rSges(3,1) * t368 + rSges(3,2) * t371;
t344 = t369 * t399 - t400;
t343 = t372 * t399 + t401;
t336 = t369 * rSges(3,3) + t372 * t393;
t335 = -t372 * rSges(3,3) + t369 * t393;
t325 = t369 * rSges(4,3) + t372 * t392;
t324 = -t372 * rSges(4,3) + t369 * t392;
t315 = -t362 * rSges(5,3) + (rSges(5,1) * t370 - rSges(5,2) * t367) * t361;
t303 = qJD(1) * t336 - t354 * t401 + t350;
t302 = -t354 * t400 + (-t335 - t357) * qJD(1);
t301 = (t335 * t369 + t336 * t372) * qJD(2);
t300 = rSges(5,1) * t342 + rSges(5,2) * t341 + rSges(5,3) * t411;
t298 = rSges(5,1) * t340 + rSges(5,2) * t339 + rSges(5,3) * t412;
t282 = qJD(1) * t325 + t369 * t391 + t394;
t281 = t363 + t372 * t391 + (-t324 + t402) * qJD(1);
t280 = (t324 * t369 + t325 * t372) * qJD(2) + t403;
t279 = t358 * t300 - t343 * t315 + t369 * t390 + t377;
t278 = -t358 * t298 + t344 * t315 + t372 * t390 + t376;
t277 = t298 * t343 - t300 * t344 + t396;
t276 = -t343 * t404 + t358 * t405 + t369 * t374 + t377;
t275 = t344 * t404 - t358 * t406 + t372 * t374 + t376;
t274 = -qJD(5) * t362 + t343 * t406 - t344 * t405 + t396;
t1 = m(3) * (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(4) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(5) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(6) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t436 * t341 + t435 * t342 + t437 * t411) * t358 + (t441 * t341 + t439 * t342 + t443 * t411) * t344 + (t440 * t341 + t438 * t342 + t442 * t411) * t343) * t343 / 0.2e1 + ((t436 * t339 + t435 * t340 + t437 * t412) * t358 + (t441 * t339 + t439 * t340 + t443 * t412) * t344 + (t440 * t339 + t438 * t340 + t442 * t412) * t343) * t344 / 0.2e1 + ((-t442 * t343 - t443 * t344 - t437 * t358) * t362 + ((-t436 * t367 + t435 * t370) * t358 + (-t441 * t367 + t439 * t370) * t344 + (-t440 * t367 + t438 * t370) * t343) * t361) * t358 / 0.2e1 + (m(2) * (t355 ^ 2 + t356 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t320 * t362 - t322 * t361 - t331 * t371 - t333 * t368) * t372 + (t321 * t362 + t323 * t361 + t332 * t371 + t334 * t368) * t369) * qJD(2) + (t346 * t362 + t347 * t361 + t371 * t352 + t368 * t353) * qJD(1)) * qJD(1) / 0.2e1 + ((t433 * t369 ^ 2 + (t429 * t372 + (t430 - t434) * t369) * t372) * qJD(2) + (t432 * t369 + t431 * t372) * qJD(1)) * t401 / 0.2e1 - ((t434 * t372 ^ 2 + (t430 * t369 + (t429 - t433) * t372) * t369) * qJD(2) + (t431 * t369 - t432 * t372) * qJD(1)) * t400 / 0.2e1;
T = t1;
