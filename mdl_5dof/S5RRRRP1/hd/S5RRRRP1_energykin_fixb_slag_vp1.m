% Calculate kinetic energy for
% S5RRRRP1
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:43
% EndTime: 2019-12-05 18:44:45
% DurationCPUTime: 1.93s
% Computational Cost: add. (1060->181), mult. (1117->284), div. (0->0), fcn. (958->8), ass. (0->117)
t447 = Icges(5,4) + Icges(6,4);
t446 = Icges(5,1) + Icges(6,1);
t445 = Icges(5,2) + Icges(6,2);
t361 = qJ(2) + qJ(3);
t358 = qJ(4) + t361;
t351 = cos(t358);
t444 = t447 * t351;
t350 = sin(t358);
t443 = t447 * t350;
t442 = Icges(5,5) + Icges(6,5);
t441 = Icges(5,6) + Icges(6,6);
t440 = -t445 * t350 + t444;
t439 = t446 * t351 - t443;
t438 = rSges(6,1) + pkin(4);
t437 = Icges(5,3) + Icges(6,3);
t363 = sin(qJ(1));
t365 = cos(qJ(1));
t436 = t440 * t363 - t441 * t365;
t435 = t441 * t363 + t440 * t365;
t434 = t439 * t363 - t442 * t365;
t433 = t442 * t363 + t439 * t365;
t432 = t445 * t351 + t443;
t431 = t446 * t350 + t444;
t430 = -t441 * t350 + t442 * t351;
t429 = rSges(6,3) + qJ(5);
t428 = -rSges(6,2) * t350 + t438 * t351;
t354 = qJD(2) * t363;
t340 = qJD(3) * t363 + t354;
t332 = qJD(4) * t363 + t340;
t401 = -qJD(2) - qJD(3);
t333 = (-qJD(4) + t401) * t365;
t427 = (-t436 * t350 + t434 * t351) * t333 + (-t435 * t350 + t433 * t351) * t332 + (-t432 * t350 + t431 * t351) * qJD(1);
t426 = (t430 * t363 - t437 * t365) * t333 + (t437 * t363 + t430 * t365) * t332 + (t442 * t350 + t441 * t351) * qJD(1);
t355 = sin(t361);
t422 = pkin(3) * t355;
t364 = cos(qJ(2));
t420 = t364 * pkin(2);
t362 = sin(qJ(2));
t418 = Icges(3,4) * t362;
t417 = Icges(3,4) * t364;
t416 = Icges(4,4) * t355;
t356 = cos(t361);
t415 = Icges(4,4) * t356;
t410 = t428 * t363 - t429 * t365;
t409 = t429 * t363 + t428 * t365;
t304 = -pkin(7) * t365 + t420 * t363;
t305 = pkin(7) * t363 + t420 * t365;
t402 = qJD(2) * t365;
t408 = t304 * t354 + t305 * t402;
t348 = pkin(1) * t363 - pkin(6) * t365;
t407 = -t304 - t348;
t405 = pkin(3) * t356;
t400 = pkin(2) * qJD(2) * t362;
t283 = -pkin(8) * t365 + t405 * t363;
t399 = -t283 + t407;
t398 = rSges(6,2) * t351 + t438 * t350;
t397 = t365 * t400;
t396 = rSges(3,1) * t364 - rSges(3,2) * t362;
t395 = rSges(4,1) * t356 - rSges(4,2) * t355;
t394 = rSges(5,1) * t351 - rSges(5,2) * t350;
t392 = Icges(3,1) * t364 - t418;
t391 = Icges(4,1) * t356 - t416;
t388 = -Icges(3,2) * t362 + t417;
t387 = -Icges(4,2) * t355 + t415;
t384 = Icges(3,5) * t364 - Icges(3,6) * t362;
t383 = Icges(4,5) * t356 - Icges(4,6) * t355;
t316 = -Icges(3,6) * t365 + t388 * t363;
t318 = -Icges(3,5) * t365 + t392 * t363;
t380 = t316 * t362 - t318 * t364;
t317 = Icges(3,6) * t363 + t388 * t365;
t319 = Icges(3,5) * t363 + t392 * t365;
t379 = -t317 * t362 + t319 * t364;
t343 = Icges(3,2) * t364 + t418;
t344 = Icges(3,1) * t362 + t417;
t378 = -t343 * t362 + t344 * t364;
t341 = t401 * t365;
t377 = t341 * t422 - t397;
t284 = pkin(8) * t363 + t405 * t365;
t376 = t340 * t283 - t284 * t341 + t408;
t339 = qJD(1) * (pkin(1) * t365 + pkin(6) * t363);
t375 = qJD(1) * t305 - t363 * t400 + t339;
t372 = (Icges(4,5) * t355 + Icges(4,6) * t356) * qJD(1) + (-Icges(4,3) * t365 + t383 * t363) * t341 + (Icges(4,3) * t363 + t383 * t365) * t340;
t371 = qJD(1) * t284 - t340 * t422 + t375;
t308 = -Icges(4,6) * t365 + t387 * t363;
t309 = Icges(4,6) * t363 + t387 * t365;
t310 = -Icges(4,5) * t365 + t391 * t363;
t311 = Icges(4,5) * t363 + t391 * t365;
t335 = Icges(4,2) * t356 + t416;
t336 = Icges(4,1) * t355 + t415;
t368 = (-t309 * t355 + t311 * t356) * t340 + (-t308 * t355 + t310 * t356) * t341 + (-t335 * t355 + t336 * t356) * qJD(1);
t347 = rSges(2,1) * t365 - rSges(2,2) * t363;
t346 = rSges(2,1) * t363 + rSges(2,2) * t365;
t345 = rSges(3,1) * t362 + rSges(3,2) * t364;
t342 = Icges(3,5) * t362 + Icges(3,6) * t364;
t337 = rSges(4,1) * t355 + rSges(4,2) * t356;
t331 = rSges(5,1) * t350 + rSges(5,2) * t351;
t322 = rSges(3,3) * t363 + t396 * t365;
t321 = -rSges(3,3) * t365 + t396 * t363;
t315 = Icges(3,3) * t363 + t384 * t365;
t314 = -Icges(3,3) * t365 + t384 * t363;
t313 = rSges(4,3) * t363 + t395 * t365;
t312 = -rSges(4,3) * t365 + t395 * t363;
t303 = rSges(5,3) * t363 + t394 * t365;
t301 = -rSges(5,3) * t365 + t394 * t363;
t281 = qJD(1) * t322 - t345 * t354 + t339;
t280 = -t345 * t402 + (-t321 - t348) * qJD(1);
t278 = (t321 * t363 + t322 * t365) * qJD(2);
t275 = qJD(1) * t313 - t337 * t340 + t375;
t274 = -t397 + t337 * t341 + (-t312 + t407) * qJD(1);
t273 = t312 * t340 - t313 * t341 + t408;
t272 = qJD(1) * t303 - t331 * t332 + t371;
t271 = t331 * t333 + (-t301 + t399) * qJD(1) + t377;
t270 = t301 * t332 - t303 * t333 + t376;
t269 = t409 * qJD(1) - qJD(5) * t365 - t398 * t332 + t371;
t268 = qJD(5) * t363 + t398 * t333 + (t399 - t410) * qJD(1) + t377;
t267 = t410 * t332 - t409 * t333 + t376;
t1 = m(3) * (t278 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + ((t363 * t342 + t378 * t365) * qJD(1) + (t363 ^ 2 * t315 + (t380 * t365 + (-t314 + t379) * t363) * t365) * qJD(2)) * t354 / 0.2e1 - ((-t365 * t342 + t378 * t363) * qJD(1) + (t365 ^ 2 * t314 + (t379 * t363 + (-t315 + t380) * t365) * t363) * qJD(2)) * t402 / 0.2e1 + m(4) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t340 * (t372 * t363 + t368 * t365) / 0.2e1 + t341 * (t368 * t363 - t372 * t365) / 0.2e1 + m(5) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + m(6) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + (t426 * t363 + t427 * t365) * t332 / 0.2e1 + (t427 * t363 - t426 * t365) * t333 / 0.2e1 + (m(2) * (t346 ^ 2 + t347 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t317 * t364 + t319 * t362) * t363 - (t316 * t364 + t318 * t362) * t365) * qJD(2) + (t309 * t356 + t311 * t355) * t340 + (t308 * t356 + t310 * t355) * t341 + (t434 * t350 + t436 * t351) * t333 + (t433 * t350 + t435 * t351) * t332 + (t356 * t335 + t355 * t336 + t364 * t343 + t362 * t344 + t431 * t350 + t432 * t351) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
