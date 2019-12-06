% Calculate kinetic energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:24
% EndTime: 2019-12-05 16:02:26
% DurationCPUTime: 2.22s
% Computational Cost: add. (1287->249), mult. (3317->395), div. (0->0), fcn. (3920->10), ass. (0->119)
t450 = Icges(3,1) + Icges(4,2);
t449 = Icges(3,4) + Icges(4,6);
t448 = Icges(3,5) - Icges(4,4);
t447 = Icges(3,2) + Icges(4,3);
t446 = Icges(3,6) - Icges(4,5);
t445 = Icges(3,3) + Icges(4,1);
t405 = cos(pkin(5));
t403 = sin(pkin(5));
t404 = cos(pkin(9));
t430 = t403 * t404;
t402 = sin(pkin(9));
t431 = t402 * t403;
t408 = sin(qJ(2));
t410 = cos(qJ(2));
t437 = t445 * t405 + (t448 * t408 + t446 * t410) * t403;
t426 = t405 * t410;
t388 = t402 * t408 - t404 * t426;
t427 = t405 * t408;
t389 = t402 * t410 + t404 * t427;
t438 = -t446 * t388 + t448 * t389 - t445 * t430;
t390 = t402 * t426 + t404 * t408;
t391 = -t402 * t427 + t404 * t410;
t439 = -t446 * t390 + t448 * t391 + t445 * t431;
t444 = -t437 * t405 + t438 * t430 - t439 * t431;
t443 = t447 * t390 - t449 * t391 - t446 * t431;
t442 = t447 * t388 - t449 * t389 + t446 * t430;
t441 = -t449 * t390 + t450 * t391 + t448 * t431;
t440 = t449 * t388 - t450 * t389 + t448 * t430;
t436 = t446 * t405 + (t449 * t408 + t447 * t410) * t403;
t435 = t448 * t405 + (t450 * t408 + t449 * t410) * t403;
t433 = qJD(2) ^ 2;
t432 = cos(qJ(4));
t407 = sin(qJ(4));
t429 = t403 * t407;
t428 = t403 * t408;
t359 = t391 * pkin(2) + t390 * qJ(3);
t401 = qJD(2) * t405;
t425 = qJD(3) * t388 + t359 * t401;
t424 = qJD(2) * t403;
t399 = t402 * t424;
t371 = qJD(4) * t391 + t399;
t395 = qJD(4) * t428 + t401;
t423 = qJD(3) * t410;
t421 = t403 * t432;
t420 = t404 * t424;
t358 = t389 * pkin(2) + t388 * qJ(3);
t419 = t358 * t399 + t359 * t420 + qJD(1);
t394 = (pkin(2) * t408 - qJ(3) * t410) * t403;
t417 = (-t405 * rSges(4,1) - (-rSges(4,2) * t408 - rSges(4,3) * t410) * t403 - t394) * t403;
t416 = (-t405 * pkin(3) - pkin(7) * t428 - t394) * t403;
t372 = qJD(4) * t389 - t420;
t373 = pkin(3) * t431 + t391 * pkin(7);
t374 = -pkin(3) * t430 + t389 * pkin(7);
t413 = t373 * t420 + t374 * t399 - t403 * t423 + t419;
t412 = qJD(2) * t402 * t416 + t373 * t401 + t425;
t387 = qJD(3) * t390;
t411 = t387 + ((-t358 - t374) * t405 + t404 * t416) * qJD(2);
t409 = cos(qJ(5));
t406 = sin(qJ(5));
t393 = t405 * t432 - t410 * t429;
t392 = t405 * t407 + t410 * t421;
t381 = t405 * rSges(3,3) + (rSges(3,1) * t408 + rSges(3,2) * t410) * t403;
t370 = t393 * t409 + t406 * t428;
t369 = -t393 * t406 + t409 * t428;
t368 = t388 * t407 - t404 * t421;
t367 = t388 * t432 + t404 * t429;
t366 = t390 * t407 + t402 * t421;
t365 = -t390 * t432 + t402 * t429;
t363 = qJD(5) * t392 + t395;
t360 = t393 * pkin(4) + t392 * pkin(8);
t356 = t393 * rSges(5,1) - t392 * rSges(5,2) + rSges(5,3) * t428;
t355 = Icges(5,1) * t393 - Icges(5,4) * t392 + Icges(5,5) * t428;
t354 = Icges(5,4) * t393 - Icges(5,2) * t392 + Icges(5,6) * t428;
t353 = Icges(5,5) * t393 - Icges(5,6) * t392 + Icges(5,3) * t428;
t352 = t391 * rSges(3,1) - t390 * rSges(3,2) + rSges(3,3) * t431;
t351 = t389 * rSges(3,1) - t388 * rSges(3,2) - rSges(3,3) * t430;
t350 = -rSges(4,1) * t430 - t389 * rSges(4,2) + t388 * rSges(4,3);
t349 = rSges(4,1) * t431 - t391 * rSges(4,2) + t390 * rSges(4,3);
t334 = t368 * t409 + t389 * t406;
t333 = -t368 * t406 + t389 * t409;
t332 = t366 * t409 + t391 * t406;
t331 = -t366 * t406 + t391 * t409;
t330 = -qJD(5) * t367 + t372;
t329 = qJD(5) * t365 + t371;
t328 = t368 * pkin(4) - t367 * pkin(8);
t327 = t366 * pkin(4) + t365 * pkin(8);
t326 = (-t351 * t405 - t381 * t430) * qJD(2);
t325 = (t352 * t405 - t381 * t431) * qJD(2);
t324 = t370 * rSges(6,1) + t369 * rSges(6,2) + t392 * rSges(6,3);
t323 = Icges(6,1) * t370 + Icges(6,4) * t369 + Icges(6,5) * t392;
t322 = Icges(6,4) * t370 + Icges(6,2) * t369 + Icges(6,6) * t392;
t321 = Icges(6,5) * t370 + Icges(6,6) * t369 + Icges(6,3) * t392;
t320 = t368 * rSges(5,1) + t367 * rSges(5,2) + t389 * rSges(5,3);
t319 = t366 * rSges(5,1) - t365 * rSges(5,2) + t391 * rSges(5,3);
t318 = Icges(5,1) * t368 + Icges(5,4) * t367 + Icges(5,5) * t389;
t317 = Icges(5,1) * t366 - Icges(5,4) * t365 + Icges(5,5) * t391;
t316 = Icges(5,4) * t368 + Icges(5,2) * t367 + Icges(5,6) * t389;
t315 = Icges(5,4) * t366 - Icges(5,2) * t365 + Icges(5,6) * t391;
t314 = Icges(5,5) * t368 + Icges(5,6) * t367 + Icges(5,3) * t389;
t313 = Icges(5,5) * t366 - Icges(5,6) * t365 + Icges(5,3) * t391;
t312 = qJD(1) + (t351 * t402 + t352 * t404) * t424;
t311 = t334 * rSges(6,1) + t333 * rSges(6,2) - t367 * rSges(6,3);
t310 = t332 * rSges(6,1) + t331 * rSges(6,2) + t365 * rSges(6,3);
t309 = Icges(6,1) * t334 + Icges(6,4) * t333 - Icges(6,5) * t367;
t308 = Icges(6,1) * t332 + Icges(6,4) * t331 + Icges(6,5) * t365;
t307 = Icges(6,4) * t334 + Icges(6,2) * t333 - Icges(6,6) * t367;
t306 = Icges(6,4) * t332 + Icges(6,2) * t331 + Icges(6,6) * t365;
t305 = Icges(6,5) * t334 + Icges(6,6) * t333 - Icges(6,3) * t367;
t304 = Icges(6,5) * t332 + Icges(6,6) * t331 + Icges(6,3) * t365;
t303 = t387 + ((-t350 - t358) * t405 + t404 * t417) * qJD(2);
t302 = (t349 * t405 + t402 * t417) * qJD(2) + t425;
t301 = (-t423 + (t349 * t404 + t350 * t402) * qJD(2)) * t403 + t419;
t300 = -t395 * t320 + t372 * t356 + t411;
t299 = t395 * t319 - t371 * t356 + t412;
t298 = -t372 * t319 + t371 * t320 + t413;
t297 = -t363 * t311 + t330 * t324 - t395 * t328 + t372 * t360 + t411;
t296 = t363 * t310 - t329 * t324 + t395 * t327 - t371 * t360 + t412;
t295 = -t330 * t310 + t329 * t311 - t372 * t327 + t371 * t328 + t413;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t312 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(5) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + t371 * ((t391 * t313 - t365 * t315 + t366 * t317) * t371 + (t391 * t314 - t365 * t316 + t366 * t318) * t372 + (t391 * t353 - t365 * t354 + t366 * t355) * t395) / 0.2e1 + t372 * ((t389 * t313 + t367 * t315 + t368 * t317) * t371 + (t389 * t314 + t367 * t316 + t368 * t318) * t372 + (t389 * t353 + t367 * t354 + t368 * t355) * t395) / 0.2e1 + t395 * ((t313 * t428 - t392 * t315 + t393 * t317) * t371 + (t314 * t428 - t392 * t316 + t393 * t318) * t372 + (t353 * t428 - t392 * t354 + t393 * t355) * t395) / 0.2e1 + m(6) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + t329 * ((t365 * t304 + t331 * t306 + t332 * t308) * t329 + (t365 * t305 + t331 * t307 + t332 * t309) * t330 + (t365 * t321 + t331 * t322 + t332 * t323) * t363) / 0.2e1 + t330 * ((-t367 * t304 + t333 * t306 + t334 * t308) * t329 + (-t367 * t305 + t333 * t307 + t334 * t309) * t330 + (-t367 * t321 + t333 * t322 + t334 * t323) * t363) / 0.2e1 + t363 * ((t392 * t304 + t369 * t306 + t370 * t308) * t329 + (t392 * t305 + t369 * t307 + t370 * t309) * t330 + (t392 * t321 + t369 * t322 + t370 * t323) * t363) / 0.2e1 - ((t388 * t443 + t441 * t389) * t431 + (-t388 * t436 + t389 * t435) * t405 + (-t388 * t442 + t389 * t440 + t444) * t430) * t433 * t430 / 0.2e1 + ((t437 * t405 ^ 2 + (((t408 * t440 + t410 * t442) * t404 + (t441 * t408 - t443 * t410) * t402) * t403 + (t402 * t439 - t404 * t438 + t408 * t435 + t410 * t436) * t405) * t403) * t405 + ((-t390 * t442 + t391 * t440) * t430 + (-t390 * t436 + t391 * t435) * t405 + (t390 * t443 + t441 * t391 - t444) * t431) * t431) * t433 / 0.2e1;
T = t1;
