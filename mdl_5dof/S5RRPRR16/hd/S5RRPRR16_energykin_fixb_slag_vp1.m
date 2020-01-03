% Calculate kinetic energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR16_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR16_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR16_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:41
% EndTime: 2019-12-31 20:44:43
% DurationCPUTime: 2.22s
% Computational Cost: add. (1339->258), mult. (3361->401), div. (0->0), fcn. (3948->10), ass. (0->122)
t458 = Icges(3,1) + Icges(4,2);
t457 = Icges(3,4) + Icges(4,6);
t456 = Icges(3,5) - Icges(4,4);
t455 = Icges(3,2) + Icges(4,3);
t454 = Icges(3,6) - Icges(4,5);
t453 = Icges(3,3) + Icges(4,1);
t407 = sin(pkin(5));
t408 = cos(pkin(5));
t414 = cos(qJ(2));
t415 = cos(qJ(1));
t434 = t414 * t415;
t411 = sin(qJ(2));
t412 = sin(qJ(1));
t437 = t411 * t412;
t390 = -t408 * t434 + t437;
t435 = t412 * t414;
t436 = t411 * t415;
t391 = t408 * t436 + t435;
t392 = t408 * t435 + t436;
t393 = -t408 * t437 + t434;
t438 = t407 * t415;
t439 = t407 * t412;
t444 = (-t454 * t390 + t456 * t391 - t453 * t438) * t415 + (t454 * t392 - t456 * t393 - t453 * t439) * t412;
t452 = t444 * t407;
t451 = t455 * t392 - t457 * t393 - t454 * t439;
t450 = t455 * t390 - t457 * t391 + t454 * t438;
t449 = -t457 * t392 + t458 * t393 + t456 * t439;
t448 = t457 * t390 - t458 * t391 + t456 * t438;
t447 = t453 * t408 + (t456 * t411 + t454 * t414) * t407;
t446 = t454 * t408 + (t457 * t411 + t455 * t414) * t407;
t445 = t456 * t408 + (t458 * t411 + t457 * t414) * t407;
t441 = cos(qJ(4));
t440 = t407 * t411;
t359 = pkin(2) * t391 + qJ(3) * t390;
t360 = pkin(2) * t393 + qJ(3) * t392;
t431 = qJD(2) * t407;
t404 = t412 * t431;
t427 = t415 * t431;
t433 = t359 * t404 + t360 * t427;
t371 = qJD(4) * t393 + t404;
t432 = qJD(1) * (pkin(1) * t412 - pkin(7) * t438);
t430 = qJD(3) * t414;
t405 = qJD(2) * t408 + qJD(1);
t396 = qJD(1) * (pkin(1) * t415 + pkin(7) * t439);
t429 = qJD(3) * t390 + t405 * t360 + t396;
t428 = t407 * t441;
t395 = qJD(4) * t440 + t405;
t426 = qJD(3) * t392 - t432;
t394 = (pkin(2) * t411 - qJ(3) * t414) * t407;
t423 = (-rSges(4,1) * t408 - (-rSges(4,2) * t411 - rSges(4,3) * t414) * t407 - t394) * t431;
t422 = (-pkin(3) * t408 - pkin(8) * t440 - t394) * t431;
t372 = qJD(4) * t391 - t427;
t373 = pkin(3) * t439 + pkin(8) * t393;
t374 = -pkin(3) * t438 + pkin(8) * t391;
t419 = t373 * t427 + t374 * t404 - t407 * t430 + t433;
t418 = t405 * t373 + t412 * t422 + t429;
t417 = (-t359 - t374) * t405 + t415 * t422 + t426;
t413 = cos(qJ(5));
t410 = sin(qJ(4));
t409 = sin(qJ(5));
t400 = rSges(2,1) * t415 - rSges(2,2) * t412;
t399 = rSges(2,1) * t412 + rSges(2,2) * t415;
t389 = -t407 * t414 * t410 + t408 * t441;
t388 = t408 * t410 + t414 * t428;
t381 = rSges(3,3) * t408 + (rSges(3,1) * t411 + rSges(3,2) * t414) * t407;
t370 = t390 * t410 - t415 * t428;
t369 = t390 * t441 + t410 * t438;
t368 = t392 * t410 + t412 * t428;
t367 = -t392 * t441 + t410 * t439;
t366 = t389 * t413 + t409 * t440;
t365 = -t389 * t409 + t413 * t440;
t361 = qJD(5) * t388 + t395;
t358 = pkin(4) * t389 + pkin(9) * t388;
t357 = rSges(3,1) * t393 - rSges(3,2) * t392 + rSges(3,3) * t439;
t356 = rSges(3,1) * t391 - rSges(3,2) * t390 - rSges(3,3) * t438;
t355 = -rSges(4,1) * t438 - rSges(4,2) * t391 + rSges(4,3) * t390;
t354 = rSges(4,1) * t439 - rSges(4,2) * t393 + rSges(4,3) * t392;
t338 = rSges(5,1) * t389 - rSges(5,2) * t388 + rSges(5,3) * t440;
t337 = Icges(5,1) * t389 - Icges(5,4) * t388 + Icges(5,5) * t440;
t336 = Icges(5,4) * t389 - Icges(5,2) * t388 + Icges(5,6) * t440;
t335 = Icges(5,5) * t389 - Icges(5,6) * t388 + Icges(5,3) * t440;
t334 = t370 * t413 + t391 * t409;
t333 = -t370 * t409 + t391 * t413;
t332 = t368 * t413 + t393 * t409;
t331 = -t368 * t409 + t393 * t413;
t330 = -qJD(5) * t369 + t372;
t329 = qJD(5) * t367 + t371;
t328 = pkin(4) * t370 - pkin(9) * t369;
t327 = pkin(4) * t368 + pkin(9) * t367;
t326 = rSges(5,1) * t370 + rSges(5,2) * t369 + rSges(5,3) * t391;
t325 = rSges(5,1) * t368 - rSges(5,2) * t367 + rSges(5,3) * t393;
t324 = Icges(5,1) * t370 + Icges(5,4) * t369 + Icges(5,5) * t391;
t323 = Icges(5,1) * t368 - Icges(5,4) * t367 + Icges(5,5) * t393;
t322 = Icges(5,4) * t370 + Icges(5,2) * t369 + Icges(5,6) * t391;
t321 = Icges(5,4) * t368 - Icges(5,2) * t367 + Icges(5,6) * t393;
t320 = Icges(5,5) * t370 + Icges(5,6) * t369 + Icges(5,3) * t391;
t319 = Icges(5,5) * t368 - Icges(5,6) * t367 + Icges(5,3) * t393;
t318 = rSges(6,1) * t366 + rSges(6,2) * t365 + rSges(6,3) * t388;
t317 = Icges(6,1) * t366 + Icges(6,4) * t365 + Icges(6,5) * t388;
t316 = Icges(6,4) * t366 + Icges(6,2) * t365 + Icges(6,6) * t388;
t315 = Icges(6,5) * t366 + Icges(6,6) * t365 + Icges(6,3) * t388;
t314 = t357 * t405 - t381 * t404 + t396;
t313 = -t356 * t405 - t381 * t427 - t432;
t312 = (t356 * t412 + t357 * t415) * t431;
t311 = rSges(6,1) * t334 + rSges(6,2) * t333 - rSges(6,3) * t369;
t310 = rSges(6,1) * t332 + rSges(6,2) * t331 + rSges(6,3) * t367;
t309 = Icges(6,1) * t334 + Icges(6,4) * t333 - Icges(6,5) * t369;
t308 = Icges(6,1) * t332 + Icges(6,4) * t331 + Icges(6,5) * t367;
t307 = Icges(6,4) * t334 + Icges(6,2) * t333 - Icges(6,6) * t369;
t306 = Icges(6,4) * t332 + Icges(6,2) * t331 + Icges(6,6) * t367;
t305 = Icges(6,5) * t334 + Icges(6,6) * t333 - Icges(6,3) * t369;
t304 = Icges(6,5) * t332 + Icges(6,6) * t331 + Icges(6,3) * t367;
t303 = t354 * t405 + t412 * t423 + t429;
t302 = (-t355 - t359) * t405 + t415 * t423 + t426;
t301 = (-t430 + (t354 * t415 + t355 * t412) * qJD(2)) * t407 + t433;
t300 = t325 * t395 - t338 * t371 + t418;
t299 = -t326 * t395 + t338 * t372 + t417;
t298 = -t325 * t372 + t326 * t371 + t419;
t297 = t310 * t361 - t318 * t329 + t327 * t395 - t358 * t371 + t418;
t296 = -t311 * t361 + t318 * t330 - t328 * t395 + t358 * t372 + t417;
t295 = -t310 * t330 + t311 * t329 - t327 * t372 + t328 * t371 + t419;
t1 = m(3) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(4) * (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(5) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + t371 * ((t393 * t319 - t367 * t321 + t368 * t323) * t371 + (t320 * t393 - t322 * t367 + t324 * t368) * t372 + (t335 * t393 - t336 * t367 + t337 * t368) * t395) / 0.2e1 + t372 * ((t319 * t391 + t321 * t369 + t323 * t370) * t371 + (t391 * t320 + t369 * t322 + t370 * t324) * t372 + (t335 * t391 + t336 * t369 + t337 * t370) * t395) / 0.2e1 + t395 * ((t319 * t440 - t321 * t388 + t323 * t389) * t371 + (t320 * t440 - t322 * t388 + t324 * t389) * t372 + (t335 * t440 - t388 * t336 + t389 * t337) * t395) / 0.2e1 + m(6) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + t329 * ((t367 * t304 + t331 * t306 + t332 * t308) * t329 + (t305 * t367 + t307 * t331 + t309 * t332) * t330 + (t315 * t367 + t316 * t331 + t317 * t332) * t361) / 0.2e1 + t330 * ((-t304 * t369 + t306 * t333 + t308 * t334) * t329 + (-t369 * t305 + t333 * t307 + t334 * t309) * t330 + (-t315 * t369 + t316 * t333 + t317 * t334) * t361) / 0.2e1 + t361 * ((t304 * t388 + t306 * t365 + t308 * t366) * t329 + (t305 * t388 + t307 * t365 + t309 * t366) * t330 + (t388 * t315 + t365 * t316 + t366 * t317) * t361) / 0.2e1 + ((-t444 * t408 + ((t448 * t411 + t450 * t414) * t415 + (t449 * t411 - t451 * t414) * t412) * t407) * t431 + (t447 * t408 + (t445 * t411 + t446 * t414) * t407) * t405) * t405 / 0.2e1 + (m(2) * (t399 ^ 2 + t400 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t450 * t392 + t448 * t393) * t415 + (t451 * t392 + t449 * t393 - t452) * t412) * t431 + (-t446 * t392 + t445 * t393 + t447 * t439) * t405) * t404 / 0.2e1 - (((-t450 * t390 + t448 * t391 + t452) * t415 + (t451 * t390 + t449 * t391) * t412) * t431 + (-t446 * t390 + t445 * t391 - t447 * t438) * t405) * t427 / 0.2e1;
T = t1;
