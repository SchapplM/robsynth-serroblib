% Calculate kinetic energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:05
% EndTime: 2019-12-05 16:54:07
% DurationCPUTime: 1.94s
% Computational Cost: add. (1774->247), mult. (4507->387), div. (0->0), fcn. (5519->10), ass. (0->119)
t463 = Icges(5,1) + Icges(6,1);
t462 = Icges(5,4) + Icges(6,4);
t461 = Icges(5,5) + Icges(6,5);
t460 = Icges(5,2) + Icges(6,2);
t459 = Icges(5,6) + Icges(6,6);
t458 = Icges(5,3) + Icges(6,3);
t457 = rSges(6,3) + qJ(5);
t409 = sin(pkin(9));
t411 = cos(pkin(9));
t418 = cos(qJ(2));
t412 = cos(pkin(5));
t416 = sin(qJ(2));
t434 = t412 * t416;
t396 = t409 * t418 + t411 * t434;
t415 = sin(qJ(3));
t410 = sin(pkin(5));
t435 = t411 * t410;
t443 = cos(qJ(3));
t382 = t396 * t443 - t415 * t435;
t433 = t412 * t418;
t395 = t409 * t416 - t411 * t433;
t414 = sin(qJ(4));
t417 = cos(qJ(4));
t358 = -t382 * t414 + t395 * t417;
t440 = t395 * t414;
t359 = t382 * t417 + t440;
t426 = t410 * t443;
t381 = t396 * t415 + t411 * t426;
t456 = t459 * t358 + t461 * t359 + t458 * t381;
t398 = -t409 * t434 + t411 * t418;
t437 = t410 * t415;
t384 = t398 * t443 + t409 * t437;
t397 = t409 * t433 + t411 * t416;
t360 = -t384 * t414 + t397 * t417;
t439 = t397 * t414;
t361 = t384 * t417 + t439;
t383 = t398 * t415 - t409 * t426;
t455 = t459 * t360 + t461 * t361 + t458 * t383;
t454 = t460 * t358 + t462 * t359 + t459 * t381;
t453 = t460 * t360 + t462 * t361 + t459 * t383;
t452 = t462 * t358 + t463 * t359 + t461 * t381;
t451 = t462 * t360 + t463 * t361 + t461 * t383;
t400 = t412 * t415 + t416 * t426;
t436 = t410 * t418;
t385 = -t400 * t414 - t417 * t436;
t427 = t414 * t436;
t386 = t400 * t417 - t427;
t399 = -t412 * t443 + t416 * t437;
t450 = t459 * t385 + t461 * t386 + t458 * t399;
t449 = t460 * t385 + t462 * t386 + t459 * t399;
t448 = t462 * t385 + t463 * t386 + t461 * t399;
t447 = qJD(2) ^ 2;
t442 = t417 * pkin(4);
t438 = t409 * t410;
t432 = t359 * rSges(6,1) + t358 * rSges(6,2) + pkin(4) * t440 + t457 * t381 + t442 * t382;
t431 = t361 * rSges(6,1) + t360 * rSges(6,2) + pkin(4) * t439 + t457 * t383 + t442 * t384;
t430 = t386 * rSges(6,1) + t385 * rSges(6,2) - pkin(4) * t427 + t457 * t399 + t442 * t400;
t429 = qJD(2) * t410;
t405 = t409 * t429;
t387 = qJD(3) * t397 + t405;
t408 = qJD(2) * t412;
t425 = t411 * t429;
t377 = t396 * pkin(2) + t395 * pkin(7);
t378 = t398 * pkin(2) + t397 * pkin(7);
t424 = t377 * t405 + t378 * t425 + qJD(1);
t388 = qJD(3) * t395 - t425;
t402 = -qJD(3) * t436 + t408;
t401 = (pkin(2) * t416 - pkin(7) * t418) * t410;
t423 = t378 * t408 - t401 * t405;
t354 = t382 * pkin(3) + t381 * pkin(8);
t355 = t384 * pkin(3) + t383 * pkin(8);
t422 = t387 * t354 - t388 * t355 + t424;
t421 = (-t377 * t412 - t401 * t435) * qJD(2);
t379 = t400 * pkin(3) + t399 * pkin(8);
t420 = t402 * t355 - t387 * t379 + t423;
t419 = -t402 * t354 + t388 * t379 + t421;
t392 = t412 * rSges(3,3) + (rSges(3,1) * t416 + rSges(3,2) * t418) * t410;
t391 = Icges(3,5) * t412 + (Icges(3,1) * t416 + Icges(3,4) * t418) * t410;
t390 = Icges(3,6) * t412 + (Icges(3,4) * t416 + Icges(3,2) * t418) * t410;
t389 = Icges(3,3) * t412 + (Icges(3,5) * t416 + Icges(3,6) * t418) * t410;
t380 = qJD(4) * t399 + t402;
t375 = t400 * rSges(4,1) - t399 * rSges(4,2) - rSges(4,3) * t436;
t374 = Icges(4,1) * t400 - Icges(4,4) * t399 - Icges(4,5) * t436;
t373 = Icges(4,4) * t400 - Icges(4,2) * t399 - Icges(4,6) * t436;
t372 = Icges(4,5) * t400 - Icges(4,6) * t399 - Icges(4,3) * t436;
t369 = t398 * rSges(3,1) - t397 * rSges(3,2) + rSges(3,3) * t438;
t368 = t396 * rSges(3,1) - t395 * rSges(3,2) - rSges(3,3) * t435;
t367 = Icges(3,1) * t398 - Icges(3,4) * t397 + Icges(3,5) * t438;
t366 = Icges(3,1) * t396 - Icges(3,4) * t395 - Icges(3,5) * t435;
t365 = Icges(3,4) * t398 - Icges(3,2) * t397 + Icges(3,6) * t438;
t364 = Icges(3,4) * t396 - Icges(3,2) * t395 - Icges(3,6) * t435;
t363 = Icges(3,5) * t398 - Icges(3,6) * t397 + Icges(3,3) * t438;
t362 = Icges(3,5) * t396 - Icges(3,6) * t395 - Icges(3,3) * t435;
t357 = qJD(4) * t381 + t388;
t356 = qJD(4) * t383 + t387;
t351 = (-t368 * t412 - t392 * t435) * qJD(2);
t350 = (t369 * t412 - t392 * t438) * qJD(2);
t349 = t386 * rSges(5,1) + t385 * rSges(5,2) + t399 * rSges(5,3);
t341 = t384 * rSges(4,1) - t383 * rSges(4,2) + t397 * rSges(4,3);
t340 = t382 * rSges(4,1) - t381 * rSges(4,2) + t395 * rSges(4,3);
t339 = Icges(4,1) * t384 - Icges(4,4) * t383 + Icges(4,5) * t397;
t338 = Icges(4,1) * t382 - Icges(4,4) * t381 + Icges(4,5) * t395;
t337 = Icges(4,4) * t384 - Icges(4,2) * t383 + Icges(4,6) * t397;
t336 = Icges(4,4) * t382 - Icges(4,2) * t381 + Icges(4,6) * t395;
t335 = Icges(4,5) * t384 - Icges(4,6) * t383 + Icges(4,3) * t397;
t334 = Icges(4,5) * t382 - Icges(4,6) * t381 + Icges(4,3) * t395;
t331 = qJD(1) + (t368 * t409 + t369 * t411) * t429;
t330 = t361 * rSges(5,1) + t360 * rSges(5,2) + t383 * rSges(5,3);
t328 = t359 * rSges(5,1) + t358 * rSges(5,2) + t381 * rSges(5,3);
t312 = -t402 * t340 + t388 * t375 + t421;
t311 = t402 * t341 - t387 * t375 + t423;
t310 = t387 * t340 - t388 * t341 + t424;
t309 = -t380 * t328 + t357 * t349 + t419;
t308 = t380 * t330 - t356 * t349 + t420;
t307 = t356 * t328 - t357 * t330 + t422;
t306 = qJD(5) * t383 + t430 * t357 - t432 * t380 + t419;
t305 = qJD(5) * t381 - t430 * t356 + t431 * t380 + t420;
t304 = qJD(5) * t399 + t432 * t356 - t431 * t357 + t422;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t331 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 - t447 * ((-t363 * t435 - t395 * t365 + t396 * t367) * t438 - (-t362 * t435 - t395 * t364 + t396 * t366) * t435 + (-t389 * t435 - t395 * t390 + t396 * t391) * t412) * t435 / 0.2e1 + m(4) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + t387 * ((t397 * t335 - t383 * t337 + t384 * t339) * t387 + (t397 * t334 - t383 * t336 + t384 * t338) * t388 + (t397 * t372 - t383 * t373 + t384 * t374) * t402) / 0.2e1 + t388 * ((t395 * t335 - t381 * t337 + t382 * t339) * t387 + (t395 * t334 - t381 * t336 + t382 * t338) * t388 + (t395 * t372 - t381 * t373 + t382 * t374) * t402) / 0.2e1 + t402 * ((-t335 * t436 - t399 * t337 + t400 * t339) * t387 + (-t334 * t436 - t399 * t336 + t400 * t338) * t388 + (-t372 * t436 - t399 * t373 + t400 * t374) * t402) / 0.2e1 + m(5) * (t307 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + m(6) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + ((t449 * t360 + t448 * t361 + t450 * t383) * t380 + (t454 * t360 + t452 * t361 + t456 * t383) * t357 + (t453 * t360 + t451 * t361 + t455 * t383) * t356) * t356 / 0.2e1 + ((t449 * t358 + t448 * t359 + t450 * t381) * t380 + (t454 * t358 + t452 * t359 + t456 * t381) * t357 + (t453 * t358 + t451 * t359 + t455 * t381) * t356) * t357 / 0.2e1 + ((t449 * t385 + t448 * t386 + t450 * t399) * t380 + (t454 * t385 + t452 * t386 + t456 * t399) * t357 + (t453 * t385 + t451 * t386 + t455 * t399) * t356) * t380 / 0.2e1 + (((t363 * t438 - t397 * t365 + t398 * t367) * t438 - (t362 * t438 - t397 * t364 + t398 * t366) * t435 + (t389 * t438 - t397 * t390 + t398 * t391) * t412) * t438 + t412 * (t412 ^ 2 * t389 + (((t365 * t418 + t367 * t416) * t409 - (t364 * t418 + t366 * t416) * t411) * t410 + (-t362 * t411 + t363 * t409 + t390 * t418 + t391 * t416) * t412) * t410)) * t447 / 0.2e1;
T = t1;
