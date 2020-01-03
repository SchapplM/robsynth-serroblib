% Calculate kinetic energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:57
% EndTime: 2019-12-31 22:31:59
% DurationCPUTime: 1.85s
% Computational Cost: add. (2120->307), mult. (4009->487), div. (0->0), fcn. (4774->12), ass. (0->140)
t427 = sin(qJ(2));
t428 = sin(qJ(1));
t431 = cos(qJ(2));
t432 = cos(qJ(1));
t457 = cos(pkin(5));
t443 = t432 * t457;
t406 = t427 * t428 - t431 * t443;
t407 = t427 * t443 + t428 * t431;
t424 = sin(pkin(5));
t453 = t432 * t424;
t367 = Icges(3,5) * t407 - Icges(3,6) * t406 - Icges(3,3) * t453;
t444 = t428 * t457;
t408 = t432 * t427 + t431 * t444;
t409 = -t427 * t444 + t432 * t431;
t454 = t428 * t424;
t368 = Icges(3,5) * t409 - Icges(3,6) * t408 + Icges(3,3) * t454;
t461 = (t367 * t432 - t368 * t428) * t424;
t430 = cos(qJ(3));
t459 = pkin(3) * t430;
t456 = t424 * t427;
t455 = t424 * t431;
t452 = qJ(3) + qJ(4);
t379 = pkin(2) * t407 + pkin(8) * t406;
t380 = pkin(2) * t409 + pkin(8) * t408;
t449 = qJD(2) * t424;
t419 = t428 * t449;
t446 = t432 * t449;
t451 = t379 * t419 + t380 * t446;
t391 = qJD(3) * t408 + t419;
t450 = qJD(1) * (pkin(1) * t428 - pkin(7) * t453);
t420 = qJD(2) * t457 + qJD(1);
t426 = sin(qJ(3));
t448 = t426 * t454;
t447 = t426 * t453;
t363 = qJD(4) * t408 + t391;
t445 = cos(t452);
t442 = t457 * t426;
t441 = t424 * t445;
t392 = qJD(3) * t406 - t446;
t336 = -pkin(3) * t447 + pkin(9) * t406 + t407 * t459;
t337 = pkin(3) * t448 + pkin(9) * t408 + t409 * t459;
t439 = t391 * t336 - t337 * t392 + t451;
t364 = qJD(4) * t406 + t392;
t410 = (pkin(2) * t427 - pkin(8) * t431) * t424;
t412 = qJD(1) * (pkin(1) * t432 + pkin(7) * t454);
t438 = t420 * t380 - t410 * t419 + t412;
t393 = (-qJD(3) - qJD(4)) * t455 + t420;
t437 = -t379 * t420 - t410 * t446 - t450;
t373 = pkin(3) * t442 + (-pkin(9) * t431 + t427 * t459) * t424;
t411 = -qJD(3) * t455 + t420;
t436 = t411 * t337 - t373 * t391 + t438;
t435 = -t336 * t411 + t392 * t373 + t437;
t429 = cos(qJ(5));
t425 = sin(qJ(5));
t423 = sin(t452);
t416 = rSges(2,1) * t432 - rSges(2,2) * t428;
t415 = rSges(2,1) * t428 + rSges(2,2) * t432;
t405 = t430 * t456 + t442;
t404 = -t426 * t456 + t430 * t457;
t399 = t423 * t457 + t427 * t441;
t398 = t423 * t456 - t445 * t457;
t397 = t457 * rSges(3,3) + (rSges(3,1) * t427 + rSges(3,2) * t431) * t424;
t396 = Icges(3,5) * t457 + (Icges(3,1) * t427 + Icges(3,4) * t431) * t424;
t395 = Icges(3,6) * t457 + (Icges(3,4) * t427 + Icges(3,2) * t431) * t424;
t394 = Icges(3,3) * t457 + (Icges(3,5) * t427 + Icges(3,6) * t431) * t424;
t390 = t409 * t430 + t448;
t389 = -t409 * t426 + t430 * t454;
t388 = t407 * t430 - t447;
t387 = -t407 * t426 - t430 * t453;
t386 = t409 * t445 + t423 * t454;
t385 = t409 * t423 - t428 * t441;
t384 = t407 * t445 - t423 * t453;
t383 = t407 * t423 + t432 * t441;
t382 = t399 * t429 - t425 * t455;
t381 = -t399 * t425 - t429 * t455;
t376 = rSges(3,1) * t409 - rSges(3,2) * t408 + rSges(3,3) * t454;
t375 = rSges(3,1) * t407 - rSges(3,2) * t406 - rSges(3,3) * t453;
t372 = Icges(3,1) * t409 - Icges(3,4) * t408 + Icges(3,5) * t454;
t371 = Icges(3,1) * t407 - Icges(3,4) * t406 - Icges(3,5) * t453;
t370 = Icges(3,4) * t409 - Icges(3,2) * t408 + Icges(3,6) * t454;
t369 = Icges(3,4) * t407 - Icges(3,2) * t406 - Icges(3,6) * t453;
t366 = pkin(4) * t399 + pkin(10) * t398;
t365 = rSges(4,1) * t405 + rSges(4,2) * t404 - rSges(4,3) * t455;
t362 = Icges(4,1) * t405 + Icges(4,4) * t404 - Icges(4,5) * t455;
t361 = Icges(4,4) * t405 + Icges(4,2) * t404 - Icges(4,6) * t455;
t360 = Icges(4,5) * t405 + Icges(4,6) * t404 - Icges(4,3) * t455;
t359 = qJD(5) * t398 + t393;
t358 = rSges(5,1) * t399 - rSges(5,2) * t398 - rSges(5,3) * t455;
t357 = Icges(5,1) * t399 - Icges(5,4) * t398 - Icges(5,5) * t455;
t356 = Icges(5,4) * t399 - Icges(5,2) * t398 - Icges(5,6) * t455;
t355 = Icges(5,5) * t399 - Icges(5,6) * t398 - Icges(5,3) * t455;
t354 = t386 * t429 + t408 * t425;
t353 = -t386 * t425 + t408 * t429;
t352 = t384 * t429 + t406 * t425;
t351 = -t384 * t425 + t406 * t429;
t350 = pkin(4) * t386 + pkin(10) * t385;
t349 = pkin(4) * t384 + pkin(10) * t383;
t347 = qJD(5) * t383 + t364;
t346 = qJD(5) * t385 + t363;
t345 = rSges(4,1) * t390 + rSges(4,2) * t389 + rSges(4,3) * t408;
t344 = rSges(4,1) * t388 + rSges(4,2) * t387 + rSges(4,3) * t406;
t343 = Icges(4,1) * t390 + Icges(4,4) * t389 + Icges(4,5) * t408;
t342 = Icges(4,1) * t388 + Icges(4,4) * t387 + Icges(4,5) * t406;
t341 = Icges(4,4) * t390 + Icges(4,2) * t389 + Icges(4,6) * t408;
t340 = Icges(4,4) * t388 + Icges(4,2) * t387 + Icges(4,6) * t406;
t339 = Icges(4,5) * t390 + Icges(4,6) * t389 + Icges(4,3) * t408;
t338 = Icges(4,5) * t388 + Icges(4,6) * t387 + Icges(4,3) * t406;
t335 = rSges(5,1) * t386 - rSges(5,2) * t385 + rSges(5,3) * t408;
t334 = rSges(5,1) * t384 - rSges(5,2) * t383 + rSges(5,3) * t406;
t333 = Icges(5,1) * t386 - Icges(5,4) * t385 + Icges(5,5) * t408;
t332 = Icges(5,1) * t384 - Icges(5,4) * t383 + Icges(5,5) * t406;
t331 = Icges(5,4) * t386 - Icges(5,2) * t385 + Icges(5,6) * t408;
t330 = Icges(5,4) * t384 - Icges(5,2) * t383 + Icges(5,6) * t406;
t329 = Icges(5,5) * t386 - Icges(5,6) * t385 + Icges(5,3) * t408;
t328 = Icges(5,5) * t384 - Icges(5,6) * t383 + Icges(5,3) * t406;
t327 = rSges(6,1) * t382 + rSges(6,2) * t381 + rSges(6,3) * t398;
t326 = Icges(6,1) * t382 + Icges(6,4) * t381 + Icges(6,5) * t398;
t325 = Icges(6,4) * t382 + Icges(6,2) * t381 + Icges(6,6) * t398;
t324 = Icges(6,5) * t382 + Icges(6,6) * t381 + Icges(6,3) * t398;
t323 = t376 * t420 - t397 * t419 + t412;
t322 = -t375 * t420 - t397 * t446 - t450;
t320 = (t375 * t428 + t376 * t432) * t449;
t318 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t385;
t317 = rSges(6,1) * t352 + rSges(6,2) * t351 + rSges(6,3) * t383;
t316 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t385;
t315 = Icges(6,1) * t352 + Icges(6,4) * t351 + Icges(6,5) * t383;
t314 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t385;
t313 = Icges(6,4) * t352 + Icges(6,2) * t351 + Icges(6,6) * t383;
t312 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t385;
t311 = Icges(6,5) * t352 + Icges(6,6) * t351 + Icges(6,3) * t383;
t310 = t345 * t411 - t365 * t391 + t438;
t309 = -t344 * t411 + t365 * t392 + t437;
t308 = t344 * t391 - t345 * t392 + t451;
t307 = t335 * t393 - t358 * t363 + t436;
t306 = -t334 * t393 + t358 * t364 + t435;
t305 = t334 * t363 - t335 * t364 + t439;
t304 = t318 * t359 - t327 * t346 + t350 * t393 - t363 * t366 + t436;
t303 = -t317 * t359 + t327 * t347 - t349 * t393 + t364 * t366 + t435;
t302 = t317 * t346 - t318 * t347 + t349 * t363 - t350 * t364 + t439;
t1 = m(3) * (t320 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + ((t394 * t454 - t395 * t408 + t396 * t409) * t420 + (-(-t369 * t408 + t371 * t409) * t432 + (-t370 * t408 + t372 * t409 - t461) * t428) * t449) * t419 / 0.2e1 - ((-t394 * t453 - t395 * t406 + t396 * t407) * t420 + ((-t370 * t406 + t372 * t407) * t428 + (t369 * t406 - t371 * t407 + t461) * t432) * t449) * t446 / 0.2e1 + t420 * ((t457 * t368 + (t370 * t431 + t372 * t427) * t424) * t419 - (t457 * t367 + (t369 * t431 + t371 * t427) * t424) * t446 + (t457 * t394 + (t395 * t431 + t396 * t427) * t424) * t420) / 0.2e1 + m(4) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t391 * ((t339 * t408 + t341 * t389 + t343 * t390) * t391 + (t338 * t408 + t340 * t389 + t342 * t390) * t392 + (t360 * t408 + t361 * t389 + t362 * t390) * t411) / 0.2e1 + t392 * ((t339 * t406 + t341 * t387 + t343 * t388) * t391 + (t338 * t406 + t340 * t387 + t342 * t388) * t392 + (t360 * t406 + t361 * t387 + t362 * t388) * t411) / 0.2e1 + t411 * ((-t339 * t455 + t341 * t404 + t343 * t405) * t391 + (-t338 * t455 + t340 * t404 + t342 * t405) * t392 + (-t360 * t455 + t361 * t404 + t362 * t405) * t411) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + t363 * ((t329 * t408 - t331 * t385 + t333 * t386) * t363 + (t328 * t408 - t330 * t385 + t332 * t386) * t364 + (t355 * t408 - t356 * t385 + t357 * t386) * t393) / 0.2e1 + t364 * ((t329 * t406 - t331 * t383 + t333 * t384) * t363 + (t328 * t406 - t330 * t383 + t332 * t384) * t364 + (t355 * t406 - t356 * t383 + t357 * t384) * t393) / 0.2e1 + t393 * ((-t329 * t455 - t331 * t398 + t333 * t399) * t363 + (-t328 * t455 - t330 * t398 + t332 * t399) * t364 + (-t355 * t455 - t356 * t398 + t357 * t399) * t393) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + t346 * ((t312 * t385 + t314 * t353 + t316 * t354) * t346 + (t311 * t385 + t313 * t353 + t315 * t354) * t347 + (t324 * t385 + t325 * t353 + t326 * t354) * t359) / 0.2e1 + t347 * ((t312 * t383 + t314 * t351 + t316 * t352) * t346 + (t311 * t383 + t313 * t351 + t315 * t352) * t347 + (t324 * t383 + t325 * t351 + t326 * t352) * t359) / 0.2e1 + t359 * ((t312 * t398 + t314 * t381 + t316 * t382) * t346 + (t311 * t398 + t313 * t381 + t315 * t382) * t347 + (t324 * t398 + t325 * t381 + t326 * t382) * t359) / 0.2e1 + (m(2) * (t415 ^ 2 + t416 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
