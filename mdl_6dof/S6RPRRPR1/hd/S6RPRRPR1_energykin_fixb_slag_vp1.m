% Calculate kinetic energy for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:30
% EndTime: 2019-03-09 04:57:31
% DurationCPUTime: 1.76s
% Computational Cost: add. (1855->244), mult. (1331->375), div. (0->0), fcn. (1208->12), ass. (0->139)
t477 = Icges(5,3) + Icges(6,3);
t410 = qJ(3) + qJ(4);
t404 = pkin(11) + t410;
t398 = sin(t404);
t399 = cos(t404);
t405 = sin(t410);
t406 = cos(t410);
t476 = Icges(5,5) * t406 + Icges(6,5) * t399 - Icges(5,6) * t405 - Icges(6,6) * t398;
t409 = qJ(1) + pkin(10);
t402 = sin(t409);
t403 = cos(t409);
t460 = Icges(6,4) * t399;
t433 = -Icges(6,2) * t398 + t460;
t335 = -Icges(6,6) * t403 + t433 * t402;
t336 = Icges(6,6) * t402 + t433 * t403;
t461 = Icges(6,4) * t398;
t436 = Icges(6,1) * t399 - t461;
t337 = -Icges(6,5) * t403 + t436 * t402;
t338 = Icges(6,5) * t402 + t436 * t403;
t462 = Icges(5,4) * t406;
t434 = -Icges(5,2) * t405 + t462;
t345 = -Icges(5,6) * t403 + t434 * t402;
t346 = Icges(5,6) * t402 + t434 * t403;
t463 = Icges(5,4) * t405;
t437 = Icges(5,1) * t406 - t463;
t347 = -Icges(5,5) * t403 + t437 * t402;
t348 = Icges(5,5) * t402 + t437 * t403;
t375 = Icges(6,2) * t399 + t461;
t376 = Icges(6,1) * t398 + t460;
t397 = qJD(3) * t402;
t380 = qJD(4) * t402 + t397;
t381 = (-qJD(3) - qJD(4)) * t403;
t384 = Icges(5,2) * t406 + t463;
t385 = Icges(5,1) * t405 + t462;
t475 = (-t335 * t398 + t337 * t399 - t345 * t405 + t347 * t406) * t381 + (-t336 * t398 + t338 * t399 - t346 * t405 + t348 * t406) * t380 + (-t375 * t398 + t376 * t399 - t384 * t405 + t385 * t406) * qJD(1);
t474 = (t476 * t402 - t477 * t403) * t381 + (t477 * t402 + t476 * t403) * t380 + (Icges(5,5) * t405 + Icges(6,5) * t398 + Icges(5,6) * t406 + Icges(6,6) * t399) * qJD(1);
t413 = sin(qJ(1));
t470 = pkin(1) * t413;
t469 = pkin(4) * t405;
t415 = cos(qJ(3));
t467 = t415 * pkin(3);
t412 = sin(qJ(3));
t465 = Icges(4,4) * t412;
t464 = Icges(4,4) * t415;
t459 = t398 * t402;
t458 = t398 * t403;
t411 = sin(qJ(6));
t457 = t402 * t411;
t414 = cos(qJ(6));
t456 = t402 * t414;
t455 = t403 * t411;
t454 = t403 * t414;
t416 = cos(qJ(1));
t401 = qJD(1) * t416 * pkin(1);
t453 = qJD(1) * (pkin(2) * t403 + pkin(7) * t402) + t401;
t452 = pkin(4) * t406;
t450 = qJD(3) * t403;
t449 = qJD(6) * t398;
t448 = pkin(3) * qJD(3) * t412;
t339 = -pkin(8) * t403 + t467 * t402;
t340 = pkin(8) * t402 + t467 * t403;
t447 = t339 * t397 + t340 * t450 + qJD(2);
t446 = -pkin(2) * t402 + pkin(7) * t403 - t470;
t445 = t403 * t448;
t328 = -qJ(5) * t403 + t452 * t402;
t444 = t380 * t328 + t447;
t443 = -t339 + t446;
t442 = pkin(5) * t399 + pkin(9) * t398;
t441 = rSges(4,1) * t415 - rSges(4,2) * t412;
t440 = rSges(5,1) * t406 - rSges(5,2) * t405;
t439 = rSges(6,1) * t399 - rSges(6,2) * t398;
t438 = Icges(4,1) * t415 - t465;
t435 = -Icges(4,2) * t412 + t464;
t432 = Icges(4,5) * t415 - Icges(4,6) * t412;
t359 = -Icges(4,6) * t403 + t435 * t402;
t361 = -Icges(4,5) * t403 + t438 * t402;
t429 = t359 * t412 - t361 * t415;
t360 = Icges(4,6) * t402 + t435 * t403;
t362 = Icges(4,5) * t402 + t438 * t403;
t428 = -t360 * t412 + t362 * t415;
t390 = Icges(4,2) * t415 + t465;
t391 = Icges(4,1) * t412 + t464;
t427 = -t390 * t412 + t391 * t415;
t426 = -t328 + t443;
t425 = qJD(5) * t402 + t381 * t469 - t445;
t424 = qJD(1) * t340 - t402 * t448 + t453;
t329 = qJ(5) * t402 + t452 * t403;
t421 = qJD(1) * t329 - qJD(5) * t403 + t424;
t394 = rSges(2,1) * t416 - rSges(2,2) * t413;
t393 = rSges(2,1) * t413 + rSges(2,2) * t416;
t392 = rSges(4,1) * t412 + rSges(4,2) * t415;
t389 = Icges(4,5) * t412 + Icges(4,6) * t415;
t388 = -qJD(6) * t399 + qJD(1);
t386 = rSges(5,1) * t405 + rSges(5,2) * t406;
t378 = pkin(5) * t398 - pkin(9) * t399;
t377 = rSges(6,1) * t398 + rSges(6,2) * t399;
t372 = t401 + qJD(1) * (rSges(3,1) * t403 - rSges(3,2) * t402);
t371 = (-rSges(3,1) * t402 - rSges(3,2) * t403 - t470) * qJD(1);
t370 = t399 * t454 + t457;
t369 = -t399 * t455 + t456;
t368 = t399 * t456 - t455;
t367 = -t399 * t457 - t454;
t366 = t442 * t403;
t365 = t442 * t402;
t364 = rSges(4,3) * t402 + t441 * t403;
t363 = -rSges(4,3) * t403 + t441 * t402;
t358 = Icges(4,3) * t402 + t432 * t403;
t357 = -Icges(4,3) * t403 + t432 * t402;
t356 = t402 * t449 + t381;
t355 = t403 * t449 + t380;
t354 = -rSges(7,3) * t399 + (rSges(7,1) * t414 - rSges(7,2) * t411) * t398;
t353 = rSges(5,3) * t402 + t440 * t403;
t352 = -rSges(5,3) * t403 + t440 * t402;
t351 = -Icges(7,5) * t399 + (Icges(7,1) * t414 - Icges(7,4) * t411) * t398;
t350 = -Icges(7,6) * t399 + (Icges(7,4) * t414 - Icges(7,2) * t411) * t398;
t349 = -Icges(7,3) * t399 + (Icges(7,5) * t414 - Icges(7,6) * t411) * t398;
t342 = rSges(6,3) * t402 + t439 * t403;
t341 = -rSges(6,3) * t403 + t439 * t402;
t325 = rSges(7,1) * t370 + rSges(7,2) * t369 + rSges(7,3) * t458;
t324 = rSges(7,1) * t368 + rSges(7,2) * t367 + rSges(7,3) * t459;
t323 = Icges(7,1) * t370 + Icges(7,4) * t369 + Icges(7,5) * t458;
t322 = Icges(7,1) * t368 + Icges(7,4) * t367 + Icges(7,5) * t459;
t321 = Icges(7,4) * t370 + Icges(7,2) * t369 + Icges(7,6) * t458;
t320 = Icges(7,4) * t368 + Icges(7,2) * t367 + Icges(7,6) * t459;
t319 = Icges(7,5) * t370 + Icges(7,6) * t369 + Icges(7,3) * t458;
t318 = Icges(7,5) * t368 + Icges(7,6) * t367 + Icges(7,3) * t459;
t317 = qJD(1) * t364 - t392 * t397 + t453;
t316 = -t392 * t450 + (-t363 + t446) * qJD(1);
t315 = qJD(2) + (t363 * t402 + t364 * t403) * qJD(3);
t314 = qJD(1) * t353 - t380 * t386 + t424;
t313 = -t445 + t381 * t386 + (-t352 + t443) * qJD(1);
t312 = t352 * t380 - t353 * t381 + t447;
t311 = qJD(1) * t342 + (-t377 - t469) * t380 + t421;
t310 = t377 * t381 + (-t341 + t426) * qJD(1) + t425;
t309 = t341 * t380 + (-t329 - t342) * t381 + t444;
t308 = qJD(1) * t366 + t325 * t388 - t354 * t355 + (-t378 - t469) * t380 + t421;
t307 = -t324 * t388 + t354 * t356 + t378 * t381 + (-t365 + t426) * qJD(1) + t425;
t306 = t324 * t355 - t325 * t356 + t365 * t380 + (-t329 - t366) * t381 + t444;
t1 = ((t402 * t389 + t427 * t403) * qJD(1) + (t402 ^ 2 * t358 + (t429 * t403 + (-t357 + t428) * t402) * t403) * qJD(3)) * t397 / 0.2e1 + m(6) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(7) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + m(5) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(4) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + t356 * ((t319 * t459 + t321 * t367 + t323 * t368) * t355 + (t318 * t459 + t367 * t320 + t368 * t322) * t356 + (t349 * t459 + t350 * t367 + t351 * t368) * t388) / 0.2e1 + t388 * ((-t318 * t356 - t319 * t355 - t349 * t388) * t399 + ((-t321 * t411 + t323 * t414) * t355 + (-t320 * t411 + t322 * t414) * t356 + (-t350 * t411 + t351 * t414) * t388) * t398) / 0.2e1 + t355 * ((t319 * t458 + t369 * t321 + t370 * t323) * t355 + (t318 * t458 + t320 * t369 + t322 * t370) * t356 + (t349 * t458 + t350 * t369 + t351 * t370) * t388) / 0.2e1 - ((-t403 * t389 + t427 * t402) * qJD(1) + (t403 ^ 2 * t357 + (t428 * t402 + (-t358 + t429) * t403) * t402) * qJD(3)) * t450 / 0.2e1 + (t474 * t402 + t475 * t403) * t380 / 0.2e1 + (t475 * t402 - t474 * t403) * t381 / 0.2e1 + (Icges(2,3) + m(2) * (t393 ^ 2 + t394 ^ 2) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + (((t360 * t415 + t362 * t412) * t402 - (t359 * t415 + t361 * t412) * t403) * qJD(3) + (t335 * t399 + t337 * t398 + t345 * t406 + t347 * t405) * t381 + (t336 * t399 + t338 * t398 + t346 * t406 + t348 * t405) * t380 + (t399 * t375 + t398 * t376 + t406 * t384 + t405 * t385 + t415 * t390 + t412 * t391) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
