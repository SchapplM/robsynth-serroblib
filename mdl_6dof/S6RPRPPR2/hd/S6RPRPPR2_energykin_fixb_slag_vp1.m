% Calculate kinetic energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:13
% EndTime: 2019-03-09 02:41:15
% DurationCPUTime: 2.29s
% Computational Cost: add. (1473->218), mult. (1269->330), div. (0->0), fcn. (1146->10), ass. (0->126)
t499 = Icges(5,4) + Icges(6,6);
t498 = Icges(5,1) + Icges(6,2);
t497 = -Icges(5,2) - Icges(6,3);
t406 = qJ(3) + pkin(10);
t404 = cos(t406);
t496 = t499 * t404;
t402 = sin(t406);
t495 = t499 * t402;
t494 = -Icges(6,4) + Icges(5,5);
t493 = Icges(6,5) - Icges(5,6);
t492 = t497 * t402 + t496;
t491 = -t498 * t404 + t495;
t407 = qJ(1) + pkin(9);
t403 = sin(t407);
t405 = cos(t407);
t490 = t492 * t403 + t493 * t405;
t489 = -t493 * t403 + t492 * t405;
t488 = t491 * t403 + t494 * t405;
t487 = t494 * t403 - t491 * t405;
t486 = t497 * t404 - t495;
t485 = t498 * t402 + t496;
t484 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t410 = sin(qJ(3));
t413 = cos(qJ(3));
t483 = Icges(4,5) * t413 - Icges(4,6) * t410 + t493 * t402 + t494 * t404;
t482 = t484 * t403 + t483 * t405;
t481 = t483 * t403 - t484 * t405;
t480 = Icges(4,5) * t410 + Icges(4,6) * t413 + t494 * t402 - t493 * t404;
t468 = Icges(4,4) * t410;
t393 = Icges(4,2) * t413 + t468;
t467 = Icges(4,4) * t413;
t394 = Icges(4,1) * t410 + t467;
t479 = -t393 * t410 + t394 * t413 + t486 * t402 + t485 * t404;
t437 = -Icges(4,2) * t410 + t467;
t359 = Icges(4,6) * t403 + t437 * t405;
t439 = Icges(4,1) * t413 - t468;
t362 = Icges(4,5) * t403 + t439 * t405;
t478 = -t359 * t410 + t362 * t413 - t489 * t402 + t487 * t404;
t358 = -Icges(4,6) * t405 + t437 * t403;
t361 = -Icges(4,5) * t405 + t439 * t403;
t477 = t358 * t410 - t361 * t413 + t490 * t402 + t488 * t404;
t473 = pkin(3) * t410;
t411 = sin(qJ(1));
t472 = t411 * pkin(1);
t470 = t413 * pkin(3);
t462 = t403 * t404;
t409 = sin(qJ(6));
t461 = t403 * t409;
t412 = cos(qJ(6));
t460 = t403 * t412;
t459 = t404 * t405;
t458 = t405 * t409;
t457 = t405 * t412;
t414 = cos(qJ(1));
t401 = qJD(1) * t414 * pkin(1);
t456 = qJD(1) * (t405 * pkin(2) + t403 * pkin(7)) + t401;
t399 = qJD(4) * t403;
t452 = qJD(5) * t402;
t455 = t405 * t452 + t399;
t454 = qJD(3) * t403;
t453 = qJD(3) * t405;
t451 = qJD(6) * t404;
t334 = -qJ(4) * t405 + t470 * t403;
t335 = qJ(4) * t403 + t470 * t405;
t450 = t334 * t454 + t335 * t453 + qJD(2);
t447 = -t403 * pkin(2) + t405 * pkin(7) - t472;
t446 = -t402 * pkin(4) + t404 * qJ(5) - t473;
t445 = -t334 + t447;
t444 = rSges(4,1) * t413 - rSges(4,2) * t410;
t443 = rSges(5,1) * t404 - rSges(5,2) * t402;
t442 = -rSges(6,2) * t404 + rSges(6,3) * t402;
t441 = pkin(4) * t404 + qJ(5) * t402;
t440 = qJD(3) * (-t402 * rSges(5,1) - t404 * rSges(5,2) - t473);
t367 = t441 * t403;
t421 = -t367 + t445;
t420 = qJD(1) * t335 - qJD(4) * t405 + t456;
t419 = qJD(3) * (t402 * rSges(6,2) + t404 * rSges(6,3) + t446);
t418 = qJD(3) * (-pkin(8) * t402 + t446);
t368 = t441 * t405;
t417 = qJD(1) * t368 + t403 * t452 + t420;
t416 = -qJD(5) * t404 + t367 * t454 + t368 * t453 + t450;
t398 = qJD(6) * t402 + qJD(1);
t397 = t414 * rSges(2,1) - t411 * rSges(2,2);
t396 = t411 * rSges(2,1) + t414 * rSges(2,2);
t395 = t410 * rSges(4,1) + t413 * rSges(4,2);
t378 = -t405 * pkin(5) + pkin(8) * t462;
t377 = t403 * pkin(5) + pkin(8) * t459;
t376 = t403 * t451 - t453;
t375 = t405 * t451 + t454;
t374 = t402 * t461 - t457;
t373 = t402 * t460 + t458;
t372 = t402 * t458 + t460;
t371 = t402 * t457 - t461;
t370 = t401 + qJD(1) * (t405 * rSges(3,1) - t403 * rSges(3,2));
t369 = (-t403 * rSges(3,1) - t405 * rSges(3,2) - t472) * qJD(1);
t365 = t403 * rSges(4,3) + t444 * t405;
t364 = t402 * rSges(7,3) + (-rSges(7,1) * t409 - rSges(7,2) * t412) * t404;
t363 = -t405 * rSges(4,3) + t444 * t403;
t360 = Icges(7,5) * t402 + (-Icges(7,1) * t409 - Icges(7,4) * t412) * t404;
t357 = Icges(7,6) * t402 + (-Icges(7,4) * t409 - Icges(7,2) * t412) * t404;
t354 = Icges(7,3) * t402 + (-Icges(7,5) * t409 - Icges(7,6) * t412) * t404;
t353 = -t405 * rSges(6,1) + t442 * t403;
t352 = t403 * rSges(6,1) + t442 * t405;
t351 = t403 * rSges(5,3) + t443 * t405;
t350 = -t405 * rSges(5,3) + t443 * t403;
t330 = t374 * rSges(7,1) + t373 * rSges(7,2) + rSges(7,3) * t462;
t329 = t372 * rSges(7,1) + t371 * rSges(7,2) + rSges(7,3) * t459;
t328 = Icges(7,1) * t374 + Icges(7,4) * t373 + Icges(7,5) * t462;
t327 = Icges(7,1) * t372 + Icges(7,4) * t371 + Icges(7,5) * t459;
t326 = Icges(7,4) * t374 + Icges(7,2) * t373 + Icges(7,6) * t462;
t325 = Icges(7,4) * t372 + Icges(7,2) * t371 + Icges(7,6) * t459;
t324 = Icges(7,5) * t374 + Icges(7,6) * t373 + Icges(7,3) * t462;
t323 = Icges(7,5) * t372 + Icges(7,6) * t371 + Icges(7,3) * t459;
t322 = qJD(1) * t365 - t395 * t454 + t456;
t321 = -t395 * t453 + (-t363 + t447) * qJD(1);
t320 = qJD(2) + (t363 * t403 + t365 * t405) * qJD(3);
t319 = qJD(1) * t351 + t403 * t440 + t420;
t318 = t399 + t405 * t440 + (-t350 + t445) * qJD(1);
t317 = (t350 * t403 + t351 * t405) * qJD(3) + t450;
t316 = qJD(1) * t352 + t403 * t419 + t417;
t315 = t405 * t419 + (-t353 + t421) * qJD(1) + t455;
t314 = (t352 * t405 + t353 * t403) * qJD(3) + t416;
t313 = qJD(1) * t377 + t398 * t329 - t375 * t364 + t403 * t418 + t417;
t312 = -t398 * t330 + t376 * t364 + t405 * t418 + (-t378 + t421) * qJD(1) + t455;
t311 = -t376 * t329 + t375 * t330 + (t377 * t405 + t378 * t403) * qJD(3) + t416;
t1 = m(7) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + t398 * ((t323 * t375 + t324 * t376 + t354 * t398) * t402 + ((-t325 * t412 - t327 * t409) * t375 + (-t326 * t412 - t328 * t409) * t376 + (-t357 * t412 - t360 * t409) * t398) * t404) / 0.2e1 + t375 * ((t323 * t459 + t371 * t325 + t372 * t327) * t375 + (t324 * t459 + t371 * t326 + t372 * t328) * t376 + (t354 * t459 + t371 * t357 + t372 * t360) * t398) / 0.2e1 + t376 * ((t323 * t462 + t373 * t325 + t374 * t327) * t375 + (t324 * t462 + t373 * t326 + t374 * t328) * t376 + (t354 * t462 + t373 * t357 + t374 * t360) * t398) / 0.2e1 + m(5) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + m(6) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(4) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t396 ^ 2 + t397 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t413 * t358 - t410 * t361 + t488 * t402 - t490 * t404) * t405 + (t413 * t359 + t410 * t362 + t487 * t402 + t489 * t404) * t403) * qJD(3) + (t413 * t393 + t410 * t394 + t485 * t402 - t486 * t404) * qJD(1)) * qJD(1) / 0.2e1 + ((t482 * t403 ^ 2 + (t477 * t405 + (t478 - t481) * t403) * t405) * qJD(3) + (t403 * t480 + t405 * t479) * qJD(1)) * t454 / 0.2e1 - ((t481 * t405 ^ 2 + (t478 * t403 + (t477 - t482) * t405) * t403) * qJD(3) + (t403 * t479 - t405 * t480) * qJD(1)) * t453 / 0.2e1;
T  = t1;
