% Calculate kinetic energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:52
% EndTime: 2019-03-09 04:53:53
% DurationCPUTime: 1.88s
% Computational Cost: add. (771->197), mult. (1822->298), div. (0->0), fcn. (1838->6), ass. (0->105)
t457 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t456 = Icges(5,4) + Icges(6,6) - Icges(7,6);
t455 = Icges(5,5) + Icges(7,5) - Icges(6,4);
t454 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t453 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t452 = Icges(5,3) + Icges(7,1) + Icges(6,1);
t451 = rSges(7,1) + pkin(5);
t450 = rSges(7,3) + qJ(6);
t400 = sin(qJ(3));
t402 = cos(qJ(4));
t404 = cos(qJ(1));
t428 = t404 * t402;
t399 = sin(qJ(4));
t401 = sin(qJ(1));
t433 = t401 * t399;
t373 = t400 * t433 - t428;
t432 = t401 * t402;
t434 = t399 * t404;
t374 = t400 * t432 + t434;
t403 = cos(qJ(3));
t431 = t401 * t403;
t449 = -t456 * t373 + t457 * t374 - t455 * t431;
t375 = t400 * t434 + t432;
t376 = t400 * t428 - t433;
t429 = t403 * t404;
t448 = t456 * t375 - t457 * t376 + t455 * t429;
t447 = t454 * t373 - t456 * t374 + t453 * t431;
t446 = -t454 * t375 + t456 * t376 - t453 * t429;
t445 = -t453 * t373 + t455 * t374 - t452 * t431;
t444 = t453 * t375 - t455 * t376 + t452 * t429;
t443 = (-t453 * t399 + t455 * t402) * t403 + t452 * t400;
t442 = (-t454 * t399 + t456 * t402) * t403 + t453 * t400;
t441 = (-t456 * t399 + t457 * t402) * t403 + t455 * t400;
t436 = Icges(4,4) * t400;
t435 = Icges(4,4) * t403;
t430 = t402 * t403;
t427 = rSges(7,2) * t373 + t450 * t374 - t431 * t451;
t426 = -rSges(7,2) * t375 - t450 * t376 + t429 * t451;
t425 = (rSges(7,2) * t399 + rSges(7,3) * t402) * t403 + qJ(6) * t430 + t451 * t400;
t384 = qJD(1) * (pkin(1) * t404 + qJ(2) * t401);
t424 = qJD(1) * t404 * pkin(7) + t384;
t423 = qJD(3) * t401;
t422 = qJD(3) * t404;
t421 = qJD(4) * t403;
t388 = pkin(1) * t401 - qJ(2) * t404;
t420 = -pkin(7) * t401 - t388;
t419 = pkin(3) * t400 - pkin(8) * t403;
t378 = t419 * t401;
t379 = t419 * t404;
t418 = -t378 * t423 - t379 * t422;
t417 = rSges(4,1) * t400 + rSges(4,2) * t403;
t416 = Icges(4,1) * t400 + t435;
t415 = Icges(4,2) * t403 + t436;
t414 = Icges(4,5) * t400 + Icges(4,6) * t403;
t353 = Icges(4,6) * t404 + t401 * t415;
t356 = Icges(4,5) * t404 + t401 * t416;
t413 = -t353 * t403 - t356 * t400;
t354 = Icges(4,6) * t401 - t404 * t415;
t357 = Icges(4,5) * t401 - t404 * t416;
t412 = t354 * t403 + t357 * t400;
t386 = -Icges(4,2) * t400 + t435;
t387 = Icges(4,1) * t403 - t436;
t411 = t386 * t403 + t387 * t400;
t343 = -pkin(4) * t376 - qJ(5) * t375;
t382 = -t401 * t421 + t422;
t410 = qJD(5) * t403 * t399 + t382 * t343 + t418;
t392 = pkin(3) * t403 + pkin(8) * t400;
t398 = qJD(2) * t401;
t409 = t392 * t423 + t398 + (t379 + t420) * qJD(1);
t408 = qJD(1) * t378 + (-qJD(3) * t392 - qJD(2)) * t404 + t424;
t377 = (pkin(4) * t402 + qJ(5) * t399) * t403;
t381 = t404 * t421 + t423;
t407 = qJD(5) * t373 + t381 * t377 + t409;
t342 = pkin(4) * t374 + qJ(5) * t373;
t394 = qJD(4) * t400 + qJD(1);
t406 = -qJD(5) * t375 + t394 * t342 + t408;
t391 = rSges(2,1) * t404 - rSges(2,2) * t401;
t390 = rSges(4,1) * t403 - rSges(4,2) * t400;
t389 = rSges(2,1) * t401 + rSges(2,2) * t404;
t385 = Icges(4,5) * t403 - Icges(4,6) * t400;
t368 = rSges(6,1) * t400 + (-rSges(6,2) * t402 + rSges(6,3) * t399) * t403;
t366 = rSges(4,3) * t401 - t404 * t417;
t365 = rSges(5,3) * t400 + (rSges(5,1) * t402 - rSges(5,2) * t399) * t403;
t364 = rSges(4,3) * t404 + t401 * t417;
t351 = Icges(4,3) * t401 - t404 * t414;
t350 = Icges(4,3) * t404 + t401 * t414;
t345 = t384 - qJD(2) * t404 + qJD(1) * (-rSges(3,2) * t404 + rSges(3,3) * t401);
t344 = t398 + (rSges(3,2) * t401 + rSges(3,3) * t404 - t388) * qJD(1);
t341 = -rSges(5,1) * t376 + rSges(5,2) * t375 + rSges(5,3) * t429;
t340 = rSges(5,1) * t374 - rSges(5,2) * t373 - rSges(5,3) * t431;
t339 = rSges(6,1) * t429 + rSges(6,2) * t376 - rSges(6,3) * t375;
t337 = -rSges(6,1) * t431 - rSges(6,2) * t374 + rSges(6,3) * t373;
t315 = (-t364 * t401 + t366 * t404) * qJD(3);
t314 = qJD(1) * t364 + (-qJD(3) * t390 - qJD(2)) * t404 + t424;
t313 = t390 * t423 + t398 + (-t366 + t420) * qJD(1);
t312 = t340 * t394 - t365 * t382 + t408;
t311 = -t341 * t394 + t365 * t381 + t409;
t310 = -t340 * t381 + t341 * t382 + t418;
t309 = t337 * t394 + (-t368 - t377) * t382 + t406;
t308 = t368 * t381 + (-t339 - t343) * t394 + t407;
t307 = t339 * t382 + (-t337 - t342) * t381 + t410;
t306 = -qJD(6) * t376 + t427 * t394 + (-t377 - t425) * t382 + t406;
t305 = qJD(6) * t374 + t425 * t381 + (-t343 - t426) * t394 + t407;
t304 = qJD(6) * t430 + t426 * t382 + (-t342 - t427) * t381 + t410;
t1 = ((t401 * t385 - t404 * t411) * qJD(1) + (t401 ^ 2 * t351 + (t413 * t404 + (t350 - t412) * t401) * t404) * qJD(3)) * t423 / 0.2e1 + qJD(1) * ((-t400 * t386 + t403 * t387) * qJD(1) + ((-t353 * t400 + t356 * t403) * t404 + (-t354 * t400 + t357 * t403) * t401) * qJD(3)) / 0.2e1 + m(4) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(3) * (t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(6) * (t307 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + m(7) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + m(5) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + ((t404 * t385 + t401 * t411) * qJD(1) + (t404 ^ 2 * t350 + (t412 * t401 + (t351 - t413) * t404) * t401) * qJD(3)) * t422 / 0.2e1 + ((t442 * t375 - t441 * t376 + t443 * t429) * t394 + (-t447 * t375 - t449 * t376 + t445 * t429) * t382 + (-t446 * t375 - t448 * t376 + t444 * t429) * t381) * t381 / 0.2e1 + ((-t442 * t373 + t441 * t374 - t443 * t431) * t394 + (t447 * t373 + t449 * t374 - t445 * t431) * t382 + (t446 * t373 + t448 * t374 - t444 * t431) * t381) * t382 / 0.2e1 + (((-t442 * t399 + t441 * t402) * t394 + (t447 * t399 + t449 * t402) * t382 + (t446 * t399 + t448 * t402) * t381) * t403 + (t444 * t381 + t445 * t382 + t443 * t394) * t400) * t394 / 0.2e1 + (m(2) * (t389 ^ 2 + t391 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
