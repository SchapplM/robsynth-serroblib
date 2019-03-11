% Calculate kinetic energy for
% S6RPRRPP7
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:04
% EndTime: 2019-03-09 04:50:06
% DurationCPUTime: 1.78s
% Computational Cost: add. (769->198), mult. (1817->297), div. (0->0), fcn. (1831->6), ass. (0->107)
t457 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t456 = Icges(5,4) - Icges(7,4) - Icges(6,5);
t455 = -Icges(7,5) + Icges(6,4) + Icges(5,5);
t454 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t453 = Icges(6,6) - Icges(7,6) - Icges(5,6);
t452 = -Icges(7,3) - Icges(5,3) - Icges(6,2);
t451 = rSges(7,1) + pkin(5);
t450 = rSges(7,3) + qJ(6);
t398 = sin(qJ(3));
t400 = cos(qJ(4));
t402 = cos(qJ(1));
t428 = t402 * t400;
t397 = sin(qJ(4));
t399 = sin(qJ(1));
t433 = t397 * t399;
t372 = t398 * t433 - t428;
t431 = t399 * t400;
t432 = t397 * t402;
t373 = t398 * t431 + t432;
t341 = pkin(4) * t373 + qJ(5) * t372;
t374 = t398 * t432 + t431;
t393 = qJD(4) * t398 + qJD(1);
t449 = -qJD(5) * t374 + t393 * t341;
t401 = cos(qJ(3));
t430 = t399 * t401;
t448 = -t453 * t372 - t455 * t373 - t452 * t430;
t375 = -t398 * t428 + t433;
t429 = t401 * t402;
t447 = t453 * t374 - t455 * t375 + t452 * t429;
t446 = t454 * t372 - t456 * t373 - t453 * t430;
t445 = -t454 * t374 - t456 * t375 + t453 * t429;
t444 = -t456 * t372 + t457 * t373 - t455 * t430;
t443 = t456 * t374 + t457 * t375 + t455 * t429;
t442 = (-t453 * t397 - t455 * t400) * t401 + t452 * t398;
t441 = (t454 * t397 - t456 * t400) * t401 + t453 * t398;
t440 = (-t456 * t397 + t457 * t400) * t401 + t455 * t398;
t435 = Icges(4,4) * t398;
t434 = Icges(4,4) * t401;
t427 = rSges(7,2) * t372 + t451 * t373 + t450 * t430;
t426 = -rSges(7,2) * t374 + t451 * t375 - t450 * t429;
t425 = (rSges(7,2) * t397 + t451 * t400) * t401 - t450 * t398;
t383 = qJD(1) * (pkin(1) * t402 + qJ(2) * t399);
t424 = qJD(1) * t402 * pkin(7) + t383;
t423 = qJD(3) * t399;
t422 = qJD(3) * t402;
t421 = qJD(4) * t401;
t420 = qJD(6) * t401;
t416 = pkin(3) * t398 - pkin(8) * t401;
t377 = t416 * t399;
t419 = qJD(1) * t377 + t424;
t387 = pkin(1) * t399 - qJ(2) * t402;
t418 = -pkin(7) * t399 - t387;
t391 = pkin(3) * t401 + pkin(8) * t398;
t417 = -qJD(3) * t391 - qJD(2);
t378 = t416 * t402;
t415 = -t377 * t423 - t378 * t422;
t414 = rSges(4,1) * t398 + rSges(4,2) * t401;
t413 = Icges(4,1) * t398 + t434;
t412 = Icges(4,2) * t401 + t435;
t411 = Icges(4,5) * t398 + Icges(4,6) * t401;
t356 = Icges(4,6) * t402 + t399 * t412;
t361 = Icges(4,5) * t402 + t399 * t413;
t410 = -t356 * t401 - t361 * t398;
t357 = Icges(4,6) * t399 - t402 * t412;
t362 = Icges(4,5) * t399 - t402 * t413;
t409 = t357 * t401 + t362 * t398;
t385 = -Icges(4,2) * t398 + t434;
t386 = Icges(4,1) * t401 - t435;
t408 = t385 * t401 + t386 * t398;
t342 = pkin(4) * t375 - qJ(5) * t374;
t381 = -t399 * t421 + t422;
t407 = qJD(5) * t401 * t397 + t381 * t342 + t415;
t396 = qJD(2) * t399;
t406 = t391 * t423 + t396 + (t378 + t418) * qJD(1);
t405 = t402 * t417 + t419;
t376 = (pkin(4) * t400 + qJ(5) * t397) * t401;
t380 = t402 * t421 + t423;
t404 = qJD(5) * t372 + t380 * t376 + t406;
t390 = rSges(2,1) * t402 - rSges(2,2) * t399;
t389 = rSges(4,1) * t401 - rSges(4,2) * t398;
t388 = rSges(2,1) * t399 + rSges(2,2) * t402;
t384 = Icges(4,5) * t401 - Icges(4,6) * t398;
t367 = rSges(4,3) * t399 - t402 * t414;
t366 = rSges(5,3) * t398 + (rSges(5,1) * t400 - rSges(5,2) * t397) * t401;
t365 = rSges(6,2) * t398 + (rSges(6,1) * t400 + rSges(6,3) * t397) * t401;
t363 = rSges(4,3) * t402 + t399 * t414;
t352 = Icges(4,3) * t399 - t402 * t411;
t351 = Icges(4,3) * t402 + t399 * t411;
t344 = t383 - qJD(2) * t402 + qJD(1) * (-rSges(3,2) * t402 + rSges(3,3) * t399);
t343 = t396 + (rSges(3,2) * t399 + rSges(3,3) * t402 - t387) * qJD(1);
t340 = rSges(5,1) * t375 + rSges(5,2) * t374 + rSges(5,3) * t429;
t339 = rSges(6,1) * t375 + rSges(6,2) * t429 - rSges(6,3) * t374;
t337 = rSges(5,1) * t373 - rSges(5,2) * t372 - rSges(5,3) * t430;
t336 = rSges(6,1) * t373 - rSges(6,2) * t430 + rSges(6,3) * t372;
t314 = (-t363 * t399 + t367 * t402) * qJD(3);
t313 = qJD(1) * t363 + (-qJD(3) * t389 - qJD(2)) * t402 + t424;
t312 = t389 * t423 + t396 + (-t367 + t418) * qJD(1);
t311 = t337 * t393 - t366 * t381 + t405;
t310 = -t340 * t393 + t366 * t380 + t406;
t309 = -t337 * t380 + t340 * t381 + t415;
t308 = t336 * t393 + (-t365 - t376) * t381 + t405 + t449;
t307 = t365 * t380 + (-t339 - t342) * t393 + t404;
t306 = t339 * t381 + (-t336 - t341) * t380 + t407;
t305 = t427 * t393 + (t417 - t420) * t402 + (-t376 - t425) * t381 + t419 + t449;
t304 = t399 * t420 + t425 * t380 + (-t342 - t426) * t393 + t404;
t303 = -qJD(6) * t398 + t426 * t381 + (-t341 - t427) * t380 + t407;
t1 = ((t402 * t384 + t399 * t408) * qJD(1) + (t402 ^ 2 * t351 + (t409 * t399 + (t352 - t410) * t402) * t399) * qJD(3)) * t422 / 0.2e1 + ((t399 * t384 - t408 * t402) * qJD(1) + (t399 ^ 2 * t352 + (t410 * t402 + (t351 - t409) * t399) * t402) * qJD(3)) * t423 / 0.2e1 + qJD(1) * ((-t398 * t385 + t401 * t386) * qJD(1) + ((-t356 * t398 + t361 * t401) * t402 + (-t357 * t398 + t362 * t401) * t399) * qJD(3)) / 0.2e1 + m(6) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + m(7) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(5) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(4) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(3) * (t343 ^ 2 + t344 ^ 2) / 0.2e1 + ((-t441 * t374 + t440 * t375 - t442 * t429) * t393 + (-t446 * t374 + t444 * t375 - t448 * t429) * t381 + (-t445 * t374 + t443 * t375 - t447 * t429) * t380) * t380 / 0.2e1 + ((t441 * t372 + t440 * t373 + t442 * t430) * t393 + (t446 * t372 + t444 * t373 + t448 * t430) * t381 + (t445 * t372 + t443 * t373 + t447 * t430) * t380) * t381 / 0.2e1 + (((t441 * t397 + t440 * t400) * t393 + (t446 * t397 + t444 * t400) * t381 + (t445 * t397 + t443 * t400) * t380) * t401 + (-t447 * t380 - t448 * t381 - t442 * t393) * t398) * t393 / 0.2e1 + (Icges(2,3) + Icges(3,1) + m(2) * (t388 ^ 2 + t390 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
