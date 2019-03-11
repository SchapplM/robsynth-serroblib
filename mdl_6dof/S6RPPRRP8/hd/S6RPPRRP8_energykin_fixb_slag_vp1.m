% Calculate kinetic energy for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:14
% EndTime: 2019-03-09 02:15:15
% DurationCPUTime: 1.35s
% Computational Cost: add. (981->186), mult. (1307->287), div. (0->0), fcn. (1280->8), ass. (0->102)
t466 = Icges(6,1) + Icges(7,1);
t465 = Icges(6,4) - Icges(7,5);
t464 = Icges(7,4) + Icges(6,5);
t463 = Icges(6,2) + Icges(7,3);
t462 = Icges(7,6) - Icges(6,6);
t461 = Icges(6,3) + Icges(7,2);
t460 = rSges(7,1) + pkin(5);
t459 = rSges(7,3) + qJ(6);
t404 = pkin(9) + qJ(4);
t399 = sin(t404);
t410 = cos(qJ(5));
t411 = cos(qJ(1));
t436 = t411 * t410;
t408 = sin(qJ(5));
t409 = sin(qJ(1));
t439 = t408 * t409;
t381 = t399 * t439 - t436;
t437 = t409 * t410;
t438 = t408 * t411;
t382 = t399 * t437 + t438;
t400 = cos(t404);
t441 = t400 * t409;
t458 = t463 * t381 - t465 * t382 - t462 * t441;
t383 = t399 * t438 + t437;
t384 = -t399 * t436 + t439;
t440 = t400 * t411;
t457 = -t463 * t383 - t465 * t384 + t462 * t440;
t456 = t462 * t381 + t464 * t382 - t461 * t441;
t455 = -t462 * t383 + t464 * t384 + t461 * t440;
t454 = -t465 * t381 + t466 * t382 - t464 * t441;
t453 = t465 * t383 + t466 * t384 + t464 * t440;
t452 = (t463 * t408 - t465 * t410) * t400 + t462 * t399;
t451 = (t462 * t408 + t464 * t410) * t400 + t461 * t399;
t450 = (-t465 * t408 + t466 * t410) * t400 + t464 * t399;
t449 = qJD(1) * t411 * qJ(3) + qJD(3) * t409;
t405 = sin(pkin(9));
t444 = pkin(3) * t405;
t443 = Icges(5,4) * t399;
t442 = Icges(5,4) * t400;
t434 = -rSges(7,2) * t441 + t459 * t381 + t460 * t382;
t433 = rSges(7,2) * t440 - t459 * t383 + t460 * t384;
t432 = rSges(7,2) * t399 + (t459 * t408 + t460 * t410) * t400;
t403 = qJD(2) * t409;
t431 = qJD(3) * t411 + t403;
t430 = qJD(4) * t409;
t429 = qJD(4) * t411;
t428 = qJD(5) * t400;
t392 = qJD(1) * (pkin(1) * t411 + qJ(2) * t409);
t427 = -qJD(2) * t411 + t392;
t426 = qJD(1) * (pkin(7) * t411 + t409 * t444) + t392 + t449;
t425 = pkin(4) * t399 - pkin(8) * t400;
t393 = pkin(1) * t409 - qJ(2) * t411;
t424 = t411 * t444 - t393 + (-pkin(7) - qJ(3)) * t409;
t377 = t425 * t409;
t378 = t425 * t411;
t423 = -t377 * t430 - t378 * t429;
t406 = cos(pkin(9));
t422 = rSges(4,1) * t405 + rSges(4,2) * t406;
t421 = rSges(5,1) * t399 + rSges(5,2) * t400;
t420 = Icges(5,1) * t399 + t442;
t419 = Icges(5,2) * t400 + t443;
t418 = Icges(5,5) * t399 + Icges(5,6) * t400;
t367 = Icges(5,6) * t411 + t419 * t409;
t369 = Icges(5,5) * t411 + t420 * t409;
t417 = -t367 * t400 - t369 * t399;
t368 = Icges(5,6) * t409 - t419 * t411;
t370 = Icges(5,5) * t409 - t420 * t411;
t416 = t368 * t400 + t370 * t399;
t388 = -Icges(5,2) * t399 + t442;
t389 = Icges(5,1) * t400 - t443;
t415 = t388 * t400 + t389 * t399;
t391 = pkin(4) * t400 + pkin(8) * t399;
t414 = t391 * t430 + (t378 + t424) * qJD(1) + t431;
t413 = qJD(1) * t377 + (-qJD(4) * t391 - qJD(2)) * t411 + t426;
t396 = qJD(5) * t399 + qJD(1);
t395 = rSges(2,1) * t411 - rSges(2,2) * t409;
t394 = rSges(2,1) * t409 + rSges(2,2) * t411;
t390 = rSges(5,1) * t400 - rSges(5,2) * t399;
t387 = Icges(5,5) * t400 - Icges(5,6) * t399;
t386 = -t409 * t428 + t429;
t385 = t411 * t428 + t430;
t372 = rSges(5,3) * t409 - t421 * t411;
t371 = rSges(5,3) * t411 + t421 * t409;
t366 = Icges(5,3) * t409 - t418 * t411;
t365 = Icges(5,3) * t411 + t418 * t409;
t364 = rSges(6,3) * t399 + (rSges(6,1) * t410 - rSges(6,2) * t408) * t400;
t356 = qJD(1) * (-rSges(3,2) * t411 + rSges(3,3) * t409) + t427;
t355 = t403 + (rSges(3,2) * t409 + rSges(3,3) * t411 - t393) * qJD(1);
t352 = rSges(6,1) * t384 + rSges(6,2) * t383 + rSges(6,3) * t440;
t350 = rSges(6,1) * t382 - rSges(6,2) * t381 - rSges(6,3) * t441;
t336 = qJD(1) * (rSges(4,3) * t411 + t422 * t409) + t427 + t449;
t335 = (-t393 + t422 * t411 + (-rSges(4,3) - qJ(3)) * t409) * qJD(1) + t431;
t334 = (-t371 * t409 + t372 * t411) * qJD(4);
t333 = qJD(1) * t371 + (-qJD(4) * t390 - qJD(2)) * t411 + t426;
t332 = t390 * t430 + (-t372 + t424) * qJD(1) + t431;
t331 = -t350 * t385 + t352 * t386 + t423;
t330 = t350 * t396 - t364 * t386 + t413;
t329 = -t352 * t396 + t364 * t385 + t414;
t328 = -qJD(6) * t383 - t432 * t386 + t434 * t396 + t413;
t327 = qJD(6) * t381 + t432 * t385 - t433 * t396 + t414;
t326 = qJD(6) * t400 * t408 - t434 * t385 + t433 * t386 + t423;
t1 = ((t411 * t387 + t415 * t409) * qJD(1) + (t411 ^ 2 * t365 + (t416 * t409 + (t366 - t417) * t411) * t409) * qJD(4)) * t429 / 0.2e1 + m(3) * (t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(5) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(6) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(7) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + ((t409 * t387 - t415 * t411) * qJD(1) + (t409 ^ 2 * t366 + (t417 * t411 + (t365 - t416) * t409) * t411) * qJD(4)) * t430 / 0.2e1 + qJD(1) * ((-t399 * t388 + t400 * t389) * qJD(1) + ((-t367 * t399 + t369 * t400) * t411 + (-t368 * t399 + t370 * t400) * t409) * qJD(4)) / 0.2e1 + ((-t383 * t452 + t384 * t450 + t440 * t451) * t396 + (-t383 * t458 + t454 * t384 + t456 * t440) * t386 + (-t457 * t383 + t453 * t384 + t455 * t440) * t385) * t385 / 0.2e1 + ((t381 * t452 + t382 * t450 - t441 * t451) * t396 + (t458 * t381 + t454 * t382 - t456 * t441) * t386 + (t381 * t457 + t382 * t453 - t441 * t455) * t385) * t386 / 0.2e1 + (((t408 * t452 + t410 * t450) * t396 + (t408 * t458 + t454 * t410) * t386 + (t408 * t457 + t410 * t453) * t385) * t400 + (t385 * t455 + t386 * t456 + t396 * t451) * t399) * t396 / 0.2e1 + (Icges(2,3) + Icges(3,1) + Icges(4,1) * t406 ^ 2 + (-0.2e1 * Icges(4,4) * t406 + Icges(4,2) * t405) * t405 + m(2) * (t394 ^ 2 + t395 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
