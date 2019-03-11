% Calculate kinetic energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:10
% EndTime: 2019-03-09 03:21:12
% DurationCPUTime: 1.83s
% Computational Cost: add. (1109->225), mult. (1537->331), div. (0->0), fcn. (1471->8), ass. (0->123)
t506 = Icges(6,1) + Icges(7,1);
t505 = Icges(6,4) + Icges(7,4);
t504 = Icges(7,5) + Icges(6,5);
t503 = Icges(6,2) + Icges(7,2);
t502 = Icges(7,6) + Icges(6,6);
t501 = Icges(4,3) + Icges(5,3);
t500 = Icges(7,3) + Icges(6,3);
t418 = qJ(3) + pkin(9);
t413 = sin(t418);
t414 = cos(t418);
t422 = sin(qJ(3));
t425 = cos(qJ(3));
t499 = Icges(4,5) * t422 + Icges(5,5) * t413 + Icges(4,6) * t425 + Icges(5,6) * t414;
t424 = cos(qJ(5));
t426 = cos(qJ(1));
t461 = t424 * t426;
t421 = sin(qJ(5));
t423 = sin(qJ(1));
t464 = t421 * t423;
t390 = -t413 * t464 + t461;
t462 = t423 * t424;
t463 = t421 * t426;
t391 = t413 * t462 + t463;
t466 = t414 * t423;
t498 = t502 * t390 + t504 * t391 - t500 * t466;
t392 = t413 * t463 + t462;
t393 = -t413 * t461 + t464;
t465 = t414 * t426;
t497 = t502 * t392 + t504 * t393 + t500 * t465;
t496 = t503 * t390 + t505 * t391 - t502 * t466;
t495 = t503 * t392 + t505 * t393 + t502 * t465;
t494 = t505 * t390 + t506 * t391 - t504 * t466;
t493 = t505 * t392 + t506 * t393 + t504 * t465;
t492 = (-t502 * t421 + t504 * t424) * t414 + t500 * t413;
t491 = (-t503 * t421 + t505 * t424) * t414 + t502 * t413;
t490 = (-t505 * t421 + t506 * t424) * t414 + t504 * t413;
t489 = t499 * t423 + t501 * t426;
t488 = t501 * t423 - t499 * t426;
t487 = Icges(4,5) * t425 + Icges(5,5) * t414 - Icges(4,6) * t422 - Icges(5,6) * t413;
t467 = Icges(5,4) * t414;
t397 = -Icges(5,2) * t413 + t467;
t468 = Icges(5,4) * t413;
t398 = Icges(5,1) * t414 - t468;
t469 = Icges(4,4) * t425;
t403 = -Icges(4,2) * t422 + t469;
t470 = Icges(4,4) * t422;
t404 = Icges(4,1) * t425 - t470;
t486 = t397 * t414 + t398 * t413 + t403 * t425 + t404 * t422;
t439 = Icges(5,2) * t414 + t468;
t368 = Icges(5,6) * t423 - t426 * t439;
t441 = Icges(5,1) * t413 + t467;
t370 = Icges(5,5) * t423 - t426 * t441;
t440 = Icges(4,2) * t425 + t470;
t380 = Icges(4,6) * t423 - t426 * t440;
t442 = Icges(4,1) * t422 + t469;
t382 = Icges(4,5) * t423 - t426 * t442;
t485 = t368 * t414 + t370 * t413 + t380 * t425 + t382 * t422;
t367 = Icges(5,6) * t426 + t423 * t439;
t369 = Icges(5,5) * t426 + t423 * t441;
t379 = Icges(4,6) * t426 + t423 * t440;
t381 = Icges(4,5) * t426 + t423 * t442;
t484 = -t367 * t414 - t369 * t413 - t379 * t425 - t381 * t422;
t473 = pkin(5) * t424;
t483 = -qJ(6) * t414 + t413 * t473;
t476 = pkin(3) * t422;
t475 = pkin(3) * t425;
t460 = rSges(7,1) * t391 + rSges(7,2) * t390 - rSges(7,3) * t466 + pkin(5) * t463 + t423 * t483;
t459 = rSges(7,1) * t393 + rSges(7,2) * t392 + rSges(7,3) * t465 + pkin(5) * t464 - t426 * t483;
t458 = (rSges(7,1) * t424 - rSges(7,2) * t421 + t473) * t414 + (qJ(6) + rSges(7,3)) * t413;
t401 = qJD(1) * (pkin(1) * t426 + qJ(2) * t423);
t457 = qJD(1) * t426 * pkin(7) + t401;
t456 = qJD(3) * t423;
t455 = qJD(3) * t426;
t454 = qJD(5) * t414;
t453 = qJD(6) * t414;
t417 = qJD(2) * t423;
t452 = qJD(4) * t426 + t456 * t475 + t417;
t405 = pkin(1) * t423 - qJ(2) * t426;
t449 = -pkin(7) * t423 - t405;
t389 = qJ(4) * t426 + t423 * t476;
t448 = qJD(1) * t389 + qJD(4) * t423 + t457;
t388 = qJ(4) * t423 - t426 * t476;
t447 = -t388 + t449;
t446 = pkin(4) * t413 - pkin(8) * t414;
t445 = rSges(4,1) * t422 + rSges(4,2) * t425;
t444 = rSges(5,1) * t413 + rSges(5,2) * t414;
t385 = t446 * t423;
t443 = qJD(1) * t385 + t448;
t400 = pkin(4) * t414 + pkin(8) * t413;
t430 = -qJD(2) + (-t400 - t475) * qJD(3);
t374 = t388 * t455;
t386 = t446 * t426;
t429 = -t386 * t455 + t374 + (-t385 - t389) * t456;
t428 = t400 * t456 + (t386 + t447) * qJD(1) + t452;
t409 = qJD(5) * t413 + qJD(1);
t408 = rSges(2,1) * t426 - rSges(2,2) * t423;
t407 = rSges(4,1) * t425 - rSges(4,2) * t422;
t406 = rSges(2,1) * t423 + rSges(2,2) * t426;
t399 = rSges(5,1) * t414 - rSges(5,2) * t413;
t395 = -t423 * t454 + t455;
t394 = t426 * t454 + t456;
t384 = rSges(4,3) * t423 - t426 * t445;
t383 = rSges(4,3) * t426 + t423 * t445;
t372 = rSges(5,3) * t423 - t426 * t444;
t371 = rSges(5,3) * t426 + t423 * t444;
t364 = rSges(6,3) * t413 + (rSges(6,1) * t424 - rSges(6,2) * t421) * t414;
t356 = t401 - qJD(2) * t426 + qJD(1) * (-rSges(3,2) * t426 + rSges(3,3) * t423);
t355 = t417 + (rSges(3,2) * t423 + rSges(3,3) * t426 - t405) * qJD(1);
t353 = (-t383 * t423 + t384 * t426) * qJD(3);
t352 = rSges(6,1) * t393 + rSges(6,2) * t392 + rSges(6,3) * t465;
t350 = rSges(6,1) * t391 + rSges(6,2) * t390 - rSges(6,3) * t466;
t334 = qJD(1) * t383 + (-qJD(3) * t407 - qJD(2)) * t426 + t457;
t333 = t407 * t456 + t417 + (-t384 + t449) * qJD(1);
t332 = qJD(1) * t371 + (-qJD(2) + (-t399 - t475) * qJD(3)) * t426 + t448;
t331 = t399 * t456 + (-t372 + t447) * qJD(1) + t452;
t330 = t374 + (t372 * t426 + (-t371 - t389) * t423) * qJD(3);
t329 = t350 * t409 - t364 * t395 + t426 * t430 + t443;
t328 = -t352 * t409 + t364 * t394 + t428;
t327 = -t350 * t394 + t352 * t395 + t429;
t326 = t460 * t409 - t458 * t395 + (t430 + t453) * t426 + t443;
t325 = t394 * t458 - t409 * t459 - t423 * t453 + t428;
t324 = qJD(6) * t413 - t394 * t460 + t395 * t459 + t429;
t1 = m(7) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(6) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(5) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(4) * (t333 ^ 2 + t334 ^ 2 + t353 ^ 2) / 0.2e1 + m(3) * (t355 ^ 2 + t356 ^ 2) / 0.2e1 + ((t392 * t491 + t393 * t490 + t465 * t492) * t409 + (t496 * t392 + t494 * t393 + t465 * t498) * t395 + (t495 * t392 + t493 * t393 + t497 * t465) * t394) * t394 / 0.2e1 + ((t390 * t491 + t391 * t490 - t466 * t492) * t409 + (t496 * t390 + t494 * t391 - t498 * t466) * t395 + (t390 * t495 + t391 * t493 - t466 * t497) * t394) * t395 / 0.2e1 + (((-t421 * t491 + t424 * t490) * t409 + (-t421 * t496 + t424 * t494) * t395 + (-t421 * t495 + t424 * t493) * t394) * t414 + (t497 * t394 + t395 * t498 + t492 * t409) * t413) * t409 / 0.2e1 + (((-t367 * t413 + t414 * t369 - t379 * t422 + t381 * t425) * t426 + (-t368 * t413 + t370 * t414 - t380 * t422 + t425 * t382) * t423) * qJD(3) + (-t413 * t397 + t414 * t398 - t422 * t403 + t425 * t404) * qJD(1)) * qJD(1) / 0.2e1 + ((t488 * t423 ^ 2 + (t484 * t426 + (-t485 + t489) * t423) * t426) * qJD(3) + (t423 * t487 - t426 * t486) * qJD(1)) * t456 / 0.2e1 + ((t489 * t426 ^ 2 + (t485 * t423 + (-t484 + t488) * t426) * t423) * qJD(3) + (t423 * t486 + t426 * t487) * qJD(1)) * t455 / 0.2e1 + (Icges(3,1) + m(2) * (t406 ^ 2 + t408 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
