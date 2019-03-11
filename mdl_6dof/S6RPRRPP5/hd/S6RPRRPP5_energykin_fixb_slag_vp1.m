% Calculate kinetic energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:35
% EndTime: 2019-03-09 04:42:37
% DurationCPUTime: 1.84s
% Computational Cost: add. (1423->204), mult. (1861->308), div. (0->0), fcn. (1875->8), ass. (0->112)
t477 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t476 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t475 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t474 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t473 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t472 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t471 = rSges(7,1) + pkin(5);
t470 = rSges(7,3) + qJ(6);
t412 = pkin(9) + qJ(3);
t410 = cos(t412);
t417 = sin(qJ(1));
t418 = cos(qJ(4));
t447 = t417 * t418;
t416 = sin(qJ(4));
t419 = cos(qJ(1));
t448 = t416 * t419;
t391 = t410 * t448 - t447;
t446 = t418 * t419;
t449 = t416 * t417;
t392 = t410 * t446 + t449;
t353 = pkin(4) * t392 + qJ(5) * t391;
t389 = t410 * t449 + t446;
t406 = -qJD(4) * t410 + qJD(1);
t469 = qJD(5) * t389 + t406 * t353;
t409 = sin(t412);
t386 = (pkin(4) * t418 + qJ(5) * t416) * t409;
t437 = qJD(4) * t409;
t438 = qJD(3) * t419;
t395 = t417 * t437 - t438;
t468 = qJD(5) * t391 + t395 * t386;
t390 = t410 * t447 - t448;
t451 = t409 * t417;
t467 = t473 * t389 + t475 * t390 - t472 * t451;
t450 = t409 * t419;
t466 = t473 * t391 + t475 * t392 - t472 * t450;
t465 = t474 * t389 + t476 * t390 - t473 * t451;
t464 = t474 * t391 + t476 * t392 - t473 * t450;
t463 = t476 * t389 + t477 * t390 - t475 * t451;
t462 = t476 * t391 + t477 * t392 - t475 * t450;
t461 = t472 * t410 + (t473 * t416 + t475 * t418) * t409;
t460 = t473 * t410 + (t474 * t416 + t476 * t418) * t409;
t459 = t475 * t410 + (t476 * t416 + t477 * t418) * t409;
t414 = cos(pkin(9));
t454 = pkin(2) * t414;
t453 = Icges(4,4) * t409;
t452 = Icges(4,4) * t410;
t444 = rSges(7,2) * t389 + t471 * t390 - t470 * t451;
t443 = rSges(7,2) * t391 + t471 * t392 - t470 * t450;
t442 = t470 * t410 + (rSges(7,2) * t416 + t471 * t418) * t409;
t403 = pkin(1) * t417 - qJ(2) * t419;
t441 = pkin(7) * t419 - t417 * t454 - t403;
t435 = pkin(3) * t410 + pkin(8) * t409;
t387 = t435 * t417;
t388 = t435 * t419;
t439 = qJD(3) * t417;
t440 = t387 * t439 + t388 * t438;
t352 = pkin(4) * t390 + qJ(5) * t389;
t394 = t419 * t437 + t439;
t436 = qJD(5) * t409 * t416 + t394 * t352 + t440;
t401 = qJD(1) * (pkin(1) * t419 + qJ(2) * t417);
t434 = -qJD(2) * t419 + qJD(1) * (pkin(7) * t417 + t419 * t454) + t401;
t413 = sin(pkin(9));
t433 = rSges(3,1) * t414 - rSges(3,2) * t413;
t432 = rSges(4,1) * t410 - rSges(4,2) * t409;
t431 = Icges(4,1) * t410 - t453;
t430 = -Icges(4,2) * t409 + t452;
t429 = Icges(4,5) * t410 - Icges(4,6) * t409;
t375 = -Icges(4,6) * t419 + t417 * t430;
t377 = -Icges(4,5) * t419 + t417 * t431;
t428 = t375 * t409 - t377 * t410;
t376 = Icges(4,6) * t417 + t419 * t430;
t378 = Icges(4,5) * t417 + t419 * t431;
t427 = -t376 * t409 + t378 * t410;
t397 = Icges(4,2) * t410 + t453;
t398 = Icges(4,1) * t409 + t452;
t426 = -t397 * t409 + t398 * t410;
t400 = pkin(3) * t409 - pkin(8) * t410;
t425 = -qJD(3) * t400 - qJD(6) * t409;
t424 = qJD(1) * t388 + t434;
t411 = qJD(2) * t417;
t423 = t411 + (-t387 + t441) * qJD(1);
t422 = -t400 * t439 + t424;
t421 = -t400 * t438 + t423;
t405 = rSges(2,1) * t419 - rSges(2,2) * t417;
t404 = rSges(2,1) * t417 + rSges(2,2) * t419;
t399 = rSges(4,1) * t409 + rSges(4,2) * t410;
t396 = Icges(4,5) * t409 + Icges(4,6) * t410;
t380 = rSges(4,3) * t417 + t419 * t432;
t379 = -rSges(4,3) * t419 + t417 * t432;
t374 = Icges(4,3) * t417 + t419 * t429;
t373 = -Icges(4,3) * t419 + t417 * t429;
t371 = -rSges(5,3) * t410 + (rSges(5,1) * t418 - rSges(5,2) * t416) * t409;
t370 = -rSges(6,2) * t410 + (rSges(6,1) * t418 + rSges(6,3) * t416) * t409;
t355 = qJD(1) * t417 * rSges(3,3) + t401 + (qJD(1) * t433 - qJD(2)) * t419;
t354 = t411 + (t419 * rSges(3,3) - t417 * t433 - t403) * qJD(1);
t351 = rSges(5,1) * t392 - rSges(5,2) * t391 + rSges(5,3) * t450;
t350 = rSges(6,1) * t392 + rSges(6,2) * t450 + rSges(6,3) * t391;
t348 = rSges(5,1) * t390 - rSges(5,2) * t389 + rSges(5,3) * t451;
t347 = rSges(6,1) * t390 + rSges(6,2) * t451 + rSges(6,3) * t389;
t325 = (t379 * t417 + t380 * t419) * qJD(3);
t324 = qJD(1) * t380 - t399 * t439 + t434;
t323 = -t399 * t438 + t411 + (-t379 + t441) * qJD(1);
t322 = t348 * t394 - t351 * t395 + t440;
t321 = t351 * t406 - t371 * t394 + t422;
t320 = -t348 * t406 + t371 * t395 + t421;
t319 = t350 * t406 + (-t370 - t386) * t394 + t422 + t469;
t318 = t370 * t395 + (-t347 - t352) * t406 + t421 + t468;
t317 = t347 * t394 + (-t350 - t353) * t395 + t436;
t316 = t425 * t417 + t443 * t406 + (-t386 - t442) * t394 + t424 + t469;
t315 = t425 * t419 + t442 * t395 + (-t352 - t444) * t406 + t423 + t468;
t314 = qJD(6) * t410 + t444 * t394 + (-t353 - t443) * t395 + t436;
t1 = ((t417 * t396 + t419 * t426) * qJD(1) + (t417 ^ 2 * t374 + (t428 * t419 + (-t373 + t427) * t417) * t419) * qJD(3)) * t439 / 0.2e1 - ((-t419 * t396 + t417 * t426) * qJD(1) + (t419 ^ 2 * t373 + (t427 * t417 + (-t374 + t428) * t419) * t417) * qJD(3)) * t438 / 0.2e1 + qJD(1) * ((t410 * t397 + t409 * t398) * qJD(1) + ((t376 * t410 + t378 * t409) * t417 - (t375 * t410 + t377 * t409) * t419) * qJD(3)) / 0.2e1 + m(6) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + m(7) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(5) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(4) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(3) * (t354 ^ 2 + t355 ^ 2) / 0.2e1 + ((t460 * t391 + t459 * t392 - t461 * t450) * t406 + (t465 * t391 + t463 * t392 - t467 * t450) * t395 + (t464 * t391 + t462 * t392 - t466 * t450) * t394) * t394 / 0.2e1 + ((t460 * t389 + t459 * t390 - t461 * t451) * t406 + (t465 * t389 + t463 * t390 - t467 * t451) * t395 + (t464 * t389 + t462 * t390 - t466 * t451) * t394) * t395 / 0.2e1 + ((t466 * t394 + t467 * t395 + t461 * t406) * t410 + ((t460 * t416 + t459 * t418) * t406 + (t465 * t416 + t463 * t418) * t395 + (t464 * t416 + t462 * t418) * t394) * t409) * t406 / 0.2e1 + (Icges(2,3) + Icges(3,2) * t414 ^ 2 + (Icges(3,1) * t413 + 0.2e1 * Icges(3,4) * t414) * t413 + m(2) * (t404 ^ 2 + t405 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
