% Calculate kinetic energy for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:45
% EndTime: 2019-03-09 21:41:48
% DurationCPUTime: 2.49s
% Computational Cost: add. (2576->280), mult. (6456->421), div. (0->0), fcn. (7943->10), ass. (0->131)
t548 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t547 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t546 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t545 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t544 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t543 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t542 = rSges(7,1) + pkin(5);
t541 = rSges(7,3) + qJ(6);
t496 = sin(qJ(2));
t497 = sin(qJ(1));
t498 = cos(qJ(2));
t499 = cos(qJ(1));
t523 = cos(pkin(6));
t510 = t499 * t523;
t477 = t496 * t497 - t498 * t510;
t478 = t496 * t510 + t497 * t498;
t494 = sin(pkin(6));
t521 = t494 * t499;
t443 = Icges(3,5) * t478 - Icges(3,6) * t477 - Icges(3,3) * t521;
t511 = t497 * t523;
t479 = t499 * t496 + t498 * t511;
t480 = -t496 * t511 + t499 * t498;
t520 = t497 * t494;
t444 = Icges(3,5) * t480 - Icges(3,6) * t479 + Icges(3,3) * t520;
t540 = (t443 * t499 - t444 * t497) * t494;
t495 = sin(qJ(3));
t526 = cos(qJ(3));
t462 = t478 * t526 - t495 * t521;
t524 = sin(qJ(4));
t525 = cos(qJ(4));
t434 = t462 * t524 - t477 * t525;
t435 = t462 * t525 + t477 * t524;
t513 = t494 * t526;
t461 = t478 * t495 + t499 * t513;
t539 = t546 * t434 + t435 * t548 + t545 * t461;
t464 = t480 * t526 + t495 * t520;
t436 = t464 * t524 - t479 * t525;
t437 = t464 * t525 + t479 * t524;
t463 = t480 * t495 - t497 * t513;
t538 = t546 * t436 + t437 * t548 + t545 * t463;
t537 = t434 * t543 + t435 * t546 + t461 * t544;
t536 = t436 * t543 + t437 * t546 + t463 * t544;
t535 = t434 * t544 + t435 * t545 + t461 * t547;
t534 = t436 * t544 + t437 * t545 + t463 * t547;
t476 = t495 * t523 + t496 * t513;
t522 = t494 * t498;
t459 = t476 * t524 + t522 * t525;
t460 = t476 * t525 - t522 * t524;
t475 = t494 * t495 * t496 - t523 * t526;
t533 = t546 * t459 + t460 * t548 + t545 * t475;
t532 = t459 * t543 + t460 * t546 + t475 * t544;
t531 = t459 * t544 + t460 * t545 + t475 * t547;
t519 = rSges(7,2) * t434 + t541 * t435 + t542 * t461;
t518 = rSges(7,2) * t436 + t541 * t437 + t542 * t463;
t517 = rSges(7,2) * t459 + t541 * t460 + t542 * t475;
t455 = pkin(2) * t478 + pkin(9) * t477;
t456 = pkin(2) * t480 + pkin(9) * t479;
t514 = qJD(2) * t494;
t490 = t497 * t514;
t512 = t499 * t514;
t516 = t455 * t490 + t456 * t512;
t465 = qJD(3) * t479 + t490;
t515 = qJD(1) * (pkin(1) * t497 - pkin(8) * t521);
t491 = qJD(2) * t523 + qJD(1);
t466 = qJD(3) * t477 - t512;
t428 = pkin(3) * t462 + pkin(10) * t461;
t429 = pkin(3) * t464 + pkin(10) * t463;
t508 = t465 * t428 - t429 * t466 + t516;
t482 = -qJD(3) * t522 + t491;
t481 = (pkin(2) * t496 - pkin(9) * t498) * t494;
t483 = qJD(1) * (pkin(1) * t499 + pkin(8) * t520);
t507 = t491 * t456 - t481 * t490 + t483;
t397 = pkin(4) * t435 + qJ(5) * t434;
t430 = qJD(4) * t463 + t465;
t506 = qJD(5) * t459 + t430 * t397 + t508;
t505 = -t455 * t491 - t481 * t512 - t515;
t454 = pkin(3) * t476 + pkin(10) * t475;
t504 = t482 * t429 - t454 * t465 + t507;
t398 = pkin(4) * t437 + qJ(5) * t436;
t457 = qJD(4) * t475 + t482;
t503 = qJD(5) * t434 + t457 * t398 + t504;
t502 = -t428 * t482 + t466 * t454 + t505;
t427 = pkin(4) * t460 + qJ(5) * t459;
t431 = qJD(4) * t461 + t466;
t501 = qJD(5) * t436 + t431 * t427 + t502;
t486 = rSges(2,1) * t499 - rSges(2,2) * t497;
t485 = rSges(2,1) * t497 + rSges(2,2) * t499;
t470 = t523 * rSges(3,3) + (rSges(3,1) * t496 + rSges(3,2) * t498) * t494;
t469 = Icges(3,5) * t523 + (Icges(3,1) * t496 + Icges(3,4) * t498) * t494;
t468 = Icges(3,6) * t523 + (Icges(3,4) * t496 + Icges(3,2) * t498) * t494;
t467 = Icges(3,3) * t523 + (Icges(3,5) * t496 + Icges(3,6) * t498) * t494;
t451 = rSges(3,1) * t480 - rSges(3,2) * t479 + rSges(3,3) * t520;
t450 = rSges(3,1) * t478 - rSges(3,2) * t477 - rSges(3,3) * t521;
t448 = Icges(3,1) * t480 - Icges(3,4) * t479 + Icges(3,5) * t520;
t447 = Icges(3,1) * t478 - Icges(3,4) * t477 - Icges(3,5) * t521;
t446 = Icges(3,4) * t480 - Icges(3,2) * t479 + Icges(3,6) * t520;
t445 = Icges(3,4) * t478 - Icges(3,2) * t477 - Icges(3,6) * t521;
t442 = rSges(4,1) * t476 - rSges(4,2) * t475 - rSges(4,3) * t522;
t441 = Icges(4,1) * t476 - Icges(4,4) * t475 - Icges(4,5) * t522;
t440 = Icges(4,4) * t476 - Icges(4,2) * t475 - Icges(4,6) * t522;
t439 = Icges(4,5) * t476 - Icges(4,6) * t475 - Icges(4,3) * t522;
t424 = rSges(4,1) * t464 - rSges(4,2) * t463 + rSges(4,3) * t479;
t423 = rSges(4,1) * t462 - rSges(4,2) * t461 + rSges(4,3) * t477;
t422 = Icges(4,1) * t464 - Icges(4,4) * t463 + Icges(4,5) * t479;
t421 = Icges(4,1) * t462 - Icges(4,4) * t461 + Icges(4,5) * t477;
t420 = Icges(4,4) * t464 - Icges(4,2) * t463 + Icges(4,6) * t479;
t419 = Icges(4,4) * t462 - Icges(4,2) * t461 + Icges(4,6) * t477;
t418 = Icges(4,5) * t464 - Icges(4,6) * t463 + Icges(4,3) * t479;
t417 = Icges(4,5) * t462 - Icges(4,6) * t461 + Icges(4,3) * t477;
t416 = rSges(5,1) * t460 - rSges(5,2) * t459 + rSges(5,3) * t475;
t415 = rSges(6,1) * t475 - rSges(6,2) * t460 + rSges(6,3) * t459;
t401 = t451 * t491 - t470 * t490 + t483;
t400 = -t450 * t491 - t470 * t512 - t515;
t399 = (t450 * t497 + t451 * t499) * t514;
t395 = rSges(5,1) * t437 - rSges(5,2) * t436 + rSges(5,3) * t463;
t394 = rSges(5,1) * t435 - rSges(5,2) * t434 + rSges(5,3) * t461;
t393 = rSges(6,1) * t463 - rSges(6,2) * t437 + rSges(6,3) * t436;
t391 = rSges(6,1) * t461 - rSges(6,2) * t435 + rSges(6,3) * t434;
t369 = t424 * t482 - t442 * t465 + t507;
t368 = -t423 * t482 + t442 * t466 + t505;
t367 = t423 * t465 - t424 * t466 + t516;
t366 = t395 * t457 - t416 * t430 + t504;
t365 = -t394 * t457 + t416 * t431 + t502;
t364 = t394 * t430 - t395 * t431 + t508;
t363 = t393 * t457 + (-t415 - t427) * t430 + t503;
t362 = t415 * t431 + (-t391 - t397) * t457 + t501;
t361 = t391 * t430 + (-t393 - t398) * t431 + t506;
t360 = qJD(6) * t435 + t518 * t457 + (-t427 - t517) * t430 + t503;
t359 = qJD(6) * t437 + t517 * t431 + (-t397 - t519) * t457 + t501;
t358 = qJD(6) * t460 + t519 * t430 + (-t398 - t518) * t431 + t506;
t1 = ((t467 * t520 - t468 * t479 + t469 * t480) * t491 + (-(-t445 * t479 + t447 * t480) * t499 + (-t479 * t446 + t480 * t448 - t540) * t497) * t514) * t490 / 0.2e1 - ((-t467 * t521 - t468 * t477 + t469 * t478) * t491 + ((-t446 * t477 + t448 * t478) * t497 + (t477 * t445 - t478 * t447 + t540) * t499) * t514) * t512 / 0.2e1 + t465 * ((t479 * t418 - t463 * t420 + t464 * t422) * t465 + (t417 * t479 - t419 * t463 + t421 * t464) * t466 + (t439 * t479 - t440 * t463 + t441 * t464) * t482) / 0.2e1 + t466 * ((t418 * t477 - t420 * t461 + t422 * t462) * t465 + (t477 * t417 - t461 * t419 + t462 * t421) * t466 + (t439 * t477 - t440 * t461 + t441 * t462) * t482) / 0.2e1 + t482 * ((-t418 * t522 - t420 * t475 + t422 * t476) * t465 + (-t417 * t522 - t419 * t475 + t421 * t476) * t466 + (-t439 * t522 - t440 * t475 + t441 * t476) * t482) / 0.2e1 + t491 * ((t523 * t444 + (t446 * t498 + t448 * t496) * t494) * t490 - (t523 * t443 + (t445 * t498 + t447 * t496) * t494) * t512 + (t523 * t467 + (t468 * t498 + t469 * t496) * t494) * t491) / 0.2e1 + m(7) * (t358 ^ 2 + t359 ^ 2 + t360 ^ 2) / 0.2e1 + m(6) * (t361 ^ 2 + t362 ^ 2 + t363 ^ 2) / 0.2e1 + m(5) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + m(4) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(3) * (t399 ^ 2 + t400 ^ 2 + t401 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t485 ^ 2 + t486 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t532 * t436 + t533 * t437 + t531 * t463) * t457 + (t537 * t436 + t539 * t437 + t535 * t463) * t431 + (t536 * t436 + t538 * t437 + t534 * t463) * t430) * t430 / 0.2e1 + ((t532 * t434 + t533 * t435 + t531 * t461) * t457 + (t537 * t434 + t539 * t435 + t535 * t461) * t431 + (t536 * t434 + t538 * t435 + t534 * t461) * t430) * t431 / 0.2e1 + ((t532 * t459 + t533 * t460 + t531 * t475) * t457 + (t537 * t459 + t539 * t460 + t535 * t475) * t431 + (t536 * t459 + t538 * t460 + t534 * t475) * t430) * t457 / 0.2e1;
T  = t1;
