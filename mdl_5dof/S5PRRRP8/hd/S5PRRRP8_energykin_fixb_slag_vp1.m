% Calculate kinetic energy for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:12
% EndTime: 2019-12-05 16:58:14
% DurationCPUTime: 2.08s
% Computational Cost: add. (1745->242), mult. (4468->380), div. (0->0), fcn. (5481->10), ass. (0->116)
t471 = Icges(5,1) + Icges(6,1);
t470 = -Icges(5,4) + Icges(6,5);
t469 = Icges(6,4) + Icges(5,5);
t468 = Icges(5,2) + Icges(6,3);
t467 = Icges(6,2) + Icges(5,3);
t466 = -Icges(5,6) + Icges(6,6);
t465 = rSges(6,1) + pkin(4);
t464 = rSges(6,3) + qJ(5);
t422 = sin(pkin(9));
t424 = cos(pkin(9));
t429 = cos(qJ(2));
t425 = cos(pkin(5));
t428 = sin(qJ(2));
t444 = t425 * t428;
t410 = t422 * t429 + t424 * t444;
t427 = sin(qJ(3));
t423 = sin(pkin(5));
t445 = t424 * t423;
t450 = cos(qJ(3));
t394 = t410 * t450 - t427 * t445;
t443 = t425 * t429;
t409 = t422 * t428 - t424 * t443;
t426 = sin(qJ(4));
t449 = cos(qJ(4));
t370 = t394 * t426 - t409 * t449;
t371 = t394 * t449 + t409 * t426;
t437 = t423 * t450;
t393 = t410 * t427 + t424 * t437;
t463 = t468 * t370 + t470 * t371 + t466 * t393;
t412 = -t422 * t444 + t424 * t429;
t447 = t423 * t427;
t396 = t412 * t450 + t422 * t447;
t411 = t422 * t443 + t424 * t428;
t372 = t396 * t426 - t411 * t449;
t373 = t396 * t449 + t411 * t426;
t395 = t412 * t427 - t422 * t437;
t462 = t468 * t372 + t470 * t373 + t466 * t395;
t461 = t466 * t370 + t469 * t371 + t467 * t393;
t460 = t466 * t372 + t469 * t373 + t467 * t395;
t459 = t470 * t370 + t471 * t371 + t469 * t393;
t458 = t470 * t372 + t471 * t373 + t469 * t395;
t414 = t425 * t427 + t428 * t437;
t446 = t423 * t429;
t397 = t414 * t426 + t446 * t449;
t398 = t414 * t449 - t426 * t446;
t413 = -t425 * t450 + t428 * t447;
t457 = t468 * t397 + t470 * t398 + t466 * t413;
t456 = t466 * t397 + t469 * t398 + t467 * t413;
t455 = t470 * t397 + t471 * t398 + t469 * t413;
t454 = qJD(2) ^ 2;
t448 = t422 * t423;
t442 = rSges(6,2) * t393 + t464 * t370 + t465 * t371;
t441 = rSges(6,2) * t395 + t464 * t372 + t465 * t373;
t440 = rSges(6,2) * t413 + t464 * t397 + t465 * t398;
t439 = qJD(2) * t423;
t419 = t422 * t439;
t399 = qJD(3) * t411 + t419;
t421 = qJD(2) * t425;
t436 = t424 * t439;
t389 = pkin(2) * t410 + pkin(7) * t409;
t390 = pkin(2) * t412 + pkin(7) * t411;
t435 = t389 * t419 + t390 * t436 + qJD(1);
t400 = qJD(3) * t409 - t436;
t416 = -qJD(3) * t446 + t421;
t415 = (pkin(2) * t428 - pkin(7) * t429) * t423;
t434 = t390 * t421 - t415 * t419;
t365 = pkin(3) * t394 + pkin(8) * t393;
t366 = pkin(3) * t396 + pkin(8) * t395;
t433 = t399 * t365 - t366 * t400 + t435;
t432 = (-t389 * t425 - t415 * t445) * qJD(2);
t391 = pkin(3) * t414 + pkin(8) * t413;
t431 = t416 * t366 - t391 * t399 + t434;
t430 = -t365 * t416 + t400 * t391 + t432;
t404 = t425 * rSges(3,3) + (rSges(3,1) * t428 + rSges(3,2) * t429) * t423;
t403 = Icges(3,5) * t425 + (Icges(3,1) * t428 + Icges(3,4) * t429) * t423;
t402 = Icges(3,6) * t425 + (Icges(3,4) * t428 + Icges(3,2) * t429) * t423;
t401 = Icges(3,3) * t425 + (Icges(3,5) * t428 + Icges(3,6) * t429) * t423;
t392 = qJD(4) * t413 + t416;
t387 = t414 * rSges(4,1) - t413 * rSges(4,2) - rSges(4,3) * t446;
t386 = Icges(4,1) * t414 - Icges(4,4) * t413 - Icges(4,5) * t446;
t385 = Icges(4,4) * t414 - Icges(4,2) * t413 - Icges(4,6) * t446;
t384 = Icges(4,5) * t414 - Icges(4,6) * t413 - Icges(4,3) * t446;
t381 = rSges(3,1) * t412 - rSges(3,2) * t411 + rSges(3,3) * t448;
t380 = rSges(3,1) * t410 - rSges(3,2) * t409 - rSges(3,3) * t445;
t379 = Icges(3,1) * t412 - Icges(3,4) * t411 + Icges(3,5) * t448;
t378 = Icges(3,1) * t410 - Icges(3,4) * t409 - Icges(3,5) * t445;
t377 = Icges(3,4) * t412 - Icges(3,2) * t411 + Icges(3,6) * t448;
t376 = Icges(3,4) * t410 - Icges(3,2) * t409 - Icges(3,6) * t445;
t375 = Icges(3,5) * t412 - Icges(3,6) * t411 + Icges(3,3) * t448;
t374 = Icges(3,5) * t410 - Icges(3,6) * t409 - Icges(3,3) * t445;
t369 = qJD(4) * t393 + t400;
t368 = qJD(4) * t395 + t399;
t362 = (-t380 * t425 - t404 * t445) * qJD(2);
t361 = (t381 * t425 - t404 * t448) * qJD(2);
t360 = rSges(5,1) * t398 - rSges(5,2) * t397 + rSges(5,3) * t413;
t352 = rSges(4,1) * t396 - rSges(4,2) * t395 + rSges(4,3) * t411;
t351 = rSges(4,1) * t394 - rSges(4,2) * t393 + rSges(4,3) * t409;
t350 = Icges(4,1) * t396 - Icges(4,4) * t395 + Icges(4,5) * t411;
t349 = Icges(4,1) * t394 - Icges(4,4) * t393 + Icges(4,5) * t409;
t348 = Icges(4,4) * t396 - Icges(4,2) * t395 + Icges(4,6) * t411;
t347 = Icges(4,4) * t394 - Icges(4,2) * t393 + Icges(4,6) * t409;
t346 = Icges(4,5) * t396 - Icges(4,6) * t395 + Icges(4,3) * t411;
t345 = Icges(4,5) * t394 - Icges(4,6) * t393 + Icges(4,3) * t409;
t343 = qJD(1) + (t380 * t422 + t381 * t424) * t439;
t340 = rSges(5,1) * t373 - rSges(5,2) * t372 + rSges(5,3) * t395;
t338 = rSges(5,1) * t371 - rSges(5,2) * t370 + rSges(5,3) * t393;
t324 = -t351 * t416 + t387 * t400 + t432;
t323 = t352 * t416 - t387 * t399 + t434;
t322 = t351 * t399 - t352 * t400 + t435;
t321 = -t338 * t392 + t360 * t369 + t430;
t320 = t340 * t392 - t360 * t368 + t431;
t319 = t338 * t368 - t340 * t369 + t433;
t318 = qJD(5) * t372 + t369 * t440 - t392 * t442 + t430;
t317 = qJD(5) * t370 - t368 * t440 + t392 * t441 + t431;
t316 = qJD(5) * t397 + t368 * t442 - t369 * t441 + t433;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t343 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 - t454 * ((-t375 * t445 - t377 * t409 + t379 * t410) * t448 - (-t374 * t445 - t376 * t409 + t378 * t410) * t445 + (-t401 * t445 - t402 * t409 + t403 * t410) * t425) * t445 / 0.2e1 + m(4) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + t399 * ((t411 * t346 - t395 * t348 + t396 * t350) * t399 + (t345 * t411 - t347 * t395 + t349 * t396) * t400 + (t384 * t411 - t385 * t395 + t386 * t396) * t416) / 0.2e1 + t400 * ((t346 * t409 - t348 * t393 + t350 * t394) * t399 + (t409 * t345 - t393 * t347 + t394 * t349) * t400 + (t384 * t409 - t393 * t385 + t394 * t386) * t416) / 0.2e1 + t416 * ((-t346 * t446 - t413 * t348 + t414 * t350) * t399 + (-t345 * t446 - t413 * t347 + t414 * t349) * t400 + (-t384 * t446 - t413 * t385 + t414 * t386) * t416) / 0.2e1 + m(5) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(6) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + ((t457 * t372 + t455 * t373 + t456 * t395) * t392 + (t463 * t372 + t459 * t373 + t461 * t395) * t369 + (t462 * t372 + t458 * t373 + t460 * t395) * t368) * t368 / 0.2e1 + ((t457 * t370 + t455 * t371 + t456 * t393) * t392 + (t463 * t370 + t459 * t371 + t461 * t393) * t369 + (t462 * t370 + t458 * t371 + t460 * t393) * t368) * t369 / 0.2e1 + ((t457 * t397 + t455 * t398 + t456 * t413) * t392 + (t463 * t397 + t459 * t398 + t461 * t413) * t369 + (t462 * t397 + t458 * t398 + t460 * t413) * t368) * t392 / 0.2e1 + (((t375 * t448 - t377 * t411 + t379 * t412) * t448 - (t374 * t448 - t376 * t411 + t378 * t412) * t445 + (t401 * t448 - t402 * t411 + t403 * t412) * t425) * t448 + t425 * (t425 ^ 2 * t401 + (((t377 * t429 + t379 * t428) * t422 - (t376 * t429 + t378 * t428) * t424) * t423 + (-t374 * t424 + t375 * t422 + t402 * t429 + t403 * t428) * t425) * t423)) * t454 / 0.2e1;
T = t1;
