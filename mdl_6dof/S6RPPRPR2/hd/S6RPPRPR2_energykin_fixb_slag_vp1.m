% Calculate kinetic energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:32
% EndTime: 2019-03-09 01:41:34
% DurationCPUTime: 1.99s
% Computational Cost: add. (1245->190), mult. (1043->295), div. (0->0), fcn. (948->10), ass. (0->109)
t472 = Icges(5,4) + Icges(6,6);
t471 = Icges(5,1) + Icges(6,2);
t470 = -Icges(5,2) - Icges(6,3);
t390 = pkin(10) + qJ(4);
t388 = cos(t390);
t469 = t472 * t388;
t386 = sin(t390);
t468 = t472 * t386;
t467 = -Icges(6,4) + Icges(5,5);
t466 = Icges(6,5) - Icges(5,6);
t465 = t470 * t386 + t469;
t464 = -t471 * t388 + t468;
t463 = Icges(6,1) + Icges(5,3);
t391 = qJ(1) + pkin(9);
t387 = sin(t391);
t389 = cos(t391);
t462 = t465 * t387 + t466 * t389;
t461 = -t466 * t387 + t465 * t389;
t460 = t464 * t387 + t467 * t389;
t459 = t467 * t387 - t464 * t389;
t458 = t470 * t388 - t468;
t457 = t471 * t386 + t469;
t456 = t466 * t386 + t467 * t388;
t455 = t463 * t387 + t456 * t389;
t454 = t456 * t387 - t463 * t389;
t453 = t467 * t386 - t466 * t388;
t452 = t458 * t386 + t457 * t388;
t451 = -t461 * t386 + t459 * t388;
t450 = t462 * t386 + t460 * t388;
t396 = sin(qJ(1));
t446 = t396 * pkin(1);
t393 = cos(pkin(10));
t444 = pkin(3) * t393;
t395 = sin(qJ(6));
t439 = t387 * t395;
t397 = cos(qJ(6));
t438 = t387 * t397;
t437 = t388 * t387;
t436 = t388 * t389;
t435 = t389 * t395;
t434 = t389 * t397;
t398 = cos(qJ(1));
t385 = qJD(1) * t398 * pkin(1);
t432 = qJD(1) * (pkin(2) * t389 + qJ(3) * t387) + t385;
t383 = qJD(3) * t387;
t428 = qJD(5) * t386;
t431 = t389 * t428 + t383;
t430 = qJD(4) * t387;
t429 = qJD(4) * t389;
t427 = qJD(6) * t388;
t424 = -pkin(2) * t387 + qJ(3) * t389 - t446;
t374 = pkin(4) * t386 - qJ(5) * t388;
t423 = qJD(4) * (rSges(6,2) * t386 + rSges(6,3) * t388 - t374);
t422 = pkin(7) * t389 - t444 * t387 + t424;
t392 = sin(pkin(10));
t421 = rSges(4,1) * t393 - rSges(4,2) * t392;
t420 = rSges(5,1) * t388 - rSges(5,2) * t386;
t419 = -rSges(6,2) * t388 + rSges(6,3) * t386;
t418 = pkin(4) * t388 + qJ(5) * t386;
t417 = qJD(4) * (-pkin(8) * t386 - t374);
t355 = t418 * t387;
t404 = -t355 + t422;
t403 = -qJD(3) * t389 + qJD(1) * (pkin(7) * t387 + t444 * t389) + t432;
t356 = t418 * t389;
t402 = -qJD(5) * t388 + t355 * t430 + t356 * t429 + qJD(2);
t401 = qJD(1) * t356 + t387 * t428 + t403;
t399 = qJD(2) ^ 2;
t382 = qJD(6) * t386 + qJD(1);
t381 = rSges(2,1) * t398 - rSges(2,2) * t396;
t380 = rSges(2,1) * t396 + rSges(2,2) * t398;
t376 = rSges(5,1) * t386 + rSges(5,2) * t388;
t366 = -pkin(5) * t389 + pkin(8) * t437;
t365 = pkin(5) * t387 + pkin(8) * t436;
t364 = t387 * t427 - t429;
t363 = t389 * t427 + t430;
t362 = t386 * t439 - t434;
t361 = t386 * t438 + t435;
t360 = t386 * t435 + t438;
t359 = t386 * t434 - t439;
t358 = t385 + qJD(1) * (rSges(3,1) * t389 - rSges(3,2) * t387);
t357 = (-rSges(3,1) * t387 - rSges(3,2) * t389 - t446) * qJD(1);
t353 = rSges(7,3) * t386 + (-rSges(7,1) * t395 - rSges(7,2) * t397) * t388;
t352 = Icges(7,5) * t386 + (-Icges(7,1) * t395 - Icges(7,4) * t397) * t388;
t351 = Icges(7,6) * t386 + (-Icges(7,4) * t395 - Icges(7,2) * t397) * t388;
t350 = Icges(7,3) * t386 + (-Icges(7,5) * t395 - Icges(7,6) * t397) * t388;
t349 = -rSges(6,1) * t389 + t419 * t387;
t348 = rSges(6,1) * t387 + t419 * t389;
t347 = rSges(5,3) * t387 + t420 * t389;
t346 = -rSges(5,3) * t389 + t420 * t387;
t329 = qJD(1) * t387 * rSges(4,3) + (qJD(1) * t421 - qJD(3)) * t389 + t432;
t328 = t383 + (t389 * rSges(4,3) - t421 * t387 + t424) * qJD(1);
t327 = rSges(7,1) * t362 + rSges(7,2) * t361 + rSges(7,3) * t437;
t326 = rSges(7,1) * t360 + rSges(7,2) * t359 + rSges(7,3) * t436;
t325 = Icges(7,1) * t362 + Icges(7,4) * t361 + Icges(7,5) * t437;
t324 = Icges(7,1) * t360 + Icges(7,4) * t359 + Icges(7,5) * t436;
t323 = Icges(7,4) * t362 + Icges(7,2) * t361 + Icges(7,6) * t437;
t322 = Icges(7,4) * t360 + Icges(7,2) * t359 + Icges(7,6) * t436;
t321 = Icges(7,5) * t362 + Icges(7,6) * t361 + Icges(7,3) * t437;
t320 = Icges(7,5) * t360 + Icges(7,6) * t359 + Icges(7,3) * t436;
t319 = qJD(2) + (t346 * t387 + t347 * t389) * qJD(4);
t318 = qJD(1) * t347 - t376 * t430 + t403;
t317 = -t376 * t429 + t383 + (-t346 + t422) * qJD(1);
t316 = (t348 * t389 + t349 * t387) * qJD(4) + t402;
t315 = qJD(1) * t348 + t387 * t423 + t401;
t314 = t389 * t423 + (-t349 + t404) * qJD(1) + t431;
t313 = qJD(1) * t365 + t326 * t382 - t353 * t363 + t387 * t417 + t401;
t312 = -t327 * t382 + t353 * t364 + t389 * t417 + (-t366 + t404) * qJD(1) + t431;
t311 = -t326 * t364 + t327 * t363 + (t365 * t389 + t366 * t387) * qJD(4) + t402;
t1 = m(4) * (t328 ^ 2 + t329 ^ 2 + t399) / 0.2e1 + t363 * ((t320 * t436 + t359 * t322 + t360 * t324) * t363 + (t321 * t436 + t323 * t359 + t325 * t360) * t364 + (t350 * t436 + t351 * t359 + t352 * t360) * t382) / 0.2e1 + t364 * ((t320 * t437 + t322 * t361 + t324 * t362) * t363 + (t321 * t437 + t361 * t323 + t362 * t325) * t364 + (t350 * t437 + t351 * t361 + t352 * t362) * t382) / 0.2e1 + t382 * ((t320 * t363 + t321 * t364 + t350 * t382) * t386 + ((-t322 * t397 - t324 * t395) * t363 + (-t323 * t397 - t325 * t395) * t364 + (-t351 * t397 - t352 * t395) * t382) * t388) / 0.2e1 + m(7) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(6) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(5) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + m(3) * (t357 ^ 2 + t358 ^ 2 + t399) / 0.2e1 + (((t460 * t386 - t462 * t388) * t389 + (t459 * t386 + t461 * t388) * t387) * qJD(4) + (t457 * t386 - t458 * t388) * qJD(1)) * qJD(1) / 0.2e1 + ((t455 * t387 ^ 2 + (t450 * t389 + (t451 - t454) * t387) * t389) * qJD(4) + (t453 * t387 + t452 * t389) * qJD(1)) * t430 / 0.2e1 - ((t454 * t389 ^ 2 + (t451 * t387 + (t450 - t455) * t389) * t387) * qJD(4) + (t452 * t387 - t453 * t389) * qJD(1)) * t429 / 0.2e1 + (m(2) * (t380 ^ 2 + t381 ^ 2) + Icges(3,3) + Icges(4,2) * t393 ^ 2 + (Icges(4,1) * t392 + 0.2e1 * Icges(4,4) * t393) * t392 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
