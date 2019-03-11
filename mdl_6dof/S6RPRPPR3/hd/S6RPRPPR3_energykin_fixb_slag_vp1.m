% Calculate kinetic energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:05
% EndTime: 2019-03-09 02:44:07
% DurationCPUTime: 2.06s
% Computational Cost: add. (1070->193), mult. (1260->302), div. (0->0), fcn. (1137->8), ass. (0->111)
t479 = Icges(4,4) + Icges(6,4) - Icges(5,5);
t478 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t477 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t394 = sin(qJ(3));
t476 = t479 * t394;
t397 = cos(qJ(3));
t475 = t479 * t397;
t474 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t473 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t472 = t477 * t394 - t475;
t471 = -t478 * t397 + t476;
t470 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t392 = qJ(1) + pkin(9);
t390 = sin(t392);
t391 = cos(t392);
t469 = t390 * t472 + t391 * t473;
t468 = -t390 * t473 + t391 * t472;
t467 = t471 * t390 + t391 * t474;
t466 = t390 * t474 - t471 * t391;
t465 = -t477 * t397 - t476;
t464 = t478 * t394 + t475;
t463 = -t473 * t394 + t397 * t474;
t462 = t470 * t390 + t463 * t391;
t461 = t463 * t390 - t470 * t391;
t460 = -t394 * t474 - t473 * t397;
t459 = t465 * t394 + t464 * t397;
t458 = -t469 * t394 + t467 * t397;
t457 = t468 * t394 + t466 * t397;
t395 = sin(qJ(1));
t453 = pkin(1) * t395;
t445 = t390 * t397;
t444 = t391 * t397;
t393 = sin(qJ(6));
t443 = t393 * t394;
t396 = cos(qJ(6));
t442 = t394 * t396;
t398 = cos(qJ(1));
t389 = qJD(1) * t398 * pkin(1);
t441 = qJD(1) * (pkin(2) * t391 + pkin(7) * t390) + t389;
t440 = qJD(3) * t390;
t439 = qJD(3) * t391;
t438 = qJD(4) * t394;
t437 = qJD(6) * t397;
t434 = -pkin(2) * t390 + pkin(7) * t391 - t453;
t380 = pkin(3) * t394 - qJ(4) * t397;
t433 = -pkin(4) * t394 - t380;
t432 = qJD(3) * (-rSges(5,1) * t394 + rSges(5,3) * t397 - t380);
t379 = t391 * t438;
t431 = -qJD(5) * t390 + t379;
t424 = pkin(3) * t397 + qJ(4) * t394;
t358 = t424 * t391;
t430 = qJD(1) * t358 + t390 * t438 + t441;
t357 = t424 * t390;
t429 = -t357 + t434;
t428 = pkin(5) * t394 + pkin(8) * t397;
t427 = rSges(4,1) * t397 - rSges(4,2) * t394;
t426 = rSges(5,1) * t397 + rSges(5,3) * t394;
t425 = rSges(6,1) * t394 - rSges(6,2) * t397;
t365 = pkin(4) * t445 + qJ(5) * t391;
t405 = -t365 + t429;
t404 = -qJD(4) * t397 + t357 * t440 + t358 * t439 + qJD(2);
t403 = qJD(3) * (rSges(6,1) * t397 + rSges(6,2) * t394 + t433);
t402 = qJD(3) * (pkin(5) * t397 - pkin(8) * t394 + t433);
t366 = pkin(4) * t444 - qJ(5) * t390;
t401 = qJD(1) * t366 + qJD(5) * t391 + t430;
t400 = t365 * t440 + t366 * t439 + t404;
t388 = qJD(6) * t394 + qJD(1);
t385 = rSges(2,1) * t398 - rSges(2,2) * t395;
t383 = rSges(2,1) * t395 + rSges(2,2) * t398;
t382 = rSges(4,1) * t394 + rSges(4,2) * t397;
t364 = t390 * t437 - t439;
t363 = t391 * t437 + t440;
t361 = t428 * t391;
t360 = t428 * t390;
t359 = rSges(7,3) * t394 + (-rSges(7,1) * t396 + rSges(7,2) * t393) * t397;
t356 = Icges(7,5) * t394 + (-Icges(7,1) * t396 + Icges(7,4) * t393) * t397;
t355 = Icges(7,6) * t394 + (-Icges(7,4) * t396 + Icges(7,2) * t393) * t397;
t354 = Icges(7,3) * t394 + (-Icges(7,5) * t396 + Icges(7,6) * t393) * t397;
t353 = -t390 * t393 + t391 * t442;
t352 = -t390 * t396 - t391 * t443;
t351 = t390 * t442 + t391 * t393;
t350 = -t390 * t443 + t391 * t396;
t348 = t389 + qJD(1) * (rSges(3,1) * t391 - rSges(3,2) * t390);
t347 = (-rSges(3,1) * t390 - rSges(3,2) * t391 - t453) * qJD(1);
t344 = rSges(4,3) * t390 + t391 * t427;
t343 = rSges(5,2) * t390 + t391 * t426;
t342 = -rSges(6,3) * t390 + t391 * t425;
t341 = -rSges(4,3) * t391 + t390 * t427;
t340 = -rSges(5,2) * t391 + t390 * t426;
t339 = rSges(6,3) * t391 + t390 * t425;
t318 = rSges(7,1) * t353 + rSges(7,2) * t352 + rSges(7,3) * t444;
t317 = rSges(7,1) * t351 + rSges(7,2) * t350 + rSges(7,3) * t445;
t316 = Icges(7,1) * t353 + Icges(7,4) * t352 + Icges(7,5) * t444;
t315 = Icges(7,1) * t351 + Icges(7,4) * t350 + Icges(7,5) * t445;
t314 = Icges(7,4) * t353 + Icges(7,2) * t352 + Icges(7,6) * t444;
t313 = Icges(7,4) * t351 + Icges(7,2) * t350 + Icges(7,6) * t445;
t312 = Icges(7,5) * t353 + Icges(7,6) * t352 + Icges(7,3) * t444;
t311 = Icges(7,5) * t351 + Icges(7,6) * t350 + Icges(7,3) * t445;
t310 = qJD(1) * t344 - t382 * t440 + t441;
t309 = -t382 * t439 + (-t341 + t434) * qJD(1);
t308 = qJD(2) + (t341 * t390 + t344 * t391) * qJD(3);
t307 = qJD(1) * t343 + t390 * t432 + t430;
t306 = t379 + t391 * t432 + (-t340 + t429) * qJD(1);
t305 = (t340 * t390 + t343 * t391) * qJD(3) + t404;
t304 = qJD(1) * t342 + t390 * t403 + t401;
t303 = t391 * t403 + (-t339 + t405) * qJD(1) + t431;
t302 = (t339 * t390 + t342 * t391) * qJD(3) + t400;
t301 = qJD(1) * t361 + t318 * t388 - t359 * t363 + t390 * t402 + t401;
t300 = -t317 * t388 + t359 * t364 + t391 * t402 + (-t360 + t405) * qJD(1) + t431;
t299 = t317 * t363 - t318 * t364 + (t360 * t390 + t361 * t391) * qJD(3) + t400;
t1 = t363 * ((t312 * t444 + t352 * t314 + t353 * t316) * t363 + (t311 * t444 + t313 * t352 + t315 * t353) * t364 + (t352 * t355 + t353 * t356 + t354 * t444) * t388) / 0.2e1 + t388 * ((t311 * t364 + t312 * t363 + t354 * t388) * t394 + ((t314 * t393 - t316 * t396) * t363 + (t313 * t393 - t315 * t396) * t364 + (t355 * t393 - t356 * t396) * t388) * t397) / 0.2e1 + t364 * ((t312 * t445 + t314 * t350 + t316 * t351) * t363 + (t311 * t445 + t350 * t313 + t351 * t315) * t364 + (t350 * t355 + t351 * t356 + t354 * t445) * t388) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(7) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + m(4) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t383 ^ 2 + t385 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t467 * t394 + t469 * t397) * t391 + (t466 * t394 - t468 * t397) * t390) * qJD(3) + (t464 * t394 - t465 * t397) * qJD(1)) * qJD(1) / 0.2e1 + ((t462 * t390 ^ 2 + (t458 * t391 + (t457 - t461) * t390) * t391) * qJD(3) + (-t460 * t390 + t459 * t391) * qJD(1)) * t440 / 0.2e1 - ((t461 * t391 ^ 2 + (t457 * t390 + (t458 - t462) * t391) * t390) * qJD(3) + (t459 * t390 + t460 * t391) * qJD(1)) * t439 / 0.2e1;
T  = t1;
