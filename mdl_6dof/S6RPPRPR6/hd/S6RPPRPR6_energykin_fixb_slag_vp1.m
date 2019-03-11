% Calculate kinetic energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:35
% EndTime: 2019-03-09 01:50:37
% DurationCPUTime: 1.70s
% Computational Cost: add. (439->186), mult. (1017->285), div. (0->0), fcn. (924->6), ass. (0->97)
t427 = -Icges(5,4) - Icges(6,6);
t426 = Icges(5,1) + Icges(6,2);
t425 = Icges(5,2) + Icges(6,3);
t360 = cos(qJ(4));
t424 = t427 * t360;
t357 = sin(qJ(4));
t423 = t427 * t357;
t422 = -Icges(6,4) + Icges(5,5);
t421 = Icges(6,5) - Icges(5,6);
t420 = t360 * t425 - t423;
t419 = -t357 * t426 + t424;
t418 = Icges(6,1) + Icges(5,3);
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t417 = t420 * t358 - t421 * t361;
t416 = t421 * t358 + t420 * t361;
t415 = -t419 * t358 + t422 * t361;
t414 = t422 * t358 + t419 * t361;
t413 = t357 * t425 + t424;
t412 = t360 * t426 + t423;
t411 = t422 * t357 - t421 * t360;
t410 = t411 * t358 + t418 * t361;
t409 = -t418 * t358 + t411 * t361;
t408 = t421 * t357 + t422 * t360;
t407 = t412 * t357 - t413 * t360;
t406 = t415 * t357 + t417 * t360;
t405 = t414 * t357 - t416 * t360;
t400 = pkin(7) * qJD(1);
t395 = t357 * t358;
t394 = t357 * t361;
t393 = t358 * t360;
t392 = t360 * t361;
t355 = qJD(2) * t358;
t391 = qJD(3) * t361 + t355;
t390 = qJD(4) * t358;
t389 = qJD(4) * t361;
t388 = qJD(6) * t357;
t346 = pkin(4) * t360 + qJ(5) * t357;
t387 = t346 * t389 + t391;
t384 = -qJD(2) * t361 + qJD(1) * (pkin(1) * t361 + qJ(2) * t358);
t383 = rSges(5,1) * t357 + rSges(5,2) * t360;
t382 = -rSges(6,2) * t357 - rSges(6,3) * t360;
t381 = pkin(4) * t357 - qJ(5) * t360;
t380 = (pkin(8) * qJD(4) - qJD(5)) * t360;
t367 = qJD(4) * (-rSges(6,2) * t360 + rSges(6,3) * t357) - qJD(5) * t360;
t366 = qJD(1) * t361 * qJ(3) + qJD(3) * t358 + t384;
t344 = pkin(1) * t358 - qJ(2) * t361;
t365 = -pkin(7) * t361 - qJ(3) * t358 - t344;
t329 = t381 * t358;
t364 = -t329 + t365;
t330 = t381 * t361;
t363 = qJD(1) * t330 + t346 * t390 + t366;
t359 = cos(qJ(6));
t356 = sin(qJ(6));
t352 = qJD(5) * t357;
t350 = qJD(6) * t360 + qJD(1);
t349 = rSges(2,1) * t361 - rSges(2,2) * t358;
t348 = rSges(5,1) * t360 - rSges(5,2) * t357;
t345 = rSges(2,1) * t358 + rSges(2,2) * t361;
t336 = pkin(5) * t361 + pkin(8) * t395;
t335 = -pkin(5) * t358 + pkin(8) * t394;
t334 = t358 * t388 + t389;
t333 = t361 * t388 - t390;
t328 = -t356 * t393 + t359 * t361;
t327 = -t356 * t361 - t359 * t393;
t326 = -t356 * t392 - t358 * t359;
t325 = t356 * t358 - t359 * t392;
t323 = rSges(6,1) * t361 + t382 * t358;
t322 = -rSges(6,1) * t358 + t382 * t361;
t321 = -rSges(5,3) * t358 + t383 * t361;
t320 = rSges(5,3) * t361 + t383 * t358;
t319 = rSges(7,3) * t360 + (rSges(7,1) * t356 + rSges(7,2) * t359) * t357;
t310 = Icges(7,5) * t360 + (Icges(7,1) * t356 + Icges(7,4) * t359) * t357;
t307 = Icges(7,6) * t360 + (Icges(7,4) * t356 + Icges(7,2) * t359) * t357;
t304 = Icges(7,3) * t360 + (Icges(7,5) * t356 + Icges(7,6) * t359) * t357;
t303 = qJD(1) * (-rSges(3,2) * t361 + rSges(3,3) * t358) + t384;
t302 = t355 + (rSges(3,2) * t358 + rSges(3,3) * t361 - t344) * qJD(1);
t301 = qJD(1) * (rSges(4,2) * t358 + rSges(4,3) * t361) + t366;
t300 = (t361 * rSges(4,2) - t344 + (-rSges(4,3) - qJ(3)) * t358) * qJD(1) + t391;
t299 = rSges(7,1) * t328 + rSges(7,2) * t327 + rSges(7,3) * t395;
t298 = rSges(7,1) * t326 + rSges(7,2) * t325 + rSges(7,3) * t394;
t297 = Icges(7,1) * t328 + Icges(7,4) * t327 + Icges(7,5) * t395;
t296 = Icges(7,1) * t326 + Icges(7,4) * t325 + Icges(7,5) * t394;
t295 = Icges(7,4) * t328 + Icges(7,2) * t327 + Icges(7,6) * t395;
t294 = Icges(7,4) * t326 + Icges(7,2) * t325 + Icges(7,6) * t394;
t293 = Icges(7,5) * t328 + Icges(7,6) * t327 + Icges(7,3) * t395;
t292 = Icges(7,5) * t326 + Icges(7,6) * t325 + Icges(7,3) * t394;
t291 = (-t320 * t358 - t321 * t361) * qJD(4);
t290 = t348 * t390 + (-pkin(7) * t358 + t321) * qJD(1) + t366;
t289 = t348 * t389 + (-t320 + t365) * qJD(1) + t391;
t288 = t352 + ((-t322 - t330) * t361 + (-t323 - t329) * t358) * qJD(4);
t287 = qJD(1) * t322 + (t367 - t400) * t358 + t363;
t286 = t367 * t361 + (-t323 + t364) * qJD(1) + t387;
t285 = qJD(1) * t335 + t298 * t350 - t319 * t333 + (t380 - t400) * t358 + t363;
t284 = -t299 * t350 + t319 * t334 + t361 * t380 + (-t336 + t364) * qJD(1) + t387;
t283 = -t298 * t334 + t299 * t333 + t352 + ((-t330 - t335) * t361 + (-t329 - t336) * t358) * qJD(4);
t1 = t333 * ((t292 * t394 + t294 * t325 + t296 * t326) * t333 + (t293 * t394 + t295 * t325 + t297 * t326) * t334 + (t304 * t394 + t307 * t325 + t310 * t326) * t350) / 0.2e1 + t334 * ((t292 * t395 + t294 * t327 + t296 * t328) * t333 + (t293 * t395 + t295 * t327 + t297 * t328) * t334 + (t304 * t395 + t307 * t327 + t310 * t328) * t350) / 0.2e1 + t350 * ((t292 * t333 + t293 * t334 + t304 * t350) * t360 + ((t294 * t359 + t296 * t356) * t333 + (t295 * t359 + t297 * t356) * t334 + (t307 * t359 + t310 * t356) * t350) * t357) / 0.2e1 + m(6) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(7) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(5) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(3) * (t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(4) * (t300 ^ 2 + t301 ^ 2) / 0.2e1 + (((-t417 * t357 + t415 * t360) * t361 + (t416 * t357 + t414 * t360) * t358) * qJD(4) + (t413 * t357 + t412 * t360) * qJD(1)) * qJD(1) / 0.2e1 - ((t409 * t358 ^ 2 + (t406 * t361 + (t405 - t410) * t358) * t361) * qJD(4) + (-t408 * t358 + t407 * t361) * qJD(1)) * t390 / 0.2e1 + ((t410 * t361 ^ 2 + (t405 * t358 + (t406 - t409) * t361) * t358) * qJD(4) + (t407 * t358 + t408 * t361) * qJD(1)) * t389 / 0.2e1 + (Icges(2,3) + Icges(3,1) + Icges(4,1) + m(2) * (t345 ^ 2 + t349 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
