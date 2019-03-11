% Calculate kinetic energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:27
% EndTime: 2019-03-09 01:33:27
% DurationCPUTime: 0.81s
% Computational Cost: add. (917->174), mult. (1434->284), div. (0->0), fcn. (1674->10), ass. (0->92)
t393 = cos(qJ(1));
t392 = sin(qJ(1));
t356 = cos(pkin(10));
t391 = pkin(4) * t356;
t390 = cos(pkin(9));
t389 = sin(pkin(9));
t354 = pkin(10) + qJ(5);
t351 = sin(t354);
t388 = Icges(6,4) * t351;
t352 = cos(t354);
t387 = Icges(6,4) * t352;
t339 = -t392 * t389 - t393 * t390;
t386 = t339 * t351;
t340 = t393 * t389 - t392 * t390;
t385 = t340 * t351;
t358 = sin(qJ(6));
t384 = t352 * t358;
t359 = cos(qJ(6));
t383 = t352 * t359;
t353 = qJD(2) * t392;
t381 = qJD(4) * t340 + t353;
t380 = qJD(5) * t339;
t379 = qJD(5) * t340;
t378 = qJD(6) * t351;
t342 = t392 * pkin(1) - t393 * qJ(2);
t377 = -t392 * pkin(2) - t342;
t376 = -pkin(5) * t352 - pkin(8) * t351;
t375 = -qJD(2) * t393 + qJD(1) * (t393 * pkin(1) + t392 * qJ(2));
t355 = sin(pkin(10));
t374 = rSges(5,1) * t356 - rSges(5,2) * t355;
t373 = -rSges(6,1) * t352 + rSges(6,2) * t351;
t372 = -Icges(6,1) * t352 + t388;
t371 = Icges(6,2) * t351 - t387;
t370 = -Icges(6,5) * t352 + Icges(6,6) * t351;
t308 = -Icges(6,6) * t339 + t371 * t340;
t310 = -Icges(6,5) * t339 + t372 * t340;
t369 = -t308 * t351 + t310 * t352;
t309 = Icges(6,6) * t340 + t371 * t339;
t311 = Icges(6,5) * t340 + t372 * t339;
t368 = t309 * t351 - t311 * t352;
t335 = -Icges(6,2) * t352 - t388;
t336 = -Icges(6,1) * t351 - t387;
t367 = t335 * t351 - t336 * t352;
t366 = pkin(3) * t340 + qJ(4) * t339 + t377;
t365 = qJD(1) * t393 * pkin(2) + t375;
t364 = pkin(7) * t339 + t391 * t340 + t366;
t363 = qJD(1) * (-pkin(3) * t339 + qJ(4) * t340) - qJD(4) * t339 + t365;
t362 = qJD(1) * (pkin(7) * t340 - t391 * t339) + t363;
t360 = qJD(3) ^ 2;
t345 = qJD(6) * t352 + qJD(1);
t344 = t393 * rSges(2,1) - t392 * rSges(2,2);
t343 = t392 * rSges(2,1) + t393 * rSges(2,2);
t338 = -pkin(5) * t351 + pkin(8) * t352;
t337 = -rSges(6,1) * t351 - rSges(6,2) * t352;
t334 = -Icges(6,5) * t351 - Icges(6,6) * t352;
t331 = rSges(7,3) * t352 + (-rSges(7,1) * t359 + rSges(7,2) * t358) * t351;
t330 = Icges(7,5) * t352 + (-Icges(7,1) * t359 + Icges(7,4) * t358) * t351;
t329 = Icges(7,6) * t352 + (-Icges(7,4) * t359 + Icges(7,2) * t358) * t351;
t328 = Icges(7,3) * t352 + (-Icges(7,5) * t359 + Icges(7,6) * t358) * t351;
t327 = qJD(1) * (t393 * rSges(3,1) + t392 * rSges(3,3)) + t375;
t326 = t353 + (-t392 * rSges(3,1) + t393 * rSges(3,3) - t342) * qJD(1);
t323 = -t340 * t378 - t380;
t322 = -t339 * t378 + t379;
t321 = -t339 * t383 + t340 * t358;
t320 = t339 * t384 + t340 * t359;
t319 = -t339 * t358 - t340 * t383;
t318 = -t339 * t359 + t340 * t384;
t317 = t376 * t339;
t316 = t376 * t340;
t315 = qJD(1) * (-rSges(4,1) * t339 - rSges(4,2) * t340) + t365;
t314 = t353 + (t340 * rSges(4,1) - t339 * rSges(4,2) + t377) * qJD(1);
t313 = rSges(6,3) * t340 + t373 * t339;
t312 = -rSges(6,3) * t339 + t373 * t340;
t307 = Icges(6,3) * t340 + t370 * t339;
t306 = -Icges(6,3) * t339 + t370 * t340;
t303 = rSges(7,1) * t321 + rSges(7,2) * t320 - rSges(7,3) * t386;
t302 = rSges(7,1) * t319 + rSges(7,2) * t318 - rSges(7,3) * t385;
t301 = Icges(7,1) * t321 + Icges(7,4) * t320 - Icges(7,5) * t386;
t300 = Icges(7,1) * t319 + Icges(7,4) * t318 - Icges(7,5) * t385;
t299 = Icges(7,4) * t321 + Icges(7,2) * t320 - Icges(7,6) * t386;
t298 = Icges(7,4) * t319 + Icges(7,2) * t318 - Icges(7,6) * t385;
t297 = Icges(7,5) * t321 + Icges(7,6) * t320 - Icges(7,3) * t386;
t296 = Icges(7,5) * t319 + Icges(7,6) * t318 - Icges(7,3) * t385;
t295 = qJD(1) * (rSges(5,3) * t340 - t374 * t339) + t363;
t294 = (t339 * rSges(5,3) + t374 * t340 + t366) * qJD(1) + t381;
t293 = -qJD(3) + (t312 * t340 + t313 * t339) * qJD(5);
t292 = qJD(1) * t313 - t337 * t379 + t362;
t291 = -t337 * t380 + (-t312 + t364) * qJD(1) + t381;
t290 = qJD(1) * t317 + t345 * t303 - t322 * t331 - t338 * t379 + t362;
t289 = -t338 * t380 - t345 * t302 + t323 * t331 + (-t316 + t364) * qJD(1) + t381;
t288 = t302 * t322 - t303 * t323 - qJD(3) + (t316 * t340 + t317 * t339) * qJD(5);
t1 = t345 * ((t296 * t323 + t297 * t322 + t328 * t345) * t352 + ((t299 * t358 - t301 * t359) * t322 + (t298 * t358 - t300 * t359) * t323 + (t329 * t358 - t330 * t359) * t345) * t351) / 0.2e1 + t322 * ((-t297 * t386 + t320 * t299 + t321 * t301) * t322 + (-t296 * t386 + t298 * t320 + t300 * t321) * t323 + (t320 * t329 + t321 * t330 - t328 * t386) * t345) / 0.2e1 + t323 * ((-t297 * t385 + t299 * t318 + t301 * t319) * t322 + (-t296 * t385 + t318 * t298 + t319 * t300) * t323 + (t318 * t329 + t319 * t330 - t328 * t385) * t345) / 0.2e1 + ((t340 * t334 + t367 * t339) * qJD(1) + (t340 ^ 2 * t307 + (t369 * t339 + (-t306 + t368) * t340) * t339) * qJD(5)) * t379 / 0.2e1 - ((-t339 * t334 + t367 * t340) * qJD(1) + (t339 ^ 2 * t306 + (t368 * t340 + (-t307 + t369) * t339) * t340) * qJD(5)) * t380 / 0.2e1 + qJD(1) * ((-t352 * t335 - t351 * t336) * qJD(1) + ((-t309 * t352 - t311 * t351) * t340 - (-t308 * t352 - t310 * t351) * t339) * qJD(5)) / 0.2e1 + m(7) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + m(6) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(3) * (t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(4) * (t314 ^ 2 + t315 ^ 2 + t360) / 0.2e1 + m(5) * (t294 ^ 2 + t295 ^ 2 + t360) / 0.2e1 + (Icges(3,2) + Icges(2,3) + t356 ^ 2 * Icges(5,2) + (Icges(5,1) * t355 + 0.2e1 * Icges(5,4) * t356) * t355 + Icges(4,3) + m(2) * (t343 ^ 2 + t344 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
