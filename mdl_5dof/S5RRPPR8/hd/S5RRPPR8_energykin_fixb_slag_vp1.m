% Calculate kinetic energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:05
% EndTime: 2019-12-31 19:38:07
% DurationCPUTime: 2.25s
% Computational Cost: add. (769->209), mult. (1460->317), div. (0->0), fcn. (1457->8), ass. (0->119)
t442 = Icges(3,4) - Icges(4,5);
t441 = Icges(3,1) + Icges(4,1);
t440 = Icges(3,2) + Icges(4,3);
t366 = cos(qJ(2));
t439 = t442 * t366;
t364 = sin(qJ(2));
t438 = t442 * t364;
t437 = Icges(4,4) + Icges(3,5);
t436 = Icges(3,6) - Icges(4,6);
t435 = t440 * t364 - t439;
t434 = t441 * t366 - t438;
t433 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t365 = sin(qJ(1));
t367 = cos(qJ(1));
t431 = t435 * t365 + t436 * t367;
t430 = -t436 * t365 + t435 * t367;
t429 = -t434 * t365 + t437 * t367;
t428 = t437 * t365 + t434 * t367;
t427 = -t440 * t366 - t438;
t426 = t441 * t364 + t439;
t425 = -t436 * t364 + t437 * t366;
t362 = cos(pkin(8));
t361 = sin(pkin(8));
t405 = t366 * t361;
t336 = t364 * t362 - t405;
t323 = t336 * t365;
t406 = t364 * t361;
t373 = t366 * t362 + t406;
t324 = t373 * t365;
t424 = -Icges(5,5) * t324 - Icges(5,6) * t323 + t425 * t365 - t433 * t367;
t325 = t336 * t367;
t326 = t373 * t367;
t423 = -Icges(5,5) * t326 - Icges(5,6) * t325 + t433 * t365 + t425 * t367;
t420 = t364 * t427 + t366 * t426;
t419 = t364 * t430 + t366 * t428;
t418 = -t364 * t431 + t366 * t429;
t417 = Icges(5,5) * t336 - Icges(5,6) * t373 - t437 * t364 - t436 * t366;
t413 = pkin(3) * t366;
t411 = t362 * pkin(4);
t388 = pkin(2) * t366 + qJ(3) * t364;
t332 = t388 * t365;
t353 = t365 * pkin(1) - t367 * pkin(6);
t403 = -t332 - t353;
t402 = qJD(2) * t365;
t401 = qJD(2) * t367;
t400 = qJD(3) * t364;
t399 = qJD(2) - qJD(5);
t333 = t388 * t367;
t339 = qJD(1) * (t367 * pkin(1) + t365 * pkin(6));
t398 = qJD(1) * t333 + t365 * t400 + t339;
t337 = t367 * qJ(4) + t365 * t413;
t397 = -t337 + t403;
t348 = t364 * pkin(2) - t366 * qJ(3);
t394 = -t364 * pkin(3) - t348;
t393 = qJD(2) * (-t364 * rSges(4,1) + t366 * rSges(4,3) - t348);
t355 = t367 * t400;
t392 = -qJD(4) * t365 + t355;
t391 = -qJD(3) * t366 + t332 * t402 + t333 * t401;
t390 = rSges(3,1) * t366 - rSges(3,2) * t364;
t389 = rSges(4,1) * t366 + rSges(4,3) * t364;
t338 = -t365 * qJ(4) + t367 * t413;
t387 = qJD(1) * t338 + qJD(4) * t367 + t398;
t360 = pkin(8) + qJ(5);
t357 = sin(t360);
t358 = cos(t360);
t328 = -t366 * t357 + t364 * t358;
t374 = t364 * t357 + t366 * t358;
t372 = qJD(2) * (-t336 * rSges(5,1) + rSges(5,2) * t373 + t394);
t371 = qJD(2) * (pkin(4) * t405 - t411 * t364 + t394);
t370 = t337 * t402 + t338 * t401 + t391;
t369 = pkin(4) * t406 + t411 * t366;
t352 = t367 * rSges(2,1) - t365 * rSges(2,2);
t351 = t365 * rSges(2,1) + t367 * rSges(2,2);
t350 = t364 * rSges(3,1) + t366 * rSges(3,2);
t341 = t399 * t367;
t340 = t399 * t365;
t322 = t365 * rSges(3,3) + t390 * t367;
t321 = t365 * rSges(4,2) + t389 * t367;
t320 = -t367 * rSges(3,3) + t390 * t365;
t319 = -t367 * rSges(4,2) + t389 * t365;
t303 = t374 * t367;
t302 = t328 * t367;
t301 = t374 * t365;
t300 = t328 * t365;
t298 = Icges(5,1) * t336 - Icges(5,4) * t373;
t297 = Icges(5,4) * t336 - Icges(5,2) * t373;
t295 = t328 * rSges(6,1) - rSges(6,2) * t374;
t294 = Icges(6,1) * t328 - Icges(6,4) * t374;
t293 = Icges(6,4) * t328 - Icges(6,2) * t374;
t292 = Icges(6,5) * t328 - Icges(6,6) * t374;
t291 = -pkin(7) * t365 + t369 * t367;
t290 = pkin(7) * t367 + t369 * t365;
t289 = t326 * rSges(5,1) + t325 * rSges(5,2) - t365 * rSges(5,3);
t288 = t324 * rSges(5,1) + t323 * rSges(5,2) + t367 * rSges(5,3);
t287 = Icges(5,1) * t326 + Icges(5,4) * t325 - Icges(5,5) * t365;
t286 = Icges(5,1) * t324 + Icges(5,4) * t323 + Icges(5,5) * t367;
t285 = Icges(5,4) * t326 + Icges(5,2) * t325 - Icges(5,6) * t365;
t284 = Icges(5,4) * t324 + Icges(5,2) * t323 + Icges(5,6) * t367;
t281 = qJD(1) * t322 - t350 * t402 + t339;
t280 = -t350 * t401 + (-t320 - t353) * qJD(1);
t279 = (t320 * t365 + t322 * t367) * qJD(2);
t278 = t303 * rSges(6,1) + t302 * rSges(6,2) - t365 * rSges(6,3);
t277 = t301 * rSges(6,1) + t300 * rSges(6,2) + t367 * rSges(6,3);
t276 = Icges(6,1) * t303 + Icges(6,4) * t302 - Icges(6,5) * t365;
t275 = Icges(6,1) * t301 + Icges(6,4) * t300 + Icges(6,5) * t367;
t274 = Icges(6,4) * t303 + Icges(6,2) * t302 - Icges(6,6) * t365;
t273 = Icges(6,4) * t301 + Icges(6,2) * t300 + Icges(6,6) * t367;
t272 = Icges(6,5) * t303 + Icges(6,6) * t302 - Icges(6,3) * t365;
t271 = Icges(6,5) * t301 + Icges(6,6) * t300 + Icges(6,3) * t367;
t270 = qJD(1) * t321 + t365 * t393 + t398;
t269 = t355 + t367 * t393 + (-t319 + t403) * qJD(1);
t268 = (t319 * t365 + t321 * t367) * qJD(2) + t391;
t267 = qJD(1) * t289 + t365 * t372 + t387;
t266 = t367 * t372 + (-t288 + t397) * qJD(1) + t392;
t265 = (t288 * t365 + t289 * t367) * qJD(2) + t370;
t264 = -t340 * t295 + (t278 + t291) * qJD(1) + t365 * t371 + t387;
t263 = -t341 * t295 + t367 * t371 + (-t277 - t290 + t397) * qJD(1) + t392;
t262 = t340 * t277 + t341 * t278 + (t290 * t365 + t291 * t367) * qJD(2) + t370;
t1 = m(3) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(4) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t340 * ((-t365 * t272 + t302 * t274 + t303 * t276) * t340 - (-t365 * t271 + t302 * t273 + t303 * t275) * t341 + (-t365 * t292 + t302 * t293 + t303 * t294) * qJD(1)) / 0.2e1 - t341 * ((t367 * t272 + t300 * t274 + t301 * t276) * t340 - (t367 * t271 + t300 * t273 + t301 * t275) * t341 + (t367 * t292 + t300 * t293 + t301 * t294) * qJD(1)) / 0.2e1 + (m(2) * (t351 ^ 2 + t352 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t325 * t284 - t326 * t286 + t418 * t367) * t367 + (t325 * t285 + t326 * t287 + (t419 - t424) * t367 + t423 * t365) * t365) * qJD(2) + (t325 * t297 + t326 * t298 - t417 * t365 + t420 * t367) * qJD(1)) * t402 / 0.2e1 - (((t323 * t285 + t324 * t287 + t419 * t365) * t365 + (-t323 * t284 - t324 * t286 + (t418 - t423) * t365 + t424 * t367) * t367) * qJD(2) + (t323 * t297 + t324 * t298 + t420 * t365 + t417 * t367) * qJD(1)) * t401 / 0.2e1 + ((-t274 * t374 + t328 * t276) * t340 - (-t273 * t374 + t328 * t275) * t341 + ((t373 * t284 - t336 * t286 + t364 * t429 + t366 * t431) * t367 + (-t373 * t285 + t336 * t287 + t364 * t428 - t366 * t430) * t365) * qJD(2) + (-t374 * t293 + t328 * t294 - t373 * t297 + t336 * t298 + t426 * t364 - t427 * t366) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
