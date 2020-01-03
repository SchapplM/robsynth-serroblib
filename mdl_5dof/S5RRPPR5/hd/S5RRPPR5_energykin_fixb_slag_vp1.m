% Calculate kinetic energy for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:51
% EndTime: 2019-12-31 19:28:53
% DurationCPUTime: 2.00s
% Computational Cost: add. (902->194), mult. (1231->296), div. (0->0), fcn. (1150->8), ass. (0->115)
t436 = Icges(4,4) - Icges(5,5);
t435 = Icges(4,1) + Icges(5,1);
t434 = Icges(4,2) + Icges(5,3);
t351 = qJ(2) + pkin(8);
t349 = cos(t351);
t433 = t436 * t349;
t348 = sin(t351);
t432 = t436 * t348;
t431 = Icges(5,4) + Icges(4,5);
t430 = Icges(4,6) - Icges(5,6);
t429 = t434 * t348 - t433;
t428 = t435 * t349 - t432;
t355 = sin(qJ(1));
t358 = cos(qJ(1));
t427 = t429 * t355 + t430 * t358;
t426 = -t430 * t355 + t429 * t358;
t425 = -t428 * t355 + t431 * t358;
t424 = t431 * t355 + t428 * t358;
t423 = -t434 * t349 - t432;
t422 = t435 * t348 + t433;
t421 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t354 = sin(qJ(2));
t357 = cos(qJ(2));
t420 = Icges(3,5) * t357 - Icges(3,6) * t354 - t430 * t348 + t431 * t349;
t419 = t355 * t420 - t358 * t421;
t418 = t355 * t421 + t358 * t420;
t417 = Icges(3,5) * t354 + Icges(3,6) * t357 + t431 * t348 + t430 * t349;
t405 = Icges(3,4) * t354;
t339 = Icges(3,2) * t357 + t405;
t404 = Icges(3,4) * t357;
t340 = Icges(3,1) * t354 + t404;
t416 = -t339 * t354 + t340 * t357 + t348 * t423 + t349 * t422;
t379 = -Icges(3,2) * t354 + t404;
t315 = Icges(3,6) * t355 + t358 * t379;
t382 = Icges(3,1) * t357 - t405;
t317 = Icges(3,5) * t355 + t358 * t382;
t415 = -t315 * t354 + t317 * t357 + t348 * t426 + t349 * t424;
t314 = -Icges(3,6) * t358 + t355 * t379;
t316 = -Icges(3,5) * t358 + t355 * t382;
t414 = t314 * t354 - t316 * t357 - t348 * t427 + t425 * t349;
t410 = pkin(2) * t354;
t409 = pkin(4) * t349;
t407 = t357 * pkin(2);
t287 = -qJ(3) * t358 + t355 * t407;
t288 = qJ(3) * t355 + t358 * t407;
t395 = qJD(2) * t358;
t396 = qJD(2) * t355;
t399 = t287 * t396 + t288 * t395;
t346 = t355 * pkin(1) - t358 * pkin(6);
t398 = -t287 - t346;
t350 = qJD(3) * t355;
t394 = qJD(4) * t348;
t397 = t358 * t394 + t350;
t393 = qJD(2) - qJD(5);
t384 = pkin(3) * t349 + qJ(4) * t348;
t318 = t384 * t355;
t392 = -t318 + t398;
t389 = -t348 * pkin(3) + t349 * qJ(4) - t410;
t335 = qJD(1) * (t358 * pkin(1) + t355 * pkin(6));
t388 = qJD(1) * t288 - qJD(3) * t358 + t335;
t387 = rSges(3,1) * t357 - rSges(3,2) * t354;
t386 = rSges(4,1) * t349 - rSges(4,2) * t348;
t385 = rSges(5,1) * t349 + rSges(5,3) * t348;
t383 = qJD(2) * (-t348 * rSges(4,1) - t349 * rSges(4,2) - t410);
t353 = sin(qJ(5));
t356 = cos(qJ(5));
t323 = t348 * t356 - t349 * t353;
t364 = t348 * t353 + t349 * t356;
t363 = qJD(2) * (-t348 * rSges(5,1) + t349 * rSges(5,3) + t389);
t319 = t384 * t358;
t362 = qJD(1) * t319 + t355 * t394 + t388;
t361 = -qJD(4) * t349 + t318 * t396 + t319 * t395 + t399;
t360 = qJD(2) * (-pkin(4) * t348 + t389);
t345 = t358 * rSges(2,1) - t355 * rSges(2,2);
t344 = t355 * rSges(2,1) + t358 * rSges(2,2);
t343 = t354 * rSges(3,1) + t357 * rSges(3,2);
t337 = t393 * t358;
t336 = t393 * t355;
t331 = -t355 * pkin(7) + t358 * t409;
t330 = t358 * pkin(7) + t355 * t409;
t321 = t355 * rSges(3,3) + t358 * t387;
t320 = -t358 * rSges(3,3) + t355 * t387;
t310 = t364 * t358;
t309 = t323 * t358;
t308 = t364 * t355;
t307 = t323 * t355;
t306 = t355 * rSges(4,3) + t358 * t386;
t305 = t355 * rSges(5,2) + t358 * t385;
t304 = -t358 * rSges(4,3) + t355 * t386;
t303 = -t358 * rSges(5,2) + t355 * t385;
t283 = t323 * rSges(6,1) - rSges(6,2) * t364;
t282 = Icges(6,1) * t323 - Icges(6,4) * t364;
t281 = Icges(6,4) * t323 - Icges(6,2) * t364;
t280 = Icges(6,5) * t323 - Icges(6,6) * t364;
t279 = qJD(1) * t321 - t343 * t396 + t335;
t278 = -t343 * t395 + (-t320 - t346) * qJD(1);
t277 = (t320 * t355 + t321 * t358) * qJD(2);
t276 = t310 * rSges(6,1) + t309 * rSges(6,2) - t355 * rSges(6,3);
t275 = t308 * rSges(6,1) + t307 * rSges(6,2) + t358 * rSges(6,3);
t274 = Icges(6,1) * t310 + Icges(6,4) * t309 - Icges(6,5) * t355;
t273 = Icges(6,1) * t308 + Icges(6,4) * t307 + Icges(6,5) * t358;
t272 = Icges(6,4) * t310 + Icges(6,2) * t309 - Icges(6,6) * t355;
t271 = Icges(6,4) * t308 + Icges(6,2) * t307 + Icges(6,6) * t358;
t270 = Icges(6,5) * t310 + Icges(6,6) * t309 - Icges(6,3) * t355;
t269 = Icges(6,5) * t308 + Icges(6,6) * t307 + Icges(6,3) * t358;
t268 = qJD(1) * t306 + t355 * t383 + t388;
t267 = t350 + t358 * t383 + (-t304 + t398) * qJD(1);
t266 = (t304 * t355 + t306 * t358) * qJD(2) + t399;
t265 = qJD(1) * t305 + t355 * t363 + t362;
t264 = t358 * t363 + (-t303 + t392) * qJD(1) + t397;
t263 = (t303 * t355 + t305 * t358) * qJD(2) + t361;
t262 = -t336 * t283 + (t276 + t331) * qJD(1) + t355 * t360 + t362;
t261 = -t337 * t283 + t358 * t360 + (-t275 - t330 + t392) * qJD(1) + t397;
t260 = t336 * t275 + t337 * t276 + (t330 * t355 + t331 * t358) * qJD(2) + t361;
t1 = m(3) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(4) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(5) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t336 * ((-t355 * t270 + t309 * t272 + t310 * t274) * t336 - (-t355 * t269 + t309 * t271 + t310 * t273) * t337 + (-t355 * t280 + t309 * t281 + t310 * t282) * qJD(1)) / 0.2e1 - t337 * ((t358 * t270 + t307 * t272 + t308 * t274) * t336 - (t358 * t269 + t307 * t271 + t308 * t273) * t337 + (t358 * t280 + t307 * t281 + t308 * t282) * qJD(1)) / 0.2e1 + (m(2) * (t344 ^ 2 + t345 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t418 * t355 ^ 2 + (t414 * t358 + (t415 - t419) * t355) * t358) * qJD(2) + (t417 * t355 + t416 * t358) * qJD(1)) * t396 / 0.2e1 - ((t419 * t358 ^ 2 + (t415 * t355 + (t414 - t418) * t358) * t355) * qJD(2) + (t416 * t355 - t417 * t358) * qJD(1)) * t395 / 0.2e1 + ((-t272 * t364 + t323 * t274) * t336 - (-t271 * t364 + t323 * t273) * t337 + ((-t357 * t314 - t354 * t316 + t425 * t348 + t349 * t427) * t358 + (t357 * t315 + t354 * t317 + t348 * t424 - t349 * t426) * t355) * qJD(2) + (-t364 * t281 + t323 * t282 + t357 * t339 + t354 * t340 + t422 * t348 - t423 * t349) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
