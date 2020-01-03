% Calculate kinetic energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:39
% EndTime: 2019-12-31 19:40:41
% DurationCPUTime: 1.99s
% Computational Cost: add. (514->180), mult. (1226->287), div. (0->0), fcn. (1123->6), ass. (0->106)
t433 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t432 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t431 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t349 = sin(qJ(2));
t430 = t433 * t349;
t352 = cos(qJ(2));
t429 = t433 * t352;
t428 = Icges(4,4) + Icges(3,5) + Icges(5,6);
t427 = Icges(5,5) + Icges(3,6) - Icges(4,6);
t426 = t431 * t349 - t429;
t425 = -t432 * t352 + t430;
t424 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t350 = sin(qJ(1));
t353 = cos(qJ(1));
t423 = t426 * t350 + t427 * t353;
t422 = -t427 * t350 + t426 * t353;
t421 = t425 * t350 + t428 * t353;
t420 = t428 * t350 - t425 * t353;
t419 = -t431 * t352 - t430;
t418 = t432 * t349 + t429;
t417 = -t427 * t349 + t428 * t352;
t416 = t417 * t350 - t424 * t353;
t415 = t424 * t350 + t417 * t353;
t414 = -t428 * t349 - t427 * t352;
t413 = t419 * t349 + t418 * t352;
t412 = t422 * t349 + t420 * t352;
t411 = -t423 * t349 + t421 * t352;
t348 = sin(qJ(5));
t400 = t348 * t350;
t399 = t348 * t353;
t351 = cos(qJ(5));
t398 = t350 * t351;
t397 = t350 * t352;
t396 = t351 * t353;
t395 = t352 * t353;
t377 = pkin(2) * t352 + qJ(3) * t349;
t317 = t377 * t350;
t342 = pkin(1) * t350 - pkin(6) * t353;
t394 = -t317 - t342;
t393 = qJD(2) * t350;
t392 = qJD(2) * t353;
t391 = qJD(3) * t349;
t390 = qJD(5) * t352;
t318 = t377 * t353;
t326 = qJD(1) * (pkin(1) * t353 + pkin(6) * t350);
t389 = qJD(1) * t318 + t350 * t391 + t326;
t324 = pkin(3) * t397 + qJ(4) * t353;
t388 = -t324 + t394;
t336 = pkin(2) * t349 - qJ(3) * t352;
t385 = -pkin(3) * t349 - t336;
t384 = qJD(2) * (-rSges(4,1) * t349 + rSges(4,3) * t352 - t336);
t345 = t353 * t391;
t383 = -qJD(4) * t350 + t345;
t382 = pkin(4) * t349 + pkin(7) * t352;
t381 = -qJD(3) * t352 + t317 * t393 + t318 * t392;
t380 = rSges(3,1) * t352 - rSges(3,2) * t349;
t379 = rSges(4,1) * t352 + rSges(4,3) * t349;
t378 = rSges(5,1) * t349 - rSges(5,2) * t352;
t325 = pkin(3) * t395 - qJ(4) * t350;
t376 = qJD(1) * t325 + qJD(4) * t353 + t389;
t357 = qJD(2) * (rSges(5,1) * t352 + rSges(5,2) * t349 + t385);
t356 = qJD(2) * (pkin(4) * t352 - pkin(7) * t349 + t385);
t355 = t324 * t393 + t325 * t392 + t381;
t346 = qJD(5) * t349 + qJD(1);
t341 = rSges(2,1) * t353 - rSges(2,2) * t350;
t339 = rSges(2,1) * t350 + rSges(2,2) * t353;
t338 = rSges(3,1) * t349 + rSges(3,2) * t352;
t323 = t350 * t390 - t392;
t322 = t353 * t390 + t393;
t320 = t382 * t353;
t319 = t382 * t350;
t316 = t349 * t396 - t400;
t315 = -t349 * t399 - t398;
t314 = t349 * t398 + t399;
t313 = -t349 * t400 + t396;
t309 = rSges(3,3) * t350 + t380 * t353;
t308 = rSges(4,2) * t350 + t379 * t353;
t307 = -rSges(5,3) * t350 + t378 * t353;
t306 = rSges(6,3) * t349 + (-rSges(6,1) * t351 + rSges(6,2) * t348) * t352;
t305 = -rSges(3,3) * t353 + t380 * t350;
t304 = -rSges(4,2) * t353 + t379 * t350;
t303 = rSges(5,3) * t353 + t378 * t350;
t294 = Icges(6,5) * t349 + (-Icges(6,1) * t351 + Icges(6,4) * t348) * t352;
t287 = Icges(6,6) * t349 + (-Icges(6,4) * t351 + Icges(6,2) * t348) * t352;
t280 = Icges(6,3) * t349 + (-Icges(6,5) * t351 + Icges(6,6) * t348) * t352;
t279 = rSges(6,1) * t316 + rSges(6,2) * t315 + rSges(6,3) * t395;
t278 = rSges(6,1) * t314 + rSges(6,2) * t313 + rSges(6,3) * t397;
t277 = Icges(6,1) * t316 + Icges(6,4) * t315 + Icges(6,5) * t395;
t276 = Icges(6,1) * t314 + Icges(6,4) * t313 + Icges(6,5) * t397;
t275 = Icges(6,4) * t316 + Icges(6,2) * t315 + Icges(6,6) * t395;
t274 = Icges(6,4) * t314 + Icges(6,2) * t313 + Icges(6,6) * t397;
t273 = Icges(6,5) * t316 + Icges(6,6) * t315 + Icges(6,3) * t395;
t272 = Icges(6,5) * t314 + Icges(6,6) * t313 + Icges(6,3) * t397;
t271 = qJD(1) * t309 - t338 * t393 + t326;
t270 = -t338 * t392 + (-t305 - t342) * qJD(1);
t269 = (t305 * t350 + t309 * t353) * qJD(2);
t268 = qJD(1) * t308 + t350 * t384 + t389;
t267 = t345 + t353 * t384 + (-t304 + t394) * qJD(1);
t266 = (t304 * t350 + t308 * t353) * qJD(2) + t381;
t265 = qJD(1) * t307 + t350 * t357 + t376;
t264 = t353 * t357 + (-t303 + t388) * qJD(1) + t383;
t263 = (t303 * t350 + t307 * t353) * qJD(2) + t355;
t262 = qJD(1) * t320 + t279 * t346 - t306 * t322 + t350 * t356 + t376;
t261 = -t278 * t346 + t306 * t323 + t353 * t356 + (-t319 + t388) * qJD(1) + t383;
t260 = t278 * t322 - t279 * t323 + (t319 * t350 + t320 * t353) * qJD(2) + t355;
t1 = m(3) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(4) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(5) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t322 * ((t273 * t395 + t315 * t275 + t316 * t277) * t322 + (t272 * t395 + t274 * t315 + t276 * t316) * t323 + (t280 * t395 + t287 * t315 + t294 * t316) * t346) / 0.2e1 + t323 * ((t273 * t397 + t275 * t313 + t277 * t314) * t322 + (t272 * t397 + t313 * t274 + t314 * t276) * t323 + (t280 * t397 + t287 * t313 + t294 * t314) * t346) / 0.2e1 + t346 * ((t272 * t323 + t273 * t322 + t280 * t346) * t349 + ((t275 * t348 - t277 * t351) * t322 + (t274 * t348 - t276 * t351) * t323 + (t287 * t348 - t294 * t351) * t346) * t352) / 0.2e1 + (m(2) * (t339 ^ 2 + t341 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t421 * t349 + t423 * t352) * t353 + (t420 * t349 - t422 * t352) * t350) * qJD(2) + (t418 * t349 - t419 * t352) * qJD(1)) * qJD(1) / 0.2e1 + ((t415 * t350 ^ 2 + (t411 * t353 + (t412 - t416) * t350) * t353) * qJD(2) + (-t350 * t414 + t353 * t413) * qJD(1)) * t393 / 0.2e1 - ((t416 * t353 ^ 2 + (t412 * t350 + (t411 - t415) * t353) * t350) * qJD(2) + (t350 * t413 + t353 * t414) * qJD(1)) * t392 / 0.2e1;
T = t1;
