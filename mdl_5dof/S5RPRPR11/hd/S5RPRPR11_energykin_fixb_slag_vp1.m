% Calculate kinetic energy for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:09
% EndTime: 2019-12-31 18:27:10
% DurationCPUTime: 1.61s
% Computational Cost: add. (792->165), mult. (1004->261), div. (0->0), fcn. (952->8), ass. (0->97)
t407 = Icges(4,4) - Icges(5,5);
t406 = Icges(4,1) + Icges(5,1);
t405 = Icges(4,2) + Icges(5,3);
t334 = pkin(8) + qJ(3);
t332 = cos(t334);
t404 = t407 * t332;
t331 = sin(t334);
t403 = t407 * t331;
t402 = Icges(5,4) + Icges(4,5);
t401 = Icges(4,6) - Icges(5,6);
t400 = t405 * t331 - t404;
t399 = t406 * t332 - t403;
t398 = Icges(5,2) + Icges(4,3);
t339 = sin(qJ(1));
t341 = cos(qJ(1));
t397 = t400 * t339 + t401 * t341;
t396 = -t401 * t339 + t400 * t341;
t395 = -t399 * t339 + t402 * t341;
t394 = t402 * t339 + t399 * t341;
t393 = -t405 * t332 - t403;
t392 = t406 * t331 + t404;
t391 = -t401 * t331 + t402 * t332;
t390 = t391 * t339 - t341 * t398;
t389 = t339 * t398 + t391 * t341;
t388 = t402 * t331 + t401 * t332;
t387 = t331 * t393 + t332 * t392;
t386 = t331 * t396 + t332 * t394;
t385 = -t331 * t397 + t332 * t395;
t381 = pkin(4) * t332;
t336 = cos(pkin(8));
t379 = pkin(2) * t336;
t327 = pkin(1) * t339 - qJ(2) * t341;
t373 = pkin(6) * t341 - t339 * t379 - t327;
t333 = qJD(2) * t339;
t369 = qJD(4) * t331;
t372 = t341 * t369 + t333;
t371 = qJD(3) * t339;
t370 = qJD(3) * t341;
t368 = qJD(3) - qJD(5);
t358 = pkin(3) * t332 + qJ(4) * t331;
t307 = t358 * t339;
t367 = -t307 + t373;
t319 = pkin(3) * t331 - qJ(4) * t332;
t364 = qJD(3) * (-rSges(5,1) * t331 + rSges(5,3) * t332 - t319);
t322 = qJD(1) * (pkin(1) * t341 + qJ(2) * t339);
t363 = -qJD(2) * t341 + qJD(1) * (pkin(6) * t339 + t341 * t379) + t322;
t308 = t358 * t341;
t362 = -qJD(4) * t332 + t307 * t371 + t308 * t370;
t335 = sin(pkin(8));
t361 = rSges(3,1) * t336 - rSges(3,2) * t335;
t360 = rSges(4,1) * t332 - rSges(4,2) * t331;
t359 = rSges(5,1) * t332 + rSges(5,3) * t331;
t357 = qJD(3) * (-pkin(4) * t331 - t319);
t338 = sin(qJ(5));
t340 = cos(qJ(5));
t310 = t331 * t340 - t332 * t338;
t344 = t331 * t338 + t332 * t340;
t343 = qJD(1) * t308 + t339 * t369 + t363;
t329 = rSges(2,1) * t341 - rSges(2,2) * t339;
t328 = rSges(2,1) * t339 + rSges(2,2) * t341;
t324 = t368 * t341;
t323 = t368 * t339;
t321 = rSges(4,1) * t331 + rSges(4,2) * t332;
t318 = -pkin(7) * t339 + t341 * t381;
t317 = pkin(7) * t341 + t339 * t381;
t305 = t344 * t341;
t304 = t310 * t341;
t303 = t344 * t339;
t302 = t310 * t339;
t301 = rSges(4,3) * t339 + t341 * t360;
t300 = rSges(5,2) * t339 + t341 * t359;
t299 = -rSges(4,3) * t341 + t339 * t360;
t298 = -rSges(5,2) * t341 + t339 * t359;
t281 = rSges(6,1) * t310 - rSges(6,2) * t344;
t280 = Icges(6,1) * t310 - Icges(6,4) * t344;
t279 = Icges(6,4) * t310 - Icges(6,2) * t344;
t278 = Icges(6,5) * t310 - Icges(6,6) * t344;
t277 = qJD(1) * t339 * rSges(3,3) + t322 + (qJD(1) * t361 - qJD(2)) * t341;
t276 = t333 + (t341 * rSges(3,3) - t339 * t361 - t327) * qJD(1);
t275 = rSges(6,1) * t305 + rSges(6,2) * t304 - rSges(6,3) * t339;
t274 = rSges(6,1) * t303 + rSges(6,2) * t302 + rSges(6,3) * t341;
t273 = Icges(6,1) * t305 + Icges(6,4) * t304 - Icges(6,5) * t339;
t272 = Icges(6,1) * t303 + Icges(6,4) * t302 + Icges(6,5) * t341;
t271 = Icges(6,4) * t305 + Icges(6,2) * t304 - Icges(6,6) * t339;
t270 = Icges(6,4) * t303 + Icges(6,2) * t302 + Icges(6,6) * t341;
t269 = Icges(6,5) * t305 + Icges(6,6) * t304 - Icges(6,3) * t339;
t268 = Icges(6,5) * t303 + Icges(6,6) * t302 + Icges(6,3) * t341;
t267 = (t299 * t339 + t301 * t341) * qJD(3);
t266 = qJD(1) * t301 - t321 * t371 + t363;
t265 = -t321 * t370 + t333 + (-t299 + t373) * qJD(1);
t264 = (t298 * t339 + t300 * t341) * qJD(3) + t362;
t263 = qJD(1) * t300 + t339 * t364 + t343;
t262 = t341 * t364 + (-t298 + t367) * qJD(1) + t372;
t261 = t274 * t323 + t275 * t324 + (t317 * t339 + t318 * t341) * qJD(3) + t362;
t260 = -t281 * t323 + t339 * t357 + (t275 + t318) * qJD(1) + t343;
t259 = -t281 * t324 + t341 * t357 + (-t274 - t317 + t367) * qJD(1) + t372;
t1 = m(3) * (t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(4) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(5) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + m(6) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t323 * ((-t269 * t339 + t271 * t304 + t273 * t305) * t323 - (-t268 * t339 + t270 * t304 + t272 * t305) * t324 + (-t278 * t339 + t279 * t304 + t280 * t305) * qJD(1)) / 0.2e1 - t324 * ((t269 * t341 + t271 * t302 + t273 * t303) * t323 - (t268 * t341 + t270 * t302 + t272 * t303) * t324 + (t278 * t341 + t279 * t302 + t280 * t303) * qJD(1)) / 0.2e1 + ((t389 * t339 ^ 2 + (t385 * t341 + (t386 - t390) * t339) * t341) * qJD(3) + (t388 * t339 + t387 * t341) * qJD(1)) * t371 / 0.2e1 - ((t390 * t341 ^ 2 + (t386 * t339 + (t385 - t389) * t341) * t339) * qJD(3) + (t387 * t339 - t388 * t341) * qJD(1)) * t370 / 0.2e1 + (m(2) * (t328 ^ 2 + t329 ^ 2) + Icges(2,3) + Icges(3,2) * t336 ^ 2 + (Icges(3,1) * t335 + 0.2e1 * Icges(3,4) * t336) * t335) * qJD(1) ^ 2 / 0.2e1 + ((-t271 * t344 + t273 * t310) * t323 - (-t270 * t344 + t272 * t310) * t324 + ((t331 * t395 + t332 * t397) * t341 + (t331 * t394 - t332 * t396) * t339) * qJD(3) + (-t344 * t279 + t310 * t280 + t392 * t331 - t393 * t332) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
