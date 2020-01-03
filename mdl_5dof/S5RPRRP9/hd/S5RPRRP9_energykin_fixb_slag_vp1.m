% Calculate kinetic energy for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:34
% EndTime: 2019-12-31 18:48:36
% DurationCPUTime: 1.41s
% Computational Cost: add. (873->148), mult. (833->231), div. (0->0), fcn. (703->8), ass. (0->94)
t417 = Icges(5,4) - Icges(6,5);
t416 = Icges(5,1) + Icges(6,1);
t415 = Icges(5,2) + Icges(6,3);
t341 = pkin(8) + qJ(3);
t337 = qJ(4) + t341;
t333 = cos(t337);
t414 = t417 * t333;
t332 = sin(t337);
t413 = t417 * t332;
t412 = Icges(6,4) + Icges(5,5);
t411 = Icges(5,6) - Icges(6,6);
t410 = t415 * t332 - t414;
t409 = t416 * t333 - t413;
t408 = rSges(6,1) + pkin(4);
t407 = rSges(6,3) + qJ(5);
t406 = Icges(6,2) + Icges(5,3);
t345 = sin(qJ(1));
t346 = cos(qJ(1));
t405 = -t410 * t345 - t411 * t346;
t404 = -t411 * t345 + t410 * t346;
t403 = t409 * t345 - t412 * t346;
t402 = t412 * t345 + t409 * t346;
t401 = -t415 * t333 - t413;
t400 = t416 * t332 + t414;
t399 = -t411 * t332 + t412 * t333;
t398 = t407 * t332 + t408 * t333;
t374 = qJD(3) + qJD(4);
t327 = t374 * t345;
t328 = t374 * t346;
t397 = (t405 * t332 - t403 * t333) * t328 + (t404 * t332 + t402 * t333) * t327 + (t401 * t332 + t400 * t333) * qJD(1);
t396 = (-t399 * t345 + t406 * t346) * t328 + (t406 * t345 + t399 * t346) * t327 + (t412 * t332 + t411 * t333) * qJD(1);
t343 = cos(pkin(8));
t391 = t343 * pkin(2);
t335 = sin(t341);
t390 = Icges(4,4) * t335;
t336 = cos(t341);
t389 = Icges(4,4) * t336;
t378 = pkin(3) * t336;
t280 = -pkin(7) * t346 + t378 * t345;
t281 = pkin(7) * t345 + t378 * t346;
t375 = qJD(3) * t346;
t376 = qJD(3) * t345;
t383 = t280 * t376 + t281 * t375;
t382 = -rSges(6,2) * t346 + t398 * t345;
t381 = rSges(6,2) * t345 + t398 * t346;
t329 = pkin(1) * t345 - qJ(2) * t346;
t380 = pkin(6) * t346 - t391 * t345 - t329;
t379 = t408 * t332 - t407 * t333;
t373 = pkin(3) * qJD(3) * t335;
t372 = -t280 + t380;
t326 = qJD(1) * (pkin(1) * t346 + qJ(2) * t345);
t371 = -qJD(2) * t346 + qJD(1) * (pkin(6) * t345 + t391 * t346) + t326;
t342 = sin(pkin(8));
t370 = rSges(3,1) * t343 - rSges(3,2) * t342;
t369 = rSges(4,1) * t336 - rSges(4,2) * t335;
t368 = rSges(5,1) * t333 - rSges(5,2) * t332;
t365 = Icges(4,1) * t336 - t390;
t362 = -Icges(4,2) * t335 + t389;
t359 = Icges(4,5) * t336 - Icges(4,6) * t335;
t304 = -Icges(4,6) * t346 + t362 * t345;
t306 = -Icges(4,5) * t346 + t365 * t345;
t356 = t304 * t335 - t306 * t336;
t305 = Icges(4,6) * t345 + t362 * t346;
t307 = Icges(4,5) * t345 + t365 * t346;
t355 = -t305 * t335 + t307 * t336;
t322 = Icges(4,2) * t336 + t390;
t323 = Icges(4,1) * t335 + t389;
t354 = -t322 * t335 + t323 * t336;
t353 = qJD(1) * t281 + t371;
t352 = qJD(5) * t332 - t373;
t338 = qJD(2) * t345;
t331 = rSges(2,1) * t346 - rSges(2,2) * t345;
t330 = rSges(2,1) * t345 + rSges(2,2) * t346;
t324 = rSges(4,1) * t335 + rSges(4,2) * t336;
t321 = Icges(4,5) * t335 + Icges(4,6) * t336;
t320 = rSges(5,1) * t332 + rSges(5,2) * t333;
t309 = rSges(4,3) * t345 + t369 * t346;
t308 = -rSges(4,3) * t346 + t369 * t345;
t303 = Icges(4,3) * t345 + t359 * t346;
t302 = -Icges(4,3) * t346 + t359 * t345;
t300 = rSges(5,3) * t345 + t368 * t346;
t298 = -rSges(5,3) * t346 + t368 * t345;
t283 = qJD(1) * t345 * rSges(3,3) + t326 + (qJD(1) * t370 - qJD(2)) * t346;
t282 = t338 + (t346 * rSges(3,3) - t370 * t345 - t329) * qJD(1);
t276 = (t308 * t345 + t309 * t346) * qJD(3);
t275 = qJD(1) * t309 - t324 * t376 + t371;
t274 = -t324 * t375 + t338 + (-t308 + t380) * qJD(1);
t273 = qJD(1) * t300 - t320 * t327 - t345 * t373 + t353;
t272 = -t346 * t373 - t320 * t328 + t338 + (-t298 + t372) * qJD(1);
t271 = t298 * t327 + t300 * t328 + t383;
t270 = t381 * qJD(1) - t379 * t327 + t352 * t345 + t353;
t269 = t338 + t352 * t346 - t379 * t328 + (t372 - t382) * qJD(1);
t268 = -qJD(5) * t333 + t382 * t327 + t381 * t328 + t383;
t1 = m(3) * (t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(4) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t345 * t321 + t354 * t346) * qJD(1) + (t345 ^ 2 * t303 + (t356 * t346 + (-t302 + t355) * t345) * t346) * qJD(3)) * t376 / 0.2e1 - ((-t346 * t321 + t354 * t345) * qJD(1) + (t346 ^ 2 * t302 + (t355 * t345 + (-t303 + t356) * t346) * t345) * qJD(3)) * t375 / 0.2e1 + m(5) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + m(6) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + (t396 * t345 + t397 * t346) * t327 / 0.2e1 - (t397 * t345 - t396 * t346) * t328 / 0.2e1 + (m(2) * (t330 ^ 2 + t331 ^ 2) + Icges(2,3) + Icges(3,2) * t343 ^ 2 + (Icges(3,1) * t342 + 0.2e1 * Icges(3,4) * t343) * t342) * qJD(1) ^ 2 / 0.2e1 + (((t305 * t336 + t307 * t335) * t345 - (t304 * t336 + t306 * t335) * t346) * qJD(3) - (t403 * t332 + t405 * t333) * t328 + (t402 * t332 - t404 * t333) * t327 + (t336 * t322 + t335 * t323 + t400 * t332 - t401 * t333) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
