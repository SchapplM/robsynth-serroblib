% Calculate kinetic energy for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:17
% EndTime: 2019-12-31 17:42:18
% DurationCPUTime: 1.23s
% Computational Cost: add. (883->155), mult. (964->264), div. (0->0), fcn. (912->10), ass. (0->86)
t392 = Icges(4,3) + Icges(5,3);
t334 = qJ(2) + qJ(3);
t329 = pkin(9) + t334;
t323 = sin(t329);
t324 = cos(t329);
t330 = sin(t334);
t331 = cos(t334);
t391 = Icges(4,5) * t331 + Icges(5,5) * t324 - Icges(4,6) * t330 - Icges(5,6) * t323;
t335 = sin(pkin(8));
t336 = cos(pkin(8));
t390 = t335 * t336;
t328 = qJD(2) * t335;
t319 = qJD(3) * t335 + t328;
t320 = (-qJD(2) - qJD(3)) * t336;
t355 = Icges(5,4) * t324 - Icges(5,2) * t323;
t356 = Icges(4,4) * t331 - Icges(4,2) * t330;
t358 = Icges(5,1) * t324 - Icges(5,4) * t323;
t359 = Icges(4,1) * t331 - Icges(4,4) * t330;
t389 = (-(-Icges(4,6) * t336 + t356 * t335) * t330 + (-Icges(4,5) * t336 + t359 * t335) * t331 - (-Icges(5,6) * t336 + t355 * t335) * t323 + (-Icges(5,5) * t336 + t358 * t335) * t324) * t320 + (-(Icges(4,6) * t335 + t356 * t336) * t330 + (Icges(4,5) * t335 + t359 * t336) * t331 - (Icges(5,6) * t335 + t355 * t336) * t323 + (Icges(5,5) * t335 + t358 * t336) * t324) * t319;
t388 = (t391 * t335 - t392 * t336) * t320 + (t392 * t335 + t391 * t336) * t319;
t385 = t336 ^ 2;
t386 = t335 ^ 2;
t387 = t385 + t386;
t384 = qJD(2) ^ 2;
t381 = pkin(3) * t330;
t340 = cos(qJ(2));
t380 = t340 * pkin(2);
t378 = t323 * t335;
t377 = t323 * t336;
t337 = sin(qJ(5));
t376 = t335 * t337;
t339 = cos(qJ(5));
t375 = t335 * t339;
t374 = t336 * t337;
t373 = t336 * t339;
t372 = pkin(3) * t331;
t370 = qJD(5) * t323;
t369 = qJD(5) * t324;
t338 = sin(qJ(2));
t368 = pkin(2) * qJD(2) * t338;
t367 = qJD(1) + (-pkin(6) * t336 + t380 * t335) * t328 + qJD(2) * t336 * (pkin(6) * t335 + t380 * t336);
t366 = t335 * t368;
t365 = t336 * t368;
t364 = t319 * (-qJ(4) * t336 + t372 * t335) + t367;
t363 = pkin(4) * t324 + pkin(7) * t323;
t362 = rSges(4,1) * t331 - rSges(4,2) * t330;
t361 = rSges(5,1) * t324 - rSges(5,2) * t323;
t354 = Icges(3,5) * t340 - Icges(3,6) * t338;
t347 = qJD(4) * t335 + t320 * t381 - t365;
t345 = -qJD(4) * t336 - t366;
t322 = t338 * rSges(3,1) + t340 * rSges(3,2);
t317 = t330 * rSges(4,1) + t331 * rSges(4,2);
t316 = t323 * pkin(4) - t324 * pkin(7);
t315 = t323 * rSges(5,1) + t324 * rSges(5,2);
t313 = t324 * t373 + t376;
t312 = -t324 * t374 + t375;
t311 = t324 * t375 - t374;
t310 = -t324 * t376 - t373;
t305 = Icges(3,3) * t335 + t354 * t336;
t304 = -Icges(3,3) * t336 + t354 * t335;
t303 = t335 * t370 + t320;
t302 = t336 * t370 + t319;
t289 = -t324 * rSges(6,3) + (rSges(6,1) * t339 - rSges(6,2) * t337) * t323;
t288 = -Icges(6,5) * t324 + (Icges(6,1) * t339 - Icges(6,4) * t337) * t323;
t287 = -Icges(6,6) * t324 + (Icges(6,4) * t339 - Icges(6,2) * t337) * t323;
t286 = -Icges(6,3) * t324 + (Icges(6,5) * t339 - Icges(6,6) * t337) * t323;
t283 = t320 * t317 - t365;
t282 = -t319 * t317 - t366;
t281 = qJ(4) * t335 + t372 * t336;
t279 = qJD(1) + t387 * qJD(2) * (rSges(3,1) * t340 - rSges(3,2) * t338);
t278 = t313 * rSges(6,1) + t312 * rSges(6,2) + rSges(6,3) * t377;
t277 = t311 * rSges(6,1) + t310 * rSges(6,2) + rSges(6,3) * t378;
t276 = Icges(6,1) * t313 + Icges(6,4) * t312 + Icges(6,5) * t377;
t275 = Icges(6,1) * t311 + Icges(6,4) * t310 + Icges(6,5) * t378;
t274 = Icges(6,4) * t313 + Icges(6,2) * t312 + Icges(6,6) * t377;
t273 = Icges(6,4) * t311 + Icges(6,2) * t310 + Icges(6,6) * t378;
t272 = Icges(6,5) * t313 + Icges(6,6) * t312 + Icges(6,3) * t377;
t271 = Icges(6,5) * t311 + Icges(6,6) * t310 + Icges(6,3) * t378;
t270 = t320 * t315 + t347;
t269 = (-t315 - t381) * t319 + t345;
t268 = t319 * (-t336 * rSges(4,3) + t362 * t335) - t320 * (t335 * rSges(4,3) + t362 * t336) + t367;
t267 = t277 * t369 + t303 * t289 + t320 * t316 + t347;
t266 = -t278 * t369 - t302 * t289 + (-t316 - t381) * t319 + t345;
t265 = t319 * (-t336 * rSges(5,3) + t361 * t335) + (-t335 * rSges(5,3) - t361 * t336 - t281) * t320 + t364;
t264 = t302 * t277 - t303 * t278 + t319 * t363 * t335 + (-t363 * t336 - t281) * t320 + t364;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t387 * t384 * t322 ^ 2 + t279 ^ 2) / 0.2e1 + t384 * t335 * (-t304 * t390 + t386 * t305) / 0.2e1 - t384 * t336 * (t385 * t304 - t305 * t390) / 0.2e1 + m(4) * (t268 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(6) * (t264 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t302 * ((t272 * t377 + t312 * t274 + t313 * t276) * t302 + (t271 * t377 + t312 * t273 + t313 * t275) * t303 - (t286 * t377 + t312 * t287 + t313 * t288) * t369) / 0.2e1 + t303 * ((t272 * t378 + t310 * t274 + t311 * t276) * t302 + (t271 * t378 + t310 * t273 + t311 * t275) * t303 - (t286 * t378 + t310 * t287 + t311 * t288) * t369) / 0.2e1 - ((-t271 * t303 - t272 * t302 + t286 * t369) * t324 + ((-t274 * t337 + t276 * t339) * t302 + (-t273 * t337 + t275 * t339) * t303 - (-t287 * t337 + t288 * t339) * t369) * t323) * t369 / 0.2e1 + (t388 * t335 + t389 * t336) * t319 / 0.2e1 + (t389 * t335 - t388 * t336) * t320 / 0.2e1;
T = t1;
