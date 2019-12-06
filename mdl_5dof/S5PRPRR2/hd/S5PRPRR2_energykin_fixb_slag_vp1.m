% Calculate kinetic energy for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:45
% EndTime: 2019-12-05 15:44:46
% DurationCPUTime: 1.02s
% Computational Cost: add. (867->140), mult. (948->244), div. (0->0), fcn. (896->10), ass. (0->83)
t395 = Icges(3,3) + Icges(4,3);
t333 = qJ(2) + pkin(9);
t328 = sin(t333);
t329 = cos(t333);
t338 = sin(qJ(2));
t340 = cos(qJ(2));
t394 = Icges(3,5) * t340 + Icges(4,5) * t329 - Icges(3,6) * t338 - Icges(4,6) * t328;
t334 = sin(pkin(8));
t335 = cos(pkin(8));
t393 = t335 * t334;
t392 = t394 * t334 - t395 * t335;
t391 = t395 * t334 + t394 * t335;
t385 = t335 ^ 2;
t386 = t334 ^ 2;
t388 = t385 + t386;
t387 = qJD(2) * t388;
t384 = qJD(2) ^ 2;
t383 = pkin(2) * t338;
t382 = t340 * pkin(2);
t330 = qJ(4) + t333;
t322 = sin(t330);
t380 = t322 * t334;
t379 = t322 * t335;
t337 = sin(qJ(5));
t378 = t334 * t337;
t339 = cos(qJ(5));
t377 = t334 * t339;
t376 = t335 * t337;
t375 = t335 * t339;
t374 = pkin(3) * t329;
t327 = qJD(2) * t334;
t318 = qJD(4) * t334 + t327;
t372 = qJD(2) * t335;
t371 = qJD(3) * t335;
t370 = qJD(5) * t322;
t323 = cos(t330);
t369 = qJD(5) * t323;
t368 = qJD(1) + (-qJ(3) * t335 + t382 * t334) * t327 + (qJ(3) * t334 + t382 * t335) * t372;
t319 = (-qJD(2) - qJD(4)) * t335;
t364 = rSges(5,1) * t323 - rSges(5,2) * t322;
t363 = qJD(2) * (-t328 * rSges(4,1) - t329 * rSges(4,2) - t383);
t362 = (-pkin(6) * t335 + t374 * t334) * t327 + (pkin(6) * t334 + t374 * t335) * t372 + t368;
t359 = Icges(5,1) * t323 - Icges(5,4) * t322;
t356 = Icges(5,4) * t323 - Icges(5,2) * t322;
t353 = Icges(5,5) * t323 - Icges(5,6) * t322;
t352 = (-Icges(5,3) * t335 + t353 * t334) * t319 + (Icges(5,3) * t334 + t353 * t335) * t318;
t347 = qJD(2) * (-pkin(3) * t328 - t383);
t326 = qJD(3) * t334;
t344 = t335 * t347 + t326;
t343 = t334 * t347 - t371;
t342 = (-(Icges(5,6) * t334 + t356 * t335) * t322 + (Icges(5,5) * t334 + t359 * t335) * t323) * t318 + (-(-Icges(5,6) * t335 + t356 * t334) * t322 + (-Icges(5,5) * t335 + t359 * t334) * t323) * t319;
t321 = t338 * rSges(3,1) + t340 * rSges(3,2);
t315 = t322 * pkin(4) - t323 * pkin(7);
t314 = t322 * rSges(5,1) + t323 * rSges(5,2);
t313 = t323 * t375 + t378;
t312 = -t323 * t376 + t377;
t311 = t323 * t377 - t376;
t310 = -t323 * t378 - t375;
t303 = t334 * t370 + t319;
t302 = t335 * t370 + t318;
t289 = -t323 * rSges(6,3) + (rSges(6,1) * t339 - rSges(6,2) * t337) * t322;
t288 = -Icges(6,5) * t323 + (Icges(6,1) * t339 - Icges(6,4) * t337) * t322;
t287 = -Icges(6,6) * t323 + (Icges(6,4) * t339 - Icges(6,2) * t337) * t322;
t286 = -Icges(6,3) * t323 + (Icges(6,5) * t339 - Icges(6,6) * t337) * t322;
t283 = t335 * t363 + t326;
t282 = t334 * t363 - t371;
t279 = qJD(1) + (rSges(3,1) * t340 - rSges(3,2) * t338) * t387;
t278 = t319 * t314 + t344;
t277 = -t318 * t314 + t343;
t276 = t313 * rSges(6,1) + t312 * rSges(6,2) + rSges(6,3) * t379;
t275 = t311 * rSges(6,1) + t310 * rSges(6,2) + rSges(6,3) * t380;
t274 = Icges(6,1) * t313 + Icges(6,4) * t312 + Icges(6,5) * t379;
t273 = Icges(6,1) * t311 + Icges(6,4) * t310 + Icges(6,5) * t380;
t272 = Icges(6,4) * t313 + Icges(6,2) * t312 + Icges(6,6) * t379;
t271 = Icges(6,4) * t311 + Icges(6,2) * t310 + Icges(6,6) * t380;
t270 = Icges(6,5) * t313 + Icges(6,6) * t312 + Icges(6,3) * t379;
t269 = Icges(6,5) * t311 + Icges(6,6) * t310 + Icges(6,3) * t380;
t268 = t368 + (rSges(4,1) * t329 - rSges(4,2) * t328) * t387;
t267 = t275 * t369 + t303 * t289 + t319 * t315 + t344;
t266 = -t276 * t369 - t302 * t289 - t318 * t315 + t343;
t265 = t318 * (-t335 * rSges(5,3) + t364 * t334) - t319 * (t334 * rSges(5,3) + t364 * t335) + t362;
t264 = t302 * t275 - t303 * t276 + t362 + (t334 * t318 - t335 * t319) * (pkin(4) * t323 + pkin(7) * t322);
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t321 ^ 2 * t384 * t388 + t279 ^ 2) / 0.2e1 + m(4) * (t268 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t318 * (t352 * t334 + t342 * t335) / 0.2e1 + t319 * (t342 * t334 - t352 * t335) / 0.2e1 + m(6) * (t264 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t302 * ((t270 * t379 + t312 * t272 + t313 * t274) * t302 + (t269 * t379 + t312 * t271 + t313 * t273) * t303 - (t286 * t379 + t312 * t287 + t313 * t288) * t369) / 0.2e1 + t303 * ((t270 * t380 + t310 * t272 + t311 * t274) * t302 + (t269 * t380 + t310 * t271 + t311 * t273) * t303 - (t286 * t380 + t310 * t287 + t311 * t288) * t369) / 0.2e1 - ((-t269 * t303 - t270 * t302 + t286 * t369) * t323 + ((-t272 * t337 + t274 * t339) * t302 + (-t271 * t337 + t273 * t339) * t303 - (-t287 * t337 + t288 * t339) * t369) * t322) * t369 / 0.2e1 + (t391 * t386 - t392 * t393) * t334 * t384 / 0.2e1 - (t392 * t385 - t391 * t393) * t335 * t384 / 0.2e1;
T = t1;
