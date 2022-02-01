% Calculate kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:02
% EndTime: 2022-01-23 08:59:03
% DurationCPUTime: 0.70s
% Computational Cost: add. (710->170), mult. (1628->262), div. (0->0), fcn. (1947->10), ass. (0->90)
t369 = sin(qJ(1));
t371 = cos(qJ(1));
t352 = t369 * pkin(1) - t371 * qJ(2);
t367 = cos(pkin(7));
t364 = sin(pkin(7));
t399 = t364 * qJ(3) + pkin(1);
t350 = pkin(2) * t367 + t399;
t400 = -pkin(1) + t350;
t383 = -t400 * t369 - t352;
t366 = cos(pkin(8));
t363 = sin(pkin(8));
t398 = t363 * qJ(4) + pkin(2);
t334 = (t366 * pkin(3) + t398) * t367 + t399;
t384 = t334 - t350;
t377 = qJ(4) * t366 - qJ(2);
t349 = -t363 * pkin(3) + t377;
t386 = qJ(2) + t349;
t403 = -t384 * t369 - t386 * t371 + t383;
t402 = t367 ^ 2;
t397 = t363 * t364;
t368 = sin(qJ(5));
t396 = t363 * t368;
t370 = cos(qJ(5));
t395 = t363 * t370;
t394 = t364 * t366;
t393 = t364 * t369;
t392 = t364 * t371;
t365 = cos(pkin(9));
t391 = t367 * t365;
t390 = t369 * t363;
t389 = t369 * t366;
t388 = t371 * t363;
t387 = t371 * t366;
t362 = sin(pkin(9));
t351 = t365 * pkin(4) + t362 * pkin(6) + pkin(3);
t385 = -t334 + (t351 * t366 + t398) * t367 + pkin(1) + (t362 * pkin(4) - t365 * pkin(6) + qJ(3)) * t364;
t360 = qJD(2) * t369;
t381 = qJD(3) * t364;
t382 = t371 * t381 + t360;
t344 = t367 * t389 - t388;
t331 = t344 * t362 - t365 * t393;
t380 = qJD(5) * t331;
t346 = t367 * t387 + t390;
t332 = t346 * t362 - t365 * t392;
t379 = qJD(5) * t332;
t345 = t367 * t388 - t389;
t378 = qJD(4) * t345 + t382;
t347 = -qJD(3) * t367 + qJD(4) * t397;
t376 = rSges(3,1) * t367 - rSges(3,2) * t364;
t348 = qJD(1) * (t371 * pkin(1) + t369 * qJ(2));
t375 = t369 * t381 + t348 + (qJD(1) * t400 - qJD(2)) * t371;
t374 = -t351 * t363 - t349 + t377;
t343 = t367 * t390 + t387;
t373 = qJD(1) * (-t386 * t369 + t384 * t371) + qJD(4) * t343 + t375;
t354 = t371 * rSges(2,1) - t369 * rSges(2,2);
t353 = t369 * rSges(2,1) + t371 * rSges(2,2);
t342 = t365 * t396 + t370 * t366;
t341 = t364 * t362 + t366 * t391;
t340 = -t367 * t362 + t365 * t394;
t339 = t362 * t394 + t391;
t335 = qJD(5) * t339 + qJD(1);
t330 = -t341 * t368 + t367 * t395;
t329 = -t340 * t368 + t364 * t395;
t328 = t340 * t370 + t364 * t396;
t327 = qJD(1) * t369 * rSges(3,3) + t348 + (qJD(1) * t376 - qJD(2)) * t371;
t326 = t360 + (t371 * rSges(3,3) - t376 * t369 - t352) * qJD(1);
t324 = t330 * t371 - t369 * t342;
t323 = t330 * t369 + t371 * t342;
t322 = (t341 * t370 + t367 * t396) * t371 + t369 * (t365 * t395 - t368 * t366);
t321 = (t341 * t369 - t365 * t388) * t370 + t343 * t368;
t318 = t328 * rSges(6,1) + t329 * rSges(6,2) + t339 * rSges(6,3);
t317 = Icges(6,1) * t328 + Icges(6,4) * t329 + Icges(6,5) * t339;
t316 = Icges(6,4) * t328 + Icges(6,2) * t329 + Icges(6,6) * t339;
t315 = Icges(6,5) * t328 + Icges(6,6) * t329 + Icges(6,3) * t339;
t314 = qJD(1) * (t346 * rSges(4,1) - t345 * rSges(4,2) + rSges(4,3) * t392) + t375;
t313 = (-t344 * rSges(4,1) + t343 * rSges(4,2) - rSges(4,3) * t393 + t383) * qJD(1) + t382;
t312 = t322 * rSges(6,1) + t324 * rSges(6,2) + t332 * rSges(6,3);
t311 = t321 * rSges(6,1) + t323 * rSges(6,2) + t331 * rSges(6,3);
t310 = Icges(6,1) * t322 + Icges(6,4) * t324 + Icges(6,5) * t332;
t309 = Icges(6,1) * t321 + Icges(6,4) * t323 + Icges(6,5) * t331;
t308 = Icges(6,4) * t322 + Icges(6,2) * t324 + Icges(6,6) * t332;
t307 = Icges(6,4) * t321 + Icges(6,2) * t323 + Icges(6,6) * t331;
t306 = Icges(6,5) * t322 + Icges(6,6) * t324 + Icges(6,3) * t332;
t305 = Icges(6,5) * t321 + Icges(6,6) * t323 + Icges(6,3) * t331;
t304 = qJD(1) * ((t346 * t365 + t362 * t392) * rSges(5,1) - t332 * rSges(5,2) + t345 * rSges(5,3)) + t373;
t303 = (-(t344 * t365 + t362 * t393) * rSges(5,1) + t331 * rSges(5,2) - t343 * rSges(5,3) + t403) * qJD(1) + t378;
t302 = (t311 * t332 - t312 * t331) * qJD(5) + t347;
t301 = qJD(1) * (-t374 * t369 + t385 * t371) - t318 * t379 + t335 * t312 + t373;
t300 = -t335 * t311 + t318 * t380 + t378 + (-t385 * t369 - t374 * t371 + t403) * qJD(1);
t1 = m(3) * (t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t402 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(5) * (t303 ^ 2 + t304 ^ 2 + t347 ^ 2) / 0.2e1 + m(6) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + ((t332 * t315 + t324 * t316 + t322 * t317) * t335 + ((t332 * t306 + t324 * t308 + t322 * t310) * t332 + (t332 * t305 + t324 * t307 + t322 * t309) * t331) * qJD(5)) * t379 / 0.2e1 + ((t331 * t315 + t323 * t316 + t321 * t317) * t335 + ((t331 * t306 + t323 * t308 + t321 * t310) * t332 + (t331 * t305 + t323 * t307 + t321 * t309) * t331) * qJD(5)) * t380 / 0.2e1 + t335 * ((t339 * t315 + t329 * t316 + t328 * t317) * t335 + ((t339 * t306 + t329 * t308 + t328 * t310) * t332 + (t339 * t305 + t329 * t307 + t328 * t309) * t331) * qJD(5)) / 0.2e1 + (m(2) * (t353 ^ 2 + t354 ^ 2) + Icges(2,3) + (Icges(5,1) * t340 + 0.2e1 * Icges(5,5) * t397) * t340 + (-0.2e1 * Icges(5,4) * t340 + Icges(5,2) * t339 - 0.2e1 * Icges(5,6) * t397) * t339 + (Icges(3,2) + Icges(4,3)) * t402 + ((Icges(4,1) * t366 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t366 + (Icges(4,2) + Icges(5,3)) * t363) * t363) * t364 + 0.2e1 * (-Icges(4,5) * t366 + Icges(4,6) * t363 + Icges(3,4)) * t367) * t364) * qJD(1) ^ 2 / 0.2e1;
T = t1;
