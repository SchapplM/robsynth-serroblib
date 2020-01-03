% Calculate kinetic energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:19
% EndTime: 2019-12-31 18:12:20
% DurationCPUTime: 1.22s
% Computational Cost: add. (630->118), mult. (789->183), div. (0->0), fcn. (659->6), ass. (0->76)
t405 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t404 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t403 = -Icges(4,2) - Icges(6,2) - Icges(5,3);
t322 = pkin(7) + qJ(3);
t320 = cos(t322);
t402 = t405 * t320;
t319 = sin(t322);
t401 = t405 * t319;
t400 = -Icges(5,4) + Icges(4,5) + Icges(6,5);
t399 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t398 = t404 * t320 - t401;
t397 = t403 * t319 + t402;
t396 = rSges(6,3) + qJ(5);
t395 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t326 = sin(qJ(1));
t327 = cos(qJ(1));
t394 = t397 * t326 + t399 * t327;
t393 = -t399 * t326 + t397 * t327;
t392 = -t398 * t326 + t400 * t327;
t391 = t400 * t326 + t398 * t327;
t390 = t404 * t319 + t402;
t389 = t403 * t320 - t401;
t388 = t399 * t319 + t400 * t320;
t387 = rSges(6,1) + pkin(4);
t386 = rSges(6,2) * t319 + t396 * t320;
t385 = t388 * t326 - t395 * t327;
t384 = t395 * t326 + t388 * t327;
t383 = t400 * t319 - t399 * t320;
t382 = t389 * t319 + t390 * t320;
t381 = -t393 * t319 + t391 * t320;
t380 = t394 * t319 + t392 * t320;
t324 = cos(pkin(7));
t375 = t324 * pkin(2);
t315 = t326 * pkin(1) - t327 * qJ(2);
t366 = pkin(6) * t327 - t375 * t326 - t315;
t365 = t387 * t326 + t386 * t327;
t364 = t386 * t326 - t387 * t327;
t321 = qJD(2) * t326;
t360 = qJD(4) * t319;
t363 = t327 * t360 + t321;
t362 = qJD(3) * t326;
t361 = qJD(3) * t327;
t349 = pkin(3) * t320 + qJ(4) * t319;
t295 = t349 * t326;
t359 = -t295 + t366;
t308 = t319 * pkin(3) - t320 * qJ(4);
t356 = qJD(3) * (t319 * rSges(5,2) + t320 * rSges(5,3) - t308);
t312 = qJD(1) * (t327 * pkin(1) + t326 * qJ(2));
t355 = -qJD(2) * t327 + qJD(1) * (pkin(6) * t326 + t375 * t327) + t312;
t296 = t349 * t327;
t354 = -qJD(4) * t320 + t295 * t362 + t296 * t361;
t323 = sin(pkin(7));
t353 = rSges(3,1) * t324 - rSges(3,2) * t323;
t352 = rSges(4,1) * t320 - rSges(4,2) * t319;
t351 = -rSges(5,2) * t320 + rSges(5,3) * t319;
t330 = qJD(1) * t296 + t326 * t360 + t355;
t329 = qJD(5) * t320 + (t320 * rSges(6,2) - t396 * t319 - t308) * qJD(3);
t317 = t327 * rSges(2,1) - t326 * rSges(2,2);
t316 = t326 * rSges(2,1) + t327 * rSges(2,2);
t310 = t319 * rSges(4,1) + t320 * rSges(4,2);
t293 = -t327 * rSges(5,1) + t351 * t326;
t291 = t326 * rSges(5,1) + t351 * t327;
t289 = t326 * rSges(4,3) + t352 * t327;
t288 = -t327 * rSges(4,3) + t352 * t326;
t265 = qJD(1) * t326 * rSges(3,3) + t312 + (qJD(1) * t353 - qJD(2)) * t327;
t264 = t321 + (t327 * rSges(3,3) - t353 * t326 - t315) * qJD(1);
t263 = (t288 * t326 + t289 * t327) * qJD(3);
t262 = qJD(1) * t289 - t310 * t362 + t355;
t261 = -t310 * t361 + t321 + (-t288 + t366) * qJD(1);
t260 = (t291 * t327 + t293 * t326) * qJD(3) + t354;
t259 = qJD(1) * t291 + t326 * t356 + t330;
t258 = t327 * t356 + (-t293 + t359) * qJD(1) + t363;
t257 = qJD(5) * t319 + (t364 * t326 + t365 * t327) * qJD(3) + t354;
t256 = t365 * qJD(1) + t329 * t326 + t330;
t255 = t329 * t327 + (t359 - t364) * qJD(1) + t363;
t1 = m(3) * (t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(4) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(5) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(6) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + (m(2) * (t316 ^ 2 + t317 ^ 2) + Icges(2,3) + Icges(3,2) * t324 ^ 2 + (Icges(3,1) * t323 + 0.2e1 * Icges(3,4) * t324) * t323) * qJD(1) ^ 2 / 0.2e1 + (((t392 * t319 - t394 * t320) * t327 + (t391 * t319 + t393 * t320) * t326) * qJD(3) + (t390 * t319 - t389 * t320) * qJD(1)) * qJD(1) / 0.2e1 + ((t384 * t326 ^ 2 + (t380 * t327 + (t381 - t385) * t326) * t327) * qJD(3) + (t383 * t326 + t382 * t327) * qJD(1)) * t362 / 0.2e1 - ((t385 * t327 ^ 2 + (t381 * t326 + (t380 - t384) * t327) * t326) * qJD(3) + (t382 * t326 - t383 * t327) * qJD(1)) * t361 / 0.2e1;
T = t1;
