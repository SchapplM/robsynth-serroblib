% Calculate kinetic energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:16
% EndTime: 2019-12-31 18:14:17
% DurationCPUTime: 1.18s
% Computational Cost: add. (505->141), mult. (749->201), div. (0->0), fcn. (619->6), ass. (0->84)
t405 = -Icges(5,4) + Icges(6,5);
t404 = Icges(5,1) + Icges(6,1);
t403 = Icges(5,2) + Icges(6,3);
t325 = qJ(3) + pkin(7);
t321 = cos(t325);
t402 = t405 * t321;
t320 = sin(t325);
t401 = t405 * t320;
t400 = Icges(6,4) + Icges(5,5);
t399 = Icges(5,6) - Icges(6,6);
t398 = -t403 * t321 + t401;
t397 = t404 * t320 - t402;
t396 = rSges(6,1) + pkin(4);
t395 = rSges(6,3) + qJ(5);
t328 = sin(qJ(1));
t330 = cos(qJ(1));
t394 = t398 * t328 - t399 * t330;
t393 = -t399 * t328 - t398 * t330;
t392 = t397 * t328 + t400 * t330;
t391 = t400 * t328 - t397 * t330;
t390 = t403 * t320 + t402;
t389 = t404 * t321 + t401;
t388 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t327 = sin(qJ(3));
t329 = cos(qJ(3));
t387 = Icges(4,5) * t327 + Icges(4,6) * t329 + t400 * t320 + t399 * t321;
t386 = t396 * t320 - t395 * t321;
t385 = t387 * t328 + t388 * t330;
t384 = t388 * t328 - t387 * t330;
t383 = Icges(4,5) * t329 - Icges(4,6) * t327 - t399 * t320 + t400 * t321;
t371 = Icges(4,4) * t329;
t312 = -Icges(4,2) * t327 + t371;
t372 = Icges(4,4) * t327;
t313 = Icges(4,1) * t329 - t372;
t382 = t312 * t329 + t313 * t327 + t389 * t320 - t390 * t321;
t346 = Icges(4,2) * t329 + t372;
t292 = Icges(4,6) * t328 - t346 * t330;
t349 = Icges(4,1) * t327 + t371;
t294 = Icges(4,5) * t328 - t349 * t330;
t381 = t292 * t329 + t294 * t327 + t391 * t320 - t393 * t321;
t291 = Icges(4,6) * t330 + t346 * t328;
t293 = Icges(4,5) * t330 + t349 * t328;
t380 = -t291 * t329 - t293 * t327 - t392 * t320 + t394 * t321;
t376 = pkin(3) * t327;
t375 = pkin(3) * t329;
t366 = rSges(6,2) * t330 + t386 * t328;
t365 = rSges(6,2) * t328 - t386 * t330;
t364 = t395 * t320 + t396 * t321;
t310 = qJD(1) * (pkin(1) * t330 + qJ(2) * t328);
t363 = qJD(1) * t330 * pkin(6) + t310;
t362 = qJD(3) * t328;
t361 = qJD(3) * t330;
t360 = qJD(5) * t321;
t324 = qJD(2) * t328;
t359 = qJD(4) * t330 + t362 * t375 + t324;
t314 = pkin(1) * t328 - qJ(2) * t330;
t356 = -pkin(6) * t328 - t314;
t300 = qJ(4) * t330 + t328 * t376;
t355 = qJD(1) * t300 + qJD(4) * t328 + t363;
t299 = qJ(4) * t328 - t330 * t376;
t354 = -t299 + t356;
t353 = rSges(4,1) * t327 + rSges(4,2) * t329;
t352 = rSges(5,1) * t320 + rSges(5,2) * t321;
t317 = rSges(2,1) * t330 - rSges(2,2) * t328;
t316 = rSges(4,1) * t329 - rSges(4,2) * t327;
t315 = rSges(2,1) * t328 + rSges(2,2) * t330;
t309 = rSges(5,1) * t321 - rSges(5,2) * t320;
t298 = rSges(4,3) * t328 - t353 * t330;
t297 = rSges(4,3) * t330 + t353 * t328;
t287 = t299 * t361;
t286 = rSges(5,3) * t328 - t352 * t330;
t284 = rSges(5,3) * t330 + t352 * t328;
t270 = t310 - qJD(2) * t330 + qJD(1) * (-rSges(3,2) * t330 + rSges(3,3) * t328);
t269 = t324 + (rSges(3,2) * t328 + rSges(3,3) * t330 - t314) * qJD(1);
t268 = (-t297 * t328 + t298 * t330) * qJD(3);
t267 = qJD(1) * t297 + (-qJD(3) * t316 - qJD(2)) * t330 + t363;
t266 = t316 * t362 + t324 + (-t298 + t356) * qJD(1);
t265 = qJD(1) * t284 + (-qJD(2) + (-t309 - t375) * qJD(3)) * t330 + t355;
t264 = t309 * t362 + (-t286 + t354) * qJD(1) + t359;
t263 = t287 + (t286 * t330 + (-t284 - t300) * t328) * qJD(3);
t262 = t366 * qJD(1) + (t360 - qJD(2) + (-t364 - t375) * qJD(3)) * t330 + t355;
t261 = (t364 * qJD(3) - t360) * t328 + (t354 - t365) * qJD(1) + t359;
t260 = qJD(5) * t320 + t287 + (t365 * t330 + (-t300 - t366) * t328) * qJD(3);
t1 = m(3) * (t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(4) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(5) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + (m(2) * (t315 ^ 2 + t317 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1 + (((-t291 * t327 + t293 * t329 + t394 * t320 + t392 * t321) * t330 + (-t292 * t327 + t294 * t329 + t393 * t320 + t391 * t321) * t328) * qJD(3) + (-t327 * t312 + t329 * t313 + t390 * t320 + t389 * t321) * qJD(1)) * qJD(1) / 0.2e1 + ((t384 * t328 ^ 2 + (t380 * t330 + (-t381 + t385) * t328) * t330) * qJD(3) + (t383 * t328 - t382 * t330) * qJD(1)) * t362 / 0.2e1 + ((t385 * t330 ^ 2 + (t381 * t328 + (-t380 + t384) * t330) * t328) * qJD(3) + (t382 * t328 + t383 * t330) * qJD(1)) * t361 / 0.2e1;
T = t1;
