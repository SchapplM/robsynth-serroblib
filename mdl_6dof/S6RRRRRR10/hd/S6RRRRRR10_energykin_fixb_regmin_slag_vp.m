% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:16:47
% EndTime: 2019-03-10 06:16:48
% DurationCPUTime: 0.42s
% Computational Cost: add. (2094->71), mult. (5637->147), div. (0->0), fcn. (4907->16), ass. (0->66)
t355 = cos(pkin(6)) * qJD(1);
t325 = qJD(2) + t355;
t336 = sin(qJ(3));
t331 = cos(pkin(7));
t341 = cos(qJ(2));
t329 = sin(pkin(6));
t356 = qJD(1) * t329;
t348 = t341 * t356;
t347 = t331 * t348;
t337 = sin(qJ(2));
t349 = t337 * t356;
t328 = sin(pkin(7));
t367 = cos(qJ(3));
t350 = t328 * t367;
t310 = -t325 * t350 + t336 * t349 - t367 * t347;
t318 = -t331 * t325 + t328 * t348 - qJD(3);
t327 = sin(pkin(8));
t330 = cos(pkin(8));
t370 = t310 * t330 + t318 * t327;
t354 = pkin(1) * t355;
t357 = pkin(10) * t348 + t337 * t354;
t308 = (t325 * t328 + t347) * pkin(11) + t357;
t324 = t341 * t354;
t309 = t325 * pkin(2) + t324 + (-pkin(11) * t331 - pkin(10)) * t349;
t317 = (-pkin(11) * t328 * t337 - pkin(2) * t341 - pkin(1)) * t356;
t346 = t331 * t367 * t309 - t336 * t308 + t317 * t350;
t359 = t331 * t336;
t360 = t328 * t336;
t311 = t325 * t360 + (t367 * t337 + t341 * t359) * t356;
t366 = pkin(12) * t311;
t292 = -t318 * pkin(3) - t330 * t366 + t346;
t301 = -t328 * t309 + t331 * t317;
t295 = t310 * pkin(3) - t327 * t366 + t301;
t369 = t292 * t330 + t295 * t327;
t351 = t367 * t308 + t309 * t359 + t317 * t360;
t291 = -pkin(12) * t370 + t351;
t335 = sin(qJ(4));
t340 = cos(qJ(4));
t368 = -t335 * t291 + t369 * t340;
t342 = qJD(1) ^ 2;
t361 = t329 ^ 2 * t342;
t302 = -t327 * t310 + t330 * t318 - qJD(4);
t352 = t340 * t291 + t369 * t335;
t279 = -t302 * pkin(13) + t352;
t282 = -t327 * t292 + t330 * t295;
t298 = t335 * t311 + t370 * t340;
t299 = t340 * t311 - t335 * t370;
t281 = t298 * pkin(4) - t299 * pkin(13) + t282;
t334 = sin(qJ(5));
t339 = cos(qJ(5));
t358 = t339 * t279 + t334 * t281;
t353 = t341 * t361;
t289 = t334 * t299 + t339 * t302;
t345 = -t334 * t279 + t339 * t281;
t278 = t302 * pkin(4) - t368;
t338 = cos(qJ(6));
t333 = sin(qJ(6));
t297 = qJD(5) + t298;
t290 = t339 * t299 - t334 * t302;
t288 = qJD(6) + t289;
t284 = t338 * t290 + t333 * t297;
t283 = t333 * t290 - t338 * t297;
t276 = t289 * pkin(5) - t290 * pkin(14) + t278;
t275 = t297 * pkin(14) + t358;
t274 = -t297 * pkin(5) - t345;
t1 = [t342 / 0.2e1, 0, 0, t337 ^ 2 * t361 / 0.2e1, t337 * t353, t325 * t349, t325 * t348, t325 ^ 2 / 0.2e1 (-pkin(10) * t349 + t324) * t325 + pkin(1) * t353, -pkin(1) * t337 * t361 - t357 * t325, t311 ^ 2 / 0.2e1, -t310 * t311, -t318 * t311, t318 * t310, t318 ^ 2 / 0.2e1, t301 * t310 - t346 * t318, t301 * t311 + t351 * t318, t299 ^ 2 / 0.2e1, -t298 * t299, -t302 * t299, t302 * t298, t302 ^ 2 / 0.2e1, t282 * t298 - t302 * t368, t282 * t299 + t352 * t302, t290 ^ 2 / 0.2e1, -t290 * t289, t297 * t290, -t297 * t289, t297 ^ 2 / 0.2e1, t278 * t289 + t345 * t297, t278 * t290 - t358 * t297, t284 ^ 2 / 0.2e1, -t284 * t283, t288 * t284, -t288 * t283, t288 ^ 2 / 0.2e1 (-t333 * t275 + t338 * t276) * t288 + t274 * t283 -(t338 * t275 + t333 * t276) * t288 + t274 * t284;];
T_reg  = t1;
