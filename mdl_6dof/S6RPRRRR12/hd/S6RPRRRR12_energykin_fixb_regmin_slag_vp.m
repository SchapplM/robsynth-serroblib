% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:40
% EndTime: 2019-03-09 08:00:41
% DurationCPUTime: 0.39s
% Computational Cost: add. (1699->74), mult. (5653->151), div. (0->0), fcn. (4918->16), ass. (0->64)
t344 = sin(qJ(3));
t333 = sin(pkin(14));
t336 = sin(pkin(6));
t358 = qJD(1) * t336;
t353 = t333 * t358;
t370 = cos(qJ(3));
t335 = sin(pkin(7));
t339 = cos(pkin(7));
t340 = cos(pkin(6));
t337 = cos(pkin(14));
t361 = t336 * t337;
t373 = qJD(1) * (t335 * t340 + t339 * t361);
t315 = t344 * t353 - t370 * t373;
t323 = -qJD(3) - (-t335 * t361 + t339 * t340) * qJD(1);
t334 = sin(pkin(8));
t338 = cos(pkin(8));
t375 = t315 * t338 + t323 * t334;
t357 = pkin(1) * qJD(1) * t340;
t326 = t337 * qJ(2) * t358 + t333 * t357;
t313 = pkin(10) * t373 + t326;
t331 = t337 * t357;
t364 = t333 * t336;
t314 = t331 + (pkin(2) * t340 + (-pkin(10) * t339 - qJ(2)) * t364) * qJD(1);
t320 = qJD(2) + (-pkin(10) * t333 * t335 - pkin(2) * t337 - pkin(1)) * t358;
t351 = -t344 * t313 + (t314 * t339 + t320 * t335) * t370;
t360 = t339 * t344;
t362 = t335 * t344;
t316 = (t340 * t362 + (t333 * t370 + t337 * t360) * t336) * qJD(1);
t369 = pkin(11) * t316;
t294 = -t323 * pkin(3) - t338 * t369 + t351;
t305 = -t335 * t314 + t339 * t320;
t300 = t315 * pkin(3) - t334 * t369 + t305;
t374 = t294 * t338 + t300 * t334;
t354 = t313 * t370 + t314 * t360 + t320 * t362;
t293 = -pkin(11) * t375 + t354;
t343 = sin(qJ(4));
t347 = cos(qJ(4));
t371 = -t343 * t293 + t347 * t374;
t307 = -t334 * t315 + t338 * t323 - qJD(4);
t355 = t347 * t293 + t374 * t343;
t284 = -t307 * pkin(12) + t355;
t287 = -t334 * t294 + t338 * t300;
t303 = t343 * t316 + t375 * t347;
t304 = t347 * t316 - t343 * t375;
t286 = t303 * pkin(4) - t304 * pkin(12) + t287;
t342 = sin(qJ(5));
t346 = cos(qJ(5));
t359 = t346 * t284 + t342 * t286;
t296 = t342 * t304 + t346 * t307;
t350 = -t342 * t284 + t346 * t286;
t283 = t307 * pkin(4) - t371;
t345 = cos(qJ(6));
t341 = sin(qJ(6));
t332 = -pkin(1) * t358 + qJD(2);
t325 = -qJ(2) * t353 + t331;
t302 = qJD(5) + t303;
t297 = t346 * t304 - t342 * t307;
t295 = qJD(6) + t296;
t289 = t345 * t297 + t341 * t302;
t288 = t341 * t297 - t345 * t302;
t281 = t296 * pkin(5) - t297 * pkin(13) + t283;
t280 = t302 * pkin(13) + t359;
t279 = -t302 * pkin(5) - t350;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t325 * t340 - t332 * t361) * qJD(1) (-t326 * t340 + t332 * t364) * qJD(1) (-t325 * t333 + t326 * t337) * t358, t326 ^ 2 / 0.2e1 + t325 ^ 2 / 0.2e1 + t332 ^ 2 / 0.2e1, t316 ^ 2 / 0.2e1, -t316 * t315, -t323 * t316, t323 * t315, t323 ^ 2 / 0.2e1, t305 * t315 - t323 * t351, t305 * t316 + t323 * t354, t304 ^ 2 / 0.2e1, -t303 * t304, -t307 * t304, t307 * t303, t307 ^ 2 / 0.2e1, t287 * t303 - t371 * t307, t287 * t304 + t307 * t355, t297 ^ 2 / 0.2e1, -t297 * t296, t297 * t302, -t296 * t302, t302 ^ 2 / 0.2e1, t283 * t296 + t302 * t350, t283 * t297 - t302 * t359, t289 ^ 2 / 0.2e1, -t289 * t288, t289 * t295, -t288 * t295, t295 ^ 2 / 0.2e1 (-t341 * t280 + t345 * t281) * t295 + t279 * t288 -(t345 * t280 + t341 * t281) * t295 + t279 * t289;];
T_reg  = t1;
