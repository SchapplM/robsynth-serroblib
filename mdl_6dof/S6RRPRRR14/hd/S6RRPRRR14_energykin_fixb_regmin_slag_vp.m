% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:08:37
% EndTime: 2019-03-09 15:08:38
% DurationCPUTime: 0.42s
% Computational Cost: add. (2027->72), mult. (5772->148), div. (0->0), fcn. (4990->16), ass. (0->66)
t363 = cos(pkin(6)) * qJD(1);
t333 = qJD(2) + t363;
t335 = sin(pkin(14));
t340 = cos(pkin(7));
t349 = cos(qJ(2));
t338 = sin(pkin(6));
t364 = qJD(1) * t338;
t358 = t349 * t364;
t356 = t340 * t358;
t337 = sin(pkin(7));
t373 = cos(pkin(14));
t357 = t337 * t373;
t345 = sin(qJ(2));
t359 = t345 * t364;
t319 = -t333 * t357 + t335 * t359 - t356 * t373;
t336 = sin(pkin(8));
t339 = cos(pkin(8));
t353 = -t340 * t333 + t337 * t358;
t377 = t339 * t319 + t336 * t353;
t362 = pkin(1) * t363;
t365 = pkin(10) * t358 + t345 * t362;
t318 = (t333 * t337 + t356) * qJ(3) + t365;
t332 = t349 * t362;
t321 = t333 * pkin(2) + t332 + (-qJ(3) * t340 - pkin(10)) * t359;
t326 = (-qJ(3) * t337 * t345 - pkin(2) * t349 - pkin(1)) * t364;
t304 = t340 * t373 * t321 - t335 * t318 + t326 * t357;
t368 = t335 * t340;
t369 = t335 * t337;
t320 = t333 * t369 + (t345 * t373 + t349 * t368) * t364;
t374 = pkin(11) * t320;
t300 = -pkin(3) * t353 - t339 * t374 + t304;
t310 = -t337 * t321 + t340 * t326 + qJD(3);
t303 = t319 * pkin(3) - t336 * t374 + t310;
t376 = t300 * t339 + t303 * t336;
t305 = t373 * t318 + t321 * t368 + t326 * t369;
t297 = -pkin(11) * t377 + t305;
t344 = sin(qJ(4));
t348 = cos(qJ(4));
t375 = -t344 * t297 + t348 * t376;
t350 = qJD(1) ^ 2;
t370 = t338 ^ 2 * t350;
t312 = -t336 * t319 + t339 * t353 - qJD(4);
t360 = t348 * t297 + t376 * t344;
t287 = -t312 * pkin(12) + t360;
t290 = -t336 * t300 + t339 * t303;
t308 = t344 * t320 + t377 * t348;
t309 = t348 * t320 - t344 * t377;
t289 = t308 * pkin(4) - t309 * pkin(12) + t290;
t343 = sin(qJ(5));
t347 = cos(qJ(5));
t366 = t347 * t287 + t343 * t289;
t361 = t349 * t370;
t298 = t343 * t309 + t347 * t312;
t355 = -t343 * t287 + t347 * t289;
t286 = t312 * pkin(4) - t375;
t346 = cos(qJ(6));
t342 = sin(qJ(6));
t307 = qJD(5) + t308;
t299 = t347 * t309 - t343 * t312;
t296 = qJD(6) + t298;
t292 = t346 * t299 + t342 * t307;
t291 = t342 * t299 - t346 * t307;
t284 = t298 * pkin(5) - t299 * pkin(13) + t286;
t283 = t307 * pkin(13) + t366;
t282 = -t307 * pkin(5) - t355;
t1 = [t350 / 0.2e1, 0, 0, t345 ^ 2 * t370 / 0.2e1, t345 * t361, t333 * t359, t333 * t358, t333 ^ 2 / 0.2e1 (-pkin(10) * t359 + t332) * t333 + pkin(1) * t361, -pkin(1) * t345 * t370 - t333 * t365, -t304 * t353 + t310 * t319, t305 * t353 + t310 * t320, -t304 * t320 - t305 * t319, t305 ^ 2 / 0.2e1 + t304 ^ 2 / 0.2e1 + t310 ^ 2 / 0.2e1, t309 ^ 2 / 0.2e1, -t308 * t309, -t312 * t309, t312 * t308, t312 ^ 2 / 0.2e1, t290 * t308 - t375 * t312, t290 * t309 + t312 * t360, t299 ^ 2 / 0.2e1, -t299 * t298, t299 * t307, -t298 * t307, t307 ^ 2 / 0.2e1, t286 * t298 + t307 * t355, t286 * t299 - t307 * t366, t292 ^ 2 / 0.2e1, -t292 * t291, t292 * t296, -t291 * t296, t296 ^ 2 / 0.2e1 (-t342 * t283 + t346 * t284) * t296 + t282 * t291 -(t346 * t283 + t342 * t284) * t296 + t282 * t292;];
T_reg  = t1;
