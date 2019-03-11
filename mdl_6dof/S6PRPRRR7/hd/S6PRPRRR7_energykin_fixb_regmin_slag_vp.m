% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:53
% EndTime: 2019-03-08 20:55:53
% DurationCPUTime: 0.23s
% Computational Cost: add. (650->61), mult. (1993->133), div. (0->0), fcn. (1709->16), ass. (0->56)
t297 = cos(qJ(2));
t308 = qJD(1) * sin(pkin(6));
t276 = qJD(2) * pkin(2) + t297 * t308;
t284 = sin(pkin(7));
t288 = cos(pkin(7));
t307 = qJD(1) * cos(pkin(6));
t317 = t276 * t288 + t284 * t307;
t293 = sin(qJ(2));
t306 = qJD(2) * t284;
t274 = qJ(3) * t306 + t293 * t308;
t282 = sin(pkin(14));
t286 = cos(pkin(14));
t261 = t286 * t274 + t317 * t282;
t283 = sin(pkin(8));
t287 = cos(pkin(8));
t311 = t284 * t287;
t253 = (t283 * t288 + t286 * t311) * qJD(2) * pkin(10) + t261;
t260 = -t282 * t274 + t317 * t286;
t315 = pkin(10) * t282;
t254 = (pkin(3) * t288 - t311 * t315) * qJD(2) + t260;
t304 = t288 * t307 + qJD(3);
t262 = (-t276 + (-pkin(3) * t286 - t283 * t315) * qJD(2)) * t284 + t304;
t292 = sin(qJ(4));
t296 = cos(qJ(4));
t316 = -t292 * t253 + (t254 * t287 + t262 * t283) * t296;
t267 = -t284 * t276 + t304;
t314 = t267 * t284;
t312 = t283 * t292;
t310 = t287 * t292;
t301 = t286 * t306;
t305 = qJD(2) * t288;
t269 = t283 * t301 - t287 * t305 - qJD(4);
t303 = t296 * t253 + t254 * t310 + t262 * t312;
t244 = -t269 * pkin(11) + t303;
t247 = -t283 * t254 + t287 * t262;
t265 = t292 * t282 * t306 + (-t283 * t305 - t287 * t301) * t296;
t266 = (t288 * t312 + (t282 * t296 + t286 * t310) * t284) * qJD(2);
t246 = t265 * pkin(4) - t266 * pkin(11) + t247;
t291 = sin(qJ(5));
t295 = cos(qJ(5));
t309 = t295 * t244 + t291 * t246;
t300 = qJD(2) * t308;
t256 = t291 * t266 + t295 * t269;
t299 = -t291 * t244 + t295 * t246;
t243 = t269 * pkin(4) - t316;
t294 = cos(qJ(6));
t290 = sin(qJ(6));
t264 = qJD(5) + t265;
t257 = t295 * t266 - t291 * t269;
t255 = qJD(6) + t256;
t249 = t294 * t257 + t290 * t264;
t248 = t290 * t257 - t294 * t264;
t241 = t256 * pkin(5) - t257 * pkin(12) + t243;
t240 = t264 * pkin(12) + t309;
t239 = -t264 * pkin(5) - t299;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t297 * t300, -t293 * t300 (t260 * t288 - t286 * t314) * qJD(2) (-t261 * t288 + t282 * t314) * qJD(2) (-t260 * t282 + t261 * t286) * t306, t261 ^ 2 / 0.2e1 + t260 ^ 2 / 0.2e1 + t267 ^ 2 / 0.2e1, t266 ^ 2 / 0.2e1, -t265 * t266, -t269 * t266, t269 * t265, t269 ^ 2 / 0.2e1, t247 * t265 - t316 * t269, t247 * t266 + t303 * t269, t257 ^ 2 / 0.2e1, -t257 * t256, t257 * t264, -t256 * t264, t264 ^ 2 / 0.2e1, t243 * t256 + t299 * t264, t243 * t257 - t309 * t264, t249 ^ 2 / 0.2e1, -t249 * t248, t249 * t255, -t248 * t255, t255 ^ 2 / 0.2e1 (-t290 * t240 + t294 * t241) * t255 + t239 * t248 -(t294 * t240 + t290 * t241) * t255 + t239 * t249;];
T_reg  = t1;
