% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:43
% EndTime: 2019-03-10 05:06:43
% DurationCPUTime: 0.29s
% Computational Cost: add. (1315->65), mult. (3420->135), div. (0->0), fcn. (2920->14), ass. (0->59)
t299 = cos(pkin(6)) * qJD(1);
t274 = qJD(2) + t299;
t276 = sin(pkin(7));
t278 = cos(pkin(7));
t288 = cos(qJ(2));
t277 = sin(pkin(6));
t300 = qJD(1) * t277;
t294 = t288 * t300;
t310 = t274 * t276 + t278 * t294;
t284 = sin(qJ(2));
t298 = pkin(1) * t299;
t301 = pkin(9) * t294 + t284 * t298;
t256 = t310 * pkin(10) + t301;
t273 = t288 * t298;
t295 = t284 * t300;
t259 = t274 * pkin(2) + t273 + (-pkin(10) * t278 - pkin(9)) * t295;
t265 = (-pkin(10) * t276 * t284 - pkin(2) * t288 - pkin(1)) * t300;
t283 = sin(qJ(3));
t287 = cos(qJ(3));
t309 = -t283 * t256 + (t259 * t278 + t265 * t276) * t287;
t308 = cos(qJ(4));
t289 = qJD(1) ^ 2;
t306 = t277 ^ 2 * t289;
t305 = t276 * t283;
t304 = t278 * t283;
t261 = t274 * t305 + (t284 * t287 + t288 * t304) * t300;
t266 = -t278 * t274 + t276 * t294 - qJD(3);
t282 = sin(qJ(4));
t251 = t308 * t261 - t282 * t266;
t260 = t283 * t295 - t310 * t287;
t258 = qJD(4) + t260;
t248 = -t276 * t259 + t278 * t265;
t241 = t260 * pkin(3) - t261 * pkin(11) + t248;
t296 = t287 * t256 + t259 * t304 + t265 * t305;
t244 = -t266 * pkin(11) + t296;
t293 = t308 * t241 - t282 * t244;
t233 = t258 * pkin(4) - t251 * pkin(12) + t293;
t250 = t282 * t261 + t308 * t266;
t302 = t282 * t241 + t308 * t244;
t235 = -t250 * pkin(12) + t302;
t281 = sin(qJ(5));
t286 = cos(qJ(5));
t303 = t281 * t233 + t286 * t235;
t297 = t288 * t306;
t246 = t286 * t250 + t281 * t251;
t291 = t286 * t233 - t281 * t235;
t243 = t266 * pkin(3) - t309;
t236 = t250 * pkin(4) + t243;
t285 = cos(qJ(6));
t280 = sin(qJ(6));
t257 = qJD(5) + t258;
t247 = -t281 * t250 + t286 * t251;
t245 = qJD(6) + t246;
t238 = t285 * t247 + t280 * t257;
t237 = t280 * t247 - t285 * t257;
t231 = t246 * pkin(5) - t247 * pkin(13) + t236;
t230 = t257 * pkin(13) + t303;
t229 = -t257 * pkin(5) - t291;
t1 = [t289 / 0.2e1, 0, 0, t284 ^ 2 * t306 / 0.2e1, t284 * t297, t274 * t295, t274 * t294, t274 ^ 2 / 0.2e1, pkin(1) * t297 + (-pkin(9) * t295 + t273) * t274, -pkin(1) * t284 * t306 - t301 * t274, t261 ^ 2 / 0.2e1, -t261 * t260, -t261 * t266, t260 * t266, t266 ^ 2 / 0.2e1, t248 * t260 - t266 * t309, t248 * t261 + t296 * t266, t251 ^ 2 / 0.2e1, -t251 * t250, t251 * t258, -t250 * t258, t258 ^ 2 / 0.2e1, t243 * t250 + t293 * t258, t243 * t251 - t302 * t258, t247 ^ 2 / 0.2e1, -t247 * t246, t247 * t257, -t246 * t257, t257 ^ 2 / 0.2e1, t236 * t246 + t291 * t257, t236 * t247 - t303 * t257, t238 ^ 2 / 0.2e1, -t238 * t237, t238 * t245, -t237 * t245, t245 ^ 2 / 0.2e1 (-t280 * t230 + t285 * t231) * t245 + t229 * t237 -(t285 * t230 + t280 * t231) * t245 + t229 * t238;];
T_reg  = t1;
