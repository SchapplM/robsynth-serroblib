% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:08
% EndTime: 2019-03-09 07:34:08
% DurationCPUTime: 0.25s
% Computational Cost: add. (1073->70), mult. (3436->143), div. (0->0), fcn. (2931->14), ass. (0->59)
t278 = sin(pkin(13));
t281 = cos(pkin(13));
t280 = sin(pkin(6));
t299 = qJD(1) * t280;
t294 = t281 * t299;
t283 = cos(pkin(6));
t298 = qJD(1) * t283;
t297 = pkin(1) * t298;
t270 = qJ(2) * t294 + t278 * t297;
t279 = sin(pkin(7));
t282 = cos(pkin(7));
t303 = t280 * t281;
t257 = (t279 * t283 + t282 * t303) * qJD(1) * pkin(9) + t270;
t275 = t281 * t297;
t305 = t278 * t280;
t260 = t275 + (pkin(2) * t283 + (-pkin(9) * t282 - qJ(2)) * t305) * qJD(1);
t265 = qJD(2) + (-pkin(9) * t278 * t279 - pkin(2) * t281 - pkin(1)) * t299;
t287 = sin(qJ(3));
t290 = cos(qJ(3));
t307 = -t287 * t257 + (t260 * t282 + t265 * t279) * t290;
t306 = cos(qJ(4));
t304 = t279 * t287;
t302 = t282 * t287;
t262 = (t283 * t304 + (t278 * t290 + t281 * t302) * t280) * qJD(1);
t267 = t279 * t294 - t282 * t298 - qJD(3);
t286 = sin(qJ(4));
t252 = t262 * t306 - t286 * t267;
t295 = t278 * t299;
t261 = t287 * t295 + (-t279 * t298 - t282 * t294) * t290;
t259 = qJD(4) + t261;
t249 = -t279 * t260 + t282 * t265;
t242 = t261 * pkin(3) - t262 * pkin(10) + t249;
t296 = t290 * t257 + t260 * t302 + t265 * t304;
t245 = -t267 * pkin(10) + t296;
t293 = t242 * t306 - t286 * t245;
t234 = t259 * pkin(4) - t252 * pkin(11) + t293;
t251 = t286 * t262 + t267 * t306;
t300 = t286 * t242 + t245 * t306;
t236 = -t251 * pkin(11) + t300;
t285 = sin(qJ(5));
t289 = cos(qJ(5));
t301 = t285 * t234 + t289 * t236;
t247 = t289 * t251 + t285 * t252;
t292 = t289 * t234 - t285 * t236;
t244 = t267 * pkin(3) - t307;
t237 = t251 * pkin(4) + t244;
t288 = cos(qJ(6));
t284 = sin(qJ(6));
t276 = -pkin(1) * t299 + qJD(2);
t269 = -qJ(2) * t295 + t275;
t258 = qJD(5) + t259;
t248 = -t285 * t251 + t289 * t252;
t246 = qJD(6) + t247;
t239 = t288 * t248 + t284 * t258;
t238 = t284 * t248 - t288 * t258;
t232 = t247 * pkin(5) - t248 * pkin(12) + t237;
t231 = t258 * pkin(12) + t301;
t230 = -t258 * pkin(5) - t292;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t269 * t283 - t276 * t303) * qJD(1) (-t270 * t283 + t276 * t305) * qJD(1) (-t269 * t278 + t270 * t281) * t299, t270 ^ 2 / 0.2e1 + t269 ^ 2 / 0.2e1 + t276 ^ 2 / 0.2e1, t262 ^ 2 / 0.2e1, -t262 * t261, -t262 * t267, t261 * t267, t267 ^ 2 / 0.2e1, t249 * t261 - t307 * t267, t249 * t262 + t267 * t296, t252 ^ 2 / 0.2e1, -t252 * t251, t252 * t259, -t251 * t259, t259 ^ 2 / 0.2e1, t244 * t251 + t259 * t293, t244 * t252 - t259 * t300, t248 ^ 2 / 0.2e1, -t248 * t247, t248 * t258, -t247 * t258, t258 ^ 2 / 0.2e1, t237 * t247 + t258 * t292, t237 * t248 - t258 * t301, t239 ^ 2 / 0.2e1, -t239 * t238, t239 * t246, -t238 * t246, t246 ^ 2 / 0.2e1 (-t284 * t231 + t288 * t232) * t246 + t230 * t238 -(t288 * t231 + t284 * t232) * t246 + t230 * t239;];
T_reg  = t1;
