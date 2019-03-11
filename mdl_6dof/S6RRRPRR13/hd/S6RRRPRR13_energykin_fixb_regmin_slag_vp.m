% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:08
% EndTime: 2019-03-09 20:05:08
% DurationCPUTime: 0.30s
% Computational Cost: add. (1412->66), mult. (3753->136), div. (0->0), fcn. (3156->14), ass. (0->58)
t299 = cos(pkin(6)) * qJD(1);
t275 = qJD(2) + t299;
t278 = sin(pkin(7));
t280 = cos(pkin(7));
t289 = cos(qJ(2));
t279 = sin(pkin(6));
t300 = qJD(1) * t279;
t294 = t289 * t300;
t309 = t275 * t278 + t280 * t294;
t285 = sin(qJ(2));
t298 = pkin(1) * t299;
t301 = pkin(9) * t294 + t285 * t298;
t258 = t309 * pkin(10) + t301;
t274 = t289 * t298;
t295 = t285 * t300;
t260 = t275 * pkin(2) + t274 + (-pkin(10) * t280 - pkin(9)) * t295;
t266 = (-pkin(10) * t278 * t285 - pkin(2) * t289 - pkin(1)) * t300;
t284 = sin(qJ(3));
t288 = cos(qJ(3));
t308 = -t284 * t258 + (t260 * t280 + t266 * t278) * t288;
t307 = cos(pkin(13));
t290 = qJD(1) ^ 2;
t305 = t279 ^ 2 * t290;
t304 = t278 * t284;
t303 = t280 * t284;
t250 = -t278 * t260 + t280 * t266;
t261 = t284 * t295 - t309 * t288;
t262 = t275 * t304 + (t285 * t288 + t289 * t303) * t300;
t243 = t261 * pkin(3) - t262 * qJ(4) + t250;
t267 = -t280 * t275 + t278 * t294 - qJD(3);
t296 = t288 * t258 + t260 * t303 + t266 * t304;
t246 = -t267 * qJ(4) + t296;
t277 = sin(pkin(13));
t236 = t307 * t243 - t277 * t246;
t253 = t307 * t262 - t277 * t267;
t233 = t261 * pkin(4) - t253 * pkin(11) + t236;
t237 = t277 * t243 + t307 * t246;
t252 = t277 * t262 + t307 * t267;
t235 = -t252 * pkin(11) + t237;
t283 = sin(qJ(5));
t287 = cos(qJ(5));
t302 = t283 * t233 + t287 * t235;
t297 = t289 * t305;
t248 = t287 * t252 + t283 * t253;
t292 = t287 * t233 - t283 * t235;
t245 = t267 * pkin(3) + qJD(4) - t308;
t238 = t252 * pkin(4) + t245;
t286 = cos(qJ(6));
t282 = sin(qJ(6));
t259 = qJD(5) + t261;
t249 = -t283 * t252 + t287 * t253;
t247 = qJD(6) + t248;
t240 = t286 * t249 + t282 * t259;
t239 = t282 * t249 - t286 * t259;
t231 = t248 * pkin(5) - t249 * pkin(12) + t238;
t230 = t259 * pkin(12) + t302;
t229 = -t259 * pkin(5) - t292;
t1 = [t290 / 0.2e1, 0, 0, t285 ^ 2 * t305 / 0.2e1, t285 * t297, t275 * t295, t275 * t294, t275 ^ 2 / 0.2e1, pkin(1) * t297 + (-pkin(9) * t295 + t274) * t275, -pkin(1) * t285 * t305 - t301 * t275, t262 ^ 2 / 0.2e1, -t262 * t261, -t262 * t267, t261 * t267, t267 ^ 2 / 0.2e1, t250 * t261 - t267 * t308, t250 * t262 + t296 * t267, t236 * t261 + t245 * t252, -t237 * t261 + t245 * t253, -t236 * t253 - t237 * t252, t237 ^ 2 / 0.2e1 + t236 ^ 2 / 0.2e1 + t245 ^ 2 / 0.2e1, t249 ^ 2 / 0.2e1, -t249 * t248, t249 * t259, -t248 * t259, t259 ^ 2 / 0.2e1, t238 * t248 + t292 * t259, t238 * t249 - t302 * t259, t240 ^ 2 / 0.2e1, -t240 * t239, t240 * t247, -t239 * t247, t247 ^ 2 / 0.2e1 (-t282 * t230 + t286 * t231) * t247 + t229 * t239 -(t286 * t230 + t282 * t231) * t247 + t229 * t240;];
T_reg  = t1;
