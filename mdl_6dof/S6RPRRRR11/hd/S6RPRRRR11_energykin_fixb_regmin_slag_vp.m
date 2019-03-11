% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR11
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:46:03
% EndTime: 2019-03-09 07:46:04
% DurationCPUTime: 0.27s
% Computational Cost: add. (1035->70), mult. (3313->142), div. (0->0), fcn. (2826->14), ass. (0->60)
t275 = sin(pkin(13));
t303 = cos(pkin(6));
t291 = qJD(1) * t303;
t288 = pkin(1) * t291;
t278 = cos(pkin(13));
t277 = sin(pkin(6));
t296 = qJD(1) * t277;
t293 = t278 * t296;
t267 = qJ(2) * t293 + t275 * t288;
t279 = cos(pkin(7));
t276 = sin(pkin(7));
t292 = t303 * t276;
t300 = t277 * t278;
t255 = (t279 * t300 + t292) * qJD(1) * pkin(9) + t267;
t272 = t278 * t288;
t301 = t275 * t277;
t257 = t272 + (t303 * pkin(2) + (-pkin(9) * t279 - qJ(2)) * t301) * qJD(1);
t283 = sin(qJ(3));
t286 = cos(qJ(3));
t262 = qJD(2) + (-pkin(9) * t275 * t276 - pkin(2) * t278 - pkin(1)) * t296;
t302 = t262 * t276;
t305 = -t283 * t255 + (t257 * t279 + t302) * t286;
t304 = cos(qJ(5));
t299 = t279 * t283;
t294 = t275 * t296;
t258 = t283 * t294 + (-t276 * t291 - t279 * t293) * t286;
t256 = qJD(4) + t258;
t246 = -t276 * t257 + t279 * t262;
t259 = (t283 * t292 + (t275 * t286 + t278 * t299) * t277) * qJD(1);
t238 = t258 * pkin(3) - t259 * pkin(10) + t246;
t264 = t276 * t293 - t279 * t291 - qJD(3);
t295 = t286 * t255 + t257 * t299 + t283 * t302;
t242 = -t264 * pkin(10) + t295;
t282 = sin(qJ(4));
t285 = cos(qJ(4));
t297 = t282 * t238 + t285 * t242;
t231 = t256 * pkin(11) + t297;
t241 = t264 * pkin(3) - t305;
t249 = t282 * t259 + t285 * t264;
t250 = t285 * t259 - t282 * t264;
t234 = t249 * pkin(4) - t250 * pkin(11) + t241;
t281 = sin(qJ(5));
t298 = t304 * t231 + t281 * t234;
t290 = -t281 * t231 + t304 * t234;
t289 = t285 * t238 - t282 * t242;
t248 = qJD(5) + t249;
t230 = -t256 * pkin(4) - t289;
t284 = cos(qJ(6));
t280 = sin(qJ(6));
t273 = -pkin(1) * t296 + qJD(2);
t266 = -qJ(2) * t294 + t272;
t247 = qJD(6) + t248;
t245 = t304 * t250 + t281 * t256;
t244 = t281 * t250 - t304 * t256;
t236 = -t280 * t244 + t284 * t245;
t235 = t284 * t244 + t280 * t245;
t228 = t244 * pkin(5) + t230;
t227 = -t244 * pkin(12) + t298;
t226 = t248 * pkin(5) - t245 * pkin(12) + t290;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t303 * t266 - t273 * t300) * qJD(1) (-t303 * t267 + t273 * t301) * qJD(1) (-t266 * t275 + t267 * t278) * t296, t267 ^ 2 / 0.2e1 + t266 ^ 2 / 0.2e1 + t273 ^ 2 / 0.2e1, t259 ^ 2 / 0.2e1, -t259 * t258, -t259 * t264, t258 * t264, t264 ^ 2 / 0.2e1, t246 * t258 - t264 * t305, t246 * t259 + t295 * t264, t250 ^ 2 / 0.2e1, -t250 * t249, t250 * t256, -t249 * t256, t256 ^ 2 / 0.2e1, t241 * t249 + t289 * t256, t241 * t250 - t297 * t256, t245 ^ 2 / 0.2e1, -t245 * t244, t245 * t248, -t244 * t248, t248 ^ 2 / 0.2e1, t230 * t244 + t290 * t248, t230 * t245 - t298 * t248, t236 ^ 2 / 0.2e1, -t236 * t235, t236 * t247, -t235 * t247, t247 ^ 2 / 0.2e1 (t284 * t226 - t280 * t227) * t247 + t228 * t235 -(t280 * t226 + t284 * t227) * t247 + t228 * t236;];
T_reg  = t1;
