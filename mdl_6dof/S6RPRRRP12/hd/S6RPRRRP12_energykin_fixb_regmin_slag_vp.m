% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:38
% EndTime: 2019-03-09 06:50:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (1097->67), mult. (3502->136), div. (0->0), fcn. (2916->12), ass. (0->54)
t265 = sin(pkin(12));
t268 = cos(pkin(12));
t267 = sin(pkin(6));
t285 = qJD(1) * t267;
t280 = t268 * t285;
t270 = cos(pkin(6));
t284 = qJD(1) * t270;
t283 = pkin(1) * t284;
t257 = qJ(2) * t280 + t265 * t283;
t266 = sin(pkin(7));
t269 = cos(pkin(7));
t289 = t267 * t268;
t245 = (t266 * t270 + t269 * t289) * qJD(1) * pkin(9) + t257;
t262 = t268 * t283;
t291 = t265 * t267;
t247 = t262 + (pkin(2) * t270 + (-pkin(9) * t269 - qJ(2)) * t291) * qJD(1);
t252 = qJD(2) + (-pkin(9) * t265 * t266 - pkin(2) * t268 - pkin(1)) * t285;
t273 = sin(qJ(3));
t276 = cos(qJ(3));
t292 = -t273 * t245 + (t247 * t269 + t252 * t266) * t276;
t290 = t266 * t273;
t288 = t269 * t273;
t281 = t265 * t285;
t248 = t273 * t281 + (-t266 * t284 - t269 * t280) * t276;
t246 = qJD(4) + t248;
t237 = -t266 * t247 + t269 * t252;
t249 = (t270 * t290 + (t265 * t276 + t268 * t288) * t267) * qJD(1);
t230 = t248 * pkin(3) - t249 * pkin(10) + t237;
t254 = t266 * t280 - t269 * t284 - qJD(3);
t282 = t276 * t245 + t247 * t288 + t252 * t290;
t234 = -t254 * pkin(10) + t282;
t272 = sin(qJ(4));
t275 = cos(qJ(4));
t286 = t272 * t230 + t275 * t234;
t226 = t246 * pkin(11) + t286;
t233 = t254 * pkin(3) - t292;
t239 = t272 * t249 + t275 * t254;
t240 = t275 * t249 - t272 * t254;
t228 = t239 * pkin(4) - t240 * pkin(11) + t233;
t271 = sin(qJ(5));
t274 = cos(qJ(5));
t287 = t274 * t226 + t271 * t228;
t279 = t275 * t230 - t272 * t234;
t278 = -t271 * t226 + t274 * t228;
t225 = -t246 * pkin(4) - t279;
t263 = -pkin(1) * t285 + qJD(2);
t256 = -qJ(2) * t281 + t262;
t238 = qJD(5) + t239;
t236 = t274 * t240 + t271 * t246;
t235 = t271 * t240 - t274 * t246;
t223 = t235 * pkin(5) - t236 * qJ(6) + t225;
t222 = t238 * qJ(6) + t287;
t221 = -t238 * pkin(5) + qJD(6) - t278;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t256 * t270 - t263 * t289) * qJD(1) (-t257 * t270 + t263 * t291) * qJD(1) (-t256 * t265 + t257 * t268) * t285, t257 ^ 2 / 0.2e1 + t256 ^ 2 / 0.2e1 + t263 ^ 2 / 0.2e1, t249 ^ 2 / 0.2e1, -t249 * t248, -t249 * t254, t248 * t254, t254 ^ 2 / 0.2e1, t237 * t248 - t292 * t254, t237 * t249 + t282 * t254, t240 ^ 2 / 0.2e1, -t240 * t239, t240 * t246, -t239 * t246, t246 ^ 2 / 0.2e1, t233 * t239 + t279 * t246, t233 * t240 - t286 * t246, t236 ^ 2 / 0.2e1, -t236 * t235, t236 * t238, -t235 * t238, t238 ^ 2 / 0.2e1, t225 * t235 + t278 * t238, t225 * t236 - t287 * t238, -t221 * t238 + t223 * t235, t221 * t236 - t222 * t235, t222 * t238 - t223 * t236, t222 ^ 2 / 0.2e1 + t223 ^ 2 / 0.2e1 + t221 ^ 2 / 0.2e1;];
T_reg  = t1;
