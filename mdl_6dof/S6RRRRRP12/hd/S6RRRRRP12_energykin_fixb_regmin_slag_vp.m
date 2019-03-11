% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:36
% EndTime: 2019-03-10 03:20:37
% DurationCPUTime: 0.28s
% Computational Cost: add. (1327->62), mult. (3486->128), div. (0->0), fcn. (2905->12), ass. (0->54)
t286 = cos(pkin(6)) * qJD(1);
t262 = qJD(2) + t286;
t264 = sin(pkin(7));
t266 = cos(pkin(7));
t275 = cos(qJ(2));
t265 = sin(pkin(6));
t287 = qJD(1) * t265;
t281 = t275 * t287;
t296 = t262 * t264 + t266 * t281;
t271 = sin(qJ(2));
t285 = pkin(1) * t286;
t288 = pkin(9) * t281 + t271 * t285;
t245 = t296 * pkin(10) + t288;
t261 = t275 * t285;
t282 = t271 * t287;
t247 = t262 * pkin(2) + t261 + (-pkin(10) * t266 - pkin(9)) * t282;
t253 = (-pkin(10) * t264 * t271 - pkin(2) * t275 - pkin(1)) * t287;
t270 = sin(qJ(3));
t274 = cos(qJ(3));
t295 = -t270 * t245 + (t247 * t266 + t253 * t264) * t274;
t276 = qJD(1) ^ 2;
t293 = t265 ^ 2 * t276;
t292 = t264 * t270;
t291 = t266 * t270;
t248 = t270 * t282 - t296 * t274;
t246 = qJD(4) + t248;
t237 = -t264 * t247 + t266 * t253;
t249 = t262 * t292 + (t271 * t274 + t275 * t291) * t287;
t230 = t248 * pkin(3) - t249 * pkin(11) + t237;
t254 = -t266 * t262 + t264 * t281 - qJD(3);
t283 = t274 * t245 + t247 * t291 + t253 * t292;
t234 = -t254 * pkin(11) + t283;
t269 = sin(qJ(4));
t273 = cos(qJ(4));
t289 = t269 * t230 + t273 * t234;
t226 = t246 * pkin(12) + t289;
t233 = t254 * pkin(3) - t295;
t239 = t269 * t249 + t273 * t254;
t240 = t273 * t249 - t269 * t254;
t228 = t239 * pkin(4) - t240 * pkin(12) + t233;
t268 = sin(qJ(5));
t272 = cos(qJ(5));
t290 = t272 * t226 + t268 * t228;
t284 = t275 * t293;
t280 = t273 * t230 - t269 * t234;
t278 = -t268 * t226 + t272 * t228;
t225 = -t246 * pkin(4) - t280;
t238 = qJD(5) + t239;
t236 = t272 * t240 + t268 * t246;
t235 = t268 * t240 - t272 * t246;
t223 = t235 * pkin(5) - t236 * qJ(6) + t225;
t222 = t238 * qJ(6) + t290;
t221 = -t238 * pkin(5) + qJD(6) - t278;
t1 = [t276 / 0.2e1, 0, 0, t271 ^ 2 * t293 / 0.2e1, t271 * t284, t262 * t282, t262 * t281, t262 ^ 2 / 0.2e1, pkin(1) * t284 + (-pkin(9) * t282 + t261) * t262, -pkin(1) * t271 * t293 - t288 * t262, t249 ^ 2 / 0.2e1, -t249 * t248, -t249 * t254, t248 * t254, t254 ^ 2 / 0.2e1, t237 * t248 - t295 * t254, t237 * t249 + t283 * t254, t240 ^ 2 / 0.2e1, -t240 * t239, t240 * t246, -t239 * t246, t246 ^ 2 / 0.2e1, t233 * t239 + t280 * t246, t233 * t240 - t289 * t246, t236 ^ 2 / 0.2e1, -t236 * t235, t236 * t238, -t235 * t238, t238 ^ 2 / 0.2e1, t225 * t235 + t238 * t278, t225 * t236 - t290 * t238, -t221 * t238 + t223 * t235, t221 * t236 - t222 * t235, t222 * t238 - t223 * t236, t222 ^ 2 / 0.2e1 + t223 ^ 2 / 0.2e1 + t221 ^ 2 / 0.2e1;];
T_reg  = t1;
