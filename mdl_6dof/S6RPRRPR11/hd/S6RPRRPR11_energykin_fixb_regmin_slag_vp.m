% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:31
% EndTime: 2019-03-09 05:44:31
% DurationCPUTime: 0.26s
% Computational Cost: add. (1257->71), mult. (4023->144), div. (0->0), fcn. (3369->14), ass. (0->58)
t276 = sin(pkin(12));
t279 = cos(pkin(12));
t278 = sin(pkin(6));
t295 = qJD(1) * t278;
t290 = t279 * t295;
t281 = cos(pkin(6));
t294 = qJD(1) * t281;
t293 = pkin(1) * t294;
t267 = qJ(2) * t290 + t276 * t293;
t277 = sin(pkin(7));
t280 = cos(pkin(7));
t298 = t278 * t279;
t255 = (t277 * t281 + t280 * t298) * qJD(1) * pkin(9) + t267;
t272 = t279 * t293;
t300 = t276 * t278;
t257 = t272 + (pkin(2) * t281 + (-pkin(9) * t280 - qJ(2)) * t300) * qJD(1);
t262 = qJD(2) + (-pkin(9) * t276 * t277 - pkin(2) * t279 - pkin(1)) * t295;
t284 = sin(qJ(3));
t287 = cos(qJ(3));
t302 = -t284 * t255 + (t257 * t280 + t262 * t277) * t287;
t301 = cos(pkin(13));
t299 = t277 * t284;
t297 = t280 * t284;
t291 = t276 * t295;
t258 = t284 * t291 + (-t277 * t294 - t280 * t290) * t287;
t256 = qJD(4) + t258;
t247 = -t277 * t257 + t280 * t262;
t259 = (t281 * t299 + (t276 * t287 + t279 * t297) * t278) * qJD(1);
t239 = t258 * pkin(3) - t259 * pkin(10) + t247;
t264 = t277 * t290 - t280 * t294 - qJD(3);
t292 = t287 * t255 + t257 * t297 + t262 * t299;
t243 = -t264 * pkin(10) + t292;
t283 = sin(qJ(4));
t286 = cos(qJ(4));
t296 = t283 * t239 + t286 * t243;
t232 = t256 * qJ(5) + t296;
t242 = t264 * pkin(3) - t302;
t249 = t283 * t259 + t286 * t264;
t250 = t286 * t259 - t283 * t264;
t235 = t249 * pkin(4) - t250 * qJ(5) + t242;
t275 = sin(pkin(13));
t228 = t301 * t232 + t275 * t235;
t227 = -t275 * t232 + t301 * t235;
t289 = t286 * t239 - t283 * t243;
t231 = -t256 * pkin(4) + qJD(5) - t289;
t285 = cos(qJ(6));
t282 = sin(qJ(6));
t273 = -pkin(1) * t295 + qJD(2);
t266 = -qJ(2) * t291 + t272;
t248 = qJD(6) + t249;
t246 = t301 * t250 + t275 * t256;
t245 = t275 * t250 - t301 * t256;
t237 = -t282 * t245 + t285 * t246;
t236 = t285 * t245 + t282 * t246;
t229 = t245 * pkin(5) + t231;
t226 = -t245 * pkin(11) + t228;
t225 = t249 * pkin(5) - t246 * pkin(11) + t227;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t266 * t281 - t273 * t298) * qJD(1) (-t267 * t281 + t273 * t300) * qJD(1) (-t266 * t276 + t267 * t279) * t295, t267 ^ 2 / 0.2e1 + t266 ^ 2 / 0.2e1 + t273 ^ 2 / 0.2e1, t259 ^ 2 / 0.2e1, -t259 * t258, -t259 * t264, t258 * t264, t264 ^ 2 / 0.2e1, t247 * t258 - t302 * t264, t247 * t259 + t292 * t264, t250 ^ 2 / 0.2e1, -t250 * t249, t250 * t256, -t249 * t256, t256 ^ 2 / 0.2e1, t242 * t249 + t289 * t256, t242 * t250 - t296 * t256, t227 * t249 + t231 * t245, -t228 * t249 + t231 * t246, -t227 * t246 - t228 * t245, t228 ^ 2 / 0.2e1 + t227 ^ 2 / 0.2e1 + t231 ^ 2 / 0.2e1, t237 ^ 2 / 0.2e1, -t237 * t236, t237 * t248, -t236 * t248, t248 ^ 2 / 0.2e1 (t285 * t225 - t282 * t226) * t248 + t229 * t236 -(t282 * t225 + t285 * t226) * t248 + t229 * t237;];
T_reg  = t1;
