% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:25:29
% EndTime: 2019-03-10 00:25:29
% DurationCPUTime: 0.30s
% Computational Cost: add. (1522->66), mult. (4007->136), div. (0->0), fcn. (3358->14), ass. (0->58)
t295 = cos(pkin(6)) * qJD(1);
t271 = qJD(2) + t295;
t274 = sin(pkin(7));
t276 = cos(pkin(7));
t285 = cos(qJ(2));
t275 = sin(pkin(6));
t296 = qJD(1) * t275;
t290 = t285 * t296;
t305 = t271 * t274 + t276 * t290;
t281 = sin(qJ(2));
t294 = pkin(1) * t295;
t297 = pkin(9) * t290 + t281 * t294;
t254 = t305 * pkin(10) + t297;
t270 = t285 * t294;
t291 = t281 * t296;
t256 = t271 * pkin(2) + t270 + (-pkin(10) * t276 - pkin(9)) * t291;
t262 = (-pkin(10) * t274 * t281 - pkin(2) * t285 - pkin(1)) * t296;
t280 = sin(qJ(3));
t284 = cos(qJ(3));
t304 = -t280 * t254 + (t256 * t276 + t262 * t274) * t284;
t303 = cos(pkin(13));
t286 = qJD(1) ^ 2;
t301 = t275 ^ 2 * t286;
t300 = t274 * t280;
t299 = t276 * t280;
t257 = t280 * t291 - t305 * t284;
t255 = qJD(4) + t257;
t246 = -t274 * t256 + t276 * t262;
t258 = t271 * t300 + (t281 * t284 + t285 * t299) * t296;
t238 = t257 * pkin(3) - t258 * pkin(11) + t246;
t263 = -t276 * t271 + t274 * t290 - qJD(3);
t292 = t284 * t254 + t256 * t299 + t262 * t300;
t242 = -t263 * pkin(11) + t292;
t279 = sin(qJ(4));
t283 = cos(qJ(4));
t298 = t279 * t238 + t283 * t242;
t231 = t255 * qJ(5) + t298;
t241 = t263 * pkin(3) - t304;
t248 = t279 * t258 + t283 * t263;
t249 = t283 * t258 - t279 * t263;
t234 = t248 * pkin(4) - t249 * qJ(5) + t241;
t273 = sin(pkin(13));
t227 = t303 * t231 + t273 * t234;
t293 = t285 * t301;
t226 = -t273 * t231 + t303 * t234;
t289 = t283 * t238 - t279 * t242;
t230 = -t255 * pkin(4) + qJD(5) - t289;
t282 = cos(qJ(6));
t278 = sin(qJ(6));
t247 = qJD(6) + t248;
t245 = t303 * t249 + t273 * t255;
t244 = t273 * t249 - t303 * t255;
t236 = -t278 * t244 + t282 * t245;
t235 = t282 * t244 + t278 * t245;
t228 = t244 * pkin(5) + t230;
t225 = -t244 * pkin(12) + t227;
t224 = t248 * pkin(5) - t245 * pkin(12) + t226;
t1 = [t286 / 0.2e1, 0, 0, t281 ^ 2 * t301 / 0.2e1, t281 * t293, t271 * t291, t271 * t290, t271 ^ 2 / 0.2e1 (-pkin(9) * t291 + t270) * t271 + pkin(1) * t293, -pkin(1) * t281 * t301 - t297 * t271, t258 ^ 2 / 0.2e1, -t258 * t257, -t258 * t263, t257 * t263, t263 ^ 2 / 0.2e1, t246 * t257 - t263 * t304, t246 * t258 + t292 * t263, t249 ^ 2 / 0.2e1, -t249 * t248, t249 * t255, -t248 * t255, t255 ^ 2 / 0.2e1, t241 * t248 + t289 * t255, t241 * t249 - t298 * t255, t226 * t248 + t230 * t244, -t227 * t248 + t230 * t245, -t226 * t245 - t227 * t244, t227 ^ 2 / 0.2e1 + t226 ^ 2 / 0.2e1 + t230 ^ 2 / 0.2e1, t236 ^ 2 / 0.2e1, -t236 * t235, t236 * t247, -t235 * t247, t247 ^ 2 / 0.2e1 (t282 * t224 - t278 * t225) * t247 + t228 * t235 -(t278 * t224 + t282 * t225) * t247 + t228 * t236;];
T_reg  = t1;
