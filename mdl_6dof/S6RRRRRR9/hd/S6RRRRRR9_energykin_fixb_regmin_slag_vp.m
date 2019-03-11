% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR9
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:12
% EndTime: 2019-03-10 05:39:12
% DurationCPUTime: 0.29s
% Computational Cost: add. (1266->65), mult. (3297->135), div. (0->0), fcn. (2815->14), ass. (0->59)
t296 = cos(pkin(6)) * qJD(1);
t271 = qJD(2) + t296;
t273 = sin(pkin(7));
t275 = cos(pkin(7));
t285 = cos(qJ(2));
t274 = sin(pkin(6));
t297 = qJD(1) * t274;
t291 = t285 * t297;
t307 = t271 * t273 + t275 * t291;
t281 = sin(qJ(2));
t295 = pkin(1) * t296;
t298 = pkin(9) * t291 + t281 * t295;
t254 = pkin(10) * t307 + t298;
t270 = t285 * t295;
t292 = t281 * t297;
t256 = t271 * pkin(2) + t270 + (-pkin(10) * t275 - pkin(9)) * t292;
t262 = (-pkin(10) * t273 * t281 - pkin(2) * t285 - pkin(1)) * t297;
t280 = sin(qJ(3));
t284 = cos(qJ(3));
t306 = -t280 * t254 + (t256 * t275 + t262 * t273) * t284;
t305 = cos(qJ(5));
t286 = qJD(1) ^ 2;
t303 = t274 ^ 2 * t286;
t302 = t273 * t280;
t301 = t275 * t280;
t257 = t280 * t292 - t284 * t307;
t255 = qJD(4) + t257;
t245 = -t273 * t256 + t275 * t262;
t258 = t271 * t302 + (t281 * t284 + t285 * t301) * t297;
t237 = t257 * pkin(3) - t258 * pkin(11) + t245;
t263 = -t275 * t271 + t273 * t291 - qJD(3);
t293 = t284 * t254 + t256 * t301 + t262 * t302;
t241 = -t263 * pkin(11) + t293;
t279 = sin(qJ(4));
t283 = cos(qJ(4));
t299 = t279 * t237 + t283 * t241;
t230 = t255 * pkin(12) + t299;
t240 = t263 * pkin(3) - t306;
t248 = t279 * t258 + t283 * t263;
t249 = t283 * t258 - t279 * t263;
t233 = t248 * pkin(4) - t249 * pkin(12) + t240;
t278 = sin(qJ(5));
t300 = t305 * t230 + t278 * t233;
t294 = t285 * t303;
t290 = -t278 * t230 + t305 * t233;
t289 = t283 * t237 - t279 * t241;
t247 = qJD(5) + t248;
t229 = -t255 * pkin(4) - t289;
t282 = cos(qJ(6));
t277 = sin(qJ(6));
t246 = qJD(6) + t247;
t244 = t305 * t249 + t278 * t255;
t243 = t278 * t249 - t305 * t255;
t235 = -t277 * t243 + t282 * t244;
t234 = t282 * t243 + t277 * t244;
t227 = t243 * pkin(5) + t229;
t226 = -t243 * pkin(13) + t300;
t225 = t247 * pkin(5) - t244 * pkin(13) + t290;
t1 = [t286 / 0.2e1, 0, 0, t281 ^ 2 * t303 / 0.2e1, t281 * t294, t271 * t292, t271 * t291, t271 ^ 2 / 0.2e1 (-pkin(9) * t292 + t270) * t271 + pkin(1) * t294, -pkin(1) * t281 * t303 - t298 * t271, t258 ^ 2 / 0.2e1, -t258 * t257, -t258 * t263, t257 * t263, t263 ^ 2 / 0.2e1, t245 * t257 - t306 * t263, t245 * t258 + t293 * t263, t249 ^ 2 / 0.2e1, -t248 * t249, t255 * t249, -t248 * t255, t255 ^ 2 / 0.2e1, t240 * t248 + t289 * t255, t240 * t249 - t299 * t255, t244 ^ 2 / 0.2e1, -t244 * t243, t247 * t244, -t247 * t243, t247 ^ 2 / 0.2e1, t229 * t243 + t290 * t247, t229 * t244 - t300 * t247, t235 ^ 2 / 0.2e1, -t235 * t234, t235 * t246, -t234 * t246, t246 ^ 2 / 0.2e1 (t282 * t225 - t277 * t226) * t246 + t227 * t234 -(t277 * t225 + t282 * t226) * t246 + t227 * t235;];
T_reg  = t1;
