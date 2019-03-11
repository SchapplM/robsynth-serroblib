% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:05
% EndTime: 2019-03-09 04:05:05
% DurationCPUTime: 0.25s
% Computational Cost: add. (1028->69), mult. (3513->140), div. (0->0), fcn. (2971->14), ass. (0->60)
t312 = cos(qJ(3));
t311 = cos(pkin(6));
t310 = cos(pkin(7));
t285 = sin(pkin(12));
t287 = sin(pkin(6));
t309 = t285 * t287;
t289 = cos(pkin(12));
t308 = t287 * t289;
t292 = sin(qJ(3));
t286 = sin(pkin(7));
t300 = t311 * t286;
t301 = t292 * t310;
t267 = (t292 * t300 + (t312 * t285 + t289 * t301) * t287) * qJD(1);
t299 = qJD(1) * t311;
t306 = qJD(1) * t287;
t302 = t289 * t306;
t273 = t286 * t302 - t310 * t299 - qJD(3);
t298 = pkin(1) * t299;
t276 = qJ(2) * t302 + t285 * t298;
t263 = (t310 * t308 + t300) * qJD(1) * pkin(9) + t276;
t281 = t289 * t298;
t265 = t281 + (t311 * pkin(2) + (-t310 * pkin(9) - qJ(2)) * t309) * qJD(1);
t271 = qJD(2) + (-pkin(9) * t285 * t286 - pkin(2) * t289 - pkin(1)) * t306;
t297 = t310 * t312;
t304 = t286 * t312;
t296 = -t292 * t263 + t265 * t297 + t271 * t304;
t247 = -t273 * pkin(3) - t267 * qJ(4) + t296;
t303 = t285 * t306;
t266 = t292 * t303 - t297 * t302 - t299 * t304;
t305 = t286 * t292 * t271 + t312 * t263 + t265 * t301;
t250 = -t266 * qJ(4) + t305;
t284 = sin(pkin(13));
t288 = cos(pkin(13));
t241 = t284 * t247 + t288 * t250;
t239 = -t273 * pkin(10) + t241;
t259 = -t286 * t265 + t310 * t271;
t251 = t266 * pkin(3) + qJD(4) + t259;
t257 = -t288 * t266 - t284 * t267;
t258 = -t284 * t266 + t288 * t267;
t243 = -t257 * pkin(4) - t258 * pkin(10) + t251;
t291 = sin(qJ(5));
t294 = cos(qJ(5));
t307 = t294 * t239 + t291 * t243;
t240 = t288 * t247 - t284 * t250;
t253 = t291 * t258 + t294 * t273;
t295 = -t291 * t239 + t294 * t243;
t238 = t273 * pkin(4) - t240;
t293 = cos(qJ(6));
t290 = sin(qJ(6));
t282 = -pkin(1) * t306 + qJD(2);
t275 = -qJ(2) * t303 + t281;
t256 = qJD(5) - t257;
t254 = t294 * t258 - t291 * t273;
t252 = qJD(6) + t253;
t245 = t293 * t254 + t290 * t256;
t244 = t290 * t254 - t293 * t256;
t236 = t253 * pkin(5) - t254 * pkin(11) + t238;
t235 = t256 * pkin(11) + t307;
t234 = -t256 * pkin(5) - t295;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t311 * t275 - t282 * t308) * qJD(1) (-t311 * t276 + t282 * t309) * qJD(1) (-t275 * t285 + t276 * t289) * t306, t276 ^ 2 / 0.2e1 + t275 ^ 2 / 0.2e1 + t282 ^ 2 / 0.2e1, t267 ^ 2 / 0.2e1, -t267 * t266, -t267 * t273, t266 * t273, t273 ^ 2 / 0.2e1, t259 * t266 - t296 * t273, t259 * t267 + t305 * t273, -t240 * t258 + t241 * t257, t241 ^ 2 / 0.2e1 + t240 ^ 2 / 0.2e1 + t251 ^ 2 / 0.2e1, t254 ^ 2 / 0.2e1, -t254 * t253, t254 * t256, -t253 * t256, t256 ^ 2 / 0.2e1, t238 * t253 + t295 * t256, t238 * t254 - t307 * t256, t245 ^ 2 / 0.2e1, -t245 * t244, t245 * t252, -t244 * t252, t252 ^ 2 / 0.2e1 (-t290 * t235 + t293 * t236) * t252 + t234 * t244 -(t293 * t235 + t290 * t236) * t252 + t234 * t245;];
T_reg  = t1;
