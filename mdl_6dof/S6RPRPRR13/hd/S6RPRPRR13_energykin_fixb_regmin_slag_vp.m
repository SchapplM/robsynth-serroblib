% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:20
% EndTime: 2019-03-09 04:25:20
% DurationCPUTime: 0.26s
% Computational Cost: add. (704->67), mult. (2306->134), div. (0->0), fcn. (1893->12), ass. (0->56)
t266 = cos(pkin(12));
t291 = cos(pkin(6));
t279 = qJD(1) * t291;
t278 = pkin(1) * t279;
t260 = t266 * t278;
t263 = sin(pkin(12));
t265 = sin(pkin(6));
t288 = t263 * t265;
t290 = cos(pkin(7));
t245 = t260 + (t291 * pkin(2) + (-t290 * pkin(9) - qJ(2)) * t288) * qJD(1);
t264 = sin(pkin(7));
t286 = qJD(1) * t265;
t250 = qJD(2) + (-pkin(9) * t263 * t264 - pkin(2) * t266 - pkin(1)) * t286;
t294 = t290 * t245 + t250 * t264;
t283 = t266 * t286;
t255 = qJ(2) * t283 + t263 * t278;
t282 = t266 * t290;
t277 = t265 * t282;
t281 = t291 * t264;
t242 = (t277 + t281) * qJD(1) * pkin(9) + t255;
t269 = sin(qJ(3));
t272 = cos(qJ(3));
t293 = -t269 * t242 + t294 * t272;
t292 = pkin(3) + pkin(10);
t247 = (t269 * t281 + (t263 * t272 + t269 * t282) * t265) * qJD(1);
t252 = t264 * t283 - t290 * t279 - qJD(3);
t273 = qJD(4) - t293;
t225 = t247 * pkin(4) + t292 * t252 + t273;
t284 = t263 * t286;
t246 = t269 * t284 + (-qJD(1) * t277 - t264 * t279) * t272;
t234 = -t264 * t245 + t290 * t250;
t275 = -t247 * qJ(4) + t234;
t227 = t292 * t246 + t275;
t268 = sin(qJ(5));
t271 = cos(qJ(5));
t287 = t268 * t225 + t271 * t227;
t285 = t272 * t242 + t294 * t269;
t236 = -t271 * t246 - t268 * t252;
t231 = t252 * qJ(4) - t285;
t276 = t271 * t225 - t268 * t227;
t228 = -t246 * pkin(4) - t231;
t270 = cos(qJ(6));
t267 = sin(qJ(6));
t261 = -pkin(1) * t286 + qJD(2);
t254 = -qJ(2) * t284 + t260;
t244 = qJD(5) + t247;
t237 = t268 * t246 - t271 * t252;
t235 = qJD(6) + t236;
t233 = t270 * t237 + t267 * t244;
t232 = t267 * t237 - t270 * t244;
t230 = t252 * pkin(3) + t273;
t229 = t246 * pkin(3) + t275;
t223 = t236 * pkin(5) - t237 * pkin(11) + t228;
t222 = t244 * pkin(11) + t287;
t221 = -t244 * pkin(5) - t276;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (-t261 * t265 * t266 + t291 * t254) * qJD(1) (-t291 * t255 + t261 * t288) * qJD(1) (-t254 * t263 + t255 * t266) * t286, t255 ^ 2 / 0.2e1 + t254 ^ 2 / 0.2e1 + t261 ^ 2 / 0.2e1, t247 ^ 2 / 0.2e1, -t247 * t246, -t247 * t252, t246 * t252, t252 ^ 2 / 0.2e1, t234 * t246 - t293 * t252, t234 * t247 + t285 * t252, t230 * t247 + t231 * t246, -t229 * t246 - t230 * t252, -t229 * t247 + t231 * t252, t229 ^ 2 / 0.2e1 + t231 ^ 2 / 0.2e1 + t230 ^ 2 / 0.2e1, t237 ^ 2 / 0.2e1, -t237 * t236, t237 * t244, -t236 * t244, t244 ^ 2 / 0.2e1, t228 * t236 + t276 * t244, t228 * t237 - t287 * t244, t233 ^ 2 / 0.2e1, -t233 * t232, t233 * t235, -t232 * t235, t235 ^ 2 / 0.2e1 (-t267 * t222 + t270 * t223) * t235 + t221 * t232 -(t270 * t222 + t267 * t223) * t235 + t221 * t233;];
T_reg  = t1;
