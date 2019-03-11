% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR11
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
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:15:11
% EndTime: 2019-03-09 04:15:11
% DurationCPUTime: 0.28s
% Computational Cost: add. (1161->71), mult. (3769->143), div. (0->0), fcn. (3167->14), ass. (0->59)
t283 = sin(pkin(12));
t309 = cos(pkin(6));
t297 = qJD(1) * t309;
t296 = pkin(1) * t297;
t286 = cos(pkin(12));
t285 = sin(pkin(6));
t302 = qJD(1) * t285;
t299 = t286 * t302;
t274 = qJ(2) * t299 + t283 * t296;
t287 = cos(pkin(7));
t284 = sin(pkin(7));
t298 = t309 * t284;
t305 = t285 * t286;
t262 = (t287 * t305 + t298) * qJD(1) * pkin(9) + t274;
t279 = t286 * t296;
t306 = t283 * t285;
t264 = t279 + (t309 * pkin(2) + (-pkin(9) * t287 - qJ(2)) * t306) * qJD(1);
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t269 = qJD(2) + (-pkin(9) * t283 * t284 - pkin(2) * t286 - pkin(1)) * t302;
t307 = t269 * t284;
t310 = -t290 * t262 + (t264 * t287 + t307) * t293;
t308 = cos(pkin(13));
t304 = t287 * t290;
t254 = -t284 * t264 + t287 * t269;
t300 = t283 * t302;
t265 = t290 * t300 + (-t284 * t297 - t287 * t299) * t293;
t266 = (t290 * t298 + (t283 * t293 + t286 * t304) * t285) * qJD(1);
t248 = t265 * pkin(3) - t266 * qJ(4) + t254;
t271 = t284 * t299 - t287 * t297 - qJD(3);
t301 = t293 * t262 + t264 * t304 + t290 * t307;
t250 = -t271 * qJ(4) + t301;
t282 = sin(pkin(13));
t240 = t308 * t248 - t282 * t250;
t257 = t308 * t266 - t282 * t271;
t237 = t265 * pkin(4) - t257 * pkin(10) + t240;
t241 = t282 * t248 + t308 * t250;
t256 = t282 * t266 + t308 * t271;
t239 = -t256 * pkin(10) + t241;
t289 = sin(qJ(5));
t292 = cos(qJ(5));
t303 = t289 * t237 + t292 * t239;
t252 = t292 * t256 + t289 * t257;
t295 = t292 * t237 - t289 * t239;
t249 = t271 * pkin(3) + qJD(4) - t310;
t242 = t256 * pkin(4) + t249;
t291 = cos(qJ(6));
t288 = sin(qJ(6));
t280 = -pkin(1) * t302 + qJD(2);
t273 = -qJ(2) * t300 + t279;
t263 = qJD(5) + t265;
t253 = -t289 * t256 + t292 * t257;
t251 = qJD(6) + t252;
t244 = t291 * t253 + t288 * t263;
t243 = t288 * t253 - t291 * t263;
t235 = t252 * pkin(5) - t253 * pkin(11) + t242;
t234 = t263 * pkin(11) + t303;
t233 = -t263 * pkin(5) - t295;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t309 * t273 - t280 * t305) * qJD(1) (-t309 * t274 + t280 * t306) * qJD(1) (-t273 * t283 + t274 * t286) * t302, t274 ^ 2 / 0.2e1 + t273 ^ 2 / 0.2e1 + t280 ^ 2 / 0.2e1, t266 ^ 2 / 0.2e1, -t266 * t265, -t266 * t271, t265 * t271, t271 ^ 2 / 0.2e1, t254 * t265 - t271 * t310, t254 * t266 + t301 * t271, t240 * t265 + t249 * t256, -t241 * t265 + t249 * t257, -t240 * t257 - t241 * t256, t241 ^ 2 / 0.2e1 + t240 ^ 2 / 0.2e1 + t249 ^ 2 / 0.2e1, t253 ^ 2 / 0.2e1, -t253 * t252, t253 * t263, -t252 * t263, t263 ^ 2 / 0.2e1, t242 * t252 + t295 * t263, t242 * t253 - t303 * t263, t244 ^ 2 / 0.2e1, -t244 * t243, t244 * t251, -t243 * t251, t251 ^ 2 / 0.2e1 (-t288 * t234 + t291 * t235) * t251 + t233 * t243 -(t291 * t234 + t288 * t235) * t251 + t233 * t244;];
T_reg  = t1;
