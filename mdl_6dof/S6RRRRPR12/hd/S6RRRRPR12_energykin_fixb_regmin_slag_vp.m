% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR12
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
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:39
% EndTime: 2019-03-09 23:44:39
% DurationCPUTime: 0.28s
% Computational Cost: add. (1356->64), mult. (3590->132), div. (0->0), fcn. (3016->14), ass. (0->58)
t304 = cos(pkin(6)) * qJD(1);
t280 = qJD(2) + t304;
t283 = sin(pkin(7));
t286 = cos(pkin(7));
t294 = cos(qJ(2));
t284 = sin(pkin(6));
t305 = qJD(1) * t284;
t299 = t294 * t305;
t314 = t280 * t283 + t286 * t299;
t291 = sin(qJ(2));
t303 = pkin(1) * t304;
t306 = pkin(9) * t299 + t291 * t303;
t263 = t314 * pkin(10) + t306;
t279 = t294 * t303;
t300 = t291 * t305;
t265 = t280 * pkin(2) + t279 + (-pkin(10) * t286 - pkin(9)) * t300;
t271 = (-pkin(10) * t283 * t291 - pkin(2) * t294 - pkin(1)) * t305;
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t313 = -t290 * t263 + (t265 * t286 + t271 * t283) * t293;
t312 = cos(qJ(4));
t295 = qJD(1) ^ 2;
t310 = t284 ^ 2 * t295;
t309 = t283 * t290;
t308 = t286 * t290;
t267 = t280 * t309 + (t291 * t293 + t294 * t308) * t305;
t272 = -t286 * t280 + t283 * t299 - qJD(3);
t289 = sin(qJ(4));
t258 = t267 * t312 - t289 * t272;
t266 = t290 * t300 - t314 * t293;
t264 = qJD(4) + t266;
t255 = -t283 * t265 + t286 * t271;
t248 = t266 * pkin(3) - t267 * pkin(11) + t255;
t301 = t293 * t263 + t265 * t308 + t271 * t309;
t251 = -t272 * pkin(11) + t301;
t298 = t248 * t312 - t289 * t251;
t240 = t264 * pkin(4) - t258 * qJ(5) + t298;
t257 = t289 * t267 + t272 * t312;
t307 = t289 * t248 + t312 * t251;
t242 = -t257 * qJ(5) + t307;
t282 = sin(pkin(13));
t285 = cos(pkin(13));
t237 = t282 * t240 + t285 * t242;
t302 = t294 * t310;
t253 = -t285 * t257 - t282 * t258;
t236 = t285 * t240 - t282 * t242;
t250 = t272 * pkin(3) - t313;
t243 = t257 * pkin(4) + qJD(5) + t250;
t292 = cos(qJ(6));
t288 = sin(qJ(6));
t254 = -t282 * t257 + t285 * t258;
t252 = qJD(6) - t253;
t245 = t292 * t254 + t288 * t264;
t244 = t288 * t254 - t292 * t264;
t238 = -t253 * pkin(5) - t254 * pkin(12) + t243;
t235 = t264 * pkin(12) + t237;
t234 = -t264 * pkin(5) - t236;
t1 = [t295 / 0.2e1, 0, 0, t291 ^ 2 * t310 / 0.2e1, t291 * t302, t280 * t300, t280 * t299, t280 ^ 2 / 0.2e1, pkin(1) * t302 + (-pkin(9) * t300 + t279) * t280, -pkin(1) * t291 * t310 - t280 * t306, t267 ^ 2 / 0.2e1, -t267 * t266, -t267 * t272, t266 * t272, t272 ^ 2 / 0.2e1, t255 * t266 - t313 * t272, t255 * t267 + t272 * t301, t258 ^ 2 / 0.2e1, -t258 * t257, t258 * t264, -t257 * t264, t264 ^ 2 / 0.2e1, t250 * t257 + t264 * t298, t250 * t258 - t264 * t307, -t236 * t254 + t237 * t253, t237 ^ 2 / 0.2e1 + t236 ^ 2 / 0.2e1 + t243 ^ 2 / 0.2e1, t245 ^ 2 / 0.2e1, -t245 * t244, t245 * t252, -t244 * t252, t252 ^ 2 / 0.2e1 (-t288 * t235 + t292 * t238) * t252 + t234 * t244 -(t292 * t235 + t288 * t238) * t252 + t234 * t245;];
T_reg  = t1;
