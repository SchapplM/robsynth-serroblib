% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:52:56
% EndTime: 2019-03-09 05:52:56
% DurationCPUTime: 0.23s
% Computational Cost: add. (839->68), mult. (2718->135), div. (0->0), fcn. (2244->12), ass. (0->56)
t260 = sin(pkin(12));
t287 = cos(pkin(6));
t276 = qJD(1) * t287;
t274 = pkin(1) * t276;
t263 = cos(pkin(12));
t262 = sin(pkin(6));
t281 = qJD(1) * t262;
t278 = t263 * t281;
t252 = qJ(2) * t278 + t260 * t274;
t264 = cos(pkin(7));
t261 = sin(pkin(7));
t277 = t287 * t261;
t284 = t262 * t263;
t240 = (t264 * t284 + t277) * qJD(1) * pkin(9) + t252;
t257 = t263 * t274;
t285 = t260 * t262;
t242 = t257 + (t287 * pkin(2) + (-pkin(9) * t264 - qJ(2)) * t285) * qJD(1);
t267 = sin(qJ(3));
t270 = cos(qJ(3));
t247 = qJD(2) + (-pkin(9) * t260 * t261 - pkin(2) * t263 - pkin(1)) * t281;
t286 = t247 * t261;
t289 = -t267 * t240 + (t242 * t264 + t286) * t270;
t288 = pkin(4) + pkin(11);
t283 = t264 * t267;
t232 = -t261 * t242 + t264 * t247;
t279 = t260 * t281;
t243 = t267 * t279 + (-t261 * t276 - t264 * t278) * t270;
t244 = (t267 * t277 + (t260 * t270 + t263 * t283) * t262) * qJD(1);
t225 = t243 * pkin(3) - t244 * pkin(10) + t232;
t249 = t261 * t278 - t264 * t276 - qJD(3);
t280 = t270 * t240 + t242 * t283 + t267 * t286;
t229 = -t249 * pkin(10) + t280;
t266 = sin(qJ(4));
t269 = cos(qJ(4));
t282 = t266 * t225 + t269 * t229;
t275 = t269 * t225 - t266 * t229;
t241 = qJD(4) + t243;
t222 = -t241 * qJ(5) - t282;
t273 = qJD(5) - t275;
t236 = t269 * t244 - t266 * t249;
t228 = t249 * pkin(3) - t289;
t271 = -t236 * qJ(5) + t228;
t268 = cos(qJ(6));
t265 = sin(qJ(6));
t258 = -pkin(1) * t281 + qJD(2);
t251 = -qJ(2) * t279 + t257;
t235 = t266 * t244 + t269 * t249;
t234 = qJD(6) + t236;
t231 = t265 * t235 + t268 * t241;
t230 = -t268 * t235 + t265 * t241;
t223 = t235 * pkin(4) + t271;
t221 = -t241 * pkin(4) + t273;
t220 = t288 * t235 + t271;
t219 = -t235 * pkin(5) - t222;
t218 = t236 * pkin(5) - t288 * t241 + t273;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t287 * t251 - t258 * t284) * qJD(1) (-t287 * t252 + t258 * t285) * qJD(1) (-t251 * t260 + t252 * t263) * t281, t252 ^ 2 / 0.2e1 + t251 ^ 2 / 0.2e1 + t258 ^ 2 / 0.2e1, t244 ^ 2 / 0.2e1, -t244 * t243, -t244 * t249, t243 * t249, t249 ^ 2 / 0.2e1, t232 * t243 - t289 * t249, t232 * t244 + t280 * t249, t236 ^ 2 / 0.2e1, -t236 * t235, t236 * t241, -t235 * t241, t241 ^ 2 / 0.2e1, t228 * t235 + t275 * t241, t228 * t236 - t282 * t241, t221 * t236 + t222 * t235, t221 * t241 - t223 * t235, -t222 * t241 - t223 * t236, t223 ^ 2 / 0.2e1 + t222 ^ 2 / 0.2e1 + t221 ^ 2 / 0.2e1, t231 ^ 2 / 0.2e1, -t231 * t230, t231 * t234, -t230 * t234, t234 ^ 2 / 0.2e1 (t268 * t218 - t265 * t220) * t234 + t219 * t230 -(t265 * t218 + t268 * t220) * t234 + t219 * t231;];
T_reg  = t1;
