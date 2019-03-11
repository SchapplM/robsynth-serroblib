% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR15_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:39
% EndTime: 2019-03-09 20:38:39
% DurationCPUTime: 0.22s
% Computational Cost: add. (853->63), mult. (2290->128), div. (0->0), fcn. (1882->12), ass. (0->56)
t283 = cos(pkin(6)) * qJD(1);
t257 = qJD(2) + t283;
t259 = sin(pkin(7));
t269 = cos(qJ(2));
t290 = cos(pkin(7));
t277 = t269 * t290;
t260 = sin(pkin(6));
t284 = qJD(1) * t260;
t293 = t257 * t259 + t277 * t284;
t265 = sin(qJ(2));
t278 = t269 * t284;
t282 = pkin(1) * t283;
t285 = pkin(9) * t278 + t265 * t282;
t239 = t293 * pkin(10) + t285;
t247 = (-pkin(10) * t259 * t265 - pkin(2) * t269 - pkin(1)) * t284;
t264 = sin(qJ(3));
t268 = cos(qJ(3));
t256 = t269 * t282;
t279 = t265 * t284;
t242 = t257 * pkin(2) + t256 + (-t290 * pkin(10) - pkin(9)) * t279;
t276 = t290 * t242;
t292 = -t264 * t239 + (t247 * t259 + t276) * t268;
t291 = pkin(3) + pkin(11);
t270 = qJD(1) ^ 2;
t288 = t260 ^ 2 * t270;
t287 = t259 * t264;
t244 = t257 * t287 + (t264 * t277 + t265 * t268) * t284;
t249 = -t290 * t257 + t259 * t278 - qJD(3);
t271 = qJD(4) - t292;
t222 = t244 * pkin(4) + t291 * t249 + t271;
t243 = t264 * t279 - t293 * t268;
t231 = -t259 * t242 + t290 * t247;
t273 = -t244 * qJ(4) + t231;
t224 = t291 * t243 + t273;
t263 = sin(qJ(5));
t267 = cos(qJ(5));
t286 = t263 * t222 + t267 * t224;
t281 = t269 * t288;
t280 = t268 * t239 + t247 * t287 + t264 * t276;
t233 = -t267 * t243 - t263 * t249;
t228 = t249 * qJ(4) - t280;
t275 = t267 * t222 - t263 * t224;
t225 = -t243 * pkin(4) - t228;
t266 = cos(qJ(6));
t262 = sin(qJ(6));
t241 = qJD(5) + t244;
t234 = t263 * t243 - t267 * t249;
t232 = qJD(6) + t233;
t230 = t266 * t234 + t262 * t241;
t229 = t262 * t234 - t266 * t241;
t227 = t249 * pkin(3) + t271;
t226 = t243 * pkin(3) + t273;
t220 = t233 * pkin(5) - t234 * pkin(12) + t225;
t219 = t241 * pkin(12) + t286;
t218 = -t241 * pkin(5) - t275;
t1 = [t270 / 0.2e1, 0, 0, t265 ^ 2 * t288 / 0.2e1, t265 * t281, t257 * t279, t257 * t278, t257 ^ 2 / 0.2e1, pkin(1) * t281 + (-pkin(9) * t279 + t256) * t257, -pkin(1) * t265 * t288 - t285 * t257, t244 ^ 2 / 0.2e1, -t244 * t243, -t244 * t249, t243 * t249, t249 ^ 2 / 0.2e1, t231 * t243 - t292 * t249, t231 * t244 + t280 * t249, t227 * t244 + t228 * t243, -t226 * t243 - t227 * t249, -t226 * t244 + t228 * t249, t226 ^ 2 / 0.2e1 + t228 ^ 2 / 0.2e1 + t227 ^ 2 / 0.2e1, t234 ^ 2 / 0.2e1, -t234 * t233, t234 * t241, -t233 * t241, t241 ^ 2 / 0.2e1, t225 * t233 + t275 * t241, t225 * t234 - t286 * t241, t230 ^ 2 / 0.2e1, -t230 * t229, t230 * t232, -t229 * t232, t232 ^ 2 / 0.2e1 (-t262 * t219 + t266 * t220) * t232 + t218 * t229 -(t266 * t219 + t262 * t220) * t232 + t218 * t230;];
T_reg  = t1;
