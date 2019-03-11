% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:10:16
% EndTime: 2019-03-09 19:10:16
% DurationCPUTime: 0.23s
% Computational Cost: add. (1253->65), mult. (3497->134), div. (0->0), fcn. (2960->14), ass. (0->60)
t312 = cos(qJ(3));
t311 = cos(pkin(7));
t283 = sin(pkin(6));
t293 = qJD(1) ^ 2;
t310 = t283 ^ 2 * t293;
t282 = sin(pkin(7));
t288 = sin(qJ(3));
t309 = t282 * t288;
t305 = cos(pkin(6)) * qJD(1);
t279 = qJD(2) + t305;
t289 = sin(qJ(2));
t292 = cos(qJ(2));
t297 = t292 * t311;
t306 = qJD(1) * t283;
t265 = t279 * t309 + (t288 * t297 + t312 * t289) * t306;
t299 = t292 * t306;
t271 = -t311 * t279 + t282 * t299 - qJD(3);
t294 = t297 * t306;
t304 = pkin(1) * t305;
t307 = pkin(9) * t299 + t289 * t304;
t262 = (t279 * t282 + t294) * pkin(10) + t307;
t270 = (-pkin(10) * t282 * t289 - pkin(2) * t292 - pkin(1)) * t306;
t278 = t292 * t304;
t300 = t289 * t306;
t263 = t279 * pkin(2) + t278 + (-t311 * pkin(10) - pkin(9)) * t300;
t298 = t263 * t311;
t301 = t282 * t312;
t296 = -t262 * t288 + t270 * t301 + t312 * t298;
t245 = -pkin(3) * t271 - qJ(4) * t265 + t296;
t264 = -t279 * t301 + t288 * t300 - t312 * t294;
t302 = t312 * t262 + t270 * t309 + t288 * t298;
t248 = -qJ(4) * t264 + t302;
t281 = sin(pkin(13));
t284 = cos(pkin(13));
t239 = t281 * t245 + t284 * t248;
t237 = -pkin(11) * t271 + t239;
t257 = -t263 * t282 + t311 * t270;
t252 = pkin(3) * t264 + qJD(4) + t257;
t255 = -t284 * t264 - t265 * t281;
t256 = -t264 * t281 + t265 * t284;
t241 = -pkin(4) * t255 - pkin(11) * t256 + t252;
t287 = sin(qJ(5));
t291 = cos(qJ(5));
t308 = t291 * t237 + t287 * t241;
t303 = t292 * t310;
t238 = t245 * t284 - t281 * t248;
t250 = t256 * t287 + t291 * t271;
t295 = -t237 * t287 + t241 * t291;
t236 = pkin(4) * t271 - t238;
t290 = cos(qJ(6));
t286 = sin(qJ(6));
t254 = qJD(5) - t255;
t251 = t256 * t291 - t271 * t287;
t249 = qJD(6) + t250;
t243 = t251 * t290 + t254 * t286;
t242 = t251 * t286 - t290 * t254;
t234 = pkin(5) * t250 - pkin(12) * t251 + t236;
t233 = pkin(12) * t254 + t308;
t232 = -pkin(5) * t254 - t295;
t1 = [t293 / 0.2e1, 0, 0, t289 ^ 2 * t310 / 0.2e1, t289 * t303, t279 * t300, t279 * t299, t279 ^ 2 / 0.2e1, pkin(1) * t303 + (-pkin(9) * t300 + t278) * t279, -pkin(1) * t289 * t310 - t307 * t279, t265 ^ 2 / 0.2e1, -t265 * t264, -t265 * t271, t264 * t271, t271 ^ 2 / 0.2e1, t257 * t264 - t296 * t271, t257 * t265 + t302 * t271, -t238 * t256 + t239 * t255, t239 ^ 2 / 0.2e1 + t238 ^ 2 / 0.2e1 + t252 ^ 2 / 0.2e1, t251 ^ 2 / 0.2e1, -t251 * t250, t251 * t254, -t250 * t254, t254 ^ 2 / 0.2e1, t236 * t250 + t295 * t254, t236 * t251 - t308 * t254, t243 ^ 2 / 0.2e1, -t243 * t242, t243 * t249, -t242 * t249, t249 ^ 2 / 0.2e1 (-t233 * t286 + t234 * t290) * t249 + t232 * t242 -(t233 * t290 + t234 * t286) * t249 + t232 * t243;];
T_reg  = t1;
