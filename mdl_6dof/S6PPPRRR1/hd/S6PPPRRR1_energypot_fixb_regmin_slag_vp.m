% Calculate minimal parameter regressor of potential energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPPRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:34
% EndTime: 2019-03-08 18:40:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (374->67), mult. (1050->127), div. (0->0), fcn. (1385->18), ass. (0->60)
t261 = sin(pkin(13));
t265 = sin(pkin(6));
t293 = t261 * t265;
t262 = sin(pkin(12));
t271 = cos(pkin(6));
t292 = t262 * t271;
t264 = sin(pkin(7));
t291 = t264 * t265;
t290 = t264 * t271;
t268 = cos(pkin(12));
t289 = t265 * t268;
t270 = cos(pkin(7));
t288 = t265 * t270;
t267 = cos(pkin(13));
t287 = t267 * t270;
t286 = t268 * t271;
t285 = t271 * qJ(2) + qJ(1);
t284 = t262 * t265 * qJ(2) + t268 * pkin(1);
t283 = t262 * pkin(1) - qJ(2) * t289;
t253 = t261 * t286 + t262 * t267;
t260 = sin(pkin(14));
t266 = cos(pkin(14));
t252 = -t262 * t261 + t267 * t286;
t279 = t252 * t270 - t264 * t289;
t242 = -t253 * t260 + t279 * t266;
t249 = -t252 * t264 - t268 * t288;
t263 = sin(pkin(8));
t269 = cos(pkin(8));
t282 = t242 * t269 + t249 * t263;
t255 = -t261 * t292 + t268 * t267;
t254 = -t268 * t261 - t267 * t292;
t278 = t254 * t270 + t262 * t291;
t244 = -t255 * t260 + t278 * t266;
t250 = -t254 * t264 + t262 * t288;
t281 = t244 * t269 + t250 * t263;
t247 = t266 * t290 + (-t260 * t261 + t266 * t287) * t265;
t251 = -t267 * t291 + t271 * t270;
t280 = t247 * t269 + t251 * t263;
t277 = cos(qJ(4));
t276 = cos(qJ(5));
t275 = cos(qJ(6));
t274 = sin(qJ(4));
t273 = sin(qJ(5));
t272 = sin(qJ(6));
t248 = t266 * t293 + (t265 * t287 + t290) * t260;
t246 = -t247 * t263 + t251 * t269;
t245 = t255 * t266 + t278 * t260;
t243 = t253 * t266 + t279 * t260;
t241 = -t244 * t263 + t250 * t269;
t240 = -t242 * t263 + t249 * t269;
t239 = t248 * t277 + t280 * t274;
t238 = t248 * t274 - t280 * t277;
t237 = t245 * t277 + t281 * t274;
t236 = t245 * t274 - t281 * t277;
t235 = t243 * t277 + t282 * t274;
t234 = t243 * t274 - t282 * t277;
t233 = t239 * t276 + t246 * t273;
t232 = t237 * t276 + t241 * t273;
t231 = t235 * t276 + t240 * t273;
t1 = [-g(3) * qJ(1), -g(1) * t284 - g(2) * t283 - g(3) * t285, -g(1) * (t255 * pkin(2) + t284) - g(2) * (t253 * pkin(2) + t283) - g(3) * (pkin(2) * t293 + t285) + (-g(1) * t250 - g(2) * t249 - g(3) * t251) * qJ(3), 0, -g(1) * t237 - g(2) * t235 - g(3) * t239, g(1) * t236 + g(2) * t234 + g(3) * t238, 0, 0, 0, 0, 0, -g(1) * t232 - g(2) * t231 - g(3) * t233, -g(1) * (-t237 * t273 + t241 * t276) - g(2) * (-t235 * t273 + t240 * t276) - g(3) * (-t239 * t273 + t246 * t276) 0, 0, 0, 0, 0, -g(1) * (t232 * t275 + t236 * t272) - g(2) * (t231 * t275 + t234 * t272) - g(3) * (t233 * t275 + t238 * t272) -g(1) * (-t232 * t272 + t236 * t275) - g(2) * (-t231 * t272 + t234 * t275) - g(3) * (-t233 * t272 + t238 * t275);];
U_reg  = t1;
