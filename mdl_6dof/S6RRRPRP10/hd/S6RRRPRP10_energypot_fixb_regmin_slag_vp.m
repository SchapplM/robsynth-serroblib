% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:36:55
% EndTime: 2019-03-09 17:36:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (245->81), mult. (492->117), div. (0->0), fcn. (620->12), ass. (0->48)
t277 = sin(pkin(6));
t282 = sin(qJ(2));
t301 = t277 * t282;
t283 = sin(qJ(1));
t300 = t277 * t283;
t284 = cos(qJ(3));
t299 = t277 * t284;
t285 = cos(qJ(2));
t298 = t277 * t285;
t286 = cos(qJ(1));
t297 = t277 * t286;
t296 = t283 * t282;
t295 = t283 * t285;
t294 = t286 * t282;
t293 = t286 * t285;
t279 = cos(pkin(6));
t292 = pkin(2) * t301 + t279 * pkin(8) + pkin(7);
t263 = -t279 * t296 + t293;
t291 = t286 * pkin(1) + t263 * pkin(2) + pkin(8) * t300;
t276 = sin(pkin(11));
t290 = pkin(4) * t276 + pkin(9);
t261 = t279 * t294 + t295;
t289 = t283 * pkin(1) + t261 * pkin(2) - pkin(8) * t297;
t281 = sin(qJ(3));
t251 = t261 * t284 - t281 * t297;
t260 = -t279 * t293 + t296;
t275 = pkin(11) + qJ(5);
t270 = sin(t275);
t271 = cos(t275);
t244 = t251 * t270 - t260 * t271;
t253 = t263 * t284 + t281 * t300;
t262 = t279 * t295 + t294;
t246 = t253 * t270 - t262 * t271;
t259 = t279 * t281 + t282 * t299;
t248 = t259 * t270 + t271 * t298;
t288 = g(1) * t246 + g(2) * t244 + g(3) * t248;
t250 = t261 * t281 + t284 * t297;
t252 = t263 * t281 - t283 * t299;
t258 = -t279 * t284 + t281 * t301;
t287 = g(1) * t252 + g(2) * t250 + g(3) * t258;
t280 = -pkin(10) - qJ(4);
t278 = cos(pkin(11));
t269 = t278 * pkin(4) + pkin(3);
t249 = t259 * t271 - t270 * t298;
t247 = t253 * t271 + t262 * t270;
t245 = t251 * t271 + t260 * t270;
t242 = -g(1) * t247 - g(2) * t245 - g(3) * t249;
t1 = [0, -g(1) * t286 - g(2) * t283, g(1) * t283 - g(2) * t286, 0, 0, 0, 0, 0, -g(1) * t263 - g(2) * t261 - g(3) * t301, g(1) * t262 + g(2) * t260 - g(3) * t298, 0, 0, 0, 0, 0, -g(1) * t253 - g(2) * t251 - g(3) * t259, t287, -g(1) * (t253 * t278 + t262 * t276) - g(2) * (t251 * t278 + t260 * t276) - g(3) * (t259 * t278 - t276 * t298) -g(1) * (-t253 * t276 + t262 * t278) - g(2) * (-t251 * t276 + t260 * t278) - g(3) * (-t259 * t276 - t278 * t298) -t287, -g(1) * (t253 * pkin(3) + t262 * pkin(9) + t252 * qJ(4) + t291) - g(2) * (t251 * pkin(3) + t260 * pkin(9) + t250 * qJ(4) + t289) - g(3) * (t259 * pkin(3) - pkin(9) * t298 + t258 * qJ(4) + t292) 0, 0, 0, 0, 0, t242, t288, t242, -t287, -t288, -g(1) * (t247 * pkin(5) + t246 * qJ(6) - t252 * t280 + t253 * t269 + t290 * t262 + t291) - g(2) * (t245 * pkin(5) + t244 * qJ(6) - t250 * t280 + t251 * t269 + t290 * t260 + t289) - g(3) * (t249 * pkin(5) + t248 * qJ(6) - t258 * t280 + t259 * t269 - t290 * t298 + t292);];
U_reg  = t1;
