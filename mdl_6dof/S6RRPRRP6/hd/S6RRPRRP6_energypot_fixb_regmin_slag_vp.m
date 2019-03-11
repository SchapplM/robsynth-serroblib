% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:09
% EndTime: 2019-03-09 12:14:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (254->68), mult. (614->110), div. (0->0), fcn. (796->12), ass. (0->52)
t311 = pkin(8) + qJ(3);
t285 = sin(pkin(6));
t290 = sin(qJ(2));
t310 = t285 * t290;
t291 = sin(qJ(1));
t309 = t285 * t291;
t295 = cos(qJ(1));
t308 = t285 * t295;
t307 = t291 * t290;
t294 = cos(qJ(2));
t306 = t291 * t294;
t305 = t295 * t290;
t304 = t295 * t294;
t287 = cos(pkin(6));
t273 = t287 * t290 * pkin(2) - t311 * t285;
t281 = t294 * pkin(2) + pkin(1);
t303 = t295 * t273 + t291 * t281;
t302 = -t291 * t273 + t295 * t281;
t301 = pkin(2) * t310 + t311 * t287 + pkin(7);
t300 = g(1) * t291 - g(2) * t295;
t284 = sin(pkin(11));
t286 = cos(pkin(11));
t299 = t294 * t284 + t290 * t286;
t298 = t290 * t284 - t294 * t286;
t272 = t299 * t287;
t262 = t295 * t272 - t291 * t298;
t289 = sin(qJ(4));
t293 = cos(qJ(4));
t256 = t262 * t293 - t289 * t308;
t271 = t298 * t287;
t261 = -t295 * t271 - t291 * t299;
t288 = sin(qJ(5));
t292 = cos(qJ(5));
t249 = t256 * t288 + t261 * t292;
t264 = -t291 * t272 - t295 * t298;
t258 = t264 * t293 + t289 * t309;
t263 = t291 * t271 - t295 * t299;
t251 = t258 * t288 + t263 * t292;
t270 = t299 * t285;
t266 = t270 * t293 + t287 * t289;
t269 = t298 * t285;
t253 = t266 * t288 - t269 * t292;
t297 = g(1) * t251 + g(2) * t249 + g(3) * t253;
t255 = t262 * t289 + t293 * t308;
t257 = t264 * t289 - t293 * t309;
t265 = t270 * t289 - t287 * t293;
t296 = g(1) * t257 + g(2) * t255 + g(3) * t265;
t254 = t266 * t292 + t269 * t288;
t252 = t258 * t292 - t263 * t288;
t250 = t256 * t292 - t261 * t288;
t248 = -g(1) * t252 - g(2) * t250 - g(3) * t254;
t1 = [0, -g(1) * t295 - g(2) * t291, t300, 0, 0, 0, 0, 0, -g(1) * (-t287 * t307 + t304) - g(2) * (t287 * t305 + t306) - g(3) * t310, -g(1) * (-t287 * t306 - t305) - g(2) * (t287 * t304 - t307) - g(3) * t285 * t294, -g(3) * t287 - t300 * t285, -g(1) * t302 - g(2) * t303 - g(3) * t301, 0, 0, 0, 0, 0, -g(1) * t258 - g(2) * t256 - g(3) * t266, t296, 0, 0, 0, 0, 0, t248, t297, t248, -t296, -t297, -g(1) * (t264 * pkin(3) + t258 * pkin(4) + t252 * pkin(5) - t263 * pkin(9) + t257 * pkin(10) + t251 * qJ(6) + t302) - g(2) * (t262 * pkin(3) + t256 * pkin(4) + t250 * pkin(5) - t261 * pkin(9) + t255 * pkin(10) + t249 * qJ(6) + t303) - g(3) * (t270 * pkin(3) + t266 * pkin(4) + t254 * pkin(5) + t269 * pkin(9) + t265 * pkin(10) + t253 * qJ(6) + t301);];
U_reg  = t1;
