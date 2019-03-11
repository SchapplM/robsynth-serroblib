% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP9
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:33:57
% EndTime: 2019-03-09 12:33:57
% DurationCPUTime: 0.19s
% Computational Cost: add. (181->79), mult. (309->118), div. (0->0), fcn. (374->12), ass. (0->48)
t256 = sin(pkin(11));
t286 = pkin(3) * t256;
t264 = sin(qJ(1));
t285 = g(1) * t264;
t267 = cos(qJ(1));
t284 = g(2) * t267;
t259 = cos(pkin(6));
t283 = t259 * pkin(8) + pkin(7);
t266 = cos(qJ(2));
t271 = t266 * t267;
t263 = sin(qJ(2));
t274 = t263 * t264;
t240 = -t259 * t271 + t274;
t262 = sin(qJ(5));
t282 = t240 * t262;
t272 = t264 * t266;
t273 = t263 * t267;
t242 = t259 * t272 + t273;
t281 = t242 * t262;
t280 = t256 * t259;
t257 = sin(pkin(6));
t279 = t257 * t263;
t278 = t257 * t264;
t277 = t257 * t266;
t276 = t257 * t267;
t275 = t262 * t266;
t270 = t267 * pkin(1) + pkin(8) * t278;
t241 = t259 * t273 + t272;
t255 = pkin(11) + qJ(4);
t250 = sin(t255);
t251 = cos(t255);
t234 = t241 * t250 + t251 * t276;
t243 = -t259 * t274 + t271;
t236 = t243 * t250 - t251 * t278;
t238 = t250 * t279 - t259 * t251;
t269 = g(1) * t236 + g(2) * t234 + g(3) * t238;
t268 = -g(1) * t242 - g(2) * t240 + g(3) * t277;
t265 = cos(qJ(5));
t261 = -pkin(9) - qJ(3);
t260 = -qJ(6) - pkin(10);
t258 = cos(pkin(11));
t253 = t264 * pkin(1);
t249 = pkin(5) * t265 + pkin(4);
t248 = pkin(3) * t258 + pkin(2);
t239 = t250 * t259 + t251 * t279;
t237 = t243 * t251 + t250 * t278;
t235 = t241 * t251 - t250 * t276;
t1 = [0, -g(1) * t267 - g(2) * t264, -t284 + t285, 0, 0, 0, 0, 0, -g(1) * t243 - g(2) * t241 - g(3) * t279, -t268, -g(1) * (t243 * t258 + t256 * t278) - g(2) * (t241 * t258 - t256 * t276) - g(3) * (t258 * t279 + t280) -g(1) * (-t243 * t256 + t258 * t278) - g(2) * (-t241 * t256 - t258 * t276) - g(3) * (-t256 * t279 + t258 * t259) t268, -g(1) * (pkin(2) * t243 + qJ(3) * t242 + t270) - g(2) * (t241 * pkin(2) - pkin(8) * t276 + t240 * qJ(3) + t253) - g(3) * ((pkin(2) * t263 - qJ(3) * t266) * t257 + t283) 0, 0, 0, 0, 0, -g(1) * t237 - g(2) * t235 - g(3) * t239, t269, 0, 0, 0, 0, 0, -g(1) * (t237 * t265 + t281) - g(2) * (t235 * t265 + t282) - g(3) * (t239 * t265 - t257 * t275) -g(1) * (-t237 * t262 + t242 * t265) - g(2) * (-t235 * t262 + t240 * t265) - g(3) * (-t239 * t262 - t265 * t277) -t269, -g(1) * (pkin(5) * t281 - t236 * t260 + t237 * t249 - t242 * t261 + t243 * t248 + t270) - g(2) * (pkin(5) * t282 - t234 * t260 + t235 * t249 - t240 * t261 + t241 * t248 + t253) - g(3) * (pkin(3) * t280 - t238 * t260 + t239 * t249 + t283) + (-t285 * t286 - g(3) * (-pkin(5) * t275 + t248 * t263 + t261 * t266) - (-pkin(8) - t286) * t284) * t257;];
U_reg  = t1;
