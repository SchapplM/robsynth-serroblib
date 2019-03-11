% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:48
% EndTime: 2019-03-10 01:55:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (225->70), mult. (377->107), div. (0->0), fcn. (477->12), ass. (0->46)
t276 = sin(qJ(1));
t280 = cos(qJ(1));
t297 = -g(1) * t276 + g(2) * t280;
t271 = sin(pkin(6));
t275 = sin(qJ(2));
t294 = t271 * t275;
t293 = t271 * t276;
t278 = cos(qJ(3));
t292 = t271 * t278;
t279 = cos(qJ(2));
t291 = t271 * t279;
t290 = t271 * t280;
t272 = cos(pkin(6));
t274 = sin(qJ(3));
t289 = t272 * t274;
t288 = t276 * t275;
t287 = t276 * t279;
t286 = t280 * t275;
t285 = t280 * t279;
t261 = t272 * t286 + t287;
t270 = qJ(3) + qJ(4);
t268 = sin(t270);
t269 = cos(t270);
t253 = t261 * t269 - t268 * t290;
t260 = -t272 * t285 + t288;
t273 = sin(qJ(5));
t277 = cos(qJ(5));
t246 = t253 * t273 - t260 * t277;
t263 = -t272 * t288 + t285;
t255 = t263 * t269 + t268 * t293;
t262 = t272 * t287 + t286;
t248 = t255 * t273 - t262 * t277;
t257 = t272 * t268 + t269 * t294;
t250 = t257 * t273 + t277 * t291;
t283 = g(1) * t248 + g(2) * t246 + g(3) * t250;
t252 = t261 * t268 + t269 * t290;
t254 = t263 * t268 - t269 * t293;
t256 = t268 * t294 - t272 * t269;
t282 = g(1) * t254 + g(2) * t252 + g(3) * t256;
t281 = -pkin(10) - pkin(9);
t267 = t278 * pkin(3) + pkin(2);
t251 = t257 * t277 - t273 * t291;
t249 = t255 * t277 + t262 * t273;
t247 = t253 * t277 + t260 * t273;
t245 = -g(1) * t249 - g(2) * t247 - g(3) * t251;
t1 = [0, -g(1) * t280 - g(2) * t276, -t297, 0, 0, 0, 0, 0, -g(1) * t263 - g(2) * t261 - g(3) * t294, g(1) * t262 + g(2) * t260 - g(3) * t291, 0, 0, 0, 0, 0, -g(1) * (t263 * t278 + t274 * t293) - g(2) * (t261 * t278 - t274 * t290) - g(3) * (t275 * t292 + t289) -g(1) * (-t263 * t274 + t276 * t292) - g(2) * (-t261 * t274 - t278 * t290) - g(3) * (t272 * t278 - t274 * t294) 0, 0, 0, 0, 0, -g(1) * t255 - g(2) * t253 - g(3) * t257, t282, 0, 0, 0, 0, 0, t245, t283, t245, -t282, -t283, -g(1) * (t280 * pkin(1) + t255 * pkin(4) + t249 * pkin(5) + t254 * pkin(11) + t248 * qJ(6) - t262 * t281 + t263 * t267) - g(2) * (t276 * pkin(1) + t253 * pkin(4) + t247 * pkin(5) + t252 * pkin(11) + t246 * qJ(6) - t260 * t281 + t261 * t267) - g(3) * (pkin(3) * t289 + t257 * pkin(4) + t251 * pkin(5) + t272 * pkin(8) + t256 * pkin(11) + t250 * qJ(6) + pkin(7)) + (-g(3) * (t267 * t275 + t279 * t281) + t297 * (pkin(3) * t274 + pkin(8))) * t271;];
U_reg  = t1;
