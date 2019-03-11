% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:34:08
% EndTime: 2019-03-09 22:34:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (141->62), mult. (215->102), div. (0->0), fcn. (262->14), ass. (0->38)
t273 = qJ(3) + qJ(4);
t270 = sin(t273);
t277 = sin(qJ(3));
t298 = t277 * pkin(3) + pkin(4) * t270 + pkin(8);
t279 = sin(qJ(1));
t283 = cos(qJ(1));
t297 = -g(1) * t279 + g(2) * t283;
t274 = sin(pkin(6));
t278 = sin(qJ(2));
t293 = t274 * t278;
t292 = t274 * t279;
t281 = cos(qJ(3));
t291 = t274 * t281;
t282 = cos(qJ(2));
t290 = t274 * t282;
t289 = t274 * t283;
t288 = t279 * t278;
t287 = t279 * t282;
t286 = t283 * t278;
t285 = t283 * t282;
t275 = cos(pkin(6));
t260 = -t275 * t285 + t288;
t262 = t275 * t287 + t286;
t284 = -g(1) * t262 - g(2) * t260 + g(3) * t290;
t280 = cos(qJ(6));
t276 = sin(qJ(6));
t272 = -qJ(5) - pkin(10) - pkin(9);
t271 = cos(t273);
t269 = pkin(12) + t273;
t268 = cos(t269);
t267 = sin(t269);
t264 = t281 * pkin(3) + pkin(4) * t271 + pkin(2);
t263 = -t275 * t288 + t285;
t261 = t275 * t286 + t287;
t259 = t275 * t267 + t268 * t293;
t258 = t263 * t268 + t267 * t292;
t257 = t261 * t268 - t267 * t289;
t1 = [0, -g(1) * t283 - g(2) * t279, -t297, 0, 0, 0, 0, 0, -g(1) * t263 - g(2) * t261 - g(3) * t293, -t284, 0, 0, 0, 0, 0, -g(1) * (t263 * t281 + t277 * t292) - g(2) * (t261 * t281 - t277 * t289) - g(3) * (t275 * t277 + t278 * t291) -g(1) * (-t263 * t277 + t279 * t291) - g(2) * (-t261 * t277 - t281 * t289) - g(3) * (t275 * t281 - t277 * t293) 0, 0, 0, 0, 0, -g(1) * (t263 * t271 + t270 * t292) - g(2) * (t261 * t271 - t270 * t289) - g(3) * (t275 * t270 + t271 * t293) -g(1) * (-t263 * t270 + t271 * t292) - g(2) * (-t261 * t270 - t271 * t289) - g(3) * (-t270 * t293 + t275 * t271) t284, -g(1) * (t283 * pkin(1) - t262 * t272 + t263 * t264) - g(2) * (t279 * pkin(1) - t260 * t272 + t261 * t264) - g(3) * (t298 * t275 + pkin(7)) + (-g(3) * (t264 * t278 + t272 * t282) + t297 * t298) * t274, 0, 0, 0, 0, 0, -g(1) * (t258 * t280 + t262 * t276) - g(2) * (t257 * t280 + t260 * t276) - g(3) * (t259 * t280 - t276 * t290) -g(1) * (-t258 * t276 + t262 * t280) - g(2) * (-t257 * t276 + t260 * t280) - g(3) * (-t259 * t276 - t280 * t290);];
U_reg  = t1;
