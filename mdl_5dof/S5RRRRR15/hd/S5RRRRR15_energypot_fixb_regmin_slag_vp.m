% Calculate minimal parameter regressor of potential energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR15_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:18
% EndTime: 2024-09-27 22:26:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (158->64), mult. (144->106), div. (0->0), fcn. (168->22), ass. (0->53)
t283 = sin(qJ(1));
t306 = -t283 / 0.2e1;
t286 = cos(qJ(1));
t305 = t286 / 0.2e1;
t278 = sin(pkin(5));
t304 = g(3) * t278;
t276 = qJ(2) + qJ(3);
t275 = qJ(4) + t276;
t270 = cos(t275);
t279 = cos(pkin(6));
t303 = t270 * t279;
t277 = sin(pkin(6));
t302 = t277 * t278;
t280 = cos(pkin(5));
t301 = t277 * t280;
t281 = sin(qJ(5));
t300 = t283 * t281;
t282 = sin(qJ(2));
t299 = t283 * t282;
t285 = cos(qJ(2));
t298 = t283 * t285;
t284 = cos(qJ(5));
t297 = t284 * t283;
t296 = t286 * t281;
t295 = t286 * t282;
t294 = t286 * t284;
t293 = t286 * t285;
t292 = t280 * t296;
t291 = t283 * t302;
t290 = t280 * t300;
t289 = t280 * t297;
t288 = t286 * t302;
t287 = t280 * t294;
t272 = pkin(5) - t276;
t271 = pkin(5) + t276;
t274 = cos(t276);
t273 = sin(t276);
t269 = sin(t275);
t268 = -qJ(4) + t272;
t267 = qJ(4) + t271;
t266 = cos(t272);
t265 = cos(t271);
t264 = sin(t272);
t263 = sin(t271);
t262 = cos(t267);
t261 = sin(t268);
t260 = cos(t268) / 0.2e1;
t259 = sin(t267) / 0.2e1;
t258 = t266 + t265;
t257 = t263 - t264;
t256 = t260 + t262 / 0.2e1;
t255 = t259 - t261 / 0.2e1;
t1 = [0, -g(1) * t286 - g(2) * t283, g(1) * t283 - g(2) * t286, 0, 0, 0, 0, 0, -g(1) * (-t280 * t299 + t293) - g(2) * (t280 * t295 + t298) - t282 * t304, -g(1) * (-t280 * t298 - t295) - g(2) * (t280 * t293 - t299) - t285 * t304, 0, 0, 0, 0, 0, -g(1) * (t257 * t306 + t286 * t274) - g(2) * (t257 * t305 + t283 * t274) - g(3) * (t266 / 0.2e1 - t265 / 0.2e1), -g(1) * (t258 * t306 - t286 * t273) - g(2) * (t258 * t305 - t283 * t273) - g(3) * (t263 / 0.2e1 + t264 / 0.2e1), 0, 0, 0, 0, 0, -g(1) * (-t283 * t255 + t286 * t270) - g(2) * (t286 * t255 + t283 * t270) - g(3) * (t260 - t262 / 0.2e1), -g(1) * (-t283 * t256 - t286 * t269) - g(2) * (t286 * t256 - t283 * t269) - g(3) * (t259 + t261 / 0.2e1), 0, 0, 0, 0, 0, -g(1) * ((-t279 * t290 + t294) * t270 + (-t279 * t296 - t289) * t269 + t281 * t291) - g(2) * ((t279 * t292 + t297) * t270 + (-t279 * t300 + t287) * t269 - t281 * t288) - g(3) * (t281 * t301 + (t269 * t284 + t281 * t303) * t278), -g(1) * ((-t279 * t289 - t296) * t270 + (-t279 * t294 + t290) * t269 + t284 * t291) - g(2) * ((t279 * t287 - t300) * t270 + (-t279 * t297 - t292) * t269 - t284 * t288) - g(3) * (t284 * t301 + (-t269 * t281 + t284 * t303) * t278);];
U_reg = t1;
