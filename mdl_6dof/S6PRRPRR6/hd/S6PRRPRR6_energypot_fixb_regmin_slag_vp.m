% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:37
% EndTime: 2019-03-08 22:27:37
% DurationCPUTime: 0.17s
% Computational Cost: add. (259->77), mult. (636->135), div. (0->0), fcn. (824->16), ass. (0->47)
t290 = sin(pkin(12));
t292 = sin(pkin(6));
t315 = t290 * t292;
t291 = sin(pkin(7));
t296 = cos(pkin(6));
t314 = t291 * t296;
t294 = cos(pkin(12));
t313 = t292 * t294;
t295 = cos(pkin(7));
t312 = t292 * t295;
t299 = sin(qJ(2));
t311 = t292 * t299;
t301 = cos(qJ(3));
t310 = t292 * t301;
t302 = cos(qJ(2));
t309 = t292 * t302;
t308 = t295 * t301;
t307 = t295 * t302;
t306 = t296 * t299;
t305 = t296 * t302;
t304 = t291 * t310;
t279 = -t290 * t299 + t294 * t305;
t274 = -t279 * t291 - t294 * t312;
t281 = -t290 * t305 - t294 * t299;
t275 = -t281 * t291 + t290 * t312;
t278 = -t291 * t309 + t295 * t296;
t280 = t290 * t302 + t294 * t306;
t298 = sin(qJ(3));
t268 = -t279 * t308 + t280 * t298 + t294 * t304;
t282 = -t290 * t306 + t294 * t302;
t270 = -t281 * t308 + t282 * t298 - t290 * t304;
t272 = t298 * t311 - t301 * t314 - t307 * t310;
t303 = g(1) * t270 + g(2) * t268 + g(3) * t272;
t300 = cos(qJ(6));
t297 = sin(qJ(6));
t293 = cos(pkin(13));
t289 = sin(pkin(13));
t288 = pkin(13) + qJ(5);
t287 = cos(t288);
t286 = sin(t288);
t273 = t298 * t314 + (t298 * t307 + t299 * t301) * t292;
t271 = t282 * t301 + (t281 * t295 + t291 * t315) * t298;
t269 = t280 * t301 + (t279 * t295 - t291 * t313) * t298;
t267 = t273 * t287 + t278 * t286;
t266 = t271 * t287 + t275 * t286;
t265 = t269 * t287 + t274 * t286;
t1 = [-g(3) * qJ(1), 0, -g(1) * t282 - g(2) * t280 - g(3) * t311, -g(1) * t281 - g(2) * t279 - g(3) * t309, 0, 0, 0, 0, 0, -g(1) * t271 - g(2) * t269 - g(3) * t273, t303, -g(1) * (t271 * t293 + t275 * t289) - g(2) * (t269 * t293 + t274 * t289) - g(3) * (t273 * t293 + t278 * t289) -g(1) * (-t271 * t289 + t275 * t293) - g(2) * (-t269 * t289 + t274 * t293) - g(3) * (-t273 * t289 + t278 * t293) -t303, -g(1) * (t294 * pkin(1) + t282 * pkin(2) + t271 * pkin(3) + pkin(8) * t315 + t270 * qJ(4)) - g(2) * (t290 * pkin(1) + t280 * pkin(2) + t269 * pkin(3) - pkin(8) * t313 + t268 * qJ(4)) - g(3) * (pkin(2) * t311 + t273 * pkin(3) + t296 * pkin(8) + t272 * qJ(4) + qJ(1)) + (-g(1) * t275 - g(2) * t274 - g(3) * t278) * pkin(9), 0, 0, 0, 0, 0, -g(1) * t266 - g(2) * t265 - g(3) * t267, -g(1) * (-t271 * t286 + t275 * t287) - g(2) * (-t269 * t286 + t274 * t287) - g(3) * (-t273 * t286 + t278 * t287) 0, 0, 0, 0, 0, -g(1) * (t266 * t300 + t270 * t297) - g(2) * (t265 * t300 + t268 * t297) - g(3) * (t267 * t300 + t272 * t297) -g(1) * (-t266 * t297 + t270 * t300) - g(2) * (-t265 * t297 + t268 * t300) - g(3) * (-t267 * t297 + t272 * t300);];
U_reg  = t1;
