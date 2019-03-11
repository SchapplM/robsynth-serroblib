% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:44:37
% EndTime: 2019-03-08 23:44:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (323->83), mult. (841->141), div. (0->0), fcn. (1099->16), ass. (0->50)
t289 = sin(pkin(12));
t291 = sin(pkin(6));
t317 = t289 * t291;
t290 = sin(pkin(7));
t295 = cos(pkin(6));
t316 = t290 * t295;
t293 = cos(pkin(12));
t315 = t291 * t293;
t294 = cos(pkin(7));
t314 = t291 * t294;
t298 = sin(qJ(2));
t313 = t291 * t298;
t300 = cos(qJ(3));
t312 = t291 * t300;
t301 = cos(qJ(2));
t311 = t291 * t301;
t310 = t294 * t300;
t309 = t294 * t301;
t308 = t295 * t298;
t307 = t295 * t301;
t306 = t290 * t312;
t278 = -t289 * t298 + t293 * t307;
t305 = t278 * t290 + t293 * t314;
t280 = -t289 * t307 - t293 * t298;
t304 = -t280 * t290 + t289 * t314;
t303 = t290 * t311 - t295 * t294;
t279 = t289 * t301 + t293 * t308;
t297 = sin(qJ(3));
t266 = t279 * t300 + (t278 * t294 - t290 * t315) * t297;
t296 = sin(qJ(4));
t299 = cos(qJ(4));
t261 = t266 * t296 + t305 * t299;
t281 = -t289 * t308 + t293 * t301;
t268 = t281 * t300 + (t280 * t294 + t290 * t317) * t297;
t263 = t268 * t296 - t304 * t299;
t274 = t297 * t316 + (t297 * t309 + t298 * t300) * t291;
t269 = t274 * t296 + t303 * t299;
t302 = g(1) * t263 + g(2) * t261 + g(3) * t269;
t292 = cos(pkin(13));
t288 = sin(pkin(13));
t287 = pkin(13) + qJ(6);
t286 = cos(t287);
t285 = sin(t287);
t273 = t297 * t313 - t300 * t316 - t309 * t312;
t270 = t274 * t299 - t303 * t296;
t267 = -t280 * t310 + t281 * t297 - t289 * t306;
t265 = -t278 * t310 + t279 * t297 + t293 * t306;
t264 = t268 * t299 + t304 * t296;
t262 = t266 * t299 - t305 * t296;
t1 = [-g(3) * qJ(1), 0, -g(1) * t281 - g(2) * t279 - g(3) * t313, -g(1) * t280 - g(2) * t278 - g(3) * t311, 0, 0, 0, 0, 0, -g(1) * t268 - g(2) * t266 - g(3) * t274, g(1) * t267 + g(2) * t265 + g(3) * t273, 0, 0, 0, 0, 0, -g(1) * t264 - g(2) * t262 - g(3) * t270, t302, -g(1) * (t264 * t292 + t267 * t288) - g(2) * (t262 * t292 + t265 * t288) - g(3) * (t270 * t292 + t273 * t288) -g(1) * (-t264 * t288 + t267 * t292) - g(2) * (-t262 * t288 + t265 * t292) - g(3) * (-t270 * t288 + t273 * t292) -t302, -g(1) * (t293 * pkin(1) + t281 * pkin(2) + t268 * pkin(3) + t264 * pkin(4) + pkin(8) * t317 + t267 * pkin(10) + t263 * qJ(5)) - g(2) * (t289 * pkin(1) + t279 * pkin(2) + t266 * pkin(3) + t262 * pkin(4) - pkin(8) * t315 + t265 * pkin(10) + t261 * qJ(5)) - g(3) * (pkin(2) * t313 + t274 * pkin(3) + t270 * pkin(4) + t295 * pkin(8) + t273 * pkin(10) + t269 * qJ(5) + qJ(1)) + (-g(1) * t304 + g(2) * t305 + g(3) * t303) * pkin(9), 0, 0, 0, 0, 0, -g(1) * (t264 * t286 + t267 * t285) - g(2) * (t262 * t286 + t265 * t285) - g(3) * (t270 * t286 + t273 * t285) -g(1) * (-t264 * t285 + t267 * t286) - g(2) * (-t262 * t285 + t265 * t286) - g(3) * (-t270 * t285 + t273 * t286);];
U_reg  = t1;
