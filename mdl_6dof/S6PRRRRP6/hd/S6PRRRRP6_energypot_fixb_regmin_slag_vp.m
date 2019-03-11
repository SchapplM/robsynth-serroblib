% Calculate minimal parameter regressor of potential energy for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:24
% EndTime: 2019-03-09 00:34:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (379->79), mult. (1025->129), div. (0->0), fcn. (1345->14), ass. (0->55)
t286 = sin(pkin(12));
t288 = sin(pkin(6));
t316 = t286 * t288;
t287 = sin(pkin(7));
t291 = cos(pkin(6));
t315 = t287 * t291;
t289 = cos(pkin(12));
t314 = t288 * t289;
t290 = cos(pkin(7));
t313 = t288 * t290;
t295 = sin(qJ(2));
t312 = t288 * t295;
t298 = cos(qJ(3));
t311 = t288 * t298;
t299 = cos(qJ(2));
t310 = t288 * t299;
t309 = t290 * t298;
t308 = t290 * t299;
t307 = t291 * t295;
t306 = t291 * t299;
t305 = t287 * t311;
t279 = -t286 * t295 + t289 * t306;
t304 = t279 * t287 + t289 * t313;
t281 = -t286 * t306 - t289 * t295;
t303 = -t281 * t287 + t286 * t313;
t302 = t287 * t310 - t291 * t290;
t280 = t286 * t299 + t289 * t307;
t294 = sin(qJ(3));
t266 = t280 * t298 + (t279 * t290 - t287 * t314) * t294;
t293 = sin(qJ(4));
t297 = cos(qJ(4));
t260 = t266 * t297 - t304 * t293;
t265 = -t279 * t309 + t280 * t294 + t289 * t305;
t292 = sin(qJ(5));
t296 = cos(qJ(5));
t253 = t260 * t292 - t265 * t296;
t282 = -t286 * t307 + t289 * t299;
t268 = t282 * t298 + (t281 * t290 + t287 * t316) * t294;
t262 = t268 * t297 + t303 * t293;
t267 = -t281 * t309 + t282 * t294 - t286 * t305;
t255 = t262 * t292 - t267 * t296;
t275 = t294 * t315 + (t294 * t308 + t295 * t298) * t288;
t270 = t275 * t297 - t302 * t293;
t274 = t294 * t312 - t298 * t315 - t308 * t311;
t257 = t270 * t292 - t274 * t296;
t301 = g(1) * t255 + g(2) * t253 + g(3) * t257;
t259 = t266 * t293 + t304 * t297;
t261 = t268 * t293 - t303 * t297;
t269 = t275 * t293 + t302 * t297;
t300 = g(1) * t261 + g(2) * t259 + g(3) * t269;
t258 = t270 * t296 + t274 * t292;
t256 = t262 * t296 + t267 * t292;
t254 = t260 * t296 + t265 * t292;
t252 = -g(1) * t256 - g(2) * t254 - g(3) * t258;
t1 = [-g(3) * qJ(1), 0, -g(1) * t282 - g(2) * t280 - g(3) * t312, -g(1) * t281 - g(2) * t279 - g(3) * t310, 0, 0, 0, 0, 0, -g(1) * t268 - g(2) * t266 - g(3) * t275, g(1) * t267 + g(2) * t265 + g(3) * t274, 0, 0, 0, 0, 0, -g(1) * t262 - g(2) * t260 - g(3) * t270, t300, 0, 0, 0, 0, 0, t252, t301, t252, -t300, -t301, -g(1) * (t289 * pkin(1) + t282 * pkin(2) + t268 * pkin(3) + t262 * pkin(4) + t256 * pkin(5) + pkin(8) * t316 + t267 * pkin(10) + t261 * pkin(11) + t255 * qJ(6)) - g(2) * (t286 * pkin(1) + t280 * pkin(2) + t266 * pkin(3) + t260 * pkin(4) + t254 * pkin(5) - pkin(8) * t314 + t265 * pkin(10) + t259 * pkin(11) + t253 * qJ(6)) - g(3) * (pkin(2) * t312 + t275 * pkin(3) + t270 * pkin(4) + t258 * pkin(5) + t291 * pkin(8) + t274 * pkin(10) + t269 * pkin(11) + t257 * qJ(6) + qJ(1)) + (-g(1) * t303 + g(2) * t304 + g(3) * t302) * pkin(9);];
U_reg  = t1;
