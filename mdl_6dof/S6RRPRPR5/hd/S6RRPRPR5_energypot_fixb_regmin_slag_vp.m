% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:10
% EndTime: 2019-03-09 10:35:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (220->72), mult. (500->122), div. (0->0), fcn. (642->14), ass. (0->47)
t312 = pkin(8) + qJ(3);
t287 = sin(pkin(6));
t292 = sin(qJ(2));
t311 = t287 * t292;
t293 = sin(qJ(1));
t310 = t287 * t293;
t296 = cos(qJ(1));
t309 = t287 * t296;
t308 = t293 * t292;
t295 = cos(qJ(2));
t307 = t293 * t295;
t306 = t296 * t292;
t305 = t296 * t295;
t290 = cos(pkin(6));
t271 = t290 * t292 * pkin(2) - t312 * t287;
t279 = t295 * pkin(2) + pkin(1);
t304 = t296 * t271 + t293 * t279;
t303 = -t293 * t271 + t296 * t279;
t302 = pkin(2) * t311 + t312 * t290 + pkin(7);
t301 = g(1) * t293 - g(2) * t296;
t286 = sin(pkin(11));
t289 = cos(pkin(11));
t300 = t295 * t286 + t292 * t289;
t299 = t292 * t286 - t295 * t289;
t298 = t299 * t290;
t270 = t300 * t290;
t261 = t296 * t270 - t293 * t299;
t291 = sin(qJ(4));
t294 = cos(qJ(4));
t256 = t261 * t291 + t294 * t309;
t263 = -t293 * t270 - t296 * t299;
t258 = t263 * t291 - t294 * t310;
t269 = t300 * t287;
t264 = t269 * t291 - t290 * t294;
t297 = g(1) * t258 + g(2) * t256 + g(3) * t264;
t288 = cos(pkin(12));
t285 = sin(pkin(12));
t284 = pkin(12) + qJ(6);
t281 = cos(t284);
t280 = sin(t284);
t268 = t299 * t287;
t265 = t269 * t294 + t290 * t291;
t262 = t293 * t298 - t296 * t300;
t260 = -t293 * t300 - t296 * t298;
t259 = t263 * t294 + t291 * t310;
t257 = t261 * t294 - t291 * t309;
t1 = [0, -g(1) * t296 - g(2) * t293, t301, 0, 0, 0, 0, 0, -g(1) * (-t290 * t308 + t305) - g(2) * (t290 * t306 + t307) - g(3) * t311, -g(1) * (-t290 * t307 - t306) - g(2) * (t290 * t305 - t308) - g(3) * t287 * t295, -g(3) * t290 - t301 * t287, -g(1) * t303 - g(2) * t304 - g(3) * t302, 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t257 - g(3) * t265, t297, -g(1) * (t259 * t288 - t262 * t285) - g(2) * (t257 * t288 - t260 * t285) - g(3) * (t265 * t288 + t268 * t285) -g(1) * (-t259 * t285 - t262 * t288) - g(2) * (-t257 * t285 - t260 * t288) - g(3) * (-t265 * t285 + t268 * t288) -t297, -g(1) * (t263 * pkin(3) + t259 * pkin(4) - t262 * pkin(9) + t258 * qJ(5) + t303) - g(2) * (t261 * pkin(3) + t257 * pkin(4) - t260 * pkin(9) + t256 * qJ(5) + t304) - g(3) * (t269 * pkin(3) + t265 * pkin(4) + t268 * pkin(9) + t264 * qJ(5) + t302) 0, 0, 0, 0, 0, -g(1) * (t259 * t281 - t262 * t280) - g(2) * (t257 * t281 - t260 * t280) - g(3) * (t265 * t281 + t268 * t280) -g(1) * (-t259 * t280 - t262 * t281) - g(2) * (-t257 * t280 - t260 * t281) - g(3) * (-t265 * t280 + t268 * t281);];
U_reg  = t1;
