% Calculate minimal parameter regressor of potential energy for
% S6PRRRRP5
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
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:53
% EndTime: 2019-03-09 00:25:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (264->77), mult. (697->127), div. (0->0), fcn. (903->14), ass. (0->52)
t270 = sin(pkin(12));
t273 = cos(pkin(12));
t280 = sin(qJ(2));
t275 = cos(pkin(6));
t284 = cos(qJ(2));
t290 = t275 * t284;
t262 = -t270 * t280 + t273 * t290;
t291 = t275 * t280;
t263 = t270 * t284 + t273 * t291;
t279 = sin(qJ(3));
t271 = sin(pkin(7));
t272 = sin(pkin(6));
t283 = cos(qJ(3));
t295 = t272 * t283;
t289 = t271 * t295;
t274 = cos(pkin(7));
t293 = t274 * t283;
t249 = -t262 * t293 + t263 * t279 + t273 * t289;
t277 = sin(qJ(5));
t303 = t249 * t277;
t264 = -t270 * t290 - t273 * t280;
t265 = -t270 * t291 + t273 * t284;
t251 = -t264 * t293 + t265 * t279 - t270 * t289;
t302 = t251 * t277;
t292 = t274 * t284;
t296 = t272 * t280;
t299 = t271 * t275;
t257 = t279 * t296 - t283 * t299 - t292 * t295;
t301 = t257 * t277;
t300 = t270 * t272;
t298 = t272 * t273;
t297 = t272 * t274;
t294 = t272 * t284;
t288 = t262 * t271 + t273 * t297;
t287 = -t264 * t271 + t270 * t297;
t286 = t271 * t294 - t275 * t274;
t250 = t263 * t283 + (t262 * t274 - t271 * t298) * t279;
t278 = sin(qJ(4));
t282 = cos(qJ(4));
t245 = t250 * t278 + t288 * t282;
t252 = t265 * t283 + (t264 * t274 + t271 * t300) * t279;
t247 = t252 * t278 - t287 * t282;
t258 = t279 * t299 + (t279 * t292 + t280 * t283) * t272;
t253 = t258 * t278 + t286 * t282;
t285 = g(1) * t247 + g(2) * t245 + g(3) * t253;
t281 = cos(qJ(5));
t276 = -qJ(6) - pkin(11);
t269 = t281 * pkin(5) + pkin(4);
t254 = t258 * t282 - t286 * t278;
t248 = t252 * t282 + t287 * t278;
t246 = t250 * t282 - t288 * t278;
t1 = [-g(3) * qJ(1), 0, -g(1) * t265 - g(2) * t263 - g(3) * t296, -g(1) * t264 - g(2) * t262 - g(3) * t294, 0, 0, 0, 0, 0, -g(1) * t252 - g(2) * t250 - g(3) * t258, g(1) * t251 + g(2) * t249 + g(3) * t257, 0, 0, 0, 0, 0, -g(1) * t248 - g(2) * t246 - g(3) * t254, t285, 0, 0, 0, 0, 0, -g(1) * (t248 * t281 + t302) - g(2) * (t246 * t281 + t303) - g(3) * (t254 * t281 + t301) -g(1) * (-t248 * t277 + t251 * t281) - g(2) * (-t246 * t277 + t249 * t281) - g(3) * (-t254 * t277 + t257 * t281) -t285, -g(1) * (t273 * pkin(1) + t265 * pkin(2) + t252 * pkin(3) + pkin(5) * t302 + pkin(8) * t300 + t251 * pkin(10) - t247 * t276 + t248 * t269) - g(2) * (t270 * pkin(1) + t263 * pkin(2) + t250 * pkin(3) + pkin(5) * t303 - pkin(8) * t298 + t249 * pkin(10) - t245 * t276 + t246 * t269) - g(3) * (pkin(2) * t296 + t258 * pkin(3) + pkin(5) * t301 + t275 * pkin(8) + t257 * pkin(10) - t253 * t276 + t254 * t269 + qJ(1)) + (-g(1) * t287 + g(2) * t288 + g(3) * t286) * pkin(9);];
U_reg  = t1;
