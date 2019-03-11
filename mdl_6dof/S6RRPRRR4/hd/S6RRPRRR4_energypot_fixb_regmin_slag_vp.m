% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:34:48
% EndTime: 2019-03-09 13:34:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (141->57), mult. (288->107), div. (0->0), fcn. (369->14), ass. (0->40)
t300 = pkin(8) + qJ(3);
t279 = sin(pkin(6));
t284 = sin(qJ(2));
t299 = t279 * t284;
t285 = sin(qJ(1));
t298 = t279 * t285;
t289 = cos(qJ(1));
t297 = t279 * t289;
t296 = t285 * t284;
t288 = cos(qJ(2));
t295 = t285 * t288;
t294 = t289 * t284;
t293 = t289 * t288;
t292 = g(1) * t285 - g(2) * t289;
t278 = sin(pkin(12));
t280 = cos(pkin(12));
t291 = t288 * t278 + t284 * t280;
t290 = t284 * t278 - t288 * t280;
t287 = cos(qJ(4));
t286 = cos(qJ(6));
t283 = sin(qJ(4));
t282 = sin(qJ(6));
t281 = cos(pkin(6));
t277 = qJ(4) + qJ(5);
t276 = cos(t277);
t275 = sin(t277);
t274 = t288 * pkin(2) + pkin(1);
t271 = t281 * t284 * pkin(2) - t279 * t300;
t270 = t291 * t281;
t269 = t290 * t281;
t268 = t291 * t279;
t267 = t290 * t279;
t266 = t268 * t276 + t281 * t275;
t265 = -t285 * t270 - t289 * t290;
t264 = -t285 * t269 + t289 * t291;
t263 = t289 * t270 - t285 * t290;
t262 = t289 * t269 + t285 * t291;
t261 = t265 * t276 + t275 * t298;
t260 = t263 * t276 - t275 * t297;
t1 = [0, -g(1) * t289 - g(2) * t285, t292, 0, 0, 0, 0, 0, -g(1) * (-t281 * t296 + t293) - g(2) * (t281 * t294 + t295) - g(3) * t299, -g(1) * (-t281 * t295 - t294) - g(2) * (t281 * t293 - t296) - g(3) * t279 * t288, -g(3) * t281 - t279 * t292, -g(1) * (-t285 * t271 + t289 * t274) - g(2) * (t289 * t271 + t285 * t274) - g(3) * (pkin(2) * t299 + t281 * t300 + pkin(7)) 0, 0, 0, 0, 0, -g(1) * (t265 * t287 + t283 * t298) - g(2) * (t263 * t287 - t283 * t297) - g(3) * (t268 * t287 + t281 * t283) -g(1) * (-t265 * t283 + t287 * t298) - g(2) * (-t263 * t283 - t287 * t297) - g(3) * (-t268 * t283 + t281 * t287) 0, 0, 0, 0, 0, -g(1) * t261 - g(2) * t260 - g(3) * t266, -g(1) * (-t265 * t275 + t276 * t298) - g(2) * (-t263 * t275 - t276 * t297) - g(3) * (-t268 * t275 + t281 * t276) 0, 0, 0, 0, 0, -g(1) * (t261 * t286 + t264 * t282) - g(2) * (t260 * t286 + t262 * t282) - g(3) * (t266 * t286 + t267 * t282) -g(1) * (-t261 * t282 + t264 * t286) - g(2) * (-t260 * t282 + t262 * t286) - g(3) * (-t266 * t282 + t267 * t286);];
U_reg  = t1;
