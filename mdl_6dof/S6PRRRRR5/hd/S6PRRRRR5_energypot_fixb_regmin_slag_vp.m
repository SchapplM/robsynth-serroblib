% Calculate minimal parameter regressor of potential energy for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:01
% EndTime: 2019-03-09 01:10:01
% DurationCPUTime: 0.15s
% Computational Cost: add. (217->57), mult. (567->109), div. (0->0), fcn. (750->16), ass. (0->43)
t283 = sin(pkin(7));
t284 = sin(pkin(6));
t304 = t283 * t284;
t287 = cos(pkin(6));
t303 = t283 * t287;
t286 = cos(pkin(7));
t302 = t284 * t286;
t295 = cos(qJ(2));
t301 = t284 * t295;
t300 = t286 * t295;
t291 = sin(qJ(2));
t299 = t287 * t291;
t298 = t287 * t295;
t282 = sin(pkin(13));
t285 = cos(pkin(13));
t275 = -t282 * t291 + t285 * t298;
t297 = -t275 * t286 + t285 * t304;
t277 = -t282 * t298 - t285 * t291;
t296 = t277 * t286 + t282 * t304;
t294 = cos(qJ(3));
t293 = cos(qJ(4));
t292 = cos(qJ(5));
t290 = sin(qJ(3));
t289 = sin(qJ(4));
t288 = sin(qJ(5));
t281 = qJ(5) + qJ(6);
t280 = cos(t281);
t279 = sin(t281);
t278 = -t282 * t299 + t285 * t295;
t276 = t282 * t295 + t285 * t299;
t274 = -t283 * t301 + t287 * t286;
t273 = -t277 * t283 + t282 * t302;
t272 = -t275 * t283 - t285 * t302;
t271 = t290 * t303 + (t290 * t300 + t291 * t294) * t284;
t270 = -t294 * t303 + (t290 * t291 - t294 * t300) * t284;
t269 = t271 * t293 + t274 * t289;
t268 = t278 * t294 + t290 * t296;
t267 = t278 * t290 - t294 * t296;
t266 = t276 * t294 - t290 * t297;
t265 = t276 * t290 + t294 * t297;
t264 = t268 * t293 + t273 * t289;
t263 = t266 * t293 + t272 * t289;
t1 = [-g(3) * qJ(1), 0, -g(3) * t284 * t291 - g(1) * t278 - g(2) * t276, -g(1) * t277 - g(2) * t275 - g(3) * t301, 0, 0, 0, 0, 0, -g(1) * t268 - g(2) * t266 - g(3) * t271, g(1) * t267 + g(2) * t265 + g(3) * t270, 0, 0, 0, 0, 0, -g(1) * t264 - g(2) * t263 - g(3) * t269, -g(1) * (-t268 * t289 + t273 * t293) - g(2) * (-t266 * t289 + t272 * t293) - g(3) * (-t271 * t289 + t274 * t293) 0, 0, 0, 0, 0, -g(1) * (t264 * t292 + t267 * t288) - g(2) * (t263 * t292 + t265 * t288) - g(3) * (t269 * t292 + t270 * t288) -g(1) * (-t264 * t288 + t267 * t292) - g(2) * (-t263 * t288 + t265 * t292) - g(3) * (-t269 * t288 + t270 * t292) 0, 0, 0, 0, 0, -g(1) * (t264 * t280 + t267 * t279) - g(2) * (t263 * t280 + t265 * t279) - g(3) * (t269 * t280 + t270 * t279) -g(1) * (-t264 * t279 + t267 * t280) - g(2) * (-t263 * t279 + t265 * t280) - g(3) * (-t269 * t279 + t270 * t280);];
U_reg  = t1;
