% Calculate minimal parameter regressor of potential energy for
% S6PRRRRR4
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:00
% EndTime: 2019-03-09 01:00:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (203->57), mult. (493->109), div. (0->0), fcn. (650->16), ass. (0->43)
t288 = sin(pkin(7));
t289 = sin(pkin(6));
t309 = t288 * t289;
t292 = cos(pkin(6));
t308 = t288 * t292;
t291 = cos(pkin(7));
t307 = t289 * t291;
t300 = cos(qJ(2));
t306 = t289 * t300;
t305 = t291 * t300;
t296 = sin(qJ(2));
t304 = t292 * t296;
t303 = t292 * t300;
t287 = sin(pkin(13));
t290 = cos(pkin(13));
t280 = -t287 * t296 + t290 * t303;
t302 = -t280 * t291 + t290 * t309;
t282 = -t287 * t303 - t290 * t296;
t301 = t282 * t291 + t287 * t309;
t299 = cos(qJ(3));
t298 = cos(qJ(4));
t297 = cos(qJ(6));
t295 = sin(qJ(3));
t294 = sin(qJ(4));
t293 = sin(qJ(6));
t286 = qJ(4) + qJ(5);
t285 = cos(t286);
t284 = sin(t286);
t283 = -t287 * t304 + t290 * t300;
t281 = t287 * t300 + t290 * t304;
t279 = -t288 * t306 + t292 * t291;
t278 = -t282 * t288 + t287 * t307;
t277 = -t280 * t288 - t290 * t307;
t276 = t295 * t308 + (t295 * t305 + t296 * t299) * t289;
t275 = -t299 * t308 + (t295 * t296 - t299 * t305) * t289;
t274 = t283 * t299 + t295 * t301;
t273 = t283 * t295 - t299 * t301;
t272 = t281 * t299 - t295 * t302;
t271 = t281 * t295 + t299 * t302;
t270 = t276 * t285 + t279 * t284;
t269 = t274 * t285 + t278 * t284;
t268 = t272 * t285 + t277 * t284;
t1 = [-g(3) * qJ(1), 0, -g(3) * t289 * t296 - g(1) * t283 - g(2) * t281, -g(1) * t282 - g(2) * t280 - g(3) * t306, 0, 0, 0, 0, 0, -g(1) * t274 - g(2) * t272 - g(3) * t276, g(1) * t273 + g(2) * t271 + g(3) * t275, 0, 0, 0, 0, 0, -g(1) * (t274 * t298 + t278 * t294) - g(2) * (t272 * t298 + t277 * t294) - g(3) * (t276 * t298 + t279 * t294) -g(1) * (-t274 * t294 + t278 * t298) - g(2) * (-t272 * t294 + t277 * t298) - g(3) * (-t276 * t294 + t279 * t298) 0, 0, 0, 0, 0, -g(1) * t269 - g(2) * t268 - g(3) * t270, -g(1) * (-t274 * t284 + t278 * t285) - g(2) * (-t272 * t284 + t277 * t285) - g(3) * (-t276 * t284 + t279 * t285) 0, 0, 0, 0, 0, -g(1) * (t269 * t297 + t273 * t293) - g(2) * (t268 * t297 + t271 * t293) - g(3) * (t270 * t297 + t275 * t293) -g(1) * (-t269 * t293 + t273 * t297) - g(2) * (-t268 * t293 + t271 * t297) - g(3) * (-t270 * t293 + t275 * t297);];
U_reg  = t1;
