% Calculate minimal parameter regressor of potential energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:27
% EndTime: 2019-03-08 18:47:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (320->81), mult. (831->138), div. (0->0), fcn. (1082->16), ass. (0->51)
t266 = sin(pkin(12));
t269 = sin(pkin(6));
t295 = t266 * t269;
t267 = sin(pkin(11));
t294 = t267 * t269;
t274 = cos(pkin(6));
t293 = t267 * t274;
t268 = sin(pkin(7));
t278 = cos(qJ(3));
t292 = t268 * t278;
t271 = cos(pkin(12));
t291 = t269 * t271;
t272 = cos(pkin(11));
t290 = t269 * t272;
t273 = cos(pkin(7));
t289 = t269 * t273;
t288 = t272 * t274;
t287 = t273 * t278;
t286 = t274 * qJ(2) + qJ(1);
t285 = t272 * pkin(1) + qJ(2) * t294;
t284 = t269 * t292;
t283 = t267 * pkin(1) - qJ(2) * t290;
t251 = -t267 * t266 + t271 * t288;
t282 = t251 * t268 + t272 * t289;
t253 = -t272 * t266 - t271 * t293;
t281 = -t253 * t268 + t267 * t289;
t280 = t268 * t291 - t274 * t273;
t252 = t266 * t288 + t267 * t271;
t276 = sin(qJ(3));
t239 = t252 * t278 + (t251 * t273 - t268 * t290) * t276;
t275 = sin(qJ(4));
t277 = cos(qJ(4));
t234 = t239 * t275 + t282 * t277;
t254 = -t266 * t293 + t272 * t271;
t241 = t254 * t278 + (t253 * t273 + t268 * t294) * t276;
t236 = t241 * t275 - t281 * t277;
t247 = t274 * t268 * t276 + (t271 * t273 * t276 + t266 * t278) * t269;
t242 = t247 * t275 + t280 * t277;
t279 = g(1) * t236 + g(2) * t234 + g(3) * t242;
t270 = cos(pkin(13));
t265 = sin(pkin(13));
t264 = pkin(13) + qJ(6);
t260 = cos(t264);
t259 = sin(t264);
t246 = -t274 * t292 + t276 * t295 - t287 * t291;
t243 = t247 * t277 - t280 * t275;
t240 = -t253 * t287 + t254 * t276 - t267 * t284;
t238 = -t251 * t287 + t252 * t276 + t272 * t284;
t237 = t241 * t277 + t281 * t275;
t235 = t239 * t277 - t282 * t275;
t1 = [-g(3) * qJ(1), -g(1) * t285 - g(2) * t283 - g(3) * t286, 0, -g(1) * t241 - g(2) * t239 - g(3) * t247, g(1) * t240 + g(2) * t238 + g(3) * t246, 0, 0, 0, 0, 0, -g(1) * t237 - g(2) * t235 - g(3) * t243, t279, -g(1) * (t237 * t270 + t240 * t265) - g(2) * (t235 * t270 + t238 * t265) - g(3) * (t243 * t270 + t246 * t265) -g(1) * (-t237 * t265 + t240 * t270) - g(2) * (-t235 * t265 + t238 * t270) - g(3) * (-t243 * t265 + t246 * t270) -t279, -g(1) * (t254 * pkin(2) + t241 * pkin(3) + t237 * pkin(4) + t240 * pkin(9) + t236 * qJ(5) + t285) - g(2) * (t252 * pkin(2) + t239 * pkin(3) + t235 * pkin(4) + t238 * pkin(9) + t234 * qJ(5) + t283) - g(3) * (pkin(2) * t295 + t247 * pkin(3) + t243 * pkin(4) + t246 * pkin(9) + t242 * qJ(5) + t286) + (-g(1) * t281 + g(2) * t282 + g(3) * t280) * pkin(8), 0, 0, 0, 0, 0, -g(1) * (t237 * t260 + t240 * t259) - g(2) * (t235 * t260 + t238 * t259) - g(3) * (t243 * t260 + t246 * t259) -g(1) * (-t237 * t259 + t240 * t260) - g(2) * (-t235 * t259 + t238 * t260) - g(3) * (-t243 * t259 + t246 * t260);];
U_reg  = t1;
