% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:37
% EndTime: 2019-03-08 22:40:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (214->67), mult. (578->117), div. (0->0), fcn. (745->14), ass. (0->46)
t267 = sin(pkin(12));
t270 = cos(pkin(12));
t276 = sin(qJ(2));
t272 = cos(pkin(6));
t280 = cos(qJ(2));
t284 = t272 * t280;
t260 = -t267 * t276 + t270 * t284;
t268 = sin(pkin(7));
t269 = sin(pkin(6));
t271 = cos(pkin(7));
t291 = t269 * t271;
t255 = -t260 * t268 - t270 * t291;
t262 = -t267 * t284 - t270 * t276;
t256 = -t262 * t268 + t267 * t291;
t288 = t269 * t280;
t259 = -t268 * t288 + t272 * t271;
t298 = -g(1) * t256 - g(2) * t255 - g(3) * t259;
t294 = t267 * t269;
t293 = t268 * t272;
t292 = t269 * t270;
t290 = t269 * t276;
t279 = cos(qJ(3));
t289 = t269 * t279;
t287 = t271 * t279;
t286 = t271 * t280;
t285 = t272 * t276;
t283 = t268 * t289;
t261 = t267 * t280 + t270 * t285;
t275 = sin(qJ(3));
t248 = -t260 * t287 + t261 * t275 + t270 * t283;
t263 = -t267 * t285 + t270 * t280;
t250 = -t262 * t287 + t263 * t275 - t267 * t283;
t253 = t275 * t290 - t279 * t293 - t286 * t289;
t282 = g(1) * t250 + g(2) * t248 + g(3) * t253;
t249 = t261 * t279 + (t260 * t271 - t268 * t292) * t275;
t251 = t263 * t279 + (t262 * t271 + t268 * t294) * t275;
t254 = t275 * t293 + (t275 * t286 + t276 * t279) * t269;
t281 = g(1) * t251 + g(2) * t249 + g(3) * t254;
t278 = cos(qJ(5));
t277 = cos(qJ(6));
t274 = sin(qJ(5));
t273 = sin(qJ(6));
t252 = t253 * t274 + t259 * t278;
t247 = t250 * t274 + t256 * t278;
t246 = t248 * t274 + t255 * t278;
t1 = [-g(3) * qJ(1), 0, -g(1) * t263 - g(2) * t261 - g(3) * t290, -g(1) * t262 - g(2) * t260 - g(3) * t288, 0, 0, 0, 0, 0, -t281, t282, t298, t281, -t282, -g(1) * (t270 * pkin(1) + t263 * pkin(2) + t251 * pkin(3) + pkin(8) * t294 + t250 * qJ(4)) - g(2) * (t267 * pkin(1) + t261 * pkin(2) + t249 * pkin(3) - pkin(8) * t292 + t248 * qJ(4)) - g(3) * (pkin(2) * t290 + t254 * pkin(3) + t272 * pkin(8) + t253 * qJ(4) + qJ(1)) + t298 * pkin(9), 0, 0, 0, 0, 0, -g(1) * t247 - g(2) * t246 - g(3) * t252, -g(1) * (t250 * t278 - t256 * t274) - g(2) * (t248 * t278 - t255 * t274) - g(3) * (t253 * t278 - t259 * t274) 0, 0, 0, 0, 0, -g(1) * (t247 * t277 + t251 * t273) - g(2) * (t246 * t277 + t249 * t273) - g(3) * (t252 * t277 + t254 * t273) -g(1) * (-t247 * t273 + t251 * t277) - g(2) * (-t246 * t273 + t249 * t277) - g(3) * (-t252 * t273 + t254 * t277);];
U_reg  = t1;
