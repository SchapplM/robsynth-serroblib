% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:03
% EndTime: 2019-03-10 00:04:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (186->65), mult. (460->104), div. (0->0), fcn. (596->12), ass. (0->41)
t266 = sin(pkin(6));
t271 = sin(qJ(2));
t288 = t266 * t271;
t272 = sin(qJ(1));
t287 = t266 * t272;
t275 = cos(qJ(3));
t286 = t266 * t275;
t276 = cos(qJ(2));
t285 = t266 * t276;
t277 = cos(qJ(1));
t284 = t266 * t277;
t283 = t272 * t271;
t282 = t272 * t276;
t281 = t277 * t271;
t280 = t277 * t276;
t267 = cos(pkin(6));
t260 = t267 * t281 + t282;
t270 = sin(qJ(3));
t252 = t260 * t275 - t270 * t284;
t259 = -t267 * t280 + t283;
t269 = sin(qJ(4));
t274 = cos(qJ(4));
t245 = t252 * t269 - t259 * t274;
t262 = -t267 * t283 + t280;
t254 = t262 * t275 + t270 * t287;
t261 = t267 * t282 + t281;
t247 = t254 * t269 - t261 * t274;
t258 = t267 * t270 + t271 * t286;
t249 = t258 * t269 + t274 * t285;
t279 = g(1) * t247 + g(2) * t245 + g(3) * t249;
t251 = t260 * t270 + t275 * t284;
t253 = t262 * t270 - t272 * t286;
t257 = -t267 * t275 + t270 * t288;
t278 = g(1) * t253 + g(2) * t251 + g(3) * t257;
t273 = cos(qJ(6));
t268 = sin(qJ(6));
t250 = t258 * t274 - t269 * t285;
t248 = t254 * t274 + t261 * t269;
t246 = t252 * t274 + t259 * t269;
t244 = -g(1) * t248 - g(2) * t246 - g(3) * t250;
t1 = [0, -g(1) * t277 - g(2) * t272, g(1) * t272 - g(2) * t277, 0, 0, 0, 0, 0, -g(1) * t262 - g(2) * t260 - g(3) * t288, g(1) * t261 + g(2) * t259 - g(3) * t285, 0, 0, 0, 0, 0, -g(1) * t254 - g(2) * t252 - g(3) * t258, t278, 0, 0, 0, 0, 0, t244, t279, t244, -t278, -t279, -g(1) * (t277 * pkin(1) + t262 * pkin(2) + t254 * pkin(3) + t248 * pkin(4) + pkin(8) * t287 + t261 * pkin(9) + t253 * pkin(10) + t247 * qJ(5)) - g(2) * (t272 * pkin(1) + t260 * pkin(2) + t252 * pkin(3) + t246 * pkin(4) - pkin(8) * t284 + t259 * pkin(9) + t251 * pkin(10) + t245 * qJ(5)) - g(3) * (t258 * pkin(3) + t250 * pkin(4) + t267 * pkin(8) + t257 * pkin(10) + t249 * qJ(5) + pkin(7) + (pkin(2) * t271 - pkin(9) * t276) * t266) 0, 0, 0, 0, 0, -g(1) * (t247 * t268 + t248 * t273) - g(2) * (t245 * t268 + t246 * t273) - g(3) * (t249 * t268 + t250 * t273) -g(1) * (t247 * t273 - t248 * t268) - g(2) * (t245 * t273 - t246 * t268) - g(3) * (t249 * t273 - t250 * t268);];
U_reg  = t1;
