% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:02
% EndTime: 2019-03-09 09:05:02
% DurationCPUTime: 0.21s
% Computational Cost: add. (146->58), mult. (344->104), div. (0->0), fcn. (430->12), ass. (0->41)
t286 = pkin(8) + qJ(3);
t261 = sin(pkin(6));
t266 = sin(qJ(2));
t285 = t261 * t266;
t267 = sin(qJ(1));
t284 = t261 * t267;
t271 = cos(qJ(1));
t283 = t261 * t271;
t282 = t267 * t266;
t270 = cos(qJ(2));
t281 = t267 * t270;
t280 = t271 * t266;
t279 = t271 * t270;
t263 = cos(pkin(6));
t251 = t263 * t266 * pkin(2) - t286 * t261;
t257 = t270 * pkin(2) + pkin(1);
t278 = t271 * t251 + t267 * t257;
t277 = -t267 * t251 + t271 * t257;
t276 = pkin(2) * t285 + t286 * t263 + pkin(7);
t275 = g(1) * t267 - g(2) * t271;
t260 = sin(pkin(11));
t262 = cos(pkin(11));
t274 = t270 * t260 + t266 * t262;
t273 = t266 * t260 - t270 * t262;
t272 = t273 * t263;
t269 = cos(qJ(5));
t268 = cos(qJ(6));
t265 = sin(qJ(5));
t264 = sin(qJ(6));
t250 = t274 * t263;
t249 = t274 * t261;
t248 = t273 * t261;
t247 = -g(3) * t263 - t275 * t261;
t244 = t248 * t265 + t263 * t269;
t243 = -t267 * t250 - t271 * t273;
t242 = t267 * t272 - t271 * t274;
t241 = t271 * t250 - t267 * t273;
t240 = -t267 * t274 - t271 * t272;
t239 = -t240 * t265 - t269 * t283;
t238 = -t242 * t265 + t269 * t284;
t1 = [0, -g(1) * t271 - g(2) * t267, t275, 0, 0, 0, 0, 0, -g(1) * (-t263 * t282 + t279) - g(2) * (t263 * t280 + t281) - g(3) * t285, -g(1) * (-t263 * t281 - t280) - g(2) * (t263 * t279 - t282) - g(3) * t261 * t270, t247, -g(1) * t277 - g(2) * t278 - g(3) * t276, t247, g(1) * t243 + g(2) * t241 + g(3) * t249, g(1) * t242 + g(2) * t240 - g(3) * t248, -g(1) * (t243 * pkin(3) - t242 * qJ(4) + t277) - g(2) * (t241 * pkin(3) - t240 * qJ(4) + t278) - g(3) * (t249 * pkin(3) + t248 * qJ(4) + t276) 0, 0, 0, 0, 0, -g(1) * t238 - g(2) * t239 - g(3) * t244, -g(1) * (-t242 * t269 - t265 * t284) - g(2) * (-t240 * t269 + t265 * t283) - g(3) * (t248 * t269 - t263 * t265) 0, 0, 0, 0, 0, -g(1) * (t238 * t268 + t243 * t264) - g(2) * (t239 * t268 + t241 * t264) - g(3) * (t244 * t268 + t249 * t264) -g(1) * (-t238 * t264 + t243 * t268) - g(2) * (-t239 * t264 + t241 * t268) - g(3) * (-t244 * t264 + t249 * t268);];
U_reg  = t1;
