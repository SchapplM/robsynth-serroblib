% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:25
% EndTime: 2019-03-10 02:27:25
% DurationCPUTime: 0.13s
% Computational Cost: add. (207->70), mult. (408->108), div. (0->0), fcn. (520->12), ass. (0->46)
t271 = cos(pkin(6));
t278 = cos(qJ(2));
t279 = cos(qJ(1));
t284 = t279 * t278;
t274 = sin(qJ(2));
t275 = sin(qJ(1));
t287 = t275 * t274;
t259 = -t271 * t284 + t287;
t272 = sin(qJ(4));
t293 = t259 * t272;
t285 = t279 * t274;
t286 = t275 * t278;
t261 = t271 * t286 + t285;
t292 = t261 * t272;
t270 = sin(pkin(6));
t291 = t270 * t274;
t277 = cos(qJ(3));
t290 = t270 * t277;
t289 = t270 * t278;
t288 = t270 * t279;
t283 = g(1) * t275 - g(2) * t279;
t260 = t271 * t285 + t286;
t273 = sin(qJ(3));
t252 = t260 * t277 - t273 * t288;
t269 = qJ(4) + qJ(5);
t267 = sin(t269);
t268 = cos(t269);
t245 = t252 * t267 - t259 * t268;
t262 = -t271 * t287 + t284;
t254 = t275 * t270 * t273 + t262 * t277;
t247 = t254 * t267 - t261 * t268;
t258 = t271 * t273 + t274 * t290;
t249 = t258 * t267 + t268 * t289;
t282 = g(1) * t247 + g(2) * t245 + g(3) * t249;
t251 = t260 * t273 + t277 * t288;
t253 = t262 * t273 - t275 * t290;
t257 = -t271 * t277 + t273 * t291;
t281 = g(1) * t253 + g(2) * t251 + g(3) * t257;
t280 = -pkin(11) - pkin(10);
t276 = cos(qJ(4));
t266 = t276 * pkin(4) + pkin(3);
t250 = t258 * t268 - t267 * t289;
t248 = t254 * t268 + t261 * t267;
t246 = t252 * t268 + t259 * t267;
t244 = -g(1) * t248 - g(2) * t246 - g(3) * t250;
t1 = [0, -g(1) * t279 - g(2) * t275, t283, 0, 0, 0, 0, 0, -g(1) * t262 - g(2) * t260 - g(3) * t291, g(1) * t261 + g(2) * t259 - g(3) * t289, 0, 0, 0, 0, 0, -g(1) * t254 - g(2) * t252 - g(3) * t258, t281, 0, 0, 0, 0, 0, -g(1) * (t254 * t276 + t292) - g(2) * (t252 * t276 + t293) - g(3) * (t258 * t276 - t272 * t289) -g(1) * (-t254 * t272 + t261 * t276) - g(2) * (-t252 * t272 + t259 * t276) - g(3) * (-t258 * t272 - t276 * t289) 0, 0, 0, 0, 0, t244, t282, t244, -t281, -t282, -g(1) * (t279 * pkin(1) + t262 * pkin(2) + pkin(4) * t292 + t248 * pkin(5) + t261 * pkin(9) + t247 * qJ(6) - t253 * t280 + t254 * t266) - g(2) * (t275 * pkin(1) + t260 * pkin(2) + pkin(4) * t293 + t246 * pkin(5) + t259 * pkin(9) + t245 * qJ(6) - t251 * t280 + t252 * t266) - g(3) * (t250 * pkin(5) + t271 * pkin(8) + t249 * qJ(6) - t257 * t280 + t258 * t266 + pkin(7)) + (-g(3) * (pkin(2) * t274 + (-pkin(4) * t272 - pkin(9)) * t278) - t283 * pkin(8)) * t270;];
U_reg  = t1;
