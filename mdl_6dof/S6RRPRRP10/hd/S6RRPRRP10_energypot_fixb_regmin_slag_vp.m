% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:42:31
% EndTime: 2019-03-09 12:42:32
% DurationCPUTime: 0.20s
% Computational Cost: add. (244->81), mult. (417->119), div. (0->0), fcn. (520->12), ass. (0->51)
t271 = sin(pkin(11));
t298 = pkin(3) * t271;
t278 = sin(qJ(1));
t297 = g(1) * t278;
t281 = cos(qJ(1));
t296 = g(2) * t281;
t274 = cos(pkin(6));
t295 = t274 * pkin(8) + pkin(7);
t272 = sin(pkin(6));
t277 = sin(qJ(2));
t294 = t272 * t277;
t293 = t272 * t278;
t280 = cos(qJ(2));
t292 = t272 * t280;
t291 = t272 * t281;
t290 = t274 * t271;
t289 = t277 * t278;
t288 = t277 * t281;
t287 = t278 * t280;
t286 = t280 * t281;
t285 = t281 * pkin(1) + pkin(8) * t293;
t257 = t274 * t288 + t287;
t270 = pkin(11) + qJ(4);
t265 = sin(t270);
t266 = cos(t270);
t249 = t257 * t266 - t265 * t291;
t256 = -t274 * t286 + t289;
t276 = sin(qJ(5));
t279 = cos(qJ(5));
t242 = t249 * t276 - t256 * t279;
t259 = -t274 * t289 + t286;
t251 = t259 * t266 + t265 * t293;
t258 = t274 * t287 + t288;
t244 = t251 * t276 - t258 * t279;
t253 = t265 * t274 + t266 * t294;
t246 = t253 * t276 + t279 * t292;
t284 = g(1) * t244 + g(2) * t242 + g(3) * t246;
t248 = t257 * t265 + t266 * t291;
t250 = t259 * t265 - t266 * t293;
t252 = t265 * t294 - t274 * t266;
t283 = g(1) * t250 + g(2) * t248 + g(3) * t252;
t282 = -g(1) * t258 - g(2) * t256 + g(3) * t292;
t275 = -pkin(9) - qJ(3);
t273 = cos(pkin(11));
t268 = t278 * pkin(1);
t264 = pkin(3) * t273 + pkin(2);
t247 = t253 * t279 - t276 * t292;
t245 = t251 * t279 + t258 * t276;
t243 = t249 * t279 + t256 * t276;
t241 = -g(1) * t245 - g(2) * t243 - g(3) * t247;
t1 = [0, -g(1) * t281 - g(2) * t278, -t296 + t297, 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t257 - g(3) * t294, -t282, -g(1) * (t259 * t273 + t271 * t293) - g(2) * (t257 * t273 - t271 * t291) - g(3) * (t273 * t294 + t290) -g(1) * (-t259 * t271 + t273 * t293) - g(2) * (-t257 * t271 - t273 * t291) - g(3) * (-t271 * t294 + t273 * t274) t282, -g(1) * (pkin(2) * t259 + qJ(3) * t258 + t285) - g(2) * (t257 * pkin(2) - pkin(8) * t291 + t256 * qJ(3) + t268) - g(3) * ((pkin(2) * t277 - qJ(3) * t280) * t272 + t295) 0, 0, 0, 0, 0, -g(1) * t251 - g(2) * t249 - g(3) * t253, t283, 0, 0, 0, 0, 0, t241, t284, t241, -t283, -t284, -g(1) * (pkin(4) * t251 + pkin(5) * t245 + pkin(10) * t250 + qJ(6) * t244 - t258 * t275 + t259 * t264 + t285) - g(2) * (t249 * pkin(4) + t243 * pkin(5) + t248 * pkin(10) + t242 * qJ(6) - t256 * t275 + t257 * t264 + t268) - g(3) * (pkin(3) * t290 + pkin(4) * t253 + pkin(5) * t247 + pkin(10) * t252 + qJ(6) * t246 + t295) + (-t297 * t298 - g(3) * (t264 * t277 + t275 * t280) - (-pkin(8) - t298) * t296) * t272;];
U_reg  = t1;
