% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:14
% EndTime: 2019-03-09 19:34:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (143->59), mult. (354->98), div. (0->0), fcn. (453->12), ass. (0->38)
t260 = sin(pkin(6));
t265 = sin(qJ(2));
t282 = t260 * t265;
t266 = sin(qJ(1));
t281 = t260 * t266;
t269 = cos(qJ(3));
t280 = t260 * t269;
t270 = cos(qJ(2));
t279 = t260 * t270;
t271 = cos(qJ(1));
t278 = t260 * t271;
t277 = t266 * t265;
t276 = t266 * t270;
t275 = t271 * t265;
t274 = t271 * t270;
t261 = cos(pkin(6));
t254 = t261 * t275 + t276;
t264 = sin(qJ(3));
t247 = t254 * t264 + t269 * t278;
t256 = -t261 * t277 + t274;
t249 = t256 * t264 - t266 * t280;
t251 = -t261 * t269 + t264 * t282;
t273 = g(1) * t249 + g(2) * t247 + g(3) * t251;
t253 = -t261 * t274 + t277;
t255 = t261 * t276 + t275;
t272 = -g(1) * t255 - g(2) * t253 + g(3) * t279;
t268 = cos(qJ(5));
t267 = cos(qJ(6));
t263 = sin(qJ(5));
t262 = sin(qJ(6));
t252 = t261 * t264 + t265 * t280;
t250 = t256 * t269 + t264 * t281;
t248 = t254 * t269 - t264 * t278;
t246 = t251 * t263 + t252 * t268;
t245 = t249 * t263 + t250 * t268;
t244 = t247 * t263 + t248 * t268;
t243 = -g(1) * t250 - g(2) * t248 - g(3) * t252;
t1 = [0, -g(1) * t271 - g(2) * t266, g(1) * t266 - g(2) * t271, 0, 0, 0, 0, 0, -g(1) * t256 - g(2) * t254 - g(3) * t282, -t272, 0, 0, 0, 0, 0, t243, t273, t243, t272, -t273, -g(1) * (t271 * pkin(1) + t256 * pkin(2) + t250 * pkin(3) + pkin(8) * t281 + t255 * pkin(9) + t249 * qJ(4)) - g(2) * (t266 * pkin(1) + t254 * pkin(2) + t248 * pkin(3) - pkin(8) * t278 + t253 * pkin(9) + t247 * qJ(4)) - g(3) * (t252 * pkin(3) + t261 * pkin(8) + t251 * qJ(4) + pkin(7) + (pkin(2) * t265 - pkin(9) * t270) * t260) 0, 0, 0, 0, 0, -g(1) * t245 - g(2) * t244 - g(3) * t246, -g(1) * (t249 * t268 - t250 * t263) - g(2) * (t247 * t268 - t248 * t263) - g(3) * (t251 * t268 - t252 * t263) 0, 0, 0, 0, 0, -g(1) * (t245 * t267 - t255 * t262) - g(2) * (t244 * t267 - t253 * t262) - g(3) * (t246 * t267 + t262 * t279) -g(1) * (-t245 * t262 - t255 * t267) - g(2) * (-t244 * t262 - t253 * t267) - g(3) * (-t246 * t262 + t267 * t279);];
U_reg  = t1;
