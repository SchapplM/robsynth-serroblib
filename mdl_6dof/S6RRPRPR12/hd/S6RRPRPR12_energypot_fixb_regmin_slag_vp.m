% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:21:46
% EndTime: 2019-03-09 11:21:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (123->64), mult. (248->97), div. (0->0), fcn. (294->12), ass. (0->42)
t263 = sin(qJ(1));
t286 = g(1) * t263;
t267 = cos(qJ(1));
t285 = g(2) * t267;
t258 = cos(pkin(6));
t266 = cos(qJ(2));
t274 = t266 * t267;
t262 = sin(qJ(2));
t276 = t263 * t262;
t242 = -t258 * t274 + t276;
t261 = sin(qJ(4));
t284 = t242 * t261;
t275 = t263 * t266;
t277 = t262 * t267;
t244 = t258 * t275 + t277;
t283 = t244 * t261;
t257 = sin(pkin(6));
t282 = t257 * t262;
t281 = t257 * t263;
t280 = t257 * t266;
t279 = t257 * t267;
t278 = t261 * t266;
t273 = pkin(2) * t282 + t258 * pkin(8) + pkin(7);
t272 = -t285 + t286;
t243 = t258 * t277 + t275;
t271 = t263 * pkin(1) + t243 * pkin(2) + t242 * qJ(3);
t245 = -t258 * t276 + t274;
t270 = t267 * pkin(1) + t245 * pkin(2) + pkin(8) * t281 + qJ(3) * t244;
t269 = -g(1) * t244 - g(2) * t242 + g(3) * t280;
t268 = g(1) * t245 + g(2) * t243 + g(3) * t282;
t265 = cos(qJ(4));
t264 = cos(qJ(6));
t260 = sin(qJ(6));
t259 = -qJ(5) - pkin(9);
t256 = qJ(4) + pkin(11);
t252 = cos(t256);
t251 = sin(t256);
t250 = pkin(4) * t265 + pkin(3);
t239 = -t251 * t280 + t252 * t258;
t238 = t242 * t251 - t252 * t279;
t237 = t244 * t251 + t252 * t281;
t1 = [0, -g(1) * t267 - g(2) * t263, t272, 0, 0, 0, 0, 0, -t268, -t269, -g(3) * t258 - t272 * t257, t268, t269, -g(1) * t270 - g(2) * (-pkin(8) * t279 + t271) - g(3) * (-qJ(3) * t280 + t273) 0, 0, 0, 0, 0, -g(1) * (t265 * t281 + t283) - g(2) * (-t265 * t279 + t284) - g(3) * (-t257 * t278 + t258 * t265) -g(1) * (t244 * t265 - t261 * t281) - g(2) * (t242 * t265 + t261 * t279) - g(3) * (-t258 * t261 - t265 * t280) -t268, -g(1) * (pkin(4) * t283 - t245 * t259 + t270) - g(2) * (pkin(4) * t284 - t243 * t259 + t271) - g(3) * (t250 * t258 + t273) + (-t250 * t286 - g(3) * (-pkin(4) * t278 - qJ(3) * t266 - t259 * t262) - (-pkin(8) - t250) * t285) * t257, 0, 0, 0, 0, 0, -g(1) * (t237 * t264 + t245 * t260) - g(2) * (t238 * t264 + t243 * t260) - g(3) * (t239 * t264 + t260 * t282) -g(1) * (-t237 * t260 + t245 * t264) - g(2) * (-t238 * t260 + t243 * t264) - g(3) * (-t239 * t260 + t264 * t282);];
U_reg  = t1;
