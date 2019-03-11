% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:25:22
% EndTime: 2019-03-09 14:25:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (153->64), mult. (260->109), div. (0->0), fcn. (327->14), ass. (0->33)
t267 = sin(pkin(6));
t271 = sin(qJ(2));
t284 = t267 * t271;
t272 = sin(qJ(1));
t283 = t267 * t272;
t274 = cos(qJ(2));
t282 = t267 * t274;
t275 = cos(qJ(1));
t281 = t267 * t275;
t280 = t272 * t271;
t279 = t272 * t274;
t278 = t275 * t271;
t277 = t275 * t274;
t269 = cos(pkin(6));
t255 = -t269 * t277 + t280;
t257 = t269 * t279 + t278;
t276 = -g(1) * t257 - g(2) * t255 + g(3) * t282;
t273 = cos(qJ(5));
t270 = sin(qJ(5));
t268 = cos(pkin(12));
t266 = sin(pkin(12));
t265 = qJ(5) + qJ(6);
t264 = pkin(12) + qJ(4);
t263 = cos(t265);
t262 = sin(t265);
t261 = cos(t264);
t260 = sin(t264);
t258 = -t269 * t280 + t277;
t256 = t269 * t278 + t279;
t254 = t269 * t260 + t261 * t284;
t253 = t258 * t261 + t260 * t283;
t252 = t256 * t261 - t260 * t281;
t1 = [0, -g(1) * t275 - g(2) * t272, g(1) * t272 - g(2) * t275, 0, 0, 0, 0, 0, -g(1) * t258 - g(2) * t256 - g(3) * t284, -t276, -g(1) * (t258 * t268 + t266 * t283) - g(2) * (t256 * t268 - t266 * t281) - g(3) * (t269 * t266 + t268 * t284) -g(1) * (-t258 * t266 + t268 * t283) - g(2) * (-t256 * t266 - t268 * t281) - g(3) * (-t266 * t284 + t269 * t268) t276, -g(1) * (t275 * pkin(1) + t258 * pkin(2) + pkin(8) * t283 + t257 * qJ(3)) - g(2) * (t272 * pkin(1) + t256 * pkin(2) - pkin(8) * t281 + t255 * qJ(3)) - g(3) * (t269 * pkin(8) + pkin(7) + (pkin(2) * t271 - qJ(3) * t274) * t267) 0, 0, 0, 0, 0, -g(1) * t253 - g(2) * t252 - g(3) * t254, -g(1) * (-t258 * t260 + t261 * t283) - g(2) * (-t256 * t260 - t261 * t281) - g(3) * (-t260 * t284 + t269 * t261) 0, 0, 0, 0, 0, -g(1) * (t253 * t273 + t257 * t270) - g(2) * (t252 * t273 + t255 * t270) - g(3) * (t254 * t273 - t270 * t282) -g(1) * (-t253 * t270 + t257 * t273) - g(2) * (-t252 * t270 + t255 * t273) - g(3) * (-t254 * t270 - t273 * t282) 0, 0, 0, 0, 0, -g(1) * (t253 * t263 + t257 * t262) - g(2) * (t252 * t263 + t255 * t262) - g(3) * (t254 * t263 - t262 * t282) -g(1) * (-t253 * t262 + t257 * t263) - g(2) * (-t252 * t262 + t255 * t263) - g(3) * (-t254 * t262 - t263 * t282);];
U_reg  = t1;
