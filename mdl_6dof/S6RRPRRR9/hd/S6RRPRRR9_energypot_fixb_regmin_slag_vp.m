% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR9
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:13:58
% EndTime: 2019-03-09 14:13:58
% DurationCPUTime: 0.15s
% Computational Cost: add. (155->64), mult. (234->109), div. (0->0), fcn. (291->14), ass. (0->33)
t269 = sin(pkin(6));
t273 = sin(qJ(2));
t286 = t269 * t273;
t274 = sin(qJ(1));
t285 = t269 * t274;
t276 = cos(qJ(2));
t284 = t269 * t276;
t277 = cos(qJ(1));
t283 = t269 * t277;
t282 = t274 * t273;
t281 = t274 * t276;
t280 = t277 * t273;
t279 = t277 * t276;
t267 = pkin(12) + qJ(4);
t271 = cos(pkin(6));
t257 = -t271 * t279 + t282;
t259 = t271 * t281 + t280;
t278 = -g(1) * t259 - g(2) * t257 + g(3) * t284;
t275 = cos(qJ(6));
t272 = sin(qJ(6));
t270 = cos(pkin(12));
t268 = sin(pkin(12));
t266 = qJ(5) + t267;
t265 = cos(t267);
t264 = sin(t267);
t263 = cos(t266);
t262 = sin(t266);
t260 = -t271 * t282 + t279;
t258 = t271 * t280 + t281;
t256 = t271 * t262 + t263 * t286;
t255 = t260 * t263 + t262 * t285;
t254 = t258 * t263 - t262 * t283;
t1 = [0, -g(1) * t277 - g(2) * t274, g(1) * t274 - g(2) * t277, 0, 0, 0, 0, 0, -g(1) * t260 - g(2) * t258 - g(3) * t286, -t278, -g(1) * (t260 * t270 + t268 * t285) - g(2) * (t258 * t270 - t268 * t283) - g(3) * (t271 * t268 + t270 * t286) -g(1) * (-t260 * t268 + t270 * t285) - g(2) * (-t258 * t268 - t270 * t283) - g(3) * (-t268 * t286 + t271 * t270) t278, -g(1) * (t277 * pkin(1) + t260 * pkin(2) + pkin(8) * t285 + t259 * qJ(3)) - g(2) * (t274 * pkin(1) + t258 * pkin(2) - pkin(8) * t283 + t257 * qJ(3)) - g(3) * (t271 * pkin(8) + pkin(7) + (pkin(2) * t273 - qJ(3) * t276) * t269) 0, 0, 0, 0, 0, -g(1) * (t260 * t265 + t264 * t285) - g(2) * (t258 * t265 - t264 * t283) - g(3) * (t271 * t264 + t265 * t286) -g(1) * (-t260 * t264 + t265 * t285) - g(2) * (-t258 * t264 - t265 * t283) - g(3) * (-t264 * t286 + t271 * t265) 0, 0, 0, 0, 0, -g(1) * t255 - g(2) * t254 - g(3) * t256, -g(1) * (-t260 * t262 + t263 * t285) - g(2) * (-t258 * t262 - t263 * t283) - g(3) * (-t262 * t286 + t271 * t263) 0, 0, 0, 0, 0, -g(1) * (t255 * t275 + t259 * t272) - g(2) * (t254 * t275 + t257 * t272) - g(3) * (t256 * t275 - t272 * t284) -g(1) * (-t255 * t272 + t259 * t275) - g(2) * (-t254 * t272 + t257 * t275) - g(3) * (-t256 * t272 - t275 * t284);];
U_reg  = t1;
