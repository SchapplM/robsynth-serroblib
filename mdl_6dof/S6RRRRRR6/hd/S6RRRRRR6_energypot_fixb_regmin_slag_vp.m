% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:34
% EndTime: 2019-03-10 04:21:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (134->52), mult. (220->95), div. (0->0), fcn. (284->14), ass. (0->33)
t269 = sin(pkin(6));
t273 = sin(qJ(2));
t287 = t269 * t273;
t274 = sin(qJ(1));
t286 = t269 * t274;
t276 = cos(qJ(3));
t285 = t269 * t276;
t277 = cos(qJ(2));
t284 = t269 * t277;
t278 = cos(qJ(1));
t283 = t269 * t278;
t282 = t274 * t273;
t281 = t274 * t277;
t280 = t278 * t273;
t279 = t278 * t277;
t275 = cos(qJ(5));
t272 = sin(qJ(3));
t271 = sin(qJ(5));
t270 = cos(pkin(6));
t268 = qJ(3) + qJ(4);
t267 = qJ(5) + qJ(6);
t266 = cos(t268);
t265 = cos(t267);
t264 = sin(t268);
t263 = sin(t267);
t262 = -t270 * t282 + t279;
t261 = t270 * t281 + t280;
t260 = t270 * t280 + t281;
t259 = -t270 * t279 + t282;
t258 = t270 * t264 + t266 * t287;
t257 = t262 * t266 + t264 * t286;
t256 = t260 * t266 - t264 * t283;
t1 = [0, -g(1) * t278 - g(2) * t274, g(1) * t274 - g(2) * t278, 0, 0, 0, 0, 0, -g(1) * t262 - g(2) * t260 - g(3) * t287, g(1) * t261 + g(2) * t259 - g(3) * t284, 0, 0, 0, 0, 0, -g(1) * (t262 * t276 + t272 * t286) - g(2) * (t260 * t276 - t272 * t283) - g(3) * (t270 * t272 + t273 * t285) -g(1) * (-t262 * t272 + t274 * t285) - g(2) * (-t260 * t272 - t276 * t283) - g(3) * (t270 * t276 - t272 * t287) 0, 0, 0, 0, 0, -g(1) * t257 - g(2) * t256 - g(3) * t258, -g(1) * (-t262 * t264 + t266 * t286) - g(2) * (-t260 * t264 - t266 * t283) - g(3) * (-t264 * t287 + t270 * t266) 0, 0, 0, 0, 0, -g(1) * (t257 * t275 + t261 * t271) - g(2) * (t256 * t275 + t259 * t271) - g(3) * (t258 * t275 - t271 * t284) -g(1) * (-t257 * t271 + t261 * t275) - g(2) * (-t256 * t271 + t259 * t275) - g(3) * (-t258 * t271 - t275 * t284) 0, 0, 0, 0, 0, -g(1) * (t257 * t265 + t261 * t263) - g(2) * (t256 * t265 + t259 * t263) - g(3) * (t258 * t265 - t263 * t284) -g(1) * (-t257 * t263 + t261 * t265) - g(2) * (-t256 * t263 + t259 * t265) - g(3) * (-t258 * t263 - t265 * t284);];
U_reg  = t1;
