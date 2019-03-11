% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:28
% EndTime: 2019-03-09 21:38:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (245->65), mult. (587->94), div. (0->0), fcn. (748->10), ass. (0->43)
t281 = pkin(10) - qJ(6);
t258 = sin(pkin(6));
t262 = sin(qJ(2));
t280 = t258 * t262;
t263 = sin(qJ(1));
t279 = t258 * t263;
t265 = cos(qJ(3));
t278 = t258 * t265;
t266 = cos(qJ(2));
t277 = t258 * t266;
t267 = cos(qJ(1));
t276 = t258 * t267;
t275 = t263 * t262;
t274 = t263 * t266;
t273 = t267 * t262;
t272 = t267 * t266;
t259 = cos(pkin(6));
t247 = t259 * t273 + t274;
t261 = sin(qJ(3));
t236 = t247 * t265 - t261 * t276;
t246 = -t259 * t272 + t275;
t260 = sin(qJ(4));
t264 = cos(qJ(4));
t226 = t236 * t260 - t246 * t264;
t249 = -t259 * t275 + t272;
t238 = t249 * t265 + t261 * t279;
t248 = t259 * t274 + t273;
t228 = t238 * t260 - t248 * t264;
t245 = t259 * t261 + t262 * t278;
t233 = t245 * t260 + t264 * t277;
t271 = g(1) * t228 + g(2) * t226 + g(3) * t233;
t235 = t247 * t261 + t265 * t276;
t237 = t249 * t261 - t263 * t278;
t244 = -t259 * t265 + t261 * t280;
t223 = g(1) * t237 + g(2) * t235 + g(3) * t244;
t229 = t238 * t264 + t248 * t260;
t270 = t267 * pkin(1) + t249 * pkin(2) + t238 * pkin(3) + t229 * pkin(4) + pkin(8) * t279 + t248 * pkin(9) + t228 * qJ(5);
t234 = t245 * t264 - t260 * t277;
t269 = pkin(2) * t280 + t245 * pkin(3) + t234 * pkin(4) + t259 * pkin(8) - pkin(9) * t277 + t233 * qJ(5) + pkin(7);
t227 = t236 * t264 + t246 * t260;
t268 = t263 * pkin(1) + t247 * pkin(2) + t236 * pkin(3) + t227 * pkin(4) - pkin(8) * t276 + t246 * pkin(9) + t226 * qJ(5);
t222 = -g(1) * t229 - g(2) * t227 - g(3) * t234;
t1 = [0, -g(1) * t267 - g(2) * t263, g(1) * t263 - g(2) * t267, 0, 0, 0, 0, 0, -g(1) * t249 - g(2) * t247 - g(3) * t280, g(1) * t248 + g(2) * t246 - g(3) * t277, 0, 0, 0, 0, 0, -g(1) * t238 - g(2) * t236 - g(3) * t245, t223, 0, 0, 0, 0, 0, t222, t271, t222, -t223, -t271, -g(1) * (t237 * pkin(10) + t270) - g(2) * (t235 * pkin(10) + t268) - g(3) * (t244 * pkin(10) + t269) t222, -t271, t223, -g(1) * (t229 * pkin(5) + t281 * t237 + t270) - g(2) * (t227 * pkin(5) + t281 * t235 + t268) - g(3) * (t234 * pkin(5) + t281 * t244 + t269);];
U_reg  = t1;
