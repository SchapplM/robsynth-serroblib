% Calculate minimal parameter regressor of potential energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:11
% EndTime: 2019-12-31 22:55:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (140->47), mult. (386->92), div. (0->0), fcn. (508->14), ass. (0->42)
t274 = sin(pkin(6));
t277 = cos(pkin(5));
t298 = t274 * t277;
t275 = sin(pkin(5));
t282 = sin(qJ(1));
t297 = t275 * t282;
t286 = cos(qJ(2));
t296 = t275 * t286;
t287 = cos(qJ(1));
t295 = t275 * t287;
t276 = cos(pkin(6));
t294 = t276 * t286;
t281 = sin(qJ(2));
t293 = t282 * t281;
t292 = t282 * t286;
t291 = t287 * t281;
t290 = t287 * t286;
t270 = t277 * t290 - t293;
t289 = -t270 * t276 + t274 * t295;
t272 = -t277 * t292 - t291;
t288 = t272 * t276 + t274 * t297;
t285 = cos(qJ(3));
t284 = cos(qJ(4));
t283 = cos(qJ(5));
t280 = sin(qJ(3));
t279 = sin(qJ(4));
t278 = sin(qJ(5));
t273 = -t277 * t293 + t290;
t271 = t277 * t291 + t292;
t269 = -t274 * t296 + t277 * t276;
t268 = -t272 * t274 + t276 * t297;
t267 = -t270 * t274 - t276 * t295;
t266 = t280 * t298 + (t280 * t294 + t281 * t285) * t275;
t265 = -t285 * t298 + (t280 * t281 - t285 * t294) * t275;
t264 = t273 * t285 + t288 * t280;
t263 = t273 * t280 - t288 * t285;
t262 = t271 * t285 - t289 * t280;
t261 = t271 * t280 + t289 * t285;
t260 = t266 * t284 + t269 * t279;
t259 = t264 * t284 + t268 * t279;
t258 = t262 * t284 + t267 * t279;
t1 = [0, -g(1) * t287 - g(2) * t282, g(1) * t282 - g(2) * t287, 0, 0, 0, 0, 0, -g(3) * t275 * t281 - g(1) * t273 - g(2) * t271, -g(1) * t272 - g(2) * t270 - g(3) * t296, 0, 0, 0, 0, 0, -g(1) * t264 - g(2) * t262 - g(3) * t266, g(1) * t263 + g(2) * t261 + g(3) * t265, 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t258 - g(3) * t260, -g(1) * (-t264 * t279 + t268 * t284) - g(2) * (-t262 * t279 + t267 * t284) - g(3) * (-t266 * t279 + t269 * t284), 0, 0, 0, 0, 0, -g(1) * (t259 * t283 + t263 * t278) - g(2) * (t258 * t283 + t261 * t278) - g(3) * (t260 * t283 + t265 * t278), -g(1) * (-t259 * t278 + t263 * t283) - g(2) * (-t258 * t278 + t261 * t283) - g(3) * (-t260 * t278 + t265 * t283);];
U_reg = t1;
