% Calculate minimal parameter regressor of potential energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:42
% EndTime: 2019-12-31 19:19:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (147->53), mult. (401->102), div. (0->0), fcn. (520->14), ass. (0->43)
t273 = sin(pkin(6));
t277 = cos(pkin(5));
t297 = t273 * t277;
t274 = sin(pkin(5));
t275 = cos(pkin(11));
t296 = t274 * t275;
t281 = sin(qJ(1));
t295 = t274 * t281;
t285 = cos(qJ(1));
t294 = t274 * t285;
t276 = cos(pkin(6));
t293 = t275 * t276;
t272 = sin(pkin(11));
t292 = t281 * t272;
t291 = t281 * t275;
t290 = t285 * t272;
t289 = t285 * t275;
t288 = g(1) * t281 - g(2) * t285;
t268 = t277 * t289 - t292;
t287 = -t268 * t276 + t273 * t294;
t270 = -t277 * t291 - t290;
t286 = t270 * t276 + t273 * t295;
t284 = cos(qJ(3));
t283 = cos(qJ(4));
t282 = cos(qJ(5));
t280 = sin(qJ(3));
t279 = sin(qJ(4));
t278 = sin(qJ(5));
t271 = -t277 * t292 + t289;
t269 = t277 * t290 + t291;
t267 = -t273 * t296 + t277 * t276;
t266 = -t270 * t273 + t276 * t295;
t265 = -t268 * t273 - t276 * t294;
t264 = t280 * t297 + (t272 * t284 + t280 * t293) * t274;
t263 = -t284 * t297 + (t272 * t280 - t284 * t293) * t274;
t262 = t271 * t284 + t286 * t280;
t261 = t271 * t280 - t286 * t284;
t260 = t269 * t284 - t287 * t280;
t259 = t269 * t280 + t287 * t284;
t258 = t264 * t283 + t267 * t279;
t257 = t262 * t283 + t266 * t279;
t256 = t260 * t283 + t265 * t279;
t1 = [0, -g(1) * t285 - g(2) * t281, t288, -g(3) * t274 * t272 - g(1) * t271 - g(2) * t269, -g(1) * t270 - g(2) * t268 - g(3) * t296, -g(3) * t277 - t288 * t274, -g(1) * (t285 * pkin(1) + qJ(2) * t295) - g(2) * (t281 * pkin(1) - qJ(2) * t294) - g(3) * (t277 * qJ(2) + pkin(7)), 0, 0, 0, 0, 0, -g(1) * t262 - g(2) * t260 - g(3) * t264, g(1) * t261 + g(2) * t259 + g(3) * t263, 0, 0, 0, 0, 0, -g(1) * t257 - g(2) * t256 - g(3) * t258, -g(1) * (-t262 * t279 + t266 * t283) - g(2) * (-t260 * t279 + t265 * t283) - g(3) * (-t264 * t279 + t267 * t283), 0, 0, 0, 0, 0, -g(1) * (t257 * t282 + t261 * t278) - g(2) * (t256 * t282 + t259 * t278) - g(3) * (t258 * t282 + t263 * t278), -g(1) * (-t257 * t278 + t261 * t282) - g(2) * (-t256 * t278 + t259 * t282) - g(3) * (-t258 * t278 + t263 * t282);];
U_reg = t1;
