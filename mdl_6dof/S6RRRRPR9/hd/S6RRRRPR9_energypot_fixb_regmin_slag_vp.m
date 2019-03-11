% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:55:43
% EndTime: 2019-03-09 22:55:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (199->74), mult. (315->119), div. (0->0), fcn. (395->14), ass. (0->41)
t280 = sin(qJ(1));
t283 = cos(qJ(1));
t299 = -g(1) * t280 + g(2) * t283;
t275 = sin(pkin(6));
t279 = sin(qJ(2));
t296 = t275 * t279;
t295 = t275 * t280;
t281 = cos(qJ(3));
t294 = t275 * t281;
t282 = cos(qJ(2));
t293 = t275 * t282;
t292 = t275 * t283;
t277 = cos(pkin(6));
t278 = sin(qJ(3));
t291 = t277 * t278;
t290 = t280 * t279;
t289 = t280 * t282;
t288 = t283 * t279;
t287 = t283 * t282;
t261 = t277 * t288 + t289;
t273 = qJ(3) + qJ(4);
t270 = sin(t273);
t271 = cos(t273);
t254 = t261 * t270 + t271 * t292;
t263 = -t277 * t290 + t287;
t256 = t263 * t270 - t271 * t295;
t258 = t270 * t296 - t277 * t271;
t285 = g(1) * t256 + g(2) * t254 + g(3) * t258;
t284 = -pkin(10) - pkin(9);
t276 = cos(pkin(12));
t274 = sin(pkin(12));
t272 = pkin(12) + qJ(6);
t269 = cos(t272);
t268 = sin(t272);
t267 = t281 * pkin(3) + pkin(2);
t262 = t277 * t289 + t288;
t260 = -t277 * t287 + t290;
t259 = t277 * t270 + t271 * t296;
t257 = t263 * t271 + t270 * t295;
t255 = t261 * t271 - t270 * t292;
t1 = [0, -g(1) * t283 - g(2) * t280, -t299, 0, 0, 0, 0, 0, -g(1) * t263 - g(2) * t261 - g(3) * t296, g(1) * t262 + g(2) * t260 - g(3) * t293, 0, 0, 0, 0, 0, -g(1) * (t263 * t281 + t278 * t295) - g(2) * (t261 * t281 - t278 * t292) - g(3) * (t279 * t294 + t291) -g(1) * (-t263 * t278 + t280 * t294) - g(2) * (-t261 * t278 - t281 * t292) - g(3) * (t277 * t281 - t278 * t296) 0, 0, 0, 0, 0, -g(1) * t257 - g(2) * t255 - g(3) * t259, t285, -g(1) * (t257 * t276 + t262 * t274) - g(2) * (t255 * t276 + t260 * t274) - g(3) * (t259 * t276 - t274 * t293) -g(1) * (-t257 * t274 + t262 * t276) - g(2) * (-t255 * t274 + t260 * t276) - g(3) * (-t259 * t274 - t276 * t293) -t285, -g(1) * (t283 * pkin(1) + t257 * pkin(4) + t256 * qJ(5) - t262 * t284 + t263 * t267) - g(2) * (t280 * pkin(1) + t255 * pkin(4) + t254 * qJ(5) - t260 * t284 + t261 * t267) - g(3) * (pkin(3) * t291 + t259 * pkin(4) + t277 * pkin(8) + t258 * qJ(5) + pkin(7)) + (-g(3) * (t267 * t279 + t282 * t284) + t299 * (pkin(3) * t278 + pkin(8))) * t275, 0, 0, 0, 0, 0, -g(1) * (t257 * t269 + t262 * t268) - g(2) * (t255 * t269 + t260 * t268) - g(3) * (t259 * t269 - t268 * t293) -g(1) * (-t257 * t268 + t262 * t269) - g(2) * (-t255 * t268 + t260 * t269) - g(3) * (-t259 * t268 - t269 * t293);];
U_reg  = t1;
