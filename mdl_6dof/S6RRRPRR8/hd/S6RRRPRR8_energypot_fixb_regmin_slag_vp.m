% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:21
% EndTime: 2019-03-09 18:52:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (136->61), mult. (235->101), div. (0->0), fcn. (292->14), ass. (0->38)
t278 = sin(qJ(1));
t282 = cos(qJ(1));
t297 = -g(1) * t278 + g(2) * t282;
t272 = sin(pkin(6));
t277 = sin(qJ(2));
t294 = t272 * t277;
t293 = t272 * t278;
t280 = cos(qJ(3));
t292 = t272 * t280;
t281 = cos(qJ(2));
t291 = t272 * t281;
t290 = t272 * t282;
t273 = cos(pkin(6));
t276 = sin(qJ(3));
t289 = t273 * t276;
t288 = t278 * t277;
t287 = t278 * t281;
t286 = t282 * t277;
t285 = t282 * t281;
t260 = -t273 * t285 + t288;
t262 = t273 * t287 + t286;
t283 = -g(1) * t262 - g(2) * t260 + g(3) * t291;
t279 = cos(qJ(5));
t275 = sin(qJ(5));
t274 = -qJ(4) - pkin(9);
t271 = qJ(5) + qJ(6);
t270 = qJ(3) + pkin(12);
t269 = cos(t271);
t268 = sin(t271);
t267 = cos(t270);
t266 = sin(t270);
t265 = t280 * pkin(3) + pkin(2);
t263 = -t273 * t288 + t285;
t261 = t273 * t286 + t287;
t259 = t273 * t266 + t267 * t294;
t258 = t263 * t267 + t266 * t293;
t257 = t261 * t267 - t266 * t290;
t1 = [0, -g(1) * t282 - g(2) * t278, -t297, 0, 0, 0, 0, 0, -g(1) * t263 - g(2) * t261 - g(3) * t294, -t283, 0, 0, 0, 0, 0, -g(1) * (t263 * t280 + t276 * t293) - g(2) * (t261 * t280 - t276 * t290) - g(3) * (t277 * t292 + t289) -g(1) * (-t263 * t276 + t278 * t292) - g(2) * (-t261 * t276 - t280 * t290) - g(3) * (t273 * t280 - t276 * t294) t283, -g(1) * (t282 * pkin(1) - t262 * t274 + t263 * t265) - g(2) * (t278 * pkin(1) - t260 * t274 + t261 * t265) - g(3) * (pkin(3) * t289 + t273 * pkin(8) + pkin(7)) + (-g(3) * (t265 * t277 + t274 * t281) + t297 * (pkin(3) * t276 + pkin(8))) * t272, 0, 0, 0, 0, 0, -g(1) * (t258 * t279 + t262 * t275) - g(2) * (t257 * t279 + t260 * t275) - g(3) * (t259 * t279 - t275 * t291) -g(1) * (-t258 * t275 + t262 * t279) - g(2) * (-t257 * t275 + t260 * t279) - g(3) * (-t259 * t275 - t279 * t291) 0, 0, 0, 0, 0, -g(1) * (t258 * t269 + t262 * t268) - g(2) * (t257 * t269 + t260 * t268) - g(3) * (t259 * t269 - t268 * t291) -g(1) * (-t258 * t268 + t262 * t269) - g(2) * (-t257 * t268 + t260 * t269) - g(3) * (-t259 * t268 - t269 * t291);];
U_reg  = t1;
