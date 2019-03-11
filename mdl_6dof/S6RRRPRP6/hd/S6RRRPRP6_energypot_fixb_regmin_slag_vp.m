% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:08
% EndTime: 2019-03-09 17:01:08
% DurationCPUTime: 0.14s
% Computational Cost: add. (164->69), mult. (284->105), div. (0->0), fcn. (339->12), ass. (0->46)
t268 = cos(pkin(6));
t277 = cos(qJ(2));
t278 = cos(qJ(1));
t285 = t277 * t278;
t273 = sin(qJ(2));
t274 = sin(qJ(1));
t288 = t273 * t274;
t249 = -t268 * t285 + t288;
t271 = sin(qJ(5));
t296 = t249 * t271;
t286 = t274 * t277;
t287 = t273 * t278;
t251 = t268 * t286 + t287;
t295 = t251 * t271;
t267 = sin(pkin(6));
t294 = t267 * t273;
t293 = t267 * t274;
t276 = cos(qJ(3));
t292 = t267 * t276;
t291 = t267 * t277;
t290 = t267 * t278;
t272 = sin(qJ(3));
t289 = t268 * t272;
t284 = t271 * t291;
t283 = t272 * t293;
t260 = pkin(3) * t276 + pkin(2);
t270 = -qJ(4) - pkin(9);
t282 = pkin(3) * t289 + t268 * pkin(8) + t260 * t294 + t270 * t291 + pkin(7);
t252 = -t268 * t288 + t285;
t281 = t278 * pkin(1) + pkin(3) * t283 + pkin(8) * t293 - t251 * t270 + t252 * t260;
t280 = -g(1) * t251 - g(2) * t249 + g(3) * t291;
t250 = t268 * t287 + t286;
t279 = t250 * t260 - t249 * t270 + t274 * pkin(1) + (-pkin(3) * t272 - pkin(8)) * t290;
t275 = cos(qJ(5));
t269 = -qJ(6) - pkin(10);
t266 = qJ(3) + pkin(11);
t262 = cos(t266);
t261 = sin(t266);
t259 = pkin(5) * t275 + pkin(4);
t246 = t261 * t268 + t262 * t294;
t245 = t261 * t294 - t262 * t268;
t242 = t252 * t262 + t261 * t293;
t241 = t252 * t261 - t262 * t293;
t240 = t250 * t262 - t261 * t290;
t239 = t250 * t261 + t262 * t290;
t1 = [0, -g(1) * t278 - g(2) * t274, g(1) * t274 - g(2) * t278, 0, 0, 0, 0, 0, -g(1) * t252 - g(2) * t250 - g(3) * t294, -t280, 0, 0, 0, 0, 0, -g(1) * (t252 * t276 + t283) - g(2) * (t250 * t276 - t272 * t290) - g(3) * (t273 * t292 + t289) -g(1) * (-t252 * t272 + t274 * t292) - g(2) * (-t250 * t272 - t276 * t290) - g(3) * (t268 * t276 - t272 * t294) t280, -g(1) * t281 - g(2) * t279 - g(3) * t282, 0, 0, 0, 0, 0, -g(1) * (t242 * t275 + t295) - g(2) * (t240 * t275 + t296) - g(3) * (t246 * t275 - t284) -g(1) * (-t242 * t271 + t251 * t275) - g(2) * (-t240 * t271 + t249 * t275) - g(3) * (-t246 * t271 - t275 * t291) -g(1) * t241 - g(2) * t239 - g(3) * t245, -g(1) * (pkin(5) * t295 - t241 * t269 + t242 * t259 + t281) - g(2) * (pkin(5) * t296 - t239 * t269 + t240 * t259 + t279) - g(3) * (-pkin(5) * t284 - t245 * t269 + t246 * t259 + t282);];
U_reg  = t1;
