% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP7
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
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:11:06
% EndTime: 2019-03-09 17:11:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (227->71), mult. (392->107), div. (0->0), fcn. (485->12), ass. (0->49)
t282 = sin(pkin(6));
t287 = sin(qJ(2));
t308 = t282 * t287;
t288 = sin(qJ(1));
t307 = t282 * t288;
t290 = cos(qJ(3));
t306 = t282 * t290;
t291 = cos(qJ(2));
t305 = t282 * t291;
t292 = cos(qJ(1));
t304 = t282 * t292;
t283 = cos(pkin(6));
t286 = sin(qJ(3));
t303 = t283 * t286;
t302 = t288 * t287;
t301 = t288 * t291;
t300 = t292 * t287;
t299 = t292 * t291;
t298 = t286 * t307;
t275 = t290 * pkin(3) + pkin(2);
t284 = -qJ(4) - pkin(9);
t297 = pkin(3) * t303 + t283 * pkin(8) + t275 * t308 + t284 * t305 + pkin(7);
t267 = t283 * t301 + t300;
t268 = -t283 * t302 + t299;
t296 = t292 * pkin(1) + pkin(3) * t298 + pkin(8) * t307 - t267 * t284 + t268 * t275;
t266 = t283 * t300 + t301;
t281 = qJ(3) + pkin(11);
t276 = sin(t281);
t277 = cos(t281);
t254 = t266 * t277 - t276 * t304;
t265 = -t283 * t299 + t302;
t285 = sin(qJ(5));
t289 = cos(qJ(5));
t247 = t254 * t285 - t265 * t289;
t256 = t268 * t277 + t276 * t307;
t249 = t256 * t285 - t267 * t289;
t260 = t283 * t276 + t277 * t308;
t251 = t260 * t285 + t289 * t305;
t295 = g(1) * t249 + g(2) * t247 + g(3) * t251;
t294 = -g(1) * t267 - g(2) * t265 + g(3) * t305;
t293 = t266 * t275 - t265 * t284 + t288 * pkin(1) + (-pkin(3) * t286 - pkin(8)) * t304;
t259 = t276 * t308 - t283 * t277;
t255 = t268 * t276 - t277 * t307;
t253 = t266 * t276 + t277 * t304;
t252 = t260 * t289 - t285 * t305;
t250 = t256 * t289 + t267 * t285;
t248 = t254 * t289 + t265 * t285;
t246 = -g(1) * t250 - g(2) * t248 - g(3) * t252;
t1 = [0, -g(1) * t292 - g(2) * t288, g(1) * t288 - g(2) * t292, 0, 0, 0, 0, 0, -g(1) * t268 - g(2) * t266 - g(3) * t308, -t294, 0, 0, 0, 0, 0, -g(1) * (t268 * t290 + t298) - g(2) * (t266 * t290 - t286 * t304) - g(3) * (t287 * t306 + t303) -g(1) * (-t268 * t286 + t288 * t306) - g(2) * (-t266 * t286 - t290 * t304) - g(3) * (t283 * t290 - t286 * t308) t294, -g(1) * t296 - g(2) * t293 - g(3) * t297, 0, 0, 0, 0, 0, t246, t295, t246, -g(1) * t255 - g(2) * t253 - g(3) * t259, -t295, -g(1) * (t256 * pkin(4) + t250 * pkin(5) + t255 * pkin(10) + t249 * qJ(6) + t296) - g(2) * (t254 * pkin(4) + t248 * pkin(5) + t253 * pkin(10) + t247 * qJ(6) + t293) - g(3) * (t260 * pkin(4) + t252 * pkin(5) + t259 * pkin(10) + t251 * qJ(6) + t297);];
U_reg  = t1;
