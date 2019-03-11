% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:54
% EndTime: 2019-03-09 15:53:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (170->66), mult. (296->104), div. (0->0), fcn. (356->12), ass. (0->41)
t265 = sin(pkin(6));
t270 = sin(qJ(2));
t289 = t265 * t270;
t271 = sin(qJ(1));
t288 = t265 * t271;
t273 = cos(qJ(3));
t287 = t265 * t273;
t274 = cos(qJ(2));
t286 = t265 * t274;
t275 = cos(qJ(1));
t285 = t265 * t275;
t266 = cos(pkin(6));
t269 = sin(qJ(3));
t284 = t266 * t269;
t283 = t271 * t270;
t282 = t271 * t274;
t281 = t275 * t270;
t280 = t275 * t274;
t279 = t269 * t288;
t258 = t273 * pkin(3) + pkin(2);
t267 = -qJ(4) - pkin(9);
t278 = pkin(3) * t284 + t266 * pkin(8) + t258 * t289 + t267 * t286 + pkin(7);
t250 = t266 * t282 + t281;
t251 = -t266 * t283 + t280;
t277 = t275 * pkin(1) + pkin(3) * t279 + pkin(8) * t288 - t250 * t267 + t251 * t258;
t248 = -t266 * t280 + t283;
t237 = -g(1) * t250 - g(2) * t248 + g(3) * t286;
t249 = t266 * t281 + t282;
t276 = t249 * t258 - t248 * t267 + t271 * pkin(1) + (-pkin(3) * t269 - pkin(8)) * t285;
t272 = cos(qJ(6));
t268 = sin(qJ(6));
t264 = qJ(3) + pkin(11);
t260 = cos(t264);
t259 = sin(t264);
t245 = t266 * t259 + t260 * t289;
t244 = t259 * t289 - t266 * t260;
t241 = t251 * t260 + t259 * t288;
t240 = t251 * t259 - t260 * t288;
t239 = t249 * t260 - t259 * t285;
t238 = t249 * t259 + t260 * t285;
t1 = [0, -g(1) * t275 - g(2) * t271, g(1) * t271 - g(2) * t275, 0, 0, 0, 0, 0, -g(1) * t251 - g(2) * t249 - g(3) * t289, -t237, 0, 0, 0, 0, 0, -g(1) * (t251 * t273 + t279) - g(2) * (t249 * t273 - t269 * t285) - g(3) * (t270 * t287 + t284) -g(1) * (-t251 * t269 + t271 * t287) - g(2) * (-t249 * t269 - t273 * t285) - g(3) * (t266 * t273 - t269 * t289) t237, -g(1) * t277 - g(2) * t276 - g(3) * t278, t237, g(1) * t241 + g(2) * t239 + g(3) * t245, -g(1) * t240 - g(2) * t238 - g(3) * t244, -g(1) * (t241 * pkin(4) + t240 * qJ(5) + t277) - g(2) * (t239 * pkin(4) + t238 * qJ(5) + t276) - g(3) * (t245 * pkin(4) + t244 * qJ(5) + t278) 0, 0, 0, 0, 0, -g(1) * (t240 * t268 + t250 * t272) - g(2) * (t238 * t268 + t248 * t272) - g(3) * (t244 * t268 - t272 * t286) -g(1) * (t240 * t272 - t250 * t268) - g(2) * (t238 * t272 - t248 * t268) - g(3) * (t244 * t272 + t268 * t286);];
U_reg  = t1;
