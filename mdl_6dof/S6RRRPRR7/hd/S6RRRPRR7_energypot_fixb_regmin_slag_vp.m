% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR7
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:25
% EndTime: 2019-03-09 18:39:26
% DurationCPUTime: 0.16s
% Computational Cost: add. (138->57), mult. (209->95), div. (0->0), fcn. (256->12), ass. (0->35)
t272 = sin(qJ(1));
t276 = cos(qJ(1));
t291 = -g(1) * t272 + g(2) * t276;
t266 = sin(pkin(6));
t271 = sin(qJ(2));
t288 = t266 * t271;
t287 = t266 * t272;
t274 = cos(qJ(3));
t286 = t266 * t274;
t275 = cos(qJ(2));
t285 = t266 * t275;
t284 = t266 * t276;
t267 = cos(pkin(6));
t270 = sin(qJ(3));
t283 = t267 * t270;
t282 = t272 * t271;
t281 = t272 * t275;
t280 = t276 * t271;
t279 = t276 * t275;
t257 = -t267 * t279 + t282;
t259 = t267 * t281 + t280;
t277 = -g(1) * t259 - g(2) * t257 + g(3) * t285;
t273 = cos(qJ(6));
t269 = sin(qJ(6));
t268 = -qJ(4) - pkin(9);
t265 = qJ(3) + pkin(12) + qJ(5);
t264 = t274 * pkin(3) + pkin(2);
t263 = cos(t265);
t262 = sin(t265);
t260 = -t267 * t282 + t279;
t258 = t267 * t280 + t281;
t256 = t267 * t262 + t263 * t288;
t255 = t260 * t263 + t262 * t287;
t254 = t258 * t263 - t262 * t284;
t1 = [0, -g(1) * t276 - g(2) * t272, -t291, 0, 0, 0, 0, 0, -g(1) * t260 - g(2) * t258 - g(3) * t288, -t277, 0, 0, 0, 0, 0, -g(1) * (t260 * t274 + t270 * t287) - g(2) * (t258 * t274 - t270 * t284) - g(3) * (t271 * t286 + t283) -g(1) * (-t260 * t270 + t272 * t286) - g(2) * (-t258 * t270 - t274 * t284) - g(3) * (t267 * t274 - t270 * t288) t277, -g(1) * (t276 * pkin(1) - t259 * t268 + t260 * t264) - g(2) * (t272 * pkin(1) - t257 * t268 + t258 * t264) - g(3) * (pkin(3) * t283 + t267 * pkin(8) + pkin(7)) + (-g(3) * (t264 * t271 + t268 * t275) + t291 * (pkin(3) * t270 + pkin(8))) * t266, 0, 0, 0, 0, 0, -g(1) * t255 - g(2) * t254 - g(3) * t256, -g(1) * (-t260 * t262 + t263 * t287) - g(2) * (-t258 * t262 - t263 * t284) - g(3) * (-t262 * t288 + t267 * t263) 0, 0, 0, 0, 0, -g(1) * (t255 * t273 + t259 * t269) - g(2) * (t254 * t273 + t257 * t269) - g(3) * (t256 * t273 - t269 * t285) -g(1) * (-t255 * t269 + t259 * t273) - g(2) * (-t254 * t269 + t257 * t273) - g(3) * (-t256 * t269 - t273 * t285);];
U_reg  = t1;
