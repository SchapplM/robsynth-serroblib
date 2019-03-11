% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:25
% EndTime: 2019-03-09 12:05:25
% DurationCPUTime: 0.19s
% Computational Cost: add. (180->64), mult. (415->106), div. (0->0), fcn. (524->12), ass. (0->47)
t296 = pkin(8) + qJ(3);
t268 = sin(pkin(6));
t274 = sin(qJ(2));
t295 = t268 * t274;
t275 = sin(qJ(1));
t294 = t268 * t275;
t279 = cos(qJ(1));
t293 = t268 * t279;
t292 = t275 * t274;
t278 = cos(qJ(2));
t291 = t275 * t278;
t290 = t279 * t274;
t289 = t279 * t278;
t270 = cos(pkin(6));
t255 = t270 * t274 * pkin(2) - t296 * t268;
t264 = t278 * pkin(2) + pkin(1);
t288 = t279 * t255 + t275 * t264;
t272 = sin(qJ(5));
t287 = pkin(5) * t272 + pkin(9);
t286 = -t275 * t255 + t279 * t264;
t285 = pkin(2) * t295 + t296 * t270 + pkin(7);
t284 = g(1) * t275 - g(2) * t279;
t267 = sin(pkin(11));
t269 = cos(pkin(11));
t283 = t278 * t267 + t274 * t269;
t282 = t274 * t267 - t278 * t269;
t281 = t282 * t270;
t254 = t283 * t270;
t245 = t279 * t254 - t275 * t282;
t273 = sin(qJ(4));
t277 = cos(qJ(4));
t240 = t245 * t273 + t277 * t293;
t247 = -t275 * t254 - t279 * t282;
t242 = t247 * t273 - t277 * t294;
t253 = t283 * t268;
t248 = t253 * t273 - t270 * t277;
t280 = g(1) * t242 + g(2) * t240 + g(3) * t248;
t276 = cos(qJ(5));
t271 = -qJ(6) - pkin(10);
t263 = t276 * pkin(5) + pkin(4);
t252 = t282 * t268;
t249 = t253 * t277 + t270 * t273;
t246 = t275 * t281 - t279 * t283;
t244 = -t275 * t283 - t279 * t281;
t243 = t247 * t277 + t273 * t294;
t241 = t245 * t277 - t273 * t293;
t1 = [0, -g(1) * t279 - g(2) * t275, t284, 0, 0, 0, 0, 0, -g(1) * (-t270 * t292 + t289) - g(2) * (t270 * t290 + t291) - g(3) * t295, -g(1) * (-t270 * t291 - t290) - g(2) * (t270 * t289 - t292) - g(3) * t268 * t278, -g(3) * t270 - t284 * t268, -g(1) * t286 - g(2) * t288 - g(3) * t285, 0, 0, 0, 0, 0, -g(1) * t243 - g(2) * t241 - g(3) * t249, t280, 0, 0, 0, 0, 0, -g(1) * (t243 * t276 - t246 * t272) - g(2) * (t241 * t276 - t244 * t272) - g(3) * (t249 * t276 + t252 * t272) -g(1) * (-t243 * t272 - t246 * t276) - g(2) * (-t241 * t272 - t244 * t276) - g(3) * (-t249 * t272 + t252 * t276) -t280, -g(1) * (t247 * pkin(3) - t242 * t271 + t243 * t263 - t287 * t246 + t286) - g(2) * (t245 * pkin(3) - t240 * t271 + t241 * t263 - t287 * t244 + t288) - g(3) * (t253 * pkin(3) - t248 * t271 + t249 * t263 + t287 * t252 + t285);];
U_reg  = t1;
