% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:24
% EndTime: 2019-03-09 11:09:24
% DurationCPUTime: 0.14s
% Computational Cost: add. (187->75), mult. (321->113), div. (0->0), fcn. (391->12), ass. (0->44)
t255 = sin(pkin(11));
t281 = pkin(3) * t255;
t262 = sin(qJ(1));
t280 = g(1) * t262;
t265 = cos(qJ(1));
t279 = g(2) * t265;
t258 = cos(pkin(6));
t278 = t258 * pkin(8) + pkin(7);
t256 = sin(pkin(6));
t261 = sin(qJ(2));
t277 = t256 * t261;
t276 = t256 * t262;
t264 = cos(qJ(2));
t275 = t256 * t264;
t274 = t256 * t265;
t273 = t258 * t255;
t272 = t262 * t261;
t271 = t262 * t264;
t270 = t265 * t261;
t269 = t265 * t264;
t268 = t265 * pkin(1) + pkin(8) * t276;
t241 = t258 * t270 + t271;
t254 = pkin(11) + qJ(4);
t249 = sin(t254);
t250 = cos(t254);
t234 = t241 * t249 + t250 * t274;
t243 = -t258 * t272 + t269;
t236 = t243 * t249 - t250 * t276;
t238 = t249 * t277 - t258 * t250;
t267 = g(1) * t236 + g(2) * t234 + g(3) * t238;
t235 = t241 * t250 - t249 * t274;
t237 = t243 * t250 + t249 * t276;
t239 = t258 * t249 + t250 * t277;
t266 = g(1) * t237 + g(2) * t235 + g(3) * t239;
t240 = -t258 * t269 + t272;
t242 = t258 * t271 + t270;
t233 = -g(1) * t242 - g(2) * t240 + g(3) * t275;
t263 = cos(qJ(6));
t260 = sin(qJ(6));
t259 = -pkin(9) - qJ(3);
t257 = cos(pkin(11));
t252 = t262 * pkin(1);
t248 = t257 * pkin(3) + pkin(2);
t1 = [0, -g(1) * t265 - g(2) * t262, -t279 + t280, 0, 0, 0, 0, 0, -g(1) * t243 - g(2) * t241 - g(3) * t277, -t233, -g(1) * (t243 * t257 + t255 * t276) - g(2) * (t241 * t257 - t255 * t274) - g(3) * (t257 * t277 + t273) -g(1) * (-t243 * t255 + t257 * t276) - g(2) * (-t241 * t255 - t257 * t274) - g(3) * (-t255 * t277 + t258 * t257) t233, -g(1) * (t243 * pkin(2) + t242 * qJ(3) + t268) - g(2) * (t241 * pkin(2) - pkin(8) * t274 + t240 * qJ(3) + t252) - g(3) * ((pkin(2) * t261 - qJ(3) * t264) * t256 + t278) 0, 0, 0, 0, 0, -t266, t267, t233, t266, -t267, -g(1) * (t237 * pkin(4) + t236 * qJ(5) - t242 * t259 + t243 * t248 + t268) - g(2) * (t235 * pkin(4) + t234 * qJ(5) - t240 * t259 + t241 * t248 + t252) - g(3) * (pkin(3) * t273 + t239 * pkin(4) + t238 * qJ(5) + t278) + (-t280 * t281 - g(3) * (t248 * t261 + t259 * t264) - (-pkin(8) - t281) * t279) * t256, 0, 0, 0, 0, 0, -g(1) * (t236 * t260 + t242 * t263) - g(2) * (t234 * t260 + t240 * t263) - g(3) * (t238 * t260 - t263 * t275) -g(1) * (t236 * t263 - t242 * t260) - g(2) * (t234 * t263 - t240 * t260) - g(3) * (t238 * t263 + t260 * t275);];
U_reg  = t1;
