% Calculate minimal parameter regressor of potential energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:29
% EndTime: 2019-03-08 19:01:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (200->58), mult. (483->111), div. (0->0), fcn. (633->16), ass. (0->43)
t252 = sin(pkin(12));
t258 = cos(pkin(6));
t273 = t252 * t258;
t253 = sin(pkin(7));
t254 = sin(pkin(6));
t272 = t253 * t254;
t271 = t253 * t258;
t256 = cos(pkin(12));
t270 = t254 * t256;
t257 = cos(pkin(7));
t269 = t254 * t257;
t255 = cos(pkin(13));
t268 = t255 * t257;
t267 = t256 * t258;
t251 = sin(pkin(13));
t244 = -t252 * t251 + t255 * t267;
t266 = -t244 * t257 + t253 * t270;
t246 = -t256 * t251 - t255 * t273;
t265 = t246 * t257 + t252 * t272;
t264 = cos(qJ(3));
t263 = cos(qJ(4));
t262 = cos(qJ(6));
t261 = sin(qJ(3));
t260 = sin(qJ(4));
t259 = sin(qJ(6));
t250 = qJ(4) + qJ(5);
t249 = cos(t250);
t248 = sin(t250);
t247 = -t251 * t273 + t256 * t255;
t245 = t251 * t267 + t252 * t255;
t243 = -t255 * t272 + t258 * t257;
t242 = -t246 * t253 + t252 * t269;
t241 = -t244 * t253 - t256 * t269;
t240 = t261 * t271 + (t251 * t264 + t261 * t268) * t254;
t239 = -t264 * t271 + (t251 * t261 - t264 * t268) * t254;
t238 = t247 * t264 + t265 * t261;
t237 = t247 * t261 - t265 * t264;
t236 = t245 * t264 - t266 * t261;
t235 = t245 * t261 + t266 * t264;
t234 = t240 * t249 + t243 * t248;
t233 = t238 * t249 + t242 * t248;
t232 = t236 * t249 + t241 * t248;
t1 = [-g(3) * qJ(1), -g(1) * (t252 * t254 * qJ(2) + t256 * pkin(1)) - g(2) * (t252 * pkin(1) - qJ(2) * t270) - g(3) * (t258 * qJ(2) + qJ(1)) 0, -g(1) * t238 - g(2) * t236 - g(3) * t240, g(1) * t237 + g(2) * t235 + g(3) * t239, 0, 0, 0, 0, 0, -g(1) * (t238 * t263 + t242 * t260) - g(2) * (t236 * t263 + t241 * t260) - g(3) * (t240 * t263 + t243 * t260) -g(1) * (-t238 * t260 + t242 * t263) - g(2) * (-t236 * t260 + t241 * t263) - g(3) * (-t240 * t260 + t243 * t263) 0, 0, 0, 0, 0, -g(1) * t233 - g(2) * t232 - g(3) * t234, -g(1) * (-t238 * t248 + t242 * t249) - g(2) * (-t236 * t248 + t241 * t249) - g(3) * (-t240 * t248 + t243 * t249) 0, 0, 0, 0, 0, -g(1) * (t233 * t262 + t237 * t259) - g(2) * (t232 * t262 + t235 * t259) - g(3) * (t234 * t262 + t239 * t259) -g(1) * (-t233 * t259 + t237 * t262) - g(2) * (-t232 * t259 + t235 * t262) - g(3) * (-t234 * t259 + t239 * t262);];
U_reg  = t1;
