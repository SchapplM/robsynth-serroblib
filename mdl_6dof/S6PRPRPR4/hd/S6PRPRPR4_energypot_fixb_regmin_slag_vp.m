% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:17
% EndTime: 2019-03-08 19:41:17
% DurationCPUTime: 0.18s
% Computational Cost: add. (217->84), mult. (352->132), div. (0->0), fcn. (434->14), ass. (0->42)
t231 = sin(pkin(11));
t252 = pkin(3) * t231;
t237 = cos(pkin(6));
t251 = t231 * t237;
t232 = sin(pkin(10));
t233 = sin(pkin(6));
t250 = t232 * t233;
t236 = cos(pkin(10));
t249 = t233 * t236;
t239 = sin(qJ(2));
t248 = t233 * t239;
t240 = cos(qJ(2));
t247 = t233 * t240;
t246 = t237 * t239;
t245 = t237 * t240;
t244 = t237 * pkin(7) + qJ(1);
t243 = t236 * pkin(1) + pkin(7) * t250;
t213 = t232 * t240 + t236 * t246;
t229 = pkin(11) + qJ(4);
t222 = sin(t229);
t224 = cos(t229);
t206 = t213 * t222 + t224 * t249;
t215 = -t232 * t246 + t236 * t240;
t208 = t215 * t222 - t224 * t250;
t210 = t222 * t248 - t237 * t224;
t242 = g(1) * t208 + g(2) * t206 + g(3) * t210;
t212 = t232 * t239 - t236 * t245;
t214 = t232 * t245 + t236 * t239;
t241 = -g(1) * t214 - g(2) * t212 + g(3) * t247;
t238 = -pkin(8) - qJ(3);
t235 = cos(pkin(11));
t234 = cos(pkin(12));
t230 = sin(pkin(12));
t228 = pkin(12) + qJ(6);
t225 = t232 * pkin(1);
t223 = cos(t228);
t221 = sin(t228);
t220 = pkin(3) * t235 + pkin(2);
t211 = t222 * t237 + t224 * t248;
t209 = t215 * t224 + t222 * t250;
t207 = t213 * t224 - t222 * t249;
t1 = [-g(3) * qJ(1), 0, -g(1) * t215 - g(2) * t213 - g(3) * t248, -t241, -g(1) * (t215 * t235 + t231 * t250) - g(2) * (t213 * t235 - t231 * t249) - g(3) * (t235 * t248 + t251) -g(1) * (-t215 * t231 + t235 * t250) - g(2) * (-t213 * t231 - t235 * t249) - g(3) * (-t231 * t248 + t235 * t237) t241, -g(1) * (pkin(2) * t215 + qJ(3) * t214 + t243) - g(2) * (pkin(2) * t213 - pkin(7) * t249 + qJ(3) * t212 + t225) - g(3) * ((pkin(2) * t239 - qJ(3) * t240) * t233 + t244) 0, 0, 0, 0, 0, -g(1) * t209 - g(2) * t207 - g(3) * t211, t242, -g(1) * (t209 * t234 + t214 * t230) - g(2) * (t207 * t234 + t212 * t230) - g(3) * (t211 * t234 - t230 * t247) -g(1) * (-t209 * t230 + t214 * t234) - g(2) * (-t207 * t230 + t212 * t234) - g(3) * (-t211 * t230 - t234 * t247) -t242, -g(1) * (pkin(4) * t209 + qJ(5) * t208 - t214 * t238 + t215 * t220 + t243) - g(2) * (t207 * pkin(4) + t206 * qJ(5) - t212 * t238 + t213 * t220 + t225) - g(3) * (pkin(3) * t251 + t211 * pkin(4) + t210 * qJ(5) + t244) + (-g(1) * t232 * t252 - g(3) * (t220 * t239 + t238 * t240) - g(2) * (-pkin(7) - t252) * t236) * t233, 0, 0, 0, 0, 0, -g(1) * (t209 * t223 + t214 * t221) - g(2) * (t207 * t223 + t212 * t221) - g(3) * (t211 * t223 - t221 * t247) -g(1) * (-t209 * t221 + t214 * t223) - g(2) * (-t207 * t221 + t212 * t223) - g(3) * (-t211 * t221 - t223 * t247);];
U_reg  = t1;
