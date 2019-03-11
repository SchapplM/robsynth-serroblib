% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:30
% EndTime: 2019-03-08 21:06:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (169->65), mult. (293->104), div. (0->0), fcn. (352->12), ass. (0->40)
t230 = sin(pkin(10));
t231 = sin(pkin(6));
t253 = t230 * t231;
t232 = cos(pkin(10));
t252 = t231 * t232;
t236 = sin(qJ(3));
t251 = t231 * t236;
t237 = sin(qJ(2));
t250 = t231 * t237;
t239 = cos(qJ(3));
t249 = t231 * t239;
t240 = cos(qJ(2));
t248 = t231 * t240;
t233 = cos(pkin(6));
t247 = t233 * t236;
t246 = t233 * t237;
t245 = t233 * t240;
t244 = t230 * t251;
t215 = t230 * t245 + t232 * t237;
t216 = -t230 * t246 + t232 * t240;
t223 = t239 * pkin(3) + pkin(2);
t234 = -qJ(4) - pkin(8);
t243 = t232 * pkin(1) + pkin(3) * t244 + pkin(7) * t253 - t215 * t234 + t216 * t223;
t242 = pkin(3) * t247 + t233 * pkin(7) + t223 * t250 + t234 * t248 + qJ(1);
t213 = t230 * t237 - t232 * t245;
t202 = -g(1) * t215 - g(2) * t213 + g(3) * t248;
t214 = t230 * t240 + t232 * t246;
t241 = t214 * t223 - t213 * t234 + t230 * pkin(1) + (-pkin(3) * t236 - pkin(7)) * t252;
t238 = cos(qJ(6));
t235 = sin(qJ(6));
t229 = qJ(3) + pkin(11);
t225 = cos(t229);
t224 = sin(t229);
t210 = t233 * t224 + t225 * t250;
t209 = t224 * t250 - t233 * t225;
t206 = t216 * t225 + t224 * t253;
t205 = t216 * t224 - t225 * t253;
t204 = t214 * t225 - t224 * t252;
t203 = t214 * t224 + t225 * t252;
t1 = [-g(3) * qJ(1), 0, -g(1) * t216 - g(2) * t214 - g(3) * t250, -t202, 0, 0, 0, 0, 0, -g(1) * (t216 * t239 + t244) - g(2) * (t214 * t239 - t232 * t251) - g(3) * (t237 * t249 + t247) -g(1) * (-t216 * t236 + t230 * t249) - g(2) * (-t214 * t236 - t232 * t249) - g(3) * (t233 * t239 - t236 * t250) t202, -g(1) * t243 - g(2) * t241 - g(3) * t242, t202, g(1) * t206 + g(2) * t204 + g(3) * t210, -g(1) * t205 - g(2) * t203 - g(3) * t209, -g(1) * (t206 * pkin(4) + t205 * qJ(5) + t243) - g(2) * (t204 * pkin(4) + t203 * qJ(5) + t241) - g(3) * (t210 * pkin(4) + t209 * qJ(5) + t242) 0, 0, 0, 0, 0, -g(1) * (t205 * t235 + t215 * t238) - g(2) * (t203 * t235 + t213 * t238) - g(3) * (t209 * t235 - t238 * t248) -g(1) * (t205 * t238 - t215 * t235) - g(2) * (t203 * t238 - t213 * t235) - g(3) * (t209 * t238 + t235 * t248);];
U_reg  = t1;
