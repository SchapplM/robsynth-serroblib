% Calculate minimal parameter regressor of potential energy for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:19
% EndTime: 2019-03-09 00:04:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (224->69), mult. (374->109), div. (0->0), fcn. (473->12), ass. (0->44)
t229 = sin(pkin(11));
t230 = sin(pkin(6));
t251 = t229 * t230;
t231 = cos(pkin(11));
t250 = t230 * t231;
t234 = sin(qJ(3));
t249 = t230 * t234;
t235 = sin(qJ(2));
t248 = t230 * t235;
t237 = cos(qJ(3));
t247 = t230 * t237;
t238 = cos(qJ(2));
t246 = t230 * t238;
t232 = cos(pkin(6));
t245 = t232 * t234;
t244 = t232 * t235;
t243 = t232 * t238;
t219 = t229 * t238 + t231 * t244;
t228 = qJ(3) + qJ(4);
t226 = sin(t228);
t227 = cos(t228);
t209 = t219 * t227 - t226 * t250;
t218 = t229 * t235 - t231 * t243;
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t204 = t209 * t233 - t218 * t236;
t221 = -t229 * t244 + t231 * t238;
t211 = t221 * t227 + t226 * t251;
t220 = t229 * t243 + t231 * t235;
t206 = t211 * t233 - t220 * t236;
t215 = t232 * t226 + t227 * t248;
t212 = t215 * t233 + t236 * t246;
t241 = g(1) * t206 + g(2) * t204 + g(3) * t212;
t208 = t219 * t226 + t227 * t250;
t210 = t221 * t226 - t227 * t251;
t214 = t226 * t248 - t232 * t227;
t240 = g(1) * t210 + g(2) * t208 + g(3) * t214;
t239 = -pkin(9) - pkin(8);
t225 = t237 * pkin(3) + pkin(2);
t213 = t215 * t236 - t233 * t246;
t207 = t211 * t236 + t220 * t233;
t205 = t209 * t236 + t218 * t233;
t203 = -g(1) * t207 - g(2) * t205 - g(3) * t213;
t1 = [-g(3) * qJ(1), 0, -g(1) * t221 - g(2) * t219 - g(3) * t248, g(1) * t220 + g(2) * t218 - g(3) * t246, 0, 0, 0, 0, 0, -g(1) * (t221 * t237 + t229 * t249) - g(2) * (t219 * t237 - t231 * t249) - g(3) * (t235 * t247 + t245) -g(1) * (-t221 * t234 + t229 * t247) - g(2) * (-t219 * t234 - t231 * t247) - g(3) * (t232 * t237 - t234 * t248) 0, 0, 0, 0, 0, -g(1) * t211 - g(2) * t209 - g(3) * t215, t240, 0, 0, 0, 0, 0, t203, t241, t203, -t240, -t241, -g(1) * (t231 * pkin(1) + t211 * pkin(4) + t207 * pkin(5) + t210 * pkin(10) + t206 * qJ(6) - t220 * t239 + t221 * t225) - g(2) * (t229 * pkin(1) + t209 * pkin(4) + t205 * pkin(5) + t208 * pkin(10) + t204 * qJ(6) - t218 * t239 + t219 * t225) - g(3) * (pkin(3) * t245 + t215 * pkin(4) + t213 * pkin(5) + t232 * pkin(7) + t214 * pkin(10) + t212 * qJ(6) + qJ(1)) + (-g(3) * (t225 * t235 + t238 * t239) + (-g(1) * t229 + g(2) * t231) * (pkin(3) * t234 + pkin(7))) * t230;];
U_reg  = t1;
