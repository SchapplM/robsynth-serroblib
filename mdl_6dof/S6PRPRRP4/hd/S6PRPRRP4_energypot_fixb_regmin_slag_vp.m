% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:18
% EndTime: 2019-03-08 20:12:18
% DurationCPUTime: 0.17s
% Computational Cost: add. (243->80), mult. (414->120), div. (0->0), fcn. (516->12), ass. (0->47)
t227 = sin(pkin(11));
t250 = pkin(3) * t227;
t228 = sin(pkin(10));
t229 = sin(pkin(6));
t249 = t228 * t229;
t231 = cos(pkin(10));
t248 = t229 * t231;
t235 = sin(qJ(2));
t247 = t229 * t235;
t237 = cos(qJ(2));
t246 = t229 * t237;
t232 = cos(pkin(6));
t245 = t232 * t227;
t244 = t232 * t235;
t243 = t232 * t237;
t242 = t232 * pkin(7) + qJ(1);
t241 = t231 * pkin(1) + pkin(7) * t249;
t213 = t228 * t237 + t231 * t244;
t226 = pkin(11) + qJ(4);
t221 = sin(t226);
t222 = cos(t226);
t203 = t213 * t222 - t221 * t248;
t212 = t228 * t235 - t231 * t243;
t234 = sin(qJ(5));
t236 = cos(qJ(5));
t198 = t203 * t234 - t212 * t236;
t215 = -t228 * t244 + t231 * t237;
t205 = t215 * t222 + t221 * t249;
t214 = t228 * t243 + t231 * t235;
t200 = t205 * t234 - t214 * t236;
t209 = t232 * t221 + t222 * t247;
t206 = t209 * t234 + t236 * t246;
t240 = g(1) * t200 + g(2) * t198 + g(3) * t206;
t202 = t213 * t221 + t222 * t248;
t204 = t215 * t221 - t222 * t249;
t208 = t221 * t247 - t232 * t222;
t239 = g(1) * t204 + g(2) * t202 + g(3) * t208;
t238 = -g(1) * t214 - g(2) * t212 + g(3) * t246;
t233 = -pkin(8) - qJ(3);
t230 = cos(pkin(11));
t223 = t228 * pkin(1);
t220 = t230 * pkin(3) + pkin(2);
t207 = t209 * t236 - t234 * t246;
t201 = t205 * t236 + t214 * t234;
t199 = t203 * t236 + t212 * t234;
t197 = -g(1) * t201 - g(2) * t199 - g(3) * t207;
t1 = [-g(3) * qJ(1), 0, -g(1) * t215 - g(2) * t213 - g(3) * t247, -t238, -g(1) * (t215 * t230 + t227 * t249) - g(2) * (t213 * t230 - t227 * t248) - g(3) * (t230 * t247 + t245) -g(1) * (-t215 * t227 + t230 * t249) - g(2) * (-t213 * t227 - t230 * t248) - g(3) * (-t227 * t247 + t232 * t230) t238, -g(1) * (t215 * pkin(2) + t214 * qJ(3) + t241) - g(2) * (t213 * pkin(2) - pkin(7) * t248 + t212 * qJ(3) + t223) - g(3) * ((pkin(2) * t235 - qJ(3) * t237) * t229 + t242) 0, 0, 0, 0, 0, -g(1) * t205 - g(2) * t203 - g(3) * t209, t239, 0, 0, 0, 0, 0, t197, t240, t197, -t239, -t240, -g(1) * (t205 * pkin(4) + t201 * pkin(5) + t204 * pkin(9) + t200 * qJ(6) - t214 * t233 + t215 * t220 + t241) - g(2) * (t203 * pkin(4) + t199 * pkin(5) + t202 * pkin(9) + t198 * qJ(6) - t212 * t233 + t213 * t220 + t223) - g(3) * (pkin(3) * t245 + t209 * pkin(4) + t207 * pkin(5) + t208 * pkin(9) + t206 * qJ(6) + t242) + (-g(1) * t228 * t250 - g(3) * (t220 * t235 + t233 * t237) - g(2) * (-pkin(7) - t250) * t231) * t229;];
U_reg  = t1;
