% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP2
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
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:05
% EndTime: 2019-03-08 20:03:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (251->66), mult. (606->106), div. (0->0), fcn. (787->12), ass. (0->49)
t257 = pkin(7) + qJ(3);
t235 = sin(pkin(6));
t240 = sin(qJ(4));
t256 = t235 * t240;
t241 = sin(qJ(2));
t255 = t235 * t241;
t243 = cos(qJ(4));
t254 = t235 * t243;
t238 = cos(pkin(6));
t253 = t238 * t241;
t244 = cos(qJ(2));
t252 = t238 * t244;
t222 = pkin(2) * t253 - t257 * t235;
t230 = t244 * pkin(2) + pkin(1);
t234 = sin(pkin(10));
t237 = cos(pkin(10));
t251 = t237 * t222 + t234 * t230;
t250 = -t234 * t222 + t237 * t230;
t249 = pkin(2) * t255 + t257 * t238 + qJ(1);
t233 = sin(pkin(11));
t236 = cos(pkin(11));
t248 = t244 * t233 + t241 * t236;
t247 = t241 * t233 - t244 * t236;
t221 = t248 * t238;
t211 = t237 * t221 - t234 * t247;
t205 = t211 * t243 - t237 * t256;
t220 = t247 * t238;
t210 = -t237 * t220 - t234 * t248;
t239 = sin(qJ(5));
t242 = cos(qJ(5));
t198 = t205 * t239 + t210 * t242;
t213 = -t234 * t221 - t237 * t247;
t207 = t213 * t243 + t234 * t256;
t212 = t234 * t220 - t237 * t248;
t200 = t207 * t239 + t212 * t242;
t219 = t248 * t235;
t215 = t219 * t243 + t238 * t240;
t218 = t247 * t235;
t202 = t215 * t239 - t218 * t242;
t246 = g(1) * t200 + g(2) * t198 + g(3) * t202;
t204 = t211 * t240 + t237 * t254;
t206 = t213 * t240 - t234 * t254;
t214 = t219 * t240 - t238 * t243;
t245 = g(1) * t206 + g(2) * t204 + g(3) * t214;
t203 = t215 * t242 + t218 * t239;
t201 = t207 * t242 - t212 * t239;
t199 = t205 * t242 - t210 * t239;
t197 = -g(1) * t201 - g(2) * t199 - g(3) * t203;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t234 * t253 + t237 * t244) - g(2) * (t234 * t244 + t237 * t253) - g(3) * t255, -g(1) * (-t234 * t252 - t237 * t241) - g(2) * (-t234 * t241 + t237 * t252) - g(3) * t235 * t244, -g(1) * t250 - g(2) * t251 - g(3) * t249, 0, 0, 0, 0, 0, -g(1) * t207 - g(2) * t205 - g(3) * t215, t245, 0, 0, 0, 0, 0, t197, t246, t197, -t245, -t246, -g(1) * (t213 * pkin(3) + t207 * pkin(4) + t201 * pkin(5) - t212 * pkin(8) + t206 * pkin(9) + t200 * qJ(6) + t250) - g(2) * (t211 * pkin(3) + t205 * pkin(4) + t199 * pkin(5) - t210 * pkin(8) + t204 * pkin(9) + t198 * qJ(6) + t251) - g(3) * (t219 * pkin(3) + t215 * pkin(4) + t203 * pkin(5) + t218 * pkin(8) + t214 * pkin(9) + t202 * qJ(6) + t249);];
U_reg  = t1;
