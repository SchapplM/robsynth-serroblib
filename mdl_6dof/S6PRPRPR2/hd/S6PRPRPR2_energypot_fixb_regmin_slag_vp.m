% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR2
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
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:53
% EndTime: 2019-03-08 19:32:53
% DurationCPUTime: 0.22s
% Computational Cost: add. (217->70), mult. (492->118), div. (0->0), fcn. (633->14), ass. (0->44)
t257 = pkin(7) + qJ(3);
t236 = sin(pkin(6));
t241 = sin(qJ(4));
t256 = t236 * t241;
t242 = sin(qJ(2));
t255 = t236 * t242;
t243 = cos(qJ(4));
t254 = t236 * t243;
t240 = cos(pkin(6));
t253 = t240 * t242;
t244 = cos(qJ(2));
t252 = t240 * t244;
t219 = pkin(2) * t253 - t257 * t236;
t227 = t244 * pkin(2) + pkin(1);
t235 = sin(pkin(10));
t239 = cos(pkin(10));
t251 = t239 * t219 + t235 * t227;
t250 = -t235 * t219 + t239 * t227;
t249 = pkin(2) * t255 + t257 * t240 + qJ(1);
t234 = sin(pkin(11));
t238 = cos(pkin(11));
t248 = t244 * t234 + t242 * t238;
t247 = t242 * t234 - t244 * t238;
t246 = t247 * t240;
t218 = t248 * t240;
t209 = t239 * t218 - t235 * t247;
t204 = t209 * t241 + t239 * t254;
t211 = -t235 * t218 - t239 * t247;
t206 = t211 * t241 - t235 * t254;
t217 = t248 * t236;
t212 = t217 * t241 - t240 * t243;
t245 = g(1) * t206 + g(2) * t204 + g(3) * t212;
t237 = cos(pkin(12));
t233 = sin(pkin(12));
t232 = pkin(12) + qJ(6);
t229 = cos(t232);
t228 = sin(t232);
t216 = t247 * t236;
t213 = t217 * t243 + t240 * t241;
t210 = t235 * t246 - t239 * t248;
t208 = -t235 * t248 - t239 * t246;
t207 = t211 * t243 + t235 * t256;
t205 = t209 * t243 - t239 * t256;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t235 * t253 + t239 * t244) - g(2) * (t235 * t244 + t239 * t253) - g(3) * t255, -g(1) * (-t235 * t252 - t239 * t242) - g(2) * (-t235 * t242 + t239 * t252) - g(3) * t236 * t244, -g(1) * t250 - g(2) * t251 - g(3) * t249, 0, 0, 0, 0, 0, -g(1) * t207 - g(2) * t205 - g(3) * t213, t245, -g(1) * (t207 * t237 - t210 * t233) - g(2) * (t205 * t237 - t208 * t233) - g(3) * (t213 * t237 + t216 * t233) -g(1) * (-t207 * t233 - t210 * t237) - g(2) * (-t205 * t233 - t208 * t237) - g(3) * (-t213 * t233 + t216 * t237) -t245, -g(1) * (t211 * pkin(3) + t207 * pkin(4) - t210 * pkin(8) + t206 * qJ(5) + t250) - g(2) * (t209 * pkin(3) + t205 * pkin(4) - t208 * pkin(8) + t204 * qJ(5) + t251) - g(3) * (t217 * pkin(3) + t213 * pkin(4) + t216 * pkin(8) + t212 * qJ(5) + t249) 0, 0, 0, 0, 0, -g(1) * (t207 * t229 - t210 * t228) - g(2) * (t205 * t229 - t208 * t228) - g(3) * (t213 * t229 + t216 * t228) -g(1) * (-t207 * t228 - t210 * t229) - g(2) * (-t205 * t228 - t208 * t229) - g(3) * (-t213 * t228 + t216 * t229);];
U_reg  = t1;
