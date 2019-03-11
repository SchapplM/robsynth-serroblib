% Calculate minimal parameter regressor of potential energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:14
% EndTime: 2019-03-08 19:16:14
% DurationCPUTime: 0.21s
% Computational Cost: add. (180->65), mult. (370->115), div. (0->0), fcn. (468->14), ass. (0->40)
t244 = pkin(7) + qJ(3);
t223 = sin(pkin(10));
t224 = sin(pkin(6));
t243 = t223 * t224;
t227 = cos(pkin(10));
t242 = t224 * t227;
t230 = sin(qJ(2));
t241 = t224 * t230;
t228 = cos(pkin(6));
t240 = t228 * t230;
t232 = cos(qJ(2));
t239 = t228 * t232;
t209 = pkin(2) * t240 - t244 * t224;
t215 = t232 * pkin(2) + pkin(1);
t238 = t227 * t209 + t223 * t215;
t237 = -t223 * t209 + t227 * t215;
t236 = pkin(2) * t241 + t244 * t228 + qJ(1);
t222 = sin(pkin(11));
t226 = cos(pkin(11));
t235 = t232 * t222 + t230 * t226;
t234 = t230 * t222 - t232 * t226;
t233 = t234 * t228;
t231 = cos(qJ(6));
t229 = sin(qJ(6));
t225 = cos(pkin(12));
t221 = sin(pkin(12));
t220 = pkin(12) + qJ(5);
t217 = cos(t220);
t216 = sin(t220);
t208 = t235 * t228;
t207 = t235 * t224;
t206 = t234 * t224;
t203 = t207 * t217 + t228 * t216;
t202 = -t223 * t208 - t227 * t234;
t201 = t223 * t233 - t227 * t235;
t200 = t227 * t208 - t223 * t234;
t199 = -t223 * t235 - t227 * t233;
t198 = t202 * t217 + t216 * t243;
t197 = t200 * t217 - t216 * t242;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t223 * t240 + t227 * t232) - g(2) * (t223 * t232 + t227 * t240) - g(3) * t241, -g(1) * (-t223 * t239 - t227 * t230) - g(2) * (-t223 * t230 + t227 * t239) - g(3) * t224 * t232, -g(1) * t237 - g(2) * t238 - g(3) * t236, -g(1) * (t202 * t225 + t221 * t243) - g(2) * (t200 * t225 - t221 * t242) - g(3) * (t207 * t225 + t228 * t221) -g(1) * (-t202 * t221 + t225 * t243) - g(2) * (-t200 * t221 - t225 * t242) - g(3) * (-t207 * t221 + t228 * t225) g(1) * t201 + g(2) * t199 - g(3) * t206, -g(1) * (t202 * pkin(3) - t201 * qJ(4) + t237) - g(2) * (t200 * pkin(3) - t199 * qJ(4) + t238) - g(3) * (t207 * pkin(3) + t206 * qJ(4) + t236) 0, 0, 0, 0, 0, -g(1) * t198 - g(2) * t197 - g(3) * t203, -g(1) * (-t202 * t216 + t217 * t243) - g(2) * (-t200 * t216 - t217 * t242) - g(3) * (-t207 * t216 + t228 * t217) 0, 0, 0, 0, 0, -g(1) * (t198 * t231 - t201 * t229) - g(2) * (t197 * t231 - t199 * t229) - g(3) * (t203 * t231 + t206 * t229) -g(1) * (-t198 * t229 - t201 * t231) - g(2) * (-t197 * t229 - t199 * t231) - g(3) * (-t203 * t229 + t206 * t231);];
U_reg  = t1;
