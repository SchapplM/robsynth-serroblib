% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP1
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:28
% EndTime: 2019-03-08 19:58:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (177->62), mult. (407->102), div. (0->0), fcn. (515->12), ass. (0->44)
t242 = pkin(7) + qJ(3);
t218 = sin(pkin(6));
t224 = sin(qJ(4));
t241 = t218 * t224;
t225 = sin(qJ(2));
t240 = t218 * t225;
t227 = cos(qJ(4));
t239 = t218 * t227;
t221 = cos(pkin(6));
t238 = t221 * t225;
t228 = cos(qJ(2));
t237 = t221 * t228;
t204 = pkin(2) * t238 - t242 * t218;
t213 = t228 * pkin(2) + pkin(1);
t217 = sin(pkin(10));
t220 = cos(pkin(10));
t236 = t220 * t204 + t217 * t213;
t223 = sin(qJ(5));
t235 = pkin(5) * t223 + pkin(8);
t234 = -t217 * t204 + t220 * t213;
t233 = pkin(2) * t240 + t242 * t221 + qJ(1);
t216 = sin(pkin(11));
t219 = cos(pkin(11));
t232 = t228 * t216 + t225 * t219;
t231 = t225 * t216 - t228 * t219;
t230 = t231 * t221;
t203 = t232 * t221;
t194 = t220 * t203 - t217 * t231;
t189 = t194 * t224 + t220 * t239;
t196 = -t217 * t203 - t220 * t231;
t191 = t196 * t224 - t217 * t239;
t202 = t232 * t218;
t197 = t202 * t224 - t221 * t227;
t229 = g(1) * t191 + g(2) * t189 + g(3) * t197;
t226 = cos(qJ(5));
t222 = -qJ(6) - pkin(9);
t212 = t226 * pkin(5) + pkin(4);
t201 = t231 * t218;
t198 = t202 * t227 + t221 * t224;
t195 = t217 * t230 - t220 * t232;
t193 = -t217 * t232 - t220 * t230;
t192 = t196 * t227 + t217 * t241;
t190 = t194 * t227 - t220 * t241;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t217 * t238 + t220 * t228) - g(2) * (t217 * t228 + t220 * t238) - g(3) * t240, -g(1) * (-t217 * t237 - t220 * t225) - g(2) * (-t217 * t225 + t220 * t237) - g(3) * t218 * t228, -g(1) * t234 - g(2) * t236 - g(3) * t233, 0, 0, 0, 0, 0, -g(1) * t192 - g(2) * t190 - g(3) * t198, t229, 0, 0, 0, 0, 0, -g(1) * (t192 * t226 - t195 * t223) - g(2) * (t190 * t226 - t193 * t223) - g(3) * (t198 * t226 + t201 * t223) -g(1) * (-t192 * t223 - t195 * t226) - g(2) * (-t190 * t223 - t193 * t226) - g(3) * (-t198 * t223 + t201 * t226) -t229, -g(1) * (t196 * pkin(3) - t191 * t222 + t192 * t212 - t235 * t195 + t234) - g(2) * (t194 * pkin(3) - t189 * t222 + t190 * t212 - t235 * t193 + t236) - g(3) * (t202 * pkin(3) - t197 * t222 + t198 * t212 + t235 * t201 + t233);];
U_reg  = t1;
