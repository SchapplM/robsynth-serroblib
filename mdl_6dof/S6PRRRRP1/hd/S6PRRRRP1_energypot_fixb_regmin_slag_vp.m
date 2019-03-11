% Calculate minimal parameter regressor of potential energy for
% S6PRRRRP1
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
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:17
% EndTime: 2019-03-08 23:59:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (161->67), mult. (266->107), div. (0->0), fcn. (327->12), ass. (0->40)
t214 = sin(pkin(11));
t216 = cos(pkin(11));
t221 = sin(qJ(2));
t217 = cos(pkin(6));
t224 = cos(qJ(2));
t228 = t217 * t224;
t202 = t214 * t221 - t216 * t228;
t219 = sin(qJ(5));
t238 = t202 * t219;
t204 = t214 * t228 + t216 * t221;
t237 = t204 * t219;
t215 = sin(pkin(6));
t236 = t214 * t215;
t235 = t215 * t216;
t220 = sin(qJ(3));
t234 = t215 * t220;
t233 = t215 * t221;
t223 = cos(qJ(3));
t232 = t215 * t223;
t231 = t215 * t224;
t230 = t217 * t220;
t229 = t217 * t221;
t203 = t214 * t224 + t216 * t229;
t213 = qJ(3) + qJ(4);
t211 = sin(t213);
t212 = cos(t213);
t196 = t203 * t211 + t212 * t235;
t205 = -t214 * t229 + t216 * t224;
t198 = t205 * t211 - t212 * t236;
t200 = t211 * t233 - t217 * t212;
t226 = g(1) * t198 + g(2) * t196 + g(3) * t200;
t225 = -pkin(9) - pkin(8);
t222 = cos(qJ(5));
t218 = -qJ(6) - pkin(10);
t210 = t223 * pkin(3) + pkin(2);
t209 = t222 * pkin(5) + pkin(4);
t201 = t217 * t211 + t212 * t233;
t199 = t205 * t212 + t211 * t236;
t197 = t203 * t212 - t211 * t235;
t1 = [-g(3) * qJ(1), 0, -g(1) * t205 - g(2) * t203 - g(3) * t233, g(1) * t204 + g(2) * t202 - g(3) * t231, 0, 0, 0, 0, 0, -g(1) * (t205 * t223 + t214 * t234) - g(2) * (t203 * t223 - t216 * t234) - g(3) * (t221 * t232 + t230) -g(1) * (-t205 * t220 + t214 * t232) - g(2) * (-t203 * t220 - t216 * t232) - g(3) * (t217 * t223 - t220 * t233) 0, 0, 0, 0, 0, -g(1) * t199 - g(2) * t197 - g(3) * t201, t226, 0, 0, 0, 0, 0, -g(1) * (t199 * t222 + t237) - g(2) * (t197 * t222 + t238) - g(3) * (t201 * t222 - t219 * t231) -g(1) * (-t199 * t219 + t204 * t222) - g(2) * (-t197 * t219 + t202 * t222) - g(3) * (-t201 * t219 - t222 * t231) -t226, -g(1) * (t216 * pkin(1) + pkin(5) * t237 - t198 * t218 + t199 * t209 - t204 * t225 + t205 * t210) - g(2) * (t214 * pkin(1) + pkin(5) * t238 - t196 * t218 + t197 * t209 - t202 * t225 + t203 * t210) - g(3) * (pkin(3) * t230 + t217 * pkin(7) - t200 * t218 + t201 * t209 + qJ(1)) + (-g(3) * (t210 * t221 + (-pkin(5) * t219 + t225) * t224) + (-g(1) * t214 + g(2) * t216) * (pkin(3) * t220 + pkin(7))) * t215;];
U_reg  = t1;
