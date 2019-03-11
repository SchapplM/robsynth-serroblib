% Calculate minimal parameter regressor of potential energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:02
% EndTime: 2019-03-08 19:23:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (141->61), mult. (338->106), div. (0->0), fcn. (423->12), ass. (0->38)
t194 = sin(pkin(6));
t214 = pkin(7) * t194;
t199 = sin(qJ(5));
t213 = t194 * t199;
t200 = sin(qJ(2));
t212 = t194 * t200;
t202 = cos(qJ(5));
t211 = t194 * t202;
t203 = cos(qJ(2));
t210 = t194 * t203;
t197 = cos(pkin(6));
t209 = t197 * t200;
t208 = t197 * t203;
t207 = pkin(2) * t212 + t197 * pkin(7) + qJ(1);
t193 = sin(pkin(10));
t196 = cos(pkin(10));
t182 = t193 * t200 - t196 * t208;
t183 = t193 * t203 + t196 * t209;
t206 = t193 * pkin(1) + t183 * pkin(2) + t182 * qJ(3);
t184 = t193 * t208 + t196 * t200;
t185 = -t193 * t209 + t196 * t203;
t205 = t196 * pkin(1) + t185 * pkin(2) + t184 * qJ(3) + t193 * t214;
t204 = -g(1) * t184 - g(2) * t182 + g(3) * t210;
t201 = cos(qJ(6));
t198 = sin(qJ(6));
t195 = cos(pkin(11));
t192 = sin(pkin(11));
t179 = (-t192 * t203 + t195 * t200) * t194;
t178 = (t192 * t200 + t195 * t203) * t194;
t177 = t179 * t202 - t197 * t199;
t176 = t184 * t192 + t185 * t195;
t175 = -t184 * t195 + t185 * t192;
t174 = t182 * t192 + t183 * t195;
t173 = -t182 * t195 + t183 * t192;
t172 = -g(1) * t185 - g(2) * t183 - g(3) * t212;
t171 = t176 * t202 - t193 * t213;
t170 = t174 * t202 + t196 * t213;
t1 = [-g(3) * qJ(1), 0, t172, -t204, t172, t204, -g(1) * t205 - g(2) * (-t196 * t214 + t206) - g(3) * (-qJ(3) * t210 + t207) -g(1) * t176 - g(2) * t174 - g(3) * t179, g(1) * t175 + g(2) * t173 + g(3) * t178, -g(1) * (t185 * pkin(3) + t205) - g(2) * (t183 * pkin(3) + t206) - g(3) * (-t197 * qJ(4) + t207) + (g(1) * qJ(4) * t193 - g(3) * (pkin(3) * t200 - qJ(3) * t203) - g(2) * (-pkin(7) + qJ(4)) * t196) * t194, 0, 0, 0, 0, 0, -g(1) * t171 - g(2) * t170 - g(3) * t177, -g(1) * (-t176 * t199 - t193 * t211) - g(2) * (-t174 * t199 + t196 * t211) - g(3) * (-t179 * t199 - t197 * t202) 0, 0, 0, 0, 0, -g(1) * (t171 * t201 + t175 * t198) - g(2) * (t170 * t201 + t173 * t198) - g(3) * (t177 * t201 + t178 * t198) -g(1) * (-t171 * t198 + t175 * t201) - g(2) * (-t170 * t198 + t173 * t201) - g(3) * (-t177 * t198 + t178 * t201);];
U_reg  = t1;
