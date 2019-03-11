% Calculate minimal parameter regressor of potential energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:00
% EndTime: 2019-03-08 19:20:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (141->56), mult. (331->100), div. (0->0), fcn. (416->12), ass. (0->37)
t214 = pkin(7) + qJ(3);
t193 = sin(pkin(6));
t198 = sin(qJ(5));
t213 = t193 * t198;
t199 = sin(qJ(2));
t212 = t193 * t199;
t201 = cos(qJ(5));
t211 = t193 * t201;
t196 = cos(pkin(6));
t210 = t196 * t199;
t202 = cos(qJ(2));
t209 = t196 * t202;
t182 = pkin(2) * t210 - t214 * t193;
t188 = t202 * pkin(2) + pkin(1);
t192 = sin(pkin(10));
t195 = cos(pkin(10));
t208 = t195 * t182 + t192 * t188;
t207 = -t192 * t182 + t195 * t188;
t206 = pkin(2) * t212 + t214 * t196 + qJ(1);
t191 = sin(pkin(11));
t194 = cos(pkin(11));
t205 = t202 * t191 + t199 * t194;
t204 = t199 * t191 - t202 * t194;
t203 = t204 * t196;
t200 = cos(qJ(6));
t197 = sin(qJ(6));
t181 = t205 * t196;
t180 = t205 * t193;
t179 = t204 * t193;
t176 = t179 * t198 + t196 * t201;
t175 = -t192 * t181 - t195 * t204;
t174 = t192 * t203 - t195 * t205;
t173 = t195 * t181 - t192 * t204;
t172 = -t192 * t205 - t195 * t203;
t171 = -t172 * t198 - t195 * t211;
t170 = -t174 * t198 + t192 * t211;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t192 * t210 + t195 * t202) - g(2) * (t192 * t202 + t195 * t210) - g(3) * t212, -g(1) * (-t192 * t209 - t195 * t199) - g(2) * (-t192 * t199 + t195 * t209) - g(3) * t193 * t202, -g(1) * t207 - g(2) * t208 - g(3) * t206, g(1) * t175 + g(2) * t173 + g(3) * t180, g(1) * t174 + g(2) * t172 - g(3) * t179, -g(1) * (t175 * pkin(3) - t174 * qJ(4) + t207) - g(2) * (t173 * pkin(3) - t172 * qJ(4) + t208) - g(3) * (t180 * pkin(3) + t179 * qJ(4) + t206) 0, 0, 0, 0, 0, -g(1) * t170 - g(2) * t171 - g(3) * t176, -g(1) * (-t174 * t201 - t192 * t213) - g(2) * (-t172 * t201 + t195 * t213) - g(3) * (t179 * t201 - t196 * t198) 0, 0, 0, 0, 0, -g(1) * (t170 * t200 + t175 * t197) - g(2) * (t171 * t200 + t173 * t197) - g(3) * (t176 * t200 + t180 * t197) -g(1) * (-t170 * t197 + t175 * t200) - g(2) * (-t171 * t197 + t173 * t200) - g(3) * (-t176 * t197 + t180 * t200);];
U_reg  = t1;
