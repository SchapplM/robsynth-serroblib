% Calculate minimal parameter regressor of potential energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:16
% EndTime: 2019-03-09 10:06:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (107->54), mult. (216->63), div. (0->0), fcn. (229->6), ass. (0->33)
t195 = sin(qJ(2));
t218 = qJ(3) * t195 + pkin(1);
t198 = cos(qJ(2));
t217 = g(3) * t198;
t199 = cos(qJ(1));
t216 = t199 * pkin(7);
t215 = t195 * pkin(2) + pkin(6);
t194 = sin(qJ(4));
t196 = sin(qJ(1));
t213 = t196 * t194;
t197 = cos(qJ(4));
t212 = t196 * t197;
t211 = t196 * t198;
t210 = t197 * t198;
t209 = t198 * t199;
t208 = t199 * t194;
t207 = t199 * t197;
t206 = pkin(2) * t211 + t218 * t196;
t205 = t195 * pkin(8) + qJ(5) * t210 + t215;
t204 = pkin(2) * t209 + t196 * pkin(7) + t218 * t199;
t203 = g(1) * t199 + g(2) * t196;
t177 = t195 * t212 + t208;
t178 = t195 * t213 - t207;
t202 = t178 * pkin(4) + pkin(8) * t211 - t177 * qJ(5) + t206;
t175 = -t195 * t207 + t213;
t201 = g(1) * t175 - g(2) * t177 + g(3) * t210;
t176 = t195 * t208 + t212;
t200 = t196 * pkin(3) + t176 * pkin(4) + pkin(8) * t209 + t175 * qJ(5) + t204;
t179 = g(1) * t196 - g(2) * t199;
t172 = g(3) * t195 + t203 * t198;
t171 = t203 * t195 - t217;
t170 = -g(1) * t176 - g(2) * t178 + t194 * t217;
t1 = [0, -t203, t179, 0, 0, 0, 0, 0, -t172, t171, -t179, t172, -t171, -g(1) * t204 - g(2) * (t206 - t216) - g(3) * (-t198 * qJ(3) + t215) 0, 0, 0, 0, 0, t170, t201, t170, -t172, -t201, -g(1) * t200 - g(2) * ((-pkin(3) - pkin(7)) * t199 + t202) - g(3) * ((-pkin(4) * t194 - qJ(3)) * t198 + t205) t170, -t201, t172, -g(1) * (t176 * pkin(5) + t200) - g(2) * (-t199 * pkin(3) + t178 * pkin(5) + t202 - t216) - g(3) * (-t195 * qJ(6) + t205) + (-g(3) * (-qJ(3) + (-pkin(4) - pkin(5)) * t194) + t203 * qJ(6)) * t198;];
U_reg  = t1;
