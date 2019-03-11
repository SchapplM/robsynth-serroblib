% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:43
% EndTime: 2019-03-09 01:53:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (97->50), mult. (115->66), div. (0->0), fcn. (114->10), ass. (0->32)
t203 = pkin(2) + pkin(6);
t182 = pkin(9) + qJ(4);
t177 = cos(t182);
t202 = g(3) * t177;
t181 = pkin(10) + qJ(6);
t174 = sin(t181);
t188 = sin(qJ(1));
t201 = t188 * t174;
t176 = cos(t181);
t200 = t188 * t176;
t183 = sin(pkin(10));
t199 = t188 * t183;
t185 = cos(pkin(10));
t198 = t188 * t185;
t189 = cos(qJ(1));
t197 = t189 * t174;
t196 = t189 * t176;
t195 = t189 * t183;
t194 = t189 * t185;
t193 = t189 * pkin(1) + t188 * qJ(2);
t192 = g(1) * t193;
t179 = t188 * pkin(1);
t191 = -t189 * qJ(2) + t179;
t172 = g(1) * t188 - g(2) * t189;
t175 = sin(t182);
t184 = sin(pkin(9));
t190 = pkin(3) * t184 + pkin(4) * t175 - qJ(5) * t177;
t187 = -pkin(7) - qJ(3);
t186 = cos(pkin(9));
t173 = g(1) * t189 + g(2) * t188;
t171 = -g(3) * t175 + t172 * t177;
t1 = [0, -t173, t172, t173, -t172, -g(3) * pkin(6) - g(2) * t191 - t192, -g(3) * t186 - t172 * t184, g(3) * t184 - t172 * t186, -t173, -g(1) * (t189 * qJ(3) + t193) - g(2) * (t188 * qJ(3) + t191) - g(3) * t203, 0, 0, 0, 0, 0, -t172 * t175 - t202, -t171, -g(1) * (t175 * t198 + t195) - g(2) * (-t175 * t194 + t199) - t185 * t202, -g(1) * (-t175 * t199 + t194) - g(2) * (t175 * t195 + t198) + t183 * t202, t171, -t192 - g(2) * t179 - g(3) * (t186 * pkin(3) + t177 * pkin(4) + t175 * qJ(5) + t203) + (-g(1) * t190 + g(2) * t187) * t188 + (g(1) * t187 - g(2) * (-qJ(2) - t190)) * t189, 0, 0, 0, 0, 0, -g(1) * (t175 * t200 + t197) - g(2) * (-t175 * t196 + t201) - t176 * t202, -g(1) * (-t175 * t201 + t196) - g(2) * (t175 * t197 + t200) + t174 * t202;];
U_reg  = t1;
