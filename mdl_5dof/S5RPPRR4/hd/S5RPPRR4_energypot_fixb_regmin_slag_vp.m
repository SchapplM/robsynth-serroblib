% Calculate minimal parameter regressor of potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:54
% EndTime: 2019-12-05 17:44:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->44), mult. (105->66), div. (0->0), fcn. (111->10), ass. (0->33)
t181 = sin(pkin(8));
t202 = g(1) * t181;
t184 = sin(qJ(1));
t201 = g(2) * t184;
t179 = pkin(9) + qJ(4);
t175 = qJ(5) + t179;
t171 = sin(t175);
t200 = t184 * t171;
t172 = cos(t175);
t199 = t184 * t172;
t173 = sin(t179);
t198 = t184 * t173;
t174 = cos(t179);
t197 = t184 * t174;
t180 = sin(pkin(9));
t196 = t184 * t180;
t182 = cos(pkin(9));
t195 = t184 * t182;
t185 = cos(qJ(1));
t194 = t185 * t171;
t193 = t185 * t172;
t192 = t185 * t173;
t191 = t185 * t174;
t190 = t185 * t180;
t189 = t185 * t182;
t188 = t185 * pkin(1) + t184 * qJ(2);
t187 = -g(3) * t185 + t201;
t183 = cos(pkin(8));
t186 = pkin(2) * t183 + qJ(3) * t181;
t177 = t185 * qJ(2);
t170 = g(2) * t185 + g(3) * t184;
t169 = g(1) * t183 + t187 * t181;
t1 = [0, t187, t170, t187 * t183 - t202, -t169, -t170, -g(1) * pkin(5) - g(2) * (-t184 * pkin(1) + t177) - g(3) * t188, -t182 * t202 - g(2) * (-t183 * t195 + t190) - g(3) * (t183 * t189 + t196), t180 * t202 - g(2) * (t183 * t196 + t189) - g(3) * (-t183 * t190 + t195), t169, -g(1) * (t181 * pkin(2) - t183 * qJ(3) + pkin(5)) - g(2) * t177 - g(3) * (t186 * t185 + t188) - (-pkin(1) - t186) * t201, 0, 0, 0, 0, 0, -t174 * t202 - g(2) * (-t183 * t197 + t192) - g(3) * (t183 * t191 + t198), t173 * t202 - g(2) * (t183 * t198 + t191) - g(3) * (-t183 * t192 + t197), 0, 0, 0, 0, 0, -t172 * t202 - g(2) * (-t183 * t199 + t194) - g(3) * (t183 * t193 + t200), t171 * t202 - g(2) * (t183 * t200 + t193) - g(3) * (-t183 * t194 + t199);];
U_reg = t1;
