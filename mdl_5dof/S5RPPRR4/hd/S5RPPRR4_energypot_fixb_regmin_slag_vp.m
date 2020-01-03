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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:31:21
% EndTime: 2020-01-03 11:31:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->45), mult. (105->69), div. (0->0), fcn. (111->10), ass. (0->31)
t178 = sin(pkin(8));
t197 = g(1) * t178;
t176 = pkin(9) + qJ(4);
t174 = qJ(5) + t176;
t170 = sin(t174);
t181 = sin(qJ(1));
t196 = t181 * t170;
t171 = cos(t174);
t195 = t181 * t171;
t172 = sin(t176);
t194 = t181 * t172;
t173 = cos(t176);
t193 = t181 * t173;
t177 = sin(pkin(9));
t192 = t181 * t177;
t179 = cos(pkin(9));
t191 = t181 * t179;
t182 = cos(qJ(1));
t190 = t182 * t170;
t189 = t182 * t171;
t188 = t182 * t172;
t187 = t182 * t173;
t186 = t182 * t177;
t185 = t182 * t179;
t184 = -g(2) * t181 + g(3) * t182;
t180 = cos(pkin(8));
t183 = pkin(2) * t180 + qJ(3) * t178;
t175 = t181 * pkin(1);
t169 = g(2) * t182 + g(3) * t181;
t168 = g(1) * t180 + t184 * t178;
t1 = [0, t184, -t169, t184 * t180 - t197, -t168, t169, -g(1) * pkin(5) - g(2) * (-t182 * qJ(2) + t175) - g(3) * (-t182 * pkin(1) - t181 * qJ(2)), -t179 * t197 - g(2) * (t180 * t191 - t186) - g(3) * (-t180 * t185 - t192), t177 * t197 - g(2) * (-t180 * t192 - t185) - g(3) * (t180 * t186 - t191), t168, -g(1) * (t178 * pkin(2) - t180 * qJ(3) + pkin(5)) - g(2) * t175 + (-g(2) * t183 + g(3) * qJ(2)) * t181 + (g(2) * qJ(2) - g(3) * (-pkin(1) - t183)) * t182, 0, 0, 0, 0, 0, -t173 * t197 - g(2) * (t180 * t193 - t188) - g(3) * (-t180 * t187 - t194), t172 * t197 - g(2) * (-t180 * t194 - t187) - g(3) * (t180 * t188 - t193), 0, 0, 0, 0, 0, -t171 * t197 - g(2) * (t180 * t195 - t190) - g(3) * (-t180 * t189 - t196), t170 * t197 - g(2) * (-t180 * t196 - t189) - g(3) * (t180 * t190 - t195);];
U_reg = t1;
