% Calculate minimal parameter regressor of potential energy for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:46
% EndTime: 2019-12-31 22:25:46
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->25), mult. (64->40), div. (0->0), fcn. (72->10), ass. (0->23)
t173 = qJ(2) + qJ(3);
t169 = sin(t173);
t189 = g(3) * t169;
t172 = qJ(4) + qJ(5);
t168 = sin(t172);
t176 = sin(qJ(1));
t188 = t176 * t168;
t170 = cos(t172);
t187 = t176 * t170;
t174 = sin(qJ(4));
t186 = t176 * t174;
t177 = cos(qJ(4));
t185 = t176 * t177;
t179 = cos(qJ(1));
t184 = t179 * t168;
t183 = t179 * t170;
t182 = t179 * t174;
t181 = t179 * t177;
t180 = g(1) * t179 + g(2) * t176;
t178 = cos(qJ(2));
t175 = sin(qJ(2));
t171 = cos(t173);
t1 = [0, -t180, g(1) * t176 - g(2) * t179, 0, 0, 0, 0, 0, -g(3) * t175 - t180 * t178, -g(3) * t178 + t180 * t175, 0, 0, 0, 0, 0, -t180 * t171 - t189, -g(3) * t171 + t180 * t169, 0, 0, 0, 0, 0, -g(1) * (t171 * t181 + t186) - g(2) * (t171 * t185 - t182) - t177 * t189, -g(1) * (-t171 * t182 + t185) - g(2) * (-t171 * t186 - t181) + t174 * t189, 0, 0, 0, 0, 0, -g(1) * (t171 * t183 + t188) - g(2) * (t171 * t187 - t184) - t170 * t189, -g(1) * (-t171 * t184 + t187) - g(2) * (-t171 * t188 - t183) + t168 * t189;];
U_reg = t1;
