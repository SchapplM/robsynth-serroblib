% Calculate minimal parameter regressor of potential energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:34
% EndTime: 2019-12-05 15:51:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (90->44), mult. (218->85), div. (0->0), fcn. (278->12), ass. (0->34)
t191 = pkin(6) + qJ(3);
t174 = sin(pkin(5));
t179 = sin(qJ(4));
t190 = t174 * t179;
t180 = sin(qJ(2));
t189 = t174 * t180;
t182 = cos(qJ(4));
t188 = t174 * t182;
t177 = cos(pkin(5));
t187 = t177 * t180;
t183 = cos(qJ(2));
t186 = t177 * t183;
t172 = sin(pkin(10));
t175 = cos(pkin(10));
t185 = t183 * t172 + t180 * t175;
t184 = t180 * t172 - t183 * t175;
t181 = cos(qJ(5));
t178 = sin(qJ(5));
t176 = cos(pkin(9));
t173 = sin(pkin(9));
t171 = t183 * pkin(2) + pkin(1);
t168 = pkin(2) * t187 - t191 * t174;
t167 = t185 * t177;
t166 = t184 * t177;
t165 = t185 * t174;
t164 = t184 * t174;
t163 = t165 * t182 + t177 * t179;
t162 = -t173 * t167 - t176 * t184;
t161 = -t173 * t166 + t176 * t185;
t160 = t176 * t167 - t173 * t184;
t159 = t176 * t166 + t173 * t185;
t158 = t162 * t182 + t173 * t190;
t157 = t160 * t182 - t176 * t190;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t173 * t187 + t176 * t183) - g(2) * (t173 * t183 + t176 * t187) - g(3) * t189, -g(1) * (-t173 * t186 - t176 * t180) - g(2) * (-t173 * t180 + t176 * t186) - g(3) * t174 * t183, -g(1) * (-t173 * t168 + t176 * t171) - g(2) * (t176 * t168 + t173 * t171) - g(3) * (pkin(2) * t189 + t191 * t177 + qJ(1)), 0, 0, 0, 0, 0, -g(1) * t158 - g(2) * t157 - g(3) * t163, -g(1) * (-t162 * t179 + t173 * t188) - g(2) * (-t160 * t179 - t176 * t188) - g(3) * (-t165 * t179 + t177 * t182), 0, 0, 0, 0, 0, -g(1) * (t158 * t181 + t161 * t178) - g(2) * (t157 * t181 + t159 * t178) - g(3) * (t163 * t181 + t164 * t178), -g(1) * (-t158 * t178 + t161 * t181) - g(2) * (-t157 * t178 + t159 * t181) - g(3) * (-t163 * t178 + t164 * t181);];
U_reg = t1;
