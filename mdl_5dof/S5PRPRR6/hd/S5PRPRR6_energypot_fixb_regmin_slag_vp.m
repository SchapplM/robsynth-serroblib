% Calculate minimal parameter regressor of potential energy for
% S5PRPRR6
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:53
% EndTime: 2019-12-05 15:57:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (104->52), mult. (195->90), div. (0->0), fcn. (241->12), ass. (0->28)
t171 = sin(pkin(9));
t172 = sin(pkin(5));
t186 = t171 * t172;
t174 = cos(pkin(9));
t185 = t172 * t174;
t177 = sin(qJ(2));
t184 = t172 * t177;
t179 = cos(qJ(2));
t183 = t172 * t179;
t175 = cos(pkin(5));
t182 = t175 * t177;
t181 = t175 * t179;
t162 = t171 * t177 - t174 * t181;
t164 = t171 * t181 + t174 * t177;
t180 = -g(1) * t164 - g(2) * t162 + g(3) * t183;
t178 = cos(qJ(5));
t176 = sin(qJ(5));
t173 = cos(pkin(10));
t170 = sin(pkin(10));
t169 = pkin(10) + qJ(4);
t168 = cos(t169);
t167 = sin(t169);
t165 = -t171 * t182 + t174 * t179;
t163 = t171 * t179 + t174 * t182;
t161 = t175 * t167 + t168 * t184;
t160 = t165 * t168 + t167 * t186;
t159 = t163 * t168 - t167 * t185;
t1 = [-g(3) * qJ(1), 0, -g(1) * t165 - g(2) * t163 - g(3) * t184, -t180, -g(1) * (t165 * t173 + t170 * t186) - g(2) * (t163 * t173 - t170 * t185) - g(3) * (t175 * t170 + t173 * t184), -g(1) * (-t165 * t170 + t173 * t186) - g(2) * (-t163 * t170 - t173 * t185) - g(3) * (-t170 * t184 + t175 * t173), t180, -g(1) * (t174 * pkin(1) + t165 * pkin(2) + pkin(6) * t186 + t164 * qJ(3)) - g(2) * (t171 * pkin(1) + t163 * pkin(2) - pkin(6) * t185 + t162 * qJ(3)) - g(3) * (t175 * pkin(6) + qJ(1) + (pkin(2) * t177 - qJ(3) * t179) * t172), 0, 0, 0, 0, 0, -g(1) * t160 - g(2) * t159 - g(3) * t161, -g(1) * (-t165 * t167 + t168 * t186) - g(2) * (-t163 * t167 - t168 * t185) - g(3) * (-t167 * t184 + t175 * t168), 0, 0, 0, 0, 0, -g(1) * (t160 * t178 + t164 * t176) - g(2) * (t159 * t178 + t162 * t176) - g(3) * (t161 * t178 - t176 * t183), -g(1) * (-t160 * t176 + t164 * t178) - g(2) * (-t159 * t176 + t162 * t178) - g(3) * (-t161 * t176 - t178 * t183);];
U_reg = t1;
