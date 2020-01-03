% Calculate minimal parameter regressor of potential energy for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (59->30), mult. (73->47), div. (0->0), fcn. (78->10), ass. (0->24)
t168 = pkin(9) + qJ(3);
t164 = sin(t168);
t185 = g(3) * t164;
t169 = qJ(4) + qJ(5);
t166 = sin(t169);
t173 = sin(qJ(1));
t184 = t173 * t166;
t167 = cos(t169);
t183 = t173 * t167;
t172 = sin(qJ(4));
t182 = t173 * t172;
t174 = cos(qJ(4));
t181 = t173 * t174;
t175 = cos(qJ(1));
t180 = t175 * t166;
t179 = t175 * t167;
t178 = t175 * t172;
t177 = t175 * t174;
t176 = g(1) * t175 + g(2) * t173;
t171 = cos(pkin(9));
t170 = sin(pkin(9));
t165 = cos(t168);
t163 = g(1) * t173 - g(2) * t175;
t1 = [0, -t176, t163, -g(3) * t170 - t176 * t171, -g(3) * t171 + t176 * t170, -t163, -g(1) * (t175 * pkin(1) + t173 * qJ(2)) - g(2) * (t173 * pkin(1) - t175 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t176 * t165 - t185, -g(3) * t165 + t176 * t164, 0, 0, 0, 0, 0, -g(1) * (t165 * t177 + t182) - g(2) * (t165 * t181 - t178) - t174 * t185, -g(1) * (-t165 * t178 + t181) - g(2) * (-t165 * t182 - t177) + t172 * t185, 0, 0, 0, 0, 0, -g(1) * (t165 * t179 + t184) - g(2) * (t165 * t183 - t180) - t167 * t185, -g(1) * (-t165 * t180 + t183) - g(2) * (-t165 * t184 - t179) + t166 * t185;];
U_reg = t1;
