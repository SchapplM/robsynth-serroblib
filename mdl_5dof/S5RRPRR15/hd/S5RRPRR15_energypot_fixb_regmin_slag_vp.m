% Calculate minimal parameter regressor of potential energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR15_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:08
% EndTime: 2019-12-31 20:43:08
% DurationCPUTime: 0.06s
% Computational Cost: add. (47->33), mult. (83->47), div. (0->0), fcn. (88->8), ass. (0->24)
t173 = cos(qJ(2));
t185 = g(3) * t173;
t168 = qJ(4) + qJ(5);
t166 = sin(t168);
t171 = sin(qJ(1));
t184 = t171 * t166;
t167 = cos(t168);
t183 = t171 * t167;
t169 = sin(qJ(4));
t182 = t171 * t169;
t172 = cos(qJ(4));
t181 = t171 * t172;
t174 = cos(qJ(1));
t180 = t174 * t166;
t179 = t174 * t167;
t178 = t174 * t169;
t177 = t174 * t172;
t176 = g(1) * t174 + g(2) * t171;
t170 = sin(qJ(2));
t175 = pkin(2) * t173 + qJ(3) * t170 + pkin(1);
t165 = g(1) * t171 - g(2) * t174;
t164 = g(3) * t170 + t176 * t173;
t163 = t176 * t170 - t185;
t1 = [0, -t176, t165, 0, 0, 0, 0, 0, -t164, t163, -t165, t164, -t163, -g(3) * (t170 * pkin(2) - t173 * qJ(3) + pkin(5)) + (g(2) * pkin(6) - g(1) * t175) * t174 + (-g(1) * pkin(6) - g(2) * t175) * t171, 0, 0, 0, 0, 0, -g(1) * (t170 * t178 + t181) - g(2) * (t170 * t182 - t177) + t169 * t185, -g(1) * (t170 * t177 - t182) - g(2) * (t170 * t181 + t178) + t172 * t185, 0, 0, 0, 0, 0, -g(1) * (t170 * t180 + t183) - g(2) * (t170 * t184 - t179) + t166 * t185, -g(1) * (t170 * t179 - t184) - g(2) * (t170 * t183 + t180) + t167 * t185;];
U_reg = t1;
