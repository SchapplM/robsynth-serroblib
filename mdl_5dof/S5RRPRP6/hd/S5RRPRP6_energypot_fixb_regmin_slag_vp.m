% Calculate minimal parameter regressor of potential energy for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:27
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->35), mult. (78->48), div. (0->0), fcn. (76->8), ass. (0->25)
t167 = cos(qJ(4));
t156 = t167 * pkin(4) + pkin(3);
t161 = qJ(2) + pkin(8);
t158 = sin(t161);
t159 = cos(t161);
t162 = -qJ(5) - pkin(7);
t180 = t156 * t159 - t158 * t162;
t179 = g(3) * t158;
t165 = sin(qJ(2));
t178 = t165 * pkin(2) + pkin(5);
t164 = sin(qJ(4));
t166 = sin(qJ(1));
t175 = t166 * t164;
t174 = t166 * t167;
t169 = cos(qJ(1));
t173 = t169 * t164;
t172 = t169 * t167;
t168 = cos(qJ(2));
t157 = t168 * pkin(2) + pkin(1);
t163 = -pkin(6) - qJ(3);
t171 = t166 * t157 + t169 * t163;
t170 = g(1) * t169 + g(2) * t166;
t154 = t169 * t157;
t152 = g(1) * t166 - g(2) * t169;
t1 = [0, -t170, t152, 0, 0, 0, 0, 0, -g(3) * t165 - t170 * t168, -g(3) * t168 + t170 * t165, -t152, -g(1) * (-t166 * t163 + t154) - g(2) * t171 - g(3) * t178, 0, 0, 0, 0, 0, -g(1) * (t159 * t172 + t175) - g(2) * (t159 * t174 - t173) - t167 * t179, -g(1) * (-t159 * t173 + t174) - g(2) * (-t159 * t175 - t172) + t164 * t179, g(3) * t159 - t170 * t158, -g(1) * (t180 * t169 + t154) - g(2) * (-pkin(4) * t173 + t171) - g(3) * (t158 * t156 + t159 * t162 + t178) + (-g(1) * (pkin(4) * t164 - t163) - g(2) * t180) * t166;];
U_reg = t1;
