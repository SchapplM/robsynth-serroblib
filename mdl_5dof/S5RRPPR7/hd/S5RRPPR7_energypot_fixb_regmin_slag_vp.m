% Calculate minimal parameter regressor of potential energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:21
% EndTime: 2019-12-31 19:36:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (62->32), mult. (78->44), div. (0->0), fcn. (76->8), ass. (0->23)
t162 = qJ(2) + pkin(8);
t159 = cos(t162);
t179 = g(3) * t159;
t165 = sin(qJ(2));
t178 = t165 * pkin(2) + pkin(5);
t164 = sin(qJ(5));
t166 = sin(qJ(1));
t177 = t166 * t164;
t167 = cos(qJ(5));
t176 = t166 * t167;
t169 = cos(qJ(1));
t175 = t169 * t164;
t174 = t169 * t167;
t168 = cos(qJ(2));
t157 = t168 * pkin(2) + pkin(1);
t163 = -pkin(6) - qJ(3);
t173 = t166 * t157 + t169 * t163;
t172 = t169 * t157 - t166 * t163;
t171 = g(1) * t169 + g(2) * t166;
t158 = sin(t162);
t170 = pkin(3) * t159 + qJ(4) * t158;
t153 = g(1) * t166 - g(2) * t169;
t1 = [0, -t171, t153, 0, 0, 0, 0, 0, -g(3) * t165 - t171 * t168, -g(3) * t168 + t171 * t165, -t153, -g(1) * t172 - g(2) * t173 - g(3) * t178, -t153, g(3) * t158 + t171 * t159, -t171 * t158 + t179, -g(1) * (t170 * t169 + t172) - g(2) * (t170 * t166 + t173) - g(3) * (t158 * pkin(3) - t159 * qJ(4) + t178), 0, 0, 0, 0, 0, -g(1) * (t158 * t175 + t176) - g(2) * (t158 * t177 - t174) + t164 * t179, -g(1) * (t158 * t174 - t177) - g(2) * (t158 * t176 + t175) + t167 * t179;];
U_reg = t1;
