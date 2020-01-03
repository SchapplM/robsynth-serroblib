% Calculate minimal parameter regressor of potential energy for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:59
% EndTime: 2019-12-31 19:07:59
% DurationCPUTime: 0.05s
% Computational Cost: add. (57->24), mult. (63->37), div. (0->0), fcn. (64->10), ass. (0->20)
t163 = pkin(9) + qJ(3);
t162 = qJ(4) + t163;
t158 = sin(t162);
t175 = g(3) * t158;
t166 = sin(qJ(5));
t167 = sin(qJ(1));
t174 = t167 * t166;
t168 = cos(qJ(5));
t173 = t167 * t168;
t169 = cos(qJ(1));
t172 = t169 * t166;
t171 = t169 * t168;
t170 = g(1) * t169 + g(2) * t167;
t165 = cos(pkin(9));
t164 = sin(pkin(9));
t161 = cos(t163);
t160 = sin(t163);
t159 = cos(t162);
t157 = g(1) * t167 - g(2) * t169;
t1 = [0, -t170, t157, -g(3) * t164 - t170 * t165, -g(3) * t165 + t170 * t164, -t157, -g(1) * (t169 * pkin(1) + t167 * qJ(2)) - g(2) * (t167 * pkin(1) - t169 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t160 - t170 * t161, -g(3) * t161 + t170 * t160, 0, 0, 0, 0, 0, -t170 * t159 - t175, -g(3) * t159 + t170 * t158, 0, 0, 0, 0, 0, -g(1) * (t159 * t171 + t174) - g(2) * (t159 * t173 - t172) - t168 * t175, -g(1) * (-t159 * t172 + t173) - g(2) * (-t159 * t174 - t171) + t166 * t175;];
U_reg = t1;
