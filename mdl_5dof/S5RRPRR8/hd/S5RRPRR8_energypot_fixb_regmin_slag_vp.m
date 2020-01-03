% Calculate minimal parameter regressor of potential energy for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:12
% EndTime: 2019-12-31 20:18:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->24), mult. (56->35), div. (0->0), fcn. (57->8), ass. (0->19)
t167 = qJ(2) + pkin(9) + qJ(4);
t164 = sin(t167);
t180 = g(3) * t164;
t169 = sin(qJ(5));
t171 = sin(qJ(1));
t179 = t171 * t169;
t172 = cos(qJ(5));
t178 = t171 * t172;
t174 = cos(qJ(1));
t177 = t174 * t169;
t176 = t174 * t172;
t175 = g(1) * t174 + g(2) * t171;
t173 = cos(qJ(2));
t170 = sin(qJ(2));
t168 = -pkin(6) - qJ(3);
t166 = t173 * pkin(2) + pkin(1);
t165 = cos(t167);
t163 = g(1) * t171 - g(2) * t174;
t1 = [0, -t175, t163, 0, 0, 0, 0, 0, -g(3) * t170 - t175 * t173, -g(3) * t173 + t175 * t170, -t163, -g(1) * (t174 * t166 - t171 * t168) - g(2) * (t171 * t166 + t174 * t168) - g(3) * (t170 * pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -t175 * t165 - t180, -g(3) * t165 + t175 * t164, 0, 0, 0, 0, 0, -g(1) * (t165 * t176 + t179) - g(2) * (t165 * t178 - t177) - t172 * t180, -g(1) * (-t165 * t177 + t178) - g(2) * (-t165 * t179 - t176) + t169 * t180;];
U_reg = t1;
