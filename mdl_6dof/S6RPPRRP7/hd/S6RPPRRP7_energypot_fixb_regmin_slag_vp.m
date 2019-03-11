% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:53
% EndTime: 2019-03-09 02:13:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->44), mult. (102->54), div. (0->0), fcn. (97->8), ass. (0->27)
t179 = pkin(2) + pkin(6);
t160 = pkin(9) + qJ(4);
t156 = cos(t160);
t178 = g(3) * t156;
t165 = sin(qJ(5));
t166 = sin(qJ(1));
t177 = t166 * t165;
t167 = cos(qJ(5));
t176 = t166 * t167;
t168 = cos(qJ(1));
t175 = t168 * t165;
t174 = t168 * t167;
t173 = t168 * pkin(1) + t166 * qJ(2);
t172 = g(1) * t173;
t171 = pkin(5) * t165 + pkin(7) + qJ(3);
t158 = t166 * pkin(1);
t170 = -t168 * qJ(2) + t158;
t152 = g(1) * t166 - g(2) * t168;
t154 = t167 * pkin(5) + pkin(4);
t155 = sin(t160);
t161 = sin(pkin(9));
t163 = -qJ(6) - pkin(8);
t169 = pkin(3) * t161 + t154 * t155 + t156 * t163;
t162 = cos(pkin(9));
t153 = g(1) * t168 + g(2) * t166;
t151 = -g(3) * t155 + t152 * t156;
t1 = [0, -t153, t152, t153, -t152, -g(3) * pkin(6) - g(2) * t170 - t172, -g(3) * t162 - t152 * t161, g(3) * t161 - t152 * t162, -t153, -g(1) * (t168 * qJ(3) + t173) - g(2) * (t166 * qJ(3) + t170) - g(3) * t179, 0, 0, 0, 0, 0, -t152 * t155 - t178, -t151, 0, 0, 0, 0, 0, -g(1) * (t155 * t176 + t175) - g(2) * (-t155 * t174 + t177) - t167 * t178, -g(1) * (-t155 * t177 + t174) - g(2) * (t155 * t175 + t176) + t165 * t178, t151, -t172 - g(2) * t158 - g(3) * (t162 * pkin(3) + t156 * t154 - t155 * t163 + t179) + (-g(1) * t169 - g(2) * t171) * t166 + (-g(1) * t171 - g(2) * (-qJ(2) - t169)) * t168;];
U_reg  = t1;
