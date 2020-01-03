% Calculate minimal parameter regressor of potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:23
% EndTime: 2020-01-03 11:50:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (70->41), mult. (97->59), div. (0->0), fcn. (99->8), ass. (0->27)
t157 = sin(pkin(8));
t174 = g(1) * t157;
t156 = qJ(3) + qJ(4);
t152 = sin(t156);
t160 = sin(qJ(1));
t173 = t160 * t152;
t153 = cos(t156);
t172 = t160 * t153;
t159 = sin(qJ(3));
t171 = t160 * t159;
t161 = cos(qJ(3));
t170 = t160 * t161;
t162 = cos(qJ(1));
t169 = t162 * t152;
t168 = t162 * t153;
t167 = t162 * t159;
t166 = t162 * t161;
t165 = t159 * pkin(3) + pkin(4) * t152 + qJ(2);
t164 = -g(2) * t160 + g(3) * t162;
t149 = t161 * pkin(3) + pkin(4) * t153 + pkin(2);
t155 = -qJ(5) - pkin(7) - pkin(6);
t158 = cos(pkin(8));
t163 = t149 * t158 - t155 * t157;
t154 = t160 * pkin(1);
t151 = g(2) * t162 + g(3) * t160;
t148 = g(1) * t158 + t164 * t157;
t1 = [0, t164, -t151, t164 * t158 - t174, -t148, t151, -g(1) * pkin(5) - g(2) * (-t162 * qJ(2) + t154) - g(3) * (-t162 * pkin(1) - t160 * qJ(2)), 0, 0, 0, 0, 0, -t161 * t174 - g(2) * (t158 * t170 - t167) - g(3) * (-t158 * t166 - t171), t159 * t174 - g(2) * (-t158 * t171 - t166) - g(3) * (t158 * t167 - t170), 0, 0, 0, 0, 0, -t153 * t174 - g(2) * (t158 * t172 - t169) - g(3) * (-t158 * t168 - t173), t152 * t174 - g(2) * (-t158 * t173 - t168) - g(3) * (t158 * t169 - t172), t148, -g(1) * (t157 * t149 + t158 * t155 + pkin(5)) - g(2) * t154 + (-g(2) * t163 + g(3) * t165) * t160 + (g(2) * t165 - g(3) * (-pkin(1) - t163)) * t162;];
U_reg = t1;
