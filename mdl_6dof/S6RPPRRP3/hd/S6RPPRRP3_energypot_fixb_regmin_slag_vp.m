% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:47
% EndTime: 2019-03-09 02:03:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (110->44), mult. (118->58), div. (0->0), fcn. (121->8), ass. (0->28)
t163 = sin(qJ(4));
t178 = pkin(4) * t163;
t160 = qJ(1) + pkin(9);
t156 = sin(t160);
t177 = g(1) * t156;
t161 = qJ(2) + pkin(6);
t176 = g(3) * t161;
t166 = cos(qJ(4));
t175 = g(3) * t166;
t162 = sin(qJ(5));
t174 = t162 * t163;
t165 = cos(qJ(5));
t173 = t163 * t165;
t164 = sin(qJ(1));
t172 = t164 * pkin(1) + t156 * pkin(2);
t157 = cos(t160);
t167 = cos(qJ(1));
t171 = t167 * pkin(1) + t157 * pkin(2) + t156 * qJ(3);
t170 = -g(2) * t157 + t177;
t169 = -g(1) * t167 - g(2) * t164;
t148 = t156 * t174 - t157 * t165;
t150 = t156 * t165 + t157 * t174;
t168 = g(1) * t148 - g(2) * t150 + t162 * t175;
t151 = t156 * t162 - t157 * t173;
t149 = t156 * t173 + t157 * t162;
t147 = -g(3) * t163 + t170 * t166;
t146 = -g(1) * t149 - g(2) * t151 - t165 * t175;
t1 = [0, t169, g(1) * t164 - g(2) * t167, t169 * pkin(1) - t176, g(1) * t157 + g(2) * t156, -t170, -g(1) * t171 - g(2) * (-t157 * qJ(3) + t172) - t176, 0, 0, 0, 0, 0, -t170 * t163 - t175, -t147, 0, 0, 0, 0, 0, t146, t168, t146, t147, -t168, -g(1) * (t149 * pkin(5) + t148 * qJ(6) + t156 * t178 + t171) - g(2) * (t151 * pkin(5) + t156 * pkin(7) - t150 * qJ(6) + t172) - g(3) * (t163 * pkin(8) + pkin(3) + t161) + (pkin(8) * t177 - g(3) * (pkin(5) * t165 + qJ(6) * t162 + pkin(4))) * t166 + (-g(1) * pkin(7) - g(2) * (pkin(8) * t166 - qJ(3) - t178)) * t157;];
U_reg  = t1;
