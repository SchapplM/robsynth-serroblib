% Calculate minimal parameter regressor of potential energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:50
% EndTime: 2019-12-31 17:47:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (65->49), mult. (137->69), div. (0->0), fcn. (142->8), ass. (0->31)
t165 = sin(qJ(1));
t163 = cos(pkin(7));
t179 = t163 * t165;
t161 = sin(pkin(7));
t181 = qJ(3) * t161;
t183 = pkin(2) * t179 + t165 * t181;
t182 = g(3) * t163;
t160 = sin(pkin(8));
t180 = t160 * t163;
t166 = cos(qJ(5));
t178 = t163 * t166;
t167 = cos(qJ(1));
t177 = t163 * t167;
t176 = t165 * t160;
t162 = cos(pkin(8));
t175 = t165 * t162;
t174 = t167 * t160;
t173 = t167 * t162;
t172 = t167 * pkin(1) + t165 * qJ(2);
t157 = t165 * pkin(1);
t171 = -t167 * qJ(2) + t157;
t170 = pkin(2) * t177 + t167 * t181 + t172;
t169 = t161 * pkin(2) - t163 * qJ(3) + pkin(5);
t168 = g(1) * t167 + g(2) * t165;
t164 = sin(qJ(5));
t150 = g(1) * t165 - g(2) * t167;
t149 = t161 * t176 - t173;
t148 = t161 * t174 + t175;
t147 = g(3) * t161 + t168 * t163;
t146 = t168 * t161 - t182;
t1 = [0, -t168, t150, -t147, t146, -t150, -g(3) * pkin(5) - g(1) * t172 - g(2) * t171, -t150, t147, -t146, -g(1) * t170 - g(2) * (t171 + t183) - g(3) * t169, -g(1) * t148 - g(2) * t149 + g(3) * t180, -g(1) * (t161 * t173 - t176) - g(2) * (t161 * t175 + t174) + t162 * t182, -t147, -g(1) * (t165 * pkin(3) + qJ(4) * t177 + t170) - g(2) * (qJ(4) * t179 + t157 + (-pkin(3) - qJ(2)) * t167 + t183) - g(3) * (t161 * qJ(4) + t169), 0, 0, 0, 0, 0, -g(1) * (t148 * t166 + t164 * t177) - g(2) * (t149 * t166 + t164 * t179) - g(3) * (-t160 * t178 + t161 * t164), -g(1) * (-t148 * t164 + t166 * t177) - g(2) * (-t149 * t164 + t165 * t178) - g(3) * (t161 * t166 + t164 * t180);];
U_reg = t1;
