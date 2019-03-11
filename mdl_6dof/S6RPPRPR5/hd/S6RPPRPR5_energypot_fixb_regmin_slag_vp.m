% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:20
% EndTime: 2019-03-09 01:49:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (65->47), mult. (104->60), div. (0->0), fcn. (103->8), ass. (0->29)
t159 = sin(qJ(4));
t161 = cos(qJ(4));
t178 = pkin(4) * t159 - qJ(5) * t161;
t177 = pkin(2) + pkin(6);
t175 = g(3) * t161;
t156 = pkin(9) + qJ(6);
t148 = sin(t156);
t160 = sin(qJ(1));
t173 = t160 * t148;
t149 = cos(t156);
t172 = t160 * t149;
t157 = sin(pkin(9));
t171 = t160 * t157;
t158 = cos(pkin(9));
t170 = t160 * t158;
t162 = cos(qJ(1));
t169 = t162 * t148;
t168 = t162 * t149;
t167 = t162 * t157;
t166 = t162 * t158;
t165 = t162 * pkin(1) + t160 * qJ(2);
t164 = t162 * qJ(3) + t165;
t153 = t160 * pkin(1);
t163 = -t162 * qJ(2) + t153;
t147 = g(1) * t162 + g(2) * t160;
t150 = t160 * qJ(3);
t146 = g(1) * t160 - g(2) * t162;
t145 = -g(3) * t159 + t147 * t161;
t1 = [0, -t147, t146, t147, -t146, -g(3) * pkin(6) - g(1) * t165 - g(2) * t163, -t146, -t147, -g(1) * t164 - g(2) * (t150 + t163) - g(3) * t177, 0, 0, 0, 0, 0, -t147 * t159 - t175, -t145, -g(1) * (t159 * t166 - t171) - g(2) * (t159 * t170 + t167) - t158 * t175, -g(1) * (-t159 * t167 - t170) - g(2) * (-t159 * t171 + t166) + t157 * t175, t145, -g(1) * (-t160 * pkin(7) + t164) - g(2) * (t178 * t160 + t150 + t153) - g(3) * (t161 * pkin(4) + t159 * qJ(5) + pkin(3) + t177) + (-g(1) * t178 - g(2) * (pkin(7) - qJ(2))) * t162, 0, 0, 0, 0, 0, -g(1) * (t159 * t168 - t173) - g(2) * (t159 * t172 + t169) - t149 * t175, -g(1) * (-t159 * t169 - t172) - g(2) * (-t159 * t173 + t168) + t148 * t175;];
U_reg  = t1;
