% Calculate minimal parameter regressor of potential energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:31:09
% EndTime: 2024-09-27 17:31:09
% DurationCPUTime: 0.02s
% Computational Cost: add. (100->31), mult. (58->49), div. (0->0), fcn. (60->16), ass. (0->26)
t173 = g(3) * sin(pkin(5));
t166 = cos(pkin(5));
t167 = sin(qJ(4));
t172 = t166 * t167;
t169 = cos(qJ(4));
t171 = t166 * t169;
t164 = qJ(1) + qJ(2);
t163 = qJ(4) + qJ(5);
t170 = cos(qJ(1));
t168 = sin(qJ(1));
t162 = qJ(3) + t164;
t161 = cos(t164);
t160 = cos(t163);
t159 = sin(t164);
t158 = sin(t163);
t157 = pkin(5) - t163;
t156 = pkin(5) + t163;
t155 = cos(t162);
t154 = sin(t162);
t153 = cos(t156);
t152 = sin(t157);
t151 = cos(t157) / 0.2e1;
t150 = sin(t156) / 0.2e1;
t149 = t151 + t153 / 0.2e1;
t148 = t150 - t152 / 0.2e1;
t1 = [0, -g(1) * t170 - g(2) * t168, g(1) * t168 - g(2) * t170, 0, -g(1) * t161 - g(2) * t159, g(1) * t159 - g(2) * t161, 0, -g(1) * t155 - g(2) * t154, g(1) * t154 - g(2) * t155, 0, 0, 0, 0, 0, -g(1) * (-t154 * t172 + t155 * t169) - g(2) * (t154 * t169 + t155 * t172) - t167 * t173, -g(1) * (-t154 * t171 - t155 * t167) - g(2) * (-t154 * t167 + t155 * t171) - t169 * t173, 0, 0, 0, 0, 0, -g(1) * (-t154 * t148 + t155 * t160) - g(2) * (t155 * t148 + t154 * t160) - g(3) * (t151 - t153 / 0.2e1), -g(1) * (-t154 * t149 - t155 * t158) - g(2) * (t155 * t149 - t154 * t158) - g(3) * (t150 + t152 / 0.2e1);];
U_reg = t1;
