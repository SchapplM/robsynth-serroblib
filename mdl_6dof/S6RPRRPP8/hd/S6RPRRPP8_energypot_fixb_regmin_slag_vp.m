% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:52
% EndTime: 2019-03-09 04:55:52
% DurationCPUTime: 0.08s
% Computational Cost: add. (94->51), mult. (188->61), div. (0->0), fcn. (201->6), ass. (0->31)
t168 = sin(qJ(3));
t188 = pkin(3) * t168;
t172 = cos(qJ(1));
t187 = g(2) * t172;
t167 = sin(qJ(4));
t171 = cos(qJ(3));
t186 = t167 * t171;
t169 = sin(qJ(1));
t185 = t169 * t167;
t170 = cos(qJ(4));
t184 = t169 * t170;
t183 = t169 * t171;
t182 = t170 * t171;
t181 = t172 * t167;
t180 = t172 * t170;
t179 = t172 * pkin(1) + t169 * qJ(2);
t178 = -qJ(2) - t188;
t152 = g(1) * t169 - t187;
t177 = t171 * pkin(3) + pkin(4) * t182 + t168 * pkin(8) + qJ(5) * t186 + pkin(2) + pkin(6);
t150 = t168 * t181 + t184;
t151 = t168 * t180 - t185;
t163 = t169 * pkin(1);
t176 = t172 * t171 * pkin(8) - t151 * pkin(4) + t169 * pkin(7) - t150 * qJ(5) + t163;
t148 = t168 * t185 - t180;
t149 = t168 * t184 + t181;
t175 = t149 * pkin(4) + t172 * pkin(7) + t148 * qJ(5) + t169 * t188 + t179;
t174 = g(1) * t148 - g(2) * t150 + g(3) * t186;
t173 = g(1) * t149 - g(2) * t151 + g(3) * t182;
t153 = g(1) * t172 + g(2) * t169;
t145 = -g(3) * t168 + t152 * t171;
t1 = [0, -t153, t152, t153, -t152, -g(1) * t179 - g(2) * (-t172 * qJ(2) + t163) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t171 - t152 * t168, -t145, 0, 0, 0, 0, 0, -t173, t174, t145, t173, -t174, -g(1) * (-pkin(8) * t183 + t175) - g(2) * (t178 * t172 + t176) - g(3) * t177, t145, -t174, -t173, -g(1) * (t149 * qJ(6) + (-pkin(5) - pkin(8)) * t183 + t175) - g(2) * (-t151 * qJ(6) + t176) - g(3) * (t168 * pkin(5) + qJ(6) * t182 + t177) - (pkin(5) * t171 + t178) * t187;];
U_reg  = t1;
