% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP8
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
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:23
% EndTime: 2019-03-09 02:16:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (103->49), mult. (135->60), div. (0->0), fcn. (138->8), ass. (0->31)
t171 = pkin(9) + qJ(4);
t166 = sin(t171);
t172 = sin(pkin(9));
t191 = pkin(3) * t172 + pkin(4) * t166;
t190 = pkin(2) + pkin(6);
t176 = sin(qJ(1));
t187 = g(1) * t176;
t167 = cos(t171);
t186 = g(3) * t167;
t175 = sin(qJ(5));
t185 = t176 * t175;
t177 = cos(qJ(5));
t184 = t176 * t177;
t178 = cos(qJ(1));
t183 = t178 * t175;
t182 = t178 * t177;
t181 = t178 * pkin(1) + t176 * qJ(2);
t169 = t176 * pkin(1);
t180 = -t178 * qJ(2) + t169;
t163 = -g(2) * t178 + t187;
t159 = t166 * t185 - t182;
t161 = t166 * t183 + t184;
t179 = g(1) * t159 - g(2) * t161 + t175 * t186;
t174 = -pkin(7) - qJ(3);
t173 = cos(pkin(9));
t164 = g(1) * t178 + g(2) * t176;
t162 = -t166 * t182 + t185;
t160 = t166 * t184 + t183;
t158 = -g(3) * t166 + t163 * t167;
t157 = -g(1) * t160 - g(2) * t162 - t177 * t186;
t1 = [0, -t164, t163, t164, -t163, -g(3) * pkin(6) - g(1) * t181 - g(2) * t180, -g(3) * t173 - t163 * t172, g(3) * t172 - t163 * t173, -t164, -g(1) * (t178 * qJ(3) + t181) - g(2) * (t176 * qJ(3) + t180) - g(3) * t190, 0, 0, 0, 0, 0, -t163 * t166 - t186, -t158, 0, 0, 0, 0, 0, t157, t179, t157, t158, -t179, -g(1) * (t160 * pkin(5) + t159 * qJ(6) + t191 * t176 + t181) - g(2) * (t162 * pkin(5) - t161 * qJ(6) - t176 * t174 + t169) - g(3) * (t173 * pkin(3) + t166 * pkin(8) + t190) + (pkin(8) * t187 - g(3) * (pkin(5) * t177 + qJ(6) * t175 + pkin(4))) * t167 + (g(1) * t174 - g(2) * (pkin(8) * t167 - qJ(2) - t191)) * t178;];
U_reg  = t1;
