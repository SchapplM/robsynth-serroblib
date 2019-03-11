% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:48
% EndTime: 2019-03-09 02:23:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (72->32), mult. (72->54), div. (0->0), fcn. (74->10), ass. (0->21)
t182 = g(3) * (qJ(2) + pkin(6));
t173 = cos(qJ(4));
t181 = g(3) * t173;
t167 = qJ(5) + qJ(6);
t164 = sin(t167);
t170 = sin(qJ(4));
t180 = t164 * t170;
t165 = cos(t167);
t179 = t165 * t170;
t169 = sin(qJ(5));
t178 = t169 * t170;
t172 = cos(qJ(5));
t177 = t170 * t172;
t166 = qJ(1) + pkin(10);
t162 = sin(t166);
t163 = cos(t166);
t176 = -g(1) * t162 + g(2) * t163;
t171 = sin(qJ(1));
t174 = cos(qJ(1));
t175 = -g(1) * t174 - g(2) * t171;
t1 = [0, t175, g(1) * t171 - g(2) * t174, t175 * pkin(1) - t182, g(1) * t163 + g(2) * t162, t176, -g(1) * (t174 * pkin(1) + t163 * pkin(2) + t162 * qJ(3)) - g(2) * (t171 * pkin(1) + t162 * pkin(2) - t163 * qJ(3)) - t182, 0, 0, 0, 0, 0, t176 * t170 - t181, g(3) * t170 + t176 * t173, 0, 0, 0, 0, 0, -g(1) * (t162 * t177 + t163 * t169) - g(2) * (t162 * t169 - t163 * t177) - t172 * t181, -g(1) * (-t162 * t178 + t163 * t172) - g(2) * (t162 * t172 + t163 * t178) + t169 * t181, 0, 0, 0, 0, 0, -g(1) * (t162 * t179 + t163 * t164) - g(2) * (t162 * t164 - t163 * t179) - t165 * t181, -g(1) * (-t162 * t180 + t163 * t165) - g(2) * (t162 * t165 + t163 * t180) + t164 * t181;];
U_reg  = t1;
