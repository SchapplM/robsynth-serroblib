% Calculate minimal parameter regressor of potential energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:03
% EndTime: 2019-12-31 20:14:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->39), mult. (132->48), div. (0->0), fcn. (138->6), ass. (0->26)
t159 = sin(qJ(2));
t162 = cos(qJ(2));
t176 = pkin(2) * t162 + qJ(3) * t159 + pkin(1);
t174 = g(3) * t162;
t173 = t159 * pkin(2) + pkin(5);
t158 = sin(qJ(4));
t160 = sin(qJ(1));
t171 = t160 * t158;
t161 = cos(qJ(4));
t170 = t160 * t161;
t163 = cos(qJ(1));
t169 = t163 * t158;
t168 = t163 * t161;
t167 = t160 * pkin(6) + t176 * t163;
t166 = g(1) * t163 + g(2) * t160;
t165 = -t163 * pkin(6) + t176 * t160;
t144 = -t159 * t168 + t171;
t146 = t159 * t170 + t169;
t164 = g(1) * t144 - g(2) * t146 + t161 * t174;
t148 = g(1) * t160 - g(2) * t163;
t147 = t159 * t171 - t168;
t145 = t159 * t169 + t170;
t143 = g(3) * t159 + t166 * t162;
t142 = t166 * t159 - t174;
t141 = -g(1) * t145 - g(2) * t147 + t158 * t174;
t1 = [0, -t166, t148, 0, 0, 0, 0, 0, -t143, t142, -t148, t143, -t142, -g(1) * t167 - g(2) * t165 - g(3) * (-t162 * qJ(3) + t173), 0, 0, 0, 0, 0, t141, t164, t141, -t143, -t164, -g(1) * (t160 * pkin(3) + t145 * pkin(4) + t144 * qJ(5) + t167) - g(2) * (-t163 * pkin(3) + t147 * pkin(4) - t146 * qJ(5) + t165) - g(3) * (t159 * pkin(7) + t173) + (-g(3) * (-pkin(4) * t158 + qJ(5) * t161 - qJ(3)) - t166 * pkin(7)) * t162;];
U_reg = t1;
