% Calculate minimal parameter regressor of potential energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:53
% EndTime: 2019-12-31 21:05:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (82->41), mult. (173->55), div. (0->0), fcn. (189->6), ass. (0->25)
t157 = sin(qJ(3));
t158 = sin(qJ(2));
t174 = t157 * t158;
t159 = sin(qJ(1));
t173 = t158 * t159;
t160 = cos(qJ(3));
t172 = t158 * t160;
t162 = cos(qJ(1));
t171 = t158 * t162;
t161 = cos(qJ(2));
t170 = t159 * t161;
t169 = t162 * t157;
t168 = t162 * t160;
t167 = t158 * pkin(2) + pkin(3) * t172 + qJ(4) * t174 + pkin(5);
t166 = g(1) * t162 + g(2) * t159;
t144 = -t159 * t160 + t161 * t169;
t145 = t159 * t157 + t161 * t168;
t165 = t145 * pkin(3) + t159 * pkin(6) + pkin(7) * t171 + t144 * qJ(4) + (pkin(2) * t161 + pkin(1)) * t162;
t142 = t157 * t170 + t168;
t164 = g(1) * t144 + g(2) * t142 + g(3) * t174;
t143 = t160 * t170 - t169;
t163 = t159 * pkin(1) + pkin(2) * t170 + t143 * pkin(3) - t162 * pkin(6) + pkin(7) * t173 + t142 * qJ(4);
t139 = -g(3) * t161 + t166 * t158;
t138 = -g(1) * t145 - g(2) * t143 - g(3) * t172;
t1 = [0, -t166, g(1) * t159 - g(2) * t162, 0, 0, 0, 0, 0, -g(3) * t158 - t166 * t161, t139, 0, 0, 0, 0, 0, t138, t164, t138, -t139, -t164, -g(1) * t165 - g(2) * t163 - g(3) * (-t161 * pkin(7) + t167), t138, -t164, t139, -g(1) * (t145 * pkin(4) - qJ(5) * t171 + t165) - g(2) * (t143 * pkin(4) - qJ(5) * t173 + t163) - g(3) * (pkin(4) * t172 + (-pkin(7) + qJ(5)) * t161 + t167);];
U_reg = t1;
