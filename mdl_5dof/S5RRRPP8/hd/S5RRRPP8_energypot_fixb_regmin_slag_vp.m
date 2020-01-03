% Calculate minimal parameter regressor of potential energy for
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:23
% EndTime: 2019-12-31 21:09:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (82->44), mult. (173->55), div. (0->0), fcn. (189->6), ass. (0->25)
t160 = sin(qJ(3));
t161 = sin(qJ(2));
t178 = t160 * t161;
t162 = sin(qJ(1));
t177 = t161 * t162;
t163 = cos(qJ(3));
t176 = t161 * t163;
t165 = cos(qJ(1));
t175 = t161 * t165;
t164 = cos(qJ(2));
t174 = t162 * t164;
t173 = t165 * t160;
t172 = t165 * t163;
t171 = t161 * pkin(2) + pkin(3) * t176 + qJ(4) * t178 + pkin(5);
t170 = g(1) * t165 + g(2) * t162;
t145 = -t162 * t163 + t164 * t173;
t146 = t162 * t160 + t164 * t172;
t169 = t146 * pkin(3) + t162 * pkin(6) + pkin(7) * t175 + t145 * qJ(4) + (pkin(2) * t164 + pkin(1)) * t165;
t143 = t160 * t174 + t172;
t168 = g(1) * t145 + g(2) * t143 + g(3) * t178;
t144 = t163 * t174 - t173;
t167 = g(1) * t146 + g(2) * t144 + g(3) * t176;
t166 = t162 * pkin(1) + pkin(2) * t174 + t144 * pkin(3) - t165 * pkin(6) + pkin(7) * t177 + t143 * qJ(4);
t140 = -g(3) * t164 + t170 * t161;
t1 = [0, -t170, g(1) * t162 - g(2) * t165, 0, 0, 0, 0, 0, -g(3) * t161 - t170 * t164, t140, 0, 0, 0, 0, 0, -t167, t168, -t140, t167, -t168, -g(1) * t169 - g(2) * t166 - g(3) * (-t164 * pkin(7) + t171), -t140, -t168, -t167, -g(1) * (pkin(4) * t175 + t146 * qJ(5) + t169) - g(2) * (pkin(4) * t177 + t144 * qJ(5) + t166) - g(3) * (qJ(5) * t176 + (-pkin(4) - pkin(7)) * t164 + t171);];
U_reg = t1;
