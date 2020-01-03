% Calculate minimal parameter regressor of potential energy for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:26
% EndTime: 2019-12-31 22:22:26
% DurationCPUTime: 0.04s
% Computational Cost: add. (52->19), mult. (54->30), div. (0->0), fcn. (58->10), ass. (0->19)
t165 = qJ(2) + qJ(3);
t164 = qJ(4) + t165;
t160 = sin(t164);
t177 = g(3) * t160;
t166 = sin(qJ(5));
t168 = sin(qJ(1));
t176 = t168 * t166;
t169 = cos(qJ(5));
t175 = t168 * t169;
t171 = cos(qJ(1));
t174 = t171 * t166;
t173 = t171 * t169;
t172 = g(1) * t171 + g(2) * t168;
t170 = cos(qJ(2));
t167 = sin(qJ(2));
t163 = cos(t165);
t162 = sin(t165);
t161 = cos(t164);
t1 = [0, -t172, g(1) * t168 - g(2) * t171, 0, 0, 0, 0, 0, -g(3) * t167 - t172 * t170, -g(3) * t170 + t172 * t167, 0, 0, 0, 0, 0, -g(3) * t162 - t172 * t163, -g(3) * t163 + t172 * t162, 0, 0, 0, 0, 0, -t172 * t161 - t177, -g(3) * t161 + t172 * t160, 0, 0, 0, 0, 0, -g(1) * (t161 * t173 + t176) - g(2) * (t161 * t175 - t174) - t169 * t177, -g(1) * (-t161 * t174 + t175) - g(2) * (-t161 * t176 - t173) + t166 * t177;];
U_reg = t1;
