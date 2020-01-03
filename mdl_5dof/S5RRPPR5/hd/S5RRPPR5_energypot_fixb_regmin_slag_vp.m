% Calculate minimal parameter regressor of potential energy for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:57
% EndTime: 2019-12-31 19:29:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (70->28), mult. (84->38), div. (0->0), fcn. (84->8), ass. (0->20)
t171 = sin(qJ(1));
t174 = cos(qJ(1));
t178 = g(1) * t174 + g(2) * t171;
t170 = sin(qJ(2));
t181 = t170 * pkin(2) + pkin(5);
t173 = cos(qJ(2));
t162 = t173 * pkin(2) + pkin(1);
t168 = -pkin(6) - qJ(3);
t180 = t171 * t162 + t174 * t168;
t179 = t174 * t162 - t171 * t168;
t167 = qJ(2) + pkin(8);
t163 = sin(t167);
t164 = cos(t167);
t177 = pkin(3) * t164 + qJ(4) * t163;
t169 = sin(qJ(5));
t172 = cos(qJ(5));
t176 = t163 * t172 - t164 * t169;
t175 = t163 * t169 + t164 * t172;
t158 = g(1) * t171 - g(2) * t174;
t1 = [0, -t178, t158, 0, 0, 0, 0, 0, -g(3) * t170 - t178 * t173, -g(3) * t173 + t178 * t170, -t158, -g(1) * t179 - g(2) * t180 - g(3) * t181, -g(3) * t163 - t178 * t164, -t158, g(3) * t164 - t178 * t163, -g(1) * (t177 * t174 + t179) - g(2) * (t177 * t171 + t180) - g(3) * (t163 * pkin(3) - t164 * qJ(4) + t181), 0, 0, 0, 0, 0, -g(3) * t176 - t178 * t175, g(3) * t175 - t178 * t176;];
U_reg = t1;
