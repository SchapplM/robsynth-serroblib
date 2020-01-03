% Calculate minimal parameter regressor of potential energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:49
% EndTime: 2019-12-31 18:14:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (60->34), mult. (77->39), div. (0->0), fcn. (68->6), ass. (0->19)
t121 = qJ(3) + pkin(7);
t115 = sin(t121);
t116 = cos(t121);
t135 = pkin(4) * t115 - qJ(5) * t116;
t123 = sin(qJ(3));
t134 = pkin(3) * t123;
t124 = sin(qJ(1));
t126 = cos(qJ(1));
t131 = t126 * pkin(1) + t124 * qJ(2);
t125 = cos(qJ(3));
t130 = t125 * pkin(3) + pkin(2) + pkin(5);
t129 = t124 * t134 + t131;
t128 = -qJ(2) - t134;
t118 = t124 * pkin(1);
t122 = -qJ(4) - pkin(6);
t127 = -t124 * t122 + t118;
t112 = g(1) * t124 - g(2) * t126;
t113 = g(1) * t126 + g(2) * t124;
t1 = [0, -t113, t112, t113, -t112, -g(1) * t131 - g(2) * (-t126 * qJ(2) + t118) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t125 - t112 * t123, g(3) * t123 - t112 * t125, -t113, -g(1) * (-t126 * t122 + t129) - g(2) * (t128 * t126 + t127) - g(3) * t130, -g(3) * t116 - t112 * t115, -t113, -g(3) * t115 + t112 * t116, -g(1) * (t135 * t124 + t129) - g(2) * t127 - g(3) * (t116 * pkin(4) + t115 * qJ(5) + t130) + (g(1) * t122 - g(2) * (t128 - t135)) * t126;];
U_reg = t1;
