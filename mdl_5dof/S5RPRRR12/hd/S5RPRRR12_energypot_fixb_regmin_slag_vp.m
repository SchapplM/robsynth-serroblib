% Calculate minimal parameter regressor of potential energy for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:08
% DurationCPUTime: 0.04s
% Computational Cost: add. (36->21), mult. (55->33), div. (0->0), fcn. (56->8), ass. (0->17)
t128 = qJ(3) + qJ(4);
t127 = cos(t128);
t139 = g(3) * t127;
t129 = sin(qJ(5));
t131 = sin(qJ(1));
t138 = t131 * t129;
t132 = cos(qJ(5));
t137 = t131 * t132;
t134 = cos(qJ(1));
t136 = t134 * t129;
t135 = t134 * t132;
t124 = g(1) * t131 - g(2) * t134;
t133 = cos(qJ(3));
t130 = sin(qJ(3));
t126 = sin(t128);
t125 = g(1) * t134 + g(2) * t131;
t1 = [0, -t125, t124, t125, -t124, -g(1) * (t134 * pkin(1) + t131 * qJ(2)) - g(2) * (t131 * pkin(1) - t134 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t133 - t124 * t130, g(3) * t130 - t124 * t133, 0, 0, 0, 0, 0, -t124 * t126 - t139, g(3) * t126 - t124 * t127, 0, 0, 0, 0, 0, -g(1) * (t126 * t137 + t136) - g(2) * (-t126 * t135 + t138) - t132 * t139, -g(1) * (-t126 * t138 + t135) - g(2) * (t126 * t136 + t137) + t129 * t139;];
U_reg = t1;
