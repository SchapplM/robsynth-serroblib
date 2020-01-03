% Calculate minimal parameter regressor of potential energy for
% S5RRRPP3
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:40
% EndTime: 2019-12-31 20:53:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (89->34), mult. (87->36), div. (0->0), fcn. (81->6), ass. (0->18)
t123 = sin(qJ(3));
t134 = qJ(4) * t123 + pkin(2);
t122 = qJ(1) + qJ(2);
t116 = sin(t122);
t125 = cos(qJ(3));
t132 = t116 * t125;
t117 = cos(t122);
t131 = t117 * t125;
t124 = sin(qJ(1));
t130 = t124 * pkin(1) + pkin(3) * t132 + t134 * t116;
t129 = g(1) * t117 + g(2) * t116;
t126 = cos(qJ(1));
t128 = t126 * pkin(1) + pkin(3) * t131 + t116 * pkin(7) + t134 * t117;
t127 = t123 * pkin(3) - t125 * qJ(4) + pkin(5) + pkin(6);
t107 = g(1) * t116 - g(2) * t117;
t106 = g(3) * t123 + t129 * t125;
t105 = -g(3) * t125 + t129 * t123;
t1 = [0, -g(1) * t126 - g(2) * t124, g(1) * t124 - g(2) * t126, 0, -t129, t107, 0, 0, 0, 0, 0, -t106, t105, -t107, t106, -t105, -g(1) * t128 - g(2) * (-t117 * pkin(7) + t130) - g(3) * t127, -t107, -t105, -t106, -g(1) * (t116 * pkin(4) + qJ(5) * t131 + t128) - g(2) * (qJ(5) * t132 + (-pkin(4) - pkin(7)) * t117 + t130) - g(3) * (t123 * qJ(5) + t127);];
U_reg = t1;
