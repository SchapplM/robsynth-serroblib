% Calculate minimal parameter regressor of potential energy for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:13
% EndTime: 2019-12-31 19:51:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (93->34), mult. (73->40), div. (0->0), fcn. (67->8), ass. (0->20)
t133 = pkin(6) + pkin(5);
t125 = qJ(1) + qJ(2);
t120 = sin(t125);
t121 = cos(t125);
t132 = g(1) * t121 + g(2) * t120;
t124 = pkin(8) + qJ(4);
t118 = sin(t124);
t119 = cos(t124);
t127 = cos(pkin(8));
t131 = t127 * pkin(3) + pkin(4) * t119 + qJ(5) * t118 + pkin(2);
t130 = cos(qJ(1));
t129 = sin(qJ(1));
t128 = -pkin(7) - qJ(3);
t126 = sin(pkin(8));
t123 = t130 * pkin(1);
t122 = t129 * pkin(1);
t115 = g(1) * t120 - g(2) * t121;
t114 = -g(3) * t118 - t132 * t119;
t113 = -g(3) * t119 + t132 * t118;
t1 = [0, -g(1) * t130 - g(2) * t129, g(1) * t129 - g(2) * t130, 0, -t132, t115, -g(3) * t126 - t132 * t127, -g(3) * t127 + t132 * t126, -t115, -g(1) * (t121 * pkin(2) + t120 * qJ(3) + t123) - g(2) * (t120 * pkin(2) - t121 * qJ(3) + t122) - g(3) * t133, 0, 0, 0, 0, 0, t114, t113, t114, -t115, -t113, -g(1) * t123 - g(2) * t122 - g(3) * (t126 * pkin(3) + t118 * pkin(4) - t119 * qJ(5) + t133) + (-g(1) * t131 - g(2) * t128) * t121 + (g(1) * t128 - g(2) * t131) * t120;];
U_reg = t1;
