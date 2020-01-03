% Calculate minimal parameter regressor of potential energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:41
% EndTime: 2020-01-03 11:47:42
% DurationCPUTime: 0.04s
% Computational Cost: add. (54->23), mult. (46->31), div. (0->0), fcn. (40->8), ass. (0->16)
t130 = qJ(2) + pkin(5);
t122 = qJ(1) + pkin(8);
t117 = sin(t122);
t118 = cos(t122);
t129 = g(2) * t117 - g(3) * t118;
t125 = sin(qJ(1));
t127 = cos(qJ(1));
t128 = -g(2) * t125 + g(3) * t127;
t126 = cos(qJ(3));
t124 = sin(qJ(3));
t123 = qJ(3) + qJ(4);
t121 = -qJ(5) - pkin(7) - pkin(6);
t120 = cos(t123);
t119 = sin(t123);
t116 = t126 * pkin(3) + pkin(4) * t120 + pkin(2);
t1 = [0, t128, -g(2) * t127 - g(3) * t125, t128 * pkin(1) - g(1) * t130, 0, 0, 0, 0, 0, -g(1) * t124 - t129 * t126, -g(1) * t126 + t129 * t124, 0, 0, 0, 0, 0, -g(1) * t119 - t129 * t120, -g(1) * t120 + t129 * t119, g(2) * t118 + g(3) * t117, -g(1) * (t124 * pkin(3) + pkin(4) * t119 + t130) - g(2) * (t125 * pkin(1) + t117 * t116 + t118 * t121) - g(3) * (-t127 * pkin(1) - t118 * t116 + t117 * t121);];
U_reg = t1;
