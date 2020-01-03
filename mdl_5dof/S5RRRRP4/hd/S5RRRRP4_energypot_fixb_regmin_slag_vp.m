% Calculate minimal parameter regressor of potential energy for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:12
% EndTime: 2019-12-31 21:51:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (79->26), mult. (62->30), div. (0->0), fcn. (59->8), ass. (0->18)
t127 = qJ(1) + qJ(2);
t123 = sin(t127);
t125 = cos(t127);
t135 = g(1) * t125 + g(2) * t123;
t129 = sin(qJ(1));
t131 = cos(qJ(1));
t134 = -g(1) * t131 - g(2) * t129;
t126 = qJ(3) + qJ(4);
t122 = sin(t126);
t124 = cos(t126);
t130 = cos(qJ(3));
t133 = t130 * pkin(3) + pkin(4) * t124 + qJ(5) * t122 + pkin(2);
t132 = -pkin(8) - pkin(7);
t128 = sin(qJ(3));
t120 = g(1) * t123 - g(2) * t125;
t119 = -g(3) * t122 - t135 * t124;
t118 = -g(3) * t124 + t135 * t122;
t1 = [0, t134, g(1) * t129 - g(2) * t131, 0, -t135, t120, 0, 0, 0, 0, 0, -g(3) * t128 - t135 * t130, -g(3) * t130 + t135 * t128, 0, 0, 0, 0, 0, t119, t118, t119, -t120, -t118, -g(3) * (t128 * pkin(3) + t122 * pkin(4) - t124 * qJ(5) + pkin(5) + pkin(6)) + t134 * pkin(1) + (-g(1) * t133 - g(2) * t132) * t125 + (g(1) * t132 - g(2) * t133) * t123;];
U_reg = t1;
