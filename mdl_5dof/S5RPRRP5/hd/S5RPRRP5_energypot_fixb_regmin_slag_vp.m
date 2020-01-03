% Calculate minimal parameter regressor of potential energy for
% S5RPRRP5
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
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:57
% DurationCPUTime: 0.05s
% Computational Cost: add. (78->26), mult. (56->31), div. (0->0), fcn. (50->8), ass. (0->16)
t121 = qJ(2) + pkin(5);
t113 = qJ(1) + pkin(8);
t112 = qJ(3) + t113;
t110 = sin(t112);
t111 = cos(t112);
t120 = g(1) * t111 + g(2) * t110;
t115 = sin(qJ(1));
t117 = cos(qJ(1));
t119 = -g(1) * t117 - g(2) * t115;
t114 = sin(qJ(4));
t116 = cos(qJ(4));
t118 = pkin(4) * t116 + qJ(5) * t114 + pkin(3);
t109 = g(1) * t110 - g(2) * t111;
t108 = -g(3) * t114 - t120 * t116;
t107 = -g(3) * t116 + t120 * t114;
t1 = [0, t119, g(1) * t115 - g(2) * t117, t119 * pkin(1) - g(3) * t121, 0, -t120, t109, 0, 0, 0, 0, 0, t108, t107, t108, -t109, -t107, -g(1) * (pkin(2) * cos(t113) + t117 * pkin(1)) - g(2) * (pkin(2) * sin(t113) + t115 * pkin(1)) - g(3) * (t114 * pkin(4) - t116 * qJ(5) + pkin(6) + t121) + (g(2) * pkin(7) - g(1) * t118) * t111 + (-g(1) * pkin(7) - g(2) * t118) * t110;];
U_reg = t1;
