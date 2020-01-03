% Calculate minimal parameter regressor of potential energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:37
% EndTime: 2019-12-31 17:38:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (58->35), mult. (124->54), div. (0->0), fcn. (133->8), ass. (0->23)
t114 = sin(qJ(2));
t129 = qJ(3) * t114 + pkin(1);
t110 = sin(pkin(7));
t112 = cos(pkin(7));
t120 = g(1) * t112 + g(2) * t110;
t109 = sin(pkin(8));
t111 = cos(pkin(8));
t116 = cos(qJ(2));
t118 = t116 * t109 - t114 * t111;
t126 = g(3) * t118;
t124 = t110 * t116;
t123 = t112 * t116;
t122 = pkin(2) * t124 + t129 * t110;
t121 = pkin(2) * t123 + t110 * pkin(5) + t129 * t112;
t119 = t114 * pkin(2) - t116 * qJ(3) + qJ(1);
t117 = t114 * t109 + t116 * t111;
t115 = cos(qJ(5));
t113 = sin(qJ(5));
t99 = t117 * t112;
t98 = t117 * t110;
t97 = -g(3) * t114 - t120 * t116;
t96 = -g(3) * t116 + t120 * t114;
t1 = [-g(3) * qJ(1), 0, t97, t96, t97, -t96, -g(1) * t121 - g(2) * (-t112 * pkin(5) + t122) - g(3) * t119, -g(1) * t99 - g(2) * t98 + t126, g(3) * t117 + t120 * t118, -g(1) * (pkin(3) * t123 - t110 * qJ(4) + t121) - g(2) * (pkin(3) * t124 + (-pkin(5) + qJ(4)) * t112 + t122) - g(3) * (t114 * pkin(3) + t119), 0, 0, 0, 0, 0, -g(1) * (-t110 * t113 + t99 * t115) - g(2) * (t112 * t113 + t98 * t115) + t115 * t126, -g(1) * (-t110 * t115 - t99 * t113) - g(2) * (t112 * t115 - t98 * t113) - t113 * t126;];
U_reg = t1;
