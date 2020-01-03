% Calculate minimal parameter regressor of potential energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (65->34), mult. (82->38), div. (0->0), fcn. (73->6), ass. (0->19)
t120 = pkin(2) + pkin(5);
t114 = sin(qJ(1));
t115 = cos(qJ(1));
t119 = t115 * pkin(1) + t114 * qJ(2);
t118 = g(1) * t119;
t108 = t114 * pkin(1);
t117 = -t115 * qJ(2) + t108;
t103 = g(1) * t114 - g(2) * t115;
t110 = pkin(7) + qJ(4);
t105 = sin(t110);
t106 = cos(t110);
t111 = sin(pkin(7));
t116 = pkin(3) * t111 + pkin(4) * t105 - qJ(5) * t106;
t113 = -pkin(6) - qJ(3);
t112 = cos(pkin(7));
t104 = g(1) * t115 + g(2) * t114;
t102 = -g(3) * t105 + t103 * t106;
t101 = -g(3) * t106 - t103 * t105;
t1 = [0, -t104, t103, t104, -t103, -g(3) * pkin(5) - g(2) * t117 - t118, -g(3) * t112 - t103 * t111, g(3) * t111 - t103 * t112, -t104, -g(1) * (t115 * qJ(3) + t119) - g(2) * (t114 * qJ(3) + t117) - g(3) * t120, 0, 0, 0, 0, 0, t101, -t102, t101, -t104, t102, -t118 - g(2) * t108 - g(3) * (t112 * pkin(3) + t106 * pkin(4) + t105 * qJ(5) + t120) + (-g(1) * t116 + g(2) * t113) * t114 + (g(1) * t113 - g(2) * (-qJ(2) - t116)) * t115;];
U_reg = t1;
