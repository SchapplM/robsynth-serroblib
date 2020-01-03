% Calculate Gravitation load on the joints for
% S5RRRPP7
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:05:56
% DurationCPUTime: 0.52s
% Computational Cost: add. (197->80), mult. (457->110), div. (0->0), fcn. (452->6), ass. (0->39)
t132 = MDP(10) - MDP(19) + MDP(24);
t130 = MDP(16) + MDP(18) + MDP(22);
t129 = MDP(17) - MDP(20) - MDP(23);
t100 = sin(qJ(2));
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t111 = g(1) * t104 + g(2) * t101;
t131 = t111 * t100;
t128 = -pkin(3) - pkin(4);
t127 = g(1) * t101;
t94 = t100 * pkin(7);
t103 = cos(qJ(2));
t96 = t103 * pkin(2);
t99 = sin(qJ(3));
t124 = qJ(4) * t99;
t123 = t100 * t99;
t102 = cos(qJ(3));
t122 = t100 * t102;
t121 = t100 * t104;
t120 = t101 * t103;
t119 = t102 * t103;
t118 = t103 * t104;
t117 = -pkin(1) - t96;
t116 = -pkin(2) - t124;
t79 = t102 * t104 + t99 * t120;
t80 = t101 * t119 - t104 * t99;
t115 = -t79 * pkin(3) + qJ(4) * t80;
t81 = -t101 * t102 + t99 * t118;
t82 = t101 * t99 + t102 * t118;
t114 = -t81 * pkin(3) + qJ(4) * t82;
t113 = pkin(3) * t119 + t103 * t124 + t94 + t96;
t112 = -t80 * pkin(3) + t104 * pkin(6) - qJ(4) * t79;
t107 = t104 * pkin(1) + pkin(2) * t118 + t82 * pkin(3) + t101 * pkin(6) + pkin(7) * t121 + t81 * qJ(4);
t106 = g(1) * t81 + g(2) * t79 + g(3) * t123;
t73 = -g(3) * t103 + t131;
t89 = pkin(7) * t118;
t86 = pkin(7) * t120;
t84 = qJ(4) * t122;
t1 = [t111 * MDP(3) + (-g(1) * t112 - g(2) * t107 - (t117 - t94) * t127) * MDP(21) + (-g(1) * (-pkin(4) * t80 + t112) - g(2) * (t82 * pkin(4) - qJ(5) * t121 + t107) - ((-pkin(7) + qJ(5)) * t100 + t117) * t127) * MDP(25) + (t103 * MDP(9) + MDP(2)) * (-g(2) * t104 + t127) - t132 * (-g(2) * t121 + t100 * t127) + t130 * (g(1) * t80 - g(2) * t82) - t129 * (g(1) * t79 - g(2) * t81); (-g(1) * t89 - g(2) * t86 - g(3) * t113 + (pkin(3) * t102 - t116) * t131) * MDP(21) + (-g(1) * (-qJ(5) * t118 + t89) - g(2) * (-qJ(5) * t120 + t86) - g(3) * (pkin(4) * t119 + t113) + (g(3) * qJ(5) + t111 * (-t128 * t102 - t116)) * t100) * MDP(25) + t132 * (g(3) * t100 + t111 * t103) + (t130 * t102 - t129 * t99 + MDP(9)) * t73; (-g(1) * t114 - g(2) * t115 - g(3) * (-pkin(3) * t123 + t84)) * MDP(21) + (-g(1) * (-pkin(4) * t81 + t114) - g(2) * (-pkin(4) * t79 + t115) - g(3) * (t128 * t123 + t84)) * MDP(25) + t130 * t106 + t129 * (g(1) * t82 + g(2) * t80 + g(3) * t122); -(MDP(21) + MDP(25)) * t106; t73 * MDP(25);];
taug = t1;
