% Calculate Gravitation load on the joints for
% S5RRRPP8
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
%   see S5RRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:25
% EndTime: 2019-12-31 21:09:27
% DurationCPUTime: 0.60s
% Computational Cost: add. (199->80), mult. (462->108), div. (0->0), fcn. (459->6), ass. (0->39)
t133 = MDP(17) - MDP(20) - MDP(23);
t131 = MDP(16) - MDP(19) + MDP(24);
t130 = MDP(10) - MDP(18) - MDP(22);
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t108 = g(1) * t101 + g(2) * t98;
t97 = sin(qJ(2));
t129 = t108 * t97;
t126 = g(1) * t98;
t91 = t97 * pkin(7);
t100 = cos(qJ(2));
t93 = t100 * pkin(2);
t96 = sin(qJ(3));
t123 = t96 * t97;
t99 = cos(qJ(3));
t122 = t97 * t99;
t121 = -pkin(3) - qJ(5);
t120 = qJ(4) * t96;
t119 = t100 * t98;
t118 = t100 * t99;
t117 = t101 * t97;
t116 = t100 * t101;
t115 = -pkin(1) - t93;
t114 = -pkin(2) - t120;
t77 = t101 * t99 + t96 * t119;
t78 = -t101 * t96 + t98 * t118;
t113 = -t77 * pkin(3) + qJ(4) * t78;
t79 = t96 * t116 - t98 * t99;
t80 = t99 * t116 + t98 * t96;
t112 = -t79 * pkin(3) + qJ(4) * t80;
t111 = pkin(3) * t118 + t100 * t120 + t91 + t93;
t110 = -t78 * pkin(3) + t101 * pkin(6) - qJ(4) * t77;
t106 = t101 * pkin(1) + pkin(2) * t116 + t80 * pkin(3) + t98 * pkin(6) + pkin(7) * t117 + t79 * qJ(4);
t104 = g(1) * t79 + g(2) * t77 + g(3) * t123;
t103 = g(1) * t80 + g(2) * t78 + g(3) * t122;
t87 = pkin(7) * t116;
t84 = pkin(7) * t119;
t82 = qJ(4) * t122;
t1 = [t108 * MDP(3) + (-g(1) * t110 - g(2) * t106 - (t115 - t91) * t126) * MDP(21) + (-g(1) * (-qJ(5) * t78 + t110) - g(2) * (pkin(4) * t117 + t80 * qJ(5) + t106) - ((-pkin(4) - pkin(7)) * t97 + t115) * t126) * MDP(25) + t131 * (g(1) * t78 - g(2) * t80) - t133 * (g(1) * t77 - g(2) * t79) + (t100 * MDP(9) - t130 * t97 + MDP(2)) * (-g(2) * t101 + t126); (-g(1) * t87 - g(2) * t84 - g(3) * t111 + (pkin(3) * t99 - t114) * t129) * MDP(21) + (-g(1) * (pkin(4) * t116 + t87) - g(2) * (pkin(4) * t119 + t84) - g(3) * (qJ(5) * t118 + t111) + (-g(3) * pkin(4) + t108 * (-t121 * t99 - t114)) * t97) * MDP(25) + t130 * (g(3) * t97 + t108 * t100) + (t131 * t99 - t133 * t96 + MDP(9)) * (-g(3) * t100 + t129); (-g(1) * t112 - g(2) * t113 - g(3) * (-pkin(3) * t123 + t82)) * MDP(21) + (-g(1) * (-qJ(5) * t79 + t112) - g(2) * (-qJ(5) * t77 + t113) - g(3) * (t121 * t123 + t82)) * MDP(25) + t131 * t104 + t133 * t103; -(MDP(21) + MDP(25)) * t104; -t103 * MDP(25);];
taug = t1;
