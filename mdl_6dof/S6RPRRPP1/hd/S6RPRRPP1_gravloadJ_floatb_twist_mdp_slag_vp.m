% Calculate Gravitation load on the joints for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:52
% EndTime: 2019-03-09 04:29:53
% DurationCPUTime: 0.43s
% Computational Cost: add. (345->90), mult. (361->128), div. (0->0), fcn. (334->10), ass. (0->46)
t101 = cos(qJ(3));
t96 = -qJ(5) - pkin(8);
t98 = sin(qJ(3));
t122 = t98 * t96;
t100 = cos(qJ(4));
t87 = pkin(4) * t100 + pkin(3);
t112 = t101 * t87 - t122;
t129 = MDP(11) - MDP(19) - MDP(22);
t95 = qJ(1) + pkin(9);
t89 = sin(t95);
t91 = cos(t95);
t111 = g(1) * t91 + g(2) * t89;
t69 = -g(3) * t101 + t111 * t98;
t126 = g(3) * t98;
t97 = sin(qJ(4));
t124 = t89 * t97;
t123 = t91 * t97;
t121 = t100 * t89;
t120 = t100 * t91;
t119 = t101 * t89;
t118 = t101 * t91;
t117 = t101 * t96;
t116 = t101 * t97;
t115 = t100 * t101;
t114 = -MDP(20) - MDP(24);
t113 = t91 * t116;
t94 = qJ(4) + pkin(10);
t88 = sin(t94);
t90 = cos(t94);
t108 = pkin(5) * t90 + qJ(6) * t88;
t71 = t89 * t116 + t120;
t65 = t88 * t119 + t90 * t91;
t67 = t88 * t118 - t89 * t90;
t106 = g(1) * t67 + g(2) * t65 + t88 * t126;
t102 = cos(qJ(1));
t105 = t102 * pkin(1) + pkin(4) * t124 + t89 * pkin(7) + t87 * t118 + (pkin(2) - t122) * t91;
t70 = t111 * t101 + t126;
t99 = sin(qJ(1));
t103 = -pkin(1) * t99 + pkin(4) * t123 + t91 * pkin(7) + (-pkin(2) - t112) * t89;
t81 = pkin(4) * t121;
t74 = t91 * t115 + t124;
t73 = -t113 + t121;
t72 = -t89 * t115 + t123;
t68 = t90 * t118 + t88 * t89;
t66 = t90 * t119 - t91 * t88;
t1 = [(g(1) * t102 + g(2) * t99) * MDP(3) + (-g(1) * t72 - g(2) * t74) * MDP(17) + (-g(1) * t71 - g(2) * t73) * MDP(18) + (-g(1) * t103 - g(2) * t105) * MDP(20) + (g(1) * t66 - g(2) * t68) * MDP(21) + (g(1) * t65 - g(2) * t67) * MDP(23) + (-g(1) * (-pkin(5) * t66 - qJ(6) * t65 + t103) - g(2) * (pkin(5) * t68 + qJ(6) * t67 + t105)) * MDP(24) + (t101 * MDP(10) - t129 * t98) * (g(1) * t89 - g(2) * t91) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t99 - g(2) * t102); (-MDP(4) + t114) * g(3); (-g(3) * t112 + t111 * (t87 * t98 + t117)) * MDP(20) + (-g(3) * (t108 * t101 + t112) + t111 * (t117 - (-t108 - t87) * t98)) * MDP(24) + t129 * t70 + (t100 * MDP(17) - t97 * MDP(18) + t90 * MDP(21) + t88 * MDP(23) + MDP(10)) * t69; (-g(1) * t73 + g(2) * t71 + t97 * t126) * MDP(17) + (g(1) * t74 - g(2) * t72 + t100 * t126) * MDP(18) + (-g(1) * t81 + (g(2) * t120 + t70 * t97) * pkin(4)) * MDP(20) + t106 * MDP(21) + (-g(1) * t68 - g(2) * t66 - t90 * t126) * MDP(23) + (-g(1) * (-pkin(4) * t113 - pkin(5) * t67 + qJ(6) * t68 + t81) - g(2) * (-t71 * pkin(4) - pkin(5) * t65 + qJ(6) * t66) - (-pkin(4) * t97 - pkin(5) * t88 + qJ(6) * t90) * t126) * MDP(24); t114 * t69; -t106 * MDP(24);];
taug  = t1;
