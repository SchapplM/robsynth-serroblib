% Calculate Gravitation load on the joints for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:23
% EndTime: 2019-12-31 22:12:25
% DurationCPUTime: 0.53s
% Computational Cost: add. (262->94), mult. (633->155), div. (0->0), fcn. (741->10), ass. (0->44)
t84 = sin(qJ(2));
t85 = sin(qJ(1));
t88 = cos(qJ(2));
t109 = cos(qJ(1));
t99 = cos(pkin(5));
t94 = t99 * t109;
t71 = t84 * t94 + t85 * t88;
t83 = sin(qJ(3));
t87 = cos(qJ(3));
t80 = sin(pkin(5));
t97 = t80 * t109;
t63 = t71 * t87 - t83 * t97;
t70 = t85 * t84 - t88 * t94;
t82 = sin(qJ(4));
t86 = cos(qJ(4));
t117 = t63 * t82 - t70 * t86;
t116 = t63 * t86 + t70 * t82;
t114 = MDP(17) - MDP(25);
t112 = g(2) * t70;
t111 = g(2) * t71;
t110 = g(3) * t80;
t106 = t80 * t84;
t105 = t80 * t85;
t104 = t80 * t87;
t103 = t80 * t88;
t102 = t82 * t87;
t101 = t86 * t87;
t100 = t87 * t88;
t98 = pkin(4) * t82 + pkin(8);
t96 = t85 * t99;
t73 = t109 * t88 - t84 * t96;
t67 = t83 * t105 + t73 * t87;
t72 = t109 * t84 + t88 * t96;
t59 = -t67 * t82 + t72 * t86;
t79 = t86 * pkin(4) + pkin(3);
t81 = -qJ(5) - pkin(9);
t93 = t79 * t87 - t81 * t83 + pkin(2);
t62 = t71 * t83 + t87 * t97;
t66 = -t85 * t104 + t73 * t83;
t68 = t83 * t106 - t99 * t87;
t92 = g(1) * t66 + g(2) * t62 + g(3) * t68;
t69 = t84 * t104 + t99 * t83;
t60 = t67 * t86 + t72 * t82;
t1 = [(g(1) * t85 - g(2) * t109) * MDP(2) + (g(1) * t109 + g(2) * t85) * MDP(3) + (g(1) * t71 - g(2) * t73) * MDP(9) + (-g(1) * t70 + g(2) * t72) * MDP(10) + (g(1) * t63 - g(2) * t67) * MDP(16) + (g(1) * t116 - g(2) * t60) * MDP(23) + (-g(1) * t117 - g(2) * t59) * MDP(24) + (-g(1) * (-t85 * pkin(1) - t71 * pkin(2) + pkin(7) * t97 + t62 * t81 - t63 * t79 - t98 * t70) - g(2) * (t109 * pkin(1) + t73 * pkin(2) + pkin(7) * t105 - t66 * t81 + t67 * t79 + t98 * t72)) * MDP(26) + t114 * (-g(1) * t62 + g(2) * t66); (g(1) * t73 + g(3) * t106 + t111) * MDP(10) + (-g(1) * (-t72 * t101 + t73 * t82) - g(2) * (-t70 * t101 + t71 * t82) - (t86 * t100 + t82 * t84) * t110) * MDP(23) + (-g(1) * (t72 * t102 + t73 * t86) - g(2) * (t70 * t102 + t71 * t86) - (-t82 * t100 + t84 * t86) * t110) * MDP(24) + (-g(1) * (-t93 * t72 + t98 * t73) - t98 * t111 + t93 * t112 - (t98 * t84 + t93 * t88) * t110) * MDP(26) + (-t87 * MDP(16) + t114 * t83 - MDP(9)) * (-g(1) * t72 + g(3) * t103 - t112); (-g(1) * (-t66 * t79 - t67 * t81) - g(2) * (-t62 * t79 - t63 * t81) - g(3) * (-t68 * t79 - t69 * t81)) * MDP(26) + t114 * (g(1) * t67 + g(2) * t63 + g(3) * t69) + (MDP(23) * t86 - MDP(24) * t82 + MDP(16)) * t92; (g(1) * t60 + g(2) * t116 - g(3) * (t82 * t103 - t69 * t86)) * MDP(24) + (pkin(4) * MDP(26) + MDP(23)) * (g(2) * t117 - g(3) * (-t86 * t103 - t69 * t82) - g(1) * t59); -t92 * MDP(26);];
taug = t1;
