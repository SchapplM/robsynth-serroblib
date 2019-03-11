% Calculate Gravitation load on the joints for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:29
% EndTime: 2019-03-08 22:14:31
% DurationCPUTime: 0.62s
% Computational Cost: add. (327->83), mult. (860->141), div. (0->0), fcn. (1059->12), ass. (0->44)
t101 = cos(qJ(6));
t102 = cos(qJ(5));
t103 = cos(qJ(3));
t95 = sin(pkin(6));
t124 = t103 * t95;
t100 = sin(qJ(2));
t104 = cos(qJ(2));
t122 = cos(pkin(6));
t96 = cos(pkin(11));
t116 = t96 * t122;
t121 = sin(pkin(11));
t86 = t100 * t116 + t121 * t104;
t99 = sin(qJ(3));
t80 = t96 * t124 + t86 * t99;
t81 = -t96 * t95 * t99 + t103 * t86;
t98 = sin(qJ(5));
t69 = t102 * t81 + t80 * t98;
t117 = t95 * t121;
t114 = t122 * t121;
t88 = -t100 * t114 + t104 * t96;
t82 = -t103 * t117 + t88 * t99;
t83 = t88 * t103 + t99 * t117;
t72 = t102 * t83 + t82 * t98;
t125 = t100 * t95;
t89 = -t122 * t103 + t99 * t125;
t90 = t100 * t124 + t122 * t99;
t79 = t102 * t90 + t89 * t98;
t97 = sin(qJ(6));
t136 = (g(1) * t72 + g(2) * t69 + g(3) * t79) * MDP(22) + (-MDP(28) * t101 + MDP(29) * t97 - MDP(21)) * (g(1) * (t102 * t82 - t83 * t98) + g(2) * (t102 * t80 - t81 * t98) + g(3) * (t102 * t89 - t90 * t98));
t135 = MDP(10) + MDP(12);
t134 = MDP(11) - MDP(14);
t87 = -t96 * t100 - t104 * t114;
t129 = g(1) * t87;
t85 = -t121 * t100 + t104 * t116;
t127 = g(2) * t85;
t123 = t104 * t95;
t115 = -g(1) * t88 - g(2) * t86;
t112 = t102 * t103 + t98 * t99;
t111 = pkin(3) * t103 + qJ(4) * t99 + pkin(2);
t108 = g(1) * t82 + g(2) * t80 + g(3) * t89;
t84 = t112 * t123;
t76 = t112 * t87;
t75 = t112 * t85;
t1 = [(-MDP(1) - MDP(15)) * g(3); (t115 * pkin(8) - t111 * t129 - t111 * t127 - g(3) * (pkin(8) * t100 + t111 * t104) * t95) * MDP(15) + (-g(1) * t76 - g(2) * t75 - g(3) * t84) * MDP(21) + (-g(1) * (t101 * t76 - t88 * t97) - g(2) * (t101 * t75 - t86 * t97) - g(3) * (t101 * t84 - t97 * t125)) * MDP(28) + (-g(1) * (-t101 * t88 - t76 * t97) - g(2) * (-t101 * t86 - t75 * t97) - g(3) * (-t101 * t125 - t84 * t97)) * MDP(29) + (MDP(4) - MDP(13)) * (g(3) * t125 - t115) + (-MDP(3) - t135 * t103 + t134 * t99 - (t102 * t99 - t103 * t98) * MDP(22)) * (g(3) * t123 + t127 + t129); (-g(1) * (-pkin(3) * t82 + qJ(4) * t83) - g(2) * (-pkin(3) * t80 + qJ(4) * t81) - g(3) * (-pkin(3) * t89 + qJ(4) * t90)) * MDP(15) + t135 * t108 + t134 * (g(1) * t83 + g(2) * t81 + g(3) * t90) - t136; -t108 * MDP(15); t136; (-g(1) * (t101 * t87 - t72 * t97) - g(2) * (t101 * t85 - t69 * t97) - g(3) * (t101 * t123 - t79 * t97)) * MDP(28) + (-g(1) * (-t101 * t72 - t87 * t97) - g(2) * (-t101 * t69 - t85 * t97) - g(3) * (-t79 * t101 - t97 * t123)) * MDP(29);];
taug  = t1;
