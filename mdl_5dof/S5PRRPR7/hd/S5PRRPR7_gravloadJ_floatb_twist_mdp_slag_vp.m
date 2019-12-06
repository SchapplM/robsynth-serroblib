% Calculate Gravitation load on the joints for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:35
% EndTime: 2019-12-05 16:37:38
% DurationCPUTime: 0.63s
% Computational Cost: add. (230->89), mult. (608->161), div. (0->0), fcn. (737->12), ass. (0->48)
t115 = MDP(11) - MDP(14);
t85 = sin(qJ(2));
t88 = cos(qJ(2));
t100 = cos(pkin(9));
t101 = cos(pkin(5));
t94 = t101 * t100;
t99 = sin(pkin(9));
t69 = t99 * t85 - t88 * t94;
t93 = t101 * t99;
t71 = t100 * t85 + t88 * t93;
t114 = -g(1) * t71 - g(2) * t69;
t81 = sin(pkin(5));
t111 = g(3) * t81;
t80 = sin(pkin(10));
t87 = cos(qJ(3));
t110 = t80 * t87;
t109 = t81 * t85;
t108 = t81 * t88;
t82 = cos(pkin(10));
t83 = sin(qJ(5));
t107 = t82 * t83;
t86 = cos(qJ(5));
t106 = t82 * t86;
t105 = t82 * t87;
t84 = sin(qJ(3));
t104 = t83 * t84;
t103 = t84 * t86;
t102 = t87 * t88;
t98 = t84 * t108;
t97 = t81 * t100;
t96 = t81 * t99;
t70 = t85 * t94 + t99 * t88;
t72 = t100 * t88 - t85 * t93;
t95 = -g(1) * t72 - g(2) * t70;
t64 = t70 * t84 + t87 * t97;
t66 = t72 * t84 - t87 * t96;
t73 = -t101 * t87 + t84 * t109;
t91 = g(1) * t66 + g(2) * t64 + g(3) * t73;
t74 = t101 * t84 + t87 * t109;
t68 = (t82 * t102 + t80 * t85) * t81;
t67 = t72 * t87 + t84 * t96;
t65 = t70 * t87 - t84 * t97;
t63 = -t80 * t108 + t74 * t82;
t62 = -t71 * t105 + t72 * t80;
t61 = -t69 * t105 + t70 * t80;
t59 = t67 * t82 + t71 * t80;
t58 = t65 * t82 + t69 * t80;
t1 = [(-MDP(1) - MDP(15)) * g(3); (g(3) * t109 - t95) * MDP(4) + (-g(1) * t62 - g(2) * t61 - g(3) * t68) * MDP(12) + (-g(1) * (t71 * t110 + t72 * t82) - g(2) * (t69 * t110 + t70 * t82) - (-t80 * t102 + t82 * t85) * t111) * MDP(13) + ((-t85 * t111 + t95) * pkin(7) + (-t88 * t111 - t114) * (pkin(3) * t87 + qJ(4) * t84 + pkin(2))) * MDP(15) + (-g(1) * (-t71 * t104 + t62 * t86) - g(2) * (-t69 * t104 + t61 * t86) - g(3) * (t68 * t86 + t83 * t98)) * MDP(21) + (-g(1) * (-t71 * t103 - t62 * t83) - g(2) * (-t69 * t103 - t61 * t83) - g(3) * (-t68 * t83 + t86 * t98)) * MDP(22) + (-t87 * MDP(10) + t115 * t84 - MDP(3)) * (g(3) * t108 + t114); (-g(1) * (-t66 * pkin(3) + t67 * qJ(4)) - g(2) * (-t64 * pkin(3) + t65 * qJ(4)) - g(3) * (-t73 * pkin(3) + t74 * qJ(4))) * MDP(15) + (-g(1) * (-t66 * t106 + t67 * t83) - g(2) * (-t64 * t106 + t65 * t83) - g(3) * (-t73 * t106 + t74 * t83)) * MDP(21) + (-g(1) * (t66 * t107 + t67 * t86) - g(2) * (t64 * t107 + t65 * t86) - g(3) * (t73 * t107 + t74 * t86)) * MDP(22) + t115 * (g(1) * t67 + g(2) * t65 + g(3) * t74) + (t82 * MDP(12) - MDP(13) * t80 + MDP(10)) * t91; -t91 * MDP(15); (-g(1) * (-t59 * t83 + t66 * t86) - g(2) * (-t58 * t83 + t64 * t86) - g(3) * (-t63 * t83 + t73 * t86)) * MDP(21) + (-g(1) * (-t59 * t86 - t66 * t83) - g(2) * (-t58 * t86 - t64 * t83) - g(3) * (-t63 * t86 - t73 * t83)) * MDP(22);];
taug = t1;
