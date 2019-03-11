% Calculate Gravitation load on the joints for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:07
% EndTime: 2019-03-09 15:23:08
% DurationCPUTime: 0.38s
% Computational Cost: add. (341->82), mult. (316->114), div. (0->0), fcn. (273->12), ass. (0->50)
t96 = sin(qJ(1));
t98 = cos(qJ(1));
t77 = g(1) * t98 + g(2) * t96;
t92 = qJ(2) + qJ(3);
t86 = sin(t92);
t120 = t77 * t86;
t85 = pkin(10) + t92;
t80 = sin(t85);
t81 = cos(t85);
t102 = t81 * pkin(4) + t80 * qJ(5);
t119 = MDP(19) + MDP(23);
t100 = -g(3) * t81 + t77 * t80;
t118 = pkin(4) * t80;
t117 = g(3) * t80;
t87 = cos(t92);
t115 = g(3) * t87;
t91 = pkin(11) + qJ(6);
t83 = sin(t91);
t114 = t96 * t83;
t84 = cos(t91);
t113 = t96 * t84;
t93 = sin(pkin(11));
t112 = t96 * t93;
t94 = cos(pkin(11));
t111 = t96 * t94;
t110 = t98 * t83;
t109 = t98 * t84;
t108 = t98 * t93;
t107 = t98 * t94;
t82 = pkin(3) * t87;
t97 = cos(qJ(2));
t88 = t97 * pkin(2);
t106 = t82 + t88;
t105 = qJ(5) * t81;
t104 = t82 + t102;
t95 = sin(qJ(2));
t73 = -t95 * pkin(2) - pkin(3) * t86;
t103 = t73 - t118;
t76 = g(1) * t96 - g(2) * t98;
t101 = (-t77 * t81 - t117) * MDP(22) + (-t115 + t120) * MDP(16) + (g(3) * t86 + t77 * t87) * MDP(17) + (t94 * MDP(20) - t93 * MDP(21) + t84 * MDP(29) - t83 * MDP(30)) * t100;
t90 = -qJ(4) - pkin(8) - pkin(7);
t75 = t98 * t105;
t74 = t96 * t105;
t72 = pkin(1) + t106;
t71 = t98 * t72;
t70 = t81 * t109 + t114;
t69 = -t81 * t110 + t113;
t68 = -t81 * t113 + t110;
t67 = t81 * t114 + t109;
t1 = [(-g(1) * (-t96 * t72 - t98 * t90) - g(2) * (-t96 * t90 + t71)) * MDP(19) + (-g(1) * (-t81 * t111 + t108) - g(2) * (t81 * t107 + t112)) * MDP(20) + (-g(1) * (t81 * t112 + t107) - g(2) * (-t81 * t108 + t111)) * MDP(21) + (-g(2) * t71 + (g(1) * t90 - g(2) * t102) * t98 + (-g(1) * (-t102 - t72) + g(2) * t90) * t96) * MDP(23) + (-g(1) * t68 - g(2) * t70) * MDP(29) + (-g(1) * t67 - g(2) * t69) * MDP(30) + (MDP(3) - MDP(18)) * t77 + (-t95 * MDP(10) + MDP(16) * t87 - MDP(17) * t86 + t80 * MDP(22) + t97 * MDP(9) + MDP(2)) * t76; (-g(3) * t97 + t77 * t95) * MDP(9) + (g(3) * t95 + t77 * t97) * MDP(10) + (-g(3) * t106 - t77 * t73) * MDP(19) + (-g(1) * (t103 * t98 + t75) - g(2) * (t103 * t96 + t74) - g(3) * (t88 + t104)) * MDP(23) + t101; (-g(1) * (-t98 * t118 + t75) - g(2) * (-t96 * t118 + t74) - g(3) * t104) * MDP(23) + (-MDP(19) * t115 + t119 * t120) * pkin(3) + t101; -t119 * t76; -t100 * MDP(23); (-g(1) * t69 + g(2) * t67 + t83 * t117) * MDP(29) + (g(1) * t70 - g(2) * t68 + t84 * t117) * MDP(30);];
taug  = t1;
