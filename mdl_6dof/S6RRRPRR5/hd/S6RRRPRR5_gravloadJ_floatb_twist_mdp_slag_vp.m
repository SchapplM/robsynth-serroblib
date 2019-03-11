% Calculate Gravitation load on the joints for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:10
% EndTime: 2019-03-09 18:23:11
% DurationCPUTime: 0.26s
% Computational Cost: add. (283->70), mult. (329->98), div. (0->0), fcn. (300->10), ass. (0->45)
t93 = qJ(2) + qJ(3);
t88 = sin(t93);
t90 = cos(t93);
t107 = t90 * pkin(3) + t88 * qJ(4);
t98 = cos(qJ(2));
t122 = t98 * pkin(2) + t107;
t121 = MDP(16) - MDP(19);
t120 = MDP(17) - MDP(20);
t119 = pkin(3) * t88;
t117 = g(3) * t90;
t92 = qJ(5) + qJ(6);
t87 = sin(t92);
t96 = sin(qJ(1));
t116 = t96 * t87;
t89 = cos(t92);
t115 = t96 * t89;
t94 = sin(qJ(5));
t114 = t96 * t94;
t97 = cos(qJ(5));
t113 = t96 * t97;
t99 = cos(qJ(1));
t112 = t99 * t87;
t111 = t99 * t89;
t110 = t99 * t94;
t109 = t99 * t97;
t70 = t88 * t111 - t116;
t71 = t88 * t112 + t115;
t72 = t88 * t115 + t112;
t73 = -t88 * t116 + t111;
t108 = (-g(1) * t70 - g(2) * t72 + t89 * t117) * MDP(34) + (g(1) * t71 - g(2) * t73 - t87 * t117) * MDP(35);
t106 = qJ(4) * t90;
t95 = sin(qJ(2));
t105 = -pkin(2) * t95 - t119;
t80 = g(1) * t99 + g(2) * t96;
t103 = pkin(1) + t122;
t68 = t80 * t88 - t117;
t102 = t121 * t68 + (-MDP(27) * t94 - MDP(28) * t97 - MDP(34) * t87 - MDP(35) * t89 + t120) * (g(3) * t88 + t80 * t90);
t100 = -pkin(8) - pkin(7);
t82 = t99 * t106;
t81 = t96 * t106;
t79 = -t88 * t114 + t109;
t78 = t88 * t113 + t110;
t77 = t88 * t110 + t113;
t76 = t88 * t109 - t114;
t1 = [((g(1) * t100 - g(2) * t103) * t99 + (g(1) * t103 + g(2) * t100) * t96) * MDP(21) + (-g(1) * t79 - g(2) * t77) * MDP(27) + (g(1) * t78 - g(2) * t76) * MDP(28) + (-g(1) * t73 - g(2) * t71) * MDP(34) + (g(1) * t72 - g(2) * t70) * MDP(35) + (MDP(3) - MDP(18)) * t80 + (-t95 * MDP(10) + t98 * MDP(9) - t120 * t88 + t121 * t90 + MDP(2)) * (g(1) * t96 - g(2) * t99); (-g(3) * t98 + t80 * t95) * MDP(9) + (g(3) * t95 + t80 * t98) * MDP(10) + (-g(1) * (t105 * t99 + t82) - g(2) * (t105 * t96 + t81) - g(3) * t122) * MDP(21) + t102; (-g(1) * (-t99 * t119 + t82) - g(2) * (-t96 * t119 + t81) - g(3) * t107) * MDP(21) + t102; -t68 * MDP(21); (-g(1) * t76 - g(2) * t78 + t97 * t117) * MDP(27) + (g(1) * t77 - g(2) * t79 - t94 * t117) * MDP(28) + t108; t108;];
taug  = t1;
