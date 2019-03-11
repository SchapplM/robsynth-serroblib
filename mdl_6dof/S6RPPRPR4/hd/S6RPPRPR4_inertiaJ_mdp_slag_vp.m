% Calculate joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRPR4_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:10
% EndTime: 2019-03-09 01:47:11
% DurationCPUTime: 0.37s
% Computational Cost: add. (427->108), mult. (689->170), div. (0->0), fcn. (709->8), ass. (0->48)
t97 = sin(qJ(6));
t99 = cos(qJ(6));
t105 = t99 * MDP(24) - t97 * MDP(25);
t121 = cos(qJ(4));
t93 = sin(pkin(10));
t95 = cos(pkin(10));
t98 = sin(qJ(4));
t76 = t93 * t121 + t95 * t98;
t124 = (MDP(21) * t99 - MDP(22) * t97) * t76;
t104 = MDP(24) * t97 + MDP(25) * t99;
t123 = t104 * (t93 * pkin(4) + pkin(8)) - t97 * MDP(21) - t99 * MDP(22);
t122 = t76 ^ 2;
t120 = t76 * t97;
t119 = t76 * t99;
t118 = t97 * t99;
t100 = -pkin(1) - pkin(2);
t94 = sin(pkin(9));
t96 = cos(pkin(9));
t79 = t94 * qJ(2) - t96 * t100;
t81 = t96 * qJ(2) + t94 * t100;
t117 = MDP(18) * pkin(4);
t78 = pkin(3) + t79;
t71 = t121 * pkin(4) + t78;
t116 = t71 * MDP(18);
t113 = -pkin(7) + t81;
t112 = MDP(20) * t118;
t111 = t121 * MDP(15);
t103 = t95 * t121 - t93 * t98;
t110 = t105 * t103;
t109 = (qJ(5) - t113) * t98;
t108 = t121 * t113;
t102 = t98 * MDP(16) - t111;
t92 = t99 ^ 2;
t91 = t97 ^ 2;
t90 = t96 ^ 2;
t84 = -t95 * pkin(4) - pkin(5);
t72 = t103 ^ 2;
t69 = t103 * t94;
t67 = t76 * t94;
t66 = -t121 * qJ(5) + t108;
t64 = t99 * t69 - t97 * t96;
t63 = -t97 * t69 - t99 * t96;
t62 = pkin(5) * t103 + t76 * pkin(8) + t71;
t61 = t93 * t109 + t95 * t66;
t59 = -t95 * t109 + t93 * t66;
t58 = t99 * t61 + t97 * t62;
t57 = -t97 * t61 + t99 * t62;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t79 ^ 2 + t81 ^ 2) * MDP(9) + 0.2e1 * t78 * t111 + (t59 ^ 2 + t61 ^ 2 + t71 ^ 2) * MDP(18) + t72 * MDP(23) - 0.2e1 * t103 * t124 + (t92 * MDP(19) - 0.2e1 * t112) * t122 + (MDP(10) * t98 + 0.2e1 * t121 * MDP(11) - 0.2e1 * t78 * MDP(16)) * t98 + 0.2e1 * t79 * MDP(7) + 0.2e1 * t81 * MDP(8) + 0.2e1 * (-t103 * t61 - t59 * t76) * MDP(17) + 0.2e1 * (t103 * t57 - t59 * t120) * MDP(24) + 0.2e1 * (-t103 * t58 - t59 * t119) * MDP(25); -MDP(4) - pkin(1) * MDP(6) + (-t103 * t69 - t67 * t76) * MDP(17) + (t59 * t67 + t61 * t69) * MDP(18) + (t103 * t63 - t67 * t120) * MDP(24) + (-t103 * t64 - t67 * t119) * MDP(25) + (t81 * MDP(9) + MDP(8)) * t94 + (-t79 * MDP(9) - MDP(7) + t102 - t116) * t96; MDP(6) + (t94 ^ 2 + t90) * MDP(9) + (t67 ^ 2 + t69 ^ 2 + t90) * MDP(18); (-t103 * t59 + t61 * t76) * MDP(18); (-t103 * t67 + t69 * t76) * MDP(18); MDP(9) + (t72 + t122) * MDP(18); -t121 * MDP(13) - MDP(16) * t108 + (-t113 * MDP(15) - MDP(12)) * t98 - t105 * t59 - t123 * t103 + (-MDP(19) * t118 + (t91 - t92) * MDP(20) - t104 * t84) * t76 + ((-t103 * t93 + t76 * t95) * MDP(17) + (-t59 * t95 + t61 * t93) * MDP(18)) * pkin(4); (-t98 * MDP(15) - t121 * MDP(16)) * t94 - t105 * t67 + (-t67 * t95 + t69 * t93) * t117; (t103 * t95 + t76 * t93) * t117 - t102 + t110; 0.2e1 * t112 + t91 * MDP(19) + MDP(14) - 0.2e1 * t105 * t84 + (t93 ^ 2 + t95 ^ 2) * MDP(18) * pkin(4) ^ 2; t110 + t116; -t96 * MDP(18); 0; 0; MDP(18); MDP(23) * t103 + t57 * MDP(24) - t58 * MDP(25) - t124; t63 * MDP(24) - t64 * MDP(25); -t104 * t76; -t123; t105; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
