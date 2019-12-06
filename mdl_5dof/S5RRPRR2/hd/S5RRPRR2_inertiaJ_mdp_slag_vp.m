% Calculate joint inertia matrix for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:28
% EndTime: 2019-12-05 18:28:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (510->82), mult. (936->120), div. (0->0), fcn. (1071->8), ass. (0->48)
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t82 = -t92 * t96 + t93 * t99;
t83 = t92 * t99 + t93 * t96;
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t72 = -t98 * t82 + t95 * t83;
t91 = -t99 * pkin(2) - pkin(1);
t77 = -t82 * pkin(3) + t91;
t118 = 0.2e1 * t72 * pkin(4) + 0.2e1 * t77;
t117 = 0.2e1 * t77;
t116 = 0.2e1 * t99;
t115 = pkin(2) * t92;
t90 = t93 * pkin(2) + pkin(3);
t80 = t98 * t115 + t95 * t90;
t97 = cos(qJ(5));
t114 = t97 * t80;
t113 = -qJ(3) - pkin(6);
t87 = t113 * t96;
t88 = t113 * t99;
t75 = t92 * t87 - t93 * t88;
t73 = t95 * t82 + t98 * t83;
t94 = sin(qJ(5));
t61 = t97 * t72 + t94 * t73;
t112 = t61 * MDP(25);
t79 = -t95 * t115 + t98 * t90;
t78 = pkin(4) + t79;
t65 = t97 * t78 - t94 * t80;
t111 = t65 * MDP(25);
t110 = (-t94 * t78 - t114) * MDP(26);
t109 = t72 * MDP(18);
t108 = t79 * MDP(18);
t107 = t80 * MDP(19);
t106 = MDP(17) + MDP(24);
t74 = t93 * t87 + t92 * t88;
t67 = -t83 * pkin(7) + t74;
t68 = t82 * pkin(7) + t75;
t105 = t98 * t67 - t95 * t68;
t55 = -t73 * pkin(8) + t105;
t103 = -t95 * t67 - t98 * t68;
t56 = -t72 * pkin(8) - t103;
t62 = -t94 * t72 + t97 * t73;
t104 = t62 * MDP(22) - t61 * MDP(23) + (t97 * t55 - t94 * t56) * MDP(25) + (-t94 * t55 - t97 * t56) * MDP(26);
t102 = (t97 * MDP(25) - t94 * MDP(26)) * pkin(4);
t101 = t73 * MDP(15) - t72 * MDP(16) + t105 * MDP(18) + t103 * MDP(19) + t104;
t1 = [MDP(1) + pkin(1) * MDP(9) * t116 + 0.2e1 * (-t74 * t83 + t75 * t82) * MDP(11) + (t74 ^ 2 + t75 ^ 2 + t91 ^ 2) * MDP(12) + t109 * t117 + t112 * t118 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t96 + MDP(5) * t116) * t96 + (MDP(13) * t73 - 0.2e1 * t72 * MDP(14) + MDP(19) * t117) * t73 + (MDP(20) * t62 - 0.2e1 * t61 * MDP(21) + MDP(26) * t118) * t62; t96 * MDP(6) + t99 * MDP(7) + (-t99 * MDP(10) - t96 * MDP(9)) * pkin(6) + ((t82 * t92 - t83 * t93) * MDP(11) + (t74 * t93 + t75 * t92) * MDP(12)) * pkin(2) + t101; MDP(8) + (t92 ^ 2 + t93 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t108 - 0.2e1 * t107 + 0.2e1 * t111 + 0.2e1 * t110 + t106; t91 * MDP(12) + t73 * MDP(19) + t62 * MDP(26) + t109 + t112; 0; MDP(12); t101; t108 - t107 + (t97 * pkin(4) + t65) * MDP(25) + (-t114 + (-pkin(4) - t78) * t94) * MDP(26) + t106; 0; 0.2e1 * t102 + t106; t104; MDP(24) + t110 + t111; 0; MDP(24) + t102; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
