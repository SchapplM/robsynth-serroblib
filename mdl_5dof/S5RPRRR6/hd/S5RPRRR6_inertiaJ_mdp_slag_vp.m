% Calculate joint inertia matrix for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR6_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:41
% EndTime: 2019-12-31 19:01:42
% DurationCPUTime: 0.22s
% Computational Cost: add. (261->76), mult. (487->111), div. (0->0), fcn. (481->8), ass. (0->44)
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t105 = t85 * MDP(21) + t88 * MDP(22);
t95 = t88 * MDP(24) - t85 * MDP(25);
t84 = cos(pkin(9));
t74 = -t84 * pkin(1) - pkin(2);
t89 = cos(qJ(3));
t71 = -t89 * pkin(3) + t74;
t111 = 0.2e1 * t71;
t83 = sin(pkin(9));
t73 = t83 * pkin(1) + pkin(6);
t109 = pkin(7) + t73;
t108 = cos(qJ(4));
t87 = sin(qJ(3));
t66 = t109 * t87;
t67 = t109 * t89;
t86 = sin(qJ(4));
t55 = t108 * t66 + t86 * t67;
t51 = t55 * t85;
t107 = t55 * t88;
t106 = t85 * t88;
t69 = -t108 * t89 + t86 * t87;
t104 = t69 * MDP(23);
t70 = t108 * t87 + t86 * t89;
t65 = t70 * MDP(18);
t101 = t89 * MDP(10);
t100 = MDP(20) * t106;
t81 = t85 ^ 2;
t99 = t81 * MDP(19) + MDP(16) + 0.2e1 * t100;
t98 = -pkin(4) * t70 - pkin(8) * t69;
t76 = t86 * pkin(3) + pkin(8);
t77 = -t108 * pkin(3) - pkin(4);
t97 = -t69 * t76 + t70 * t77;
t96 = MDP(21) * t88 - MDP(22) * t85;
t94 = -MDP(24) * t85 - MDP(25) * t88;
t93 = -t65 + (-MDP(17) - t95) * t69;
t56 = t108 * t67 - t86 * t66;
t82 = t88 ^ 2;
t92 = -t55 * MDP(17) - t56 * MDP(18) + ((-t81 + t82) * MDP(20) + MDP(19) * t106 + MDP(14)) * t70 + (-MDP(15) + t105) * t69;
t91 = (t108 * MDP(17) - t86 * MDP(18)) * pkin(3);
t52 = t69 * pkin(4) - t70 * pkin(8) + t71;
t50 = t85 * t52 + t88 * t56;
t49 = t88 * t52 - t85 * t56;
t1 = [-0.2e1 * t74 * t101 + t65 * t111 + MDP(1) + (t83 ^ 2 + t84 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(17) * t111 + t104 + 0.2e1 * (-MDP(13) + t96) * t70) * t69 + 0.2e1 * (t49 * t69 + t70 * t51) * MDP(24) + 0.2e1 * (t70 * t107 - t50 * t69) * MDP(25) + (0.2e1 * t74 * MDP(11) + MDP(5) * t87 + 0.2e1 * t89 * MDP(6)) * t87 + (t82 * MDP(19) + MDP(12) - 0.2e1 * t100) * t70 ^ 2; 0; MDP(4); t87 * MDP(7) + t89 * MDP(8) + (t97 * t85 - t107) * MDP(24) + (t97 * t88 + t51) * MDP(25) + (-t87 * MDP(10) - t89 * MDP(11)) * t73 + t92; -t87 * MDP(11) + t101 + t93; -0.2e1 * t77 * t95 + MDP(9) + 0.2e1 * t91 + t99; (t98 * t85 - t107) * MDP(24) + (t98 * t88 + t51) * MDP(25) + t92; t93; t91 + t99 + t95 * (pkin(4) - t77); 0.2e1 * pkin(4) * t95 + t99; t49 * MDP(24) - t50 * MDP(25) + t96 * t70 + t104; t94 * t70; t94 * t76 + t105; t94 * pkin(8) + t105; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
