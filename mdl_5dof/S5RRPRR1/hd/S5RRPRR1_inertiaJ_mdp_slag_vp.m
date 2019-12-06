% Calculate joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:34
% EndTime: 2019-12-05 18:25:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (202->68), mult. (401->102), div. (0->0), fcn. (375->6), ass. (0->43)
t77 = sin(qJ(5));
t80 = cos(qJ(5));
t98 = t77 * MDP(22) + t80 * MDP(23);
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t83 = pkin(1) + pkin(2);
t91 = t80 * MDP(25) - t77 * MDP(26);
t88 = MDP(18) + t91;
t104 = (-MDP(19) * t78 + t88 * t81) * t83;
t82 = cos(qJ(2));
t66 = t83 * t82;
t103 = -0.2e1 * t66;
t79 = sin(qJ(2));
t62 = t78 * t79 - t81 * t82;
t102 = pkin(4) * t62;
t99 = pkin(3) + qJ(3);
t64 = t99 * t79;
t65 = t99 * t82;
t53 = t81 * t64 + t78 * t65;
t50 = t53 * t77;
t101 = t53 * t80;
t100 = t77 * t80;
t96 = t62 * MDP(24);
t63 = t78 * t82 + t81 * t79;
t95 = t63 * MDP(19);
t94 = MDP(21) * t100;
t73 = t77 ^ 2;
t93 = t73 * MDP(20) + MDP(17) + 0.2e1 * t94;
t92 = MDP(22) * t80 - MDP(23) * t77;
t90 = -MDP(25) * t77 - MDP(26) * t80;
t68 = t78 * t83 + pkin(4);
t89 = -t63 * t81 * t83 - t62 * t68;
t54 = -t78 * t64 + t81 * t65;
t75 = t80 ^ 2;
t87 = -t53 * MDP(18) - t54 * MDP(19) + ((-t73 + t75) * MDP(21) + MDP(20) * t100 + MDP(15)) * t63 + (-MDP(16) + t98) * t62;
t85 = pkin(1) ^ 2;
t84 = qJ(3) ^ 2;
t76 = t82 ^ 2;
t74 = t79 ^ 2;
t55 = -t63 * pkin(4) - t66;
t48 = t80 * t54 + t77 * t55;
t47 = -t77 * t54 + t80 * t55;
t1 = [MDP(1) + t74 * MDP(4) + 0.2e1 * t79 * t82 * MDP(5) + (t74 * t84 + (t84 + t85) * t76) * MDP(12) + t95 * t103 + (t75 * MDP(20) + MDP(13) - 0.2e1 * t94) * t63 ^ 2 + (MDP(18) * t103 + t96 + 0.2e1 * (-MDP(14) + t92) * t63) * t62 + 0.2e1 * (t47 * t62 + t63 * t50) * MDP(25) + 0.2e1 * (t63 * t101 - t48 * t62) * MDP(26) + 0.2e1 * (t74 + t76) * MDP(11) * qJ(3); t82 * MDP(7) + (t89 * t77 - t101) * MDP(25) + (t89 * t80 + t50) * MDP(26) + (MDP(6) + (-MDP(12) * qJ(3) - MDP(11)) * pkin(1)) * t79 + t87; t85 * MDP(12) + MDP(8) + 0.2e1 * t104 + t93; -t82 * pkin(1) * MDP(12) + t88 * t62 + t95; 0; MDP(12); (-t77 * t102 - t101) * MDP(25) + (-t80 * t102 + t50) * MDP(26) + t87; t104 + t93; 0; t93; t47 * MDP(25) - t48 * MDP(26) + t92 * t63 + t96; t90 * t68 + t98; t91; t90 * pkin(4) + t98; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
