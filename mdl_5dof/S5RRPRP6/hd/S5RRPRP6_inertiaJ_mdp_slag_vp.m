% Calculate joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP6_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:28
% EndTime: 2019-12-31 19:58:29
% DurationCPUTime: 0.28s
% Computational Cost: add. (370->92), mult. (687->142), div. (0->0), fcn. (688->6), ass. (0->45)
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t75 = sin(qJ(2));
t77 = cos(qJ(2));
t62 = t72 * t77 + t73 * t75;
t98 = 0.2e1 * t62;
t97 = 0.2e1 * t77;
t74 = sin(qJ(4));
t76 = cos(qJ(4));
t96 = t74 * t76;
t94 = -qJ(3) - pkin(6);
t65 = t94 * t77;
t85 = t94 * t75;
t57 = -t73 * t65 + t72 * t85;
t95 = t76 * t57;
t70 = t74 ^ 2;
t71 = t76 ^ 2;
t93 = t70 + t71;
t92 = MDP(21) * pkin(4);
t91 = qJ(5) * t62;
t67 = t72 * pkin(2) + pkin(7);
t90 = qJ(5) + t67;
t89 = MDP(16) * t74;
t61 = t72 * t75 - t73 * t77;
t88 = t61 * MDP(17);
t87 = t74 * MDP(19);
t86 = MDP(14) * t96;
t69 = -t77 * pkin(2) - pkin(1);
t68 = -t73 * pkin(2) - pkin(3);
t54 = t61 * pkin(3) - t62 * pkin(7) + t69;
t50 = t76 * t54 - t74 * t57;
t55 = -t72 * t65 - t73 * t85;
t84 = -MDP(20) * pkin(4) + MDP(15);
t48 = t61 * pkin(4) - t76 * t91 + t50;
t49 = t95 + (t54 - t91) * t74;
t83 = t48 * t76 + t49 * t74;
t58 = t90 * t74;
t59 = t90 * t76;
t82 = -t58 * t76 + t59 * t74;
t81 = t50 * MDP(18) - (t74 * t54 + t95) * MDP(19);
t80 = t76 * MDP(18) - t87;
t79 = MDP(18) * t74 + MDP(19) * t76;
t64 = -t76 * pkin(4) + t68;
t52 = t74 * t62 * pkin(4) + t55;
t1 = [MDP(1) + pkin(1) * MDP(9) * t97 + (t55 ^ 2 + t57 ^ 2 + t69 ^ 2) * MDP(12) + (t48 ^ 2 + t49 ^ 2 + t52 ^ 2) * MDP(21) + (t71 * MDP(13) - 0.2e1 * t86) * t62 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t75 + MDP(5) * t97) * t75 + (t88 + (MDP(15) * t76 - t89) * t98) * t61 + 0.2e1 * (-t57 * MDP(11) + t81) * t61 + (-t83 * MDP(20) + (MDP(11) + t79) * t55) * t98; t75 * MDP(6) + t77 * MDP(7) + (-t48 * t74 + t49 * t76) * MDP(20) + (-t48 * t58 + t49 * t59 + t52 * t64) * MDP(21) - t80 * t55 + (t74 * MDP(15) + t76 * MDP(16) - t79 * t67) * t61 + (-t77 * MDP(10) - t75 * MDP(9)) * pkin(6) + (MDP(13) * t96 + (-t70 + t71) * MDP(14) - t82 * MDP(20) + t79 * t68) * t62 + ((-t61 * t72 - t62 * t73) * MDP(11) + (-t55 * t73 + t57 * t72) * MDP(12)) * pkin(2); MDP(8) + t70 * MDP(13) + 0.2e1 * t86 + 0.2e1 * (t58 * t74 + t59 * t76) * MDP(20) + (t58 ^ 2 + t59 ^ 2 + t64 ^ 2) * MDP(21) + (t72 ^ 2 + t73 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t80 * t68; -t93 * MDP(20) * t62 + t69 * MDP(12) + t83 * MDP(21) + t80 * t61; t82 * MDP(21); t93 * MDP(21) + MDP(12); t48 * t92 + t88 + (t84 * t76 - t89) * t62 + t81; -t58 * t92 + (-MDP(19) * t67 + MDP(16)) * t76 + (-MDP(18) * t67 + t84) * t74; -t87 + (MDP(18) + t92) * t76; MDP(21) * pkin(4) ^ 2 + MDP(17); t52 * MDP(21); t64 * MDP(21); 0; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
