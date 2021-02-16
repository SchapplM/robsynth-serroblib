% Calculate joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP7_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:25
% EndTime: 2021-01-15 16:45:28
% DurationCPUTime: 0.33s
% Computational Cost: add. (277->118), mult. (599->168), div. (0->0), fcn. (569->8), ass. (0->49)
t71 = sin(qJ(3));
t101 = 0.2e1 * t71;
t70 = sin(qJ(4));
t100 = pkin(7) * t70;
t74 = cos(qJ(3));
t99 = pkin(7) * t74;
t69 = cos(pkin(5));
t68 = sin(pkin(5));
t72 = sin(qJ(2));
t97 = t68 * t72;
t56 = -t69 * t74 + t71 * t97;
t98 = t56 * t71;
t75 = cos(qJ(2));
t96 = t68 * t75;
t73 = cos(qJ(4));
t95 = t70 * t73;
t94 = -qJ(5) - pkin(8);
t93 = MDP(22) * pkin(4);
t92 = qJ(5) * t71;
t91 = MDP(11) * t71;
t62 = t94 * t73;
t90 = t62 * MDP(20);
t64 = -t73 * pkin(4) - pkin(3);
t89 = t64 * MDP(22);
t88 = t70 * MDP(15);
t87 = t70 * MDP(20);
t86 = t73 * MDP(19);
t85 = MDP(17) + MDP(19);
t84 = MDP(18) + MDP(20);
t83 = t73 * t99;
t82 = MDP(13) * t95;
t81 = -MDP(21) * pkin(4) + MDP(14);
t80 = MDP(19) + t93;
t79 = MDP(17) * t70 + MDP(18) * t73;
t78 = t70 * MDP(19) + t73 * MDP(20);
t60 = -t74 * pkin(3) - t71 * pkin(8) - pkin(2);
t53 = t83 + (t60 - t92) * t70;
t58 = t73 * t60;
t77 = (-t70 * t99 + t58) * MDP(17) - (t70 * t60 + t83) * MDP(18) - t53 * MDP(20);
t76 = -t86 + t87 + t89;
t67 = t73 ^ 2;
t65 = t70 ^ 2;
t61 = t94 * t70;
t59 = (pkin(4) * t70 + pkin(7)) * t71;
t57 = t69 * t71 + t74 * t97;
t52 = t57 * t73 - t70 * t96;
t51 = -t57 * t70 - t73 * t96;
t50 = -t73 * t92 + t58 + (-pkin(4) - t100) * t74;
t1 = [MDP(1) + (t51 ^ 2 + t52 ^ 2 + t56 ^ 2) * MDP(22); (t51 * t50 + t52 * t53 + t56 * t59) * MDP(22) + (-t51 * t73 - t52 * t70) * MDP(21) * t71 + t84 * (t52 * t74 + t73 * t98) + t85 * (-t51 * t74 + t70 * t98) + (-t72 * MDP(4) + (MDP(10) * t74 + MDP(3) - t91) * t75) * t68; MDP(2) - 0.2e1 * pkin(2) * t91 + (t50 ^ 2 + t53 ^ 2 + t59 ^ 2) * MDP(22) + (0.2e1 * pkin(2) * MDP(10) + t74 * MDP(16) + (-t73 * MDP(14) + MDP(6) + t88) * t101) * t74 + 0.2e1 * (-t50 * MDP(19) - t77) * t74 + ((-t50 * t73 - t53 * t70) * MDP(21) + t78 * t59) * t101 + (t67 * MDP(12) + 0.2e1 * t79 * pkin(7) + MDP(5) - 0.2e1 * t82) * t71 ^ 2; -t57 * MDP(11) + (-t51 * t70 + t52 * t73) * MDP(21) + (t51 * t61 - t52 * t62) * MDP(22) + (t84 * t70 - t73 * t85 - MDP(10) + t89) * t56; (-t50 * t70 + t53 * t73) * MDP(21) + (t50 * t61 - t53 * t62) * MDP(22) + t76 * t59 + (-pkin(7) * MDP(11) - t70 * MDP(14) - t73 * MDP(15) - t61 * MDP(19) + pkin(8) * t79 + MDP(8) - t90) * t74 + (MDP(7) - pkin(7) * MDP(10) + MDP(12) * t95 + (-t65 + t67) * MDP(13) + (-pkin(3) * t70 - pkin(7) * t73) * MDP(17) + (-pkin(3) * t73 + t100) * MDP(18) + (-t61 * t73 + t62 * t70) * MDP(21) + t78 * t64) * t71; MDP(9) + t65 * MDP(12) + 0.2e1 * t82 + 0.2e1 * (-t61 * t70 - t62 * t73) * MDP(21) + (t61 ^ 2 + t62 ^ 2) * MDP(22) + (-0.2e1 * t86 + 0.2e1 * t87 + t89) * t64 + 0.2e1 * (t73 * MDP(17) - t70 * MDP(18)) * pkin(3); -t84 * t52 + (MDP(17) + t80) * t51; t50 * t93 + t58 * MDP(19) + (-MDP(16) + (-0.2e1 * pkin(4) - t100) * MDP(19)) * t74 + (-t88 + (-MDP(19) * qJ(5) + t81) * t73) * t71 + t77; t90 + (-MDP(18) * pkin(8) + MDP(15)) * t73 + t80 * t61 + (-MDP(17) * pkin(8) + t81) * t70; MDP(16) + (0.2e1 * MDP(19) + t93) * pkin(4); t56 * MDP(22); t59 * MDP(22) + t71 * t78; t76; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
